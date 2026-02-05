import gradio as gr
from lasso_workbench.utils.sequence_io import write_fasta_entries
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any, Optional, Set
import json
import tempfile
from datetime import datetime, timezone

from lasso_workbench.services.pipeline_service import PipelineService, BGCSegment
from lasso_workbench.services.rules_service import RulesService
from lasso_workbench.core.prediction import CorePredictor, RuleEngine
from lasso_workbench import DATA_DIR, __version__
from lasso_workbench.schemas.pipeline import RankingConfig
from lasso_workbench.visualization.pipeline_viz import generate_gradio_iframe_html

# Constants
DEFAULT_PRECURSORS_PATH = DATA_DIR / "precursors" / "precursor_proteins_verified.faa"
DEFAULT_PRECURSOR_SET_NAME = "precursor_proteins_verified.faa"
UPLOAD_CUSTOM_LABEL = "Upload Custom..."


def get_available_precursor_sets() -> List[str]:
    """Get list of available precursor reference sets."""
    sets: List[str] = []
    precursors_dir = DATA_DIR / "precursors"
    if precursors_dir.exists():
        for faa_file in sorted(precursors_dir.glob("*.faa")):
            sets.append(faa_file.name)
    if DEFAULT_PRECURSOR_SET_NAME not in sets:
        sets.insert(0, DEFAULT_PRECURSOR_SET_NAME)
    sets.append(UPLOAD_CUSTOM_LABEL)
    return sets


def create_pipeline_tab(pipeline_service: PipelineService, rules_service: RulesService):
    """Create the Semantic Pipeline tab with visualization and batch core prediction."""

    def get_configured_predictor() -> CorePredictor:
        """Get a CorePredictor configured with current rules."""
        rules = rules_service.get_rules()
        rule_engine = RuleEngine(rules)
        return CorePredictor(rule_engine=rule_engine)

    with gr.TabItem("Pipeline"):
        gr.Markdown(
            """
### Semantic discovery (GBK -> ORFs -> ESM-2 -> ranked candidates)

This pipeline:
1) parses your GBK(s),
2) harvests **all ORFs** via 6-frame translation,
3) embeds ORFs with the **ESM-2 model**,
4) scores similarity vs your chosen **validated precursor set**.
"""
        )

        # Section 1: Input Configuration
        gr.Markdown("---")
        gr.Markdown("#### 1. Input Files")

        pipeline_gbk_files = gr.File(
            label="Upload GBK files (MiBIG / antiSMASH / GenBank)",
            file_count="multiple",
            file_types=[".gbk", ".gb", ".genbank"],
        )

        pipeline_parse_btn = gr.Button("Preview Files")
        pipeline_parse_status = gr.Markdown()
        pipeline_segments_state = gr.State(value=[])

        pipeline_segments_preview = gr.DataFrame(
            label="BGC Segments Preview",
            interactive=False,
        )

        # Section 2: Reference Set Selection
        gr.Markdown("---")
        gr.Markdown("#### 2. Precursor Reference Set")

        with gr.Row():
            pipeline_precursor_set = gr.Dropdown(
                label="Reference Precursors",
                choices=get_available_precursor_sets(),
                value=DEFAULT_PRECURSOR_SET_NAME,
                scale=2,
                info="Pick a .faa file from data/precursors (default is precursor_proteins_verified.faa)"
            )
            pipeline_custom_fasta = gr.File(
                label="Custom FASTA (optional)",
                file_types=[".faa", ".fasta", ".fa"],
                scale=2,
            )

        # Section 3: Parameters
        gr.Markdown("---")
        gr.Markdown("#### 3. Analysis Parameters")

        with gr.Row():
            pipeline_min_aa = gr.Number(
                label="Min ORF Length (aa)",
                value=20,
                minimum=10,
                maximum=100,
                step=5,
            )
            pipeline_max_aa = gr.Number(
                label="Max ORF Length (aa)",
                value=120,
                minimum=30,
                maximum=200,
                step=10,
            )
            pipeline_device = gr.Dropdown(
                label="ESM-2 Device",
                choices=["auto", "cpu", "cuda", "mps"],
                value="auto",
                info="auto = detect best available"
            )
            pipeline_model = gr.Dropdown(
                label="ESM-2 Model",
                choices=[
                    "facebook/esm2_t6_8M_UR50D",
                    "facebook/esm2_t12_35M_UR50D",
                    "facebook/esm2_t30_150M_UR50D",
                    "facebook/esm2_t33_650M_UR50D",
                ],
                value="facebook/esm2_t6_8M_UR50D",
                info="Smaller models use much less RAM (recommended on Mac)"
            )

        # Section 3b: Ranking Options
        gr.Markdown(
            """
**Embedding ordering + rule cutoff**:
- Candidates are always sorted by embedding score.
- Optionally drop candidates whose **rule score** is below a cutoff.
- Set cutoff to **0** to disable rule filtering.
"""
        )
        with gr.Row():
            ranking_score_mode = gr.Dropdown(
                label="Embedding Score Mode",
                choices=["top_n_mean", "best_similarity"],
                value="top_n_mean",
                info="Embedding score used for sorting"
            )
            ranking_rule_cutoff = gr.Slider(
                label="Rule Score Cutoff",
                minimum=0.0,
                maximum=10.0,
                value=0.0,
                step=0.1,
                info="Drop candidates with rule score below this threshold (0 disables)"
            )
            ranking_allow_inverted = gr.Checkbox(
                label="Allow core+leader (experimental)",
                value=False,
                info="Enable inverted core+leader scoring (unreliable; needs dedicated training data)",
            )
        filter_lasso_only = gr.Checkbox(
            label="Filter to lasso BGCs only",
            value=False,
            info="Only show candidates from GBKs with lasso-annotated regions",
        )

        # Section 4: Run
        gr.Markdown("---")
        pipeline_run_btn = gr.Button("Run Pipeline", variant="primary", size="lg")
        pipeline_status = gr.Markdown("**Live progress:** idle")
        pipeline_metadata = gr.Markdown()

        # State for pipeline results (JSON format for visualization)
        pipeline_results_json_state = gr.State(value=[])
        # State for full results dataframe (for export)
        pipeline_full_df_state = gr.State(value=None)
        # State for selected candidates (set of candidate_ids)
        selected_candidates_state = gr.State(value=[])
        # State for core predictions {candidate_id: prediction_result}
        core_predictions_state = gr.State(value={})
        # State for selected cores to export {candidate_id: [core_indices]}
        selected_cores_state = gr.State(value={})

        # Section 5: Results with Sub-Tabs
        gr.Markdown("---")
        gr.Markdown("#### 4. Results")

        # Selection controls
        with gr.Row():
            selection_display = gr.Markdown("**Selected:** 0 candidates")
            select_all_btn = gr.Button("Select All Visible", size="sm")
            clear_selection_btn = gr.Button("Clear Selection", size="sm")

        # Sub-tabs for Table and Visualization views
        with gr.Tabs() as results_tabs:
            with gr.TabItem("Table View"):
                pipeline_results_df = gr.DataFrame(
                    label="Scored Candidates (click rows to select/deselect)",
                    interactive=False,
                    wrap=True,
                )

            with gr.TabItem("Visualization"):
                pipeline_viz_html = gr.HTML(
                    value="<div style='padding:40px;text-align:center;color:#64748b;'>Run the pipeline to see visualization</div>"
                )

        # Hidden components for selection sync (iframe <-> Gradio)
        viz_selection_receiver = gr.Textbox(
            visible=True,
            elem_id="viz-selection-receiver",
            show_label=False,
        )
        viz_selection_sender = gr.Textbox(
            visible=True,
            elem_id="viz-selection-sender",
            show_label=False,
        )
        gr.HTML(
            """
            <style>
              #viz-selection-receiver, #viz-selection-sender { display: none !important; }
            </style>
            <script>
            (function() {
              function getReceiverInput() {
                const wrapper = document.getElementById("viz-selection-receiver");
                if (!wrapper) return null;
                if (wrapper.tagName === "TEXTAREA" || wrapper.tagName === "INPUT") {
                  return wrapper;
                }
                return wrapper.querySelector("textarea, input");
              }
              function getSenderInput() {
                const wrapper = document.getElementById("viz-selection-sender");
                if (!wrapper) return null;
                if (wrapper.tagName === "TEXTAREA" || wrapper.tagName === "INPUT") {
                  return wrapper;
                }
                return wrapper.querySelector("textarea, input");
              }
              function handleMessage(event) {
                if (!event || !event.data || event.data.type !== "candidate_selection") return;
                const receiver = getReceiverInput();
                if (!receiver) return;
                const payload = JSON.stringify(event.data.candidates || []);
                receiver.value = payload;
                receiver.dispatchEvent(new Event("input", { bubbles: true }));
                receiver.dispatchEvent(new Event("change", { bubbles: true }));
              }
              function sendSelection() {
                const sender = getSenderInput();
                const iframe = document.getElementById("viz-iframe");
                if (!sender || !iframe || !iframe.contentWindow) return;
                let selection = [];
                try {
                  selection = JSON.parse(sender.value || "[]");
                } catch (e) {
                  selection = [];
                }
                iframe.contentWindow.postMessage(
                  { type: "set_selection", candidates: selection },
                  "*"
                );
              }
              function attachSender() {
                const sender = getSenderInput();
                if (!sender) return;
                sender.addEventListener("input", sendSelection);
                sender.addEventListener("change", sendSelection);
              }
              window.addEventListener("message", handleMessage);
              if (document.readyState === "loading") {
                document.addEventListener("DOMContentLoaded", attachSender);
              } else {
                attachSender();
              }
            })();
            </script>
            """
        )

        # Section 6: Core Prediction
        gr.Markdown("---")
        gr.Markdown("#### 5. Core Prediction (Selected Candidates)")

        with gr.Row():
            core_top_n = gr.Slider(
                label="Top-N Cores per Candidate",
                minimum=1,
                maximum=10,
                value=3,
                step=1,
            )
            core_predict_btn = gr.Button("Predict Cores for Selected", variant="primary")

        core_prediction_status = gr.Markdown()

        # Accordion-based display for predicted cores
        core_predictions_display = gr.HTML(
            value="<div style='padding:20px;text-align:center;color:#64748b;'>Select candidates and click 'Predict Cores' to see predictions</div>"
        )

        # Section 7: Export
        gr.Markdown("---")
        gr.Markdown("#### 6. Export")

        with gr.Row():
            export_json_btn = gr.Button("Export Selected (Candidates + Cores) JSON")
            export_fasta_btn = gr.Button("Export Selected (Candidates + Cores) FASTA")
            download_full_csv_btn = gr.Button("Download Full CSV")
            download_full_json_btn = gr.Button("Download Full JSON")

        export_output = gr.File(label="Export Download", interactive=False)
        pipeline_download_csv = gr.File(label="Full CSV", interactive=False)
        pipeline_download_json = gr.File(label="Full JSON", interactive=False)

        # JSON report viewer (collapsible)
        with gr.Accordion("JSON Report", open=False):
            pipeline_json_output = gr.Code(
                label="Pipeline Report JSON",
                language="json",
                interactive=False,
            )

        # ========== Event Handlers ==========

        def parse_files(files):
            if not files:
                return "No files", [], pd.DataFrame()
            file_paths = [f.name for f in files]
            segments, df = pipeline_service.parse_uploaded_gbk_files(file_paths)
            return f"Parsed {len(segments)} segments", segments, df

        pipeline_parse_btn.click(
            parse_files,
            inputs=[pipeline_gbk_files],
            outputs=[pipeline_parse_status, pipeline_segments_state, pipeline_segments_preview]
        )

        def _build_display_df(full_df: pd.DataFrame) -> pd.DataFrame:
            """Build a UI-friendly table while preserving lossless machine columns."""
            if full_df is None or full_df.empty:
                return pd.DataFrame()

            display_df = full_df.copy()

            if "record_id" in display_df.columns and "Record" not in display_df.columns:
                display_df["Record"] = display_df["record_id"]
            if "candidate_id" in display_df.columns and "Candidate" not in display_df.columns:
                display_df["Candidate"] = display_df["candidate_id"]
            if "rule_score_raw" in display_df.columns and "Rule Score" not in display_df.columns:
                display_df["Rule Score"] = display_df["rule_score_raw"]

            if "Score" not in display_df.columns:
                score_series = pd.Series(index=display_df.index, dtype=float)
                for score_col in ("embedding_score", "top_n_mean_similarity", "best_similarity"):
                    if score_col in display_df.columns:
                        score_series = score_series.fillna(display_df[score_col])
                display_df["Score"] = score_series

            leading = ["Record", "Candidate", "Score", "Rule Score"]
            ordered = [c for c in leading if c in display_df.columns]
            trailing = [c for c in display_df.columns if c not in ordered]
            return display_df.loc[:, ordered + trailing]

        def _candidate_id_from_row(row: pd.Series) -> str:
            for key in ("Candidate", "candidate_id"):
                value = row.get(key)
                if pd.isna(value):
                    continue
                candidate_id = str(value).strip()
                if candidate_id:
                    return candidate_id
            return ""

        def run_analysis(
            segments,
            precursor_set,
            custom_fasta,
            min_aa,
            max_aa,
            device,
            model,
            ranking_score_mode,
            ranking_rule_cutoff,
            ranking_allow_inverted,
            filter_lasso_only,
            progress=gr.Progress(),
        ):
            if not segments:
                empty_viz = "<div style='padding:40px;text-align:center;color:#64748b;'>Please parse files first</div>"
                return "Please parse files first", "", pd.DataFrame(), "", [], pd.DataFrame(), empty_viz

            # Resolve precursor path
            validated_faa = DEFAULT_PRECURSORS_PATH
            if custom_fasta:
                validated_faa = Path(custom_fasta.name)
            elif precursor_set and precursor_set != UPLOAD_CUSTOM_LABEL:
                candidate = DATA_DIR / "precursors" / precursor_set
                if candidate.exists():
                    validated_faa = candidate

            ranking_config = RankingConfig(
                score_mode=str(ranking_score_mode),
                min_rule_score=float(ranking_rule_cutoff),
                allow_inverted_rules=bool(ranking_allow_inverted),
            )
            ranker_predictor = None
            if ranking_config.min_rule_score > 0:
                ranker_predictor = get_configured_predictor()

            def p_callback(msg):
                progress(0.5, desc=msg)

            try:
                results, _ = pipeline_service.run_pipeline(
                    segments=segments,
                    validated_faa_path=validated_faa,
                    min_aa=int(min_aa),
                    max_aa=int(max_aa),
                    device=device,
                    model_name=model,
                    progress_callback=p_callback,
                    ranking_config=ranking_config,
                    ranker_predictor=ranker_predictor,
                    filter_lasso_only=bool(filter_lasso_only),
                )

                # Create JSON report + full dataframe
                from lasso_workbench.pipeline.semantic_pipeline import results_to_json, results_to_dataframe
                json_data = results_to_json(results)
                full_df = results_to_dataframe(results)
                display_df = _build_display_df(full_df)
                json_str = json.dumps(json_data, indent=2)

                metadata_text = ""
                if results:
                    meta = results[0].pipeline_metadata
                    metadata_text = (
                        f"**Pipeline metadata**  \n"
                        f"- Model: `{meta.model_name}`  \n"
                        f"- Device: `{meta.device}`  \n"
                        f"- Validated FASTA: `{meta.validated_faa}`  \n"
                        f"- ORF length: `{meta.min_aa}-{meta.max_aa}`"
                    )

                # Generate visualization HTML
                viz_html = generate_gradio_iframe_html(json_data, initial_selection=[], title="Pipeline Results")

                status = f"Completed! Found {len(display_df)} candidates across {len(results)} BGCs."
                return status, metadata_text, display_df, json_str, json_data, full_df, viz_html

            except Exception as e:
                import traceback
                traceback.print_exc()
                empty_viz = f"<div style='padding:40px;text-align:center;color:#dc2626;'>Error: {e}</div>"
                return f"Error: {e}", "", pd.DataFrame(), "", [], pd.DataFrame(), empty_viz

        pipeline_run_btn.click(
            run_analysis,
            inputs=[
                pipeline_segments_state,
                pipeline_precursor_set,
                pipeline_custom_fasta,
                pipeline_min_aa,
                pipeline_max_aa,
                pipeline_device,
                pipeline_model,
                ranking_score_mode,
                ranking_rule_cutoff,
                ranking_allow_inverted,
                filter_lasso_only,
            ],
            outputs=[
                pipeline_status,
                pipeline_metadata,
                pipeline_results_df,
                pipeline_json_output,
                pipeline_results_json_state,
                pipeline_full_df_state,
                pipeline_viz_html,
            ]
        )

        def export_full_csv(df: pd.DataFrame):
            if df is None or df.empty:
                return None
            with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
                df.to_csv(tmp.name, index=False)
                return tmp.name

        def export_full_json(json_data: List[Dict]):
            if not json_data:
                return None
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as tmp:
                json.dump(json_data, tmp, indent=2)
                return tmp.name

        download_full_csv_btn.click(
            export_full_csv,
            inputs=[pipeline_full_df_state],
            outputs=[pipeline_download_csv],
        )

        download_full_json_btn.click(
            export_full_json,
            inputs=[pipeline_results_json_state],
            outputs=[pipeline_download_json],
        )

        def on_table_row_select(evt: gr.SelectData, df: pd.DataFrame, current_selection: List[str]):
            """Toggle selection when a row is clicked in the table."""
            if evt is None or df is None or df.empty:
                return current_selection, f"**Selected:** {len(current_selection)} candidates"

            row_idx = evt.index[0] if isinstance(evt.index, (list, tuple)) else int(evt.index)
            if row_idx >= len(df):
                return current_selection, f"**Selected:** {len(current_selection)} candidates"

            row = df.iloc[row_idx]
            candidate_id = _candidate_id_from_row(row)
            if not candidate_id:
                return current_selection, f"**Selected:** {len(current_selection)} candidates"

            # Toggle selection
            selection_set = set(current_selection)
            if candidate_id in selection_set:
                selection_set.remove(candidate_id)
            else:
                selection_set.add(candidate_id)

            new_selection = list(selection_set)
            return new_selection, f"**Selected:** {len(new_selection)} candidates"

        pipeline_results_df.select(
            on_table_row_select,
            inputs=[pipeline_results_df, selected_candidates_state],
            outputs=[selected_candidates_state, selection_display],
        )

        def sync_selection_from_viz(selection_json: str, current_selection: List[str]):
            if not selection_json:
                return current_selection, f"**Selected:** {len(current_selection)} candidates"
            try:
                selection = json.loads(selection_json)
            except json.JSONDecodeError:
                return current_selection, f"**Selected:** {len(current_selection)} candidates"
            if not isinstance(selection, list):
                return current_selection, f"**Selected:** {len(current_selection)} candidates"
            if set(selection) == set(current_selection):
                return current_selection, f"**Selected:** {len(current_selection)} candidates"
            return selection, f"**Selected:** {len(selection)} candidates"

        viz_selection_receiver.change(
            sync_selection_from_viz,
            inputs=[viz_selection_receiver, selected_candidates_state],
            outputs=[selected_candidates_state, selection_display],
        )

        def select_all_visible(df: pd.DataFrame, json_data: List[Dict]):
            """Select all visible candidates from the results."""
            if df is None or df.empty:
                return [], "**Selected:** 0 candidates"

            if "Candidate" in df.columns:
                all_ids = [str(v) for v in df["Candidate"].tolist() if str(v).strip()]
            elif "candidate_id" in df.columns:
                all_ids = [str(v) for v in df["candidate_id"].tolist() if str(v).strip()]
            else:
                all_ids = []
            return all_ids, f"**Selected:** {len(all_ids)} candidates"

        select_all_btn.click(
            select_all_visible,
            inputs=[pipeline_results_df, pipeline_results_json_state],
            outputs=[selected_candidates_state, selection_display],
        )

        def clear_selection():
            """Clear all selections."""
            return [], "**Selected:** 0 candidates"

        clear_selection_btn.click(
            clear_selection,
            inputs=[],
            outputs=[selected_candidates_state, selection_display],
        )

        def selection_to_sender(selection: List[str]):
            """Send selection updates to the visualization iframe."""
            return json.dumps(selection or [])

        selected_candidates_state.change(
            selection_to_sender,
            inputs=[selected_candidates_state],
            outputs=[viz_selection_sender],
        )

        def predict_cores_batch(
            selected_ids: List[str],
            json_data: List[Dict],
            top_n: int,
            progress=gr.Progress(),
        ):
            """Predict cores for all selected candidates using existing CorePredictor."""
            if not selected_ids:
                empty_html = "<div style='padding:20px;text-align:center;color:#64748b;'>No candidates selected</div>"
                return "No candidates selected.", {}, empty_html

            # Build lookup from candidate_id to sequence
            candidate_lookup = {}
            for bgc in json_data:
                for cand in bgc.get("candidates", []):
                    candidate_lookup[cand["candidate_id"]] = {
                        "sequence": cand.get("protein_sequence", ""),
                        "record_id": bgc.get("record_id", ""),
                        "score": cand.get("embedding_score") or cand.get("top_n_mean_similarity") or cand.get("best_similarity", 0),
                        "aa_length": cand.get("aa_length", 0),
                    }

            # Get predictor (reuses existing logic)
            predictor = get_configured_predictor()

            predictions = {}
            processed = 0
            total = len(selected_ids)

            for candidate_id in selected_ids:
                if candidate_id not in candidate_lookup:
                    continue

                cand_data = candidate_lookup[candidate_id]
                seq = cand_data["sequence"].upper()
                seq = "".join(c for c in seq if c.isalpha())

                if not seq:
                    continue

                progress((processed / total), desc=f"Predicting cores for {candidate_id[:20]}...")

                # Use existing CorePredictor.predict()
                result = predictor.predict(
                    seq,
                    top_n=int(top_n),
                    allow_inverted=True,
                )

                predictions[candidate_id] = {
                    "record_id": cand_data["record_id"],
                    "sequence": seq,
                    "score": cand_data["score"],
                    "aa_length": cand_data["aa_length"],
                    "cores": [
                        {
                            "rank": i + 1,
                            "core": pred.core,
                            "leader": pred.leader,
                            "cleavage_site": pred.cleavage_site,
                            "score": pred.score,
                            "confidence": pred.confidence,
                            "reasons": pred.reasons,
                            "orientation": pred.orientation,
                        }
                        for i, pred in enumerate(result.predictions)
                    ]
                }
                processed += 1

            # Generate accordion HTML for display
            html = generate_cores_accordion_html(predictions)
            status = f"Predicted cores for {len(predictions)} candidates."
            return status, predictions, html

        core_predict_btn.click(
            predict_cores_batch,
            inputs=[
                selected_candidates_state,
                pipeline_results_json_state,
                core_top_n,
            ],
            outputs=[
                core_prediction_status,
                core_predictions_state,
                core_predictions_display,
            ],
        )

        def export_json(
            selected_ids: List[str],
            json_data: List[Dict],
            core_predictions: Dict,
        ):
            """Export selected candidates with predicted cores as JSON."""
            if not selected_ids:
                return None

            # Build export structure
            export_data = {
                "candidates": [],
                "export_timestamp": datetime.now(timezone.utc).isoformat(),
                "workbench_version": __version__,
            }

            # Build lookup
            candidate_lookup = {}
            for bgc in json_data:
                for cand in bgc.get("candidates", []):
                    candidate_lookup[cand["candidate_id"]] = {
                        **cand,
                        "record_id": bgc.get("record_id", ""),
                    }

            for candidate_id in selected_ids:
                if candidate_id not in candidate_lookup:
                    continue

                cand = candidate_lookup[candidate_id]
                entry = {
                    "candidate_id": candidate_id,
                    "record_id": cand.get("record_id", ""),
                    "protein_sequence": cand.get("protein_sequence", ""),
                    "similarity_score": cand.get("embedding_score") or cand.get("top_n_mean_similarity") or cand.get("best_similarity", 0),
                    "aa_length": cand.get("aa_length", 0),
                    "genomic_start": cand.get("genomic_start", 0),
                    "genomic_end": cand.get("genomic_end", 0),
                    "strand": cand.get("strand", 1),
                    "rule_orientation": cand.get("rule_orientation"),
                    "predicted_cores": [],
                }

                # Add core predictions if available
                if candidate_id in core_predictions:
                    entry["predicted_cores"] = core_predictions[candidate_id].get("cores", [])

                export_data["candidates"].append(entry)

            # Save to temp file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as tmp:
                json.dump(export_data, tmp, indent=2)
                return tmp.name

        export_json_btn.click(
            export_json,
            inputs=[selected_candidates_state, pipeline_results_json_state, core_predictions_state],
            outputs=[export_output],
        )

        def export_fasta(
            selected_ids: List[str],
            json_data: List[Dict],
            core_predictions: Dict,
        ):
            """Export selected candidates with predicted cores as FASTA."""
            if not selected_ids:
                return None

            # Build lookup
            candidate_lookup = {}
            for bgc in json_data:
                for cand in bgc.get("candidates", []):
                    candidate_lookup[cand["candidate_id"]] = {
                        **cand,
                        "record_id": bgc.get("record_id", ""),
                    }

            entries = []
            for candidate_id in selected_ids:
                if candidate_id not in candidate_lookup:
                    continue

                cand = candidate_lookup[candidate_id]
                score = cand.get("embedding_score") or cand.get("top_n_mean_similarity") or cand.get("best_similarity", 0)

                # Export the full precursor sequence
                seq = cand.get("protein_sequence", "")
                desc = f"record={cand.get('record_id', '')} score={score:.4f} length={len(seq)}aa"
                entries.append((candidate_id, seq, desc))

                # Export predicted cores if available
                if candidate_id in core_predictions:
                    for core_data in core_predictions[candidate_id].get("cores", []):
                        core_desc = (
                            f"score={core_data['score']:.2f} "
                            f"confidence={core_data['confidence']} "
                            f"cleavage={core_data['cleavage_site']}"
                        )
                        core_id = f"{candidate_id}|core_{core_data['rank']}"
                        entries.append((core_id, core_data["core"], core_desc))

            # Save to temp file
            tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
            tmp_path = Path(tmp.name)
            tmp.close()
            write_fasta_entries(entries, tmp_path)
            return str(tmp_path)

        export_fasta_btn.click(
            export_fasta,
            inputs=[selected_candidates_state, pipeline_results_json_state, core_predictions_state],
            outputs=[export_output],
        )

    return pipeline_gbk_files


def generate_cores_accordion_html(predictions: Dict[str, Dict]) -> str:
    """Generate accordion HTML for displaying core predictions."""
    if not predictions:
        return "<div style='padding:20px;text-align:center;color:#64748b;'>No predictions available</div>"

    html_parts = ['<div style="font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;">']

    for candidate_id, data in predictions.items():
        record_id = data.get("record_id", "")
        score = data.get("score", 0)
        aa_length = data.get("aa_length", 0)
        cores = data.get("cores", [])

        # Accordion header
        html_parts.append(f'''
        <details style="margin-bottom:8px;border:1px solid #e5e7eb;border-radius:6px;overflow:hidden;">
            <summary style="padding:12px 16px;background:#f8fafc;cursor:pointer;font-weight:600;font-size:0.9rem;color:#374151;display:flex;justify-content:space-between;align-items:center;">
                <span>{_escape(candidate_id[:40])}</span>
                <span style="font-weight:normal;font-size:0.8rem;color:#64748b;">{aa_length} aa | score {score:.4f} | {len(cores)} cores</span>
            </summary>
            <div style="padding:12px 16px;background:white;">
                <div style="font-size:0.75rem;color:#64748b;margin-bottom:8px;">Record: {_escape(record_id)}</div>
        ''')

        if not cores:
            html_parts.append('<div style="color:#94a3b8;font-style:italic;">No core predictions</div>')
        else:
            for core_data in cores:
                confidence = core_data.get("confidence", "unknown")
                conf_color = "#059669" if confidence == "high" else "#d97706" if confidence == "medium" else "#dc2626"
                orientation = core_data.get("orientation", "leader_core")
                orientation_label = "leader+core" if orientation == "leader_core" else "core+leader"
                orientation_color = "#2563eb" if orientation == "leader_core" else "#f97316"
                leader_seq = _escape(core_data.get("leader", ""))
                core_seq = _escape(core_data.get("core", ""))
                if orientation == "core_leader":
                    sequence_block = (
                        "<div style='margin-bottom:4px;'>"
                        "<span style='display:inline-block;background:#2563eb;color:white;border-radius:3px;padding:1px 6px;font-size:10px;margin-right:6px;'>Core</span>"
                        f"<span style='color:#0f172a;font-weight:600;'>{core_seq}</span>"
                        "</div>"
                        "<div>"
                        "<span style='display:inline-block;background:#e2e8f0;color:#0f172a;border-radius:3px;padding:1px 6px;font-size:10px;margin-right:6px;'>Leader</span>"
                        f"<span style='color:#334155;'>{leader_seq}</span>"
                        "</div>"
                    )
                else:
                    sequence_block = (
                        "<div style='margin-bottom:4px;'>"
                        "<span style='display:inline-block;background:#e2e8f0;color:#0f172a;border-radius:3px;padding:1px 6px;font-size:10px;margin-right:6px;'>Leader</span>"
                        f"<span style='color:#334155;'>{leader_seq}</span>"
                        "</div>"
                        "<div>"
                        "<span style='display:inline-block;background:#2563eb;color:white;border-radius:3px;padding:1px 6px;font-size:10px;margin-right:6px;'>Core</span>"
                        f"<span style='color:#0f172a;font-weight:600;'>{core_seq}</span>"
                        "</div>"
                    )

                html_parts.append(f'''
                <div style="padding:10px;margin-bottom:6px;background:#f8fafc;border-radius:4px;border-left:3px solid {conf_color};">
                    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;">
                        <span style="font-weight:600;font-size:0.85rem;">Core {core_data.get("rank", "?")}</span>
                        <span style="display:flex;gap:6px;align-items:center;">
                            <span style="font-size:0.72rem;padding:2px 6px;border-radius:3px;background:{orientation_color};color:white;">{orientation_label}</span>
                            <span style="font-size:0.75rem;padding:2px 6px;border-radius:3px;background:{conf_color};color:white;">{confidence}</span>
                        </span>
                    </div>
                    <div style="font-family:SF Mono,Monaco,Consolas,monospace;font-size:11px;background:white;padding:6px 8px;border-radius:3px;word-break:break-all;border:1px solid #e5e7eb;">
                        {sequence_block}
                    </div>
                    <div style="margin-top:6px;font-size:0.75rem;color:#64748b;">
                        Score: {core_data.get("score", 0):.2f} | Cleavage: {core_data.get("cleavage_site", "?")} | Leader: {len(core_data.get("leader", ""))} aa
                    </div>
                    {f"<div style='margin-top:4px;font-size:0.72rem;color:#475569;'>Contributors: {_escape(', '.join(core_data.get('reasons', [])))}</div>" if core_data.get("reasons") else ""}
                </div>
                ''')

        html_parts.append('</div></details>')

    html_parts.append('</div>')
    return "".join(html_parts)


def _escape(s: str) -> str:
    """HTML-escape a string."""
    import html
    return html.escape(str(s) if s else "")
