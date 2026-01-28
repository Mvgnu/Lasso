import gradio as gr
import pandas as pd
from pathlib import Path
from lasso_workbench import DATA_DIR
from lasso_workbench.services.data_service import DataService

def get_available_datasets() -> list:
    """Get list of available dataset files."""
    datasets = []
    
    # Training datasets
    training_dir = DATA_DIR / "training"
    if training_dir.exists():
        for f in training_dir.glob("*.json"):
            datasets.append(f"training/{f.name}")
    
    # Lab datasets
    lab_dir = DATA_DIR / "lab_dataset"
    if lab_dir.exists():
        for f in lab_dir.glob("*.csv"):
            datasets.append(f"lab/{f.name}")
        for f in lab_dir.glob("*.json"):
            datasets.append(f"lab/{f.name}")
    
    # Precursor datasets (metadata, cases, candidates)
    precursor_dir = DATA_DIR / "precursors"
    if precursor_dir.exists():
        protected = {
            "lab_core_candidates_multi_strict.tsv",
            "lab_core_loci_multi_strict.tsv",
            "lab_core_cases_summary.tsv",
        }
        for f in precursor_dir.glob("*.tsv"):
            if f.name in protected:
                continue
            datasets.append(f"precursors/{f.name}")
    
    return datasets

def create_dataset_tab(data_service: DataService):
    """Create the Dataset Manager tab."""

    def resolve_dataset_path(name: str) -> Path:
        if not name:
            return DATA_DIR / "training" / "train.json"
        parts = name.split("/")
        if len(parts) == 2:
            cat, fname = parts
            if cat == "training": return DATA_DIR / "training" / fname
            if cat == "lab": return DATA_DIR / "lab_dataset" / fname
            if cat == "precursors": return DATA_DIR / "precursors" / fname
        return DATA_DIR / "training" / name

    def load_handler(dataset_name):
        if not dataset_name:
            return "No dataset selected.", pd.DataFrame()
        path = resolve_dataset_path(dataset_name)
        if not path.exists():
            return f"Dataset not found: {path}", pd.DataFrame()
        
        data_service.set_dataset_path(path)
        df = data_service.load_dataset()
        suffix = " (read-only)" if data_service.is_read_only() else ""
        return f"Loaded **{dataset_name}** ({len(df)} entries){suffix}", df

    def filter_handler(query, source_filter):
        # Filtering happens in memory on loaded df
        df = data_service.load_dataset()
        if df.empty: return df
        
        if source_filter != "all" and "source" in df.columns:
            df = df[df["source"] == source_filter]
            
        if query:
            q = query.lower()
            mask = df.apply(lambda row: any(q in str(val).lower() for val in row), axis=1)
            df = df[mask]
            
        return df

    def refresh_handler():
        return data_service.load_dataset()

    def add_handler(name, precursor, leader, core, source):
        if not name or not precursor:
            return "Name and Precursor are required.", data_service.load_dataset()
        if data_service.is_read_only():
            return "Dataset is read-only (use scripts/add_validated_precursors.py).", data_service.load_dataset()
        
        entry = {
            "name": name, 
            "full_precursor": precursor,  # Maintain compatibility with logic expecting full_precursor
            "sequence": precursor,        # Also save as sequence
            "leader": leader, 
            "core": core, 
            "source": source
        }
        try:
            data_service.add_entry(entry)
            return f"Added entry: {name}", data_service.load_dataset()
        except Exception as e:
            return f"Error adding entry: {e}", data_service.load_dataset()

    def update_handler(idx_val, name, precursor, leader, core, source):
        if idx_val is None:
            return "No entry selected.", data_service.load_dataset()
        if data_service.is_read_only():
            return "Dataset is read-only (use scripts/add_validated_precursors.py).", data_service.load_dataset()
        
        try:
            # idx_val comes from State, should be int
            idx = int(idx_val)
            updates = {
                "name": name,
                "full_precursor": precursor,
                "sequence": precursor,
                "leader": leader,
                "core": core,
                "source": source
            }
            success = data_service.update_entry(idx, updates)
            if success:
                return f"Updated entry #{idx}", data_service.load_dataset()
            else:
                return f"Failed to update entry #{idx}", data_service.load_dataset()
        except Exception as e:
            return f"Error updating: {e}", data_service.load_dataset()

    def delete_handler(idx_val):
        if idx_val is None:
             return "No entry selected.", data_service.load_dataset()
        if data_service.is_read_only():
            return "Dataset is read-only (use scripts/add_validated_precursors.py).", data_service.load_dataset()
        try:
            idx = int(idx_val)
            success = data_service.delete_entry(idx)
            if success:
                return f"Deleted entry #{idx}", data_service.load_dataset()
            else:
                return f"Failed to delete entry #{idx}", data_service.load_dataset()
        except Exception as e:
            return f"Error deleting: {e}", data_service.load_dataset()

    def import_handler(file_obj):
        if file_obj is None:
            return "No file selected.", data_service.load_dataset()
        if data_service.is_read_only():
            return "Dataset is read-only (use scripts/add_validated_precursors.py).", data_service.load_dataset()
        try:
            count = data_service.import_file(Path(file_obj.name))
            return f"Imported {count} entries.", data_service.load_dataset()
        except Exception as e:
            return f"Import error: {e}", data_service.load_dataset()

    def export_handler(fmt):
        # DataService saves to current path. Exporting to a download file requires copying or saving temp.
        # Minimal implementation: return current file path if matches format, or convert.
        # Actually simplest is to write current DF to a temp file requested format
        import tempfile
        df = data_service.load_dataset()
        suffix = ".json" if fmt == "JSON" else ".csv" if fmt == "CSV" else ".tsv"
        sep = "\t" if fmt == "TSV" else ","
        
        with tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as tmp:
            if fmt == "JSON":
                df.to_json(tmp.name, orient="records", indent=2)
            else:
                df.to_csv(tmp.name, sep=sep, index=False)
            return tmp.name, f"Exported to {fmt}"

    with gr.TabItem("Dataset"):
        gr.Markdown("### Dataset")
        
        # Dataset selector
        gr.Markdown("---")
        dataset_choices = get_available_datasets()
        default_dataset = None
        verified_name = "precursors/precursor_proteins_verified.tsv"
        if verified_name in dataset_choices:
            default_dataset = verified_name
        elif dataset_choices:
            default_dataset = dataset_choices[0]

        with gr.Row():
            dataset_selector = gr.Dropdown(
                label="Select Dataset",
                choices=dataset_choices,
                value=default_dataset,
                scale=2,
                info="Import a file to create a dataset" if not dataset_choices else None,
            )
            load_dataset_btn = gr.Button("Load Dataset", scale=1)
        
        dataset_load_status = gr.Markdown()
        
        gr.Markdown("---")
        gr.Markdown("### Current Dataset")

        # Search and filter
        with gr.Row():
            search_input = gr.Textbox(label="Search", placeholder="Search by name, core, or sequence...", scale=2)
            source_filter = gr.Dropdown(
                label="Filter by source",
                choices=["all", "lab", "mibig", "prediction", "literature", "import", "genome"],
                value="all",
                scale=1,
            )
            filter_btn = gr.Button("Filter", scale=1)
        
        refresh_btn = gr.Button("🔄 Refresh Dataset")
        dataset_view = gr.DataFrame(
            label="Dataset (click row to select)",
            interactive=False,
            wrap=False,
            elem_id="dataset-table",
        )
        
        # Hidden state for selected index
        selected_idx = gr.State(value=None)
        
        # Selection display
        selected_display = gr.Markdown("**Selected:** None - click a row to select")
        
        # Dataset actions with selected row
        gr.Markdown("---")
        
        with gr.Row():
            with gr.Column(scale=2):
                gr.Markdown("### Edit Selected Entry")
                entry_name = gr.Textbox(label="Name", placeholder="e.g., Microcin J25")
                entry_source = gr.Dropdown(
                    label="Source",
                    choices=["lab", "mibig", "prediction", "literature", "import", "genome"],
                    value="lab",
                )
                entry_precursor = gr.Textbox(label="Full Precursor", placeholder="Full amino acid sequence", lines=3)
                with gr.Row():
                    entry_leader = gr.Textbox(label="Leader (optional)")
                    entry_core = gr.Textbox(label="Core (optional)")
                
                with gr.Row():
                    add_btn = gr.Button("Add New Entry", variant="secondary")
                    update_btn = gr.Button("Update Selected", variant="primary")
                    delete_btn = gr.Button("Delete Selected", variant="stop")
                
                action_status = gr.Markdown()
            
            with gr.Column(scale=1):
                gr.Markdown("### Import / Export")
                
                gr.Markdown("**Import**")
                import_file = gr.File(label="Import file (CSV, TSV, JSON, FASTA)", file_types=[".csv", ".tsv", ".json", ".fasta", ".fa", ".faa"])
                import_btn = gr.Button("Import File")
                
                gr.Markdown("**Export**")
                export_format = gr.Radio(["JSON", "CSV", "TSV"], label="Format", value="JSON")
                export_btn = gr.Button("Export Dataset")
                export_file = gr.File(label="Download Export")
                export_status = gr.Markdown()
        
        # Listeners
        load_dataset_btn.click(load_handler, inputs=[dataset_selector], outputs=[dataset_load_status, dataset_view])
        dataset_selector.change(load_handler, inputs=[dataset_selector], outputs=[dataset_load_status, dataset_view])
        
        filter_btn.click(filter_handler, inputs=[search_input, source_filter], outputs=[dataset_view])
        refresh_btn.click(refresh_handler, outputs=[dataset_view])
        
        def on_select(evt: gr.SelectData, df):
            if evt is None or df is None or df.empty:
                return None, "**Selected:** None", "", "lab", "", "", ""
            row_idx = evt.index[0]
            if row_idx >= len(df): return None, "**Selected:** None", "", "lab", "", "", ""
            
            # Use 'idx' column if available (stable id), else row_idx
            try:
                actual_idx = int(df.iloc[row_idx]["idx"])
            except:
                actual_idx = row_idx
            
            # Fetch from DF is hard if we filtered. Better fetch from current full DF?
            # Or just use the row data from the view if columns match.
            # But the view might be filtered.
            # For simplicity, if we filter, we might not get correct index unless we preserve 'idx'.
            # DataService 'idx' is preserved.
            
            # We need to look up in full dataset by actual_idx
            full_df = data_service.load_dataset()
            if "idx" in full_df.columns:
                 row = full_df[full_df["idx"] == actual_idx]
                 if row.empty: return None, "**Selected:** None", "", "lab", "", "", ""
                 row = row.iloc[0]
            else:
                 # fallback
                 if actual_idx < len(full_df): row = full_df.iloc[actual_idx]
                 else: return None, "**Selected:** None", "", "lab", "", "", ""
            
            name = row.get("name", "")
            source = row.get("source", "lab")
            precursor = row.get("full_precursor", "") or row.get("sequence", "")
            leader = row.get("leader", "")
            core = row.get("core", "")
            
            return actual_idx, f"**Selected:** #{actual_idx} - {name}", name, source, precursor, leader, core

        dataset_view.select(
            on_select,
            inputs=[dataset_view],
            outputs=[selected_idx, selected_display, entry_name, entry_source, entry_precursor, entry_leader, entry_core]
        )
        
        add_btn.click(add_handler, inputs=[entry_name, entry_precursor, entry_leader, entry_core, entry_source], outputs=[action_status, dataset_view])
        update_btn.click(update_handler, inputs=[selected_idx, entry_name, entry_precursor, entry_leader, entry_core, entry_source], outputs=[action_status, dataset_view])
        delete_btn.click(delete_handler, inputs=[selected_idx], outputs=[action_status, dataset_view])
        
        import_btn.click(import_handler, inputs=[import_file], outputs=[action_status, dataset_view])
        export_btn.click(export_handler, inputs=[export_format], outputs=[export_file, export_status])

    return dataset_view
