import gradio as gr
from pathlib import Path
from lasso_workbench.services.rules_service import RulesService
from lasso_workbench import DATA_DIR

def create_rules_tab(rules_service: RulesService) -> gr.TabItem:
    with gr.TabItem("Rules") as rules_tab:
        gr.Markdown("### Configure prediction rules")
        gr.Markdown("Adjust weights and parameters for core prediction. Click **Save Rules** to persist changes.")
        
        # Preset section
        preset_choices = ["(keep current)"] + rules_service.list_presets()
        with gr.Row():
            preset_select = gr.Dropdown(
                label="Load Preset",
                choices=preset_choices,
                value="(keep current)",
                scale=3,
            )
            load_preset_btn = gr.Button("Apply Preset", scale=1)
            reset_btn = gr.Button("Reset to Defaults", scale=1)
        
        gr.Markdown("---")
        
        with gr.Accordion("Core Structure Rules", open=True):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("#### Core Start Residue")
                    cs_enabled = gr.Checkbox(label="Enabled", value=True)
                    cs_weight = gr.Slider(label="Weight", minimum=0.0, maximum=5.0, step=0.1, value=1.2)
                    cs_scores_df = gr.Dataframe(
                        headers=["AA", "Score"],
                        datatype=["str", "number"],
                        column_count=(2, "fixed"),
                        type="array",
                        label="AA Scores",
                        interactive=True,
                        value=[["G", 1.5], ["C", 1.3], ["A", 1.1], ["S", 1.1]]
                    )
                    cs_penalty_other = gr.Number(label="Penalty (other residues)", value=-0.1)
                    cs_hard_forbidden = gr.Textbox(label="Hard-forbidden residues", value="P")
                    cs_penalty_forbidden = gr.Number(label="Penalty (forbidden)", value=-2.0)
                
                with gr.Column():
                    gr.Markdown("#### Ring Closure")
                    ring_enabled = gr.Checkbox(label="Enabled", value=True)
                    ring_weight = gr.Slider(label="Weight", minimum=0.0, maximum=5.0, step=0.1, value=4.5)
                    ring_allowed_positions = gr.Textbox(label="Allowed positions", value="6,7,8,9,10")
                    ring_residues = gr.Textbox(label="Ring residues", value="D,E")
                    ring_position_weights = gr.Dataframe(
                        headers=["Position", "Weight"],
                        datatype=["number", "number"],
                        column_count=(2, "fixed"),
                        type="array",
                        label="Position Weights",
                        interactive=True,
                        value=[[6, 0.6], [7, 1.0], [8, 1.1], [9, 1.0], [10, 0.6]],
                    )
                    ring_missing_penalty = gr.Number(label="Penalty if missing in window", value=-5.0)

            with gr.Row():
                with gr.Column():
                    gr.Markdown("#### Steric Lock Region")
                    steric_enabled = gr.Checkbox(label="Enabled", value=True)
                    steric_weight = gr.Slider(label="Weight", minimum=0.0, maximum=5.0, step=0.1, value=1.0)
                    steric_residues = gr.Textbox(
                        label="Bulky residues",
                        value="F,W,Y,H,I,L,V,R,K",
                    )
                    with gr.Row():
                        steric_start_offset = gr.Number(label="Start offset (from acceptor)", value=4, precision=0)
                        steric_end_offset = gr.Number(label="End offset (from acceptor)", value=12, precision=0)
                    steric_require_at_least = gr.Number(label="Require at least N", value=1, precision=0)
        
        with gr.Accordion("Length Constraints", open=True):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("#### Core Length")
                    cl_enabled = gr.Checkbox(label="Enabled", value=True)
                    cl_weight = gr.Slider(label="Weight", minimum=0.0, maximum=5.0, step=0.1, value=1.0)
                    with gr.Row():
                        cl_min = gr.Number(label="Min length", value=10)
                        cl_max = gr.Number(label="Max length", value=40)
                    with gr.Row():
                        cl_opt_min = gr.Number(label="Optimal Min", value=14)
                        cl_opt_max = gr.Number(label="Optimal Max", value=26)

        with gr.Accordion("Sequence Motifs", open=True):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("#### Cleavage Site Motifs")
                    cleave_enabled = gr.Checkbox(label="Enabled", value=True)
                    cleave_weight = gr.Slider(label="Base Weight", minimum=0.0, maximum=5.0, step=0.1, value=0.8)
                    gr.Markdown("Edit motif scores (multiplier for base weight):")
                    # Dataframe for editing motifs
                    cleave_motifs_df = gr.Dataframe(
                        headers=["Motif", "Score"],
                        datatype=["str", "number"],
                        column_count=(2, "fixed"),
                        type="array",
                        label="Cleavage Motifs",
                        interactive=True,
                    )
                
                with gr.Column():
                    gr.Markdown("#### Leader Signature Motifs")
                    ls_enabled = gr.Checkbox(label="Enabled", value=True)
                    ls_weight = gr.Slider(label="Weight", minimum=0.0, maximum=5.0, step=0.1, value=1.2)
                    ls_allow_regex = gr.Checkbox(label="Treat patterns as regex", value=True)
                    ls_patterns_df = gr.Dataframe(
                        headers=["Pattern", "Weight"],
                        datatype=["str", "number"],
                        column_count=(2, "fixed"),
                        type="array",
                        label="Leader Signature Patterns",
                        interactive=True,
                        value=[["Y..P.L", 1.0], ["W..P.L", 0.8], ["LI.LG.A...T.", 0.8]],
                    )

        with gr.Accordion("Leader Residue Prior", open=False):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("#### Penultimate Leader Residue")
                    lp_enabled = gr.Checkbox(label="Enabled", value=True)
                    lp_weight = gr.Slider(label="Weight", minimum=0.0, maximum=5.0, step=0.1, value=1.0)
                    lp_penalty_other = gr.Number(label="Penalty (other residues)", value=-0.5)
                with gr.Column():
                    lp_preferred_df = gr.Dataframe(
                        headers=["Residue", "Weight"],
                        datatype=["str", "number"],
                        column_count=(2, "fixed"),
                        type="array",
                        label="Preferred residues",
                        interactive=True,
                        value=[["T", 1.0]],
                    )
                    lp_allowed_df = gr.Dataframe(
                        headers=["Residue", "Weight"],
                        datatype=["str", "number"],
                        column_count=(2, "fixed"),
                        type="array",
                        label="Allowed rare residues",
                        interactive=True,
                        value=[["I", 0.4], ["V", 0.4]],
                    )

        with gr.Accordion("Advanced Rules", open=False):
            gr.Markdown("*These rules are based on biochemical properties of lasso peptides.*")
            
            with gr.Row():
                with gr.Column():
                    gr.Markdown("#### Hydrophobicity Transition")
                    gr.Markdown("*Core peptides are typically more hydrophobic than leaders.*")
                    hydro_enabled = gr.Checkbox(label="Enabled", value=True)
                    hydro_weight = gr.Slider(label="Weight", minimum=0.0, maximum=5.0, step=0.1, value=0.5)
                    hydro_min_diff = gr.Number(label="Min Difference", value=0.3, minimum=0.0, maximum=2.0)
                
                with gr.Column():
                    gr.Markdown("#### Charge Transition")
                    gr.Markdown("*Leader sequences are typically more charged than cores.*")
                    charge_enabled = gr.Checkbox(label="Enabled", value=True)
                    charge_weight = gr.Slider(label="Weight", minimum=0.0, maximum=5.0, step=0.1, value=0.4)
            
            with gr.Row():
                with gr.Column():
                    gr.Markdown("#### PTM Context Hints")
                    gr.Markdown("*Optional hints from tailoring enzymes; disabled by default.*")
                    ptm_enabled = gr.Checkbox(label="Enabled", value=False)
                    ptm_weight = gr.Slider(label="Weight", minimum=0.0, maximum=5.0, step=0.1, value=0.6)
                    ptm_kinase = gr.Textbox(label="Kinase prefer residues", value="S,T")
                    ptm_gt = gr.Textbox(label="Glycosyltransferase prefer residues", value="S,T,Y")
                    ptm_acyl = gr.Textbox(label="Acyltransferase prefer residues", value="K,S,T")

        with gr.Row():
            save_rules_btn = gr.Button("Save Rules", variant="primary")
            export_rules_btn = gr.Button("Export to JSON")
            import_rules_file = gr.File(label="Import rules JSON", file_types=[".json"])
        
        rules_status = gr.Markdown()
        rules_json_output = gr.Code(label="Rules JSON", language="json", visible=False)

        # --- Helper Functions ---

        def parse_list(s):
            return [x.strip() for x in s.split(",") if x.strip()]

        def parse_int_list(s):
            return [int(x.strip()) for x in s.split(",") if x.strip() and x.strip().isdigit()]

        def format_list(l):
            return ",".join(str(x) for x in l)

        def load_ui_values():
            """Load values from service into UI components."""
            rule_set = rules_service.get_rules()
            # Access underlying rules dict (Dict[str, RuleParameter])
            r = rule_set.rules
            
            # Helper to get rule safely
            def get_r(name):
                rule = r.get(name)
                if rule is None:
                    raise KeyError(f"Missing rule configuration: {name}")
                return rule

            # Core Start
            cs = get_r("core_start_residue")
            scores_dict = cs.parameters.scores
            scores_data = [[k, v] for k, v in scores_dict.items()]
            cs_vals = [
                cs.enabled,
                cs.weight,
                scores_data,
                cs.parameters.penalty_other,
                format_list(cs.parameters.hard_forbidden),
                cs.parameters.penalty_forbidden,
            ]
            
            # Ring Closure
            ring = get_r("ring_closure_residue")
            ring_weights_data = [[int(k), v] for k, v in ring.parameters.position_weights.items()]
            ring_vals = [
                ring.enabled,
                ring.weight,
                format_list(ring.parameters.allowed_positions),
                format_list(ring.parameters.residues),
                ring_weights_data,
                ring.parameters.penalty_if_missing_in_allowed_window,
            ]

            # Steric Lock Region
            steric = get_r("steric_lock_region")
            steric_vals = [
                steric.enabled,
                steric.weight,
                format_list(steric.parameters.bulky_residues),
                steric.parameters.relative_to_acceptor.get("start_offset", 0),
                steric.parameters.relative_to_acceptor.get("end_offset", 0),
                steric.parameters.require_at_least,
            ]

            # Core Length
            cl = get_r("core_length")
            cl_vals = [
                cl.enabled, cl.weight,
                cl.parameters.min, cl.parameters.max,
                cl.parameters.optimal_min, cl.parameters.optimal_max
            ]
            
            # Cleavage Motifs
            cm = get_r("cleavage_motif")
            motifs_dict = cm.parameters.motifs
            motifs_data = [[k, v] for k, v in motifs_dict.items()] if isinstance(motifs_dict, dict) else []
            cm_vals = [cm.enabled, cm.weight, motifs_data]
            
            # Leader Signature Motifs
            ls = get_r("leader_signature_motifs")
            ls_patterns = [[k, v] for k, v in ls.parameters.patterns.items()]
            ls_vals = [ls.enabled, ls.weight, ls.parameters.allow_regex, ls_patterns]

            # Leader Penultimate Residue
            lp = get_r("leader_penultimate_residue")
            lp_preferred = [[k, v] for k, v in lp.parameters.preferred.items()]
            lp_allowed = [[k, v] for k, v in lp.parameters.allowed_rare.items()]
            lp_vals = [lp.enabled, lp.weight, lp.parameters.penalty_other, lp_preferred, lp_allowed]
            
            # Advanced Rules
            hydro = get_r("hydrophobicity_transition")
            hydro_vals = [hydro.enabled, hydro.weight, hydro.parameters.min_difference]
            
            charge = get_r("charge_transition")
            charge_vals = [charge.enabled, charge.weight]

            ptm = get_r("ptm_context_hints")
            ptm_vals = [
                ptm.enabled,
                ptm.weight,
                format_list(ptm.parameters.if_kinase_present.get("prefer_residues", [])),
                format_list(ptm.parameters.if_glycosyltransferase_present.get("prefer_residues", [])),
                format_list(ptm.parameters.if_acyltransferase_present.get("prefer_residues", [])),
            ]

            return (
                cs_vals
                + ring_vals
                + steric_vals
                + cl_vals
                + cm_vals
                + ls_vals
                + lp_vals
                + hydro_vals
                + charge_vals
                + ptm_vals
            )

        def save_ui_values(
            cs_en, cs_w, cs_scores, cs_penalty_other, cs_hard_forbidden, cs_penalty_forbidden,
            ring_en, ring_w, ring_allowed_pos, ring_res, ring_pos_weights, ring_missing_penalty,
            steric_en, steric_w, steric_residues, steric_start, steric_end, steric_require,
            cl_en, cl_w, cl_min, cl_max, cl_opt_min, cl_opt_max,
            cm_en, cm_w, cm_motifs,
            ls_en, ls_w, ls_allow_regex, ls_patterns,
            lp_en, lp_w, lp_penalty_other, lp_preferred, lp_allowed,
            hydro_en, hydro_w, hydro_min_diff,
            charge_en, charge_w,
            ptm_en, ptm_w, ptm_kinase, ptm_gt, ptm_acyl
        ):
            """Collect UI values and save to service."""
            
            # Helper to convert dataframe to dict
            motifs_dict = {}
            if cm_motifs is not None and hasattr(cm_motifs, "__iter__"):
                for row in cm_motifs:
                    if len(row) >= 2 and row[0] and str(row[1]):
                        motifs_dict[str(row[0])] = float(row[1])

            cs_scores_dict = {}
            if cs_scores is not None and hasattr(cs_scores, "__iter__"):
                for row in cs_scores:
                    if len(row) >= 2 and row[0] and str(row[1]):
                        cs_scores_dict[str(row[0])] = float(row[1])

            ring_weights = {}
            if ring_pos_weights is not None and hasattr(ring_pos_weights, "__iter__"):
                for row in ring_pos_weights:
                    if len(row) >= 2 and str(row[0]) and str(row[1]):
                        ring_weights[str(int(float(row[0])))] = float(row[1])

            ls_patterns_dict = {}
            if ls_patterns is not None and hasattr(ls_patterns, "__iter__"):
                for row in ls_patterns:
                    if len(row) >= 2 and row[0] and str(row[1]):
                        ls_patterns_dict[str(row[0])] = float(row[1])

            lp_preferred_dict = {}
            if lp_preferred is not None and hasattr(lp_preferred, "__iter__"):
                for row in lp_preferred:
                    if len(row) >= 2 and row[0] and str(row[1]):
                        lp_preferred_dict[str(row[0])] = float(row[1])

            lp_allowed_dict = {}
            if lp_allowed is not None and hasattr(lp_allowed, "__iter__"):
                for row in lp_allowed:
                    if len(row) >= 2 and row[0] and str(row[1]):
                        lp_allowed_dict[str(row[0])] = float(row[1])

            # Construct updates for each rule
            updates = {
                "core_start_residue": {
                    "enabled": cs_en, "weight": cs_w,
                    "parameters": {
                        "scores": cs_scores_dict,
                        "penalty_other": float(cs_penalty_other),
                        "hard_forbidden": parse_list(cs_hard_forbidden),
                        "penalty_forbidden": float(cs_penalty_forbidden),
                    }
                },
                "ring_closure_residue": {
                    "enabled": ring_en, "weight": ring_w,
                    "parameters": {
                        "allowed_positions": parse_int_list(ring_allowed_pos),
                        "residues": parse_list(ring_res),
                        "position_weights": ring_weights,
                        "penalty_if_missing_in_allowed_window": float(ring_missing_penalty),
                    }
                },
                "steric_lock_region": {
                    "enabled": steric_en,
                    "weight": steric_w,
                    "parameters": {
                        "bulky_residues": parse_list(steric_residues),
                        "relative_to_acceptor": {
                            "start_offset": int(steric_start),
                            "end_offset": int(steric_end),
                        },
                        "require_at_least": int(steric_require),
                    },
                },
                "core_length": {
                    "enabled": cl_en, "weight": cl_w,
                    "parameters": {
                        "min": int(cl_min), "max": int(cl_max),
                        "optimal_min": int(cl_opt_min), "optimal_max": int(cl_opt_max)
                    }
                },
                "cleavage_motif": {
                    "enabled": cm_en, "weight": cm_w,
                    "parameters": {"motifs": motifs_dict}
                },
                "leader_signature_motifs": {
                    "enabled": ls_en, "weight": ls_w,
                    "parameters": {
                        "patterns": ls_patterns_dict,
                        "allow_regex": bool(ls_allow_regex),
                    }
                },
                "leader_penultimate_residue": {
                    "enabled": lp_en,
                    "weight": lp_w,
                    "parameters": {
                        "preferred": lp_preferred_dict,
                        "allowed_rare": lp_allowed_dict,
                        "penalty_other": float(lp_penalty_other),
                    }
                },
                "hydrophobicity_transition": {
                    "enabled": hydro_en, "weight": hydro_w,
                    "parameters": {"min_difference": float(hydro_min_diff)}
                },
                "charge_transition": {
                    "enabled": charge_en, "weight": charge_w,
                    "parameters": {"leader_more_charged": True}
                },
                "ptm_context_hints": {
                    "enabled": ptm_en,
                    "weight": ptm_w,
                    "parameters": {
                        "if_kinase_present": {"prefer_residues": parse_list(ptm_kinase)},
                        "if_glycosyltransferase_present": {"prefer_residues": parse_list(ptm_gt)},
                        "if_acyltransferase_present": {"prefer_residues": parse_list(ptm_acyl)},
                    },
                },
            }
            
            # Apply updates
            for name, data in updates.items():
                rules_service.update_rule(name, data)
                
            rules_service.save_rules()
            return "Rules saved successfully."

        def apply_preset_handler(preset_name):
            if preset_name == "(keep current)":
                return load_ui_values()
            try:
                rules_service.apply_preset(preset_name)
            except FileNotFoundError:
                return load_ui_values()
            return load_ui_values()

        def reset_handler():
            rules_service.reset_to_default()
            return load_ui_values()

        def export_handler():
            return gr.update(value=rules_service.get_rules().model_dump_json(indent=2), visible=True)

        def import_handler(file_obj):
            if file_obj is None:
                return "No rules file selected.", gr.update(visible=False), *load_ui_values()
            try:
                rules_service.import_rules(Path(file_obj.name))
            except Exception as exc:
                return f"Failed to import rules: {exc}", gr.update(visible=False), *load_ui_values()
            return "Rules imported successfully.", gr.update(visible=False), *load_ui_values()

        # --- Event Wiring ---
        
        all_inputs = [
            cs_enabled, cs_weight, cs_scores_df, cs_penalty_other, cs_hard_forbidden, cs_penalty_forbidden,
            ring_enabled, ring_weight, ring_allowed_positions, ring_residues, ring_position_weights, ring_missing_penalty,
            steric_enabled, steric_weight, steric_residues, steric_start_offset, steric_end_offset, steric_require_at_least,
            cl_enabled, cl_weight, cl_min, cl_max, cl_opt_min, cl_opt_max,
            cleave_enabled, cleave_weight, cleave_motifs_df,
            ls_enabled, ls_weight, ls_allow_regex, ls_patterns_df,
            lp_enabled, lp_weight, lp_penalty_other, lp_preferred_df, lp_allowed_df,
            hydro_enabled, hydro_weight, hydro_min_diff,
            charge_enabled, charge_weight,
            ptm_enabled, ptm_weight, ptm_kinase, ptm_gt, ptm_acyl
        ]

        save_rules_btn.click(save_ui_values, inputs=all_inputs, outputs=[rules_status])
        load_preset_btn.click(apply_preset_handler, inputs=[preset_select], outputs=all_inputs)
        reset_btn.click(reset_handler, outputs=all_inputs)
        export_rules_btn.click(export_handler, outputs=[rules_json_output])
        import_rules_file.change(
            import_handler,
            inputs=[import_rules_file],
            outputs=[rules_status, rules_json_output] + all_inputs,
        )
        rules_tab.select(load_ui_values, outputs=all_inputs)
        # Return dict for tab context
        return {
            "tab": None,
            "loader": load_ui_values,
            "inputs": all_inputs
        }
