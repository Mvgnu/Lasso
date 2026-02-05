"""
Lasso Workbench - Main Application

Gradio-based GUI for lasso peptide discovery and dataset curation.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import gradio as gr
import pandas as pd

from lasso_workbench import DATA_DIR, __version__
from lasso_workbench.services.rules_service import RulesService
from lasso_workbench.services.data_service import DataService
from lasso_workbench.services.pipeline_service import PipelineService
# UI Tabs
from lasso_workbench.ui.tabs.rules import create_rules_tab
from lasso_workbench.ui.tabs.dataset import create_dataset_tab
from lasso_workbench.ui.tabs.pipeline import create_pipeline_tab

logger = logging.getLogger(__name__)

# Custom CSS for clean system fonts
CUSTOM_CSS = """
* {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif !important;
}
code, pre, .code, textarea.code {
    font-family: "SF Mono", Monaco, Consolas, "Liberation Mono", "Courier New", monospace !important;
}
"""

APP_THEME = gr.themes.Base(
    primary_hue="slate",
    secondary_hue="slate",
    neutral_hue="gray",
)

def create_app() -> gr.Blocks:
    """Create the main Gradio application."""
    
    # Initialize Services
    rules_service = RulesService()
    data_service = DataService()
    pipeline_service = PipelineService()
    
    with gr.Blocks(
        title="Lasso Workbench",
    ) as app:
        
        gr.Markdown(f"""
        # Lasso Workbench v{__version__}
        
        Lasso peptide discovery, prediction, and dataset curation.
        """)
        
        with gr.Tabs():
            
            # ================================================================
            # Tab 1: Semantic Pipeline (Refactored)
            # ================================================================
            create_pipeline_tab(pipeline_service, rules_service)
            
            # ================================================================
            # Tab 2: Rule Editor
            # ================================================================
            rules_tab_data = create_rules_tab(rules_service)
            
            
            # ================================================================
            # Tab 3: Dataset Manager
            # ================================================================
            dataset_view = create_dataset_tab(data_service)
            

            # ================================================================
            # Tab 8: About
            # ================================================================
            with gr.TabItem("About"):
                gr.Markdown(f"""
                # Lasso Workbench v{__version__}
                
                A GUI for lasso peptide discovery and dataset curation.
                
                ## Features
                
                - **Pipeline:** ESM-2 semantic discovery (6-frame ORFs + similarity)
                - **Rules:** Configure prediction parameters
                - **Dataset:** View, edit, and import/export datasets
                
                ## Tabs Overview
                
                | Tab | Purpose |
                |-----|---------|
                | Pipeline | Upload GBK, score candidates |
                | Rules | Configure prediction parameters |
                | Dataset | Manage datasets |
                
                ## Pipeline Overview
                
                The **Pipeline** tab runs the semantic precursor discovery workflow:
                
                1. **Upload GBK files** (MiBIG, antiSMASH, or custom GenBank)
                2. **6-frame ORF extraction** finds ALL potential coding sequences
                3. **ESM-2 embedding** computes semantic similarity to validated precursors
                4. **All candidates scored** (not just top N)
                5. **Export results** as TSV or add high-scoring candidates to dataset
                
                ## Getting Started
                
                1. Use **Pipeline** for ESM-2 discovery from GBK files
                2. Adjust rules in **Rules**
                3. Manage entries in **Dataset**
                
                ## Documentation
                
                See README.md and specs/*.md for technical details.
                """)
         
        # Load initial data
        app.load(data_service.load_dataset, outputs=[dataset_view])
        app.load(rules_tab_data["loader"], outputs=rules_tab_data["inputs"])
    
    return app


def launch_app():
    """Launch the Lasso Workbench application."""
    app = create_app()
    app.launch(
        server_name="127.0.0.1",
        server_port=7860,
        share=False,
        inbrowser=True,
        theme=APP_THEME,
        css=CUSTOM_CSS,
    )


if __name__ == "__main__":
    launch_app()
