"""
Main entry point for Lasso Workbench.

Usage:
    python -m lasso_workbench           # Launch UI
"""

import sys
import argparse


def main():
    """Launch the Lasso Workbench application."""
    parser = argparse.ArgumentParser(description="Lasso Workbench - AI-powered lasso peptide discovery")
    parser.add_argument(
        "--port",
        type=int,
        default=7860,
        help="Port to run the server on (default: 7860)",
    )
    args = parser.parse_args()
    
    from lasso_workbench.app import APP_THEME, CUSTOM_CSS, create_app
    app = create_app()
    
    app.launch(
        server_name="127.0.0.1",
        server_port=args.port,
        share=False,
        theme=APP_THEME,
        css=CUSTOM_CSS,
    )


if __name__ == "__main__":
    sys.exit(main() or 0)
