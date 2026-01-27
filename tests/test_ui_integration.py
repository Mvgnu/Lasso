"""
UI Integration Tests for Lasso Workbench.

Tests that verify complete workflows through UI components.
"""

class TestAppImports:
    """Tests that main app imports successfully."""
    
    def test_create_app_importable(self):
        """create_app can be imported without errors."""
        from lasso_workbench.app import create_app
        assert create_app is not None
    
    def test_all_tabs_importable(self):
        """All tab creation functions can be imported."""
        from lasso_workbench.ui.tabs.dataset import create_dataset_tab
        from lasso_workbench.ui.tabs.rules import create_rules_tab
        from lasso_workbench.ui.tabs.pipeline import create_pipeline_tab
                
        assert all([
            create_dataset_tab,
            create_rules_tab,
            create_pipeline_tab,
        ])
