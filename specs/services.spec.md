---
name: services
description: "Business logic services coordinating between core and UI"
paths:
  types: lasso_workbench/services/__init__.py
  modules:
    - lasso_workbench/services/data_service.py
    - lasso_workbench/services/pipeline_service.py
    - lasso_workbench/services/rules_service.py
  tests: tests/test_services.py
exports:
  - DataService
  - PipelineService
  - RulesService
consumes:
  - core.prediction
  - core.dataset
  - pipeline.semantic_pipeline
  - schemas.config
verification:
  typecheck: "python -m py_compile lasso_workbench/services/*.py"
  test: "pytest tests/test_services.py -v"
---

# Services Domain

The **services** domain owns business logic that coordinates between core algorithms and UI.

## Architecture

```
UI (Gradio tabs)
       ↓
   Services  ← Orchestration layer
       ↓
Core / Pipeline
```

## Service Files

| Service | Purpose |
|---------|---------|
| `data_service.py` | Data loading/saving, import/export formats |
| `pipeline_service.py` | Pipeline execution wrapper for UI |
| `rules_service.py` | Rules persistence, presets, validation |

## Responsibilities

Services are responsible for:
- **Orchestration**: Coordinating multiple core components
- **State management**: Maintaining application state
- **I/O handling**: File reading/writing, format conversion
- **Error handling**: User-friendly error messages

## Data Boundaries

- `data/precursors/*` is read-only from the UI except for `precursor_proteins_curated.tsv`.
- Curated additions are done via `scripts/add_validated_precursors.py`.
- Default precursor reference set in the UI is `data/precursors/precursor_proteins_verified.faa`.

## Constraints

> [!WARNING]
> Services should NOT contain core algorithms. If you're writing complex logic, it belongs in `core/`.
