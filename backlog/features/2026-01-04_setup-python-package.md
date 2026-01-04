---
id: "2026-01-04_setup-python-package"
title: "Set up Python package structure for Kappenball"
status: "Completed"
priority: "High"
created: "2026-01-04"
last_updated: "2026-01-04"
category: "infrastructure"
related_cips: ["0001"]
owner: ""
dependencies: []
tags:
- backlog
- python
- infrastructure
- setup
---

# Task: Set up Python package structure for Kappenball

> **Note**: Backlog tasks are DOING the work defined in CIPs (HOW).  
> Use `related_cips` to link to CIPs. Don't link directly to requirements (bottom-up pattern).

## Description

Create the foundational Python package structure for Kappenball, including directory layout, dependency management, and basic package configuration. This establishes the framework for porting the MATLAB physics simulations to Python.

This task implements the first step of CIP-0001's migration plan.

## Acceptance Criteria

- [x] Create `kappenball/` package directory with `__init__.py`
- [x] Create subdirectories: `notebooks/`, `tests/`
- [x] Create `requirements.txt` with minimal core dependencies (numpy, matplotlib, jupyter, ipywidgets, ipympl)
- [x] Create `pyproject.toml` for package installation
- [x] Create `.gitignore` for Python-specific files (*.pyc, __pycache__, .pytest_cache, etc.)
- [x] Verify package can be installed in development mode: `pip install -e .`
- [x] Create basic README.md with installation instructions

## Implementation Notes

**Package structure:**
```
kappenball/
├── kappenball/           # Python package
│   ├── __init__.py
│   ├── billiards.py      (empty, for next task)
│   ├── falling_ball.py   (empty, for next task)
│   ├── physics.py        (empty, for next task)
│   └── visualization.py  (empty, for next task)
├── notebooks/            # Jupyter notebooks
├── tests/                # Unit tests
├── requirements.txt      # Dependencies
├── setup.py             # Package configuration
└── README.md            # Installation & usage
```

**Minimal dependencies** (aligned with CIP-0001):
- numpy >= 1.21
- matplotlib >= 3.5
- jupyter >= 1.0
- ipywidgets >= 8.0
- ipympl >= 0.9

**Development dependencies** (requirements-dev.txt):
- pytest >= 7.0
- pytest-cov
- black
- ruff

## Related

- CIP: 0001
- Implements: Step 1 of CIP-0001 implementation plan

## Progress Updates

### 2026-01-04 (Initial)
Task created as first step of Python migration.

### 2026-01-04 (Completion)
✅ Package structure complete:
- Python package with 4 modules (billiards, falling_ball, physics, visualisation - British spelling)
- pyproject.toml with minimal dependencies
- requirements.txt and requirements-dev.txt
- tests/ directory with passing tests (3/3 passed, 63% coverage)
- Updated .gitignore for Python artifacts
- Updated README with installation instructions
- Package successfully installs with `pip install -e .`

