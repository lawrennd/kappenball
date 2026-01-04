---
id: "2026-01-04_setup-python-package"
title: "Set up Python package structure for Kappenball"
status: "Proposed"
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

- [ ] Create `kappenball/` package directory with `__init__.py`
- [ ] Create subdirectories: `notebooks/`, `tests/`
- [ ] Create `requirements.txt` with minimal core dependencies (numpy, matplotlib, jupyter, ipywidgets, ipympl)
- [ ] Create `pyproject.toml` for package installation
- [ ] Create `.gitignore` for Python-specific files (*.pyc, __pycache__, .pytest_cache, etc.)
- [ ] Verify package can be installed in development mode: `pip install -e .`
- [ ] Create basic README.md with installation instructions

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

### 2026-01-04
Task created as first step of Python migration.

