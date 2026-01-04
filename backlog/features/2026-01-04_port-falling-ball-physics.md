---
id: "2026-01-04_port-falling-ball-physics"
title: "Port falling ball physics simulation from MATLAB to Python"
status: "Completed"
priority: "High"
created: "2026-01-04"
last_updated: "2026-01-04"
category: "features"
related_cips: ["0001"]
owner: ""
dependencies: ["2026-01-04_setup-python-package"]
tags:
- backlog
- python
- physics
- porting
---

# Task: Port falling ball physics simulation from MATLAB to Python

> **Note**: Backlog tasks are DOING the work defined in CIPs (HOW).  
> Use `related_cips` to link to CIPs. Don't link directly to requirements (bottom-up pattern).

## Description

Port the falling ball physics simulation and interaction logic from MATLAB to Python. This includes the ball dynamics, keyboard interaction callbacks, and probabilistic outcome demonstration.

This task implements Step 3 of CIP-0001's migration plan.

## Acceptance Criteria

- [ ] Port `simulateFallingBall.m` → `kappenball/falling_ball.py::simulate()`
- [ ] Port keyboard interaction callbacks (`fallPressKey.m`, `fallReleaseKey.m`)
- [ ] Implement ball physics (gravity, initial conditions, trajectory)
- [ ] Create unit tests for physics correctness
- [ ] Document function signatures with docstrings
- [ ] Verify probabilistic behavior matches expected outcomes

## Implementation Notes

**MATLAB → Python translation:**
- MATLAB's keyboard callbacks → Python event handling (will integrate with Jupyter widgets later)
- MATLAB's graphics callbacks → Python equivalents

**Key physics to preserve:**
- Falling ball trajectory under gravity
- Keyboard influence on ball dynamics (if applicable)
- Probabilistic outcomes demonstration

**Interaction strategy:**
- Initially implement physics without UI (pure functions)
- UI/interaction layer will be added in notebook task
- Design for testability: separate physics from visualization/interaction

**Files to reference:**
- `matlab/simulateFallingBall.m`
- `matlab/fallPressKey.m`
- `matlab/fallReleaseKey.m`
- `matlab/demKappenBall.m` (for context)

## Related

- CIP: 0001
- Implements: Step 3 of CIP-0001 implementation plan
- Depends on: 2026-01-04_setup-python-package

## Progress Updates

### 2026-01-04
Task created. Lower priority than billiards since entropy demo is more central.

