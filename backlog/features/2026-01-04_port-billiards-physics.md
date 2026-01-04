---
id: "2026-01-04_port-billiards-physics"
title: "Port billiards physics simulation from MATLAB to Python"
status: "Proposed"
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

# Task: Port billiards physics simulation from MATLAB to Python

> **Note**: Backlog tasks are DOING the work defined in CIPs (HOW).  
> Use `related_cips` to link to CIPs. Don't link directly to requirements (bottom-up pattern).

## Description

Port the core billiards physics simulation functions from MATLAB to Python. This includes initialization, simulation loop, collision detection, and vector computations. The Python implementation should maintain numerical accuracy compared to the original MATLAB/Octave code.

This task implements Step 2 of CIP-0001's migration plan and is foundational for the entropy demonstrations.

## Acceptance Criteria

- [ ] Port `initializeBilliards.m` → `kappenball/billiards.py::initialize()`
- [ ] Port `simulateBilliards.m` → `kappenball/billiards.py::simulate()`
- [ ] Port `vectorBilliards.m` → `kappenball/billiards.py::compute_vectors()`
- [ ] Port `Billiards2D.m` and `Billiards3D.m` functions as needed
- [ ] Create unit tests comparing Python output to Octave reference (if available)
- [ ] Document function signatures with docstrings
- [ ] Verify billiards bounce physics is correct (energy conservation, reflection angles)

## Implementation Notes

**MATLAB → Python translation patterns:**
- MATLAB arrays → NumPy arrays
- MATLAB functions → Python functions with numpy operations
- MATLAB's `rand()` → `numpy.random.rand()`
- MATLAB's matrix operations → NumPy equivalents

**Key physics to preserve:**
- 2D/3D collision detection
- Elastic collisions (momentum and energy conservation)
- Boundary reflections
- Particle initialization (positions, velocities)

**Testing approach:**
- If Octave is working, generate reference outputs from MATLAB code
- Compare Python outputs with tolerance (e.g., `numpy.testing.assert_allclose`)
- If Octave doesn't work, verify physics principles (energy conservation, etc.)

**Files to reference:**
- `matlab/initializeBilliards.m`
- `matlab/simulateBilliards.m`
- `matlab/vectorBilliards.m`
- `matlab/Billiards2D.m`
- `matlab/Billiards3D.m`

## Related

- CIP: 0001
- Implements: Step 2 of CIP-0001 implementation plan
- Depends on: 2026-01-04_setup-python-package

## Progress Updates

### 2026-01-04
Task created. Blocked by package setup task.

