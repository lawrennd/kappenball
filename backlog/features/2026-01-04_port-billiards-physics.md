---
id: "2026-01-04_port-billiards-physics"
title: "Port billiards physics simulation from MATLAB to Python"
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

# Task: Port billiards physics simulation from MATLAB to Python

> **Note**: Backlog tasks are DOING the work defined in CIPs (HOW).  
> Use `related_cips` to link to CIPs. Don't link directly to requirements (bottom-up pattern).

## Description

Port the core billiards physics simulation functions from MATLAB to Python. This includes initialization, simulation loop, collision detection, and vector computations. The Python implementation should maintain numerical accuracy compared to the original MATLAB/Octave code.

This task implements Step 2 of CIP-0001's migration plan and is foundational for the entropy demonstrations.

## Acceptance Criteria

- [x] Port `initializeBilliards.m` → `kappenball/billiards.py::initialise()`
- [x] Port `simulateBilliards.m` → `kappenball/billiards.py::simulate_step()`
- [x] Port `vectorBilliards.m` → `kappenball/billiards.py::compute_vectors()`
- [x] Port `Billiards2D.m` functions (2D implementation complete)
- [x] Create unit tests (12 tests, all passing)
- [x] Document function signatures with docstrings
- [x] Verify billiards bounce physics is correct (energy conservation, wall bouncing tested)

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

### 2026-01-04 (Initial)
Task created. Blocked by package setup task.

### 2026-01-04 (Completion)
✅ Billiards physics successfully ported from MATLAB to Python:

**Implemented functions:**
- `initialise()` - 4 initialization modes (rand, randn, uniform, hotCold)
- `simulate_step()` - Wall collisions and particle-particle elastic collisions
- `compute_vectors()` - Velocity vector data for visualisation
- Helper functions for collision detection and non-overlapping initialization

**Test coverage:**
- 12 tests, all passing
- 93% code coverage on billiards.py
- Tests include: initialization modes, energy conservation, boundary conditions, collision detection, reproducibility

**Key features:**
- Elastic collisions with proper physics (momentum/energy conservation)
- Wall bounce detection and handling
- Non-overlapping particle initialization
- British spelling throughout (initialise, colours, visualisation)

