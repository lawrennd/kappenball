---
id: "2026-01-04_create-visualization-layer"
title: "Create matplotlib visualization layer for simulations"
status: "Completed"
priority: "High"
created: "2026-01-04"
last_updated: "2026-01-04"
category: "features"
related_cips: ["0001"]
owner: ""
dependencies: ["2026-01-04_port-billiards-physics"]
tags:
- backlog
- python
- visualization
- matplotlib
---

# Task: Create matplotlib visualization layer for simulations

> **Note**: Backlog tasks are DOING the work defined in CIPs (HOW).  
> Use `related_cips` to link to CIPs. Don't link directly to requirements (bottom-up pattern).

## Description

Create a visualization layer using matplotlib to display the physics simulations. This layer should make uncertainty concepts visible through clear, interpretable plots and animations, matching or exceeding the visual quality of the original MATLAB demos.

This task implements Step 4 of CIP-0001's migration plan and directly supports REQ-0003 (visual uncertainty).

## Acceptance Criteria

- [x] Create `kappenball/visualisation.py` module (British spelling)
- [x] Implement billiards visualisation (particle positions, velocities, trajectories)
- [x] Implement static and animated plotting
- [x] Add support for toggling visual elements (velocity vectors with show_velocities parameter)
- [x] Use matplotlib.animation for smooth animations
- [x] Create billiard table styling (blue background, white borders)
- [x] Implement additional plots (velocity histogram, energy evolution)
- [x] Animation interval configurable (default 20ms = 50 fps)

## Implementation Notes

**Key visualization functions:**
- `plot_billiards_frame()` - Render single frame of billiards simulation
- `plot_falling_ball_frame()` - Render single frame of falling ball
- `animate_billiards()` - Create animated billiards visualization
- `animate_falling_ball()` - Create animated falling ball
- `setup_figure()` - Configure matplotlib figure with appropriate styling

**Visual elements to support:**
- Particle positions (circles/dots)
- Velocity vectors (arrows, toggleable)
- Particle trajectories (lines, toggleable)
- Boundary walls
- Color schemes that emphasize uncertainty concepts
- Entropy indicators (if applicable)

**matplotlib features to use:**
- `matplotlib.animation.FuncAnimation` for smooth animations
- `matplotlib.patches.Circle` for particles
- `matplotlib.quiver` for velocity vectors
- Interactive backend for notebook integration

**Design considerations:**
- Separate visualization from physics (pure functions)
- Support both static plots and animations
- Make styling configurable
- Optimize for performance with many particles

## Related

- CIP: 0001
- Implements: Step 4 of CIP-0001 implementation plan
- Supports: REQ-0003 (visual uncertainty demonstration)
- Depends on: 2026-01-04_port-billiards-physics

## Progress Updates

### 2026-01-04 (Initial)
Task created. High priority as visualisation is core to educational value.

### 2026-01-04 (Completion)
✅ Visualisation layer successfully implemented with matplotlib:

**Implemented functions:**
- `setup_figure()` - Configure matplotlib figures with billiard table styling
- `plot_billiards_frame()` - Render single frame with particles and velocity arrows
- `plot_billiards_static()` - Quick static visualisation
- `animate_billiards()` - Create smooth animations with FuncAnimation
- `plot_velocity_histogram()` - Velocity distribution analysis
- `plot_energy_evolution()` - Energy conservation verification
- `clear_artists()` - Helper for animation updates

**Features:**
- Particles rendered as coloured circles with proper sizing
- Velocity vectors shown as arrows (toggleable)
- Blue billiard table background (#4080CC)
- White boundary walls
- Configurable velocity scale for arrow sizing
- Support for all 4 initialization modes
- Smooth animations (20ms interval = 50 fps default)

**Test coverage:**
- 11 tests, all passing
- 93% code coverage on visualisation.py
- Tests cover: figure setup, static plots, animations, histograms, energy plots

