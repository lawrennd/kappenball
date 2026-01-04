---
id: "2026-01-04_entropy-billiards-notebook"
title: "Create interactive entropy billiards Jupyter notebook"
status: "Proposed"
priority: "High"
created: "2026-01-04"
last_updated: "2026-01-04"
category: "features"
related_cips: ["0001"]
owner: ""
dependencies: ["2026-01-04_create-visualization-layer"]
tags:
- backlog
- jupyter
- notebook
- interactive
---

# Task: Create interactive entropy billiards Jupyter notebook

> **Note**: Backlog tasks are DOING the work defined in CIPs (HOW).  
> Use `related_cips` to link to CIPs. Don't link directly to requirements (bottom-up pattern).

## Description

Create an interactive Jupyter notebook that demonstrates entropy increase through billiards simulation. The notebook should use ipywidgets for interactive controls, allowing students to experiment with parameters and observe uncertainty concepts in real-time.

This task implements Step 5 of CIP-0001 and delivers the primary educational artifact for entropy demonstrations.

## Acceptance Criteria

- [ ] Create `notebooks/entropy_billiards.ipynb`
- [ ] Add markdown cells explaining entropy concepts for lecture context
- [ ] Implement ipywidgets controls (particle count, reset, start/stop, demo mode)
- [ ] Add toggle controls for velocities, trajectories, and other visual elements
- [ ] Integrate visualization layer with interactive matplotlib (%matplotlib widget)
- [ ] Provide example parameters that demonstrate entropy increase clearly
- [ ] Test that notebook runs from fresh Python environment
- [ ] Ensure "from running to modifying in minutes" tenet is met

## Implementation Notes

**Notebook structure:**
1. **Introduction** (markdown) - Explain entropy and uncertainty concepts
2. **Setup** (code) - Import kappenball package, configure visualization
3. **Interactive Demo** (code + widgets) - Main simulation with controls
4. **Exploration** (markdown + code) - Guided experiments for students
5. **Extension** (markdown) - Ideas for student modifications

**ipywidgets controls:**
- Slider: Number of particles (5-50)
- Slider: Initial configuration (ordered ↔ random)
- Buttons: Start/Stop, Reset, Demo Mode
- Toggles: Show velocities, Show trajectories, Show entropy plot
- Dropdown: Predefined scenarios (collision cascade, entropy demo, etc.)

**Interactivity approach:**
- Use `ipympl` (`%matplotlib widget`) for interactive matplotlib
- Update visualization in response to widget changes
- Real-time animation with FuncAnimation
- Responsive feedback (target < 100ms for control interactions)

**Educational content:**
- Explain entropy as measure of disorder
- Show how ordered initial states evolve to disordered states
- Connect visual observations to thermodynamic principles
- Provide "what if" prompts for exploration

**Files to reference:**
- `matlab/demEntropyBilliards.m` (original demo structure)
- `matlab/entropyBilliards.m` (entropy calculations)

## Related

- CIP: 0001
- Implements: Step 5 of CIP-0001 implementation plan
- Supports: REQ-0002 (real-time interactivity), REQ-0003 (visual uncertainty)
- Depends on: 2026-01-04_create-visualization-layer

## Progress Updates

### 2026-01-04
Task created. This is a key deliverable - the primary educational artifact.

