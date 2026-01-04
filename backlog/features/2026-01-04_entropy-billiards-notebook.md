---
id: "2026-01-04_entropy-billiards-notebook"
title: "Create interactive entropy billiards Jupyter notebook"
status: "Completed"
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

- [x] Create `notebooks/entropy_billiards.ipynb`
- [x] Add markdown cells explaining entropy concepts for lecture context
- [x] Demonstrate entropy increase with before/after comparisons
- [x] Create animated visualisations showing real-time evolution
- [x] Integrate visualisation layer with interactive matplotlib (%matplotlib widget)
- [x] Provide multiple examples (static plots, animations, velocity analysis, energy conservation)
- [x] Include exploration challenges for students
- [x] Educational content explaining Second Law, Maxwell-Boltzmann distribution, thermalization

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

### 2026-01-04 (Initial)
Task created. This is a key deliverable - the primary educational artifact.

### 2026-01-04 (Completion)
✅ Interactive entropy billiards Jupyter notebook successfully created!

**Notebook structure (14 cells):**
1. **Introduction** - Explains entropy, Second Law, Maxwell-Boltzmann distribution
2. **Setup** - Import modules, configure %matplotlib widget
3. **Quick Static Example** - Hot/cold initial state visualisation
4. **Entropy Increase Demo** - Before/after comparison showing disorder increase
5. **Interactive Animation** - Real-time visualisation of entropy increase
6. **Velocity Distribution** - Shows thermalization and equilibrium
7. **Energy Conservation** - Verifies physics accuracy
8. **Exploration Challenges** - Student exercises and experiments

**Educational features:**
- Clear explanations of uncertainty and entropy concepts
- Visual demonstrations make abstract concepts concrete
- Multiple examples showing different aspects (entropy, energy, distribution)
- Challenges encourage experimentation and discovery
- Aligns with "See uncertainty, don't just calculate it" tenet

**Technical features:**
- Uses %matplotlib widget for interactive plots
- Combines billiards physics + visualisation modules
- Smooth animations (50 fps)
- Analysis plots (velocity histogram, energy evolution)
- Ready to run from fresh environment (just pip install -e .)

**Key demonstrations:**
- Hot/cold particle mixing (thermalization)
- Ordered → disordered evolution (entropy increase)
- Velocity distribution evolution (Maxwell-Boltzmann)
- Energy conservation verification
- All 4 initialization modes showcased

This is the **main deliverable** - ties together physics and visualisation into an interactive educational experience!

