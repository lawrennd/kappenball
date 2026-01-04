---
id: "0003"
title: "Uncertainty concepts must be visually demonstrable"
status: "Proposed"
priority: "High"
created: "2026-01-04"
last_updated: "2026-01-04"
related_tenets: ["uncertainty-accessible"]
stakeholders: ["students", "educators"]
tags: ["visualization", "pedagogy"]
---

# REQ-0003: Uncertainty concepts must be visually demonstrable

> **Remember**: Requirements describe **WHAT** should be true (outcomes), not HOW to achieve it.

## Description

Complex probabilistic and entropy concepts must be made tangible through clear, intuitive visualizations. Students should be able to "see uncertainty" develop in real-time through visual representations of particle trajectories, velocity vectors, state evolution, and entropy changes. The visualizations should make abstract mathematical concepts concrete without requiring deep mathematical background.

**Why this matters**: This is the core of the **Uncertainty Accessible** tenet - "See uncertainty, don't just calculate it." Visual demonstrations bridge the gap between intuition and theory, making abstract probabilistic concepts accessible to students at all levels.

**Who benefits**: 
- Students learning probability and statistical mechanics concepts
- Educators demonstrating entropy increase, probabilistic behavior, and uncertainty quantification
- Visual learners who understand concepts better through observation than equations

## Acceptance Criteria

What does "done" look like? Be specific about outcomes, not implementation:

- [ ] Entropy increase is visually observable as ordered states evolve to disordered states
- [ ] Probabilistic outcomes in the falling ball demo are clearly visible
- [ ] Velocity vectors and particle trajectories can be toggled on/off for clarity
- [ ] Visualizations remain clear and interpretable even with many particles
- [ ] Color schemes and visual design emphasize the uncertainty concepts being taught
- [ ] Animations are smooth enough to observe continuous evolution (target: 30+ fps)

## Notes (Optional)

The original billiards demo excels at showing entropy increase visually - balls starting in an ordered configuration scatter into disorder. The falling ball demo shows probabilistic outcomes accumulating. These visual strengths must be preserved and potentially enhanced in the modernization.

## References

- **Related Tenets**: uncertainty-accessible
- **Related CIP**: CIP-0001 (proposes matplotlib for visualization, with animation support)
- **Educational Context**: These demos were created for a 2012 lecture on uncertainty

## Progress Updates

### 2026-01-04
Requirement created. CIP-0001 proposes matplotlib with animation capabilities, which should satisfy these visualization needs.

