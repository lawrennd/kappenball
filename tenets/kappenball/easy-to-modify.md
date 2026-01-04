---
id: "easy-to-modify"
title: "Easy to Modify"
status: "Active"
created: "2026-01-04"
last_reviewed: "2026-01-04"
review_frequency: "Annual"
conflicts_with: []
tags:
- tenet
- accessibility
- maintainability
---

# Tenet: Easy to Modify

## Tenet

**Description**: Students and educators should be able to run and modify Kappenball with minimal friction. The code should be readable, well-structured, and documented so that users can understand how it works and adapt it for their own purposes. Dependencies should be minimal, installation should be straightforward, and the barrier to experimentation should be as low as possible. The project should encourage tinkering and extension.

**Quote**: *"From running to modifying in minutes, not hours."*

**Examples**:
- MATLAB code with clear function names that describe their purpose (`simulateBilliards`, `entropyBilliards`)
- Self-contained simulation files that can be run independently
- Minimal external dependencies (using standard MATLAB functions)
- Clear separation between physics simulation, visualization, and interaction logic
- Comments and documentation that explain both the "what" and the "why"
- Example scripts (`demKappenBall.m`, `demEntropyBilliards.m`) that show how to use the components

**Counter-examples**:
- Requiring installation of obscure toolboxes or libraries
- Monolithic code files mixing physics, visualization, and UI in hard-to-separate ways
- Variable names like `x1`, `tmp`, `data` without context
- Complex build processes or compilation steps before running
- Assumptions about system configuration without clear documentation
- Tightly coupled code that breaks when users try to modify one part

**Conflicts**:
- **Performance Optimization**: Highly optimized code might sacrifice readability
- Resolution: Prioritize clarity over micro-optimizations; document any necessary performance tricks
- **Feature Richness**: More features can make the codebase harder to understand
- Resolution: Keep core demonstrations simple; provide advanced features as optional modules

