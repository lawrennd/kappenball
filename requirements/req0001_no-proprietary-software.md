---
id: "0001"
title: "Simulations runnable without proprietary software"
status: "Proposed"
priority: "High"
created: "2026-01-04"
last_updated: "2026-01-04"
related_tenets: ["easy-to-modify"]
stakeholders: ["students", "educators", "researchers"]
tags: ["accessibility", "licensing"]
---

# REQ-0001: Simulations runnable without proprietary software

> **Remember**: Requirements describe **WHAT** should be true (outcomes), not HOW to achieve it.

## Description

Students and educators must be able to run Kappenball simulations without purchasing proprietary software licenses. The barrier to entry should be eliminated by using freely available, open-source tools that can be installed on common platforms (macOS, Linux, Windows).

**Why this matters**: This directly supports the **Easy to Modify** tenet - proprietary software creates a financial and accessibility barrier that prevents students from going "from running to modifying in minutes, not hours."

**Who benefits**: 
- Students who cannot afford MATLAB licenses
- Educators at institutions with limited software budgets
- Self-learners and researchers worldwide

## Acceptance Criteria

What does "done" look like? Be specific about outcomes, not implementation:

- [ ] Simulations can be executed using only free, open-source software
- [ ] Installation instructions work on macOS, Linux, and Windows
- [ ] No proprietary licenses are required to run or modify the code
- [ ] Installation can be completed in under 15 minutes on a typical system

## Notes (Optional)

The original 2012 MATLAB implementation created a significant barrier. Modern alternatives (Python, JavaScript) provide equivalent or superior capabilities for educational physics simulations without licensing costs.

## References

- **Related Tenets**: easy-to-modify
- **Related CIP**: CIP-0001 (Python + Jupyter migration)

## Progress Updates

### 2026-01-04
Requirement created. CIP-0001 proposes Python + Jupyter as the solution, which fully satisfies this requirement.

