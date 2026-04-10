# Kappenball × VibeSafe: Live Demo Script

**Purpose**: Demonstrate the VibeSafe WHY → WHAT → HOW → DO → BUILD process.  
**Starting state**: Original 2012 MATLAB code only (`matlab/` + `README.md`).  
**End state**: A JavaScript web app, built live through structured AI-assisted development.  
**Estimated time**: 45–60 minutes.

---

## Resetting for a Fresh Demo

This script lives on the `master` branch, tagged `demo-start`. After a demo run, VibeSafe will have created directories (`tenets/`, `cip/`, `backlog/`, `requirements/`, `scripts/`, `whats-next`, etc.) and you may have committed files. Run these three commands to wipe everything and return to the clean starting state:

```bash
git checkout master
git reset --hard demo-start
git clean -fd
```

- `reset --hard` rewinds `master` to the `demo-start` tag, undoing any commits made during the demo
- `clean -fd` removes all untracked files and directories (the VibeSafe scaffold, `.venv-vibesafe/`, etc.)

After running these, `ls` should show only: `README.md  DEMO.md  matlab/`

---

## Pre-Demo Checklist

- [ ] Repo at the pre-VibeSafe commit — only `README.md` and `matlab/` visible
- [ ] VibeSafe NOT yet installed (no `cip/`, `tenets/`, `backlog/`, `requirements/`)
- [ ] This script open in a separate window
- [ ] Cursor IDE open on the repo, terminal at `/Users/neil/lawrennd/kappenball`

```bash
ls   # Should show: README.md  matlab/
```

---

## The Story (Tell This First)

> "Around 2012 I wrote this MATLAB code to demonstrate a concept I'd been thinking about — a falling ball whose trajectory is uncertain, and where you have to control it to land in the right place. The idea connects to decision-making under uncertainty and the mathematics of procrastination — when is it optimal to wait rather than act? I ended up calling it Kappenball.
>
> The problem is it requires MATLAB. It's 2026 and nobody teaches with MATLAB anymore. I want to modernise it — but more importantly, I want to show you how to use AI assistance *with structure*, so you're not just vibe-coding into the void.
>
> That structure is called VibeSafe. Let's install it."

---

## STAGE 1: WHY — Install VibeSafe, Define Tenets

**Concept**: Before writing a line of new code, we articulate *why* this project exists. **Tenets** are the non-negotiable principles that guide every later decision.

### Install

```bash
bash <(curl -s https://raw.githubusercontent.com/lawrennd/vibesafe/main/scripts/install-minimal.sh)
ls   # Now shows: matlab/  cip/  tenets/  backlog/  requirements/  scripts/  whats-next
```

> "VibeSafe added scaffolding for structured development. It didn't touch the MATLAB code."

### Run What's Next

```bash
./whats-next
```

> "Empty project. No tenets, no CIPs, no requirements. The system is telling us: start with WHY."

### Create the Three Tenets

**Tell the audience**: "Tenets are *not* goals. They're principles — the things you refuse to compromise on."

For each tenet, use this Cursor prompt (adapt the tenet name/description each time):

```
Create a tenet file for the Kappenball project at 
tenets/kappenball/uncertainty-accessible.md

The tenet is: Kappenball exists to make uncertainty concepts tangible and 
understandable. Complex probabilistic and entropy concepts should be 
demonstrated through simple, visual, interactive simulations that reveal 
the underlying principles without requiring deep mathematical background.

Quote: "See uncertainty, don't just calculate it."

Include examples from the MATLAB code (billiards showing entropy increase, 
falling ball showing probabilistic outcomes), counter-examples, and 
conflict resolutions. Use the VibeSafe tenet template format.
```

Repeat for:
- **Interactive Engagement** — *"Play with uncertainty, learn through exploration."* Learning should feel like discovery, not study. Keyboard controls, real-time feedback, experiment with parameters.
- **Easy to Modify** — *"From running to modifying in minutes, not hours."* Minimal dependencies, readable code, students can go from running to tinkering quickly.

### Check in with What's Next

```bash
./whats-next
```

> "Three tenets now registered. These are automatically in Cursor's context — every AI prompt from now on is guided by these principles without us repeating ourselves."

---

## STAGE 2: WHAT — Define Requirements

**Concept**: Requirements answer "what must be true?" — outcomes, not implementations. They flow *from* the tenets. The same requirement could be satisfied by many different HOWs.

### The Pivot Moment

> "The obvious move would be Python + Jupyter. But look at these tenets: 'run in minutes', 'interactive', 'accessible'. What's the most accessible platform in 2026?
>
> A web browser. Every student has one. Zero installation. A JavaScript app means you could put a link in a lecture slide and students play it on their phones. Let's write that as a requirement."

### Create Requirements

Use Cursor prompts like:

```
Based on the Kappenball tenets, create a requirements file at
requirements/req0001_no-install-required.md

The requirement: Students and educators must be able to run Kappenball 
by visiting a URL — no software installation, no accounts, no downloads.

Link it to the easy-to-modify and uncertainty-accessible tenets.
Include clear acceptance criteria (runs in any browser, works on mobile, 
embeddable in course pages).
```

Create four requirements this way:
1. **REQ-0001**: Zero-install — runs in a browser, no setup
2. **REQ-0002**: Real-time interactive controls — keyboard/touch, <100ms response
3. **REQ-0003**: Uncertainty visually demonstrable — see probability clouds, entropy in real-time
4. **REQ-0004**: Shareable and embeddable — permanent URL, iframe-embeddable, GitHub Pages

### Check in with What's Next

```bash
./whats-next
```

> "Four requirements. They describe *outcomes*. A Python notebook could satisfy REQ-0002. Only a web app can satisfy REQ-0001 and REQ-0004. Requirements give you the freedom to choose the right HOW."

---

## STAGE 3: HOW — Write a CIP

**Concept**: CIPs (Code Improvement Plans) say *how* we'll achieve the requirements. Written before code — a contract between intent and implementation.

> "Single CIP: build a JavaScript web app, HTML5 Canvas, no npm, no build step. A single HTML file students can open locally or share as a URL."

### Create CIP-0001

```bash
cp cip/cip_template.md cip/cip0001.md
```

Cursor prompt:

```
Fill out cip/cip0001.md: "Build Kappenball as a JavaScript Web App"

This CIP proposes building Kappenball as a self-contained HTML5 Canvas 
JavaScript app — no framework, no build step. A single kappenball.html 
that runs in any browser.

It addresses all four requirements (REQ-0001 through REQ-0004) and 
aligns with all three tenets.

Architecture: single HTML file (Phase 1), optionally refactored into 
web/js/ modules (Phase 2). Core components: physics engine ported from 
MATLAB, Canvas 2D renderer, requestAnimationFrame game loop, keyboard 
and touch controls.

Implementation order: falling ball physics → Canvas renderer → controls 
→ billiards physics → entropy visualisation → GitHub Pages deploy.

Reference the MATLAB source files: matlab/simulateFallingBall.m, 
matlab/demKappenBall.m, matlab/simulateBilliards.m, 
matlab/entropyBilliards.m
```

### Check in with What's Next

```bash
./whats-next
```

> "CIP is in context. When I ask the AI to implement it, it has a complete specification — not just 'write a game', but why, what outcomes it must achieve, and in what order."

---

## STAGE 4: DO — Create Backlog Tasks

**Concept**: Backlog tasks are the individual chunks of work. Each has acceptance criteria — you know exactly when you're done.

### Create Tasks

Cursor prompt:

```
Based on CIP-0001, create backlog tasks in backlog/features/ for:

1. 2026-04-10_falling-ball-physics-js.md — port simulateFallingBall.m 
   to JavaScript. Acceptance criteria: gravity with configurable noise σ, 
   left/right controls, wall bouncing, landing detection, energy metric.

2. 2026-04-10_canvas-renderer-falling-ball.md — HTML5 Canvas renderer 
   at 60fps. Acceptance criteria: uncertainty halo around ball (radius ∝ σ), 
   fading trajectory trail, target buckets at bottom, score overlay.

3. 2026-04-10_billiards-physics-js.md — port simulateBilliards.m. 
   Acceptance criteria: ordered grid init, elastic collisions, Shannon 
   entropy computed per frame, visible order-to-disorder transition.

4. 2026-04-10_github-pages-deploy.md — GitHub Pages deployment. 
   Acceptance criteria: live at lawrennd.github.io/kappenball, 
   iframe-embeddable, works on mobile.

Link each task to CIP-0001.
```

### Check in with What's Next

```bash
./whats-next
```

> "Now look: WHY (tenets), WHAT (requirements), HOW (CIP), DO (tasks). Every piece visible. When I ask the AI to implement a task, it has full context all the way back to why this project exists."

---

## STAGE 5: BUILD — Implement with AI

**Concept**: Now we build. The structure means prompts are grounded — we're not asking "write me a game", we're asking "implement this task, per this CIP, guided by these tenets."

### Prompt 1: Falling Ball Physics

```
Looking at CIP-0001 and the backlog task for "Port falling ball physics 
to JavaScript", implement the falling ball physics.

Start with a single self-contained file kappenball.html. Port the physics 
from matlab/simulateFallingBall.m. Ball falls under gravity with Gaussian 
noise, responds to arrow key controls, bounces off walls.

Keep the code readable and commented — this is educational software.
```

**While it generates**: 

> "I didn't say 'make it readable' or 'make it educational'. That's in the tenets. The AI already knows. I only had to specify the task — which was already written in the backlog."

### Prompt 2: Canvas Renderer

```
Add the Canvas renderer per the backlog task "Build Canvas renderer for 
falling ball". The ball should have an uncertainty halo whose radius 
reflects the noise level σ. Trajectory leaves a fading trail. Target 
buckets shown at the bottom.
```

### Prompt 3: Billiards (if time permits)

```
Add the entropy billiards simulation from matlab/simulateBilliards.m 
and matlab/entropyBilliards.m. N balls start in an ordered grid, then 
scatter. Display Shannon entropy in real-time. The visual transition 
from order to disorder is the whole point — make it clear.
```

### Test It

Open `kappenball.html` in a browser. Play it.

> "No npm install. No build step. Open in browser. Done."

Check against tenets:
- **Uncertainty Accessible**: Can you *see* the uncertainty?
- **Interactive Engagement**: Is it fun to play?
- **Easy to Modify**: Can you open the file and understand it?

### Commit

```bash
git add kappenball.html
git commit -m "CIP-0001: Implement falling ball physics and Canvas renderer"
```

---

## STAGE 6: REFLECT

> "What did we do?
>
> 1. **WHY**: Wrote down why this project exists — as tenets, in the repo, in Cursor's context
> 2. **WHAT**: Derived requirements from those tenets — outcomes, not implementations
> 3. **HOW**: Designed a specific approach before touching code
> 4. **DO**: Broke the design into tasks with clear acceptance criteria
> 5. **BUILD**: Used AI to implement — with all that context already loaded
>
> You can hand this repo to a student and they'll understand not just *what* the code does, but *why* every decision was made. The structure is the documentation.
>
> That's VibeSafe."

---

## Appendix: Troubleshooting

**VibeSafe install fails**
```bash
bash <(curl -s https://raw.githubusercontent.com/lawrennd/vibesafe/main/scripts/install-minimal.sh)
```

**`./whats-next` fails**
```bash
source .venv-vibesafe/bin/activate && python scripts/whats_next.py
```

**AI generates wrong physics**
> Teaching moment: open `matlab/simulateFallingBall.m`, read it with the audience, then re-prompt with the MATLAB source pasted in. "The reference implementation is right here — let's use it."

**Demo runs long** — cut billiards, skip GitHub Pages, pre-create backlog tasks  
**Demo runs short** — linger on `./whats-next` output at each stage, walk through the MATLAB source, ask audience to suggest requirements
