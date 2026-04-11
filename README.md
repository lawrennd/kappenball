# Kappenball

Kappenball is an interactive simulation game for exploring **uncertainty, decision-making, and the mathematics of procrastination**. A ball falls through a noisy channel and you must steer it to land in the right place — but the uncertainty in its trajectory means acting too early can be worse than waiting.

Created by [Neil Lawrence](https://inverseprobability.com) around 2012 as part of his inaugural lecture. It's a way of showing how uncertainty influences decision making.

## What's in this repo

The original 2012 MATLAB code lives in `matlab/`:

| File | Description |
|------|-------------|
| `demKappenBall.m` | Interactive falling ball demo |
| `simulateFallingBall.m` | Ball physics with Gaussian noise |
| `demEntropyBilliards.m` | Entropy billiards demo |
| `simulateBilliards.m` | Billiards physics engine |
| `entropyBilliards.m` | Shannon entropy computation |

## Running the MATLAB demos

Open MATLAB (or GNU Octave — `brew install octave`) and run:

```matlab
cd matlab
demKappenBall          % falling ball game
demEntropyBilliards    % entropy billiards
```

## Demo: Modernising with VibeSafe

This repo is also used to demonstrate [VibeSafe](https://github.com/lawrennd/vibesafe) — a structured approach to AI-assisted development using the **WHY → WHAT → HOW → DO** process. See `DEMO.md` for a step-by-step walkthrough of installing VibeSafe and building a JavaScript web app version of Kappenball from scratch, live.


## Running a Demo

Before starting, create a dated branch so that `main` stays clean throughout:

```bash
git checkout main
git pull
git checkout -b demo/$(date +%Y-%m-%d)
```

All work during the demo (VibeSafe scaffolding, CIPs, backlog tasks, code, etc.) is committed to this branch. See `DEMO.md` for the full step-by-step script.

## Resetting for a Fresh Demo

### If you worked on a demo branch (recommended)

Simply switch back to `main` and delete the branch:

```bash
git checkout main
git branch -D demo/YYYY-MM-DD   # replace with the actual branch name
git clean -fd                   # remove any leftover untracked files
```

### If changes landed on main

`main` is tagged at `demo-start`. Reset hard to that tag:

```bash
git checkout main
git reset --hard demo-start
git clean -fd
```

- `reset --hard` rewinds `main` to the `demo-start` tag, undoing any commits made during the demo
- `clean -fd` removes untracked files and directories (`.cursor/`, `.venv-vibesafe/`, `tenets/`, `cip/`, etc.)

Files listed in `.gitignore` — such as `*.code-workspace` and IDE settings — are preserved across resets.

After either reset, `ls` should show only: `DEMO.md  README.md  matlab/`

