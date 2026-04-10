# Kappenball

Kappenball is an interactive simulation game for exploring **uncertainty, decision-making, and the mathematics of procrastination**. A ball falls through a noisy channel and you must steer it to land in the right place — but the uncertainty in its trajectory means acting too early can be worse than waiting.

Created by [Neil Lawrence](https://inverseprobability.com) around 2012 as a playful way to build intuition for probabilistic inference and optimal control under uncertainty. It is discussed in his book [The Atomic Human](https://atomichuman.com).

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

To reset to the clean demo starting point at any time:

```bash
git checkout main
git reset --hard demo-start
git clean -fdx
```
