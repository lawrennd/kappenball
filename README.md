# Kappenball

Interactive physics simulations for uncertainty education. Originally created in MATLAB (2012), modernized for Python + Jupyter (2026).

## Project Tenets

Kappenball is guided by three core principles:

1. **[Uncertainty Accessible](tenets/kappenball/uncertainty-accessible.md)** - *"See uncertainty, don't just calculate it."*
   
   Make complex probabilistic and entropy concepts tangible through visual, interactive simulations.

2. **[Interactive Engagement](tenets/kappenball/interactive-engagement.md)** - *"Play with uncertainty, learn through exploration."*
   
   Learning should feel like discovery. Fun, responsive controls invite experimentation and deeper understanding.

3. **[Easy to Modify](tenets/kappenball/easy-to-modify.md)** - *"From running to modifying in minutes, not hours."*
   
   Students and educators should be able to run and adapt the code with minimal friction.

These tenets guide all decisions about features, code structure, and documentation.

## Installation

**Requirements:** Python 3.8 or higher

### Quick Start

```bash
# Clone the repository
git clone https://github.com/lawrennd/kappenball.git
cd kappenball

# Install in development mode
pip install -e .

# Optional: Install development dependencies
pip install -r requirements-dev.txt
```

### Using with Jupyter

```bash
# Install with dependencies
pip install -e .

# Launch Jupyter
jupyter notebook

# Open notebooks in the notebooks/ directory
```

## Project Structure

```
kappenball/
├── kappenball/           # Python package
│   ├── billiards.py      # Billiards physics simulation
│   ├── falling_ball.py   # Falling ball simulation
│   ├── physics.py        # Core physics utilities
│   └── visualization.py  # Matplotlib visualization
├── notebooks/            # Jupyter notebooks (demos)
├── matlab/              # Original MATLAB code (2012, reference)
├── tests/               # Unit tests
├── requirements.txt     # Core dependencies
└── README.md           # This file
```

## Status

🚧 **Migration in progress** from MATLAB to Python + Jupyter (see [CIP-0001](cip/cip0001.md))

- ✅ Package structure created
- 🔄 Physics porting in progress
- 🔄 Visualization layer in development
- 🔄 Interactive notebooks coming soon

## Development

```bash
# Run tests
pytest

# Format code
black kappenball/ tests/

# Lint code
ruff check kappenball/ tests/
```

## Original MATLAB Code

The original MATLAB code (2012) is preserved in `matlab/` directory:
- `demEntropyBilliards.m` - Entropy demonstration with billiards
- `demKappenBall.m` - Falling ball demonstration
- Run with GNU Octave: `brew install octave`

## License

MIT License (see matlab/license.txt for original code)
