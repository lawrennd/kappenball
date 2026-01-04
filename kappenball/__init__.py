"""
Kappenball: Interactive Physics Simulations for Uncertainty Education

A Python package for demonstrating uncertainty concepts through interactive
billiards and falling ball simulations. Originally created in MATLAB (2012),
modernized for Python + Jupyter (2026).

Modules:
    billiards: 2D/3D billiards physics simulation
    falling_ball: Falling ball with probabilistic outcomes
    physics: Core physics utilities
    visualisation: Matplotlib-based visualisation tools
"""

__version__ = "0.1.0"
__author__ = "Neil Lawrence"

from . import billiards
from . import falling_ball
from . import physics
from . import visualisation

__all__ = ["billiards", "falling_ball", "physics", "visualisation"]

