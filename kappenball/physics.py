"""
Core physics utilities module.

Provides common physics functions used across simulations, including
collision detection, energy calculations, and helper utilities.

Functions:
    detect_collision: Check for particle-particle or particle-wall collisions
    resolve_collision: Compute post-collision velocities (elastic)
    compute_energy: Calculate kinetic energy of system
    compute_entropy: Calculate entropy measure for particle distribution
"""

import numpy as np


def detect_collision(positions, radii, box_size):
    """
    Detect collisions between particles and with walls.
    
    Args:
        positions: (N, D) array of particle positions (N particles, D dimensions)
        radii: Particle radii (scalar or array of length N)
        box_size: Size of bounding box
        
    Returns:
        list: Collision events [(particle_i, particle_j or 'wall', ...)]
    """
    # To be implemented
    raise NotImplementedError("Collision detection logic")


def resolve_collision(velocities, positions, collision_events, masses=None):
    """
    Resolve elastic collisions, updating velocities.
    
    Args:
        velocities: (N, D) array of particle velocities
        positions: (N, D) array of particle positions
        collision_events: List of collision events from detect_collision()
        masses: Particle masses (defaults to equal mass)
        
    Returns:
        array: Updated velocities after collision resolution
    """
    # To be implemented: Elastic collision physics
    raise NotImplementedError("Elastic collision resolution")


def compute_energy(velocities, masses=None):
    """
    Calculate total kinetic energy of the system.
    
    Args:
        velocities: (N, D) array of particle velocities
        masses: Particle masses (defaults to equal mass)
        
    Returns:
        float: Total kinetic energy
    """
    if masses is None:
        masses = np.ones(velocities.shape[0])
    return 0.5 * np.sum(masses[:, np.newaxis] * velocities**2)


def compute_entropy(positions, box_size, n_bins=10):
    """
    Calculate entropy measure based on spatial distribution.
    
    Args:
        positions: (N, D) array of particle positions
        box_size: Size of bounding box
        n_bins: Number of bins per dimension for histogram
        
    Returns:
        float: Entropy measure (higher = more disordered)
    """
    # To be implemented: Entropy calculation (e.g., from spatial histogram)
    raise NotImplementedError("Port from MATLAB: entropyBilliards.m")

