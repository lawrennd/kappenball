"""
Billiards physics simulation module.

Provides functions for simulating 2D billiards with elastic collisions,
demonstrating entropy increase and uncertainty concepts.

Ported from MATLAB (2012) to Python (2026).

Functions:
    initialise: Initialise billiards simulation state
    simulate_step: Run single billiards simulation step
    compute_vectors: Compute particle velocities and trajectories
"""

import numpy as np
from typing import Optional, Literal


def initialise(
    n_particles: int = 10,
    init_type: Literal['rand', 'randn', 'uniform', 'hotCold'] = 'rand',
    box_xlim: tuple[float, float] = (0.0, 10.0),
    box_ylim: tuple[float, float] = (0.0, 10.0),
    random_state: Optional[int] = None
) -> dict:
    """
    Initialise billiards simulation with particles.
    
    Ported from initializeBilliards.m
    
    Args:
        n_particles: Number of particles to simulate
        init_type: Initialisation mode:
            - 'rand': Random positions/velocities, varied sizes
            - 'randn': Random normal distribution (physics-based)
            - 'uniform': Uniform vertical spacing, opposite velocities
            - 'hotCold': Two temperature regions (hot/cold)
        box_xlim: (xmin, xmax) for simulation box
        box_ylim: (ymin, ymax) for simulation box
        random_state: Random seed for reproducibility
        
    Returns:
        dict: Simulation state containing:
            - X: (n, 2) array of positions
            - V: (n, 2) array of velocities  
            - r: (n,) array of radii
            - colors: (n, 3) array of RGB colours
            - mass: (n,) array of particle masses
            - box_xlim, box_ylim: Box boundaries
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    # Physical constants for physics-based initialisation
    T = 1e-2  # Temperature in K
    BOLTZ_K = 1.3806503e-23  # Boltzmann constant
    ATOMIC_MASS = 1.66053886e-27  # 1 AMU in kg
    
    if init_type == 'hotCold':
        # Two temperature regions: hot (fast) and cold (slow)
        mass = 2 * ATOMIC_MASS * np.ones(n_particles)
        var_v = BOLTZ_K * T / mass
        v_std = np.sqrt(var_v)
        
        # Hot particles (first half) move faster
        v_x = np.concatenate([
            2 * np.ones(int(np.ceil(n_particles / 2))),
            0.5 * np.ones(int(np.floor(n_particles / 2)))
        ])
        v_y = np.random.randn(n_particles) * 0.004
        V = np.column_stack([v_x, v_y]) * v_std[:, np.newaxis]
        
        # Uniform vertical spacing
        spacing = (box_ylim[1] - box_ylim[0]) / n_particles
        r = spacing * 0.475 * np.ones(n_particles)
        x_center = (box_xlim[1] - box_xlim[0]) * 0.5 + box_xlim[0]
        X = np.column_stack([
            x_center * np.ones(n_particles),
            np.arange(box_ylim[0] + 0.5 * spacing, box_ylim[1], spacing)
        ])
        
        # Hot particles yellow, cold particles cyan
        colors = np.zeros((n_particles, 3))
        colors[:int(np.ceil(n_particles / 2))] = [1, 1, 0]  # Yellow
        colors[int(np.ceil(n_particles / 2)):] = [0, 1, 1]  # Cyan
        
    elif init_type == 'uniform':
        # Uniform spacing with opposite velocities (demonstrates mixing)
        mass = 2 * ATOMIC_MASS * np.ones(n_particles)
        var_v = BOLTZ_K * T / mass
        v_std = np.sqrt(var_v)
        
        # Alternating velocities
        v_x = np.concatenate([
            np.ones(int(np.ceil(n_particles / 2))),
            -np.ones(int(np.floor(n_particles / 2)))
        ])
        v_y = np.random.randn(n_particles) * 0.001
        V = np.column_stack([v_x, v_y]) * v_std[:, np.newaxis]
        
        # Uniform vertical spacing
        spacing = (box_ylim[1] - box_ylim[0]) / n_particles
        r = spacing * 0.475 * np.ones(n_particles)
        x_center = (box_xlim[1] - box_xlim[0]) * 0.5 + box_xlim[0]
        X = np.column_stack([
            x_center * np.ones(n_particles),
            np.arange(box_ylim[0] + 0.5 * spacing, box_ylim[1], spacing)
        ])
        
        colors = np.tile([1, 1, 0], (n_particles, 1))  # All yellow
        
    elif init_type == 'randn':
        # Random normal distribution (physics-based, equal-sized)
        r = 0.0625 * np.ones(n_particles)
        mass = 2 * ATOMIC_MASS * np.ones(n_particles)
        var_v = BOLTZ_K * T / mass
        v_std = np.sqrt(var_v)
        
        V = np.random.randn(n_particles, 2) * v_std[:, np.newaxis]
        
        # Generate non-overlapping random positions
        X = _generate_non_overlapping_positions(
            n_particles, r, box_xlim, box_ylim
        )
        
        colors = np.tile([1, 1, 0], (n_particles, 1))  # All yellow
        
    else:  # 'rand' (default)
        # Random sizes, positions, and velocities
        r = 0.125 + 0.125 * np.random.rand(n_particles)
        V = 10 * (-1 + 2 * np.random.rand(n_particles, 2))
        
        # Generate non-overlapping random positions
        X = _generate_non_overlapping_positions(
            n_particles, r, box_xlim, box_ylim
        )
        
        # Colourful balls
        ball_colours = np.array([
            [1, 0, 0], [1, 0, 0.5], [1, 0.5, 0], [0, 1, 0], [0, 0, 1],
            [1, 1, 0], [1, 1, 1], [0, 0.3, 0], [0, 0, 0], [0.65, 0.65, 0.65],
            [0, 0.75, 0.75], [0.3, 0, 0.6], [0.95, 0.65, 0.75], [0.5, 0.25, 0],
            [0, 0.2, 0.4], [0.9, 0.4, 0.7], [0.4, 0.2, 0.3], [0.65, 0.55, 0.15],
            [0.25, 0.35, 0.25], [0.5, 0, 0]
        ])
        colors = ball_colours[np.arange(n_particles) % len(ball_colours)]
        mass = np.pi * r**2  # Mass proportional to area
    
    return {
        'X': X,
        'V': V,
        'r': r,
        'colors': colors,
        'mass': mass if 'mass' in locals() else np.pi * r**2,
        'box_xlim': box_xlim,
        'box_ylim': box_ylim,
        'n_particles': n_particles
    }


def _generate_non_overlapping_positions(
    n_particles: int,
    r: np.ndarray,
    box_xlim: tuple[float, float],
    box_ylim: tuple[float, float],
    max_attempts: int = 1000
) -> np.ndarray:
    """
    Generate random non-overlapping particle positions.
    
    Args:
        n_particles: Number of particles
        r: Array of particle radii
        box_xlim: (xmin, xmax) boundaries
        box_ylim: (ymin, ymax) boundaries
        max_attempts: Maximum attempts to find valid configuration
        
    Returns:
        (n, 2) array of particle positions
    """
    r_max = np.max(r)
    
    for attempt in range(max_attempts):
        # Random positions within box (accounting for radii)
        X = np.column_stack([
            (box_xlim[1] - box_xlim[0] - 2 * r_max) * np.random.rand(n_particles) + box_xlim[0] + r_max,
            (box_ylim[1] - box_ylim[0] - 2 * r_max) * np.random.rand(n_particles) + box_ylim[0] + r_max
        ])
        
        # Check for overlaps
        if not _check_overlaps(X, r):
            return X
    
    # If we couldn't find non-overlapping positions, return last attempt
    print(f"Warning: Could not generate non-overlapping positions after {max_attempts} attempts")
    return X


def _check_overlaps(X: np.ndarray, r: np.ndarray) -> bool:
    """Check if any particles overlap."""
    n = len(X)
    for i in range(n):
        for j in range(i + 1, n):
            dist = np.linalg.norm(X[i] - X[j])
            if dist < (r[i] + r[j]):
                return True
    return False


def simulate_step(state: dict, dt: float = 2e-3) -> dict:
    """
    Simulate one step of billiards physics.
    
    Ported from simulateBilliards.m
    
    Args:
        state: Simulation state from initialise()
        dt: Time step for simulation
        
    Returns:
        dict: Updated simulation state
    """
    X = state['X'].copy()
    V = state['V'].copy()
    r = state['r']
    mass = state['mass']
    box_xlim = state['box_xlim']
    box_ylim = state['box_ylim']
    n_particles = state['n_particles']
    
    # Wall collision detection and resolution
    # Positive edge (right/top walls)
    d_pos = X + r[:, np.newaxis] - np.array([[box_xlim[1], box_ylim[1]]])
    collision_pos = d_pos >= 0
    dt_pos = np.where(collision_pos & (V != 0), d_pos / V, 0)
    X = X - V * dt_pos
    V = V * (2 * (~collision_pos).astype(float) - 1)
    
    # Negative edge (left/bottom walls)
    d_neg = X - r[:, np.newaxis] - np.array([[box_xlim[0], box_ylim[0]]])
    collision_neg = d_neg <= 0
    dt_neg = np.where(collision_neg & (V != 0), d_neg / V, 0)
    X = X - V * dt_neg
    V = V * (2 * (~collision_neg).astype(float) - 1)
    
    # Particle-particle collision detection
    collisions = _detect_collisions(X, r)
    
    if len(collisions) > 0:
        # Resolve each collision
        for i, j in collisions:
            # Normal direction between particles
            norm_dist = X[i] - X[j]
            norm_dist = norm_dist / np.linalg.norm(norm_dist)
            
            # Velocity components along collision normal
            va_i = np.dot(V[i], norm_dist)
            va_j = np.dot(V[j], norm_dist)
            
            # Time to collision
            dist = np.linalg.norm(X[i] - X[j])
            dt_collision = abs(r[i] + r[j] - dist) / (abs(va_i) + abs(va_j) + 1e-10)
            
            # Transformation matrix (rotate to collision frame)
            M = np.array([[norm_dist[0], norm_dist[1]],
                         [-norm_dist[1], norm_dist[0]]])
            
            # Transform velocities to collision frame
            v_old = np.concatenate([V[i], V[j]])
            v_new = np.zeros(4)
            v_new[:2] = M.T @ V[i]
            v_new[2:] = M.T @ V[j]
            
            # Elastic collision (conservation of momentum and energy)
            f = 1 + mass[i] / mass[j]
            g = 1 + mass[j] / mass[i]
            collision_matrix = np.array([
                [1 - 2/f, 0, 2/f, 0],
                [0, 1, 0, 0],
                [2/g, 0, 1 - 2/g, 0],
                [0, 0, 0, 1]
            ])
            v_new_col = collision_matrix @ v_new
            
            # Transform back to original frame
            V[i] = M @ v_new_col[:2]
            V[j] = M @ v_new_col[2:]
            
            # Update positions after collision
            X[i] = X[i] + V[i] * dt_collision
            X[j] = X[j] + V[j] * dt_collision
    
    # Propagate positions
    X = X + V * dt
    
    # Update state
    state['X'] = X
    state['V'] = V
    
    return state


def _detect_collisions(X: np.ndarray, r: np.ndarray) -> list[tuple[int, int]]:
    """
    Detect particle-particle collisions.
    
    Args:
        X: (n, 2) array of positions
        r: (n,) array of radii
        
    Returns:
        List of (i, j) collision pairs
    """
    n = len(X)
    collisions = []
    
    for i in range(n):
        for j in range(i + 1, n):
            dist = np.linalg.norm(X[i] - X[j])
            if dist < (r[i] + r[j]):
                collisions.append((i, j))
    
    return collisions


def compute_vectors(state: dict) -> dict:
    """
    Compute velocity vectors and trajectories for visualisation.
    
    Args:
        state: Simulation state dict
        
    Returns:
        dict: Vector data for visualisation containing:
            - arrow_starts: (n, 2) start points for velocity arrows
            - arrow_vectors: (n, 2) velocity direction vectors
    """
    X = state['X']
    V = state['V']
    r = state['r']
    
    # Velocity arrows start at particle edge in velocity direction
    vel_norm = V / (np.linalg.norm(V, axis=1, keepdims=True) + 1e-10)
    arrow_starts = X + r[:, np.newaxis] * vel_norm
    
    # Arrow length proportional to velocity magnitude
    vel_factor = 0.025
    arrow_vectors = V * vel_factor
    
    return {
        'arrow_starts': arrow_starts,
        'arrow_vectors': arrow_vectors
    }
