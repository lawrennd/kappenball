"""
Falling ball physics simulation module (Kappenball game).

Provides functions for simulating a falling ball with probabilistic outcomes
and keyboard control, demonstrating uncertainty and stochastic processes.

Ported from MATLAB (2012) to Python (2026).

The Kappenball game: Guide a falling ball through holes in a pin array
using keyboard controls while managing uncertainty in the motion.

Functions:
    initialise: Initialise falling ball simulation state
    simulate_step: Run falling ball simulation step
    apply_control: Apply keyboard/control input to ball dynamics
    check_collision: Check if ball hit pins or went through hole
"""

import numpy as np
from typing import Optional, Literal


def initialise(
    initial_position: Optional[tuple[float, float]] = None,
    initial_velocity: Optional[tuple[float, float]] = None,
    box_xlim: tuple[float, float] = (-12.0, 12.0),
    box_ylim: tuple[float, float] = (-2.0, 10.0),
    hole_center: float = 6.0,
    hole_width: float = 4.0,
    pin_height: float = 0.5,
    velocity_variance: float = 0.0,
    random_state: Optional[int] = None
) -> dict:
    """
    Initialise falling ball simulation (Kappenball game).
    
    Ported from demKappenBall.m
    
    Args:
        initial_position: Initial (x, y) position, defaults to (0, 10)
        initial_velocity: Initial (vx, vy) velocity, defaults to (0, -1)
        box_xlim: (xmin, xmax) for simulation box
        box_ylim: (ymin, ymax) for simulation box
        hole_center: X-coordinate of hole center
        hole_width: Width of the hole
        pin_height: Height of pin array above floor
        velocity_variance: Variance of random horizontal velocity noise
        random_state: Random seed for reproducibility
        
    Returns:
        dict: Simulation state containing:
            - x: Current position [x, y]
            - v: Current velocity [vx, vy]
            - r: Ball radius
            - box_xlim, box_ylim: Box boundaries
            - hole_center, hole_width: Hole geometry
            - pin_height: Pin array height
            - v_var: Velocity variance (uncertainty level)
            - score: Number of successful drops
            - energy_count: Number of control inputs used
            - gravity: Gravitational acceleration
    """
    if initial_position is None:
        initial_position = (0.0, 10.0)
    if initial_velocity is None:
        initial_velocity = (0.0, -1.0)
    
    # Create persistent random number generator
    rng = np.random.RandomState(random_state)
    
    return {
        'x': np.array(initial_position, dtype=float),
        'v': np.array(initial_velocity, dtype=float),
        'r': 0.25,  # Ball radius
        'box_xlim': box_xlim,
        'box_ylim': box_ylim,
        'hole_center': hole_center,
        'hole_width': hole_width,
        'pin_height': pin_height,
        'v_var': velocity_variance,
        'score': 0,
        'energy_count': 0,
        'gravity': -1.0,  # Downward acceleration
        'last_collision': None,  # 'success', 'bang', or None
        'rng': rng  # Persistent random number generator
    }


def simulate_step(
    state: dict,
    dt: float = 0.01,
    control_input: Optional[Literal['left', 'right']] = None
) -> dict:
    """
    Simulate one step of falling ball physics (Kappenball game).
    
    Ported from simulateFallingBall.m
    
    Args:
        state: Simulation state from initialise()
        dt: Time step for simulation
        control_input: Keyboard control ('left', 'right', or None)
        
    Returns:
        dict: Updated simulation state with 'last_collision' field indicating
              what happened ('success' if through hole, 'bang' if hit pins, None otherwise)
    """
    state = state.copy()
    
    x = state['x'].copy()
    v = state['v'].copy()
    r = state['r']
    box_xlim = state['box_xlim']
    box_ylim = state['box_ylim']
    hole_center = state['hole_center']
    hole_width = state['hole_width']
    pin_height = state['pin_height']
    v_var = state['v_var']
    
    # Store old position for collision recovery
    old_x = x.copy()
    
    # Update position
    x = x + v * dt
    
    # Wall collision (left/right)
    if x[0] + r > box_xlim[1] or x[0] - r < box_xlim[0]:
        v[0] = -v[0]
        x = old_x  # Restore position
    
    # Add random horizontal velocity noise (uncertainty)
    v[0] = state['rng'].randn() * v_var + v[0] * 0.95
    
    # Apply control input (keyboard)
    if control_input == 'left':
        v[0] -= 0.5
        state['energy_count'] += 1
    elif control_input == 'right':
        v[0] += 0.5
        state['energy_count'] += 1
    
    # Reset collision status
    state['last_collision'] = None
    
    # Check for collision with pins or success through hole
    # Pin array has TWO holes at -hole_center and +hole_center
    # Pin regions: [xlim[0], -cent-w/2], [-cent+w/2, cent-w/2], [cent+w/2, xlim[1]]
    # Hole regions: [-cent-w/2, -cent+w/2], [cent-w/2, cent+w/2]
    
    if x[1] - r < 0 and x[1] > 0:  # Near pin level
        ball_center = x[0]
        ball_left = x[0] - r
        ball_right = x[0] + r
        
        # Define pin and hole boundaries
        left_pin_right = -hole_center - hole_width / 2
        left_hole_left = -hole_center - hole_width / 2
        left_hole_right = -hole_center + hole_width / 2
        center_pin_left = -hole_center + hole_width / 2
        center_pin_right = hole_center - hole_width / 2
        right_hole_left = hole_center - hole_width / 2
        right_hole_right = hole_center + hole_width / 2
        right_pin_left = hole_center + hole_width / 2
        
        # Check if ball is in a hole (completely through)
        in_left_hole = ball_left > left_hole_left and ball_right < left_hole_right
        in_right_hole = ball_left > right_hole_left and ball_right < right_hole_right
        
        # Check if ball hits pins
        hits_left_pins = ball_left < left_pin_right
        hits_center_pins = ball_right > center_pin_left and ball_left < center_pin_right
        hits_right_pins = ball_right > right_pin_left
        
        if not (in_left_hole or in_right_hole):
            # Not cleanly in a hole - check for pin collision
            if hits_left_pins or hits_center_pins or hits_right_pins:
                # BANG! Hit the pins
                state['last_collision'] = 'bang'
                # Reset ball
                x = np.array([0.0, 10.0])
                v = np.array([0.0, -1.0])
    
    # Check if ball fell below floor (success if it got past pins!)
    if x[1] + r < box_ylim[0]:
        # If we got here, we passed through a hole successfully
        # (otherwise we would have hit pins above)
        state['score'] += 1
        state['last_collision'] = 'success'
        
        # Reset ball
        x = np.array([0.0, 10.0])
        v = np.array([0.0, -1.0])
    
    # Update state
    state['x'] = x
    state['v'] = v
    
    return state


def apply_control(state: dict, control_input: Literal['left', 'right']) -> dict:
    """
    Apply control input to ball (convenience function).
    
    This is a wrapper around simulate_step for easier control application.
    
    Args:
        state: Simulation state dict
        control_input: Control signal ('left' or 'right')
        
    Returns:
        dict: Updated simulation state
    """
    return simulate_step(state, dt=0.01, control_input=control_input)


def get_pin_positions(state: dict, pin_interval: float = 0.4) -> tuple[np.ndarray, np.ndarray]:
    """
    Get positions of pins for visualisation.
    
    Args:
        state: Simulation state dict
        pin_interval: Spacing between pins
        
    Returns:
        tuple: (pin_x, pin_y) arrays for plotting pins
    """
    box_xlim = state['box_xlim']
    hole_center = state['hole_center']
    hole_width = state['hole_width']
    pin_height = state['pin_height']
    
    # Three sections of pins with two holes
    left_wall_right = -hole_center - hole_width / 2
    hole_left = -hole_center + hole_width / 2
    hole_right = hole_center - hole_width / 2
    right_wall_left = hole_center + hole_width / 2
    
    num_left_pins = int(round((left_wall_right - box_xlim[0]) / pin_interval))
    num_center_pins = int(round((hole_right - hole_left) / pin_interval))
    num_right_pins = int(round((box_xlim[1] - right_wall_left) / pin_interval))
    
    # Create pin positions
    pins_x = []
    pins_y = []
    
    # Left section
    for i in range(num_left_pins):
        px = box_xlim[0] + i * pin_interval
        pins_x.extend([px, px, np.nan])
        pins_y.extend([-pin_height, 0, np.nan])
    
    # Center section
    for i in range(num_center_pins):
        px = hole_left + i * pin_interval
        pins_x.extend([px, px, np.nan])
        pins_y.extend([-pin_height, 0, np.nan])
    
    # Right section
    for i in range(num_right_pins):
        px = right_wall_left + i * pin_interval
        pins_x.extend([px, px, np.nan])
        pins_y.extend([-pin_height, 0, np.nan])
    
    return np.array(pins_x), np.array(pins_y)


def get_average_energy(state: dict) -> float:
    """
    Calculate average energy expenditure per successful drop.
    
    Args:
        state: Simulation state dict
        
    Returns:
        float: Average energy, or np.nan if no successes yet
    """
    if state['score'] == 0:
        return np.nan
    return state['energy_count'] / state['score']
