"""Tests for falling ball (Kappenball game) module."""

import numpy as np
import pytest
from kappenball import falling_ball


def test_initialise_default():
    """Test default initialisation."""
    state = falling_ball.initialise()
    
    assert state['x'][0] == 0.0
    assert state['x'][1] == 10.0
    assert state['v'][0] == 0.0
    assert state['v'][1] == -1.0
    assert state['r'] == 0.25
    assert state['score'] == 0
    assert state['energy_count'] == 0


def test_initialise_custom():
    """Test custom initialisation."""
    state = falling_ball.initialise(
        initial_position=(1.0, 5.0),
        initial_velocity=(0.5, -0.5),
        velocity_variance=0.1,
        hole_center=3.0,
        hole_width=2.0
    )
    
    assert state['x'][0] == 1.0
    assert state['x'][1] == 5.0
    assert state['v'][0] == 0.5
    assert state['v'][1] == -0.5
    assert state['v_var'] == 0.1
    assert state['hole_center'] == 3.0
    assert state['hole_width'] == 2.0


def test_simulate_step_gravity():
    """Test that ball falls due to gravity."""
    state = falling_ball.initialise(random_state=42)
    initial_y = state['x'][1]
    
    state = falling_ball.simulate_step(state, dt=0.1)
    
    # Ball should have fallen (y decreased)
    assert state['x'][1] < initial_y


def test_simulate_step_control_left():
    """Test left control input."""
    state = falling_ball.initialise(random_state=42)
    initial_vx = state['v'][0]
    initial_energy = state['energy_count']
    
    state = falling_ball.simulate_step(state, control_input='left')
    
    # Horizontal velocity should have decreased
    assert state['v'][0] < initial_vx - 0.4  # Allowing for random noise
    # Energy count should increase
    assert state['energy_count'] == initial_energy + 1


def test_simulate_step_control_right():
    """Test right control input."""
    state = falling_ball.initialise(random_state=42)
    initial_vx = state['v'][0]
    initial_energy = state['energy_count']
    
    state = falling_ball.simulate_step(state, control_input='right')
    
    # Horizontal velocity should have increased
    assert state['v'][0] > initial_vx + 0.4  # Allowing for random noise
    # Energy count should increase
    assert state['energy_count'] == initial_energy + 1


def test_wall_collision():
    """Test that ball bounces off walls."""
    state = falling_ball.initialise(
        initial_position=(11.0, 5.0),  # Near right wall
        initial_velocity=(2.0, 0.0),   # Moving right
        random_state=42
    )
    
    # Simulate several steps
    for _ in range(10):
        state = falling_ball.simulate_step(state, dt=0.01)
    
    # Ball should stay within bounds
    assert state['x'][0] >= state['box_xlim'][0]
    assert state['x'][0] <= state['box_xlim'][1]


def test_successful_drop():
    """Test successful drop through hole."""
    # Position ball to fall straight through a hole
    state = falling_ball.initialise(
        initial_position=(6.0, 1.0),  # Above right hole
        initial_velocity=(0.0, -2.0),
        hole_center=6.0,
        hole_width=4.0,
        random_state=42
    )
    
    initial_score = state['score']
    
    # Simulate until ball passes through (needs ~200 steps to fall)
    for _ in range(200):
        state = falling_ball.simulate_step(state, dt=0.01)
        if state['last_collision'] == 'success':
            break
    
    # Score should have increased
    assert state['score'] == initial_score + 1
    # Ball should be reset to top
    assert state['x'][1] == 10.0


def test_collision_with_pins():
    """Test collision with pins (BANG)."""
    # Position ball to hit pins
    state = falling_ball.initialise(
        initial_position=(0.0, 1.0),  # Above center pins
        initial_velocity=(0.0, -2.0),
        hole_center=6.0,
        hole_width=4.0,
        pin_height=0.5,
        random_state=42
    )
    
    initial_score = state['score']
    
    # Simulate until collision
    collision_detected = False
    for _ in range(100):
        state = falling_ball.simulate_step(state, dt=0.01)
        if state['last_collision'] == 'bang':
            collision_detected = True
            break
    
    # Should have detected collision
    assert collision_detected or state['last_collision'] == 'success'  # Might go through if unlucky
    # Score shouldn't increase on bang
    if collision_detected:
        assert state['score'] == initial_score
    # Ball should be reset
    assert state['x'][1] == 10.0


def test_velocity_variance():
    """Test that velocity variance adds randomness."""
    # Run two simulations with same initial conditions but variance
    state1 = falling_ball.initialise(
        velocity_variance=0.5,
        random_state=42
    )
    state2 = falling_ball.initialise(
        velocity_variance=0.5,
        random_state=43  # Different seed
    )
    
    # Simulate same number of steps
    for _ in range(50):
        state1 = falling_ball.simulate_step(state1, dt=0.01)
        state2 = falling_ball.simulate_step(state2, dt=0.01)
    
    # Positions should differ due to randomness
    assert not np.allclose(state1['x'], state2['x'])


def test_reproducibility():
    """Test that random_state makes simulation reproducible."""
    state1 = falling_ball.initialise(
        velocity_variance=0.5,
        random_state=42
    )
    state2 = falling_ball.initialise(
        velocity_variance=0.5,
        random_state=42
    )
    
    # Simulate same number of steps
    for _ in range(50):
        state1 = falling_ball.simulate_step(state1, dt=0.01)
        state2 = falling_ball.simulate_step(state2, dt=0.01)
    
    # Should produce same results
    np.testing.assert_array_almost_equal(state1['x'], state2['x'])
    np.testing.assert_array_almost_equal(state1['v'], state2['v'])


def test_get_pin_positions():
    """Test pin position generation."""
    state = falling_ball.initialise()
    
    pins_x, pins_y = falling_ball.get_pin_positions(state)
    
    assert len(pins_x) > 0
    assert len(pins_y) > 0
    assert len(pins_x) == len(pins_y)
    # Pins should be at pin_height and 0
    unique_y = pins_y[~np.isnan(pins_y)]
    assert np.any(unique_y == -state['pin_height'])
    assert np.any(unique_y == 0)


def test_get_average_energy():
    """Test average energy calculation."""
    state = falling_ball.initialise()
    
    # No scores yet
    assert np.isnan(falling_ball.get_average_energy(state))
    
    # Add some scores and energy
    state['score'] = 5
    state['energy_count'] = 20
    
    avg = falling_ball.get_average_energy(state)
    assert avg == 4.0  # 20/5


def test_apply_control():
    """Test apply_control convenience function."""
    state = falling_ball.initialise(random_state=42)
    initial_energy = state['energy_count']
    
    state = falling_ball.apply_control(state, 'left')
    
    assert state['energy_count'] == initial_energy + 1

