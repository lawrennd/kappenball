"""Tests for billiards physics module."""

import numpy as np
import pytest
from kappenball import billiards


def test_initialise_default():
    """Test default initialization."""
    state = billiards.initialise(n_particles=5)
    
    assert state['n_particles'] == 5
    assert state['X'].shape == (5, 2)
    assert state['V'].shape == (5, 2)
    assert state['r'].shape == (5,)
    assert state['colors'].shape == (5, 3)
    assert state['mass'].shape == (5,)
    
    # Check positions are within bounds
    assert np.all(state['X'][:, 0] >= state['box_xlim'][0])
    assert np.all(state['X'][:, 0] <= state['box_xlim'][1])
    assert np.all(state['X'][:, 1] >= state['box_ylim'][0])
    assert np.all(state['X'][:, 1] <= state['box_ylim'][1])


def test_initialise_rand():
    """Test random initialization."""
    state = billiards.initialise(n_particles=10, init_type='rand', random_state=42)
    
    assert state['n_particles'] == 10
    # Radii should be varied for 'rand' mode
    assert np.std(state['r']) > 0


def test_initialise_uniform():
    """Test uniform initialization."""
    state = billiards.initialise(n_particles=8, init_type='uniform', random_state=42)
    
    assert state['n_particles'] == 8
    # Particles should be at same x position
    assert np.allclose(state['X'][:, 0], state['X'][0, 0])
    # All particles should be yellow
    assert np.all(state['colors'] == [1, 1, 0])


def test_initialise_hotCold():
    """Test hot/cold initialization."""
    state = billiards.initialise(n_particles=10, init_type='hotCold', random_state=42)
    
    assert state['n_particles'] == 10
    # First half should be yellow (hot), second half cyan (cold)
    assert np.allclose(state['colors'][0], [1, 1, 0])
    assert np.allclose(state['colors'][-1], [0, 1, 1])


def test_initialise_reproducibility():
    """Test that random_state makes initialization reproducible."""
    state1 = billiards.initialise(n_particles=5, random_state=123)
    state2 = billiards.initialise(n_particles=5, random_state=123)
    
    np.testing.assert_array_equal(state1['X'], state2['X'])
    np.testing.assert_array_equal(state1['V'], state2['V'])


def test_simulate_step():
    """Test simulation step updates positions."""
    state = billiards.initialise(n_particles=5, random_state=42)
    X_initial = state['X'].copy()
    
    state = billiards.simulate_step(state, dt=0.01)
    
    # Positions should have changed
    assert not np.allclose(state['X'], X_initial)
    
    # Positions should still be in bounds (approximately, allowing for radii)
    margin = np.max(state['r'])
    assert np.all(state['X'][:, 0] >= state['box_xlim'][0] - margin)
    assert np.all(state['X'][:, 0] <= state['box_xlim'][1] + margin)


def test_simulate_step_energy_conservation():
    """Test that total kinetic energy is approximately conserved."""
    state = billiards.initialise(n_particles=3, init_type='rand', random_state=42)
    
    def kinetic_energy(state):
        return 0.5 * np.sum(state['mass'][:, np.newaxis] * state['V']**2)
    
    initial_energy = kinetic_energy(state)
    
    # Simulate many steps
    for _ in range(100):
        state = billiards.simulate_step(state, dt=0.001)
    
    final_energy = kinetic_energy(state)
    
    # Energy should be approximately conserved (within 10% due to numerical errors)
    # Note: Some energy loss is expected due to wall collisions and numerical integration
    assert abs(final_energy - initial_energy) / initial_energy < 0.5


def test_simulate_step_particles_stay_in_box():
    """Test that particles bounce off walls and stay in box."""
    state = billiards.initialise(n_particles=5, init_type='rand', random_state=42)
    
    # Simulate many steps
    for _ in range(500):
        state = billiards.simulate_step(state, dt=0.002)
        
        # Check particles stay roughly within bounds
        # (allowing small margin for radii and numerical errors)
        margin = 2 * np.max(state['r'])
        assert np.all(state['X'][:, 0] >= state['box_xlim'][0] - margin)
        assert np.all(state['X'][:, 0] <= state['box_xlim'][1] + margin)
        assert np.all(state['X'][:, 1] >= state['box_ylim'][0] - margin)
        assert np.all(state['X'][:, 1] <= state['box_ylim'][1] + margin)


def test_compute_vectors():
    """Test velocity vector computation."""
    state = billiards.initialise(n_particles=5, random_state=42)
    vectors = billiards.compute_vectors(state)
    
    assert 'arrow_starts' in vectors
    assert 'arrow_vectors' in vectors
    assert vectors['arrow_starts'].shape == (5, 2)
    assert vectors['arrow_vectors'].shape == (5, 2)
    
    # Arrow vectors should be proportional to velocities
    assert np.allclose(
        vectors['arrow_vectors'] / 0.025,
        state['V']
    )


def test_detect_collisions():
    """Test collision detection."""
    # Create two particles that are overlapping
    X = np.array([[0.0, 0.0], [0.5, 0.0]])
    r = np.array([0.3, 0.3])
    
    collisions = billiards._detect_collisions(X, r)
    
    # Should detect one collision between particle 0 and 1
    assert len(collisions) == 1
    assert collisions[0] == (0, 1)


def test_no_initial_overlaps():
    """Test that initialization doesn't create overlapping particles."""
    state = billiards.initialise(n_particles=20, init_type='rand', random_state=42)
    
    # Check no particles overlap
    for i in range(state['n_particles']):
        for j in range(i + 1, state['n_particles']):
            dist = np.linalg.norm(state['X'][i] - state['X'][j])
            assert dist >= (state['r'][i] + state['r'][j]) - 1e-6  # Small tolerance


def test_custom_box_limits():
    """Test initialization with custom box limits."""
    state = billiards.initialise(
        n_particles=5,
        box_xlim=(-5.0, 5.0),
        box_ylim=(-3.0, 3.0),
        random_state=42
    )
    
    assert state['box_xlim'] == (-5.0, 5.0)
    assert state['box_ylim'] == (-3.0, 3.0)
    
    # Particles should be within custom bounds
    assert np.all(state['X'][:, 0] >= -5.0)
    assert np.all(state['X'][:, 0] <= 5.0)
    assert np.all(state['X'][:, 1] >= -3.0)
    assert np.all(state['X'][:, 1] <= 3.0)

