"""Tests for visualisation module."""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt
import pytest

from kappenball import billiards, visualisation


def test_setup_figure_default():
    """Test figure setup with default parameters."""
    fig, ax = visualisation.setup_figure()
    
    assert fig is not None
    assert ax is not None
    assert ax.get_xlim() == (0.0, 10.0)
    assert ax.get_ylim() == (0.0, 10.0)
    
    plt.close(fig)


def test_setup_figure_custom():
    """Test figure setup with custom parameters."""
    fig, ax = visualisation.setup_figure(
        box_xlim=(-5.0, 5.0),
        box_ylim=(-3.0, 3.0),
        figsize=(8, 6)
    )
    
    assert ax.get_xlim() == (-5.0, 5.0)
    assert ax.get_ylim() == (-3.0, 3.0)
    assert fig.get_figwidth() == 8
    assert fig.get_figheight() == 6
    
    plt.close(fig)


def test_plot_billiards_frame():
    """Test plotting a single billiards frame."""
    state = billiards.initialise(n_particles=5, random_state=42)
    fig, ax = visualisation.setup_figure(
        box_xlim=state['box_xlim'],
        box_ylim=state['box_ylim']
    )
    
    artists = visualisation.plot_billiards_frame(ax, state)
    
    assert len(artists) > 0  # Should have created some artists
    
    plt.close(fig)


def test_plot_billiards_frame_with_velocities():
    """Test plotting frame with velocity vectors."""
    state = billiards.initialise(n_particles=3, random_state=42)
    fig, ax = visualisation.setup_figure(
        box_xlim=state['box_xlim'],
        box_ylim=state['box_ylim']
    )
    
    artists_with_vel = visualisation.plot_billiards_frame(
        ax, state, show_velocities=True
    )
    
    plt.cla()  # Clear axis
    
    artists_without_vel = visualisation.plot_billiards_frame(
        ax, state, show_velocities=False
    )
    
    # With velocities should have more artists (particles + arrows)
    assert len(artists_with_vel) >= len(artists_without_vel)
    
    plt.close(fig)


def test_plot_billiards_static():
    """Test static billiards plot."""
    state = billiards.initialise(n_particles=5, random_state=42)
    
    fig, ax = visualisation.plot_billiards_static(state)
    
    assert fig is not None
    assert ax is not None
    
    plt.close(fig)


def test_animate_billiards_creation():
    """Test that animation object can be created (doesn't actually run it)."""
    state = billiards.initialise(n_particles=3, random_state=42)
    
    fig, anim = visualisation.animate_billiards(
        state,
        n_frames=10,  # Just a few frames for testing
        interval=50
    )
    
    assert fig is not None
    assert anim is not None
    
    plt.close(fig)


def test_plot_velocity_histogram():
    """Test velocity histogram plotting."""
    # Create some states by simulating
    state = billiards.initialise(n_particles=10, random_state=42)
    states = [state]
    
    for _ in range(50):
        state = billiards.simulate_step(state)
        states.append(state.copy())
    
    fig, ax = visualisation.plot_velocity_histogram(states)
    
    assert fig is not None
    assert ax is not None
    
    plt.close(fig)


def test_plot_energy_evolution():
    """Test energy evolution plotting."""
    # Create some states by simulating
    state = billiards.initialise(n_particles=5, random_state=42)
    states = [state]
    
    for _ in range(100):
        state = billiards.simulate_step(state, dt=0.001)
        states.append(state.copy())
    
    fig, ax = visualisation.plot_energy_evolution(states)
    
    assert fig is not None
    assert ax is not None
    
    plt.close(fig)


def test_clear_artists():
    """Test clearing artists from plot."""
    state = billiards.initialise(n_particles=3, random_state=42)
    fig, ax = visualisation.setup_figure(
        box_xlim=state['box_xlim'],
        box_ylim=state['box_ylim']
    )
    
    artists = visualisation.plot_billiards_frame(ax, state)
    initial_num_children = len(ax.get_children())
    
    visualisation.clear_artists(artists)
    
    # After clearing, should have fewer children
    assert len(ax.get_children()) < initial_num_children
    
    plt.close(fig)


def test_plot_different_init_types():
    """Test plotting works for different initialization types."""
    init_types = ['rand', 'randn', 'uniform', 'hotCold']
    
    for init_type in init_types:
        state = billiards.initialise(
            n_particles=5,
            init_type=init_type,
            random_state=42
        )
        fig, ax = visualisation.plot_billiards_static(state)
        
        assert fig is not None
        assert ax is not None
        
        plt.close(fig)


def test_velocity_scale():
    """Test that velocity scale affects arrow size."""
    state = billiards.initialise(n_particles=3, random_state=42)
    fig, ax = visualisation.setup_figure(
        box_xlim=state['box_xlim'],
        box_ylim=state['box_ylim']
    )
    
    # Plot with different scales (just verify it doesn't crash)
    for scale in [0.1, 0.5, 1.0]:
        plt.cla()
        artists = visualisation.plot_billiards_frame(
            ax, state, velocity_scale=scale
        )
        assert len(artists) > 0
    
    plt.close(fig)

