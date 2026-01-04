"""
Visualisation module using matplotlib.

Provides functions for visualising billiards and falling ball simulations,
with support for animations and interactive controls.

Ported from MATLAB (2012) to Python (2026).

Functions:
    setup_figure: Configure matplotlib figure for simulation
    plot_billiards_frame: Render single frame of billiards
    plot_falling_ball_frame: Render single frame of falling ball
    animate_billiards: Create animated billiards visualisation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle, Rectangle
from matplotlib.collections import PatchCollection
from typing import Optional, List, Tuple


def setup_figure(
    box_xlim: Tuple[float, float] = (0.0, 10.0),
    box_ylim: Tuple[float, float] = (0.0, 10.0),
    figsize: Tuple[float, float] = (10, 10),
    title: str = "Billiards Simulation"
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Set up matplotlib figure for simulation visualisation.
    
    Args:
        box_xlim: (xmin, xmax) for simulation box
        box_ylim: (ymin, ymax) for simulation box
        figsize: Figure size in inches (width, height)
        title: Figure title
        
    Returns:
        tuple: (fig, ax) matplotlib figure and axis objects
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Set axis limits and aspect ratio
    ax.set_xlim(box_xlim)
    ax.set_ylim(box_ylim)
    ax.set_aspect('equal')
    
    # Style the plot
    ax.set_facecolor('#4080CC')  # Blue billiard table colour
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    
    # Add box boundary
    box_width = box_xlim[1] - box_xlim[0]
    box_height = box_ylim[1] - box_ylim[0]
    boundary = Rectangle(
        (box_xlim[0], box_ylim[0]),
        box_width,
        box_height,
        fill=False,
        edgecolor='white',
        linewidth=3
    )
    ax.add_patch(boundary)
    
    return fig, ax


def plot_billiards_frame(
    ax: plt.Axes,
    state: dict,
    show_velocities: bool = True,
    show_trails: bool = False,
    velocity_scale: float = 0.5
) -> List:
    """
    Plot a single frame of billiards simulation.
    
    Args:
        ax: Matplotlib axis to plot on
        state: Simulation state dict from billiards.initialise()
        show_velocities: Whether to show velocity vectors
        show_trails: Whether to show particle trails (not yet implemented)
        velocity_scale: Scale factor for velocity arrows
        
    Returns:
        list: Artist objects (circles and arrows) for animation update
    """
    X = state['X']
    V = state['V']
    r = state['r']
    colors = state['colors']
    n_particles = state['n_particles']
    
    artists = []
    
    # Plot particles as circles
    for i in range(n_particles):
        circle = Circle(
            (X[i, 0], X[i, 1]),
            r[i],
            facecolor=colors[i],
            edgecolor='black',
            linewidth=1,
            zorder=10
        )
        ax.add_patch(circle)
        artists.append(circle)
    
    # Plot velocity vectors
    if show_velocities:
        for i in range(n_particles):
            if np.linalg.norm(V[i]) > 1e-6:  # Only show if moving
                # Compute arrow start and end points
                vel_norm = V[i] / (np.linalg.norm(V[i]) + 1e-10)
                arrow_start = X[i] + r[i] * vel_norm
                arrow_vec = V[i] * velocity_scale
                
                arrow = ax.arrow(
                    arrow_start[0],
                    arrow_start[1],
                    arrow_vec[0],
                    arrow_vec[1],
                    head_width=0.15,
                    head_length=0.2,
                    fc=colors[i],
                    ec='black',
                    linewidth=1.5,
                    zorder=5,
                    alpha=0.8
                )
                artists.append(arrow)
    
    return artists


def clear_artists(artists: List) -> None:
    """
    Remove artist objects from plot.
    
    Args:
        artists: List of matplotlib artist objects to remove
    """
    for artist in artists:
        artist.remove()


def animate_billiards(
    initial_state: dict,
    n_frames: int = 500,
    dt: float = 2e-3,
    interval: int = 20,
    show_velocities: bool = True,
    velocity_scale: float = 0.5,
    title: str = "Billiards Simulation"
) -> Tuple[plt.Figure, FuncAnimation]:
    """
    Create animated visualisation of billiards simulation.
    
    Args:
        initial_state: Initial simulation state from billiards.initialise()
        n_frames: Number of frames to simulate
        dt: Time step for physics simulation
        interval: Milliseconds between frames (for display)
        show_velocities: Whether to show velocity vectors
        velocity_scale: Scale factor for velocity arrows
        title: Figure title
        
    Returns:
        tuple: (fig, anim) matplotlib figure and animation object
    """
    from kappenball.billiards import simulate_step
    
    # Set up figure
    fig, ax = setup_figure(
        box_xlim=initial_state['box_xlim'],
        box_ylim=initial_state['box_ylim'],
        title=title
    )
    
    # Store state and artists
    state = initial_state.copy()
    artists = []
    
    def init():
        """Initialize animation."""
        return []
    
    def update(frame):
        """Update function for animation."""
        nonlocal state, artists
        
        # Clear previous artists
        clear_artists(artists)
        
        # Simulate physics step
        state = simulate_step(state, dt=dt)
        
        # Plot new frame
        artists = plot_billiards_frame(
            ax,
            state,
            show_velocities=show_velocities,
            velocity_scale=velocity_scale
        )
        
        # Update title with frame number
        ax.set_title(f"{title} (frame {frame}/{n_frames})")
        
        return artists
    
    # Create animation
    anim = FuncAnimation(
        fig,
        update,
        frames=n_frames,
        init_func=init,
        interval=interval,
        blit=False,  # Set to False for better compatibility
        repeat=True
    )
    
    return fig, anim


def plot_billiards_static(
    state: dict,
    show_velocities: bool = True,
    velocity_scale: float = 0.5,
    title: str = "Billiards Simulation",
    save_path: Optional[str] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a static plot of billiards state.
    
    Useful for quick visualisation or saving snapshots.
    
    Args:
        state: Simulation state dict
        show_velocities: Whether to show velocity vectors
        velocity_scale: Scale factor for velocity arrows
        title: Figure title
        save_path: Optional path to save figure
        
    Returns:
        tuple: (fig, ax) matplotlib figure and axis objects
    """
    fig, ax = setup_figure(
        box_xlim=state['box_xlim'],
        box_ylim=state['box_ylim'],
        title=title
    )
    
    plot_billiards_frame(
        ax,
        state,
        show_velocities=show_velocities,
        velocity_scale=velocity_scale
    )
    
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig, ax


def plot_velocity_histogram(
    states: List[dict],
    bins: int = 30,
    title: str = "Velocity Distribution"
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot histogram of particle velocities over time.
    
    Useful for observing thermalization and Maxwell-Boltzmann distribution.
    
    Args:
        states: List of simulation states over time
        bins: Number of histogram bins
        title: Figure title
        
    Returns:
        tuple: (fig, ax) matplotlib figure and axis objects
    """
    # Collect all velocity magnitudes
    velocities = []
    for state in states:
        V = state['V']
        v_magnitudes = np.linalg.norm(V, axis=1)
        velocities.extend(v_magnitudes)
    
    # Create histogram
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(velocities, bins=bins, density=True, alpha=0.7, color='blue', edgecolor='black')
    ax.set_xlabel('Velocity Magnitude', fontsize=12)
    ax.set_ylabel('Probability Density', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    return fig, ax


def plot_energy_evolution(
    states: List[dict],
    title: str = "Energy Evolution"
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot total kinetic energy over time.
    
    Useful for verifying energy conservation.
    
    Args:
        states: List of simulation states over time
        title: Figure title
        
    Returns:
        tuple: (fig, ax) matplotlib figure and axis objects
    """
    # Compute energy at each time step
    energies = []
    for state in states:
        V = state['V']
        mass = state['mass']
        energy = 0.5 * np.sum(mass[:, np.newaxis] * V**2)
        energies.append(energy)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(energies, linewidth=2, color='blue')
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Total Kinetic Energy', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add mean energy line
    mean_energy = np.mean(energies)
    ax.axhline(mean_energy, color='red', linestyle='--', linewidth=2, 
               label=f'Mean: {mean_energy:.2e}')
    ax.legend()
    
    return fig, ax
