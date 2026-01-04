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
from typing import Optional, List, Tuple, Any


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


# ============================================================================
# Falling Ball (Kappenball) Visualisation Functions
# ============================================================================


def setup_kappenball_figure(
    state: dict,
    figsize: tuple[float, float] = (10, 12),
    screen_colour: tuple[float, float, float] = (0.4, 0.5, 0.8)
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Set up figure specifically for Kappenball game visualisation.
    
    Args:
        state: Simulation state dict from falling_ball.initialise()
        figsize: Figure size (width, height) in inches
        screen_colour: RGB colour for background
        
    Returns:
        tuple: (fig, ax) matplotlib figure and axis objects
    """
    from kappenball import falling_ball
    
    fig, ax = plt.subplots(figsize=figsize)
    
    xlim = state['box_xlim']
    ylim = state['box_ylim']
    hole_center = state['hole_center']
    hole_width = state['hole_width']
    pin_height = state['pin_height']
    
    # Set up axes
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect('equal')
    ax.set_facecolor(screen_colour)
    ax.axis('off')
    
    # Draw walls (three sections with two holes)
    box_width = 5
    left_wall_right = -hole_center - hole_width / 2
    hole_left = -hole_center + hole_width / 2
    hole_right = hole_center - hole_width / 2
    right_wall_left = hole_center + hole_width / 2
    
    # Left wall
    ax.plot([xlim[0], xlim[0], left_wall_right, left_wall_right],
            [ylim[1], -pin_height, -pin_height, ylim[0]],
            'k-', linewidth=box_width)
    
    # Center wall (between two holes)
    ax.plot([hole_left, hole_left, hole_right, hole_right],
            [ylim[0], -pin_height, -pin_height, ylim[0]],
            'k-', linewidth=box_width)
    
    # Right wall
    ax.plot([right_wall_left, right_wall_left, xlim[1], xlim[1]],
            [ylim[0], -pin_height, -pin_height, ylim[1]],
            'k-', linewidth=box_width)
    
    # Draw pins
    pins_x, pins_y = falling_ball.get_pin_positions(state)
    ax.plot(pins_x, pins_y, 'b-', linewidth=2)
    
    return fig, ax


def plot_kappenball_frame(
    ax: plt.Axes,
    state: dict,
    show_stats: bool = True
) -> List:
    """
    Plot one frame of Kappenball game state.
    
    Args:
        ax: Matplotlib axis to plot on
        state: Simulation state dict
        show_stats: Whether to show score/energy statistics
        
    Returns:
        list: List of artist objects for animation
    """
    from kappenball import falling_ball
    
    artists = []
    
    # Draw ball
    ball = plt.Circle(
        state['x'],
        state['r'],
        color='red',
        zorder=10
    )
    ax.add_patch(ball)
    artists.append(ball)
    
    # Add statistics text if requested
    if show_stats:
        xlim = state['box_xlim']
        ylim = state['box_ylim']
        
        # Score
        score_txt = ax.text(
            xlim[0], ylim[1] + 1,
            f"Score: {state['score']}",
            fontsize=20,
            ha='right',
            va='bottom'
        )
        artists.append(score_txt)
        
        # Energy count
        energy_txt = ax.text(
            xlim[1], ylim[1] + 1,
            f"Energy: {state['energy_count']}",
            fontsize=20,
            ha='right',
            va='bottom'
        )
        artists.append(energy_txt)
        
        # Average energy per score
        avg_energy = falling_ball.get_average_energy(state)
        avg_txt = '-' if np.isnan(avg_energy) else f'{avg_energy:.2f}'
        average_txt = ax.text(
            (xlim[1] - xlim[0]) / 2 + xlim[0], ylim[1] + 1,
            f"Average: {avg_txt}",
            fontsize=20,
            ha='right',
            va='bottom'
        )
        artists.append(average_txt)
        
        # Show collision feedback
        if state.get('last_collision') == 'bang':
            bang_txt = ax.text(
                state['x'][0], state['x'][1],
                'BANG!',
                fontsize=24,
                color='red',
                ha='center',
                va='center',
                weight='bold'
            )
            artists.append(bang_txt)
    
    return artists


def animate_kappenball(
    state: dict,
    num_steps: int = 500,
    dt: float = 0.01,
    interval: int = 20,
    figsize: tuple[float, float] = (10, 12),
    control_sequence: Optional[list] = None
) -> Tuple[plt.Figure, Any]:
    """
    Create animated visualisation of Kappenball game.
    
    Args:
        state: Initial simulation state
        num_steps: Number of simulation steps to run
        dt: Time step for simulation
        interval: Delay between frames in milliseconds
        figsize: Figure size (width, height)
        control_sequence: Optional list of control inputs ('left', 'right', None)
                         for each step
        
    Returns:
        tuple: (fig, animation) matplotlib figure and animation objects
    """
    from matplotlib import animation
    from kappenball import falling_ball
    
    # Generate trajectory
    trajectory = [state.copy()]
    current_state = state.copy()
    
    for i in range(num_steps):
        control = None
        if control_sequence and i < len(control_sequence):
            control = control_sequence[i]
        
        current_state = falling_ball.simulate_step(
            current_state,
            dt=dt,
            control_input=control
        )
        trajectory.append(current_state.copy())
    
    # Set up figure
    fig, ax = setup_kappenball_figure(state, figsize=figsize)
    
    # Animation update function
    def update(frame):
        # Clear previous artists
        for artist in ax.patches[:]:
            if isinstance(artist, plt.Circle):
                artist.remove()
        for txt in ax.texts[:]:
            txt.remove()
        
        # Plot new frame
        artists = plot_kappenball_frame(ax, trajectory[frame])
        return artists
    
    # Create animation
    anim = animation.FuncAnimation(
        fig,
        update,
        frames=len(trajectory),
        interval=interval,
        blit=False,
        repeat=True
    )
    
    return fig, anim


def plot_kappenball_static(
    state: dict,
    show_stats: bool = True,
    title: str = "Kappenball",
    save_path: Optional[str] = None,
    figsize: tuple[float, float] = (10, 12)
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a static plot of Kappenball game state.
    
    Args:
        state: Simulation state dict
        show_stats: Whether to show score/energy statistics
        title: Figure title
        save_path: Optional path to save figure
        figsize: Figure size (width, height)
        
    Returns:
        tuple: (fig, ax) matplotlib figure and axis objects
    """
    fig, ax = setup_kappenball_figure(state, figsize=figsize)
    
    if title:
        fig.suptitle(title, fontsize=16)
    
    plot_kappenball_frame(ax, state, show_stats=show_stats)
    
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig, ax
