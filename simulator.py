#!/usr/bin/env python3
"""
MAGVID Electromagnetic Particle Simulator
==========================================

Main entry point for the MAGVID simulation.

Simulates charged particle trajectories in the MAGVID configuration:
- Rotating magnetic field from 4-coil assembly (creates radial electric field)
- Static axial magnetic field for confinement
- Visualizes particle accumulation and trajectories

Based on the analysis from "GroomLake Colonel Reveals All" document
"""

import sys
import os
import time
import matplotlib
# Use non-interactive backend for headless environments
if '--no-display' in sys.argv or 'DISPLAY' not in os.environ:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config import MAGVIDConfig, ConfigurationError
from physics import MAGVIDPhysics, SimulationError
from visualization import MAGVIDVisualizer


def main():
    """
    Run the MAGVID simulation.

    Returns:
        int: Exit code (0 for success, 1 for error)
    """
    print("MAGVID Electromagnetic Particle Simulator")
    print("=========================================\n")

    # Load and validate configuration
    try:
        config = MAGVIDConfig('config.ini')
        config.display_parameters()
        print()
    except ConfigurationError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"ERROR: Unexpected error loading configuration: {e}", file=sys.stderr)
        return 1

    # Initialize physics engine
    try:
        physics = MAGVIDPhysics(config)
    except Exception as e:
        print(f"ERROR: Failed to initialize physics engine: {e}", file=sys.stderr)
        return 1

    # Initialize visualizer
    try:
        visualizer = MAGVIDVisualizer(config, physics)
    except Exception as e:
        print(f"ERROR: Failed to initialize visualizer: {e}", file=sys.stderr)
        return 1

    # Run multi-particle simulation
    print("Running multi-particle simulation...")
    start_time = time.time()

    def progress_callback(i, n):
        """Print progress for particle simulation."""
        print(f"  Particle {i}/{n} simulated")

    try:
        trajectories, particles_info = physics.run_multi_particle_simulation(
            progress_callback=progress_callback
        )
    except SimulationError as e:
        print(f"\nERROR: Simulation failed: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\n\nSimulation interrupted by user", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"\nERROR: Unexpected error during simulation: {e}", file=sys.stderr)
        return 1

    end_time = time.time()
    print(f"Simulation completed in {end_time - start_time:.2f} seconds\n")

    # Analyze results
    try:
        final_positions, charges = visualizer.analyze_accumulation(trajectories, particles_info)
    except Exception as e:
        print(f"WARNING: Analysis failed: {e}", file=sys.stderr)

    # Create visualizations
    print("\nCreating visualizations...")

    try:
        # 3D trajectory plot
        if config.plot_3d or config.plot_trajectories:
            fig1 = visualizer.visualize_3d_trajectories(trajectories, particles_info)
            if config.save_plots:
                visualizer.save_figure(fig1, 'magvid_trajectories.png')

        # Field visualization
        if config.plot_fields:
            fig2 = visualizer.plot_field_visualization()
            if config.save_plots:
                visualizer.save_figure(fig2, 'magvid_fields.png')

        # Show plots (only if not in headless mode)
        if '--no-display' not in sys.argv and 'DISPLAY' in os.environ:
            try:
                plt.show()
            except Exception as e:
                print(f"Note: Could not display plots (running in headless mode?): {e}")

    except Exception as e:
        print(f"WARNING: Visualization failed: {e}", file=sys.stderr)

    # Print summary
    print("\n" + "="*50)
    print("Simulation complete!")
    print("="*50)
    print("Red trajectories: Negatively charged particles")
    print("Blue trajectories: Positively charged particles")
    print("Circles: Starting positions")
    print("Squares: Final positions")
    print("="*50)

    return 0


if __name__ == "__main__":
    sys.exit(main())
