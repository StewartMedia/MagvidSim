#!/usr/bin/env python3
"""
Visualization Module for MAGVID Simulator
==========================================

Handles 3D trajectory plots, field visualizations, and analysis output.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class MAGVIDVisualizer:
    """Handles visualization and analysis for MAGVID simulation."""

    def __init__(self, config, physics):
        """
        Initialize visualizer.

        Args:
            config: MAGVIDConfig object
            physics: MAGVIDPhysics object
        """
        self.config = config
        self.physics = physics

    def visualize_3d_trajectories(self, trajectories, particles_info):
        """
        Create 3D visualization of particle trajectories.

        Args:
            trajectories: List of trajectory arrays
            particles_info: List of particle property dicts

        Returns:
            matplotlib.figure.Figure: The created figure
        """
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Plot trajectories
        for i, (traj, info) in enumerate(zip(trajectories, particles_info)):
            if len(traj) > 1:
                x, y, z = traj[:, 0], traj[:, 1], traj[:, 2]
                ax.plot(x, y, z, color=info['color'], alpha=0.7, linewidth=0.8)

                # Mark start and end points
                ax.scatter(x[0], y[0], z[0], color=info['color'], s=30, marker='o')
                ax.scatter(x[-1], y[-1], z[-1], color=info['color'], s=30, marker='s')

        # Draw geometry - gap region boundaries (cylindrical)
        theta = np.linspace(0, 2*np.pi, 50)

        # Inner cylinder
        x_inner = self.config.r_inner * np.cos(theta)
        y_inner = self.config.r_inner * np.sin(theta)
        z_inner_top = np.full_like(x_inner, self.config.z_gap/2)
        z_inner_bot = np.full_like(x_inner, -self.config.z_gap/2)

        ax.plot(x_inner, y_inner, z_inner_top, 'k-', alpha=0.5)
        ax.plot(x_inner, y_inner, z_inner_bot, 'k-', alpha=0.5)

        # Outer cylinder
        x_outer = self.config.r_outer * np.cos(theta)
        y_outer = self.config.r_outer * np.sin(theta)
        z_outer_top = np.full_like(x_outer, self.config.z_gap/2)
        z_outer_bot = np.full_like(x_outer, -self.config.z_gap/2)

        ax.plot(x_outer, y_outer, z_outer_top, 'k-', alpha=0.5)
        ax.plot(x_outer, y_outer, z_outer_bot, 'k-', alpha=0.5)

        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title('MAGVID Charged Particle Trajectories')

        # Legend
        ax.scatter([], [], [], color='red', label='Negative charges')
        ax.scatter([], [], [], color='blue', label='Positive charges')
        ax.legend()

        plt.tight_layout()
        return fig

    def plot_field_visualization(self):
        """
        Visualize the magnetic and electric field structure.

        Creates a 4-panel plot showing:
        - Radial magnetic field (Br)
        - Axial magnetic field (Bz)
        - Radial electric field (Er)
        - Azimuthal electric field (Eφ)

        Returns:
            matplotlib.figure.Figure: The created figure
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

        # Create grid for field visualization
        r_grid = np.linspace(0.05, 0.15, 30)
        z_grid = np.linspace(-0.03, 0.03, 20)
        R, Z = np.meshgrid(r_grid, z_grid)

        # Calculate fields at t=0, phi=0
        t = 0
        phi = 0

        Br_grid = np.zeros_like(R)
        Bz_grid = np.zeros_like(R)
        Er_grid = np.zeros_like(R)
        Ephi_grid = np.zeros_like(R)

        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                r, z = R[i,j], Z[i,j]
                Br, Bphi, Bz = self.physics.magnetic_field(r, phi, z, t)
                Er, Ephi, Ez = self.physics.electric_field(r, phi, z, t)

                Br_grid[i,j] = Br
                Bz_grid[i,j] = Bz
                Er_grid[i,j] = Er
                Ephi_grid[i,j] = Ephi

        # Plot magnetic field components
        im1 = ax1.contourf(R, Z, Br_grid, levels=20, cmap='RdBu_r')
        ax1.set_title('Radial Magnetic Field Br')
        ax1.set_xlabel('Radius (m)')
        ax1.set_ylabel('Z (m)')
        plt.colorbar(im1, ax=ax1, label='Tesla')

        im2 = ax2.contourf(R, Z, Bz_grid, levels=20, cmap='viridis')
        ax2.set_title('Axial Magnetic Field Bz')
        ax2.set_xlabel('Radius (m)')
        ax2.set_ylabel('Z (m)')
        plt.colorbar(im2, ax=ax2, label='Tesla')

        # Plot electric field components
        im3 = ax3.contourf(R, Z, Er_grid, levels=20, cmap='plasma')
        ax3.set_title('Radial Electric Field Er')
        ax3.set_xlabel('Radius (m)')
        ax3.set_ylabel('Z (m)')
        plt.colorbar(im3, ax=ax3, label='V/m')

        im4 = ax4.contourf(R, Z, Ephi_grid, levels=20, cmap='plasma')
        ax4.set_title('Azimuthal Electric Field Eφ')
        ax4.set_xlabel('Radius (m)')
        ax4.set_ylabel('Z (m)')
        plt.colorbar(im4, ax=ax4, label='V/m')

        # Draw gap boundaries on all plots
        for ax in [ax1, ax2, ax3, ax4]:
            ax.axvline(x=self.config.r_inner, color='white', linestyle='--', alpha=0.7)
            ax.axvline(x=self.config.r_outer, color='white', linestyle='--', alpha=0.7)
            ax.axhline(y=self.config.z_gap/2, color='white', linestyle='--', alpha=0.7)
            ax.axhline(y=-self.config.z_gap/2, color='white', linestyle='--', alpha=0.7)

        plt.tight_layout()
        return fig

    def analyze_accumulation(self, trajectories, particles_info):
        """
        Analyze where particles accumulate and print statistics.

        Args:
            trajectories: List of trajectory arrays
            particles_info: List of particle property dicts

        Returns:
            tuple: (final_positions, charges) as numpy arrays
        """
        final_positions = []
        charges = []

        for traj, info in zip(trajectories, particles_info):
            if len(traj) > 0:
                final_pos = traj[-1, :3]  # Final x, y, z
                final_positions.append(final_pos)
                charges.append(info['charge'])

        final_positions = np.array(final_positions)
        charges = np.array(charges)

        # Convert to cylindrical coordinates for analysis
        r_final = np.sqrt(final_positions[:, 0]**2 + final_positions[:, 1]**2)
        z_final = final_positions[:, 2]

        # Separate by charge sign
        positive_mask = charges > 0
        negative_mask = charges < 0

        print("\n=== PARTICLE ACCUMULATION ANALYSIS ===")
        print(f"Total particles: {len(trajectories)}")
        print(f"Positive charges: {np.sum(positive_mask)}")
        print(f"Negative charges: {np.sum(negative_mask)}")

        if np.sum(positive_mask) > 0:
            r_pos_mean = np.mean(r_final[positive_mask])
            r_pos_std = np.std(r_final[positive_mask])
            z_pos_mean = np.mean(z_final[positive_mask])
            print(f"Positive charges - Mean radius: {r_pos_mean:.4f}m ± {r_pos_std:.4f}m")
            print(f"Positive charges - Mean z: {z_pos_mean:.4f}m")

        if np.sum(negative_mask) > 0:
            r_neg_mean = np.mean(r_final[negative_mask])
            r_neg_std = np.std(r_final[negative_mask])
            z_neg_mean = np.mean(z_final[negative_mask])
            print(f"Negative charges - Mean radius: {r_neg_mean:.4f}m ± {r_neg_std:.4f}m")
            print(f"Negative charges - Mean z: {z_neg_mean:.4f}m")

        return final_positions, charges

    def save_figure(self, fig, filename):
        """
        Save figure to file.

        Args:
            fig: matplotlib Figure object
            filename: Output filename
        """
        try:
            fig.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Saved plot: {filename}")
        except Exception as e:
            print(f"Warning: Failed to save {filename}: {e}")
