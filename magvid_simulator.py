#!/usr/bin/env python3
"""
MAGVID Electromagnetic Particle Simulator
==========================================

Simulates charged particle trajectories in the MAGVID configuration:
- Rotating magnetic field from 4-coil assembly (creates radial electric field)
- Static axial magnetic field for confinement
- Visualizes particle accumulation and trajectories

Based on the analysis from "GroomLake Colonel Reveals All" document
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.integrate import ode
import time

class MAGVIDSimulator:
    def __init__(self):
        # Physical constants
        self.q_e = 1.602e-19  # Elementary charge [C]
        self.m_e = 9.109e-31  # Electron mass [kg]
        self.m_spore = 7.5e-14  # Lycopodium spore mass [kg]
        
        # Field parameters
        self.B_rotating = 0.75  # Tesla, rotating field strength
        self.B_static = 0.1     # Tesla, axial static field
        self.omega = 2 * np.pi * 1000  # rad/s, rotation frequency (1 kHz)
        
        # Geometry (in meters)
        self.r_inner = 0.09  # Inner radius of C-magnet gap
        self.r_outer = 0.11  # Outer radius of C-magnet gap
        self.z_gap = 0.02    # Gap height (2cm)
        
        # Simulation parameters
        self.dt = 1e-6  # Time step [s]
        self.t_max = 0.01  # Total simulation time [s]
        
    def magnetic_field(self, r, phi, z, t):
        """
        Calculate magnetic field at position (r, phi, z) and time t
        Returns: Br, Bphi, Bz components
        """
        # Rotating magnetic field in the gap region
        if self.r_inner <= r <= self.r_outer and abs(z) <= self.z_gap/2:
            # Four-coil rotating field creates primarily radial and azimuthal components
            Br = self.B_rotating * np.cos(self.omega * t - phi)
            Bphi = self.B_rotating * np.sin(self.omega * t - phi)
            Bz = 0
        else:
            # Fringing fields
            decay_factor = np.exp(-abs(r - (self.r_inner + self.r_outer)/2) / 0.02)
            Br = 0.1 * self.B_rotating * decay_factor * np.cos(self.omega * t - phi)
            Bphi = 0.1 * self.B_rotating * decay_factor * np.sin(self.omega * t - phi)
            Bz = 0
            
        # Add static axial field for confinement
        if abs(z) <= self.z_gap:
            Bz += self.B_static
        else:
            # Field falls off outside gap
            Bz += self.B_static * np.exp(-abs(z - self.z_gap/2) / 0.01)
            
        return Br, Bphi, Bz
    
    def electric_field(self, r, phi, z, t):
        """
        Calculate induced electric field from rotating magnetic field
        Uses Faraday's law: E = -dA/dt × B
        """
        # Induced electric field is primarily azimuthal
        if self.r_inner <= r <= self.r_outer and abs(z) <= self.z_gap/2:
            # E_phi = -dBr/dt * r/2 (simplified)
            Er = 0
            Ephi = -self.omega * self.B_rotating * r * np.sin(self.omega * t - phi) / 2
            Ez = 0
        else:
            # Fringing electric fields
            decay_factor = np.exp(-abs(r - (self.r_inner + self.r_outer)/2) / 0.02)
            Er = 0
            Ephi = -0.1 * self.omega * self.B_rotating * r * decay_factor * np.sin(self.omega * t - phi) / 2
            Ez = 0
            
        return Er, Ephi, Ez
    
    def lorentz_force(self, pos, vel, t, charge, mass):
        """
        Calculate Lorentz force: F = q(E + v × B)
        pos: [x, y, z] in Cartesian coordinates
        vel: [vx, vy, vz] in Cartesian coordinates
        """
        x, y, z = pos
        vx, vy, vz = vel
        
        # Convert to cylindrical coordinates
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        
        # Get fields in cylindrical coordinates
        Br, Bphi, Bz = self.magnetic_field(r, phi, z, t)
        Er, Ephi, Ez = self.electric_field(r, phi, z, t)
        
        # Convert B field to Cartesian
        Bx = Br * np.cos(phi) - Bphi * np.sin(phi)
        By = Br * np.sin(phi) + Bphi * np.cos(phi)
        
        # Convert E field to Cartesian
        Ex = Er * np.cos(phi) - Ephi * np.sin(phi)
        Ey = Er * np.sin(phi) + Ephi * np.cos(phi)
        
        # Calculate v × B
        cross_x = vy * Bz - vz * By
        cross_y = vz * Bx - vx * Bz
        cross_z = vx * By - vy * Bx
        
        # Total force
        Fx = charge * (Ex + cross_x)
        Fy = charge * (Ey + cross_y)
        Fz = charge * (Ez + cross_z)
        
        # Add gravity
        Fz -= mass * 9.81
        
        return np.array([Fx, Fy, Fz])
    
    def equations_of_motion(self, t, y, charge, mass):
        """
        ODE system for particle motion
        y = [x, y, z, vx, vy, vz]
        """
        pos = y[:3]
        vel = y[3:]
        
        force = self.lorentz_force(pos, vel, t, charge, mass)
        accel = force / mass
        
        return np.concatenate([vel, accel])
    
    def simulate_particle(self, initial_pos, initial_vel, charge, mass, color='blue'):
        """
        Simulate single particle trajectory
        """
        # Initial conditions
        y0 = np.concatenate([initial_pos, initial_vel])
        
        # Set up ODE solver
        r = ode(self.equations_of_motion).set_integrator('dopri5', rtol=1e-8, atol=1e-10)
        r.set_initial_value(y0, 0).set_f_params(charge, mass)
        
        # Storage for trajectory
        t_points = []
        trajectory = []
        
        # Integrate
        while r.successful() and r.t < self.t_max:
            r.integrate(r.t + self.dt)
            t_points.append(r.t)
            trajectory.append(r.y.copy())
            
        return np.array(t_points), np.array(trajectory)
    
    def run_multi_particle_simulation(self, n_particles=10):
        """
        Run simulation with multiple particles
        """
        trajectories = []
        particles_info = []
        
        # Generate random initial conditions
        np.random.seed(42)  # For reproducibility
        
        for i in range(n_particles):
            # Random initial positions around the gap region
            r0 = np.random.uniform(0.08, 0.12)
            phi0 = np.random.uniform(0, 2*np.pi)
            z0 = np.random.uniform(-self.z_gap/2, self.z_gap/2)
            
            x0 = r0 * np.cos(phi0)
            y0 = r0 * np.sin(phi0)
            initial_pos = np.array([x0, y0, z0])
            
            # Small random initial velocities
            initial_vel = np.random.normal(0, 10, 3)  # m/s
            
            # Random charge (electrons and positive ions)
            if np.random.random() > 0.5:
                charge = -1000 * self.q_e  # Negative
                color = 'red'
            else:
                charge = 1000 * self.q_e   # Positive
                color = 'blue'
            
            mass = self.m_spore
            
            particles_info.append({
                'charge': charge,
                'mass': mass,
                'color': color
            })
            
            # Simulate this particle
            t_traj, traj = self.simulate_particle(initial_pos, initial_vel, charge, mass)
            trajectories.append(traj)
            
            print(f"Particle {i+1}/{n_particles} simulated")
        
        return trajectories, particles_info
    
    def visualize_3d_trajectories(self, trajectories, particles_info):
        """
        Create 3D visualization of particle trajectories
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
        
        # Draw geometry
        # Gap region boundaries (cylindrical)
        theta = np.linspace(0, 2*np.pi, 50)
        
        # Inner cylinder
        x_inner = self.r_inner * np.cos(theta)
        y_inner = self.r_inner * np.sin(theta)
        z_inner_top = np.full_like(x_inner, self.z_gap/2)
        z_inner_bot = np.full_like(x_inner, -self.z_gap/2)
        
        ax.plot(x_inner, y_inner, z_inner_top, 'k-', alpha=0.5)
        ax.plot(x_inner, y_inner, z_inner_bot, 'k-', alpha=0.5)
        
        # Outer cylinder
        x_outer = self.r_outer * np.cos(theta)
        y_outer = self.r_outer * np.sin(theta)
        z_outer_top = np.full_like(x_outer, self.z_gap/2)
        z_outer_bot = np.full_like(x_outer, -self.z_gap/2)
        
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
    
    def analyze_accumulation(self, trajectories, particles_info):
        """
        Analyze where particles accumulate
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
    
    def plot_field_visualization(self):
        """
        Visualize the magnetic and electric field structure
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
                Br, Bphi, Bz = self.magnetic_field(r, phi, z, t)
                Er, Ephi, Ez = self.electric_field(r, phi, z, t)
                
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
            ax.axvline(x=self.r_inner, color='white', linestyle='--', alpha=0.7)
            ax.axvline(x=self.r_outer, color='white', linestyle='--', alpha=0.7)
            ax.axhline(y=self.z_gap/2, color='white', linestyle='--', alpha=0.7)
            ax.axhline(y=-self.z_gap/2, color='white', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        return fig


def main():
    """
    Run the MAGVID simulation
    """
    print("MAGVID Electromagnetic Particle Simulator")
    print("=========================================")
    
    # Create simulator
    sim = MAGVIDSimulator()
    
    # Display parameters
    print(f"Simulation parameters:")
    print(f"  Rotating B field: {sim.B_rotating} T at {sim.omega/(2*np.pi):.1f} Hz")
    print(f"  Static B field: {sim.B_static} T")
    print(f"  Gap region: {sim.r_inner*100:.1f} - {sim.r_outer*100:.1f} cm radius")
    print(f"  Gap height: {sim.z_gap*100:.1f} cm")
    print(f"  Simulation time: {sim.t_max*1000:.1f} ms")
    print()
    
    # Run multi-particle simulation
    print("Running multi-particle simulation...")
    start_time = time.time()
    trajectories, particles_info = sim.run_multi_particle_simulation(n_particles=20)
    end_time = time.time()
    print(f"Simulation completed in {end_time - start_time:.2f} seconds")
    print()
    
    # Analyze results
    final_positions, charges = sim.analyze_accumulation(trajectories, particles_info)
    
    # Create visualizations
    print("Creating visualizations...")
    
    # 3D trajectory plot
    fig1 = sim.visualize_3d_trajectories(trajectories, particles_info)
    
    # Field visualization
    fig2 = sim.plot_field_visualization()
    
    # Show plots
    plt.show()
    
    print("\nSimulation complete!")
    print("Red trajectories: Negatively charged particles")
    print("Blue trajectories: Positively charged particles")
    print("Circles: Starting positions, Squares: Final positions")


if __name__ == "__main__":
    main()
