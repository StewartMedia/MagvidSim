#!/usr/bin/env python3
"""
Physics Engine for MAGVID Simulator
====================================

Handles electromagnetic field calculations and particle dynamics.
"""

import numpy as np
from scipy.integrate import ode


class SimulationError(Exception):
    """Raised when simulation encounters an error."""
    pass


class MAGVIDPhysics:
    """
    Physics engine for MAGVID electromagnetic simulation.

    Calculates electromagnetic fields and particle trajectories using
    the Lorentz force law in the MAGVID configuration.
    """

    def __init__(self, config):
        """
        Initialize physics engine with configuration.

        Args:
            config: MAGVIDConfig object with simulation parameters
        """
        self.config = config

    def magnetic_field(self, r, phi, z, t):
        """
        Calculate magnetic field at position (r, phi, z) and time t.

        The field consists of:
        - Rotating field from 4-coil assembly (creates radial/azimuthal components)
        - Static axial field for confinement

        Args:
            r: Radial coordinate [m]
            phi: Azimuthal angle [rad]
            z: Axial coordinate [m]
            t: Time [s]

        Returns:
            tuple: (Br, Bphi, Bz) magnetic field components in Tesla
        """
        # Rotating magnetic field in the gap region
        if self.config.r_inner <= r <= self.config.r_outer and abs(z) <= self.config.z_gap/2:
            # Four-coil rotating field creates primarily radial and azimuthal components
            Br = self.config.B_rotating * np.cos(self.config.omega * t - phi)
            Bphi = self.config.B_rotating * np.sin(self.config.omega * t - phi)
            Bz = 0
        else:
            # Fringing fields outside gap
            r_center = (self.config.r_inner + self.config.r_outer) / 2
            decay_factor = np.exp(-abs(r - r_center) / 0.02)
            Br = 0.1 * self.config.B_rotating * decay_factor * np.cos(self.config.omega * t - phi)
            Bphi = 0.1 * self.config.B_rotating * decay_factor * np.sin(self.config.omega * t - phi)
            Bz = 0

        # Add static axial field for confinement
        if abs(z) <= self.config.z_gap:
            Bz += self.config.B_static
        else:
            # Field falls off outside gap
            Bz += self.config.B_static * np.exp(-abs(z - self.config.z_gap/2) / 0.01)

        return Br, Bphi, Bz

    def electric_field(self, r, phi, z, t):
        """
        Calculate induced electric field from rotating magnetic field.

        Uses Faraday's law: ∇ × E = -∂B/∂t
        The rotating magnetic field induces an azimuthal electric field.

        Args:
            r: Radial coordinate [m]
            phi: Azimuthal angle [rad]
            z: Axial coordinate [m]
            t: Time [s]

        Returns:
            tuple: (Er, Ephi, Ez) electric field components in V/m
        """
        # Induced electric field is primarily azimuthal
        if self.config.r_inner <= r <= self.config.r_outer and abs(z) <= self.config.z_gap/2:
            # E_phi = -dBr/dt * r/2 (simplified from Faraday's law)
            Er = 0
            Ephi = -self.config.omega * self.config.B_rotating * r * np.sin(self.config.omega * t - phi) / 2
            Ez = 0
        else:
            # Fringing electric fields
            r_center = (self.config.r_inner + self.config.r_outer) / 2
            decay_factor = np.exp(-abs(r - r_center) / 0.02)
            Er = 0
            Ephi = -0.1 * self.config.omega * self.config.B_rotating * r * decay_factor * np.sin(self.config.omega * t - phi) / 2
            Ez = 0

        return Er, Ephi, Ez

    def lorentz_force(self, pos, vel, t, charge, mass):
        """
        Calculate Lorentz force on a charged particle.

        F = q(E + v × B) + m*g (includes gravity)

        Args:
            pos: Position [x, y, z] in Cartesian coordinates [m]
            vel: Velocity [vx, vy, vz] in Cartesian coordinates [m/s]
            t: Time [s]
            charge: Particle charge [C]
            mass: Particle mass [kg]

        Returns:
            np.array: Force [Fx, Fy, Fz] in Newtons
        """
        x, y, z = pos
        vx, vy, vz = vel

        # Convert to cylindrical coordinates
        r = np.sqrt(x**2 + y**2)
        if r < 1e-10:  # Avoid division by zero at origin
            r = 1e-10
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

        # Total force: Lorentz force + gravity
        Fx = charge * (Ex + cross_x)
        Fy = charge * (Ey + cross_y)
        Fz = charge * (Ez + cross_z) - mass * self.config.GRAVITY

        return np.array([Fx, Fy, Fz])

    def equations_of_motion(self, t, y, charge, mass):
        """
        ODE system for particle motion.

        Defines the system: dy/dt = [velocity, acceleration]

        Args:
            t: Time [s]
            y: State vector [x, y, z, vx, vy, vz]
            charge: Particle charge [C]
            mass: Particle mass [kg]

        Returns:
            np.array: Derivative [vx, vy, vz, ax, ay, az]
        """
        pos = y[:3]
        vel = y[3:]

        force = self.lorentz_force(pos, vel, t, charge, mass)
        accel = force / mass

        return np.concatenate([vel, accel])

    def simulate_particle(self, initial_pos, initial_vel, charge, mass):
        """
        Simulate single particle trajectory.

        Uses Dormand-Prince 5th order Runge-Kutta method (dopri5)
        with adaptive step size.

        Args:
            initial_pos: Initial position [x, y, z] [m]
            initial_vel: Initial velocity [vx, vy, vz] [m/s]
            charge: Particle charge [C]
            mass: Particle mass [kg]

        Returns:
            tuple: (t_points, trajectory) where
                t_points: Array of time points [s]
                trajectory: Array of states [x, y, z, vx, vy, vz]

        Raises:
            SimulationError: If integration fails
        """
        # Initial conditions
        y0 = np.concatenate([initial_pos, initial_vel])

        # Set up ODE solver
        try:
            r = ode(self.equations_of_motion).set_integrator('dopri5', rtol=1e-8, atol=1e-10)
            r.set_initial_value(y0, 0).set_f_params(charge, mass)
        except Exception as e:
            raise SimulationError(f"Failed to initialize ODE solver: {e}")

        # Storage for trajectory
        t_points = []
        trajectory = []

        # Integrate
        try:
            while r.successful() and r.t < self.config.t_max:
                r.integrate(r.t + self.config.dt)
                t_points.append(r.t)
                trajectory.append(r.y.copy())

                # Check for numerical overflow
                if not np.all(np.isfinite(r.y)):
                    raise SimulationError(
                        f"Integration diverged at t={r.t:.2e}s. "
                        "Particle trajectory became non-finite. "
                        "Try reducing time step or field strengths."
                    )

        except Exception as e:
            if isinstance(e, SimulationError):
                raise
            raise SimulationError(f"Integration failed: {e}")

        if len(trajectory) == 0:
            raise SimulationError("Integration produced no trajectory points")

        return np.array(t_points), np.array(trajectory)

    def run_multi_particle_simulation(self, progress_callback=None):
        """
        Run simulation with multiple particles.

        Generates random initial conditions and simulates each particle.

        Args:
            progress_callback: Optional function(i, n) called after each particle

        Returns:
            tuple: (trajectories, particles_info) where
                trajectories: List of trajectory arrays
                particles_info: List of dicts with particle properties

        Raises:
            SimulationError: If simulation fails
        """
        trajectories = []
        particles_info = []

        # Generate random initial conditions
        np.random.seed(42)  # For reproducibility

        for i in range(self.config.n_particles):
            try:
                # Random initial positions around the gap region
                r0 = np.random.uniform(0.08, 0.12)
                phi0 = np.random.uniform(0, 2*np.pi)
                z0 = np.random.uniform(-self.config.z_gap/2, self.config.z_gap/2)

                x0 = r0 * np.cos(phi0)
                y0 = r0 * np.sin(phi0)
                initial_pos = np.array([x0, y0, z0])

                # Small random initial velocities
                initial_vel = np.random.normal(0, 10, 3)  # m/s

                # Random charge based on charge_ratio
                if np.random.random() > self.config.charge_ratio:
                    charge = -self.config.charge_multiplier * self.config.Q_E  # Negative
                    color = 'red'
                else:
                    charge = self.config.charge_multiplier * self.config.Q_E   # Positive
                    color = 'blue'

                mass = self.config.mass_spore

                particles_info.append({
                    'charge': charge,
                    'mass': mass,
                    'color': color
                })

                # Simulate this particle
                t_traj, traj = self.simulate_particle(initial_pos, initial_vel, charge, mass)
                trajectories.append(traj)

                # Progress callback
                if progress_callback:
                    progress_callback(i + 1, self.config.n_particles)

            except SimulationError as e:
                raise SimulationError(f"Failed simulating particle {i+1}: {e}")

        return trajectories, particles_info
