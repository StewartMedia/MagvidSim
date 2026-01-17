#!/usr/bin/env python3
"""
Configuration Management for MAGVID Simulator
==============================================

Handles loading, validation, and error checking of configuration parameters.
"""

import configparser
import os
import sys


class ConfigurationError(Exception):
    """Raised when configuration is invalid or missing."""
    pass


class MAGVIDConfig:
    """Manages configuration parameters for MAGVID simulation."""

    # Physical constants
    Q_E = 1.602176634e-19  # Elementary charge [C]
    M_E = 9.10938356e-31   # Electron mass [kg]
    M_SPORE = 7.5e-14      # Lycopodium spore mass [kg]
    GRAVITY = 9.81         # Gravitational acceleration [m/s^2]

    def __init__(self, config_file='config.ini'):
        """
        Initialize configuration from file.

        Args:
            config_file: Path to configuration file

        Raises:
            ConfigurationError: If config file is missing or invalid
        """
        self.config_file = config_file
        self.config = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))

        # Load configuration
        self._load_config()

        # Parse and validate parameters
        self._parse_simulation_params()
        self._parse_field_params()
        self._parse_geometry_params()
        self._parse_particle_params()
        self._parse_output_params()

        # Validate all parameters
        self._validate_config()

    def _load_config(self):
        """Load configuration file with error handling."""
        if not os.path.exists(self.config_file):
            raise ConfigurationError(
                f"Configuration file '{self.config_file}' not found.\n"
                f"Please create a config.ini file or specify a valid path."
            )

        try:
            files_read = self.config.read(self.config_file)
            if not files_read:
                raise ConfigurationError(
                    f"Failed to read configuration file '{self.config_file}'"
                )
        except configparser.Error as e:
            raise ConfigurationError(
                f"Error parsing configuration file: {e}"
            )

    def _parse_simulation_params(self):
        """Parse simulation parameters."""
        try:
            self.dt = self.config.getfloat('simulation', 'dt')
            self.t_max = self.config.getfloat('simulation', 't_max')
            self.n_particles = self.config.getint('simulation', 'n_particles')
        except (configparser.NoSectionError, configparser.NoOptionError) as e:
            raise ConfigurationError(
                f"Missing simulation parameter: {e}"
            )
        except ValueError as e:
            raise ConfigurationError(
                f"Invalid simulation parameter value: {e}"
            )

    def _parse_field_params(self):
        """Parse field parameters."""
        try:
            self.B_rotating = self.config.getfloat('fields', 'B_rotating')
            self.B_static = self.config.getfloat('fields', 'B_static')
            frequency = self.config.getfloat('fields', 'frequency')
            self.omega = 2 * 3.14159265359 * frequency  # Convert Hz to rad/s
        except (configparser.NoSectionError, configparser.NoOptionError) as e:
            raise ConfigurationError(
                f"Missing field parameter: {e}"
            )
        except ValueError as e:
            raise ConfigurationError(
                f"Invalid field parameter value: {e}"
            )

    def _parse_geometry_params(self):
        """Parse geometry parameters."""
        try:
            self.r_inner = self.config.getfloat('geometry', 'r_inner')
            self.r_outer = self.config.getfloat('geometry', 'r_outer')
            self.z_gap = self.config.getfloat('geometry', 'z_gap')
        except (configparser.NoSectionError, configparser.NoOptionError) as e:
            raise ConfigurationError(
                f"Missing geometry parameter: {e}"
            )
        except ValueError as e:
            raise ConfigurationError(
                f"Invalid geometry parameter value: {e}"
            )

    def _parse_particle_params(self):
        """Parse particle parameters."""
        try:
            self.charge_multiplier = self.config.getint('particles', 'charge_multiplier')
            self.mass_spore = self.config.getfloat('particles', 'mass_spore')
            self.charge_ratio = self.config.getfloat('particles', 'charge_ratio')
        except (configparser.NoSectionError, configparser.NoOptionError) as e:
            raise ConfigurationError(
                f"Missing particle parameter: {e}"
            )
        except ValueError as e:
            raise ConfigurationError(
                f"Invalid particle parameter value: {e}"
            )

    def _parse_output_params(self):
        """Parse output parameters."""
        try:
            self.save_plots = self.config.getboolean('output', 'save_plots')
            self.plot_3d = self.config.getboolean('output', 'plot_3d')
            self.plot_fields = self.config.getboolean('output', 'plot_fields')
            self.plot_trajectories = self.config.getboolean('output', 'plot_trajectories')
            self.animation = self.config.getboolean('output', 'animation')
        except (configparser.NoSectionError, configparser.NoOptionError) as e:
            raise ConfigurationError(
                f"Missing output parameter: {e}"
            )
        except ValueError as e:
            raise ConfigurationError(
                f"Invalid output parameter value: {e}"
            )

    def _validate_config(self):
        """Validate all configuration parameters."""
        errors = []

        # Validate simulation parameters
        if self.dt <= 0:
            errors.append("Time step (dt) must be positive")
        if self.t_max <= 0:
            errors.append("Simulation time (t_max) must be positive")
        if self.dt >= self.t_max:
            errors.append("Time step (dt) must be smaller than total time (t_max)")
        if self.n_particles <= 0:
            errors.append("Number of particles must be positive")
        if self.n_particles > 10000:
            errors.append("Number of particles exceeds reasonable limit (10000)")

        # Validate field parameters
        if self.B_rotating < 0:
            errors.append("Rotating magnetic field strength must be non-negative")
        if self.B_static < 0:
            errors.append("Static magnetic field strength must be non-negative")
        if self.omega <= 0:
            errors.append("Rotation frequency must be positive")
        if self.B_rotating > 10:
            errors.append("Rotating field strength exceeds typical laboratory capabilities (>10 T)")

        # Validate geometry parameters
        if self.r_inner <= 0:
            errors.append("Inner radius must be positive")
        if self.r_outer <= 0:
            errors.append("Outer radius must be positive")
        if self.r_inner >= self.r_outer:
            errors.append("Inner radius must be less than outer radius")
        if self.z_gap <= 0:
            errors.append("Gap height must be positive")

        # Validate particle parameters
        if self.charge_multiplier <= 0:
            errors.append("Charge multiplier must be positive")
        if self.mass_spore <= 0:
            errors.append("Particle mass must be positive")
        if not 0 <= self.charge_ratio <= 1:
            errors.append("Charge ratio must be between 0 and 1")

        # Raise error if any validation failed
        if errors:
            error_message = "Configuration validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
            raise ConfigurationError(error_message)

    def display_parameters(self):
        """Display configuration parameters in a formatted way."""
        print("MAGVID Simulator Configuration")
        print("=" * 50)
        print("\n[Simulation Parameters]")
        print(f"  Time step:          {self.dt:.2e} s")
        print(f"  Total time:         {self.t_max*1000:.2f} ms")
        print(f"  Number of particles: {self.n_particles}")

        print("\n[Field Parameters]")
        print(f"  Rotating B field:   {self.B_rotating:.3f} T at {self.omega/(2*3.14159265359):.1f} Hz")
        print(f"  Static B field:     {self.B_static:.3f} T")

        print("\n[Geometry]")
        print(f"  Inner radius:       {self.r_inner*100:.2f} cm")
        print(f"  Outer radius:       {self.r_outer*100:.2f} cm")
        print(f"  Gap height:         {self.z_gap*100:.2f} cm")

        print("\n[Particles]")
        print(f"  Charge multiplier:  {self.charge_multiplier}× elementary charge")
        print(f"  Particle mass:      {self.mass_spore:.2e} kg")
        print(f"  Positive charge %:  {self.charge_ratio*100:.1f}%")

        print("\n[Output Options]")
        print(f"  Save plots:         {self.save_plots}")
        print(f"  3D trajectories:    {self.plot_3d}")
        print(f"  Field visualization: {self.plot_fields}")
        print(f"  Animation:          {self.animation}")
        print("=" * 50)
