# MAGVID Electromagnetic Particle Simulator

## Overview

This simulator models the electromagnetic behavior of charged particles in the MAGVID (MAGnetic Vortex hyper-Ionization Device) configuration as described in the infamous "GroomLake Colonel Reveals All" USENET post from 1995.

## Theory

The MAGVID uses two magnetic field systems:
1. **Rotating magnetic field**: Created by 4 coils in quadrature (90° phase shifts) like a 2-phase AC motor
2. **Static axial field**: Central DC electromagnet providing vertical confinement

The interaction creates complex charged particle trajectories via the Lorentz force: **F = q(E + v × B)**

## Key Features

- **Particle trajectory simulation** using Runge-Kutta ODE integration
- **Multi-particle dynamics** with charge separation analysis
- **3D field visualization** showing magnetic and electric field distributions
- **Configurable parameters** via config.ini file

## Installation

Run the setup script:
```bash
cd ~/coding/magvidsim
./setup.sh
```

Or manually:
```bash
python3 -m venv magvid_env
source magvid_env/bin/activate
pip install -r requirements.txt
```

## Usage

### Basic Simulation
```bash
source magvid_env/bin/activate
python magvid_simulator.py
```

### Configuration
Edit `config.ini` to match your experimental parameters:
- Field strengths and frequencies
- Geometry dimensions  
- Particle properties
- Output options

## Simulation Results

The simulator outputs:
- **3D trajectory plots**: Particle paths colored by charge sign
- **Field visualizations**: Magnetic and electric field distributions
- **Accumulation analysis**: Charge separation statistics

## Key Findings

- **Charge separation**: Positive/negative particles accumulate at different radii
- **Stable confinement**: Particles remain trapped in gap region
- **Frequency dependence**: Higher rotation rates increase separation efficiency
- **Gap geometry critical**: Field strength depends on r_inner/r_outer ratio

## Experimental Validation

To correlate with physical build:
1. **Field mapping**: Compare simulated vs measured B-field
2. **Phase verification**: Confirm 90° coil relationships  
3. **Frequency response**: Test accumulation vs rotation rate
4. **Particle detection**: Use ionized air or charged dust

## Files Structure

- `magvid_simulator.py` - Main simulation code
- `config.ini` - Parameter configuration
- `requirements.txt` - Python dependencies
- `setup.sh` - Installation script
- `README.md` - This file

## References

- Mathias Båge: "The MAGVID-GLCRA" (2007)
- James Stephens: "GroomLake Colonel Reveals All" (1995)
- Lorentz Force Law and electromagnetic field theory