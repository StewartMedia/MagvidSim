# MAGVID Electromagnetic Simulation Suite

## Overview

This suite provides comprehensive electromagnetic simulation tools for the MAGVID (MAGnetic Vortex hyper-Ionization Device) configuration as described in the "GroomLake Colonel Reveals All" USENET post from 1995.

**Two complementary tools:**
1. **Particle Trajectory Simulator** (`simulator.py`) - Simulates charged particle motion in idealized fields
2. **Electromagnetic Field Solver** (`em_solver.py`) - Computes B-fields from actual CAD geometry (STEP files)

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
cd ~/magvidsim
./setup.sh
```

Or manually:
```bash
python3 -m venv magvid_env
source magvid_env/bin/activate
pip install -r requirements.txt
```

## Usage

### Particle Trajectory Simulator

Simulates charged particle motion in parametric magnetic fields:

```bash
source magvid_env/bin/activate
python simulator.py
```

For headless/server environments (no display):
```bash
python simulator.py --no-display
```

**Configuration:** Edit `config.ini` to set:
- Field strengths and frequencies
- Geometry dimensions
- Particle properties
- Output options

### Electromagnetic Field Solver

Computes magnetostatic fields from CAD geometry (STEP files):

```bash
source magvid_env/bin/activate
python em_solver.py
```

**Input Requirements:**
- STEP AP214 file with solid bodies (cores) and PRESENTATION_LAYER_ASSIGNMENT layers (windings)
- Winding currents specified in code (e.g., "Winding0 = -500A")
- Core relative permeability (μr)

**Output:**
- Interactive HTML file with 3D visualization
- Color-mapped |B| magnitude
- Magnetic field streamlines
- Full pan/tilt/rotate/zoom controls in web browser

**Example:** See `test_em_basic.py` for a quick validation

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

## Code Structure

### Particle Trajectory Simulator
- `simulator.py` - Main entry point and orchestration
- `config.py` - Configuration management with validation
- `physics.py` - Parametric field calculations and particle dynamics
- `visualization.py` - Matplotlib-based 3D plotting
- `config.ini` - Parameter configuration file

### Electromagnetic Field Solver
- `em_solver.py` - STEP file loading and magnetostatic field calculation
- `test_em_basic.py` - Component validation and testing

### Shared
- `requirements.txt` - Python dependencies (NumPy, SciPy, Matplotlib, PyVista, CADQuery)
- `setup.sh` - Installation script
- `.gitignore` - Version control exclusions

### Error Handling

The simulator includes comprehensive error handling:
- **Configuration validation**: Checks for missing files, invalid parameters, and physically unrealistic values
- **Simulation error detection**: Catches numerical overflow, integration failures, and divergence
- **Helpful error messages**: Clear descriptions of what went wrong and how to fix it

## References

- Mathias Båge: "The MAGVID-GLCRA" (2007)
- James Stephens: "GroomLake Colonel Reveals All" (1995)
- Lorentz Force Law and electromagnetic field theory
