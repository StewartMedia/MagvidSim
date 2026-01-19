#!/usr/bin/env python3
"""
Rotating Magnetic Field Frame Cache
====================================

Pre-computes and caches magnetic field frames for a rotating excitation pattern.
Used for efficient particle trajectory simulation with time-varying fields.

The MAGVID uses 4 windings in quadrature (90° phase shifts) to create a
unipolar rotating magnetic field, NOT like a standard 2-phase AC motor stator.
"""

import numpy as np
import pickle
from pathlib import Path
from typing import List, Tuple
from dataclasses import dataclass
from em_solver import STEPLoader, MagneticCore, CurrentSheet, MagnetostaticSolver


@dataclass
class RotatingFieldCache:
    """Cached rotating magnetic field data."""
    frames: List[Tuple[np.ndarray, np.ndarray]]  # [(points, B_vectors), ...]
    n_frames: int
    rotation_period: float  # seconds
    bounds: np.ndarray
    resolution: int

    def get_field_at_time(self, t: float) -> Tuple[np.ndarray, np.ndarray]:
        """Get interpolated field at time t."""
        # Calculate frame index with wrapping
        frame_idx = int((t / self.rotation_period) * self.n_frames) % self.n_frames
        return self.frames[frame_idx]

    def save(self, filepath: str):
        """Save cache to disk."""
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
        print(f"Saved rotating field cache: {filepath}")
        file_size_mb = Path(filepath).stat().st_size / (1024 * 1024)
        print(f"  File size: {file_size_mb:.1f} MB")

    @classmethod
    def load(cls, filepath: str):
        """Load cache from disk."""
        with open(filepath, 'rb') as f:
            cache = pickle.load(f)
        print(f"Loaded rotating field cache: {filepath}")
        print(f"  Frames: {cache.n_frames}")
        print(f"  Period: {cache.rotation_period} s")
        return cache


def generate_rotating_field_cache(
    step_file: str,
    winding_names: List[str],
    winding_amplitudes: List[float],
    winding_phases: List[float],  # radians
    core_mu_r: float,
    rotation_frequency: float,  # Hz
    n_frames: int = 128,
    field_resolution: int = 8,
    margin: float = 0.05
) -> RotatingFieldCache:
    """
    Generate cached magnetic field frames for rotating excitation.

    Args:
        step_file: Path to STEP file with core geometry
        winding_names: List of winding layer names (e.g., ["Winding0", "Winding1", ...])
        winding_amplitudes: Peak currents for each winding [A]
        winding_phases: Phase offsets for each winding [radians]
                       For quadrature: [0, π/2, π, 3π/2]
        core_mu_r: Relative permeability of core
        rotation_frequency: Rotation frequency [Hz]
        n_frames: Number of frames to cache (default 128)
        field_resolution: Grid resolution for field calculation
        margin: Spatial margin around core [m]

    Returns:
        RotatingFieldCache with n_frames pre-computed field distributions
    """
    print("="*70)
    print(f"Generating Rotating Magnetic Field Cache")
    print("="*70)
    print(f"  Rotation frequency: {rotation_frequency} Hz")
    print(f"  Period: {1/rotation_frequency:.4f} s")
    print(f"  Frames: {n_frames}")
    print(f"  Resolution: {field_resolution}³ = {field_resolution**3} points")
    print(f"  Windings: {len(winding_names)}")

    # Load geometry
    print("\n" + "-"*70)
    loader = STEPLoader(step_file)
    if not loader.load():
        raise RuntimeError("Failed to load STEP file")

    vertices, faces = loader.tessellate_shape()
    if len(vertices) == 0:
        raise RuntimeError("No geometry extracted")

    print(f"  Extracted mesh: {len(vertices)} vertices, {len(faces)} faces")

    # Subsample for performance
    if len(faces) > 500:
        target_faces = min(300, len(faces))
        step = max(1, len(faces) // target_faces)
        faces = faces[::step]
        print(f"  Subsampled to: {len(faces)} faces")

    # Create core
    core = MagneticCore(vertices=vertices, faces=faces, mu_r=core_mu_r)

    # Calculate field bounds
    mins = vertices.min(axis=0)
    maxs = vertices.max(axis=0)
    bounds = np.array([[mins[0]-margin, maxs[0]+margin],
                      [mins[1]-margin, maxs[1]+margin],
                      [mins[2]-margin, maxs[2]+margin]])

    # Time points for each frame
    rotation_period = 1.0 / rotation_frequency
    times = np.linspace(0, rotation_period, n_frames, endpoint=False)

    # Generate frames
    print("\n" + "-"*70)
    print("Generating field frames...")
    frames = []

    for frame_idx, t in enumerate(times):
        print(f"\nFrame {frame_idx+1}/{n_frames} (t={t:.6f}s):")

        # Calculate winding currents at this time
        winding_currents = {}
        for name, amplitude, phase in zip(winding_names, winding_amplitudes, winding_phases):
            # I(t) = I_peak * cos(2πft + φ)
            current = amplitude * np.cos(2 * np.pi * rotation_frequency * t + phase)
            winding_currents[name] = current
            print(f"  {name}: {current:+.1f} A")

        # Create current sheets (simplified - use core faces as winding approximation)
        windings = []
        for name, current in winding_currents.items():
            if abs(current) < 1.0:  # Skip if current is negligible
                continue

            winding = CurrentSheet(
                name=name,
                vertices=vertices,
                faces=faces[:len(faces)//len(winding_names)],  # Distribute faces
                current=current
            )
            windings.append(winding)

        # Calculate field
        solver = MagnetostaticSolver(cores=[core], windings=windings)
        points, B_vectors = solver.calculate_grid(bounds, resolution=field_resolution)

        frames.append((points.copy(), B_vectors.copy()))

        # Statistics
        B_mag = np.linalg.norm(B_vectors, axis=1)
        valid_B = B_mag[~np.isnan(B_mag)]
        if len(valid_B) > 0:
            print(f"  |B| range: [{valid_B.min():.6e}, {valid_B.max():.6e}] T")

    print("\n" + "="*70)
    print("Field cache generation complete!")

    return RotatingFieldCache(
        frames=frames,
        n_frames=n_frames,
        rotation_period=rotation_period,
        bounds=bounds,
        resolution=field_resolution
    )


def main():
    """Example: Generate 4-winding quadrature field cache for MAGVID."""

    # Configuration for 4-winding quadrature system
    config = {
        'step_file': 'out.STEP',
        'winding_names': ['Winding0', 'Winding1', 'Winding2', 'Winding3'],
        'winding_amplitudes': [500.0, 500.0, 500.0, 500.0],  # All same amplitude
        'winding_phases': [0, np.pi/2, np.pi, 3*np.pi/2],   # Quadrature: 90° apart
        'core_mu_r': 5000.0,
        'rotation_frequency': 60.0,  # 60 Hz
        'n_frames': 128,
        'field_resolution': 8,  # 8³ = 512 points per frame
    }

    # Generate cache
    cache = generate_rotating_field_cache(**config)

    # Save to disk
    cache.save('rotating_field_cache_60Hz_128frames.pkl')

    print("\nCache ready for particle simulation!")
    print("Usage in simulator:")
    print("  cache = RotatingFieldCache.load('rotating_field_cache_60Hz_128frames.pkl')")
    print("  points, B = cache.get_field_at_time(t)")


if __name__ == "__main__":
    main()
