#!/usr/bin/env python3
"""
Usage:
  python vibration_project_multi_staticbond.py GEOMETRY.xyz DIPOLE \
      --temperature 310 --timestep 0.24e-15 --nh_indices 985 986 987 988 --plot

Notes:
- GEOMETRY.xyz is used only to compute static bond unit vectors (one frame).
- DIPOLE should contain per-frame dipole vector columns (reads columns 4,5,6).
"""

import numpy as np
from scipy import signal, fftpack
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import argparse
import os

PLANCK = 1.05457180013E-34  # J s (kg m^2 s^-1)
BOLTZMANN = 1.38064852E-23  # J K^-1
SPEED_OF_LIGHT = 299792458  # m s^-1

def read_geometry(geom_file):
    """Read single-frame XYZ-like file and return Nx3 coords.
       Expects at least two header lines (standard xyz style)."""
    with open(geom_file) as f:
        lines = f.readlines()
    coords = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) >= 4:
            x, y, z = map(float, parts[1:4])
            coords.append([x, y, z])
    return np.array(coords)

def get_bond_unit_vector(coords, idx1, idx2):
    vec = coords[idx2] - coords[idx1]
    norm = np.linalg.norm(vec)
    if norm == 0.0:
        raise ValueError(f"Zero-length bond vector between atoms {idx1+1} and {idx2+1}")
    return vec / norm

def get_proj_time_series_from_static_bonds(dipole_file, bond_unit_vectors):
    """
    Read dipole file and project each dipole onto each static bond unit vector.
    Returns list of 1D arrays (projections), one per bond.
    """
    # read columns: we expect columns [time?, dx, dy, dz] or similar; adjust if needed
    data = np.loadtxt(dipole_file, usecols=(0, 4, 5, 6))
    dipoles = data[:, 1:4]  # shape (n_frames, 3)
    projections = []
    for bv in bond_unit_vectors:
        proj = np.dot(dipoles, bv)  # (n_frames,)
        projections.append(proj)
    return projections

def get_autocorr_from_proj(dipole_proj):
    """Compute (biased) autocorrelation using FFT convolution similar to original script."""
    n_frames = len(dipole_proj)
    # embed into a larger array (as in original script)
    padded_len = 2 * n_frames - (n_frames % 2)
    dipole_proj_shifted = np.zeros(padded_len)
    start = n_frames // 2
    dipole_proj_shifted[start:start + n_frames] = dipole_proj

    autocorr_full = signal.fftconvolve(dipole_proj_shifted, dipole_proj[::-1], mode='same')
    autocorr = autocorr_full[-n_frames:] / np.arange(n_frames, 0, -1)
    # remove any tiny imaginary component produced by numerical noise
    autocorr = np.real_if_close(autocorr)
    return autocorr

def calc_spectrum_from_autocorr(autocorr, temperature, timestep, wn_min=0, wn_max=4000):
    """
    Convert autocorrelation to frequency domain with a DCT (type 1),
    apply harmonic quantum correction factor and return wavenumbers (cm^-1)
    and normalized spectrum.
    """
    # DCT of the autocorrelation (type=1 replicates cosine transform behavior)
    lineshape = fftpack.dct(autocorr, type=1)[1:]  # drop the zero-frequency term
    n_points = len(lineshape)
    # frequency grid (Hz): Nyquist = 0.5 / timestep
    freq_Hz = np.linspace(0, 0.5 / timestep, n_points + 1)[1:]  # skip zero
    # convert Hz -> cm^-1 : wn_cm^-1 = freq_Hz / (c * 100)
    freq_wn = freq_Hz / (SPEED_OF_LIGHT * 100.0)

    # quantum correction factor for field (approx.)
    # keep factor in SI: PLANCK * freq (Hz) has units J, kB*T in J
    field_desc = freq_Hz * (1.0 - np.exp(-PLANCK * freq_Hz / (BOLTZMANN * temperature)))

    spectrum = lineshape * field_desc

    # smoothing
    smooth_spectrum = gaussian_filter(spectrum, 2.5)

    # normalization (avoid division by zero)
    peak_height = np.abs(np.max(smooth_spectrum) - np.min(smooth_spectrum))
    if peak_height == 0:
        norm_spectrum = smooth_spectrum
    else:
        norm_spectrum = smooth_spectrum / peak_height

    mask = (freq_wn >= wn_min) & (freq_wn <= wn_max)
    return freq_wn[mask], norm_spectrum[mask]

def save_spectrum(filename, wn, spectrum):
    with open(filename, 'w') as f:
        for w, s in zip(wn, spectrum):
            f.write(f"{w:8.2f} {s:12.6e}\n")

def main():
    parser = argparse.ArgumentParser(description="Vibrational spectrum projected on static N–H bond(s)")
    parser.add_argument('GEOMETRY', help='Path to GEOMETRY.xyz (static reference geometry)')
    parser.add_argument('DIPOLE', help='Path to DIPOLE file (per-frame dipole vectors)')
    parser.add_argument('--temperature', type=float, default=310.0, help='Temperature in K')
    parser.add_argument('--timestep', type=float, default=0.24e-15, help='Timestep in seconds (default 0.24 fs)')
    parser.add_argument('--nh_indices', type=int, nargs='+', required=True,
                        help='1-based indices of N and H atoms as pairs: e.g. 985 986 987 988')
    parser.add_argument('--output', default='PROJECTED_SPECTRUM.dat', help='Base output filename')
    parser.add_argument('--plot', action='store_true', help='Plot the spectra')
    args = parser.parse_args()

    if len(args.nh_indices) % 2 != 0:
        raise ValueError("Please provide nh_indices as pairs (even number of integers).")

    print("Reading static geometry...")
    coords = read_geometry(args.GEOMETRY)

    # build static bond unit vectors (one per pair)
    bond_unit_vectors = []
    n_bonds = len(args.nh_indices) // 2
    for i in range(n_bonds):
        idx1 = args.nh_indices[2*i] - 1
        idx2 = args.nh_indices[2*i + 1] - 1
        print(f"Bond {i+1}: atoms {idx1+1} - {idx2+1}")
        bv = get_bond_unit_vector(coords, idx1, idx2)
        bond_unit_vectors.append(bv)

    print("Reading dipole file and projecting on static bond vectors...")
    projections = get_proj_time_series_from_static_bonds(args.DIPOLE, bond_unit_vectors)

    spectra = []
    labels = []

    # compute spectrum for each bond separately
    for i, proj in enumerate(projections):
        print(f"Processing bond {i+1} projection -> autocorr -> spectrum ...")
        ac = get_autocorr_from_proj(proj)
        wn, spec = calc_spectrum_from_autocorr(ac, args.temperature, args.timestep)
        spectra.append((wn, spec))
        labels.append(f"Bond_{i+1}")

    # combined projection (sum of per-bond projections)
    combined_proj = np.sum(projections, axis=0)
    ac_comb = get_autocorr_from_proj(combined_proj)
    wn_comb, spec_comb = calc_spectrum_from_autocorr(ac_comb, args.temperature, args.timestep)
    spectra.append((wn_comb, spec_comb))
    labels.append("Combined")

    # save files
    base, ext = os.path.splitext(args.output)
    if ext == '':
        ext = '.dat'

    for (wn, spec), label in zip(spectra, labels):
        outname = f"{base}_{label}{ext}"
        print(f"Saving {label} to {outname}")
        save_spectrum(outname, wn, spec)

    print("All spectra saved.")

    if args.plot:
        plt.figure(figsize=(8,5))
        for (wn, spec), label in zip(spectra, labels):
            plt.plot(wn, spec, label=label)
        plt.xlabel('Wavenumber (cm$^{-1}$)')
        plt.ylabel('Normalized Intensity')
        plt.title('Projected Vibrational Spectra (static bond vectors)')
        plt.xlim(0, 4000)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("vibration_project_multi_staticbond.png", dpi=200)
        plt.show()

if __name__ == '__main__':
    main()
