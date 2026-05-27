#python vibration_project.py GEOMETRY.xyz DIPOLE --temperature 310 --timestep 0.24e-15 --nh_indices 985 986 986 987 --plot

#!/usr/bin/env python3
import numpy as np
from scipy import signal, fftpack
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

PLANCK = 1.05457180013E-34  # kg m^2 s^-1
BOLTZMANN = 1.38064852E-23  # m^2 kg s^-2 K^-1
SPEED_OF_LIGHT = 299792458  # m/s

def read_geometry_cpmd(filename):
    frames = []
    with open(filename) as f:
        lines = f.readlines()
    n_atoms = int(lines[0].strip())  # first line = number of atoms
    lines_per_frame = n_atoms + 2    # atom lines + 2 header lines
    for i in range(0, len(lines), lines_per_frame):
        frame_lines = lines[i+2 : i+2+n_atoms]
        frame = []
        for line in frame_lines:
            parts = line.split()
            x, y, z = map(float, parts[1:4])
            frame.append([x, y, z])
        frames.append(frame)
    return np.array(frames)

def read_dipoles(filename):
    """Reads dipole moment vectors from CPMD DIPOLE file."""
    # Adjust usecols if your DIPOLE format differs!
    data = np.loadtxt(filename, usecols=(0, 4, 5, 6))
    dipoles = data[:, 1:4]  # shape: (n_frames, 3)
    return dipoles

def bond_projection_timeseries(coords, idx_N, idx_H, dipoles):
    """Project dipole moment onto NH bond vector (updated per frame)."""
    n_frames = coords.shape[0]
    proj = np.zeros(n_frames)
    for i in range(n_frames):
        bond_vec = coords[i, idx_H] - coords[i, idx_N]
        bond_vec /= np.linalg.norm(bond_vec)
        proj[i] = np.dot(dipoles[i], bond_vec)
    return proj

def autocorr(x):
    """Compute autocorrelation using FFT convolution."""
    result = signal.fftconvolve(x, x[::-1], mode='full')
    mid = result.size // 2
    return result[mid:] / np.arange(len(x), 0, -1)

def spectrum_from_acf(acf, temperature, timestep, wn_min=0, wn_max=4000, smooth_sigma=2.5):
    """Compute vibrational spectrum from autocorrelation."""
    lineshape = fftpack.dct(acf, type=1)[1:]
    n_points = len(lineshape)
    freq_Hz = np.linspace(0, 0.5 / timestep, n_points + 1)[1:]
    freq_cm = freq_Hz / (100.0 * SPEED_OF_LIGHT)

    mask = (freq_cm >= wn_min) & (freq_cm <= wn_max)

    # Quantum correction factor
    qcf = freq_Hz * (1.0 - np.exp(-PLANCK * freq_Hz / (BOLTZMANN * temperature)))
    spectrum = lineshape * qcf

    smooth_spectrum = gaussian_filter(spectrum, smooth_sigma)
    # Normalize
    smooth_spectrum /= np.max(np.abs(smooth_spectrum))

    return freq_cm[mask], smooth_spectrum[mask]

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Projected vibrational spectra for one or more N–H bonds from CPMD DIPOLE + XYZ.")
    parser.add_argument('GEOMETRY', help='Path to GEOMETRY.xyz (all frames)')
    parser.add_argument('DIPOLE', help='Path to DIPOLE file')
    parser.add_argument('--temperature', type=float, default=310.0, help='Temperature in Kelvin (default 310 K)')
    parser.add_argument('--timestep', type=float, default=0.24e-15, help='Timestep in seconds (default 0.24 fs)')
    parser.add_argument('--nh_indices', type=int, nargs='+', required=True,
                        help='List of N,H pairs (1-based), e.g., --nh_indices 985 986 990 991 for two bonds')
    parser.add_argument('--output_prefix', default='PROJECTED', help='Prefix for output files')
    parser.add_argument('--plot', action='store_true', help='Plot the spectra')
    args = parser.parse_args()

    # Validate
    if len(args.nh_indices) % 2 != 0:
        raise ValueError("You must give an even number of indices (N,H pairs).")

    pairs = [(args.nh_indices[i] - 1, args.nh_indices[i+1] - 1)
             for i in range(0, len(args.nh_indices), 2)]
    n_bonds = len(pairs)

    print(f"Reading geometry with all frames: {args.GEOMETRY}")
    coords = read_geometry_cpmd(args.GEOMETRY)

    print(f"Reading dipole trajectory: {args.DIPOLE}")
    dipoles = read_dipoles(args.DIPOLE)

    if coords.shape[0] != dipoles.shape[0]:
        raise ValueError(f"Number of frames mismatch: coords={coords.shape[0]}, dipoles={dipoles.shape[0]}")

    print(f"Processing {n_bonds} N–H bonds...")
    projections = []
    for (idx_N, idx_H) in pairs:
        proj = bond_projection_timeseries(coords, idx_N, idx_H, dipoles)
        projections.append(proj)

    # Per-bond spectra
    spectra_data = []
    for i, proj in enumerate(projections):
        acf = autocorr(proj)
        wn, spec = spectrum_from_acf(acf, args.temperature, args.timestep)
        spectra_data.append(spec)
        np.savetxt(f"{args.output_prefix}_bond{i+1}.dat", np.column_stack((wn, spec)),
                   header="Wavenumber(cm^-1) Intensity(norm.)")
        print(f"Saved {args.output_prefix}_bond{i+1}.dat")

    # Combined without cross terms (just sum intensities)
    combined_no_cross = np.sum(spectra_data, axis=0)
    combined_no_cross /= np.max(combined_no_cross)
    np.savetxt(f"{args.output_prefix}_combined_nocross.dat", np.column_stack((wn, combined_no_cross)),
               header="Wavenumber(cm^-1) Intensity(norm.)")
    print(f"Saved combined spectrum without cross terms: {args.output_prefix}_combined_nocross.dat")

    # Combined with cross terms
    # First make combined projection time series
    combined_proj = np.sum(projections, axis=0)
    acf_combined = autocorr(combined_proj)
    wn, combined_with_cross = spectrum_from_acf(acf_combined, args.temperature, args.timestep)
    np.savetxt(f"{args.output_prefix}_combined_cross.dat", np.column_stack((wn, combined_with_cross)),
               header="Wavenumber(cm^-1) Intensity(norm.)")
    print(f"Saved combined spectrum WITH cross terms: {args.output_prefix}_combined_cross.dat")

    if args.plot:
        plt.figure(figsize=(8,6))
        for i, spec in enumerate(spectra_data):
            plt.plot(wn, spec, label=f'Bond {i+1}')
        plt.plot(wn, combined_no_cross, '--', label='Sum (no cross)')
        plt.plot(wn, combined_with_cross, '-', label='Sum (with cross)')
        plt.xlabel('Wavenumber (cm$^{-1}$)')
        plt.ylabel('Normalized Intensity')
        plt.legend()
        plt.title('Projected N–H Vibrational Spectra')
        plt.xlim(0, 4000)
        plt.grid(True)
        plt.savefig(f"{args.output_prefix}_spectra.png", dpi=300)
        plt.show()

if __name__ == '__main__':
    main()
