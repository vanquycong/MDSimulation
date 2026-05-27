### Usage: python vibration_project.py GEOMETRY.xyz DIPOLE --temperature 310 --timestep 0.24e-15 --nh_indices 984 985 --plot
### Calculate the projection on a static choosen vector, e.g. N-H as shown here. 

import numpy as np
from scipy import signal, fftpack
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

PLANCK = 1.05457180013E-34  # in kg m^2 s^-1
BOLTZMANN = 1.38064852E-23  # in m^2 kg s^-2 K^-1
SPEED_OF_LIGHT = 299792458  # in m s^-1

def read_geometry(geom_file):
    with open(geom_file) as f:
        lines = f.readlines()
    coords = []
    for line in lines[2:]:  # skip header lines
        parts = line.split()
        if len(parts) >= 4:
            x, y, z = map(float, parts[1:4])
            coords.append([x, y, z])
    return np.array(coords)

def get_bond_unit_vector(coords, idx1, idx2):
    vec = coords[idx2] - coords[idx1]
    return vec / np.linalg.norm(vec)

def get_projected_dipole_autocorr(dipole_file, bond_unit_vector):
    data = np.loadtxt(dipole_file, usecols=(0, 4, 5, 6))
    dipoles = data[:, 1:4]  # shape (n_frames, 3)
    dipole_proj = np.dot(dipoles, bond_unit_vector)  # scalar projection time series

    n_frames = len(dipole_proj)
    dipole_proj_shifted = np.zeros(2 * n_frames - n_frames % 2)
    dipole_proj_shifted[n_frames // 2:n_frames // 2 + n_frames] = dipole_proj

    autocorr_full = signal.fftconvolve(dipole_proj_shifted, dipole_proj[::-1], mode='same')
    autocorr = autocorr_full[-n_frames:] / np.arange(n_frames, 0, -1)
    return autocorr

def calc_spectrum_from_autocorr(autocorr, temperature, timestep, wn_min=0, wn_max=4000):
    lineshape = fftpack.dct(autocorr, type=1)[1:]
    n_points = len(lineshape)
    lineshape_freq = np.linspace(0, 0.5 / timestep, n_points + 1)[1:]  # frequencies in Hz
    lineshape_freq_wn = lineshape_freq / 100.0 / SPEED_OF_LIGHT  # convert to cm^-1

    mask = (lineshape_freq_wn >= wn_min) & (lineshape_freq_wn <= wn_max)

    field_desc = lineshape_freq * (1.0 - np.exp(-PLANCK * lineshape_freq / (BOLTZMANN * temperature)))
    spectrum = lineshape * field_desc
    smooth_spectrum = gaussian_filter(spectrum, 2.5)

    peak_height = np.abs(np.max(smooth_spectrum) - np.min(smooth_spectrum))
    smooth_spectrum /= peak_height

    return lineshape_freq_wn[mask], smooth_spectrum[mask]

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Calculate vibrational spectrum projected on N–H bond from CPMD DIPOLE")
    parser.add_argument('GEOMETRY', help='Path to GEOMETRY.xyz file')
    parser.add_argument('DIPOLE', help='Path to DIPOLE file')
    parser.add_argument('--temperature', type=float, default=310.0, help='Temperature in Kelvin (default 310 K)')
    parser.add_argument('--timestep', type=float, default=0.24e-15, help='Timestep in seconds (default 0.24 fs)')
    parser.add_argument('--nh_indices', type=int, nargs=2, default=[985, 986],
                        help='1-based indices of N and H atoms of the bond (default: 985 986)')
    parser.add_argument('--output', default='PROJECTED_SPECTRUM.dat', help='Output spectrum filename')
    parser.add_argument('--plot', action='store_true', help='Plot the vibrational spectrum')
    args = parser.parse_args()

    print("Reading geometry file...")
    coords = read_geometry(args.GEOMETRY)
    idx_N = args.nh_indices[0] - 1  # convert to zero-based
    idx_H = args.nh_indices[1] - 1

    print(f"Computing N–H bond vector between atoms {idx_N +1} and {idx_H +1}...")
    bond_vec = get_bond_unit_vector(coords, idx_N, idx_H)

    print("Reading dipole file and calculating projected dipole autocorrelation...")
    autocorr = get_projected_dipole_autocorr(args.DIPOLE, bond_vec)

    print("Calculating vibrational spectrum...")
    wn, spectrum = calc_spectrum_from_autocorr(autocorr, args.temperature, args.timestep)

    print(f"Saving spectrum to {args.output} ...")
    with open(args.output, 'w') as f:
        for w, s in zip(wn, spectrum):
            f.write(f"{w:8.2f} {s:12.6e}\n")

    print("Spectrum saved.")

    if args.plot:
        import matplotlib.pyplot as plt
        plt.plot(wn, spectrum, label='Projected N–H Bond Spectrum')
        plt.xlabel('Wavenumber (cm$^{-1}$)')
        plt.ylabel('Normalized Intensity')
        plt.title('Vibrational Spectrum Projected on N–H Bond')
        plt.xlim(0, 4000)
        plt.legend()
        plt.grid(True)
        plt.savefig("vibration_project.png")
        plt.show()

if __name__ == '__main__':
    main()
