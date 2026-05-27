# For CECAM MiMiC Summer School 2022 (Lausanne), adapted from
# Original code by Efrem Braun (https://github.com/EfremBraun/calc-ir-spectra-from-lammps)

import numpy as np
from scipy import fftpack, signal
from scipy.ndimage import gaussian_filter


PLANCK = 1.05457180013E-34  # in kg m^2 s^-1
BOLTZMANN = 1.38064852E-23  # in m^2 kg s^-2 K^-1
SPEED_OF_LIGHT = 299792458  # in m s^-1


def get_autocorr(filename):
    data = np.loadtxt(filename, usecols=(0, 4, 5, 6))
    frames, dipoles = data[:, 0], data[:, 1:4]
    n_frames = len(frames)

    dipoles_shifted = np.zeros((2*n_frames - n_frames % 2, 3))
    dipoles_shifted[n_frames//2:n_frames//2+n_frames,:] = dipoles

    autocorr_full = signal.fftconvolve(dipoles_shifted, dipoles[::-1], mode='same', axes=0)[-n_frames:]
    autocorr_full = np.sum(autocorr_full, axis=1) / np.arange(n_frames, 0, -1)

    return autocorr_full


def calc_spectrum(temperature, timestep, wn_min=0, wn_max=4000, filename='DIPOLE'):

    autocorr = get_autocorr(filename)

    lineshape = fftpack.dct(autocorr, type=1)[1:]
    lineshape_frequencies = np.linspace(0, 0.5/timestep, len(autocorr))[1:]
    lineshape_frequencies_wn = lineshape_frequencies/100.0/SPEED_OF_LIGHT
    mask = (lineshape_frequencies_wn >= wn_min) & (lineshape_frequencies_wn <= wn_max)

    field_description = lineshape_frequencies * (1.0-np.exp(-PLANCK*lineshape_frequencies/(BOLTZMANN*temperature)))
    spectrum = lineshape * field_description
    smooth_spectrum = gaussian_filter(spectrum, 2.5)
    peak_hight = abs(max(smooth_spectrum) - min(smooth_spectrum))
    smooth_spectrum /= peak_hight

    return lineshape_frequencies_wn[mask], smooth_spectrum[mask]


if __name__ == '__main__':

    TEMPERATURE = 310.0  # in K
    TIMESTEP = 0.24E-15  # in s

    wavenumbers, signal = calc_spectrum(TEMPERATURE, TIMESTEP)

    with open('SPECTRUM', 'w') as out:
        for w, s in list(zip(wavenumbers, signal)):
            out.write(f'{w:5.1f} {s: 1.5e}\n')
