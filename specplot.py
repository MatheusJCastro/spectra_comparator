###################################################
# Plot two-by-two Spectra for comparison          #
# Matheus J. Castro                               #
# Version 2.3                                     #
# Last Modification: 07/04/2021 (month/day/year)  #
###################################################

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os


def open_spec(fl_name):
    hdul = fits.open(fl_name)
    spec_data = hdul[0].data
    spec_header = hdul[0].header

    if spec_data.shape != (2048,):
        spec_data = spec_data[1][0]

    # CDELT1 or CD1_1
    wl = spec_header['CRVAL1'] + spec_header['CD1_1'] * np.arange(0, len(spec_data))
    hdul.close()

    return wl, spec_data, spec_header


def start_plot():
    plt.figure(figsize=(16, 9))
    plt.title("Spectra")
    plt.xlabel("Wavelength")
    plt.ylabel("Intensity")
    plt.yscale("log")


def finish_plot(show=False, save=False, fl1=None, fl2=None):
    plt.legend()
    plt.grid(True, which="both", linewidth=1)
    if save:
        plt.savefig("Plots_{}_{}".format(fl1, fl2))
    if show:
        plt.show()
    plt.close()


def plot_spectra(spec, name=None):
    plt.plot(spec[0], spec[1], label=name)


if __name__ == '__main__':
    files = []
    for i in os.listdir():
        if "calib_tha" in i:
            files.append(i)
    # files = ["tha_velha_4000_002.ms.fits", "tha_velha_3835_001.ms.fits", "tha_velha_3918_001.ms.fits"]

    for i in range(len(files[:-1])):
        start_plot()
        for j in [i, i+1]:
            spec_info = open_spec(files[j])
            plot_spectra(spec_info, name=files[j][-16:-12])
        finish_plot(save=True, fl1=files[i][-16:-12], fl2=files[i+1][-16:-12])
