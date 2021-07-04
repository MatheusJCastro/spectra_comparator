###################################################
# Plot Normalized Spectra and Non-Normalized      #
# Matheus J. Castro                               #
# Version 1.4                                     #
# Last Modification: 07/04/2021 (month/day/year)  #
###################################################

import sys

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


def finish_plot(show=False, save=False, fl1=None, fl2=None):
    if save:
        plt.savefig("Plots_{}_{}".format(fl1, fl2))
    if show:
        plt.show()
    plt.close()


def plot_spectra(spec, name=None):
    plt.plot(spec[0], spec[1], label=name)


if __name__ == '__main__':
    onlynorm = False
    files = []
    for i in os.listdir():
        if "tha_" in i and "norm" not in i and "list" not in i:
            files.append(i)
    files_norm = []
    for i in os.listdir():
        if "norm_tha_" in i:
            files_norm.append(i)

    for i in range(len(files)):
        plt.figure(figsize=(40, 18))
        if not onlynorm:
            plt.subplot(121)
            plt.title("Spectra")
            plt.xlabel("Pixel")
            plt.ylabel("Intensity")
            plt.yscale("log")

            spec_info = open_spec(files[i])
            plot_spectra(spec_info, name="standard")

            plt.legend()
            plt.grid(True, which="both", linewidth=1)

            plt.subplot(122)
        plt.title("Spectra")
        plt.xlabel("Pixel")
        plt.ylabel("Intensity")
        plt.yscale("log")

        spec_info = open_spec(files_norm[i])
        plot_spectra(spec_info, name="normalized")

        plt.legend()
        plt.grid(True, which="both", linewidth=1)

        finish_plot(save=True, fl1="comp_norm", fl2=files[i][-16:-12])
