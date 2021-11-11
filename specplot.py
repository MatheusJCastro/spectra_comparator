########################################################
# specplot.py                                          #
# Matheus J. Castro                                    #
# Version 1.1                                          #
# Last Modification: 11/11/2021 (month/day/year)       #
# https://github.com/MatheusJCastro/spectra_comparator #
# Licensed under MIT License                           #
########################################################

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os


def open_spec(fl_name):
    # Subroutine to open the .fits spectrum and read it

    hdul = fits.open(fl_name)  # open the file
    spec_data = hdul[0].data  # get the data
    spec_header = hdul[0].header  # get the header

    if spec_data.shape != (2048,):  # get only the actual spectrum (for multidimensional data)
        spec_data = spec_data[1][0]

    # Get the wavelength information from the header
    # CDELT1 or CD1_1
    wl = spec_header['CRVAL1'] + spec_header['CD1_1'] * np.arange(0, len(spec_data))

    hdul.close()  # close the file

    return wl, spec_data, spec_header


def start_plot():
    # Start the plot subroutine
    
    plt.figure(figsize=(16, 9))
    plt.title("Calibration Alignment", fontsize=28)
    plt.xlabel("Wavelength", fontsize=20)
    plt.ylabel("Intensity", fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tick_params(axis='y', which='minor', labelsize=16)
    plt.yscale("log")


def finish_plot(show=False, save=False, fl1=None, fl2=None):
    # End and save plot subroutine

    plt.tight_layout()
    plt.legend(fontsize=20)
    plt.grid(True, which="both", linewidth=1)
    if save:
        plt.savefig("Plots_{}_{}".format(fl1, fl2))
    if show:
        plt.show()
    plt.close()


def plot_spectra(spec, name=None):
    # Subroutine to plot the spectrum
    plt.plot(spec[0], spec[1], label=name)


def main():
    # Main subroutine, find and plot the spectra

    files = []
    for i in os.listdir():
        if "calib_tha" in i:  # search for the calibrated spectra
            files.append(i)
    # files = ["tha_velha_4000_002.ms.fits", "tha_velha_3835_001.ms.fits", "tha_velha_3918_001.ms.fits"]

    for i in range(len(files[:-1])):  # for each tuple of spectra
        start_plot()  # initial parameters of the plot
        for j in [i, i + 1]:
            spec_info = open_spec(files[j])  # open the spectrum
            plot_spectra(spec_info, name="File {}".format(files[j][-16:-12]))  # plot it
        finish_plot(save=True, fl1=files[i][-16:-12], fl2=files[i + 1][-16:-12])  # end the current plot


if __name__ == '__main__':
    main()
