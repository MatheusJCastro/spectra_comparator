########################################################
# plot_norms.py                                        #
# Matheus J. Castro                                    #
# Version 1.2                                          #
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


def finish_plot(show=False, save=False, fl1=None, fl2=None):
    # End and save plot subroutine

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

    onlynorm = False  # change to True to plot only the normalized spectrum

    files = []
    for i in os.listdir():  # search for all non-normalized files in the current directory
        if "tha_" in i and "norm" not in i and "list" not in i:
            files.append(i)

    files_norm = []
    for i in os.listdir():  # search for all normalized files in the current directory
        if "norm_tha_" in i:
            files_norm.append(i)

    for i in range(len(files)):  # for each tuple of files
        fig = plt.figure(figsize=(21, 9))
        fig.suptitle("Comparison of normalized and non normalized spectrum", fontsize=28)

        if not onlynorm:  # to plot non-normalized files as subplot
            plt.subplot(121)
            plt.title("Standard", fontsize=22)
            plt.xlabel("Pixel", fontsize=20)
            plt.ylabel("Intensity", fontsize=20)
            plt.yscale("log")

            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.tick_params(axis='y', which='minor', labelsize=16)

            spec_info = open_spec(files[i])  # open the current spectrum
            plot_spectra(spec_info)  # plot the spectrum

            plt.grid(True, which="both", linewidth=1)

            plt.subplot(122)

        plt.title("Normalized", fontsize=22)
        plt.xlabel("Pixel", fontsize=20)

        if onlynorm:
            plt.ylabel("Intensity", fontsize=20)

        plt.yscale("log")

        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.tick_params(axis='y', which='minor', labelsize=16)

        spec_info = open_spec(files_norm[i])  # open the current spectrum
        plot_spectra(spec_info)  # plot the spectrum

        plt.grid(True, which="both", linewidth=1)

        if files[i][-16:-12] == "3080":  # there are two spectra of the 3080A, save both without erasing one
            finish_plot(save=True, fl1="comp_norm", fl2=files[i][-16:-8])
        else:
            finish_plot(save=True, fl1="comp_norm", fl2=files[i][-16:-12])


if __name__ == '__main__':
    main()
