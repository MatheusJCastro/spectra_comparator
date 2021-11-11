########################################################
# divide_spectra.py                                    #
# Matheus J. Castro                                    #
# Version 2.0                                          #
# Last Modification: 11/11/2021 (month/day/year)       #
# https://github.com/MatheusJCastro/spectra_comparator #
# Licensed under MIT License                           #
########################################################

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import sys


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


def plot_spectra(spec1, spec2, specdiv, show=False, sharey=False):
    # Subroutine to plot the spectra

    specs = [spec1, spec2, specdiv]  # spectra data
    names = ["New Lamp", "Old Lamp", "Ratio"]  # spectra names
    colors = ["blue", "orange", "green"]  # spectra colors

    # noinspection PyTypeChecker
    fig, axes = plt.subplots(figsize=(21, 9), nrows=len(specs), sharex=True, sharey=sharey)
    axes[0].set_title("Difference between lamps", fontsize=28)

    for i, spec in enumerate(specs):
        axes[i].plot(spec[0], spec[1], label=names[i], color=colors[i])
        # axes[i].set_yscale("log")
        # axes[i].set_yticklabels([])
        axes[i].tick_params(axis='y', which='both', labelsize=16)
        axes[i].yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
        axes[i].legend(loc="upper left", fontsize=20)
        axes[i].grid(which="both")
        if i == 2 and not sharey:
            axes[i].set_ylim(0.99*min(spec[1]), 0.43*max(spec[1]))

    axes[len(specs)-1].set_xlabel("Wavelength (\u212b)", fontsize=20)
    plt.xticks(fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(hspace=.0)  # no horizontal space between the plots

    shy = ""
    if sharey:
        shy = "_sharey"

    plt.savefig("Lamps_Comparison{}".format(shy))

    if show:
        plt.show()

    plt.close()


def main(arg):
    # Main subroutine to plot the spectra

    # If one of the arguments passed is "sharey" the y-axis will be shared among the subplots
    sharey = False
    if any(i == "sharey" for i in arg):
        sharey = True

    # File names. 3 must be given
    files = ["multi_tha_novo.fits",
             "tha_velha.fits",
             "tha_divided.fits"]

    # Open the spectra
    old_spec = open_spec(files[0])
    new_spec = open_spec(files[1])
    div_spec = open_spec(files[2])

    for i in [True, False]:  # Bypass of the sharey variable
        plot_spectra(old_spec, new_spec, div_spec, sharey=i)


if __name__ == '__main__':
    args = sys.argv[1:]  # get the arguments passed in the execution
    main(args)
