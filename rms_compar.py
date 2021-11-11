########################################################
# rms_compar.py                                        #
# Matheus J. Castro                                    #
# Version 1.5                                          #
# Last Modification: 11/11/2021 (month/day/year)       #
# https://github.com/MatheusJCastro/spectra_comparator #
# Licensed under MIT License                           #
########################################################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def main():
    # Main subroutine

    data = pd.read_csv("rms.csv", delimiter=",\t", engine='python')  # open the .csv file with the rms for each file
    old = data[["RMS_Old"]].to_numpy()  # create the numpy arrays of the RMS of the Westinghouse lamp
    old = old[2:]
    new = data[["RMS_New"]].to_numpy()  # create the numpy arrays of the RMS of the Photron lamp

    # Save a table in latex form
    data.to_latex("rms.tex", index=False, escape=False,
                  header=[r"Wavelength ($\ang$)", "Westinghouse Lamp", "Photron Lamp"],
                  caption="Valores RMS para cada calibração de acordo com cada lâmpada.")

    def set_ticks():
        # Subroutine to set the axis ticks parameters
        plt.xlabel("Values", fontsize=20)
        plt.ylabel("Counts", fontsize=20)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

    def set_grid():
        # Subroutine to set the grid properties
        plt.legend(fontsize=16)
        plt.grid()

    fig = plt.figure(figsize=(16, 9))

    fig.suptitle("RMS Values obtained with spectra calibrations", fontsize=28, y=1)

    plt.subplot(1, 2, 1)  # first plot
    plt.title("Old Lamp", fontsize=22)
    set_ticks()
    plt.hist(old, bins=30, color="blue")
    plt.axvline(np.mean(old), ymin=0, ymax=1, linestyle="dashed",
                color="red", label="   Mean: {:.5f}".format(np.mean(old)))
    plt.axvline(np.median(old), ymin=0, ymax=1, linestyle="dashed",
                color="magenta", label="Median: {:.5f}".format(np.median(old)))
    set_grid()

    plt.subplot(1, 2, 2)  # second plot
    plt.title("New Lamp", fontsize=22)
    set_ticks()
    plt.hist(new, bins=30, color="green")
    plt.axvline(np.mean(new), ymin=0, ymax=1, linestyle="dashed",
                color="red", label="   Mean: {:.5f}".format(np.mean(new)))
    plt.axvline(np.median(new), ymin=0, ymax=1, linestyle="dashed",
                color="magenta", label="Median: {:.5f}".format(np.median(new)))
    set_grid()

    plt.tight_layout(pad=2.8, w_pad=0)
    plt.savefig("RMS_Compar")
    # plt.show()
    plt.close()


if __name__ == '__main__':
    main()
