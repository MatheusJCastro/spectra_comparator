#!/usr/bin/env python3

########################################################
# calib.py                                             #
# Matheus J. Castro                                    #
# Version 4.1                                          #
# Last Modification: 11/11/2021 (month/day/year)       #
# https://github.com/MatheusJCastro/spectra_comparator #
# Licensed under MIT License                           #
########################################################

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from astropy.io import fits
import matplotlib.pyplot as plt
import tkinter as tk
import numpy as np
# import sys
import os

window = tk.Tk()
width = 2560 // 4 * 3
height = width * 9 / 21 * 2
# posRight = int(window.winfo_screenwidth() / 2 - width / 2)
posRight = 2560 // 4 * 3
posDown = int(window.winfo_screenheight() / 2 - height / 2)
window.geometry("{}x{}+{}+{}".format(int(width), int(height), posRight, posDown))
window.title("Spectra Comparison")


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


def hunt_method(search_list, value, eps_desired=1E-3):
    # Subroutine to search in an ordered list
    
    invert = False
    if search_list[0] < search_list[-1]:
        invert = True

    lims = [0, len(search_list) - 1]

    eps = 1
    x = int(lims[0] + (lims[1] - lims[0]) / 2)

    while eps > eps_desired:
        x_old = x

        if search_list[x] > value:
            lims = [lims[0], x] if invert else [x, lims[1]]
        else:
            lims = [x, lims[1]] if invert else [lims[0], x]

        x = int(lims[0] + (lims[1] - lims[0]) / 2)
        eps = np.abs((x - x_old) / x)

    return x


def find_peaks(spec, half_window=15):
    # Subroutine to get all the emission lines on the spectrum

    threshold = np.mean(spec[1]) + max(spec[1])*0.025

    peaks = []
    for i in range(len(spec[1])):
        imin = i - half_window
        imax = i + half_window
        if imin < 0:
            imin = 0
        if imax >= len(spec[1]):
            imax = len(spec[1])

        if spec[1][i] > threshold and spec[1][i] == max(spec[1][imin:imax]):
            peaks.append([spec[0][i], spec[1][i]])

    return peaks


def plot_specs(spec_refer, spec_info, name, pks_labels=True, cont=1, xaxis=False, save=False):
    # Subroutine to get and plot the spectra

    # Cut the NOAO spectrum in the necessary length
    mid_x = int(name[-16:-12])
    if mid_x - 100 < spec_refer[0][0]:
        indl = 0
    else:
        indl = hunt_method(spec_refer[0], mid_x - 100)
    if mid_x + 100 > spec_refer[0][-1]:
        indr = len(spec_refer[0]) - 1
    else:
        indr = hunt_method(spec_refer[0], mid_x + 100)

    spec_refer = [spec_refer[0][indl:indr],
                  spec_refer[1][indl:indr],
                  spec_refer[2]]

    fig = plt.figure(figsize=(19.2, 9), dpi=100)

    for ind, spec_long in enumerate([spec_refer, spec_info]):
        sp_len = len(spec_long[0]) // cont
        unit_x = spec_long[2]["WAT1_001"]

        labels = {}

        for i in unit_x.split():
            item = i.split("=")
            labels[item[0]] = item[1]

        for j in range(cont):
            plt.subplot(2, cont, cont * ind + (j + 1))

            if "units" in labels:
                spec = [spec_long[0][sp_len * j:sp_len * (j + 1)],
                        spec_long[1][sp_len * j:sp_len * (j + 1)]]

                if j == 0:
                    plt.ylabel("Referência NOAO", fontsize=18) if spec_long == spec_refer \
                        else plt.ylabel(name, fontsize=18)
                plt.xlabel("{} ({})".format(labels["label"], labels["units"]), fontsize=18)

                if xaxis:
                    sp_len_ref = len(spec_refer[0]) // cont
                    xlim = (spec_refer[0][sp_len_ref * j], spec_refer[0][sp_len_ref * (j + 1) - 1])
                    plt.xlim(xlim)
                else:
                    plt.xlim(spec[0][0], spec[0][-1])

                plt.xticks(fontsize=16)
                plt.ylim(0, max(spec[1]) * 1.05)
                plt.yticks([])

                if pks_labels:
                    peaks = find_peaks(spec)
                    for i in peaks:
                        plt.annotate("{:.3f}".format(i[0]), (i[0], i[1]), rotation=90)
            else:
                j_rev = (cont-1) - j
                spec = [spec_long[0][sp_len * j_rev:sp_len * (j_rev + 1)],
                        spec_long[1][sp_len * j_rev:sp_len * (j_rev + 1)]]

                if j == 0:
                    plt.ylabel(name, fontsize=18)
                plt.xlabel(labels["label"], fontsize=18)

                if xaxis:
                    sp_len_ref = len(spec_refer[0]) // cont
                    xlim = (spec_refer[0][sp_len_ref * j], spec_refer[0][sp_len_ref * (j + 1) - 1])
                    plt.xlim(xlim)
                else:
                    plt.xlim(max(spec[0]), min(spec[0]))

                plt.xticks(fontsize=16)
                plt.yticks([])

            plt.plot(spec[0], spec[1])
            plt.grid()

    plt.tight_layout()

    if save:
        # Save the plot
        nm_pks = ""
        if pks_labels:
            nm_pks = "_withPeaks"
        plt.savefig("Calibration_{}_d{}{}".format(name[-16:-12], cont, nm_pks))
        plt.close()
    else:
        # Show the plot on the tiks window
        figFrame = tk.Frame(master=window)

        canvas = FigureCanvasTkAgg(fig, master=figFrame)  # A tk.DrawingArea.
        toolbar = NavigationToolbar2Tk(canvas, figFrame)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)

        figFrame.grid(row=1, column=0, columnspan=5)


def execute(fl_name, refer_spec_info, pks_labels=True, cont=1, xaxis=False, save=False):
    # Subroutine to plot the spectrum in the tiks window

    spec_info = open_spec(fl_name)
    plot_specs(refer_spec_info, spec_info, fl_name, pks_labels=pks_labels, cont=cont,
               xaxis=xaxis, save=save)


def _quit():
    # Subroutine to quit Tiks enviroment

    window.quit()
    window.destroy()


def main():
    # Main subroutine to create the Tiks enviroment

    reference_name = "thar.fits"  # name of the calibrated spectrum obtained in NOAO website

    files = []
    for i in os.listdir():  # search for all spectra in the directory
        if "tha_" in i:
            files.append(i)

    # Create the Tiks variables
    fl_plot = tk.StringVar()
    fl_plot.set(files[0])
    divisions = tk.IntVar()
    divisions.set(1)
    pks_lab = tk.BooleanVar()
    pks_lab.set(True)
    xaxis = tk.BooleanVar()
    xaxis.set(False)

    refer_spec_info = open_spec(reference_name)  # open the NOAO spectrum

    butFrame = tk.Frame(window, width=width, height=25)  # create the tiks window

    # Create the graphical interface
    tk.Label(butFrame, text="Select file:") \
        .place(x=0, y=0, width=100, height=25)
    tk.OptionMenu(butFrame, fl_plot, *files) \
        .place(x=100, y=0, width=300, height=25)
    tk.Button(butFrame, text="Plot",
              command=lambda: execute(fl_plot.get(), refer_spec_info,
                                      pks_labels=pks_lab.get(),
                                      cont=divisions.get(),
                                      xaxis=xaxis.get())) \
        .place(x=400, y=0, width=200, height=25)
    tk.Checkbutton(butFrame, text="Label on Peaks", variable=pks_lab, onvalue=True, offvalue=False) \
        .place(x=600, y=0, width=200, height=25)
    tk.Label(butFrame, text="Number of divisions:") \
        .place(x=800, y=0, width=150, height=25)
    tk.OptionMenu(window, divisions, *[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) \
        .place(x=950, y=0, width=50, height=25)
    tk.Checkbutton(butFrame, text="Share x-axis", variable=xaxis, onvalue=True, offvalue=False) \
        .place(x=1000, y=0, width=200, height=25)
    tk.Button(butFrame, text="Save",
              command=lambda: execute(fl_plot.get(), refer_spec_info,
                                      pks_labels=pks_lab.get(),
                                      cont=divisions.get(),
                                      xaxis=xaxis.get(),
                                      save=True)) \
        .place(x=1200, y=0, width=200, height=25)
    tk.Button(butFrame, text="Quit", command=_quit, bg='red') \
        .place(x=1400, y=0, width=200, height=25)

    butFrame.grid(row=0, column=0)

    # Tiks loop
    window.mainloop()


if __name__ == '__main__':
    main()
