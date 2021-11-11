########################################################
# resolution_degradation.py                            #
# Matheus J. Castro                                    #
# Version 1.0                                          #
# Last Modification: 11/11/2021 (month/day/year)       #
# https://github.com/MatheusJCastro/spectra_comparator #
# Licensed under MIT License                           #
########################################################

import sys

from astropy.convolution import Gaussian1DKernel, convolve
from scipy.optimize import curve_fit, minimize
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np


def open_spec(fl_name):
    # Subroutine to open the .fits spectrum and read it

    hdul = fits.open(fl_name)  # open the file
    spec_data = hdul[0].data  # get the data
    spec_header = hdul[0].header  # get the header

    # Get the wavelength information from the header
    # CDELT1 or CD1_1
    wl = spec_header['CRVAL1'] + spec_header['CD1_1'] * np.arange(0, len(spec_data))

    hdul.close()  # close the file

    return wl, spec_data, spec_header


def plot(spec1, res1, spec2, res2, save=False, show=False):
    # Plot the original and the result spectrum

    plt.figure(figsize=(16, 9))
    plt.title("Spectra")
    plt.xlabel("Wavelength")
    plt.ylabel("Intensity")

    plt.plot(spec1[0], spec1[1], label="R={:.0f}".format(res1))
    plt.plot(spec2[0], spec2[1], label="R={:.0f}".format(res2))

    plt.legend()
    plt.grid(True, which="both", linewidth=1)

    if save:
        plt.savefig("Plots")
    if show:
        plt.show()
    plt.close()


def crop_spec(spec):
    # Subroutine to cut the spectrum in the needed length

    lims = [3230, 4230]
    indr, indl = 0, 0

    for i, value in enumerate(spec[0]):
        if value >= lims[0]:
            indl = i
            break
    for i, value in enumerate(spec[0]):
        if value > lims[1]:
            indr = i
            break

    return spec[0][indl:indr], spec[1][indl:indr], spec[2]


def find_peaks(spec, threshold=None, half_window=10):
    # Subroutine to find the peaks in a spectrum

    peaks = []

    if threshold is None:
        threshold = np.mean(spec[1]) * 5

    lenspec = len(spec[0])
    for i in range(lenspec):
        imin = i - half_window
        imax = i + half_window
        if imin < 0:
            imin = 0
        if imax >= lenspec:
            imax = lenspec - 1

        if spec[1][i] == spec[1][imin:imax].max() and spec[1][i] > threshold:
            peaks.append(i)

    return peaks


def gauss_function(x, a, mean, std):
    # Gaussian function
    return a * np.exp(-1/2 * ((x-mean)/std)**2)


def find_min_max_gauss(spec, peak):
    # Subroutine to find the beginning and the end of the gaussian in the spectrum peak
    # by finding the minimum value before and after the center of the peak

    miny_old = spec[1][peak]
    j = peak
    while True:
        miny = spec[1][j]
        if miny <= miny_old:
            miny_old = miny
            j -= 1
        else:
            ind_miny = j
            break

    maxy_old = spec[1][peak]
    j = peak
    while True:
        maxy = spec[1][j]
        if maxy <= maxy_old:
            maxy_old = maxy
            j += 1
        else:
            ind_maxy = j
            break

    return ind_miny, ind_maxy


def find_resolution(spec):
    # Subroutine to find the resolution of a spectrum

    peaks = find_peaks(spec)  # find the emission lines

    resolutions = []
    # plot_fits = []

    for i in peaks:
        # For each peak, adjust a gaussian function

        ind_miny, ind_maxy = find_min_max_gauss(spec, i)  # find min and max wavelength of the peak

        x = spec[0][ind_miny:ind_maxy]  # cut the region of the wavelength of the spectrum
        y = spec[1][ind_miny:ind_maxy]  # cut the region of the intensity of the spectrum

        try:
            popt, pcov = curve_fit(gauss_function, x, y, p0=[1, spec[0][i], 0.1])  # fit a gaussian to the peak
        except RuntimeError:
            pass
        else:  # if doesn't return an error
            chisq = np.sum((y - gauss_function(x, *popt))**2 / gauss_function(x, *popt))  # find chi square
            fwhm = 2 * np.sqrt(2 * np.log(2)) * popt[2]  # find full width at half maximum
            resolutions.append(np.array([np.array(np.abs(spec[0][i]/fwhm)), chisq]))  # append the resolution and chi

            # plot_fits.append([x[0], x[-1], popt])

    resolutions = np.array(resolutions)

    del_res = []
    for i, chisq in enumerate(resolutions.T[1]):  # find bad fits points
        if chisq > np.median(resolutions.T[1]):
            del_res.append(i)
    resolutions = np.delete(resolutions, del_res, axis=0)  # remove bad fits points

    # plt.figure(figsize=(16, 9))
    # plt.plot(spec[0], spec[1])
    # for i in plot_fits:
    #     xfit = np.linspace(i[0], i[1], 200)
    #     plt.plot(xfit, gauss_function(xfit, *i[2]))
    # plt.show()

    return np.mean(resolutions.T[0])


def change_resolution(spec, res_desired=6e3, eps_desired=1e-2, std_init=2):
    # Change the spectrum resolution

    def res_find(std):
        # Function to be minimized

        g = Gaussian1DKernel(stddev=std)  # function used to change the resolution
        new_spec = convolve(spec[1], g)  # convolution between the spectrum and the function above
        new_spec = [spec[0], new_spec, spec[2]]
        res = find_resolution(new_spec)

        # Minimize the error of the actual resolution and the desired one
        eps = np.abs((res - res_desired) / res)

        return eps

    # Use the Nelder-Mead method to minimize the function and find the best parameter to convolve the spectrum
    min_std = minimize(res_find, np.array([std_init]), method='Nelder-Mead').x

    # Use the result of the minimization to convolve the spectrum
    g_best = Gaussian1DKernel(stddev=min_std)
    new_spec_best = convolve(spec[1], g_best)

    return spec[0], new_spec_best, spec[2]


def save_spec_txt(spec, res):
    # Subroutine to save the spectrum in .csv format

    head = "CRVAL_1, {}\nCD1_1, {}\nData".format(spec[0][0], spec[2]['CD1_1'])
    np.savetxt("Spectrum_R={:.0f}.csv".format(res), spec[1], header=head, comments="")


def main():
    # Main subroutine to degradate a spectrum

    fl_name = "thar_photron.fits"  # file name to degradate
    spec = open_spec(fl_name)  # open spectrum

    spec = crop_spec(spec)  # crop the spectrum in the needed length

    new_spec = change_resolution(spec)  # change the resolution of the spectrum

    res1 = find_resolution(spec)  # find the resolution of the original spectrum
    res2 = find_resolution(new_spec)  # find the resolution of the new spectrum

    save_spec_txt(spec, res1)  # save the original spectrum in csv format
    save_spec_txt(new_spec, res2)  # save the new spectrum in csv format

    plot(spec, res1, new_spec, res2, show=True, save=True)  # plot both spectra


if __name__ == '__main__':
    main()
