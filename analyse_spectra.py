########################################################
# analyse_spectra.py                                   #
# Matheus J. Castro                                    #
# Version 3.7                                          #
# Last Modification: 11/11/2021 (month/day/year)       #
# https://github.com/MatheusJCastro/spectra_comparator #
# Licensed under MIT License                           #
########################################################

from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as matplot
import pandas as pd
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


def cut_spec(spec, report=False):
    # Cut the Spectra to the needed length

    if report:
        cut_num = -2000  # leaves only the end of the spectrum to plot for example only
    else:
        cut_num = len(spec[0])//8-100  # remove the beginning of the spectra

    spec = [spec[0][cut_num:-1], spec[1][cut_num:-1], spec[2][cut_num:-1]]

    return spec


def find_peaks(spec, half_window=10):
    # Find the peaks in the Spectra

    mult_factor = 1.01  # 1.01 of the median of the intensity values
    # Set two thresholds, one for positive values and one for negative values
    threshold_pos = np.median(spec[1]) * mult_factor
    threshold_neg = np.median(spec[1]) * (2 - mult_factor)

    peaks = []
    for i in range(len(spec[1])):
        imin = i - half_window
        imax = i + half_window
        if imin < 0:
            imin = 0
        if imax >= len(spec[1]):
            imax = len(spec[1])

        if spec[1][i] > threshold_pos and spec[1][i] == max(spec[1][imin:imax]):
            peaks.append([spec[0][i], spec[1][i]])
        if spec[1][i] < threshold_neg and spec[1][i] == min(spec[1][imin:imax]):
            peaks.append([spec[0][i], spec[1][i]])

    return np.array(peaks), (threshold_neg, threshold_pos)


def hunt_method(list1, list2):
    # Find the best match between the peaks on Spectra and the Linelist provided
    # using the hunt method (searching in an ordered list)

    indexes = []
    for i in range(len(list2)):
        if list2[i] <= list1[0]:
            indexes.append(np.array([0, i]))
        elif list2[i] >= list1[-1]:
            indexes.append(np.array([len(list1)-1, i]))
        else:
            gap = [0, len(list1) - 1]
            while True:
                new_gap = gap[0]+(gap[1]-gap[0])//2
                if list2[i] < list1[new_gap]:
                    gap[1] = new_gap
                else:
                    gap[0] = new_gap
                if gap[0]+1 == gap[1]:
                    diffs = [list2[i]-list1[gap[0]], list1[gap[1]]-list2[i]]
                    if diffs[0] <= diffs[1]:
                        indexes.append(np.array([gap[0], i]))
                    else:
                        indexes.append(np.array([gap[1], i]))
                    break
    indexes = np.array(indexes)

    return indexes


def format_linelist(fl):
    # Format the NIST Linelist for the use in this program

    df = pd.read_csv(fl)

    elem = []
    wave = []
    wave_unc = []
    for i in range(len(df)):
        elem.append("{} {}".format(df.element[i], df.sp_num[i]))
        wave.append(np.float(df["obs_wl_air(A)"][i][2:-1]))
        try:
            wave_unc.append(np.float(df["unc_obs_wl"][i][2:-1]))
        except ValueError:
            wave_unc.append(None)
    df.element = elem
    df["obs_wl_air(A)"] = wave
    df["unc_obs_wl"] = wave_unc

    df = df[["element", "obs_wl_air(A)", "unc_obs_wl"]]
    df.rename(columns={"element": "Element", "obs_wl_air(A)": "Value", "unc_obs_wl": "Uncertainty"},
              inplace=True)

    return df


def remove_lines(linelist, fl):
    # Remove Lines wich iraf didn't found a great fit

    peaks_out = np.loadtxt(fl)
    matchs = hunt_method(list(linelist.Value), peaks_out)

    linelist.drop(labels=matchs.T[0], axis=0, inplace=True)
    linelist.reset_index(inplace=True)

    return linelist


def get_from_linelist(fl, peaks, line_remove=None, format_csv=True):
    # Get the peaks of the Spectra that have correspondents on the Linelist

    if format_csv:
        linelist = format_linelist(fl)
    else:
        linelist = pd.read_csv(fl, names=["Value", "Element"], delim_whitespace=True)

    linelist.Value = linelist.Value.astype(np.float)

    if line_remove is not None:
        linelist = remove_lines(linelist, line_remove)

    indexes = hunt_method(list(linelist.Value), peaks[0])

    peak_values = []
    correspondents = []
    for i in range(len(indexes)):
        peak_values.append("{} {}".format(linelist.Value[indexes[i][0]], linelist.Element[i]))
        correspondents.append([linelist.Value[indexes[i][0]], linelist.Element[i]])

    del linelist

    return peak_values, correspondents


def change_lang(fl_name):
    # Change the latex longtable environment to portuguese

    with open(fl_name, 'r', encoding="ISO-8859-1") as file:
        fl = file.read()

    fl = fl.replace("Continued on next page", "Continua na próxima página")

    with open(fl_name, "w", encoding="ISO-8859-1") as file:
        file.write(fl)


def create_database(peaks_values, peaks, median):
    # Create the database with the founded peaks

    mean_rms = 2*0.00878  # 2 sigmas do rms médio
    peaks_values = np.array(peaks_values).T
    peaks = np.array(peaks)

    positives = []
    for i in peaks[1]:
        if i > median:
            positives.append(True)
        else:
            positives.append(False)

    df = pd.DataFrame()
    df["Peak (A)"] = peaks[0]
    df["Linelist Correspondent (A)"] = peaks_values[0].astype(np.float)
    df["Difference (A)"] = np.abs(df["Peak (A)"] - df["Linelist Correspondent (A)"])
    df["Element"] = peaks_values[1]
    df["Positive"] = positives
    df["Greater RMS ({})".format(mean_rms)] = df["Difference (A)"] > mean_rms

    dfnew = df[(df["Greater RMS ({})".format(mean_rms)] == True) & (df["Positive"] == True)]
    print("{} lines founded;\n{} aren't on linelist;\nThis is {:.2f}% of total."
          .format(len(df), len(dfnew), 100*len(dfnew)/len(df)))

    # noinspection PyTypeChecker
    df.to_csv("Linelist_results.csv", sep=",", index=False)

    dftex = df.drop("Positive", axis=1)
    dftex["Greater RMS ({})".format(mean_rms)] = dftex["Greater RMS ({})".format(mean_rms)].astype(int)

    dftex.to_latex("Linelist_results.tex",
                   escape=False,
                   index=False,
                   caption=r"Resultado obtidos a partir das simulações feitas. Sendo $\lambda_E$ o centro da " +
                           r"linha de emissão, $\lambda_0$ o valor mais próximo na lista de linhas, o ``Erro'' " +
                           r"como definido na Equação \ref{eq:errlambs}, a ``Transição'' sendo a transição " +
                           r"eletrônica e ``Rejeita'' sendo valores de 1 (verdadeiro) ou 0 (falso) para a rejeição da " +
                           "correspondência de acordo com 2 vezes o RMS ({:.4f}).".format(mean_rms),
                   header=[r"$\lambda_E$ ($\ang$)", r"$\lambda_0$ ($\ang$)", r"Erro ($\ang$)", "Transição", "Rejeita"],
                   float_format="%.4f",
                   label="tab:resultscaracphotron",
                   position="H",
                   column_format="rrrcc",
                   encoding="ISO-8859-1",
                   longtable=True)

    change_lang("Linelist_results.tex")

    return df, [len(df), len(dfnew), 100*len(dfnew)/len(df)]


def plot_spectra(fl_name, spec, peaks, peak_values, threshs, df, res, show=False, report=False):
    # Plot the Spectra and the founded peaks

    # show = True
    if report:
        plt.figure(figsize=(21, 9))
    else:
        plt.figure(figsize=(42, 9))

    plt.title("{} Spectrum".format(fl_name), fontsize=32)
    plt.xlabel("Wavelength (\u212b)", fontsize=30)
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    plt.tick_params(axis='y', which='minor', labelsize=28)

    if not report:
        plt.ylim(spec[1].min()*0.7, spec[1].max()*1.1)
        # plt.ylim(spec[1].min() * 0.7, spec[1].max() * 0.5)  # usar esse no spec dividido

    plt.yscale("log")

    plt.plot(spec[0], spec[1], label="Spectrum")
    count = [0, 0]
    for i, peak in enumerate(peaks):

        df_peaks = df[(df["Peak (A)"] == peak[0])]
        color = "black"
        for j in range(len(df_peaks)):
            if df_peaks.iloc[j, -1]:
                color = "red"

        labels = None
        if color == "black" and count[0] == 0:
            labels = "Peaks with\nLinelist Match\n({} lines)".format(res[0]-res[1])
            count[0] += 1
        elif color == "red" and count[1] == 0:
            labels = "Peaks without\nLinelist Match\n({} lines)".format(res[1])
            count[1] += 1

        if peak[1] > spec[1].max()*0.5:
            plt.vlines(peak[0], np.mean(spec[1]), np.mean(spec[1]) + 0.05,
                       zorder=5, colors=color, label=labels)
            plt.annotate(peak_values[i], (peak[0], np.mean(spec[1]) + 0.1),
                         rotation=90, color=color)
        elif peak[1] > np.median(spec[1]):
            plt.vlines(peak[0], peak[1] + 0.05, peak[1] + 0.1,
                       colors=color, label=labels)
            plt.annotate(peak_values[i], (peak[0], peak[1] + 0.15),
                         rotation=90, color=color)
        else:
            plt.vlines(peak[0], peak[1] - 0.05, peak[1] - 0.1,
                       colors=color, label=labels)
            plt.annotate(peak_values[i], (peak[0], peak[1] - 0.215),
                         rotation=90, color=color)

    plt.hlines(threshs, min(spec[0]), max(spec[0]), linestyle="-.", colors="green", zorder=10, label="Peak Threshold")

    plt.legend(loc="best", bbox_to_anchor=(1, 1), fontsize=24)
    plt.tight_layout()
    plt.savefig("Spectra_{}.pdf".format(fl_name))

    if show:
        plt.show()

    plt.close()


def plot_spectra_half(fl_name, spec, peaks, peak_values, threshs, df, res, div=4):
    # Plot the Spectra and the founded peaks but with subplots for visualization

    # noinspection PyTypeChecker
    fig, axes = plt.subplots(figsize=(16, 15), nrows=div, sharey=True)
    axes[0].set_title("{} Spectrum".format(fl_name), fontsize=24)

    for i in range(div):
        cut = [i*len(spec[0])//div, (i+1)*len(spec[0])//div]
        spec_cut = [spec[0][cut[0]: cut[1]], spec[1][cut[0]: cut[1]]]

        axes[i].tick_params(axis='both', which='both', labelsize=14)
        axes[i].set_yscale("log")
        # noinspection PyUnresolvedReferences
        axes[i].get_yaxis().set_minor_formatter(matplot.ticker.ScalarFormatter())
        # noinspection PyUnresolvedReferences
        axes[i].get_yaxis().set_major_formatter(matplot.ticker.ScalarFormatter())

        axes[i].set_xlim(min(spec_cut[0]), max(spec_cut[0]))

        axes[i].plot(spec_cut[0], spec_cut[1], label="Spectrum", linewidth=0.5)
        axes[i].hlines(threshs, min(spec_cut[0]), max(spec_cut[0]), linestyle="-.",
                       colors="green", zorder=10, label="Peak Threshold", linewidth=0.5)

        # Cutting arrays
        # Create an array with true or false
        bool_arr = (spec[0][cut[0]] <= peaks.T[0]) & (peaks.T[0] <= spec[0][cut[1]-1])
        # Select only the true values of bool_arr to the peaks and peak_values array
        peaks_cut = peaks[bool_arr]
        peak_values_cut = np.array(peak_values)[bool_arr]

        count = [0, 0]
        for j, peak in enumerate(peaks_cut):

            df_peaks = df[(df["Peak (A)"] == peak[0])]
            color = "black"
            for k in range(len(df_peaks)):
                if df_peaks.iloc[k, -1]:
                    color = "red"

            labels = None
            if color == "black" and count[0] == 0:
                labels = "Peaks with\nLinelist Match\n({} lines)".format(res[0] - res[1])
                count[0] += 1
            elif color == "red" and count[1] == 0:
                labels = "Peaks without\nLinelist Match\n({} lines)".format(res[1])
                count[1] += 1

            if peak[1] > spec[1].max() * 0.5:
                axes[i].vlines(peak[0], np.mean(spec[1]), np.mean(spec[1]) + 0.05,
                               zorder=5, colors=color, label=labels)
                axes[i].annotate(peak_values_cut[j], (peak[0], np.mean(spec[1]) + 0.1),
                                 rotation=90, color=color)
            elif peak[1] > np.median(spec[1]):
                axes[i].vlines(peak[0], peak[1] + 0.05, peak[1] + 0.1,
                               colors=color, label=labels)
                axes[i].annotate(peak_values_cut[j], (peak[0], peak[1] + 0.15),
                                 rotation=90, color=color)
            else:
                axes[i].vlines(peak[0], peak[1] - 0.05, peak[1] - 0.1,
                               colors=color, label=labels)
                axes[i].annotate(peak_values_cut[j], (peak[0], peak[1] - 0.215),
                                 rotation=90, color=color)

        if i == 0:
            axes[0].legend(bbox_to_anchor=(1.2, 1.05), fontsize=14)

    plt.xlabel("Wavelength (\u212b)", fontsize=20)

    plt.tight_layout()
    plt.savefig("Spectra_{}_subplot.pdf".format(fl_name))

    plt.close()


def main():
    # Main subroutine

    report = False  # variable to save a reduced spectrum for better visualization, only for use in reports

    # File to analyze
    # file = "tha_divided.fits"
    file = "tha_novo.fits"
    # file = "tha_velha.fits"

    # Linelist to use
    linelist = "linelist.txt"
    # linelist = "lines1.csv"

    # Remove lines of the analyzes
    linelist_remove = "remove.txt"
    # linelist_remove = None

    spec = open_spec(file)
    spec = cut_spec(spec, report=report)
    peaks, threshs = find_peaks(spec)
    # Change the format_csv variable to True if you are using the NIST linelist
    peak_values_str, peak_values = get_from_linelist(linelist, peaks.T,
                                                     line_remove=linelist_remove, format_csv=False)

    df, res = create_database(peak_values, peaks.T, np.median(spec[1]))
    plot_spectra(file, spec, peaks, peak_values_str, threshs, df, res, report=report)
    plot_spectra_half(file, spec, peaks, peak_values_str, threshs, df, res)


if __name__ == '__main__':
    main()
