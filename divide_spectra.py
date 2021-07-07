from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


def open_spec(fl_name):
    hdul = fits.open(fl_name)
    spec_data = hdul[0].data
    spec_header = hdul[0].header

    # CDELT1 or CD1_1
    wl = spec_header['CRVAL1'] + spec_header['CD1_1'] * np.arange(0, len(spec_data))
    hdul.close()

    return wl, spec_data, spec_header


def plot_spectra(spec1, spec2, specdiv, show=False, sharey=False):
    specs = [spec1, spec2, specdiv]
    names = ["Old Lamp", "New Lamp", "Ratio"]
    colors = ["blue", "orange", "green"]

    # noinspection PyTypeChecker
    fig, axes = plt.subplots(figsize=(21, 9), nrows=len(specs), sharex=True, sharey=sharey)
    axes[0].set_title("Difference between lamps", fontsize=24)

    for i, spec in enumerate(specs):
        axes[i].plot(spec[0], spec[1], label=names[i], color=colors[i])
        # axes[i].set_yscale("log")
        # axes[i].set_yticklabels([])
        axes[i].legend()
        axes[i].grid()

    axes[len(specs)-1].set_xlabel("Wavelength (\u212b)", fontsize=18)

    plt.tight_layout()
    plt.subplots_adjust(hspace=.0)

    shy = ""
    if sharey:
        shy = "_sharey"
    plt.savefig("Lamps_Comparison".format(shy))

    if show:
        plt.show()

    plt.close()


def main():
    files = ["tha_velha.fits",
             "tha_novo.fits",
             "tha_divided.fits"]

    old_spec = open_spec(files[0])
    new_spec = open_spec(files[1])
    div_spec = open_spec(files[2])

    plot_spectra(old_spec, new_spec, div_spec)


if __name__ == '__main__':
    main()
