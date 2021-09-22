import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def main():
    data = pd.read_csv("rms.csv", delimiter=",\t", engine='python')
    old = data[["RMS_Old"]].to_numpy()
    old = old[2:]
    new = data[["RMS_New"]].to_numpy()

    data.to_latex("rms.tex", index=False, escape=False,
                  header=[r"Wavelength ($\ang$)", "Westinghouse Lamp", "Photron Lamp"],
                  caption="Valores RMS para cada calibração de acordo com cada lâmpada.")

    def set_ticks():
        plt.xlabel("Values", fontsize=20)
        plt.ylabel("Counts", fontsize=20)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

    def set_grid():
        plt.legend(fontsize=16)
        plt.grid()

    fig = plt.figure(figsize=(16, 9))

    fig.suptitle("Calibration RMS Values", fontsize=28, y=1)

    plt.subplot(1, 2, 1)
    plt.title("Old Lamp", fontsize=22)
    set_ticks()
    plt.hist(old, bins=30, color="blue")
    plt.axvline(np.mean(old), ymin=0, ymax=1, linestyle="dashed",
                color="red", label="   Mean: {:.5f}".format(np.mean(old)))
    plt.axvline(np.median(old), ymin=0, ymax=1, linestyle="dashed",
                color="magenta", label="Median: {:.5f}".format(np.median(old)))
    set_grid()

    plt.subplot(1, 2, 2)
    plt.title("New Lamp", fontsize=22)
    set_ticks()
    plt.hist(new, bins=30, color="green")
    plt.axvline(np.mean(new), ymin=0, ymax=1, linestyle="dashed",
                color="red", label="   Mean: {:.5f}".format(np.mean(new)))
    plt.axvline(np.median(new), ymin=0, ymax=1, linestyle="dashed",
                color="magenta", label="Median: {:.5f}".format(np.median(new)))
    set_grid()

    plt.tight_layout(pad=2, w_pad=0)
    plt.savefig("RMS_Compar")
    # plt.show()
    plt.close()


if __name__ == '__main__':
    main()
