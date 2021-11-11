# Spectra Comparator
The programs listed bellow aim to analyze and characterize the Photron Lamp spectrum and compare it with the Westinghouse Lamp, both are Thorium-Argon lamps. Lamps specifications can be found [here](http://iraf.noao.edu/specatlas/).  

The IRAF procedures used can be found [here](iraf.md)  

This project was done under grant #2020/14944-4, São Paulo Research Foundation (FAPESP) and the monograph with all the results are available [here](tese.pdf).

## plot_norms.py
Program to plot the spectra before and after the normalization. Search for two types of file names:

- `tha_` for non-normalized spectrum and;
- `norm_tha_` for normalized spectrum.

And the rest of the name needs to be the wavalength information.  
Obs: part of the name string will be used to write the plot name.  

Example:
![plot_norms example](Plots_comp_norm_3918.png)

## specplot.py
Program to visually verify the calibration. For **N** spectra, **N-1** plots will be generated, each with two spectra to compare the alignment of the emission lines.

Example:
![plot_norms example](Plots_3918_4000.png)

## calib.py
Comparison between a file named "thar.fits" (reference file) and any other *.fits* spectra with the string *tha_* in part of the name of the file.  

This program shows a interactive window created with Python3 that uses matplotlib to enable interact with the spectra while running the script. The plot created can also be saved.  

It has options to divide the spectra from 1 to 10 on horizontal axis. If the x axis is the same unit on both spetras it can also share the same limits on the axis to enable a better comparison.  

If the x axis is in Angstrons, the program has an option to show the wavelength of the peaks.  

Interface example:
![example](tiks_example.png)

## rms_compar.py

## divide_spectra.py

## analyse_spectra.py

## resolution_degradation.py

## spectrum_analysis.r

# Spectra Comparator
