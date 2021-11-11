# Steps to Reduce Lamps Spectra
## Spectra Website

[Noao](http://iraf.noao.edu/specatlas/)

## Activate IRAF
Open a terminal and type:
```bash
    cd iraf/
    conda activate iraf27
    xgterm -sb &
    cl
```

## IRAF Command File Execution
```bash
    cl < nome_do_arquivo
```

## Raw Spectra Extraction
### Folders to Enter
```bash
    noao
    imred
    ccdred
    twodspec
    longslit
    specred
```

### Appal Command
```bash
    apall name.fits background=no trace=no
```

### List of Usefull Letters
* **d** - delete appertures
* **m** - create appertures
* **n** - create appertures no peak
* **c** - show center
* **l** - define min of an apperture
* **u** - define max of an apperture
* **b** - background edition of an apperture
* **z** - delete initial backgrounds
* **s** - select background
* **f** - fit background
* **?** - help

### Steps After Apall Command
1. Use **d** to delete all automatic appertures;
2. Create a new centered apperture using **n**;
3. Using **l** and **u** define the limits of the new apperture, something arround more or less 100 pixels;
4. Use **b** to edit the background;
5. With **z**, delete all backgrounds regions previously created;
6. Fit the background with **f** and return with **q**;
7. Press **q** again to finish and **return** key answering the questions to save.


## Spectra Calibration
### Normalization
```bash
    continuum @tha_list_in @tha_list_out order=10 high_reject=2
```


### Calibrating
```bash
    identify name.ms.fits coordli=linelist.txt funct=chebyshev order=2
```
Obs: order 2 to start calibration, this can be changed durring calibration process.

### List of Usefull Letters
* **m** - slecect a peak to specify the wavelength
* **f** - to fit
* **d** - delete
* **l** - calibrate
* **q** - return or quit

### Try multiples fits:
Go to fit (**f**) page, delete with **d** the largests errors, **f** again to fit. If you see any patterns, try to change the fit order typing on the fit page:
```bash
:o ORDER_NUMBER
```
Increase the order slowly until the pattern is gone. Allways fit again the lines. Try to go back with **q** to the spetrum page and type **l** to findo others lines, go to the fit page. Do this process until your error is smaller than 0.025A.

### Steps After Identify Command
#### Mirror Spectra
1. Using **m** specify the wavelength of the borders of the spectra choosing the lower wavelength to the higher pixel count and vice versa;
2. Press **f** to fit and them **q** to return.

#### Mirror Spectra 2
Press **w** then **f** to mirror.

#### Calibrating
1. Using **m** specify all peaks that the wavelegth is known;
2. Use **f** to fit and return with **q**;
3. Calibrate usign **l** (as the first letter of laboratory);
4. Press **q** to exit and **return** key to answer the questions to save.

### Syncing the file with the Calibration
```bash
    dispcor name.ms.fits name.ms.fits
```

## Spectra Combine
```bash
    scombine @name_list_in name_out.fits
```

## Spectrum multiplication
```bash
    imarith name.fits * 1.7 new_name.fits
    imarith new_name.fits - 0.7 new_name.fits
```

## Spectra Division
```bash
    sarith name1.fits / name2.fits name_out.fits w1=3163 w2=4084
```

## Plot Spetrum
```bash
splot name.fits
```
### For multiple Spectra:
```bash
specplot spec1,spec2
```

## Zoom Iraf
- **w**: start  
- **e**: select range (x,y)  
- **a**: go back  







