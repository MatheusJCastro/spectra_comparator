########################################################
# spectrum_analysis.r                                  #
# Matheus J. Castro                                    #
# Version 2.3                                          #
# Last Modification: 11/11/2021 (month/day/year)       #
# https://github.com/MatheusJCastro/spectra_comparator #
# Licensed under MIT License                           #
########################################################

library("minpack.lm")

# Files Names
files = c("Spectrum_R=14026.csv", "Spectrum_R=6000.csv")

# Load Spectrum 1
spec_wave_1 = read.table(files[1], sep=",", skip=0, nrow=2, header=FALSE, row.names=1, check.names=FALSE)

spec_data_1 = read.table(files[1], sep="", skip=2, header=TRUE, check.names=FALSE)
spec_data_1 = spec_data_1[,1]

# Load Spectrum 2
spec_wave_2 = read.table(files[2], sep=",", skip=0, nrow=2, header=FALSE, row.names=1, check.names=FALSE)

spec_data_2 = read.table(files[2], sep="", skip=2, header=TRUE, check.names=FALSE)
spec_data_2 = spec_data_2[,1]

# Calculate wavelengths
spec_wave_1 = spec_wave_1["CRVAL_1",1] + spec_wave_1["CD1_1",1] * seq.int(0, length(spec_data_1)-1)
spec_wave_2 = spec_wave_2["CRVAL_1",1] + spec_wave_2["CD1_1",1] * seq.int(0, length(spec_data_2)-1)

# Plot Spectra
png(file="Espectra.png", width=1200, height=675)
par(mfrow=c(1, 2))

plot(spec_wave_1, spec_data_1, main="Spectra Comparison", 
     xlab="Wavelength", ylab="Intensity", type="l", lty=1,
     cex.axis=2, cex.lab=1.5, cex.main=2)
lines(spec_wave_2, spec_data_2, col="blue")
labels = c(unlist(strsplit(sub(".*_", "", files[1]), "[.]"))[1])
labels = c(labels, unlist(strsplit(sub(".*_", "", files[2]), "[.]"))[1])
legend("topright", legend=c(labels),
       col=c("black", "blue"), lty=1, cex=2)
grid()

interval = 1550:1700
plot(spec_wave_1[interval], spec_data_1[interval], main="Spectra Comparison - zoom in", 
     xlab="Wavelength", ylab="Intensity", type="l", lty=1,
     cex.axis=2, cex.lab=1.5, cex.main=2)
lines(spec_wave_2[interval], spec_data_2[interval], col="blue")
labels = c(unlist(strsplit(sub(".*_", "", files[1]), "[.]"))[1])
labels = c(labels, unlist(strsplit(sub(".*_", "", files[2]), "[.]"))[1])
legend("topright", legend=c(labels),
       col=c("black", "blue"), lty=1, cex=2)
grid()

# Find Peaks
find_peaks <- function(spec_data, threshold=NULL, half_window=10){
       peaks = c()

       if(is.null(threshold))
              threshold = mean(spec_data)

       lenspec = length(spec_data)
       for(i in 1:lenspec){
              imin = i - half_window
              imax = i + half_window
              if(imin < 1)
                     imin = 1
              if(imax > lenspec)
                     imax = lenspec

              if(spec_data[i] == max(spec_data[imin:imax]) & spec_data[i] > threshold)
                     peaks = c(peaks, i)
       }

       return(peaks)
}

peaks1 = find_peaks(spec_data_1)
# segments(x0=spec_wave_1[peaks1], y0=spec_data_1[peaks1]*1.1, 
       #   y1=spec_data_1[peaks1]*1.1+10000, col="red")

peaks2 = find_peaks(spec_data_2)
# segments(x0=spec_wave_2[peaks2], y0=spec_data_2[peaks2]*1.1, 
       #   y1=spec_data_2[peaks2]*1.1+10000, col="green")

# Find Min Max of each Peak
find_min_max_gauss <- function(spec_data, peak){
       miny_old = spec_data[peak]
       j = peak
       while(TRUE){
              miny = spec_data[j]
              if(miny <= miny_old){
                     miny_old = miny
                     j = j - 1
              } else{
                     ind_miny = j
                     break
              }
       }

       maxy_old = spec_data[peak]
       j = peak
       while(TRUE){
              maxy = spec_data[j]
              if(maxy <= maxy_old){
                     maxy_old = maxy
                     j = j + 1
              } else{
                     ind_maxy = j
                     break
              }
       }

       return(c(ind_miny, ind_maxy))
}

# Gaussian function
gauss <- function(x, a, m, s){
       G = a * exp(-1/2 * ((x-m)/s)**2)
}

# Gaussian Fit and Spectrum Analysis
gaussian_fit <- function(spec_wave, spec_data, peaks, plot=FALSE){
       stds = c()
       new_peaks = c()

       for(i in peaks){
              ind_min_max = find_min_max_gauss(spec_data, i)

              x = spec_wave[ind_min_max[1]:ind_min_max[2]]
              y = spec_data[ind_min_max[1]:ind_min_max[2]]

              if(plot){
                     png(file=paste("Fits/Fit_", format(spec_wave[i], digits=2, nsmall=2), ".png", sep=""), width=1200, height=675)
                     plot(x, y, main="Gauss Fit", xlab="Wavelength", ylab="Intensity", type="l", lty=1, lwd=5,
                          cex.axis=2, cex.lab=1.5, cex.main=2)
              }

              chi_square <- function(par){
                     y_gauss = gauss(x, par[1], par[2], par[3])
                     chi = y - y_gauss
                     if(plot) lines(x, gauss(x, par[1], par[2], par[3]), col="lightblue", lty=2)
                     return(chi)
              }

              fit = nls.lm(par=c(1, spec_wave[i], 0.1),# upper=c(280000, spec_wave[i]+100, 5), 
                                 fn=chi_square, control=nls.lm.control(maxiter=500))
              
              if(plot){
                     x = seq(min(x), max(x), 0.01)
                     lines(x, gauss(x, fit$par[1], fit$par[2], fit$par[3]), col="red", lty=2, lwd=2)
                     graphics.off()
              }
              
              if(fit$par[3] <= 1 & fit$par[3] > 0){
                     stds = c(stds, fit$par[3])
                     new_peaks = c(new_peaks, i)
              }
       }
       return(list(std=stds, peak=new_peaks))
}

# Standard Devations
fit_res = gaussian_fit(spec_wave_1, spec_data_1, peaks1)
stds1 = fit_res$std
std1p = format(mean(stds1), digits=2, nsmall=2)
peaks1 = fit_res$peak

fit_res = gaussian_fit(spec_wave_2, spec_data_2, peaks2)
stds2 = fit_res$std
std2p = format(mean(stds2), digits=2, nsmall=2)
peaks2 = fit_res$peak

# Plot Stds
png(file="Desvio_padrao.png", width=1200, height=675)
par(mfrow=c(1, 2))

hist(stds1, main=paste("Standard Deviation of Original Spectrum Peaks\nwith mean", std1p), 
     ylab="Frequency", xlab="Values", breaks=15, col="lightblue",
     cex.axis=2, cex.lab=1.5, cex.main=2)
abline(v=mean(stds1), lty=2, col="red")

hist(stds2, main=paste("Standard Deviation of Reduced Spectrum Peaks\nwith mean", std2p), 
     ylab="Frequency", xlab="Values", breaks=15, col="lightblue",
     cex.axis=2, cex.lab=1.5, cex.main=2)
abline(v=mean(stds2), lty=2, col="red")

# Plot Correlation of Stds and Wavelength
png(file="Desvio_Padrao_por_Comp_Onda.png", width=1200, height=675)
par(mfrow=c(1, 2))

plot(spec_wave_1[peaks1], stds1, main="Original Spectrum\nMean of Peaks x Standard Deviation", 
     xlab="Wavelength", ylab=expression(sigma), lty=1, col="blue",
     cex.axis=2, cex.lab=1.5, cex.main=2)

plot(spec_wave_2[peaks2], stds2, main="Reduced Spectrum\nMean of Peaks x Standard Deviation", 
     xlab="Wavelength", ylab=expression(sigma), lty=1, col="red",
     cex.axis=2, cex.lab=1.5, cex.main=2)

# Change of means for each peak
cross_match <- function(spec1, spec2, radius=10){
       matchs = matrix(, nrow=0, ncol=2)
       colnames(matchs) <- c("Spec1", "Spec2")

       if(length(spec1) >= length(spec2)){
              spec_big = spec1
              spec_few = spec2
       } else{
              spec_big = spec2
              spec_few = spec1
       }

       for(i in 1:length(spec_big)){
              val1 = spec_big[i]
              for(j in 1:length(spec_few)){
                     val2 = spec_few[j]
                     if(val2-radius <= val1 & val2+radius >= val1){
                            if(i %in% matchs[,1]){
                                   ind_rep = which(matchs[,1]==i)
                                   if(abs(val1-val2) < abs(val1-spec_few[matchs[ind_rep,2]]))
                                          matchs[ind_rep,2] = j
                            } else{
                                   matchs = rbind(matchs, c(i, j))
                            }
                     }
              }
       }

       if(length(spec1) < length(spec2)){
              order = c("Spec2", "Spec1")
              matchs = matchs[,order]
       }

       return(matchs)
}

match = cross_match(spec_wave_1[peaks1], spec_wave_2[peaks2])
matchc1 = spec_wave_1[peaks1[match[,1]]]
matchc2 = spec_wave_2[peaks2[match[,2]]]

# Plot Change of Means
png(file="Mudanca_das_Medias.png", width=1200, height=675)
par(mfrow=c(1, 2))

plot(matchc1, abs(matchc1-matchc2), 
     main="Wavelengths of Original Spectrum\nPeaks x Drifting from Reduced Spectrum", 
     xlab="Original Spectrum Wavelength", ylab="Drifting from Reduced Spectrum Peaks",
     cex.axis=2, cex.lab=1.5, cex.main=2)
abline(h=quantile(abs(matchc1-matchc2), prob=c(.75)), lty=2, col="red")

hist(abs(matchc1-matchc2), breaks=10, 
     main="Histogram of Peak Value Differences",
     xlab="Change of Means", ylab="Frequency", col="lightblue",
     cex.axis=2, cex.lab=1.5, cex.main=2)
abline(v=mean(abs(matchc1-matchc2)), lty=2, col="red")

# Quartile of Peaks Change
save_txt = "Qartile of Peaks Change"
save_txt = paste(save_txt, "\n50% =", quantile(abs(matchc1-matchc2), prob=c(.50)))
save_txt = paste(save_txt, "\n75% =", quantile(abs(matchc1-matchc2), prob=c(.75)))

# Identificacion of Blended Lines
blend_indentify <- function(spec1, spec2, radius=0.1){
       blends = matrix(, nrow=0, ncol=2)
       colnames(blends) <- c("Spec1", "Spec2")

       if(length(spec1) >= length(spec2)){
              spec_big = spec1
              spec_few = spec2
       } else{
              spec_big = spec2
              spec_few = spec1
       }

       for(i in 1:length(spec_big)){
              val1 = spec_big[i]
              actual_blend = c()
              for(j in 1:length(spec_few)){
                     val2 = spec_few[j]
                     if(val2-radius <= val1 & val2+radius >= val1){
                            actual_blend = c(actual_blend, j)
                     }

              }
              if(length(actual_blend) > 1){
                     for(j in 1:length(actual_blend))
                            blends = rbind(blends, c(i, actual_blend[j]))
              }
       }

       if(length(spec1) < length(spec2)){
              order = c("Spec2", "Spec1")
              blends = blends[,order]
       }

       return(blends)
}

blend = blend_indentify(spec_wave_1[peaks1], spec_wave_2[peaks2], radius=3*mean(abs(matchc1-matchc2)))
blends1_wave = spec_wave_1[peaks1[blend[,1]]]
blends2_wave = spec_wave_2[peaks2[blend[,2]]]
blends1_data = spec_data_1[peaks1[blend[,1]]]

# Find Correlation of Blend and Intensity
corr_blends <- function(bl1_w, bl1_d){
       intensities = c(50, 40, 25, 12.5, 6.25, 3.125, 0) * 1000

       intens_blends = list()

       for(i in 1:length(intensities)){
              blends_count = c()
              for(j in 1:length(bl1_w)){
                     if(bl1_d[j] >= intensities[i]){
                            if(i == 1)
                                   blends_count = c(blends_count, bl1_w[j])
                            else if(bl1_d[j] <= intensities[i-1])
                                   blends_count = c(blends_count, bl1_w[j])
                     }
              }
              intens_blends[[paste(intensities[i])]] = unique(blends_count)
       }
       
       return(intens_blends)
}

blends_correlation = corr_blends(blends1_wave, blends1_data)

# Plot Blended Lines in Spectra
png(file="Linhas_com_Blends.png", width=1200, height=675)
plot(spec_wave_1, spec_data_1, main="Blended Lines found for each Intensity level", 
     xlab="Wavelength", ylab="Intensity", type="l", lty=1,
     cex.axis=2, cex.lab=1.5, cex.main=2)

for(i in names(blends_correlation)){
       y_inten = rep(as.numeric(i), length(blends_correlation[[i]]))
       points(blends_correlation[[i]], y_inten, col=1+which(i==names(blends_correlation)), pch=4, lwd=3)
}

labels = names(blends_correlation)
legend("topright", legend=labels, title="Intensities levels",
       col=2:(length(names(blends_correlation))+1), lty=1, cex=2)

# Fit Correlation of Blend and Intensity
e_func <- function(x, a, b, c){
       E = b * exp(a*x) + c
}

chi_square_e_func <- function(x, y, par){
       e_teor = e_func(x, par[1], par[2], par[3])
       chi = sum(((y - e_teor)^2)/e_teor)
       return(chi)
}

intens = as.numeric(names(blends_correlation))
fit = optim(par=c(-0.0001, 50, 1), fn=chi_square_e_func, x=as.numeric(names(blends_correlation)), y=lengths(blends_correlation))
save_txt = paste(save_txt, "\nExponential of b*e^(ax)+c")
save_txt = paste(save_txt, "\na =", fit$par[1])
save_txt = paste(save_txt, "\nb =", fit$par[2])
save_txt = paste(save_txt, "\nc =", fit$par[3])

x_inten = seq(min(intens), max(intens), length=1000)
y_count = e_func(x_inten, fit$par[1], fit$par[2], fit$par[3])

# Plot Blended Lines as Intensity function
png(file="Relacao_Blend_Intensidade.png", width=1200, height=675)

plot(intens, lengths(blends_correlation), main="Blended Lines as funtction of Intensity", 
     xlab="Intensity", ylab="Blend Line Count", pch=8, cex=2, lwd=2,
     cex.axis=2, cex.lab=1.5, cex.main=2)

lines(x_inten, y_count, lty=2, col="red")

# Save Results
fl <- file("Resultados.txt")
writeLines(save_txt, fl)
close(fl)
