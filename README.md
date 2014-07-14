mwaved
===========

mwaved is a set of functions that generalise the waved package for wavelet deconvolution in the Fourier domain. These generalisations are the following extensions:

* Allow a multichannel model where the practitioner has access to multiple channels of data of a common signal of interest.
* Allow additive long memory errors to be present in the multichannel signals (independent between signals but long memory within each signal)
* Allow a data-driven resolution choice in the presence of box car blur (waved does not have this property)

The code is also written with the use of the Rcpp package and uses the FFTW library to achieve speeds around 12-13 times faster than the usual WaveD package (comparing the performance of a single channel waved code to the same code in the mwaved package)