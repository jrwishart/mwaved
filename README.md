mwaved
===========

mwaved is a set of functions that generalise the waved package for wavelet deconvolution in the Fourier domain. These generalisations are the following extensions:

* Allow a multichannel model where the practitioner has access to multiple channels of data of a common signal of interest.
* Allow additive long memory errors to be present in the multichannel signals (independent between signals but long memory within each signal)
* Allow a data-driven resolution choice in the presence of box car blur (waved does not have this property)

The code is also written with the use of the Rcpp package to help use the external C FFTW library to achieve speeds around 12-13 times faster than the usual WaveD package (comparing the performance of a single channel waved code to the same code in the mwaved package)

Installation instructions (currently only tested in linux)

* Ensure you have the FFTW3 libraries installed. For ubuntu this requires `sudo apt-get install libfftw3-dev`
* Download and install the package from your favourite CRAN repository. That is, run `install.packages('mwaved')` from the R prompt or download the tarball and run R CMD INSTALL mwaved_1.0.tar.gz from the linux terminal.

The package is being developed at http://github.com/justinrwishart/mwaved and any bug reports, comments or suggestions are welcomed at http://github.com/justinrwishart/issues
