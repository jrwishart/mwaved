mwaved
===========

mwaved is a set of functions that generalise the waved package for wavelet deconvolution in the Fourier domain. These generalisations are the following extensions:

* Allow a multichannel model where the practitioner has access to multiple channels of data of a common signal of interest.
* Allow additive long memory errors to be present in the multichannel signals (independent between signals but exhibit long memory within each signal)
* Allow a data-driven resolution choice in the presence of box car blur (waved does not have this feature)

The user is encouraged to view the embedded Shiny applet that showcases the mwaved pacakge and importantly lists the appropriate R commands to recreate the output given by Shiny applet. The embedded Shiny applet can be viewed as long as the user has the [shiny](http://cran.r-project.org/package=shiny) package installed on their machine and then using R command `> mwaved::mWaveDDemo()`. 

The code is also written with the use of the Rcpp package to help use the external C FFTW library to achieve speeds around 8-15 times faster than the usual WaveD package (comparing the performance of a single channel waved code to the same code in the mwaved package with various sample sizes). The relative performance improves as the sample size increases. 

The package is being developed at http://github.com/jrwishart/mwaved and any bug reports, comments or suggestions are welcomed at http://github.com/jrwishart/issues

Optional source compilation instructions (currently only tested in Ubuntu, Slackware Linux and Windows 8.1)

* Ensure you have the FFTW3 libraries installed. For ubuntu this requires `sudo apt-get install libfftw3-dev`. For Windows 8 this requires downloading the windows fftw3 binaries and adding the installed directories to your PATH.
* Download and install the package from your favourite CRAN repository. That is, run `install.packages('mwaved')` from the R prompt or download the tarball and run R CMD INSTALL mwaved_1.0.tar.gz from the linux terminal.
