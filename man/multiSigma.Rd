\name{multiSigma}
\alias{multiSigma}
\title{multiSigma}
\usage{
multiSigma(Y, deg = 3L)
}
\arguments{
  \item{Y}{An input signal either an n by M matrix
  containing the multichannel signal to be analysed (or
  single vector of n elements). Each of the M columns
  represents a channel of n observations.}

  \item{deg}{The degree of the auxiliary polynomial used in
  the Meyer wavelet.}
}
\value{
a vector of estimates containing the level of sigma in each
of the M channels.
}
\description{
Estimates the level of noise (standard deviation of noise)
in each channel.
}
\details{
This function estimates the noise present in each channel
by computing the Meyer wavelet transform of each channel of
the multichannel signal. In particular, the noise level is
computed by using the Median Absolute Deviation (MAD) of
the wavelet coefficients at the highest possible scale at J
= floor(log2(n) - 1).
}
\examples{
library(mwaved)
# Simulate matrix of Gaussian variables with three different noise levels
sig <- c(0.25, 0.5, 1.25)
Y <- sapply(1:3, function(i) sig[i]*rnorm(1024))
# Estimate the noise levels
multiSigma(Y, deg = 3)
}

