\name{mWaveD}
\alias{mWaveD}
\title{mWaveD}
\usage{
mWaveD(Y, G = directG(dim(as.matrix(Y))), alpha = rep(1,
  dim(as.matrix(Y))[1]), j0 = 3L, j1 = NA_integer_, blur = "direct",
  thresh = as.numeric(c()), eta = NA_real_, shrinkage = "Hard",
  deg = 3L)
}
\arguments{
  \item{Y}{The input multichannel signal in the form of an
  n by m matrix to denote the m separate channels with n
  observations in each channel. If a vector of size n is
  given, it is parsed as a single channel signal with n
  elements.}

  \item{G}{The input multichannel blur matrix of size n by
  m (needs to be the same size as the signal input). This
  matrix dictates the form of blur present in each of the
  channels.}

  \item{alpha}{A vector of length m that specifies the
  level of long memory in each of the m channels.}

  \item{blur}{A string that specifies the blurring regime
  across all channels. Choices include: "direct", "smooth"
  and "box.car".}

  \item{j0}{The coarsest resolution level for the wavelet
  expansion.}

  \item{j1}{The finest resolution level for the wavelet
  expansion. If unspecified, the function will compute all
  thresholds up to the maximum possible resolution level at
  j1 = log2(n) - 1.}

  \item{thresh}{A vector of scale level thresholds to use
  in the wavelet thresholded estimator of the true signal.
  It should have enough elements to construct the required
  expansion with all resolutions. Namely, have j1 - j0 + 2
  elements. If a single element is input, it is repeated to
  be the universal threshold across all resolutions.}

  \item{eta}{The smoothing parameter. The default level is
  2*sqrt(alpha_{l_*}).}

  \item{deg}{The degree of the auxiliary polynomial for the
  Meyer wavelet.}
}
\description{
Returns a mWaveD object that contains all the required
information for the multichannel analysis.
}

