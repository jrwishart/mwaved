\name{multiCoef}
\alias{multiCoef}
\title{multiCoef}
\usage{
multiCoef(Y, G = directG(dim(as.matrix(Y))), alpha = rep(1,
  dim(as.matrix(Y))[2]), blur = "direct", j0 = 3L, j1 = NA_integer_,
  thresh = multiThresh(as.matrix(Y), G = G, alpha = alpha, blur = blur, j0 =
  j0, j1 = j1, eta = eta, deg = 3L), eta = NA_real_, deg = 3L)
}
\arguments{
  \item{Y}{The input multichannel signal in the form of an
  \emph{n} by \emph{m} matrix to denote the \emph{m}
  separate channels with \emph{n} observations in each
  channel. If a vector of size \emph{n} is given, it is
  parsed as a single channel signal with n elements.}

  \item{G}{The input multichannel blur matrix of size
  \emph{n} by \emph{m} (needs to be the same size as the
  signal input). This matrix dictates the form of blur
  present in each of the channels.}

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
Estimates the wavelet coefficiets for the underlying signal
of interest embedded in the noisy multichannel
deconvolution model.
}
\details{
Returns an object of type \emph{waveletCoef} including a
numeric vector of size n the estimated wavelet coefficients
for the signal of interest embedded in the noisy
multichannel deconvolution model and an integer, \emph{j0},
that specifies the initial resolution for the coarse
expansion.
}
\examples{
library(mwaved)
# Simulate the multichannel doppler signal.
m <- 3
n <- 2^10
t <- (1:n)/n
signal <- (sqrt(t * (1 - t))) * sin((2 * pi * 1.05)/(t + 0.05))
signal <- signal * 2.4
signal <- matrix(rep(signal,m), nrow=n, ncol=m)
# Noise levels per channel
e <- rnorm(m*n)
# Create Gamma blur
alp <- c(0.5,0.75,1)
beta <- rep(0.25,3)
G <- sapply((1:n)/n, dgamma, shape = alp, scale = beta)
G <- t(G)
# Normalise the blur
Gmax <- matrix(rep(apply(G, 2, max), dim(G)[1]), ncol = m, nrow = n, byrow = T)
G <- G/Gmax
Gsum <- matrix(rep(apply(G, 2, sum), dim(G)[1]), ncol = m, nrow = n, byrow = T)
G <- G/Gsum
# Convolve the signal
X <- Re(mvfft(mvfft(G) * mvfft(signal), inverse = TRUE))/n
# Create error with custom signal to noise ratio
SNR <- c(10,15,20)
sigma <- sqrt(mean(X^2)) * 10^( -SNR/20 )
E <- matrix(rnorm(n*m, sd = sigma), ncol = m, nrow = n, byrow = T)
# Create noisy & blurred multichannel signal
Y <- X + E
# Estimate the wavelet coefficients
dopplerEstimate <- multiEstimate(Y, G = G, alpha = rep(1,3), blur = 'smooth')
# Plot the result and compare with truth
par(mfrow=c(2,1))
matplot(t, Y, type = 'l', main = 'Noisy multichannel signal')
plot(t, y, type = 'l', lty = 2, main = 'True Doppler signal and estimate', col = 'red')
lines(t, dopplerEstimate)
}

