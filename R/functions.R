# 
# Master file to control the output functions and their documentation.
# Copyright Justin R Wishart 2014
# Questions or comments?
# email: j.wishart@unsw.edu.au
#

# Function to create default direct inverse kernel matrix (no blur)
directG <- function(G.dim){
  y <- matrix(0,G.dim[1],G.dim[2])
  y[1,] <- 1 
  y
}

# Function to check resolution levels are feasible
feasibleResolutions <- function(n, j0, j1){
  J <- log2(n) - 1
  if ( j0 < 1 ){
    warning("j0 cannot be negative, setting j0 to default, j0 = 3")
    j0 <- 3
  }
  
  if ( !is.na(j1) ){
    if( j1 > J ){
      warning(paste('Specified j1 = ', j1, ' is too large. j1 set to ', J, sep = ""))
      j1 <- J
    } else {
      if( j1 < j0){
        warning(paste('j1 cannot be smaller than j0 = ', j0,'. j1 set to ', j0, sep = ""))
        j1 <- j0
      }
    }  
  }
  
  return(list(j0 = j0,j1 = j1))
}

#' @title multiSigma
#' 
#' @description Estimates the level of noise (standard deviation of noise) in each channel.
#' 
#' @param Y An input signal either an n by M matrix containing the multichannel signal to be analysed (or single vector of n elements). Each of the M columns represents a channel of n observations.
#' @param deg The degree of the auxiliary polynomial used in the Meyer wavelet.
#' 
#' @details This function estimates the noise present in each channel by computing the Meyer wavelet transform of each channel of the multichannel signal. In particular, the noise level is computed by using the Median Absolute Deviation (MAD) of the wavelet coefficients at the highest possible scale at J = floor(log2(n) - 1).
#' 
#' @return a vector of estimates containing the level of sigma in each of the M channels.
#'
#' @examples
#' library(mwaved)
#' # Simulate matrix of Gaussian variables with three different noise levels
#' sig <- c(0.25, 0.5, 1.25)
#' Y <- sapply(1:3, function(i) sig[i]*rnorm(1024))
#' # Estimate the noise levels
#' multiSigma(Y, deg = 3)
#'
#' @export 
multiSigma <- function(Y, deg = 3L){
  Y <- as.matrix(Y)
  .Call(mwaved_multiSigma, Y, deg)
}

#' @title multiProj
#' 
#' @description Reconstructs a function using wavelet coefficients as input.
#' 
#' @param beta A vector of wavelet coefficients to be projected to create the required output function expansion.
#' @param j0 The coarsest resolution to be used in the projection. The first 2^{j0} elements in the beta vector are used as the coarse expansion.
#' @param j1 The finest resolution to be used in the projection (specifies which resolution that the wavelet expansion is truncated).
#' @param deg The degree of the auxiliary polynomial used in the Meyer wavelet.
#' 
#' @details Function that takes an input of wavelet coefficients, \emph{beta} of length n and optionally a desired coarse resolution level and maximum resolution level, \emph{j1}, to create an inhomogeneous wavelet expansion starting from resolution \emph{j0} up to resolution \emph{j1}. If the maximum resolution level \emph{j1} is not specified, the full wavelet expansion will be given. Namely, for \emph{j0} <= j <= \emph{j1} and 0 <= k <= 2^j-1, \deqn{\sum_j \sum_k beta_{j,k} \psi_{j,k}.}
#' 
#' @return A numeric vector of size n giving the wavelet function expansion.
#' 
#' @seealso multiCoef
#' 
#' @examples
#' library(mwaved)
#' # Make a noiseless doppler function
#' n <- 2^8
#' t <- (1:n)/n
#' y <- (sqrt(t * (1 - t))) * sin((2 * pi * 1.05)/(t + 0.05))
#' y <- y * 2.4
#' # Determine the wavelet coefficients
#' beta <- multiCoef(y)
#' # plot the Multi-Resolution Analysis from j0 = 3 to j1 = 5
#' plot(beta, lowest = 3, highest = 5)
#' 
#' @export
multiProj <- function(beta, j0 = 3L, j1 = log2(length(beta)) - 1, deg = 3L) {
  n <- length(beta)
  jvals <- feasibleResolutions(n, j0, j1)
  
  .Call('mwaved_multiProj', beta, jvals$j0, jvals$j1, deg)
}

#' @title multiThresh
#' 
#' @description Computes the resolution level thresholds for the hard-thresholding estimate of the desired funnction f.
#' 
#' @param Y The input multichannel signal in the form of an n by m matrix to denote the m separate channels with n observations in each channel. If a vector of size n is given, it is parsed as a single channel signal with n elements.
#' @param G The input multichannel blur matrix of size n by m (needs to be the same size as the signal input). This matrix dictates the form of blur present in each of the channels.
#' @param alpha A numeric vector, with m elements, specifying the level of long memory for the noise process within each channel of the form alpha = 2 - 2H, where H is the Hurst parameter. If alpha is a single element, that same element is repeated across all required channels.
#' @param j0 The coarsest resolution level for the wavelet expansion.
#' @param j1 The finest resolution level for the wavelet expansion. If unspecified, the function will compute all thresholds up to the maximum possible resolution level at j1 = log2(n) - 1.
#' @param eta The smoothing parameter. The default level is 2*sqrt(alpha_{l_*}).
#' @param deg The degree of the auxiliary polynomial for the Meyer wavelet.
#' 
#' @details Given an input matrix of a multichannel signal (n rows and n columns) with m channels and n observations per channel, the function returns the required thresholds for the hard-thresholding estimator of the underlying function, f.
#' 
#' @examples 
#' library(mwaved)
#' # Simulate the multichannel doppler signal.
#' m <- 3
#' n <- 2^10
#' t <- (1:n)/n
#' signal <- (sqrt(t * (1 - t))) * sin((2 * pi * 1.05)/(t + 0.05))
#' signal <- signal * 2.4
#' signal <- matrix(rep(signal,m), nrow=n, ncol=m)
#' # Noise levels per channel
#' e <- rnorm(m*n)
#' # Create Gamma blur
#' alp <- c(0.5,0.75,1)
#' beta <- rep(0.25,3)
#' G <- sapply((1:n)/n, dgamma, shape = alp, scale = beta)
#' G <- t(G)
#' # Normalise the blur
#' Gmax <- matrix(rep(apply(G, 2, max), dim(G)[1]), ncol = m, nrow = n, byrow = T)
#' G <- G/Gmax
#' Gsum <- matrix(rep(apply(G, 2, sum), dim(G)[1]), ncol = m, nrow = n, byrow = T)
#' G <- G/Gsum
#' # Convolve the signal
#' X <- Re(mvfft(mvfft(G) * mvfft(signal), inverse = TRUE))/n
#' # Create error with custom signal to noise ratio
#' SNR <- c(10,15,20)
#' sigma <- sqrt(mean(X^2)) * 10^( -SNR/20 )
#' E <- matrix(rnorm(n*m, sd = sigma), ncol = m, nrow = n, byrow = T)
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' thresh <- multiThresh(Y, G, blur = 'smooth')
#' 
#' @return A numeric vector of the resolution level thresholds for the hard-thresholding nonlinear wavelet estimator from the multichannel model.
#' 
#' @export
multiThresh <- function(Y, G = directG(dim(as.matrix(Y))), alpha = rep(1,dim(as.matrix(Y))[2]), blur = "direct", j0 = 3L, j1 = NA_integer_, eta = NA_real_, deg = 3L) {
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  jvals <- feasibleResolutions(n, j0, j1)
  
  .Call('mwaved_multiThresh', Y, G, alpha, blur, jvals$j0, jvals$j1, eta, deg)
}

#' @title multiEstimate
#' 
#' @description Estimates the underlying signal of interest from a multichannel noisy deconvolution model.
#' 
#' @param Y The input multichannel signal in the form of an n by m matrix to denote the m separate channels with n observations in each channel. If a vector of size n is given, it is parsed as a single channel signal with n elements.
#' @param G The input multichannel blur matrix of size n by m (needs to be the same dimensions and size as the signal input). This matrix dictates the form of blur present in each of the channels.
#' @param alpha A numeric vector, with m elements, specifying the level of long memory for the noise process within each channel of the form alpha = 2 - 2H, where H is the Hurst parameter. If alpha is a single element, that same element is repeated across all required channels.
#' @param j0 The coarsest resolution level for the wavelet expansion.
#' @param j1 The finest resolution level for the wavelet expansion. If unspecified, the function will compute all thresholds up to the maximum possible resolution level at j1 = log2(n) - 1.
#' @param eta The smoothing parameter. The default level is 2*sqrt(alpha_{l_*}).
#' @param deg The degree of the auxiliary polynomial for the Meyer wavelet.
#' 
#' @details Function requires input of a noisy multichannel signal matrix, Y, which contains the information for each channel in each of the m columns. Optional inputs are a matrix, G, the same dimension as Y, that gives the multichannel blur information.
#' 
#' @return A numeric vector of the estimate of the underlying signal of interest.
#' 
#' @examples 
#' library(mwaved)
#' # Simulate the multichannel doppler signal.
#' m <- 3
#' n <- 2^10
#' t <- (1:n)/n
#' signal <- (sqrt(t * (1 - t))) * sin((2 * pi * 1.05)/(t + 0.05))
#' signal <- signal * 2.4
#' signalMatrix <- matrix(rep(signal,m), nrow=n, ncol=m)
#' # Noise levels per channel
#' e <- rnorm(m*n)
#' # Create Gamma blur
#' alp <- c(0.5,0.75,1)
#' beta <- rep(0.25,3)
#' G <- sapply((1:n)/n, dgamma, shape = alp, scale = beta)
#' G <- t(G)
#' # Normalise the blur
#' Gmax <- matrix(rep(apply(G, 2, max), dim(G)[1]), ncol = m, nrow = n, byrow = T)
#' G <- G/Gmax
#' Gsum <- matrix(rep(apply(G, 2, sum), dim(G)[1]), ncol = m, nrow = n, byrow = T)
#' G <- G/Gsum
#' # Convolve the signal
#' X <- Re(mvfft(mvfft(G) * mvfft(signal), inverse = TRUE))/n
#' # Create error with custom signal to noise ratio
#' SNR <- c(10,15,20)
#' sigma <- sqrt(mean(X^2)) * 10^( -SNR/20 )
#' E <- matrix(rnorm(n*m, sd = sigma), ncol = m, nrow = n, byrow = T)
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' # Estimate the underlying doppler signal
#' dopplerEstimate <- multiEstimate(Y, G = G, alpha = rep(1,3), blur = 'smooth')
#' # Plot the result and compare with truth
#' par(mfrow=c(2,1))
#' matplot(t, Y, type = 'l', main = 'Noisy multichannel signal')
#' plot(t, signal, type = 'l', lty = 2, main = 'True Doppler signal and estimate', col = 'red')
#' lines(t, dopplerEstimate)
#' 
#' @export

multiEstimate <- function(Y, G = directG(dim(as.matrix(Y))), alpha = rep(1,dim(as.matrix(Y))[2]), blur = "direct", sigma = multiSigma(as.matrix(Y), deg = 3L), j0 = 3L, j1 = NA_integer_, 
                          eta = NA_real_, thresh = multiThresh(as.matrix(Y), G = G, alpha = alpha, blur = blur, j0 = j0, j1 = j1, eta = eta, deg = 3L) , shrinkage = "Hard", deg = 3L) {
  Y <- as.matrix(Y)
  .Call('mwaved_multiEstimate', Y, G, alpha, blur, sigma, j0, j1, eta, thresh, shrinkage, deg)
}

#' @title multiCoef
#' 
#' @description Estimates the wavelet coefficiets for the underlying signal of interest embedded in the noisy multichannel deconvolution model. 
#' 
#' @param Y The input multichannel signal in the form of an \emph{n} by \emph{m} matrix to denote the \emph{m} separate channels with \emph{n} observations in each channel. If a vector of size \emph{n} is given, it is parsed as a single channel signal with n elements.
#' @param G The input multichannel blur matrix of size \emph{n} by \emph{m} (needs to be the same size as the signal input). This matrix dictates the form of blur present in each of the channels.
#' @param alpha A vector of length m that specifies the level of long memory in each of the m channels. 
#' @param blur A string that specifies the blurring regime across all channels. Choices include: "direct", "smooth" and "box.car".
#' @param j0 The coarsest resolution level for the wavelet expansion.
#' @param j1 The finest resolution level for the wavelet expansion. If unspecified, the function will compute all thresholds up to the maximum possible resolution level at j1 = log2(n) - 1.
#' @param thresh A vector of scale level thresholds to use in the wavelet thresholded estimator of the true signal. It should have enough elements to construct the required expansion with all resolutions. Namely, have j1 - j0 + 2 elements. If a single element is input, it is repeated to be the universal threshold across all resolutions.
#' @param eta The smoothing parameter. The default level is 2*sqrt(alpha_{l_*}).
#' @param deg The degree of the auxiliary polynomial for the Meyer wavelet.
#' 
#' @details Returns an object of type \emph{waveletCoef} including a numeric vector of size n the estimated wavelet coefficients for the signal of interest embedded in the noisy multichannel deconvolution model and an integer, \emph{j0}, that specifies the initial resolution for the coarse expansion. 
#' 
#' @examples
#' library(mwaved)
#' # Simulate the multichannel doppler signal.
#' m <- 3
#' n <- 2^10
#' t <- (1:n)/n
#' signal <- (sqrt(t * (1 - t))) * sin((2 * pi * 1.05)/(t + 0.05))
#' signal <- signal * 2.4
#' signal <- matrix(rep(signal,m), nrow=n, ncol=m)
#' # Noise levels per channel
#' e <- rnorm(m*n)
#' # Create Gamma blur
#' alp <- c(0.5,0.75,1)
#' beta <- rep(0.25,3)
#' G <- sapply((1:n)/n, dgamma, shape = alp, scale = beta)
#' G <- t(G)
#' # Normalise the blur
#' Gmax <- matrix(rep(apply(G, 2, max), dim(G)[1]), ncol = m, nrow = n, byrow = T)
#' G <- G/Gmax
#' Gsum <- matrix(rep(apply(G, 2, sum), dim(G)[1]), ncol = m, nrow = n, byrow = T)
#' G <- G/Gsum
#' # Convolve the signal
#' X <- Re(mvfft(mvfft(G) * mvfft(signal), inverse = TRUE))/n
#' # Create error with custom signal to noise ratio
#' SNR <- c(10,15,20)
#' sigma <- sqrt(mean(X^2)) * 10^( -SNR/20 )
#' E <- matrix(rnorm(n*m, sd = sigma), ncol = m, nrow = n, byrow = T)
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' # Estimate the wavelet coefficients
#' dopplerEstimate <- multiEstimate(Y, G = G, alpha = rep(1,3), blur = 'smooth')
#' # Plot the result and compare with truth
#' par(mfrow=c(2,1))
#' matplot(t, Y, type = 'l', main = 'Noisy multichannel signal')
#' plot(t, y, type = 'l', lty = 2, main = 'True Doppler signal and estimate', col = 'red')
#' lines(t, dopplerEstimate)
#' 
#' @export
multiCoef <- function(Y, G = directG(dim(as.matrix(Y))), alpha = rep(1,dim(as.matrix(Y))[2]), blur = "direct", j0 = 3L, j1 = NA_integer_, thresh = multiThresh(as.matrix(Y), G = G, alpha = alpha, blur = blur, j0 = j0, j1 = j1, eta = eta, deg = 3L), eta = NA_real_, deg = 3L) {
  Y <- as.matrix(Y)
  jvals <- feasibleResolutions(n, j0, j1)
  
  .Call('mwaved_multiCoef', Y, G, alpha, blur, jvals$j0, jvals$j1, thresh, eta, deg)
}

#' @title waveletThresh
#' 
#' @description Applies a resolution level thresholding technique to a set of wavelet coefficients,
#' embedded in a wavelet coefficient object.
#'
#' @param beta A waveletCoef object that includes a set of n wavelet coeffients and integer 
#' j0 specifying the coarsest resolution.
#' @param thresh A numeric vector containing the thresholds to be applied to the coefficients 
#' at each resolution.
#' @param shrinkType A character string that specifies which thresholding regime to use. 
#' Available choices are the 'hard', 'soft' or 'garrote'. 
#' 
#' @details Applies one of three specified wavelet thresholding regimes to a wavelet 
#' coefficient object (list with a vector of n wavelet coefficients and an integer 
#' specifying the coarse resolution level). If thresh is not specified, no thresholding 
#' is applied. The formulae applied for 'hard', 'soft' or 
#' 'garrote' are given by,\itemize{
#'  \item Hard: \eqn{
#'    \delta(x) = x 1(|x| > t)
#'  }
#'  \item Soft: \eqn{
#'    \delta(x) = (x - t) 1(x > t) + (x + t) 1(x > -t)
#'  }
#'  \item Garrote: \eqn{
#'    \delta(x) = (x - t^2/x) 1(|x| > t)
#'  }
#' }
#' where 1 represents the indicator function and \emph{t > 0} represents the threshold.
#' 
#' @export 
waveletThresh <- function(beta, thresh = 0, shrinkType = 'hard'){
  nthr <- length(thresh)
  n <- length(beta$coef)
  j0 <- beta$j0
  J <- log2(nthr) - 1
  req <- J - j0 + 1
  
  if( nthr < req && thresh != 0 ){
    if( ntrh == 1 ){
      warning("thresh input vector only has one element. Universal threshold applied on all resolutions.")
      thresh <- rep(thresh, req)
    } else {
      warning("length of thresh input too small, last element repeated in higher resolutions.")
      thresh <- c(thresh, rep(thresh[nthr], req - nthr))
    } 
  }
  return(.Call('mwaved_multiWaveD', beta, thresh, jvals$j0, jvals$j1))
}


#' @title multiWaveD
#' 
#' @description Returns a mWaveD object that contains all the required information for the multichannel analysis.
#' 
#' @param Y The input multichannel signal in the form of an n by m matrix to denote the m separate channels with n observations in each channel. If a vector of size n is given, it is parsed as a single channel signal with n elements.
#' @param G The input multichannel blur matrix of size n by m (needs to be the same size as the signal input). This matrix dictates the form of blur present in each of the channels.
#' @param alpha A vector of length m that specifies the level of long memory in each of the m channels. 
#' @param blur A string that specifies the blurring regime across all channels. Choices include: "direct", "smooth" and "box.car".
#' @param j0 The coarsest resolution level for the wavelet expansion.
#' @param j1 The finest resolution level for the wavelet expansion. If unspecified, the function will compute all thresholds up to the maximum possible resolution level at j1 = log2(n) - 1.
#' @param thresh A vector of scale level thresholds to use in the wavelet thresholded estimator of the true signal. It should have enough elements to construct the required expansion with all resolutions. Namely, have j1 - j0 + 2 elements. If a single element is input, it is repeated to be the universal threshold across all resolutions.
#' @param eta The smoothing parameter. The default level is 2*sqrt(alpha_{l_*}).
#' @param deg The degree of the auxiliary polynomial for the Meyer wavelet.
#' 
#' @export
multiWaveD <- function(Y, G = directG(dim(as.matrix(Y))), alpha = rep(1,dim(as.matrix(Y))[1]), j0 = 3L, j1 = NA_integer_, blur = "direct", thresh = as.numeric(c()), eta = NA_real_, 
                   shrinkage = "Hard", deg = 3L) {
  Y <- as.matrix(Y)
  jvals <- feasibleResolutions(n, j0, j1)
  
  return(.Call('mwaved_multiWaveD', Y, G, alpha, jvals$j0, jvals$j1, blur, thresh, eta, shrinkage, deg))
}
