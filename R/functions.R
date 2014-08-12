# 
# Master file to control the output functions and their documentation.
# Copyright Justin R Wishart 2014
# Questions or comments?
# email: j.wishart@unsw.edu.au

#' @name mwaved
#' @title Multichannel wavelet deconvolution with long memory using mwaved.
#'
#' @description \code{mwaved} computes the Wavelet deconvolution estimate of a common signal present in multiple channels that have possible different levels of blur and additive error. More information about each function can be found in its help documentation.
#'
#' @details \code{mwaved} uses the WaveD wavelet deconvolution paradigm and is based on the \code{waved} R-package given by Raimondo and Stewart (2007). It generalises the approach by allowing a multichannel signal instead of a single channel signal and allows long memory errors within each channel (independent of each channel). See Kulik, Sapatinas and Wishart (2014) for theoretical results and a short numerical investigation. The \code{mwaved} package also uses the external C FFTW library described in Frigo and Johnson, (2005) to dramatically increase the speed of the computations. Detailed information and instructions for implementation are available at \url{http://www.fftw.org}.
#' 
#' @import Rcpp
#' @references 
#' Frigo, M and Johnson, S.G. (2005) \emph{The design and implementation of FFTW3}, Proceedings of the IEEE \bold{93}, 216--231.
#' \url{http://dx.doi.org/10.1109/JPROC.2004.840301}
#' 
#' Kulik, R., Sapatinas, T. and Wishart, J.R. (2014) \emph{Multichannel wavelet deconvolution with long range dependence. Upper bounds on the L_p risk}  Appl. Comput. Harmon. Anal. (to appear in).
#' \url{http://dx.doi.org/10.1016/j.acha.2014.04.004}
#' 
#' Raimondo, M. and Stewart, M. (2007) \emph{The WaveD Transform in R: Performs Fast Translation-Invariant Wavelet Deconvolution}, Journal of Statistical Software \bold{21}, 1--28
#' \url{http://www.jstatsoft.org/v21/i02}
#' @docType package
NULL
#' @name directBlur
#' @title Direct kernel
#' @description Creates appropriately sized blur matrix for the special case when no blurring is apparent.
#' @details Function creates a matrix of dimension n by m which contains appropriate entries for the case when a direct multichannel signal is observed. That is, no blurring operator is apparent. This is the default argument for the blurring matrix to all the multichannel functions in the \code{mwaved} package.
#' @param G.dim A numeric vector with two elements specifying the values of n and m respectively
#' @export
directBlur <- function(G.dim){
  y <- matrix(0, G.dim[1], G.dim[2])
  y[1,] <- 1 
  y
}

# Function to check resolution levels are sane
feasibleResolutions <- function(n, j0, j1){
  J <- log2(n) - 1
  if (j0 < 1) {
    warning("j0 must be a positive integer, setting to default j0 = 3.")
    j0 <- 3
  }
  
  if (!is.na(j1)) {
    if (j1 > J ) {
      warning(paste('Specified j1 = ', j1, ' is too large. j1 set to ', J, sep = ""))
      j1 <- J
    } else {
      if (j1 < j0) {
        warning(paste('j1 cannot be smaller than j0 = ', j0,'. j1 set to ', j0, sep = ""))
        j1 <- j0
      }
    }  
  }
  
  return(list(j0 = j0,j1 = j1))
}

# Check the dimensions of alpha are sane
feasibleAlpha <- function(m, alpha){
  ma <- length(alpha)
  if (m != ma) {
    if ( ma == 1) {
      alpha <- rep(alpha, m)
    } else {
      if ( m > ma ) {
        warning("Dimension mismatch: Length of alpha too short (not the same as the number of columns of Y). Last element of alpha repeated for the remaining channels.")
        alpha <- c(alpha, rep(rev(alpha)[1], m - ma))  
      } else {
        warning("Dimension mismatch: Length of alpha too long (longer than the number of columns of Y). Extra values of alpha ignored.")
        alpha <- alpha[1:m]
      }
    }
  }
  return(alpha)
}

# Check signal length is long enough to avoid computation errors in C code
feasibleLength <- function(n){
  if (n < 16) {
    stop("Signal length too small. Need at least 16 observations.")
  }
  J <- log2(n)
  if ((J %% 1) > 0) {
    warning("Input signal has length, n, which is not a power of 2.")
  }
}

# Check blur identification string 
feasibleBlur <- function(blur){
  if (blur != 'direct') {
    if (blur != 'smooth') {
      if (blur != 'box.car') {
        stop("unrecognised blur input. Please choose 'direct', 'smooth' or 'box.car'.")
      }
    }
  }
}

feasibleShrinkage <- function(shrinkType){
  if (shrinkType != 'hard') {
    if (shrinkType != 'soft') {
      if (shrinkType != 'garrote') {
        stop("unrecognised shrinkType input. Please choose 'soft', 'hard' or 'garrote'.")
      }
    }
  }
}

#' @title Noise level estimation among multichannel signal
#' 
#' @description Estimates the level of noise (standard deviation of noise) in each channel.
#' 
#' @param Y An input signal either an n by m matrix containing the multichannel signal to be analysed or single vector of n elements for the single channel. In the multichannel case, each of the m columns represents a channel of n observations.
#' @param deg The degree of the auxiliary polynomial used in the Meyer wavelet.
#' 
#' @details This function estimates the noise present in each channel by computing the Meyer wavelet transform of each channel of the multichannel signal. In particular, the noise level is computed by using the Median Absolute Deviation (MAD) of the wavelet coefficients at the highest possible resolution at J = floor(log2(n) - 1).
#' 
#' @return A numeric vector of estimates of the standard deviation of the noise in each of the m channels.
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
  n <- dim(Y)[1]
  feasibleLength(n)
  .Call(mwaved_multiSigma, Y, deg)
}

#' @title Meyer wavelet projection given a set of wavelet coefficients
#' 
#' @description Reconstructs a function using wavelet coefficients (\code{waveletCoef} object) as input.
#' 
#' @param beta A \code{waveletCoef} object that contains a vector of wavelet coefficients and the coarsest resolution level, \code{j0} to create the required output function expansion.
#' @param j1 The finest resolution to be used in the projection (specifies which resolution that the wavelet expansion is truncated).
#' 
#' @details Function that takes an input of wavelet coefficients in the form of a \code{waveletCoef} object (see \code{\link{multiCoef}} for details) and optionally a desired maximum resolution level, \code{j1}, to create an inhomogeneous wavelet expansion starting from resolution \code{j0} up to resolution \code{j1}. 
#' Namely, it creates the wavelet expansion, \deqn{\sum_{k = 0}^{2^{j_0} - 1} \beta_k \phi_{j_0,k} + \sum_{j = j_0}^{j_1} \sum_{k = 0}^{2^j - 1} \beta_{j,k} \psi_{j,k}.}
#' where \eqn{(\phi,\psi)} denote the father and mother periodised Meyer wavelet functions and \eqn{\beta_{j,k}} denotes the mother wavelet coefficient at resolution j and location k and \eqn{\beta_{k}} denotes the father wavelet coefficients at resolution \eqn{j=j0} and location \eqn{k}. 
#' The coefficients \code{beta} need to be ordered so that the first \eqn{2^\code{j0}} elements correspond to father wavelet coefficients at resolution \eqn{j=\code{j0}} and the remaining elements correspond to the mother wavelet coefficients from resolution \eqn{j=\code{j0}} to \eqn{j = log_2 n - 1}. If the maximum resolution level \emph{j1} is not specified, the full wavelet expansion will be given. 

#' @return A numeric vector of size n giving the wavelet function expansion.
#' 
#' @seealso \code{\link{multiCoef}}
#' 
#' @examples
#' library(mwaved)
#' # Make a noiseless doppler function
#' n <- 2^8
#' x <- (1:n)/n
#' y <- makeDoppler(n)
#' # Determine the wavelet coefficients
#' beta <- multiCoef(y)
#' # plot three raw wavelet expansions truncating in each case at j1 = 3, 4 and 5 respectively
#' plot(x, y, type = 'l', main = 'Doppler and wavelet projections at three different truncations')
#' j0 <- 3
#' j1 <- 5
#' j <- j0:j1
#' lcols <- c(1, j - j0 + 2)
#' ltys <-  c(1, 1:length(j))
#' matlines(x, sapply(j, function(i) multiProj(beta, j1 = i)), type = 'l', col = lcols[-1])
#' legend("bottomright", legend = c("Signal", paste('j1 =', j)), col = lcols, lty =ltys)
#' @export
multiProj <- function(beta, j1 = log2(length(beta$coef)) - 1) {
  coefs <- beta$coef
  j0 <- beta$j0
  deg <- beta$deg
  n <- length(coefs)
  jvals <- feasibleResolutions(n, j0, j1)
  
  .Call('mwaved_multiProj', coefs, jvals$j0, jvals$j1, deg)
}

#' @title Resolution level thresholds for hard thresholded wavelet deconvolution estimator
#' 
#' @description Computes the estimated resolution level thresholds for the hard-thresholding wavelet deconvolution estimate of the desired signal in the multichannel signals.
#' 
#' @inheritParams multiSigma
#' @param G The input multichannel blur matrix/vector (needs to be the same dimension/length as the signal input which is a matrix or vector for the multichannel or single channel case respectively). This argument dictates the form of blur present in each of the channels.
#' @param alpha A numeric vector, with m elements, specifying the level of long memory for the noise process within each channel of the form alpha = 2 - 2H, where H is the Hurst parameter. If alpha is a single element, that same element is repeated across all required channels.
#' @param blur A character string describing which deconvolution regime is to be applied.\itemize{
#' \item 'direct' (default): No deconvolution required and direct signal observed. This is captured when G is an n by m matrix with all entries being zero except the first row where all elements are one. 
#' \item 'smooth': Indirect signal observed and the blurring kernel is of regular smooth blur type.
#' \item 'box.car': Indirect signal observed and the blurring kernel is of box car type.}
#' @param j0 The coarsest resolution level for the wavelet expansion.
#' @param j1 The finest resolution level for the wavelet expansion. If unspecified, the function will compute all thresholds up to the maximum possible resolution level at j1 = log2(n) - 1.
#' @param eta The smoothing parameter. The default level is \eqn{2\sqrt(\alpha^*)} where \eqn{\alpha^*} is an optimal level depending on the type of blur. (see Kulik, Sapatinas and Wishart (2014) for details and justification)
#' 
#' @details Given an input matrix of a multichannel signal (n rows and n columns) with m channels and n observations per channel, the function returns the required thresholds for the hard-thresholding estimator of the underlying function, f.
#' 
#' @examples 
#' library(mwaved)
#' # Simulate the multichannel doppler signal.
#' m <- 3
#' n <- 2^10
#' signal <- makeDoppler(n)
#' # Noise levels per channel
#' e <- rnorm(m*n)
#' # Create Gamma blur
#' shape <- seq(from = 0.5, to = 1, length = m)
#' scale <- rep(0.25,m)
#' G <- gammaBlur(n, shape = shape, scale = scale)
#' # Convolve the signal
#' X <- blurSignal(signal, G)
#' # Create error with custom signal to noise ratio
#' SNR <- c(10,15,20)
#' sigma <- sigmaSNR(X, SNR)
#' alpha <- c(0.75, 0.8, 1)
#' E <- multiNoise(n, sigma, alpha)
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' # Determine thresholds
#' thresh <- multiThresh(Y, G, blur = 'smooth')
#' 
#' @return A numeric vector of the resolution level thresholds for the hard-thresholding nonlinear wavelet estimator from the multichannel model.
#' @references
#' Kulik, R., Sapatinas, T. and Wishart, J.R. (2014) \emph{Multichannel wavelet deconvolution with long range dependence. Upper bounds on the L_p risk}  Appl. Comput. Harmon. Anal. (to appear in).
#' \url{http://dx.doi.org/10.1016/j.acha.2014.04.004}
#' @export
multiThresh <- function(Y, G = directBlur(dim(as.matrix(Y))), alpha = rep(1,dim(as.matrix(Y))[2]), 
                        blur = "direct", j0 = 3L, j1 = NA_integer_, eta = NA_real_, deg = 3L) {
  Y <- as.matrix(Y)
  G <- as.matrix(G)
  dimY <- dim(Y)
  if (any(dimY != dim(G))) {
    stop('Dimensions of Y and G do not match')
  }
  feasibleLength(dimY[1])
  jvals <- feasibleResolutions(dimY[1], j0, j1)
  alpha <- feasibleAlpha(dimY[2], alpha)
  feasibleBlur(blur)
  .Call('mwaved_multiThresh', Y, G, alpha, blur, jvals$j0, jvals$j1, eta, deg)
}

#' @title Wavelet coefficient estimation from a multichannel signal
#' 
#' @description Estimates the wavelet coefficients for the underlying signal of interest embedded in the noisy multichannel deconvolution model. 
#' 
#' @inheritParams multiThresh
#' @param thresh A numeric vector of resolution level thresholds to use in the wavelet thresholded estimator of the true signal. It should have enough elements to construct the required expansion with all resolutions. That is, have \code{j1} - \code{j0} + 2 elements. If a single element is input, it is repeated to be the universal threshold across all resolutions.
#' 
#' @details Returns an object of type \emph{waveletCoef} which is a list including the following three objects \itemize{
#' \item \code{coef} a numeric vector of size n giving the estimated wavelet coefficients for the signal of interest
#' \item \code{j0} an integer that specifies the initial coarse resolution for the inhomogeneous wavelet expansion.
#' \item \code{deg} an integer that specifies the degree of the Meyer wavelet used in the estimation of the wavelet coefficients.
#' }
#' 
#' @examples
#' library(mwaved)
#' # Simulate the multichannel doppler signal.
#' m <- 3
#' n <- 2^10
#' signal <- makeDoppler(n)
#' # Noise levels per channel
#' e <- rnorm(m*n)
#' # Create Gamma blur
#' shape <- seq(from = 0.5, to = 1, length = m)
#' scale <- rep(0.25,m)
#' G <- gammaBlur(n, shape = shape, scale = scale)
#' # Convolve the signal
#' X <- blurSignal(signal, G)
#' # Create error with custom signal to noise ratio
#' SNR <- c(10,15,20)
#' sigma <- sigmaSNR(X, SNR)
#' if (require(fracdiff)) {
#'   alpha <- c(0.75, 0.8, 1)
#' } else {
#'   alpha <- rep(1, m)
#' }
#' E <- multiNoise(n, sigma, alpha)
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' plot(signal, type='l', lty='dashed', main='dashed: True signal, solid: multichannel signals')
#' matlines(Y)
#' # Estimate the wavelet coefficients
#' estimatedCoefs <- multiCoef(Y, G, alpha = alpha, blur = 'smooth')
#' plot(estimatedCoefs)
#' # Compute true wavelet coefficients
#' trueCoefs <- multiCoef(signal)
#' plot(trueCoefs)
#' @export
multiCoef <- function(Y, G = directBlur(dim(as.matrix(Y))), alpha = rep(1,dim(as.matrix(Y))[2]),
                      blur = "direct", j0 = 3L, j1 = NA_integer_, thresh = multiThresh(as.matrix(Y),
                      G = G, alpha = alpha, blur = blur, j0 = j0, j1 = j1, eta = eta, deg = 3L), 
                      eta = NA_real_, deg = 3L) {
  Y <- as.matrix(Y)
  G <- as.matrix(G)
  dimY <- dim(Y)
  if (any(dimY != dim(G))) {
    stop('Dimensions of Y and G do not match')
  }
  feasibleLength(dimY[1])
  jvals <- feasibleResolutions(dimY[1], j0, j1)
  alpha <- feasibleAlpha(dimY[2], alpha)
  
  blur <- tolower(blur)
  feasibleBlur(blur)
  # Pass to C code
  .Call('mwaved_multiCoef', Y, G, alpha, jvals$j0, jvals$j1, thresh, deg)
}

#' @title Apply thresholding regime to a set of wavelet coefficients
#' 
#' @description Applies a resolution level thresholding technique to a set of wavelet coefficients,
#' embedded in a wavelet coefficient object.
#'
#' @param beta A \code{waveletCoef} object.
#' @param thresh A numeric vector containing the thresholds to be applied to the coefficients 
#' at each resolution.
#' @param shrinkType A character string that specifies which thresholding regime to use. 
#' Available choices are the 'hard', 'soft' or 'garrote'.
#' 
#' @details Applies one of three specified wavelet thresholding regimes to a \code{waveletCoef} object (wavelet 
#' coefficient object created by \code{\link{multiCoef}}). If \code{thresh} is not specified, no thresholding 
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
#' @examples
#' library(mwaved)
#' # Simulate the multichannel doppler signal.
#' m <- 3
#' n <- 2^10
#' signal <- makeDoppler(n)
#' # Noise levels per channel
#' e <- rnorm(m * n)
#' # Create Gamma blur
#' shape <- seq(from = 0.5, to = 1, length = m)
#' scale <- rep(0.25, m)
#' G <- gammaBlur(n, shape = shape, scale = scale)
#' # Convolve the signal
#' X <- blurSignal(signal, G)
#' # Create error with custom signal to noise ratio
#' SNR <- c(10, 15, 20)
#' sigma <- sigmaSNR(X, SNR)
#' if (require(fracdiff)){
#'   alpha <- c(0.75, 0.8, 1)
#' } else {
#'   alpha <- rep(1, m)
#' }
#' E <- multiNoise(n, sigma, alpha)
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' # Determine thresholds
#' thresh <- multiThresh(Y, G, blur = 'smooth')
#' beta <- multiCoef(Y, G)
#' betaShrunk <- waveletThresh(beta, thresh)
#' plot(beta, betaShrunk)
#' @export 
waveletThresh <- function(beta, thresh = 0, shrinkType = 'hard'){
  if (class(beta) != "waveletCoef") {
    stop('beta needs to be a waveletCoef object.')
  }
  nthr <- length(thresh)
  
  req <- log2(length(beta$coef)) - beta$j0
  j1 <- beta$j0j0 + req
  # convert to lower case to avoid trivial issues
  shrinkType <- tolower(shrinkType)
  feasibleShrinkage(shrinkType)
  if (any(thresh < 0)) {
    stop("Input thresholds must all be non-negative.")
  }
  
  if (nthr == 1 ) {
    if (thresh == 0) {
      thresh <- rep(0, req)  
    } else {
      warning("thresh input vector only has one element. Universal threshold applied on all resolutions.")
      thresh <- rep(thresh, req)
    }
    
  }
  
  if (nthr < req ){
    warning("Thresh input length too small, last element repeated in higher resolutions.")
    thresh <- c(thresh, rep(thresh[nthr], req - nthr))
  }
  return(.Call('mwaved_waveletThresh', beta$coef, thresh, shrinkType, beta$j0, beta$deg))
}


#' @title Full mWaveD analysis
#' 
#' @description Returns a mWaveD object that contains all the required information for the multichannel analysis.
#' @inheritParams multiCoef
#' @param shrinkType A character string that specifies which thresholding regime to use. 
#' Available choices are the 'hard', 'soft' or 'garrote'.
#' @examples
#' library(mwaved)
#' # Simulate the multichannel doppler signal.
#' m <- 3
#' n <- 2^10
#' signal <- makeDoppler(n)
#' # Noise levels per channel
#' e <- rnorm(m * n)
#' # Create Gamma blur
#' shape <- seq(from = 0.5, to = 1, length = m)
#' scale <- rep(0.25, m)
#' G <- gammaBlur(n, shape = shape, scale = scale)
#' # Convolve the signal
#' X <- blurSignal(signal, G)
#' # Create error with custom signal to noise ratio
#' SNR <- c(10, 15, 20)
#' sigma <- sigmaSNR(X, SNR)
#' if (require(fracdiff)) {
#'   alpha <- c(0.75, 0.8, 1)
#' } else {
#'   alpha <- rep(1, m)
#' }
#' E <- multiNoise(n, sigma, alpha)
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' # Compute mWaveD object
#' mWaveDObj <- multiWaveD(Y, G = G, alpha = alpha, blur = 'smooth')
#' plot(mWaveDObj)
#' summary(mWaveDObj)
#' 
#' @seealso \code{\link{plot.mWaveD}} and \code{\link{summary.mWaveD}}
#' 
#' @export
multiWaveD <- function(Y, G = directBlur(dim(as.matrix(Y))), alpha = rep(1,dim(as.matrix(Y))[2]),
                       j0 = 3L, j1 = NA_integer_, blur = "direct", thresh = as.numeric(c()), 
                       eta = NA_real_, shrinkType = "hard", deg = 3L) {
  Y <- as.matrix(Y)
  G <- as.matrix(G)
  dimY <- dim(Y)
  if (any(dimY != dim(G))) {
    stop('Dimensions of Y and G do not match')
  }
  feasibleLength(dimY[1])
  jvals <- feasibleResolutions(dimY[1], j0, j1)
  alpha <- feasibleAlpha(dimY[2], alpha)
  blur <- tolower(blur)
  feasibleBlur(blur)
  shrinkType <- tolower(shrinkType)
  feasibleShrinkage(shrinkType)
  # Pass to C code
  return(.Call('mwaved_multiWaveD', Y, G, alpha, jvals$j0, jvals$j1, blur, thresh, eta, shrinkType, deg))
}

#' @title Wavelet deconvolution signal estimate from the noisy multichannel convoluted signal
#' 
#' @description Estimates the underlying signal of interest from a multichannel noisy deconvolution model.
#' 
#' @inheritParams multiWaveD
#' @param sigma A numeric vector with m elements that specifies the level of noise (standard deviation) in each channel. The default method uses the Median Absolute Deviation of wavelet coefficients in the finest resolution (see \code{\link{multiSigma}}) for details.
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
#' x <- (1:n)/n
#' signal <- makeDoppler(n)
#' # Noise levels per channel
#' e <- rnorm(m * n)
#' # Create Gamma blur
#' shape <- seq(from = 0.5, to = 1, length = m)
#' scale <- rep(0.25, m)
#' G <- gammaBlur(n, shape = shape, scale = scale)
#' # Convolve the signal
#' X <- blurSignal(signal, G)
#' # Create error with custom signal to noise ratio
#' SNR <- c(10, 15, 20)
#' sigma <- sigmaSNR(X, SNR)
#' if (require(fracdiff)) {
#'   alpha <- c(0.75, 0.8, 1)
#' } else {
#'   alpha <- rep(1, m)
#' }
#' E <- multiNoise(n, sigma, alpha)
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' # Estimate the underlying doppler signal
#' dopplerEstimate <- multiEstimate(Y, G = G, alpha = rep(1, m), blur = 'smooth')
#' # Plot the result and compare with truth
#' par(mfrow=c(2, 1))
#' matplot(x, Y, type = 'l', main = 'Noisy multichannel signal')
#' plot(x, signal, type = 'l', lty = 2, main = 'True Doppler signal and estimate', col = 'red')
#' lines(x, dopplerEstimate)
#' 
#' @export
multiEstimate <- function(Y, G = directBlur(dim(as.matrix(Y))), alpha = rep(1,dim(as.matrix(Y))[2]), 
                          blur = "direct", sigma = as.numeric(c()), j0 = 3L, j1 = NA_integer_, 
                          eta = NA_real_, thresh = multiThresh(as.matrix(Y), G = G, alpha = alpha,
                                                               blur = blur, j0 = j0, j1 = j1, eta = eta, deg = 3L) , shrinkType = "hard", deg = 3L) {
  Y <- as.matrix(Y)
  G <- as.matrix(G)
  dimY <- dim(Y)
  if (any(dimY != dim(G))) {
    stop('Dimensions of Y and G do not match')
  }
  feasibleLength(dimY[1])
  jvals <- feasibleResolutions(dimY[1], j0, j1)
  alpha <- feasibleAlpha(dimY[2], alpha)
  blur <- tolower(blur)
  feasibleBlur(blur)
  shrinkType <- tolower(shrinkType)
  feasibleShrinkage(shrinkType)
  # Pass to C code
  .Call('mwaved_multiEstimate', Y, G, alpha, blur, sigma, j0, j1, eta, thresh, shrinkType, deg)
}
