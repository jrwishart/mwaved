#' @title Summary Output for the mWaveD object
#'
#' @param object A mWaveD object which is a list containing all the information for a multichannel 
#' deconvolution analysis produced by the \code{\link{multiWaveD}} function.
#' @param ... Arguments to be passed to methods.
#' @description Summarises the mWaveD object by giving
#' @return Text output given summary information of the input and output analysis including, \itemize{
#' \item Degree of Meyer wavelet used.
#' \item Number of observations, within each channel and number of channels present.
#' \item Resolution levels used (j0 to j1)
#' \item Blur type assumed in the analysis (direct, smooth or box.car)
#' \item Matrix summarising the noise levels in each channel (and Fourier decay information for the smooth case)
#' \item Summaries of the thresholding,
#' }
#' @seealso \code{\link{multiWaveD}}
#' @export
summary.mWaveD <- function(object, ...){
  n = length(object$estimate)
  t = (1:n)/n
  m = object$channels
  
  blurInfo <- object$blurInfo
  j0 <- blurInfo$j0
  j1 <- blurInfo$j1
  cat("Degree of Meyer wavelet =", object$degree, "  , Coarse resolution level j0 =", j0)
  cat("\n")
  cat("Sample size per channel = ", n, ", Maximum possible resolution level = ", log2(n) - 1, ".",sep='')
  cat("\n\n")
  blurInfo <- object$blurInfo
  
  if( object$blurType == "direct" ){
    cat("Number of channels: M =", m,"\n")
    cat('Blur type: ',object$blurType,'\n\n')
    cat("Estimated Channel information:\n\n")
    mat <- cbind(round(object$sigma,3), round(object$alpha,3), blurInfo$freq,rep(blurInfo$j1,m))
    colnames(mat) <- c("Sigma.hat", "Alpha", "Fourier freq cutoff", "Highest resolution")
    rownames(mat) <- paste("Channel: ", 1:m)
    print(mat)
  } else {
    # If Smooth blur is used, display the matrix of values
    if( object$blurType == "smooth"){
      cat("Number of channels: M =", m,"\n")
      cat('Blur type: ',object$blurType,'\n\n')
      cat("Estimated Channel Information:\n\n")
      mat <- cbind(round(object$sigma,3), round(object$alpha,3), blurInfo$freq,blurInfo$maxLevels)
      colnames(mat) <- c("Sigma.hat", "Alpha", "Fourier freq cutoff", "Highest resolution")
      rownames(mat) <- paste("Channel: ", 1:m)
      print(mat)
      cat("\n")
      cat("Estimated best channel = Channel", blurInfo$bestChannel)
      
    } else {
      if( object$blurType == "box.car"){
        cat("Number of channels: M =", m,"\n")
        cat('Blur type: ',object$blurType,'\n\n')
        cat("Estimated Channel Information:\n\n")
        mat <- cbind(round(object$sigma,3), round(object$alpha,3))
        colnames(mat) <- c("Sigma.hat", "Alpha")
        rownames(mat) <- paste("Channel: ", 1:m)
        print(mat)
      } else {
        warning('Unrecognised blur.type')
      }
    } 
  }
  cat("\n\n")
  cat("mWaveD optimal finest resolution level j1 =", j1)

  cat("\n\n")
  cat("Thresholding method:", object$shrinkType, "   Tuning parameter: eta =",object$eta,'\n\n')
  
  threshMatrix <- cbind(round(object$levelMax,4), round(object$thresh,4), object$percent)
  rownames(threshMatrix) <- paste("Level", j0:j1,":" )
  colnames(threshMatrix) <- c("Max|w|", "Threshold", "% Shrinkage" )
  print(threshMatrix)
} 

#' @title Multi-Resolution Analysis plot of wavelet coefficients
#'
#' @description Plots the wavelet coefficient object in the multiresolution analysis
#' 
#' @param x A list created by the waveletCoef.
#' @param ... Arguments to be passed to methods.
#' @param lowest Specifies the coarsest resolution to display in the Multi-resolution plot.
#' @param highest Specifies the finest resolution to display in the Multi-resolution plot.
#' @param coefTrim A numeric vector of trimmed wavelet coefficients to be overlayed on top of the plot for comparison with the 'x' wavelet coefficients. 
#' @param thickness An integer that specifies the thickness of the overlayed coefficients (coefTrim values) in the plot. Larger values increase the thickness.
#' @param descending A logical value to specify whether resolutions on the y-axis of the plot are increasing from top to bottom.
#' @param scaling A numeric value that acts as a graphical scaling parameter to rescale the wavelet coefficients in the plot. A larger scaling value will reduce the size of the coefficients in the plot.
#' @export
plot.waveletCoef <- function(x, ..., lowest = NULL, highest = NULL, coefTrim = NULL, thickness = 1, descending = FALSE, scaling = 1){
  j0 <- x$j0
  x <- x$coef
  J <- log2(length(x))
  if((J %% 1) > 0)
    warning("Vector of wavelet coefficients has length not a power of 2\n")
  wc <- x[-(1:(2^j0))]
  fine <- ceiling(J) - 1
  # Check resolution ranges
  if (is.null(lowest))
    lowest <- j0
  if (is.null(highest)) 
    highest <- fine
  
  reslev <- rep(j0:fine, 2^(j0:fine))
  num <- numeric()
  for (i in (min(reslev):max(reslev))) {
    num <- c(num, 2*(1:(2^i)) - 1)
  }
  M <- 2 * max(abs(wc))/scaling
  den <- 2^(reslev + 1)
  ind <- (reslev >= lowest) & (reslev <= highest)
  x <- (num/den)[ind]
  y1 <- reslev[ind]
  if (descending) 
    y1 <- (fine + lowest - y1)
  leny1 = length(y1)
  lenwcind <- length(wc[ind])
  y2 <- y1 + wc[ind]/M
  cbind(x, y1, y2)
  plot(0, type = "n", xlim = range(x), ylim = range(c(y1, y2)), yaxt = 'n', xlab = "Location",
       ylab = "Resolution Level", main = "MultiResolution Analysis of Coef.")
  axis(2, at = lowest:highest)
  for (i in 1:length(x)) {
    lines(c(x[i], x[i]), c(y1[i], y2[i]))
  }
  if( !is.null(coefTrim) ){
    wc <-  coefTrim[-(1:2^j0)]
    M <- 2 * max(abs(wc))/scaling
    den <- 2^(reslev + 1)
    ind <- (reslev >= lowest) & (reslev <= highest)
    x <- (num/den)[ind]
    y1 <- reslev[ind]
    if (descending) 
      y1 <- (fine + lowest - y1)
    leny1 = length(y1)
    lenwcind <- length(wc[ind])
    y2 <- y1 + wc[ind]/M
    for (i in 1:length(x)) {
      if( y1[i] != y2[i] ){
        lines(c(x[i], x[i]), c(y1[i], y2[i]), col = 2, lwd = thickness)
      }
    }
  }
}

#' @title Plot Output for the mWaveD object
#' 
#' @description  that summarises the multichannel input, a visualisation of the severity of the channels and the output estimate.Four plots: The input multichannel signal, the mWaveD
#' 
#' @param x A mWaveD object to be plotted (list created by \code{\link{multiWaveD}})
#' @param ... Arguments to be passed to methods.
#' @param singlePlot A logical value that controls whether all four plots appear on a single window (2 x 2 plot window) or are separated into separate plot windows.
#' @param prompt A logical value that specifies whether the user is prompted between plot outputs.
#' 
#' @details Four plots are output that summarise the multichannel input, a visualisation of the severity of the channels and the output estimate.\itemize{
#' \item Plot 1: Multichannel input signal overlayed.
#' \item Plot 2: Estimated output signal using the mWaveD approach.
#' \item Plot 3: 
#' \item Plot 4: Multi-resolution plot of the raw wavelet coefficients (black) and the trimmed wavelet coefficients (red)}
#' @references TBA ACHA and COMPSTAT
#' @seealso \code{\link{multiWaveD}}
#' @export
plot.mWaveD <- function(x, ..., singlePlot = TRUE, prompt = TRUE){
  n = length(x$estimate)
  t = (1:n)/n
  if( singlePlot ){
    par(mfrow=c(2,2))
  }

  matplot(t,x$signal,type='l',main="Input Signal",ylab='Signal',xlab='t',lty=1,cex=0.8)
  if( !singlePlot && prompt ){
    readline('Press any key to see the next plot.')
  }
  plot(t,x$estimate,type='l',main='mWaveD estimate',ylab='mWaveD Estimate',xlab='t')
  blurInfo = x$blurInfo
  blur = x$blurType
  j0 <- blurInfo$j0
  j1 <- blurInfo$j1
  if( blur == "direct" || blur == "smooth" ){
    blurf = blurInfo$decay
    cutf = blurInfo$cutoff
    revblur = blurf[(n/2):2,]
    revcut = cutf[(n/2):2,]
    tS = (1:(n/2))
    tS = c(-rev(tS[-(n/2)]),0,tS)
    xlim = min(2*max(blurInfo$freqCutoff),n/2)
    xlim = c(-xlim,xlim)
    ylim = c(min(cutf[2,]),0)
    blur = rbind(revblur,blurf)
    cut = rbind(revcut,cutf)
    if( !singlePlot && prompt ){
      readline('Press any key to see the next plot.')
    }
    matplot(tS,blur,type='l',lty=1,xlim=xlim,ylim=ylim,main="Fourier decay and cutoffs",xlab="Frequency", ylab="Fourier decay")
    matlines(tS,cut,lty=2)
  } else {
    J    <- floor(log2(n)) - 1
    j    <- j0:min(c(J - 1, 2 * j1))
    blkV <- blurInfo$blockVar[1:length(j)]
    blkc <- blurInfo$blockCutoff[1:length(j)]
    ylim <- range(c(blurInfo$blockVar, blurInfo$blockCutoff))
    if( !singlePlot && prompt ){
      readline('Press any key to see the next plot.')
    }
    plot(j,blkV,type='b',xlab='j',ylab='',main='Block wise resolution selection')
    lines(j,blkc,col=2)
    points(j1,blurInfo$blockVar[j == j1],col='blue')
    lines(c(j1,j1),c(ylim[1],blurInfo$blockVar[j == j1]),lty='dashed')
  }
  
  beta <- list(coef = x$coef, j0 = j0)
  class(beta) <- 'waveletCoef'
  if( !singlePlot && prompt ){
    readline('Press any key to see the next plot.')
  }
  plot(beta, highest = j1, coefTrim = x$shrinkCoef)
}

