#' @name summary.mWaveD
#' @title Summary Output for the mWaveD object
#'
#' @param object A mWaveD object which is a list containing all the information for a multichannel 
#' deconvolution analysis produced by the \code{\link{multiWaveD}} function.
#' @param ... Arguments to be passed to methods.
#' @description Gives some numerical summaries of a \code{mWaveD} object.
#' @return Text output giving summary information of the input and output analysis including, \itemize{
#' \item Degree of Meyer wavelet used in the analysis.
#' \item Number of observations, within each channel and number of channels present.
#' \item Resolution levels used (j0 to j1)
#' \item Blur type assumed in the analysis (direct, smooth or box.car)
#' \item Matrix summarising the noise levels in each channel (and Fourier decay information for the smooth case)
#' \item Summaries of the severity of the thresholding applied amongst the resolutions.
#' }
#' @seealso \code{\link{multiWaveD}}
#' 
#' @examples
#' library(mwaved)
#' # Simulate the multichannel doppler signal.
#' m <- 3
#' n <- 2^10
#' t <- (1:n)/n
#' signal <- makeDoppler(n)
#' # Create multichannel version with smooth blur
#' shape <- seq(from = 0.5, to = 1, length = m)
#' scale <- rep(0.25, m)
#' G <- gammaBlur(n, shape, scale)
#' X <- blurSignal(signal, G)
#' # Add noise with custom signal to noise ratio
#' SNR <- c(10,15,20)
#' E <- multiNoise(n, sigma = sigmaSNR(X, SNR), alpha = c(0.5, 0.75, 1))
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' mWaveDObject <- multiWaveD(Y, G)
#' summary(mWaveDObject)
#' @export
summary.mWaveD <- function(object, ...){
  n <- length(object$estimate)
  
  cat("Degree of Meyer wavelet =", object$degree, "  , Coarse resolution level j0 =", object$j0)
  cat("\n")
  cat("Sample size per channel = ", n, ", Maximum possible resolution level = ", log2(n) - 1, ".", sep = '')
  cat("\n\n")
  cat("Number of channels: m =", object$channels,"\n")
  cat('Blur type: ',object$blurType,'\n\n')
  cat("Estimated Channel information:\n\n")
  
  if (object$blurType == "direct") {
    mat <- cbind(round(object$sigma, 3), round(object$alpha, 3), object$blurInfo$freq, rep(object$j1, object$channels))
    colnames(mat) <- c("Sigma.hat", "Alpha", "Fourier freq cutoff", "Highest resolution")
    rownames(mat) <- paste("Channel ", 1:object$channels,':', sep='')
    print(mat, ...)
  } else {
    # If Smooth blur is used, display the matrix of values
    if (object$blurType == "smooth"){
      mat <- cbind(round(object$sigma, 3), round(object$alpha, 3), object$blurInfo$freq, object$blurInfo$maxLevels)
      colnames(mat) <- c("Sigma.hat", "Alpha", "Fourier freq cutoff", "Highest resolution")
      rownames(mat) <- paste("Channel ", 1:object$channels,':', sep='')
      print(mat, ...)
      cat("\n")
      cat("Estimated best channel = Channel", object$blurInfo$bestChannel)
      
    } else {
      if (object$blurType == "box.car"){
        mat <- cbind(round(object$sigma, 3), round(object$alpha, 3))
        colnames(mat) <- c("Sigma.hat", "Alpha")
        rownames(mat) <- paste("Channel ", 1:object$channels,':', sep='')
        print(mat, ...)
      } else {
        warning('Unrecognised blur.type')
      }
    } 
  }
  cat("\n\n")
  cat("mWaveD optimal finest resolution level j1 =", object$j1)

  cat("\n\n")
  cat("Thresholding method:", object$shrinkType, "   Tuning parameter: eta =", object$eta,'\n\n')
  
  threshMatrix <- cbind(round(object$levelMax,4), round(object$thresh,4), object$percent)
  rownames(threshMatrix) <- paste("Level", object$j0:object$j1,":" )
  colnames(threshMatrix) <- c("Max|w|", "Threshold", "% Shrinkage" )
  print(threshMatrix, ...)
} 

#' @name plot.waveletCoef
#' @title Multi-Resolution Analysis plot of wavelet coefficients
#'
#' @description Plots the wavelet coefficient object in the multiresolution analysis
#' 
#' @param x A list of class waveletCoef.
#' @param y An optional numeric vector of trimmed wavelet coefficients to be overlayed on top of the plot for comparison with the \code{x} wavelet coefficients. 
#' @param labels Optional character vector with two elements to give name labels to \code{x} and \code{y} respectively.
#' @param ... Arguments to be passed to methods.
#' @param lowest Specifies the coarsest resolution to display in the Multi-resolution plot.
#' @param highest Specifies the finest resolution to display in the Multi-resolution plot.
#' @param scaling A numeric value that acts as a graphical scaling parameter to rescale the wavelet coefficients in the plot. A larger scaling value will reduce the size of the coefficients in the plot.
#' @param ggplot A logical value to specify if the user wants to use base graphics (FALSE) or ggplot2  graphics (TRUE).
#' 
#' @seealso \code{\link{multiCoef}} for generating a list of class `waveletCoef`
#' 
#' @export
plot.waveletCoef <- function(x, y = NULL, labels = NULL,  ..., lowest = NULL, highest = NULL, scaling = 1, ggplot = TRUE){
  if (!is.null(y) && class(y) != "waveletCoef") {
    stop('y argument must be a waveletCoef object')
  }
  
  J <- floor(log2(length(x$coef)))
  fine <- ceiling(J) - 1
  # Check resolution ranges
  if (is.null(lowest)) {
    lowest <- x$j0
  } else {
    # Catch lowest level too low
    if (lowest > x$j0 )
      warning("lowest level shouldn't be smaller than j0 specified in wavelet coefficient object.")
  }
  if (is.null(highest)) {
    highest <- fine
  } else {
    if (highest > fine) {
      warning(paste('highest level too high. Resetting highest level to the maximum at j1 = ', fine))
    } else {
      if (highest < lowest) {
        warning('highest level must be higher than the lowest level.')
        highest <- lowest
      }
    }
  }
  
  js <- rep(lowest:highest, 2^(lowest:highest))
  ks <- unlist(lapply(lowest:highest, function(i) 0:(2^i-1)/2^i))
  wi <- (2^lowest + 1):2^(highest + 1)
  nw <- length(wi)
  w  <- x$coef[wi]*scaling
  wf <- 2.05 * max(abs(w))/scaling
  w  <- w/wf
  ws <- w + js
  
  # check shrink input is a waveletCoef object
  if (!is.null(y)) {
    if (length(x$coef) != length(y$coef)) {
      stop('length of y coefficients is different to the length of x coefficients')
    }
    j0Trim <- y$j0
    if (j0Trim != x$j0) {
      warning('y object has a different coarse resolution j0 than x object, interpret lower resolution levels with caution')
    }
    wShrink <- y$coef[wi]/wf
    survived <- which(wShrink != 0)
    kss <- ks[survived]
    jss <- js[survived]
    wss <- wShrink[survived] + jss
    ns <- 2
    nss <- length(survived)
  } else {
    ns <- 1
  }
  
  mraTitle <- 'Multiresolution Analysis of Coef.'
  mraLabels <- c("Location", "Resolution Level")
  if (!is.null(labels)) {
    if (length(labels) != ns) {
      warning('length of labels might not be long enough.')
      labels = c('x', 'y')
    }
  } else {
    labels = c('x', 'y')
  }
  
  if (ggplot && require(ggplot2)) {
    if (ns == 2) {
      nData <- data.frame(w = c(ws,wss), js = c(js, jss), ks = c(ks, kss), col = rep(labels, c(nw, nss)))
      mraPlot <- ggplot(nData) + geom_segment(aes(x = ks, xend = ks, y = js, yend = w, colour = col, size = col)) + labs(x = mraLabels[1], y = mraLabels[2]) + scale_size_discrete(range = c(1, 2))
      # Fix legend
      mraPlot <- mraPlot + theme(legend.position = "top", axis.text.y = element_text(angle = 90)) + guides(colour = guide_legend(title = mraTitle), size = guide_legend(title = mraTitle))
    } else {
      nData <- data.frame(w = ws, js = js, ks = ks)
      mraPlot <- ggplot(nData) + geom_segment(aes(x = ks, xend = ks, y = js, yend = w), colour = 'red') + ggtitle(mraTitle)
    }
    mraPlot + labs(x = mraLabels[1], y = mraLabels[2]) + scale_y_continuous(breaks = lowest:highest) 
  } else {
    buf <- 0.5
    plot(0, type = "n", xlim = c(0,1), ylim = c(lowest - buf, highest + buf), yaxt = 'n', xlab = mraLabels[1], ylab = mraLabels[2], main = mraTitle)
    axis(2, at = lowest:highest)
    if (!is.null(y)) {
      col = 'red'
    } else {
      col = 1
    }
    
    segments(ks, js, ks, ws, col = col)
    
    if (!is.null(y)) {
      segments(kss, jss, kss, wss, lwd = 2, col = 'blue')
    }
  }
}

#' @name plot.mWaveD
#' @title Plot Output for the mWaveD object
#' 
#' @description Creates plot output that summarises the \code{mWaveD} object produced by the \code{\link{multiWaveD}} function. 
#'  
#' @param x A mWaveD object to be plotted (list created by \code{\link{multiWaveD}})
#' @param ... Arguments to be passed to methods.
#' @param which A numeric vector that specifies which plots to output. Default value is \code{1:4} which  specifies that all four plots are to be displayed.
#' @param ask A logical value that specifies whether the user is \emph{ask}ed before each plot output.
#' @param singlePlot A logical value that controls whether all plots should appear on a single window. The plot window is resized depending the value of \code{which}.
#' @param ggplot A logical value to specify if the user wants to use base graphics (FALSE) or ggplot2 graphics (TRUE).
#' 
#' @details Four plots are output that summarise the multichannel input, a visualisation of the characteristics of the channels and the output estimate and a multi-resolution analysis plot.\itemize{
#' \item Plot 1: Multichannel input signal overlayed.
#' \item Plot 2: Estimated output signal using the mWaveD approach.
#' \item Plot 3: Plot of the log decay of Fourier coefficients against the log bounds (direct and smooth case) or the blockwise resolution levels against their limit (box car case)
#' \item Plot 4: Multi-resolution plot of the raw wavelet coefficients and the trimmed wavelet coefficients}
#' @references
#' Kulik, R., Sapatinas, T. and Wishart, J.R. (2014) \emph{Multichannel wavelet deconvolution with long range dependence. Upper bounds on the L_p risk}  Appl. Comput. Harmon. Anal. (to appear in).
#' \url{http://dx.doi.org/10.1016/j.acha.2014.04.004}
#' 
#' Wishart, J.R. (2014) \emph{Data-driven wavelet resolution choice in multichannel box-car deconvolution with long memory}, Proceedings of COMPSTAT 2014, Geneva Switzerland, Physica Verlag, Heidelberg (to appear)
#' @seealso \code{\link{multiWaveD}}
#' @export
plot.mWaveD <- function(x, ..., which = 1L:4L, singlePlot = TRUE, ask = !singlePlot, ggplot = TRUE){
  showPrompt <- function(condition) {
    if (condition) {
      readline('Press any key to see the next plot.')  
    }
  }
  
  # Check if ggplot is available if requested
  if (ggplot) {
    ggAvailable <- require(ggplot2)
    # Check optional dependency
    gridExtraAvailable <- require(gridExtra)
    if (ggAvailable) {
      hsize <- 1
      lsize <- 0.5
      asize <- 0.5
      library(ggplot2)
      # Initialise list
      ggList <- list(NULL)
      i <- 1
    } 
    if (singlePlot) {
      if (gridExtraAvailable) {
        library(gridExtra)
      } else {
        warning('gridExtra package required to create ggplot2 graphics in same window. Setting output to separate windows.')
        singlePlot <- FALSE
      }
    }
  } else {
    ggAvailable <- FALSE
  }
  
  show <- rep(FALSE, 4)
  # Make sure which argument is numeric
  if (!is.numeric(which)) {
    stop('`which` argument must be a vector containing numerical elements')
  }
  # Only consider the integer bits between 1:4
  which <- which[which %in% 1L:4L]
  # Complain if there are no values in 1:4
  if (length(which) == 0) {
    stop('`which` argument must be a vector containing elements: 1, 2, 3 or 4')
  }
  show[which] <- TRUE
  
  n  <- length(x$estimate)
  n2 <- n/2
  m  <- x$channels
  t  <- (1:n)/n
  
  blurInfo <- x$blurInfo
  blurType <- x$blurType
  j0 <- x$j0
  j1 <- x$j1
  
  estimateTitle <- 'mWaveD estimate'
  signalTitle <- 'Input Signal'
  fourierLabel <- 'Fourier freq'
  fourierTitle <- 'Kernel decay in Fourier domain'
  blockTitle <- 'Block wise resolution selection'
  mraLabels <- c('raw',paste(x$shrinkType, ' thresholded', sep = ''))
  mraTitle <- 'Multiresolution Analysis'
  
  if (show[3L]) {
    if (blurType != "box.car") {
      blurf <- blurInfo$decay
      cutf <- blurInfo$cutoff
      ymin <- min(cutf)
      revblur <- as.matrix(blurf[n2:2, ])
      revcut <- as.matrix(cutf[n2:2, ])
      iw = -(n2 - 1):n2
      blur <- rbind(revblur, blurf)
      cut <- rbind(revcut, cutf)
      ylim <- c(min(cutf[2, ]), 0)
      if (blurType == 'smooth') {
        xbest <- max(blurInfo$freqCutoffs) - 1
        ybest <- blurf[xbest, blurInfo$bestChannel]
        xlim <- min(2*max(blurInfo$freqCutoff), n/2)
        xlim <- c(-xlim, xlim)
      } else {
        xlim = c(-n2, n2)
      }
    } else {
      J    <- floor(log2(n)) - 1
      j    <- j0:min(c(J - 1, 2 * j1))
      blkV <- blurInfo$blockVar[1:length(j)]
      blkc <- blurInfo$blockCutoff[1:length(j)]
      ylim <- range(c(blurInfo$blockVar, blurInfo$blockCutoff))
    }
  }
  if (show[1L] && ggAvailable) {
    signalData <- data.frame(Y = as.vector(x$signal), x = rep(t, m), Channel = rep(LETTERS[1:m], each = n))
    signalPlot <- ggplot(signalData, aes_string(x = 'x', y = 'Y', colour = 'Channel')) + geom_line(size = lsize, alpha = asize) + ggtitle(signalTitle) + labs(x = '', y = '')
    ggList[[i]] <- signalPlot
    i <- i + 1
  }
  if (show[2L] && ggAvailable) {
    estimateData <- data.frame(Y = as.vector(x$estimate), x = t)
    estimatePlot <- ggplot(estimateData, aes_string(x = 'x', y = 'Y')) + geom_line(size = lsize, alpha = asize) + ggtitle(estimateTitle) + labs(x = '', y = '')
    ggList[[i]] <- estimatePlot
    i <- i + 1
  }
  if (show[3L] && ggAvailable) {
    if (blurType != 'box.car') {
      fourierData <- data.frame(Y = as.vector(blur), x = rep(iw,m), Ycut = as.vector(cut), Channel=rep(LETTERS[1:m],each=n),m=m)
      resolutionPlot <- ggplot(fourierData) + geom_line(aes_string(x = 'x', y = 'Y', colour = 'Channel', group = 'Channel'),size = 1) + geom_line(aes_string(x = 'x', y = 'Ycut', colour = 'Channel'), linetype='dashed', size = 1) + ggtitle(fourierTitle) + labs(x = fourierLabel, y = '')
      if (blurType == 'smooth') {
        highlightData <- data.frame(x = rep(xbest, 2), y = c(-Inf, ybest))
        pointData <- data.frame(xbest = xbest, ybest = ybest)
        resolutionPlot <- resolutionPlot + geom_line(aes_string(x = 'x', y = 'y'), linetype = 'dotted' , data = highlightData) + geom_line(aes_string(x = 'x', y = 'y'), linetype = 'dotted' , data = data.frame(x = rep(-xbest + 1, 2), y = c(-Inf, ybest))) + geom_point( aes_string(x = 'xbest', y = 'ybest'), size = 4, shape = 1, data = pointData) + geom_point( aes_string(x = 'x', y = 'y'), size = 4, shape = 1, data = data.frame(x = -xbest + 1, y = ybest)) + coord_cartesian(xlim = xlim)
      }
    } else {
      resolutionData <- data.frame(Y = c(blkV, blkc), x = rep(j,2), colour = rep(c("Resolution var.",'Resolution bounds'), each = length(j)) , Ycut = blkc)
      bestV <- blkV[j == j1]
      highlightData <- data.frame(x = c(j1, j1), y = c(ylim[1], bestV))
      pointData <- data.frame(j1 = j1, bestV = bestV)
      resolutionPlot <- ggplot(resolutionData) + geom_line(aes_string(x = 'x', y = 'Y', colour = 'colour', linetype = 'colour'), size = hsize) +  geom_line(aes_string(x = 'x', y = 'y'), linetype = 'dotted', data = highlightData) + labs(x = 'j', y = '') + geom_point( aes_string(x = 'j1', y = 'bestV'), size = 4, shape = 1, data = pointData)  + scale_color_discrete(labels= c('Resolution bounds', 'Resolution var.'), guide=guide_legend(title.position='left',title.theme = element_text(size=15,angle=0))) + scale_size(guide='none') + guides(colour = guide_legend( title='Blockwise resolution decay')) + theme(legend.position="top", legend.key = element_rect(fill = NA), axis.text.y = element_text(angle = 90)) + scale_linetype_manual(values=c(1,2), name="Blockwise resolution decay", labels=c('Resolution bounds', 'Resolution var.')) + scale_x_continuous(breaks = j)
    }
    ggList[[i]] <- resolutionPlot
    i <- i + 1
  }
  
  if (show[4L] && ggAvailable) {
      mraPlot <- plot(x$coef, x$shrinkCoef, highest = j1, labels = c('Raw', paste(x$shrinkType, ' Thresholding', sep = '')), ggplot = TRUE)
      ggList[[i]] <- mraPlot
  }

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  if (ggAvailable) {
    # Plot them
    if (singlePlot == TRUE) {
      do.call(grid.arrange, ggList)  
    } else {
      if (show[1L]) {
        print(signalPlot)
      }
      if (show[2L]) {
#         showPrompt(ask)
        print(estimatePlot)
      }
      if (show[3L]) {
#         showPrompt(ask)
        print(resolutionPlot)
      }
      if (show[4L]) {
#         showPrompt(ask)
        print(mraPlot)
      }
    }
  } else {
    if (singlePlot) {
      plotDims <- switch(sum(show), c(1, 1), c(2, 1), c(3, 1), c(2, 2))
      par(mfrow = plotDims)
    } else {
      par(mfrow = c(1,1))
    }

    if (show[1L]) {
      matplot(t, x$signal, type = 'l', main = signalTitle, ylab = '', xlab = '', lty = 1, cex = 0.8)
      grid()
    }

    if (show[2L]) {
      # Plot mWaveD estimate
      plot(t, x$estimate, type = 'l', main = estimateTitle, ylab = '', xlab = '')
      grid()
    }

    if (show[3L]) {
      # Plot resolution analysis
      if (blurType != 'box.car') {
        matplot(iw, blur, type = 'l', lty = 1, xlim = xlim, ylim = ylim, main = fourierTitle, xlab = fourierLabel, ylab = "")
        matlines(iw, cut, lty = 2)
        grid()      
        if (blurType == 'smooth') {
          points(xbest, ybest, col='blue')
          points(-xbest + 1, ybest, col = 'blue')
          xbest <- rep(xbest, 2)
          ybest <- c(ylim[1], ybest)
          lines(xbest, ybest, lty = 'dotted')
          lines(-xbest, ybest, lty = 'dotted')
        }
      } else {
  #       showPrompt(ask)
        plot(j, blkV, type = 'b', xlab = 'j', ylab = '', main = blockTitle)
        lines(j, blkc, col = 2)
        points(j1, blurInfo$blockVar[j == j1], col='blue')
        lines(c(j1, j1), c(ylim[1], blurInfo$blockVar[j == j1]), lty = 'dashed')
        grid()
      }
    }
#     showPrompt(ask)
    if (show[4L]) {
      plot(x$coef, x$shrinkCoef, highest = j1, ..., ggplot = FALSE)
    }
  }
}

#' @name mWaveDDemo
#' @title Interactive Demonstration
#' @description Interactive Demonstration
#' 
#' @export
mWaveDDemo <- function (){
  runApp(system.file('mWaveDDemo', package = 'mwaved'))
}