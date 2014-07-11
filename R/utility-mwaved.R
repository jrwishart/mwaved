#' @title Summary Output for the mWaveD object
#'
#' @description Summarises the mWaveD object
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
#' @param lowest Specifies the coarsest resolution to display in the Multi-resolution plot
#' @param highest Specifies the finest resolution to display in the Multi-resolution plot 
#' 
#' @export
plot.waveletCoef <- function(x, ..., lowest = NULL, highest = NULL, thickness = 1, coefTrim = NULL, descending = FALSE, sc = 1){
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
  M <- 2 * max(abs(wc))/sc
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
    M <- 2 * max(abs(wc))/sc
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
#' @description Four plots: The input multichannel signal, the mWaveD
#' 
#' 
#' @export
plot.mWaveD <- function(x,...){
  n = length(x$estimate)
  t = (1:n)/n
  
  par(mfrow=c(2,2))
  
  matplot(t,x$signal,type='l',main="Input Signal",ylab='Signal',xlab='t',lty=1,cex=0.8)
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
    matplot(tS,blur,type='l',lty=1,xlim=xlim,ylim=ylim,main="Fourier decay and cutoffs",xlab="Frequency", ylab="Fourier decay")
    matlines(tS,cut,lty=2)
  } else {
    J    <- floor(log2(n)) - 1
    j    <- j0:min(c(J - 1, 2 * j1))
    blkV <- blurInfo$blockVar[1:length(j)]
    blkc <- blurInfo$blockCutoff[1:length(j)]
    ylim <- range(c(blurInfo$blockVar, blurInfo$blockCutoff))
    plot(j,blkV,type='b',xlab='j',ylab='',main='Block wise resolution selection')
    lines(j,blkc,col=2)
    points(j1,blurInfo$blockVar[j == j1],col='blue')
    lines(c(j1,j1),c(ylim[1],blurInfo$blockVar[j == j1]),lty='dashed')
  }
  
  beta <- list(coef = x$coef, j0 = j0)
  class(beta) <- 'waveletCoef'
  plot(beta, highest = j1, coefTrim = x$shrinkCoef)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   require(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#     
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }

# ggplot.mWaveD <- function(mwaved.obj,...){
#   #   reduce.margins <- theme(axis.title = element_blank(),plot.margin=unit(c(0.25,0.25,0,0),"lines"));
#   lsize=0.25;
#   hsize=1;
#   n <- dim(mwaved.obj$signal)[1]
#   m <- mwaved.obj$channels
#   x <- 1:n/n
#   signal.df <- data.frame(Y=as.vector(mwaved.obj$signal),x=rep(x,m),Channel=rep(LETTERS[1:m],each=n))
#   est.df <- data.frame(Y=mwaved.obj$estimate,x=x)
#   sig.plot <- ggplot(signal.df) + geom_line(aes(x=x, y=Y, colour=Channel),size=lsize) + scale_color_discrete(labels= LETTERS[1:m]) + scale_size(guide='none') + guides(colour = guide_legend(title="Signal Input",title.position='left', title.theme = element_text(size=15,angle=0),override.aes = list(size = 10))) + theme(legend.position="top",legend.key=element_rect(fill=NA))
#   est.plot <- ggplot(est.df, aes(x=x, y=Y),colour='blue') + geom_line(size=lsize,alpha=1) + ggtitle(paste(mwaved.obj$shrinkType, ' Signal Estimate',sep='')) # + reduce.margins 
#   # Check blur type
#   blurInfo = mwaved.obj$blurInfo
#   blur = mwaved.obj$blurType
#   j0 <- blurInfo$j0
#   j1 <- blurInfo$j1
#   
#   if( blur == "direct" || blur == "smooth" ){
#     blurf  <- blurInfo$decay
#     cutf   <- blurInfo$cutoffs
#     rblurf <- blurf[(n/2):2,]
#     rcutf  <- cutf[(n/2):2,]
#     tS     <- (1:(n/2))
#     tS     <- c(-rev(tS[-(n/2)]),0,tS)
#     xlim   <- min(1.5*max(blurInfo$freqCutoff),n/2)
#     xlim   <- c(-xlim,xlim)
#     ylim   <- c(min(cutf[2,]),0.1)
#     blurf  <- rbind(rblurf,blurf)
#     cut    <- rbind(rcutf,cutf)
#     blur.df <- data.frame(Y=as.vector(blurf),x=rep(tS,m),Ycut = as.vector(cut),Channel=rep(LETTERS[1:m],each=n))
#     blur.plot <- ggplot(blur.df) + geom_line(aes(x=x, y=Y, colour=Channel,group=Channel),size=hsize) + geom_line(aes(x=x, y=Ycut, colour=Channel),linetype=2, size=hsize)  + coord_cartesian(xlim = xlim,ylim=ylim)  + scale_size(guide='none') + theme(legend.position=c(1,1),legend.key=element_rect(fill=NA)) + ggtitle("Fourier Decay")
#   } else {
#     # Blur is box car
#     J    <- floor(log2(n)) - 1
#     x    <- j0:min(c(J - 1, 2 * j1))
#     nx   <- length(x)
#     blkV <- blurInfo$blockVar[1:length(x)]
#     blkc <- blurInfo$blockCutoff[1:length(x)]
#     # Mark the selection
#     j.df <- data.frame(j1 = rep(j1,2), val = c(min(blkV), blkV[x == j1]))
# #     ylim <- 1.05*range(c(blurInfo$blockVar[1:nx], blurInfo$blockCutoff[1:nx]))
#     blur.df <- data.frame(y = c(blkV,blkc), x=rep(x,2), lty = rep(c('A','B'),each=nx), colour = rep(c('red','blue'),each=nx))
#     blur.plot <- ggplot(blur.df) + geom_line(aes(x=x, y=y, colour=colour, linetype=lty)) + geom_point(aes(x=x, y=y, linetype=lty, colour=colour)) + ggtitle('Block-wise resolution selection') + labs(x="j") + theme(legend.position="none") + geom_line(aes(x=j1, y=val), linetype='dashed', colour = 'red',data=j.df) + geom_point(aes(x=j1, y=val), colour = 'red',data=j.df)
#   }
#   
#   js <- rep(j0:j1,2^(j0:j1))
#   ks <- unlist(lapply(j0:j1,function(i) (0:(2^i-1)/2^i)))
#   wi <- (2^j0 + 1):2^(j1 + 1)
#   nw <- length(wi)
#   w  <- mwaved.obj$coef[wi]
#   wf <- 2.1*max(w)
#   w  <- w/wf
#   wc <- mwaved.obj$shrinkCoef[wi]/wf
#   
#   # size of the coefficents in MRA plot
#   raw.size = 1
#   shrink.size = 1
#   mra.df <- data.frame(w = c(w + js,wc+js),js=rep(js,2),ks=rep(unlist(lapply(j0:j1,function(i) (0:(2^i-1)/2^i))),2),MRA=rep(LETTERS[1:2],each=nw),wsize = rep(c(raw.size,shrink.size),each=nw),method=rep(c('Raw',paste(mwaved.obj$shrinkType,' Thresholding',sep='')),each=nw))
#   
#   mra.plot <- ggplot(mra.df) + geom_segment(aes(x=ks, y=js,xend=ks,yend=w,colour=MRA,size=wsize)) + labs(x="k",y='j') +scale_color_discrete(labels= c('Raw',paste(mwaved.obj$shrinkType,' Thresholding',sep='')),guide=guide_legend(title.position='left',title.theme = element_text(size=15,angle=0))) + scale_size(guide='none') + guides(colour = guide_legend(override.aes = list(size = 10), title='Multiresolution Analysis')) + theme(legend.position="top",legend.key=element_rect(fill=NA), axis.text.y = element_text(angle=90))
#   # , Raw = Blue, Red = Thresholded
#   multiplot(sig.plot,blur.plot,est.plot,mra.plot,cols=2)
# }