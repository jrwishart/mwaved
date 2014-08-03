library(shiny)
library(mwaved)
# Check if ggplot2 is available and use it
ggplot2Avail <- 'ggplot2' %in% .packages(all.available = TRUE)
# Graphical settings
shrinkLabel <- 'Shrinkage Type:'
shrinkChoices <- list("Hard" = "Hard", "Soft" = "Soft", "Garrote" = "Garrote")
fill.col <- 'light blue'
outline.col <- 'blue'
density.col <- 'dark grey'
highlight.col <- 'red'
if (ggplot2Avail){
  library(ggplot2)
  blanklabs <- labs(x = '', y = '')  
}
primenum <- c(83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199)
# reduce.margins <- theme(axis.title = element_blank(),plot.margin=unit(c(0.25,0.25,0,0),"lines"));
lsize=0.5;
hsize=1;
xcat <- function(obj){
  paste('c(',paste(as.numeric(obj), collapse = ','),')',sep='')
}
numclean <- function(obj){
  return(round(as.numeric(obj), 3))
}


shinyServer(function(input, output, session) {

  ###############################################
  ## Compute Signals and other bits
  ###############################################	
  
  currentShrinkage <- 'Hard'
  
  sigList <- reactive({
	 n     <- 2^input$J
	 m     <- input$m;
   if (m == 1){
     SNR <- sample(input$SNR[1]:input$SNR[2], 1)
   } else {
     SNR <- sample(seq(from = input$SNR[1], to = input$SNR[2], length = m), m);  
   }
	 
   alpha <- sample(seq(from = input$alpha[1], to = input$alpha[2], length = m), m);
	 shape <- sample(seq(from = input$gshape[1], to = input$gshape[2], length = m), m);
	 scale <- sample(seq(from = input$gscale[1], to = input$gscale[2], length = m), m);
	 BA    <- 1/sqrt(sample(primenum, m))
	 G     <- switch(input$blur, 
              'smooth'= gammaBlur(n, shape, scale),
              'direct' = directBlur(c(n, m)), 
              'box.car'= boxcarBlur(n, BA))
	 
	 signal   <- switch(input$sig,
              'lidar' = makeLIDAR(n), 
				      'doppler' = makeDoppler(n),
              'bumps' = makeBumps(n),
              'blocks' = makeBlocks(n))
	 signalName <- switch(input$sig,'lidar' = "LIDAR Signal",
				'doppler' = "Doppler Signal",
				'bumps' = "Bumps Signal",
				'blocks' = "Blocks Signal")
	 signalBlur <- blurSignal(signal, G)
   sigma <- sigmaSNR(signalBlur, SNR)
   eps <- multiNoise(n, sigma, alpha)
	 Y   <- signalBlur + eps;
	 x   <- (1:n)/n;
	 
   list(m = m, n = n, signal = Y, G = G, blur = input$blur, SNR = SNR, alpha = alpha, shape = shape, 
        scale = scale , G = G, trueSignal = signal, eps = eps, signalBlur = signalBlur, BA = BA, 
        x = x, signalName = signalName, sigma = sigma)
  })
  
  mWaveDList <- reactive({
    
    multiSig <- sigList()
    mlwvd <- multiWaveD(multiSig$signal, multiSig$G, alpha = multiSig$alpha, 
                        blur = multiSig$blur, shrinkType = input$shrinkage1, deg = as.integer(input$degree))
    
    j0 <- mlwvd$j0
    j1 <- mlwvd$j1
    
    list(n = multiSig$n, m = multiSig$m, mWaveD = mlwvd, j0 = mlwvd$j0, j1 = mlwvd$j1, coef = mlwvd$coef, shrinkCoef = mlwvd$shrinkCoef)
      
  })
  
  ###############################################
  ## Render WvdPlots
  ###############################################	
  
  output$reactiveSignal <- renderPlot({
    sList    <- sigList()
    
    # Base signal selected
    if (input$signalShow == 1){
      plotTitle <- sList$signalName
      if (ggplot2Avail){
        plot.df <- data.frame(x = sList$x, signal = sList$trueSignal)
        print(ggplot(data = plot.df, aes(x = x, y = signal)) + geom_line(size = lsize) + ggtitle(plotTitle) + blanklabs)
      } else {
        plot(sList$x, sList$trueSignal, type = 'l', main = plotTitle, ylab = '', xlab = '')
        grid()
      }
       # + reduce.margins
    }
    if (input$signalShow == 2){
      plotTitle <- paste('Blurred',sList$signalName)
      if (ggplot2Avail){
        plot.df <- data.frame(Y = as.vector(sList$signalBlur), x = rep(sList$x, sList$m), 
                              Channel = rep(LETTERS[1:sList$m], each = sList$n))
        print(ggplot(plot.df, aes(x = x, y = Y, colour=Channel)) + geom_line(size = lsize, alpha = 0.9) + ggtitle(paste('Blurred ', sList$signalName, sep = "")) + blanklabs ) 
      } else {
        par(oma = c(4, 1, 1, 1))
        plot.out <- matplot(sList$x, sList$signalBlur, type = 'l', main = plotTitle, xlab = '', ylab = '')
        grid()
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        legend("bottom", paste('Channel ', LETTERS[1:sList$m]), xpd = TRUE, horiz = TRUE, 
               inset = c(0, 0), bty = "n", lty = 1:sList$m, col = 1:sList$m, cex = 1)
      }
       # + reduce.margins 
      # }
    }
    if (input$signalShow == 3){
      plotTitle <- paste('Noisy blurred', sList$signalName)
      if (ggplot2Avail){
        plot.df <- data.frame(Y = as.vector(sList$signal), x = rep(sList$x, sList$m), Channel = rep(LETTERS[1:sList$m], each = sList$n))
        print(ggplot(plot.df, aes(x = x, y = Y, colour = Channel)) + geom_line(size = lsize, alpha = 0.5) + ggtitle(plotTitle) + blanklabs)  
      } else {
        par(oma = c(4, 1, 1, 1))
        matplot(sList$x, sList$signal, type = 'l', main = plotTitle, xlab = '', ylab = '')
        grid()
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        legend("bottom", paste('Channel ', LETTERS[1:sList$m]), xpd = TRUE, horiz = TRUE, 
               inset = -0.1, bty = "n", lty = 1:sList$m, col = 1:sList$m, cex = 1)
      }
    }
  },  height=function() { session$clientData$output_reactiveSignal_width * 3/4 })
  
  output$multiPlot <- renderPlot({
    mList <- mWaveDList()
    
    j0 <- mList$j0
    j1 <- mList$j1
    
    js <- rep(j0:j1, 2^(j0:j1))
    ks <- unlist(lapply(j0:j1, function(i) (0:(2^i-1)/2^i)))
    wi <- (2^j0 + 1):2^(j1 + 1)
    nw <- length(wi)
    w  <- mList$coef$coef[wi]
    wf <- 2.1 * max(w)
    w  <- w/wf
    wc <- mList$shrinkCoef$coef[wi]/wf
    
    
    raw.size = 1
    shrink.size = 1
    if (ggplot2Avail){
      mra.df <- data.frame(w = c(w + js, wc+js), js = rep(js, 2), ks = rep(ks , 2), MRA = rep(LETTERS[1:2], each = nw), wsize = rep(c(raw.size, shrink.size), each = nw), method = rep(c('Raw', paste(mList$mWaveD$shrinkType, ' Thresholding', sep = '')), each = nw))
      
      mra.plot <- ggplot(mra.df) + geom_segment(aes(x=ks, y=js,xend=ks,yend=w,colour=MRA,size=wsize)) + labs(x="k",y='j') +scale_color_discrete(labels= c('Raw',paste(mList$mWaveD$shrinkType,' Thresholding',sep='')),guide=guide_legend(title.position='left',title.theme = element_text(size=15,angle=0))) + scale_size(guide='none') + guides(colour = guide_legend(override.aes = list(size = 10), title='Multiresolution Analysis')) + theme(legend.position="top",legend.key=element_rect(fill=NA), axis.text.y = element_text(angle=90))  + scale_y_continuous(breaks = j0:j1)
      print(mra.plot)
    } else {
      buf <- 0.5
      thickness <- 4
      plot(0, type = "n", xlim = c(0,1), ylim = c(j0 - buf, j1 + buf), yaxt = 'n', xlab = "",
           ylab = "Resolution Level", main = "MultiResolution Analysis of Coef.")
      axis(2, at = j0:j1)
      ws <- w + js
      survived <- which(wc != 0)
      kss <- ks[survived]
      jss <- js[survived]
      wss <- wc[survived] + jss
      segments(ks, js, ks, ws, lwd = thickness, col = 'red')
      segments(kss, jss, kss, wss, lwd = thickness, col = 'blue')
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      legend("bottom", c('Raw',paste(mList$mWaveD$shrinkType,' Thresholding',sep='')), xpd = TRUE, horiz = TRUE, 
             inset = -0.1, bty = "n", lty = rep(1,2), col = c('red','blue'), cex = 1)
    }
    
  },  height=function() { session$clientData$output_multiPlot_width * 3/4 })

output$wvdPlot <- renderPlot({
  sList <- sigList()
	mList <- mWaveDList()
  
	mlwvd <- mList$mWaveD
	n <- mList$n
	m <- mList$m
	
	# Number of things to plot
	k = 2
  plotTitle <- paste(input$shrinkage1, ' thresholded estimate of ', sList$signalName,sep="")
	if (ggplot2Avail){
	  out <- c(sList$trueSignal,mlwvd$estimate)
	  lty <- c(rep('solid',n),rep('dashed',n))
	  
	  group <- rep(c("True Signal", "Estimate"), each = n)
    
	  if (input$wvdshow != 3){
	    k <- 3
	    extra <- rep('dashed', n)
	    if (m == 1){
	      out = c(out,multiEstimate(mlwvd$signal , G = mlwvd$G , blur = 'smooth', alpha = mlwvd$alpha, shrinkType = mlwvd$shrinkType, deg = as.integer(input$degree)))
	      if (input$wvdshow == 2){
	        plot.lab <- 'Naive mWaveD average'
	      } else {
	        plot.lab <- 'Best channel'
	      }
	      group <- c(group,rep(plot.lab,n))
	      lty <- c(lty,extra)
	    } else {
	      if (input$wvdshow == 2){
	        wvdMean = apply(sapply(1:m, function(x) 
	          multiEstimate(as.matrix(mlwvd$signal[, x]) ,as.matrix(mlwvd$G[, x]) , blur = 'smooth', alpha = mlwvd$alpha[x], shrinkType = mlwvd$shrinkType, deg = as.integer(input$degree))), 1, mean)
	        out <- c(out,wvdMean)
	        group <- c(group,rep('Naive mWaveD average',n))
	        lty <- c(lty,extra)
	      } else {
	        blur.ty <- mlwvd$blurType
	        if (blur.ty == 'smooth'){
	          best <- mlwvd$blurInfo$bestChannel
	        } else {
	          best <- which.min(mlwvd$sigmaEst)
	        }
          wvdBest <- multiEstimate(as.matrix(mlwvd$signal[, best]), as.matrix(mlwvd$G[, best]), blur = 'smooth', deg = as.integer(input$degree))
	        out <- c(out, wvdBest)
	        group <- c(group, rep("mWaveD best chan", n))
	        lty <- c(lty, extra)
	      }
	    }
	  }
	  est.df = data.frame(Y = out, x = rep(1:n/n, k), group = group, lty = lty)
	  est.plot <- ggplot(est.df, aes(x = x, y = Y, colour = group)) + geom_line(aes(linetype = lty, colour = group), size = lsize + 0.25) + scale_linetype(guide = 'none') + ggtitle(plotTitle)
	  print(est.plot)
	} else {
    if (input$wvdshow != 3){
      k <- 3
      if (input$wvdshow == 2){
        xtra = apply(sapply(1:m, function(x) 
          multiEstimate(as.matrix(mlwvd$signal[, x]) ,as.matrix(mlwvd$G[, x]) , blur = 'smooth', alpha = mlwvd$alpha[x], shrinkType = mlwvd$shrinkType, deg = as.integer(input$degree))), 1, mean)
        xtralab <- 'Naive mean'
      } else {
        blur.ty <- mlwvd$blurType
        if (blur.ty == 'smooth'){
          best = mlwvd$blurInfo$bestChannel
        } else {
          best = which.min(mlwvd$sigmaEst)
        }
        xtra <- multiEstimate(as.matrix(mlwvd$signal[, best]),as.matrix(mlwvd$G[, best]), blur = 'smooth', deg = as.integer(input$degree))
        xtralab <- 'Best channel'
      }
    } else {
      xtra <- NULL
      xtralab <- NULL
    }
    Y <- cbind(sList$trueSignal,mlwvd$estimate, xtra)
    t = (0:(n-1))/n
    labels <- c("True Signal",'mWaveD estimate',xtralab)
    par(oma = c(4, 1, 1, 1))
    matplot(t, Y, type = 'l', xlab = 'x', ylab = '', main = plotTitle)
    grid()
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    legend("bottom", labels, xpd = TRUE, horiz = TRUE, 
           inset = -0.1, bty = "n", lty = 1:3, col = 1:3, cex = 1)
	}
	
	},  height=function() { session$clientData$output_wvdPlot_width * 3/4 })
  
  output$resolutionPlot <- renderPlot({
    mList = mWaveDList()      
    mWaveD <- mList$mWaveD
    n <- mList$n
    m <- mList$m
    iw = -(n/2 - 1):(n/2)
    
    blurInfo = mWaveD$blurInfo
    blurType = mWaveD$blurType
    
    if (blurType != "box.car"){
      blurf <- blurInfo$decay
      cutf <- blurInfo$cutoff
      ymin <- min(cutf)
      revblur <- as.matrix(blurf[(n/2):2, ])
      revcut <- as.matrix(cutf[(n/2):2, ])
      tS <- (1:(n/2))
      
      tS <- c(-rev(tS[-(n/2)]), 0, tS)
      blur <- rbind(revblur, blurf)
      cut <- rbind(revcut, cutf)
      plotTitle <- 'Decay of weighted Kernel coefficients in Fourier domain'
      xlab <- 'Fourier freq'
      if (ggplot2Avail){
        fourier.df <- data.frame(Y = as.vector(blur), x = rep(iw,m), Ycut = as.vector(cut),Channel=rep(LETTERS[1:m],each=n),m=m)
        resolution.plot <- ggplot(fourier.df) + geom_line(aes(x=x, y=Y, colour=Channel,group=Channel),size=hsize) + geom_line(aes(x=x, y=Ycut, colour=Channel), linetype='dashed', size=hsize) + ggtitle(plotTitle) + labs(x = xlab, y = '')
        if (blurType != 'smooth'){
          print(resolution.plot)
        }
      }
      if (blurType == 'smooth'){
        xbest <- max(blurInfo$freqCutoffs) - 1
        ybest <- blurf[xbest, blurInfo$bestChannel]
        xlim <- min(2 * xbest, n/2)
        if (ggplot2Avail){
          high.df <- data.frame(x = rep(xbest, 2), y = c(-Inf, ybest))
          point.df <- data.frame(xbest = xbest, ybest = ybest)
          resolution.plot <- resolution.plot + geom_line(aes(x = x, y = y), linetype = 'dotted' , data = high.df) + geom_line(aes(x = -x + 1, y = y), linetype = 'dotted' , data = high.df) + geom_point( aes(x = xbest, y = ybest), size = 4, shape = 1, data = point.df) + geom_point( aes(x = -xbest + 1, y = ybest), size = 4, shape = 1, data = point.df)
        } 
        ylimlow <- 0
      } else {
        xlim <- n/2
        ylimlow <- 0.5
      }
      xlim <- c(-xlim, xlim)
      l <- n/2 + 1
      ylim <- c(min(cutf[2, ], blurf[l, ]), ylimlow)
      
      zoom = input$zoom
      if (zoom != TRUE){
        xlim = c(-n/2, n/2)
      }
      
      if (ggplot2Avail){
        resolution.plot <- resolution.plot + coord_cartesian(xlim = xlim, ylim = ylim)
        print(resolution.plot)
      } else {
        par(oma = c(4, 1, 1, 1))
        ylim <-  c(min(cutf[2, ], blurf[n/2 + 1, ]), 0.5)
        if (blurType != 'smooth'){
          xlim = c(-n/2, n/2)
        }
        matplot(iw, blur, type = 'l', main = plotTitle, xlab = 'Fourier freq', ylab = '', xlim = xlim, ylim = ylim, lty = 'solid')
        matlines(iw, cut, lty = 'dashed')
        if (blurType == 'smooth'){
          points(xbest, ybest, col='blue')
          points(-xbest + 1, ybest, col = 'blue')
          xbest <- rep(xbest, 2)
          ybest <- c(ylim[1], ybest)
          lines(xbest, ybest, lty = 'dotted')
          lines(-xbest + 1, ybest, lty = 'dotted')
        }
        grid()
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        legend("bottom", paste('Channel ', LETTERS[1:m]), xpd = TRUE, horiz = TRUE, 
               inset = -0.1, bty = "n", lty = 1:4, col = 1:4, cex = 1)
      }
    } else {
      j0   <- mWaveD$j0
      j1   <- mWaveD$j1
      J    <- floor(log2(n)) - 1
      j    <- j0:min(c(J - 1, 2 * j1))
      blkV <- blurInfo$blockVar[1:length(j)]
      blkc <- blurInfo$blockCutoff[1:length(j)]
      ylim <- range(c(blurInfo$blockVar, blurInfo$blockCutoff))
      if (ggplot2Avail){
        resolution.df <- data.frame(Y = c(blkV, blkc), x = rep(j,2), colour = rep(c("Resolution var.",'Resolution bounds'), each = length(j)) , Ycut = blkc)
        bestV <- blkV[j == j1]
        high.df <- data.frame(x = c(j1, j1), y = c(ylim[1], bestV))
        point.df <- data.frame(j1 = j1, bestV = bestV)
        resolution.plot <- ggplot(resolution.df) + geom_line(aes(x = x, y = Y, colour = colour, linetype = colour), size = hsize) +  geom_line(aes(x = x, y = y), linetype = 'dotted', data = high.df) + labs(x = 'j', y = '') + geom_point( aes(x = j1, y = bestV), size = 4, shape = 1, data = point.df)  + scale_color_discrete(labels= c('Resolution bounds', 'Resolution var.'), guide=guide_legend(title.position='left',title.theme = element_text(size=15,angle=0))) + scale_size(guide='none') + guides(colour = guide_legend( title='Blockwise resolution decay')) + theme(legend.position="top", legend.key = element_rect(fill = NA), axis.text.y = element_text(angle = 90)) + scale_linetype_manual(values=c(1,2), name="Blockwise resolution decay", labels=c('Resolution bounds', 'Resolution var.')) + scale_x_continuous(breaks = j)
        print(resolution.plot)
      } else {
        plot(j, blkV, type = 'b', xlab = 'j', ylab = '', main = 'Block wise resolution selection')
        lines(j, blkc, col = 2)
        points(j1, blurInfo$blockVar[j == j1], col='blue')
        lines(c(j1, j1), c(ylim[1], blurInfo$blockVar[j == j1]), lty = 'dashed')
        grid()
        legend("top", c('Resolution bounds', 'Resolution var.'), xpd = TRUE, horiz = TRUE, 
               inset = 0, bty = "n", pch = c(NA_integer_,21) , lty = rep(1,2), col = 2:1, cex = 1.2)
      }
      
    }
  },  height=function() { session$clientData$output_resolutionPlot_width * 3/4 })

	output$summaryOut <- renderPrint({ 
    mList <- mWaveDList()
	  print(summary(mList$mWaveD))
	})
  
  output$summaryWVD <- renderPrint({
  
    wList <- mWaveDList()
    sList <- sigList()
  
    cat("Degree of Meyer poly =", wList$mWaveD$degree, "\n")
    cat("# obs per channel =", wList$n, "\n\n")
    cat("Max possible resolution =", log2(wList$n) - 1, "\n")
    cat("Coarse resolution, j0 =", wList$j0,"\n")
    cat("mWaveD optimal res, j1 =", wList$j1)
    
    cat("\n\n")
    cat("Thresholding method:", wList$mWaveD$shrinkType, "\n")
    cat("Tuning parameter: eta =",wList$mWaveD$eta,'\n\n')
    
    threshMatrix <- cbind(round(wList$mWaveD$levelMax,4), round(wList$mWaveD$thresh,4), round(wList$mWaveD$percent,2))
    rownames(threshMatrix) <- paste("Res", wList$j0:wList$j1,":" )
    colnames(threshMatrix) <- c("max|coef|", "thresh", "%shrunk" )
    print(threshMatrix)
    
    cat('\n Measure of fit:\n')
    cat('L_1 norm =', sum(abs(sList$trueSignal - wList$mWaveD$estimate))/sList$n, '\n')
    cat('L_2 norm =', sqrt(sum((sList$trueSignal - wList$mWaveD$estimate)^2)/sList$n), '\n')
    
  })

  output$summaryResolution <- renderPrint({
    
    wList <- mWaveDList()
    
    cat("# obs per channel =", wList$n, "\n\n")
    cat("# channels: m = ", wList$m, ".\n", sep ='')
    cat("Max possible resolution =", log2(wList$n) - 1, "\n")
    cat("Coarse resolution, j0 =", wList$j0,"\n")
    cat("mWaveD optimal res, j1 =", wList$j1)
    
    cat("\n\n")
    cat("Estimated Channel information:\n\n")
    
    if (wList$mWaveD$blurType == "direct"){
      mat <- cbind(round(wList$mWaveD$sigmaEst, 3), round(wList$mWaveD$alpha, 3), wList$mWaveD$blurInfo$freq, rep(wList$mWaveD$blurInfo$j1, wList$m))
      colnames(mat) <- c("sigma", "alpha", "cutoff", "best res")
      rownames(mat) <- paste("Chan ", 1:wList$m,':', sep = '')
      print(mat)
    } else {
      # If Smooth blur is used, display the matrix of values
      if (wList$mWaveD$blurType == "smooth"){
        mat <- cbind(round(wList$mWaveD$sigmaEst, 3), round(wList$mWaveD$alpha, 3), wList$mWaveD$blurInfo$freq, wList$mWaveD$blurInfo$maxLevels)
        colnames(mat) <- c("sigma", "alpha", "cutoff", "best res")
        rownames(mat) <- paste("chan: ", 1:wList$m)
        print(mat)
        cat("\n")
        cat("Estimated best channel = ", wList$mWaveD$blurInfo$bestChannel)
        
      } else {
        if (wList$mWaveD$blurType == "box.car"){
          mat <- cbind(round(wList$mWaveD$sigma, 3), round(wList$mWaveD$alpha, 3))
          colnames(mat) <- c("sigma", "alpha")
          rownames(mat) <- paste("chan: ", 1:wList$m)
          print(mat)
        }
      }
    }
  })

  output$summaryMRA <- renderPrint({
    
    wList <- mWaveDList()
    
    cat("Degree of Meyer poly =", wList$mWaveD$degree, "\n")
    cat("# obs per channel =", wList$n, "\n\n")
    cat("Max possible resolution =", log2(wList$n) - 1, "\n")
    cat("Coarse resolution, j0 =", wList$j0,"\n")
    cat("mWaveD optimal res, j1 =", wList$j1)
    
    cat("\n\n")
    cat("Thresholding method:", wList$mWaveD$shrinkType, "\n")
    cat("Tuning parameter: eta =",wList$mWaveD$eta,'\n\n')
    
    threshMatrix <- cbind(round(wList$mWaveD$levelMax,4), round(wList$mWaveD$thresh,4), round(wList$mWaveD$percent,2))
    rownames(threshMatrix) <- paste("Res", wList$j0:wList$j1,":" )
    colnames(threshMatrix) <- c("max|coef|", "thresh", "%shrunk" )
    print(threshMatrix)
  })

  output$summarySignal <- renderPrint({
    
    sList <- sigList()
    cat("# obs per channel = ", sList$n, ".\n", sep = '')
    cat("# channels: m = ", sList$m, ".\n", sep ='')
    cat('Blur type: ',sList$blur, '\n\n')
    
    cat("Channel information:\n\n")
    mat <- cbind(round(sList$sigma, 3), round(sList$alpha, 3), round(sList$SNR,3))
    colnames(mat) <- c("sigma", "alpha", "SNR")
    rownames(mat) <- paste("Chan ", 1:sList$m," :", sep ='')
    if ( sList$blur == 'smooth'){
      mat <- cbind(mat, round(sList$shape, 3), round(sList$scale, 3))
      colnames(mat) <- c("sigma", "alpha","SNR", "shape", "scale")
    } else {
      if (sList$blur == 'box.car'){
        mat <- cbind(mat, paste(round(sList$BA, 3), " = 1/sqrt(", 1/sList$BA^2 ,")", sep = '') )
        colnames(mat) <- c("sigma", "alpha","SNR", "Box.car widths")        
      }
    }
    
    print(mat, quote = FALSE)  
  })
  
  output$signalCalls <- renderPrint({
    cat('Base R Function calls:\n\n')
    sList <- sigList()
    cat('J <-', log2(sList$n), '\n')
    cat('n <- 2^J', '\n')
    cat('m <-', sList$m, '\n')
    switch(input$sig,
           'lidar' = cat('signal <- makeLIDAR(n)', '\n'),
           'doppler' = cat('signal <- makeDoppler(n)', '\n'),
           'bumps' = cat('signal <- makeBumps(n)', '\n'),
           'blocks' = cat('signal <- makeBlocks(n)', '\n'))
    cat('t <- (0:(n-1))/n\n')
    if (input$signalShow == '1'){
      cat("plot(t, signal, type = 'l')\n")
    }
    if (input$signalShow == '2' || input$signalShow == '3'){
      if (sList$blur == 'direct'){
        cat('G <- directBlur(c(n, m))\n') 
      }
      if (sList$blur == 'smooth'){
        if (sList$m == 1){
          cat('shape <-', numclean(sList$shape),'\n')
          cat('scale <-', numclean(sList$scale),'\n')
        } else {
          cat('shape <-',xcat(numclean(sList$shape)),'\n', sep = '')
          cat('scale <-',xcat(numclean(sList$scale)),'\n', sep = '')  
        }
        cat('G <- gammaBlur(n, shape, scale)\n') 
      }
      if (sList$blur == 'box.car'){
        if (sList$m == 1){
          cat('boxWindow <- 1/sqrt(', 1/as.numeric(sList$BA)^2,')\n', sep = '')
        } else {
          cat('boxWindow <- 1/sqrt(', xcat(1/as.numeric(sList$BA)^2),')\n', sep = '')
        }
        cat("G <- boxcarBlur(n, boxWindow)\n")
      }
      cat('X <- blurSignal(signal, G)\n')
      if (input$signalShow == '2'){
        cat("matplot(t, X, type = 'l')")
      }
      if (input$signalShow == '3'){
        if (sList$m == 1){
          cat('alpha <-', numclean(sList$alpha),'\n')
          cat('SNR <-', numclean(sList$SNR), '\n')
        } else {
          cat('alpha <-', xcat(numclean(sList$alpha)),'\n')
          cat('SNR <-', xcat(numclean(sList$SNR)), '\n')  
        }
        cat('sigma <- sigmaSNR(X, SNR)\n')
        cat('E <- multiNoise(n, sigma, alpha)\n')
        cat('Y <- X + E\n')
        cat("matplot(t, Y, type = 'l')")
      }
    }
  })

  output$resolutionCalls <- renderPrint({
    cat('Base R Function calls:\n\n')
    
    wList <- mWaveDList()
    cat("mWaveD.output <- multiWaveD(Y, G = G, alpha = alpha, blur = '",wList$mWaveD$blurType,"', shrinkType = '",wList$mWaveD$shrinkType,"', deg = ",wList$mWaveD$degree,')\n', sep = '')
    cat('plot(mWaveD.output, which = 3)')
  })

  output$mWaveDCalls <- renderPrint({
    cat('Base R Function calls:\n\n')

    wList <- mWaveDList()
    cat('# Note the following commands just give the mWaveD estimate only.\n')
    cat("mWaveD.output <- multiWaveD(Y, G = G, alpha = alpha, blur = '",wList$mWaveD$blurType,"', shrinkType = '",wList$mWaveD$shrinkType,"', deg = ",wList$mWaveD$degree,')\n', sep = '')
    cat('plot(mWaveD.output, which = 2)\n')
    cat('## Or alternatively\n')
    cat('x = (0:(n-1))/n\n')
    cat("estimate <- multiEstimate(Y, G = G, alpha = alpha, blur = '",wList$mWaveD$blurType,"', shrinkType = '",wList$mWaveD$shrinkType,"', deg = ",wList$mWaveD$degree,')\n', sep = '')
    cat("plot(x, estimate, type = 'l')")
  })

  output$mraCalls <- renderPrint({
    cat('Base R Function calls:\n\n')
    
    wList <- mWaveDList()
    cat('j1 <- mWaveD.output$j1\n')
    cat("rawOutputCoefs <- mWaveD.output$coef\n")
    cat("shrinkCoefs <- mWaveD.output$shrinkCoef\n")
    cat('plot(rawOutputCoefs, shrinkCoef = shrinkCoefs, highest = j1)\n')
    cat('## Or alternatively\n')
    cat('plot(mWaveD.output, which = 4)')
  })

  observe({
    # This observer depends on shrinkage1 and updates shrinkage2 with any changes
    if (currentShrinkage != input$shrinkage1){
      # Then we assume text2 hasn't yet been updated
      updateSelectInput(session, inputId = 'shrinkage2', label = NULL, choices = shrinkChoices,
                      selected = input$shrinkage1)
      currentShrinkage <<- input$shrinkage1
    }
  })

  observe({
    # This observer depends on shrinkage2 and updates shrinkage1 with any changes
    if (currentShrinkage != input$shrinkage2){
    # Then we assume text2 hasn't yet been updated
      updateSelectInput(session, inputId = 'shrinkage1', label = NULL, choices = shrinkChoices,
                      selected = input$shrinkage2)
      currentShrinkage <<- input$shrinkage2
    }
  })
  updateSelectInput(session, inputId = 'shrinkage1', label = shrinkLabel, choices = shrinkChoices,
                  selected = currentShrinkage)
  updateSelectInput(session, inputId = 'shrinkage2', label = shrinkLabel, choices = shrinkChoices,
                  selected = currentShrinkage)
})


