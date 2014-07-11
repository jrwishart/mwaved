# Auxiliary file that contains the functions used in the simulations.  NB:
# REQUIRES fracdiff package

######### Simulate the LRD noise #########

LRD.noise.sim <- function(n, alpha,Blur,SNR) {
  # LRD.noise.sim: Computes a matrix of M LRD signals of length n using a
  # dependence vector alpha Usage g = LRD.noise.sim(n,alpha) Inputs n
  # length(y) <- 2^J alpha dependence vector of length M (values between 0
  # and 0.5), d = (1-alpha)/2 = H - 1/2
  # 
  # Outputs out M x n matrix of M LRD sequences with dependence level d.
  if(length(alpha) ==1)
  {
    J <- log2(length(Blur))
    y <- sqrt(mean(Blur^2)) * 10^(-SNR/20)
    out <- y*fracdiff.sim(n, d = (1 - alpha)/2)$series
  } else {	
    sds <- sqrt(apply(Blur^2, 2, mean)) * 10^(-SNR/20)
    M <- length(alpha)
    d = (1-alpha)/2
    out <- sapply(1:M, function(x,y) y[x]*fracdiff.sim(n,d = d[x])$series,y=sds)
  }
  out
}

######### Blur the multiple channels #########

BlurSignal <- function(f, g) {
    # BlurSignal: Blur the multiple signals using mvfft Usage g =
    # BlurSignal(f,g) Inputs f matrix of M signals with each with n values as
    # an n x M matrix g matrix of M convolving sequences with each with n
    # values as an n x M matrix
    # 
    # Outputs y M x n matrix of M blurred signals of length n.
    if(is.null(dim(g))){
	n = length(g)
	f <- as.vector(f);
	f.fft <- fft(f)
	g.fft <- fft(g)
	fg.fft <- f.fft * g.fft
	y <- Re(fft(fg.fft, inverse = TRUE))/n
    } else {
	n = dim(f)[1]
	f.fft <- mvfft(f)
	g.fft <- mvfft(g)
	fg.fft <- f.fft * g.fft
	y <- Re(mvfft(fg.fft, inverse = TRUE))/n
    }
}

ScaleNoise <- function(Blur,SNR,alpha) {
    # BlurSignal: Blur the multiple signals using mvfft Usage g =
    # BlurSignal(f,g) Inputs f matrix of M signals with each with n values as
    # an n x M matrix g matrix of M convolving sequences with each with n
    # values as an n x M matrix
    # 
    # Outputs y M x n matrix of M blurred signals of length n.
    if(is.null(dim(Blur))){
	J <- log2(length(Blur))
	y <- sqrt(mean(Blur^2)) * 10^(-SNR/20)/2^(-J*alpha/2)
    } else {
	n <- dim(Blur)[1]
	M <- length(alpha)
	J <- log2(n)
	y <- matrix(rep(sqrt(apply(Blur^2, 2, mean)) * 10^(-SNR/20)/2^(-J*alpha/2),n),ncol=M,byrow=T)
    }
    y
}

#########  Make the regular smooth convolving seq #########

gamma.mat <- function(n, alp, bet) {
    # gamma.mat: Create Gamma density Blur (regular smooth) Usage g =
    # gamma.mat(n,shape,scale) Inputs n length(y) = 2^J shape vector of shape
    # parameters of Gamma distribution (vector of length M) scale vector of
    # scale parameters of Gamma distribution (vector of length M)
    # 
    # Outputs G M x n matrix of gamma matrix for signals.
    
    if (length(alp) != length(bet)) 
        warning("Length of vectors for scale and shape parametrs do not match")
    M <- length(alp)
    if( M > 1){
	G <- sapply((1:n)/n, dgamma, shape = alp, scale = bet)
	G <- t(G)
	aux1 <- matrix(rep(apply(G, 2, max), dim(G)[1]), ncol = M, nrow = n, byrow = T)
	G <- G/aux1
	aux2 <- matrix(rep(apply(G, 2, sum), dim(G)[1]), ncol = M, nrow = n, byrow = T)
	G <- G/aux2
    } else {
	G <- dgamma((1:n)/n,shape=alp,scale=bet);
	aux1 <- max(G)
	G <- G/aux1
	aux2 <- sum(G)
	G <- G/aux2
    }
    G
}

#########  Make the super smooth convolving seq #########

gaussian.mat <- function(n, param) {
    # gaussian.mat: 
    # 
    # 
    # Outputs G M x n matrix of gamma matrix for signals.
    
    if(length(param)==1){
	dom <- seq(0,4,le=n/2+1);
	dom <- c(dom,-rev(dom[c(-1,-n/2-1)]));
	G = dnorm(dom,mean=0,sd=param);
    } else {
	M <- length(param);
        G <- matrix(0,nr=n,ncol=M);
        for(i in 1:M){
            dom <- seq(0,4,le=n/2+1);
            dom <- c(dom,-rev(dom[c(-1,-n/2-1)]));
            G[,i] = dnorm(dom,mean=0,sd=param[i]);
        }
    }
    G
}


######### Make the box.car convolving seq #########
box.car.mat <- function(n, BA) {
    # gamma.mat: Create Gamma density Blur (regular smooth) Usage g =
    # box.car.mat(n,shape,scale) Inputs n length(y) = 2^J BA vector of BA
    # numbers for Box.car (vector of length M)
    # 
    # Outputs G M x n matrix of gamma matrix for signals.
    
    M <- length(BA)
    
    box.car <- function(n, a) {
        left <- seq(from=0,to=1-1/n,le=n)
        left <- (left < a/2) * 1/(a*n)
        aux <- c(0,rev(left[-1]))
        y <- left + aux
#         y <- c(left, rev(left[c(-1, -(n/2 + 1))]))/n
    }
    
#       box.car <- function(n, a) {
#         left <- seq(from=0,to=1-1/n,le=n)
#         s <- sum(left < a/2)
#         d <- 1/(a)
#         left <- (left < a/2) * d/n
#         left <- (left < a/2) * d/n
#         y <- rep(0,n)
#         ind = c(1:s,n:(n-s+1))
#         y[ind] <- d
#         y
# #         y <- c(left, rev(left[c(-1, -(n/2 + 1))]))/n
#       }

    G <- sapply(BA, function(BA) box.car(n, BA))
    if(length(BA) == 1) G <- as.vector(G)
    G
}

######### Make the 'convolution' matrix for direct setting #########
direct.mat <- function(n, M) {
    # direct.mat: Create convolution matrix for direct setting Usage G =
    # direct.mat(n,M) Inputs n length(y) = 2^J M number of channels
    # 
    # Outputs G M x n zero matrix (with first row of ones).
    if( M == 1){
	y <- c(1,rep(0,n-1))
    } else {
	y <- matrix(0, nrow = n, ncol = M)
        y[1, ] <- 1
    }
    y
}


######### Make the signals #########

make.lidar <- function(n, M = 1) {
    # make.bumps: Create M `LIDAR' signals of length n in an n x M matrix Usage
    # f = make.lidar(n,M) Inputs n length(y) = 2^J M number of channels
    # 
    # Outputs y M x n matrix of lidar signals.
    t <- (1:n)/n
    y <- 0.7 * (1 * (t > 0.15) * (t < 0.65) + 1 * (t > 0.28) * (t < 0.48) + (133.33 * 
        t - 106.66) * (t > 0.8) * (t < 0.815) + (-133.33 * t + 110.6639) * (t > 
        0.815) * (t < 0.83) + (133.33 * t - 119.997) * (t > 0.9) * (t < 0.915) + 
        (-133.33 * t + 123.9969) * (t > 0.915) * (t < 0.93))
    if( M > 1) y <- matrix(rep(y, M), ncol = M, byrow = F)
	return(y)
}

make.bumps <- function(n, M = 1) {
    # make.bumps: Create M `bumps' signals of length n in an n x M matrix Usage
    # f = make.bumps(n,M) Inputs n length(y) = 2^J M number of channels
    # 
    # Outputs y M x n matrix of bumps signals.
    
    t <- 1:n/n
    pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
    hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
    wth <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 
        0.005)
    y <- rep(0, n)
    for (j in 1:length(pos)) {
        y <- y + hgt[j]/(1 + abs((t - pos[j])/wth[j]))^4
    }
    y <- y * 1.8
    if( M > 1) y <- matrix(rep(y, M), ncol = M, byrow = F)
	return(y)
}

make.doppler <- function(n, M = 1) {
    # make.doppler: Create M `doppler' signals of length n in an n x M matrix
    # Usage f = make.doppler(n,M) Inputs n length(y) = 2^J M number of channels
    # 
    # Outputs y M x n matrix of doppler signals.
    
    t <- (1:n)/n
    y = (sqrt(t * (1 - t))) * sin((2 * pi * 1.05)/(t + 0.05))
    y <- y * 2.4
    Y <- matrix(rep(y, M), ncol = M, byrow = F)
    if( M > 1) y <- matrix(rep(y, M), ncol = M, byrow = F)
    return(y)
}


make.cusp <- function(n, M = 1) {
    # make.cusp: Create M `cusp' signals of length n in an n x M matrix Usage f
    # = make.cusp(n,M) Inputs n length(y) = 2^J M number of channels
    # 
    # Outputs y M x n matrix of cusp signals.
    
    t <- (1:n)/n
    y <- 2.4 * sqrt(abs(t - 0.37))
    if( M > 1)
	y <- matrix(rep(y, M), ncol = M, byrow = F)
}

make.blocks <- function(n, M = 1) {
    # make.blocks: Create M `blocks' signals of length n in an n x M matrix
    # Usage f = make.blocks(n,M) Inputs n length(y) = 2^J M number of channels
    # 
    # Outputs y M x n matrix of blocks signals.
    
    t <- (1:n)/n
    pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
    hgt = c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
    y <- rep(0, n)
    for (j in 1:length(pos)) {
        y = y + (1 + sign(t - pos[j])) * (hgt[j]/2)
    }
    y
    if( M > 1) y <- matrix(rep(y, M), ncol = M, byrow = F)
    return(y)
}


####### Generate prime numbers for the BA numbers
# Code taken from http://www.r-bloggers.com/prime-number-in-r-language-cloudstat/
prime = function(n)
{
   n = as.integer(n)
   if(n > 1e8) stop("n too large")
   primes = rep(TRUE, n)
   primes[1] = FALSE
   last.prime = 2L
   fsqr = floor(sqrt(n))
   while (last.prime <= fsqr)
   {
      primes[seq.int(2L*last.prime, n, last.prime)] = FALSE
      sel = which(primes[(last.prime+1):(fsqr+1)])
      if(any(sel)){
        last.prime = last.prime + min(sel)
      }else last.prime = fsqr+1
   }
   which(primes)
}
# 3,5,7,11,13,17,19,23,29,31,37,41,43,
primenum <- c(47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199)

# make.cusp <- function(n) { cusppoint <- 0.37; f <- 0.5-cusppoint;
# t<-(1:(n/2))/(n); y <- sqrt(1/4 - t^2); t <- (1:n)/n; aux <- rot90(y); y
# <- c(y,aux); o <- min(which(t >= f)); y <- 2*c(y[o:n],y[1:(o-1)]);
# return(y); } 
