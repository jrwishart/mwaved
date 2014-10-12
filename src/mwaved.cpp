#include <fftw3.h>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

double MeyerPol(double x, int deg){
  double xi;
  switch(deg){
    case 0:  xi = x;
      break;
    case 1:  xi = pow(x, 2) * (3 - 2 * x);
      break;
    case 2:  xi = pow(x, 3) * (10 - 15 * x + 6  * pow(x, 2));
      break;
    case 3:  xi = pow(x, 4) * (35 - 84 * x + 70 * pow(x, 2) - 20 * pow(x, 3));
      break;
    case 4:  xi = pow(x, 5) * (126 - 420 * x + 540 * pow(x, 2) - 315 * pow(x, 3) + 70 * pow(x, 4));
      break;
    default: xi = pow(x, 4) * (35 - 84 * x + 70 * pow(x, 2) - 20 * pow(x, 3));
      break;
  }
  return xi;
}

// Taken from http://gallery.rcpp.org/articles/robust-estimators/
double median_rcpp(NumericVector x) {
   NumericVector y = clone(x);
   int n, half;
   double y1, y2;
   n    = y.size();
   half = n / 2;
   if (n % 2 == 1) {
      // median for odd length vector
      std::nth_element(y.begin(), y.begin() + half, y.end());
      return y[half];
   } else {
      // median for even length vector
      std::nth_element(y.begin(), y.begin() + half, y.end());
      y1 = y[half];
      std::nth_element(y.begin(), y.begin() + half - 1, y.begin() + half);
      y2 = y[half - 1];
      return (y1 + y2) / 2.0;
   }
}

double mad_rcpp(NumericVector x, double scale_factor = 1.4826) {
   // scale_factor = 1.4826; default for normal distribution consistent with R
   return median_rcpp(abs(x - median_rcpp(x))) * scale_factor;
}


NumericVector est_sigma(NumericMatrix noise){
  int j, m;
  m = noise.ncol();
  NumericVector final_out(m);
  for ( j = 0; j < m; ++j ){
    final_out[j] = mad_rcpp(noise(_,j));
  }  
  return final_out;
}

NumericVector est_sigma_from_mat(double * real_out, int n, int m){
  
  int i, j;
  
  NumericVector tmp(n);
  NumericVector final_out(m);
  
  for ( j = 0; j < m; ++j ){
    for ( i = 0; i < n; ++i ){
      tmp[i] = real_out[i + j * n];
    }
    final_out[j] = mad_rcpp(tmp); 
  }
  
  return final_out;
}

NumericMatrix est_noise(double * real_out, int n, int m){
  
  int i, j;
  
  NumericMatrix noiseMat(n,m);
  
  for ( j = 0; j < m; ++j ){
    for ( i = 0; i < n; ++i ){
      noiseMat(i,j) = real_out[i + j * n];
    }
  }
  
  return noiseMat;
}

//[[Rcpp::export]]
NumericVector multiSigma(NumericMatrix signal, int deg = 3){
  
  int i, j, m, n, n2, J, nj, w1, w2, w3;
  
  double x, xi;
  
  // FFTW specific definitions
  double       *x_m_real, p, c, s;
  fftw_complex *x_m_out, *x_m_in;
  fftw_plan     x_m_real_p, x_m_back_p;
  
  n  = signal.nrow();
  m  = signal.ncol();  
  n2 = n/2 + 1;
  J  = log2((double)n);
  
  x_m_real   = (double*)fftw_malloc(sizeof(double) * n * m);
  x_m_out    = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2 * m);
  x_m_in     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n * m);
  x_m_real_p = fftw_plan_many_dft_r2c(1, &n, m, x_m_real, NULL, 1, n, x_m_out , NULL, 1, n2, FFTW_ESTIMATE);
  x_m_back_p = fftw_plan_many_dft_c2r(1, &n, m, x_m_in, NULL, 1, n, x_m_real, NULL, 1, n, FFTW_ESTIMATE);
  
  for ( j = 0; j < m; ++j ){
    for ( i = 0; i < n; ++i ){
      x_m_real[i + j * n] = signal(i,j);
    }
  }
  
  fftw_execute(x_m_real_p);
    
  // Compute fine levels for scale estimation
  j  = J - 1;
  nj = 1 << j;
  w1 = ceil(nj/3.0);
  w2 = ceil(2.0 * nj/3.0);
  w3 = n2 - pow(2.0, j - 3) - 1;
  p  = 1.0/pow((double)n, 1.0/2)/pow(2.0, j/2.0);
  
  memset(x_m_in, 0, sizeof(fftw_complex) * n * m);
  memset(x_m_real, 0, sizeof(double) * n * m);
  // Input convolution matrix to compute FFT, also 
  for (j = 0; j < m; ++j) {
    for (i = w1; i < w2; ++i) {
      xi                       = (double)i/nj;
      x                        = 3 * xi - 1;
      x                        = p * sin(M_PI_2 * MeyerPol(x, deg));
      c                        = cos(M_PI * xi) * x;
      s                        = sin(- M_PI * xi) * x;
      x_m_in[i + j * n][0]     = c * x_m_out[i + j * n2][0] - s * x_m_out[i + j * n2][1];
      x_m_in[i + j * n][1]     = s * x_m_out[i + j * n2][0] + c * x_m_out[i + j * n2][1];
      x_m_in[n - i + j * n][0] = x_m_in[i + j * n][0];
      x_m_in[n - i + j * n][1] = - x_m_in[i + j * n][1];
    }
    for (i = w2; i < w3; ++i) {
      xi                       = (double)i/nj;
      x                        = 3 * xi/2 - 1;
      x                        = p * cos(M_PI_2 * MeyerPol(x, deg));
      c                        = cos(M_PI * xi) * x;
      s                        = sin(- M_PI * xi) * x;
      x_m_in[i + j * n][0]     = c * x_m_out[i + j * n2][0] - s * x_m_out[i + j * n2][1];
      x_m_in[i + j * n][1]     = s * x_m_out[i + j * n2][0] + c * x_m_out[i + j * n2][1];
      x_m_in[n - i + j * n][0] = x_m_in[i + j * n][0];
      x_m_in[n - i + j * n][1] = - x_m_in[i + j * n][1];
    }
    for (i = w3; i < n2; ++i) {
      x_m_in[i + j * n][0]     = - p * x_m_out[i + j * n2][0];
      x_m_in[i + j * n][1]     = - p * x_m_out[i + j * n2][1];
      x_m_in[n - i + j * n][0] = x_m_in[i + j * n][0];
      x_m_in[n - i + j * n][1] = - x_m_in[i + j * n][1];
    }
  }
  
  fftw_execute(x_m_back_p);
  
  NumericVector final_out(m);
  
  NumericVector tmp(n);
  
  for (j = 0; j < m; ++j) {
    for (i = 0; i < n; ++i) {
      tmp[i] = x_m_real[i + j * n];
    }
    final_out[j] = mad_rcpp(tmp);
  }

  fftw_free(x_m_real);
  fftw_free(x_m_out);
  fftw_free(x_m_in);
  
  fftw_destroy_plan(x_m_real_p);
  fftw_destroy_plan(x_m_back_p);
  
  return final_out;
}

// Function that computes the x_fft for the multichannel method
void mlwavedxfft(fftw_complex * x_fft, int m, int n, 
                fftw_complex * x_multi_out, fftw_complex * g_multi_out,
                NumericVector sigma, NumericVector alpha){
  
  memset(x_fft, 0, sizeof(fftw_complex) * n);
  
  int i, j, n2, k;
  
  double  tmp, denom, c, s, eps, powj, x;

  n2 = n/2 + 1;
  c  = 0;
  // First element is seperate due to division by zero issues.
  for (j = 0; j < m; ++j) {
    k = j * n2;
    // denominator calculation (Mod(g_multi_out[1,]))
    denom = g_multi_out[k][0] * g_multi_out[k][0] + g_multi_out[k][1] * g_multi_out[k][1];
    // numerator calculation mean(y_fft*conj(g_multi_out)[1,])
    //real and imaginary parts
    x     = x_multi_out[k][0] * g_multi_out[k][0] + x_multi_out[k][1] * g_multi_out[k][1];
    c    += x/denom;
  }
  x_fft[0][0] = c/m;
  x_fft[0][1] = 0;
  
  for (i = 1; i < n2; ++i) {
    denom = 0;
    c     = 0;
    s     = 0;
    for (j = 0; j < m; ++j) {
      k = i + j * n2;
      // denominator calculation
      eps    = pow((double)n, alpha[j])/(sigma[j] * sigma[j]);
      powj   = pow((double)i, 1 - alpha[j]);
      x      = g_multi_out[k][0] * g_multi_out[k][0] + g_multi_out[k][1] * g_multi_out[k][1];
      denom += eps * powj * x;
      // numerator calculation
      tmp    = powj * eps;
      //real and imaginary parts
      x      = x_multi_out[k][0] * g_multi_out[k][0] + x_multi_out[k][1] * g_multi_out[k][1];
      c     += x * tmp;
      x      = g_multi_out[k][0] * x_multi_out[k][1] - g_multi_out[k][1] * x_multi_out[k][0];
      s     += x * tmp;
    }
    x_fft[i][0] = c/denom;
    x_fft[i][1] = s/denom;
  }
  for (i = n2; i < n; ++i) {
    denom = 0;
    c     = 0;
    s     = 0;
    for (j = 0; j < m; ++j) {
      k      = n - i + j * n2;
      // denominator calculation
      eps    = pow((double)n, alpha[j])/(sigma[j] * sigma[j]);
      powj   = pow((double)i, 1 - alpha[j]);
      x      = g_multi_out[k][0] * g_multi_out[k][0] + g_multi_out[k][1] * g_multi_out[k][1];
      denom += eps * powj * x;
      // numerator calculation
      tmp    = powj * eps;
      //real and imaginary parts
      x      = x_multi_out[k][0] * g_multi_out[k][0] + x_multi_out[k][1] * g_multi_out[k][1];
      c     += x * tmp;
      x      = g_multi_out[k][0] * x_multi_out[k][1] - g_multi_out[k][1] * x_multi_out[k][0];
      s     += x * tmp;
    }
    x_fft[i][0] = c/denom;
    x_fft[i][1] = - s/denom;
  }
    
}

// Function that computes the Hard thresholding given fftw input and threshold level
void hardThreshFFTW(double * in, double * out, int n, double thr){
    int i;
    double x;
    // Hard-Threshold the fft inversion
    for (i = 0; i < n; ++i) {
      x = in[i];
      if (fabs(x) < thr) {
        out[i] = 0;
      } else {
        out[i] = x;
      }
    }
}

// Function that computes the Soft thresholding given fftw input and threshold level
void softThreshFFTW(double * in, double * out, int n, double thr){
    int i;
    double x;
    // Soft-threshold the fft inversion
    for (i = 0; i < n; ++i) {
      x = in[i];
      if (fabs(x) < thr) {
        out[i] = 0;
      } else {
          if (x > 0) {
            out[i] = x - thr;
          } else {
            out[i] = x + thr;
          }
      }
    }
}

// Function that computes the Soft thresholding given fftw input and threshold level
void garroteThreshFFTW(double * in, double * out, int n, double thr){
    int i;
    double x, thr2;
    
    thr2 = thr * thr;
    // Soft-threshold the fft inversion
    for (i = 0; i < n; ++i) {
      x = in[i];
      if (fabs(x) < thr) {
        out[i] = 0;
      } else {
          out[i] = x - thr2/x;
      }
    }
}

void ThresholdFFTW(double * in, double * out, int n, double thr, String shrinkType){
  
  if (shrinkType == "hard") {
    hardThreshFFTW(in, out, n, thr);
  } else if (shrinkType == "soft") {
     softThreshFFTW(in, out, n, thr);
    } else {
      garroteThreshFFTW(in, out, n, thr);
  }
}

// Function that computes the Hard thresholded version of wavelet coefficients
List HardThreshCoef(NumericVector beta, int j0, int j1, NumericVector thr){

  int i, nj, j, k, n, ti, per;
  double x, ax, cur;
  
  n  = beta.size();
  nj = j1 - j0 + 1;
  NumericVector beta_shrink(n);
  NumericVector percent_shrunk(nj);
  NumericVector level_max(nj);
    
  nj = 1 << j0;
  //Coarse level coefficients not shrunk.
  for (i = 0; i < nj; ++i)
    beta_shrink[i] = beta[i];
    
  nj /= 2;
  ti  = -1; 
  // Hard-Threshold the detail wavelet coefficients
  for (j = j0; j <= j1; ++j) {
    nj *= 2;
    k   = nj;
    ti += 1;
    per = 0;
    cur = 0;
    for (i = 0; i < nj; ++i) {
      x  = beta[k];
      ax = fabs(x);
      if (ax > cur) cur = ax;
      if (ax < thr[ti]){
        beta_shrink[k] = 0;
        per++;
      } else {
        beta_shrink[k] = beta[k];
      }
      k++;
    }
    percent_shrunk[j - j0] = 100*(double)per/nj;
    level_max[j - j0]      = cur;
  }
  
  List list_out = List::create(
    _["coef"]    = beta_shrink,
    _["percent"] = percent_shrunk,
    _["max"]     = level_max
          );
  return list_out;
}

// Function that computes the Soft thresholded version of wavelet coefficients
List SoftThreshCoef(NumericVector beta, int j0, int j1, NumericVector thr){

  int i, nj, j, k, n, ti, per;
  double x, ax, cur;
  
  n  = beta.size();
  nj = j1 - j0 + 1;
  NumericVector beta_shrink(n);
  NumericVector percent_shrunk(nj);
  NumericVector level_max(nj);
    
  nj = 1 << j0;
  //Coarse level coefficients not shrunk.
  for (i = 0; i < nj; ++i)
    beta_shrink[i] = beta[i];
    
  nj /= 2;
  ti  = -1; 
  // Hard-Threshold the detail wavelet coefficients
  for (j = j0; j <= j1; ++j) {
    nj *= 2;
    k   = nj;
    ti += 1;
    per = 0;
    cur = 0;
    for (i = 0; i < nj; ++i) {
      x = beta[k];
      ax = fabs(x);
      if (ax > cur) cur = ax;
      if (ax < thr[ti]) {
        beta_shrink[k] = 0;
        per++;
      } else {
        if (x > 0) {
          beta_shrink[k] = x - thr[ti]; 
        } else {
          beta_shrink[k] = x + thr[ti]; 
        }
      }
      k++;
    }
    percent_shrunk[j - j0] = (double)per/nj;
    level_max[j - j0]      = cur;
  }
  
  List list_out = List::create(
    _["coef"]    = beta_shrink,
    _["percent"] = percent_shrunk,
    _["max"]     = level_max
          );
  return list_out;
}

// Function that computes the Garrote thresholded version of wavelet coefficients
List GarroteThreshCoef(NumericVector beta, int j0, int j1, NumericVector thr){

  int i, nj, j, k, n, ti, per;
  double x, ax, thr2, cur;
  
  n = beta.size();
  nj = j1 - j0 + 1;
  NumericVector beta_shrink(n);
  NumericVector percent_shrunk(nj);
  NumericVector level_max(nj);
    
  nj = 1 << j0;
  //Coarse level coefficients not shrunk.
  for (i = 0; i < nj; ++i)
    beta_shrink[i] = beta[i];
    
  nj /= 2;
  ti  = -1; 
  // Hard-Threshold the detail wavelet coefficients
  for (j = j0; j <= j1; ++j) {
    nj *= 2;
    k   = nj;
    ti += 1;
    per = 0;
    cur = 0;
    thr2 = thr[ti] * thr[ti];
    for (i = 0; i < nj; ++i) {
      x = beta[k];
      ax = fabs(x);
      if (ax > cur) cur = ax;
      if (ax < thr[ti]) {
        beta_shrink[k] = 0;
        per++;
      } else {
        beta_shrink[k] = x - thr2/x;
      }
      k++;
    }
    percent_shrunk[j - j0] = (double)per/nj;
    level_max[j - j0]      = cur;
  }
  
  List list_out = List::create(
    _["coef"]    = beta_shrink,
    _["percent"] = percent_shrunk,
    _["max"]     = level_max
          );
  return list_out;
}

// Function that dispatches the wavelet coefficient shrinkage type 
List ThresholdCoef(NumericVector beta, int j0, int j1, NumericVector thr, String shrinkType){

  List beta_shrink(beta.size());

  if (shrinkType == "hard") {
    beta_shrink = HardThreshCoef(beta, j0, j1, thr);
  } else if (shrinkType == "soft") {
      beta_shrink = SoftThreshCoef(beta, j0, j1, thr);
      // Otherwise garrote thresholding used.
    } else {
        beta_shrink = GarroteThreshCoef(beta, j0, j1, thr);
  }
  return beta_shrink;
}

// Exported wavelet thresholding functions below
NumericVector hardThresh(NumericVector in, NumericVector thr, int j0, int j1){
    int i, j, k, n, l;
    double x, thr1;
    
    NumericVector out(in.size());
    
    n = 1 << j0;
    for (i = 0; i < n; ++i) {
      out[i] = in[i];
    }
    
    // Start thresholding at correct spot
    k = (1 << j0);
    l = 0;
    
    for (j = j0; j <= j1; ++j) {
      thr1 = thr[l];
      for (i = 0; i < n; ++i) {
          x = in[k];
        if (fabs(x) < thr1) {
          out[k] = 0;
        } else {
          out[k] = x;
        }
        ++k;
      }
      ++l;
      n *= 2;
    }
    
    return out;
}

// Function that computes the Soft thresholding given coefficients and threshold level
NumericVector softThresh(NumericVector in, NumericVector thr, int j0, int j1){
    int i, j, k, n, l;
    double x, thr1;
    
    NumericVector out(in.size());
    
    n = 1 << j0;
    for (i = 0; i < n; ++i) {
      out[i] = in[i];
    }
    
    // Start thresholding at correct spot
    k = n - 1;
    l = 0;
    
    for (j = j0; j <= j1; ++j) {
      thr1 = thr[l];
      for (i = 0; i < n; ++i) {
        x = in[k];
        if (fabs(x) < thr1) {
          out[k] = 0;
        } else {
          if (x > 0) {
            out[k] = x - thr1;
          } else {
            out[k] = x + thr1;
          }
        }
        ++k;
      }
      ++l;
      n *= 2;
    }
    
    return out;
}

// Function that computes the Garrote thresholding given coefficients and threshold level
NumericVector garroteThresh(NumericVector in, NumericVector thr, int j0, int j1){
    int i, j, k, n, l;
    double x, thr1, thr2;
    
    NumericVector out(in.size());
    
    n = 1 << j0;
    for (i = 0; i < n; ++i) {
      out[i] = in[i];
    }
    
    // Start thresholding at correct spot
    k = n - 1;
    l = 0;
    
    for (j = j0; j <= j1; ++j) {
      thr1 = thr[l];
      thr2 = thr1 * thr1;
      for (i = 0; i < n; ++i) {
        x = in[k];
        if (fabs(x) < thr1) {
          out[k] = 0;
        } else {
          out[k] = x - thr2/x;
        }
        ++k;
      }
      ++l;
      n *= 2;
    }
    
    return out;
}

//[[Rcpp::export]]
List waveletThresh(NumericVector beta, NumericVector thr, String shrinkType = "hard", int j0 = 3, int deg = 3){
  
  int n   = beta.size();
  int nt  = thr.size();
  int j1  = j0 + nt - 1;
  
  NumericVector betaShrink(n);
  
  if (shrinkType == "hard") {
    betaShrink = hardThresh(beta, thr, j0, j1);
  } else if (shrinkType == "soft") {
      betaShrink = softThresh(beta, thr, j0, j1);
    } else {
      betaShrink = garroteThresh(beta, thr, j0, j1);
  }
    
  List betaThresh =  List::create(
    _["coef"] = betaShrink,
    _["j0"]   = j0,
    _["deg"]  = deg
    );
  
  betaThresh.attr("class") = "waveletCoef";
  
  return betaThresh;
}

// Function that computes the boxcar feasible level information from fft knowledge
List BoxCarChanInfo(int m, int n, fftw_complex * g_fft, NumericVector sigma, NumericVector alpha, int j0, int deg){

  int i, j, j1, n2, w1, w2, w3, J, nj, k, jmax;
  
  double x, xi;
  double eps_cut = 0;

  n2     = n/2 + 1;
  J      = log2((double)n);
  jmax   = J - 1;
  NumericVector eps(m);
  NumericVector finfo(n2);
  NumericVector blockvar(J - j0);
  NumericVector blockcut(J - j0);
  
// Compute multi-channel threshold
  for (j = 0; j < m; ++j) {
    eps[j]   = pow((double)n, alpha[j])/(sigma[j] * sigma[j]);
    eps_cut += log(eps[j]);
  }
  
  // Compute the normalised Fourier decay levels.
  for (i = 1; i < n2; ++i) {
      for (j = 0; j < m; ++j)
        finfo[i] += ((g_fft[i + j * n2][0] * g_fft[i + j * n2][0]) + (g_fft[i + j * n2][1] * g_fft[i + j * n2][1])) * eps[j] * pow((double)i, 1 - alpha[j]);
  }
  
  nj = pow(2.0, j0 - 1);
  k  = -1;  
  for (j = j0; j < jmax; ++j) {
    ++k;
    nj *= 2;
    w1  = ceil(nj/3.0);
    w2  = 2 * w1 + j % 2 -1;
    w3  = w1 + nj;
    for (i = w1; i < w2; ++i) {
      xi           = (double)i/nj;
      x            = 3 * xi - 1;
      xi           = sin(M_PI_2 * MeyerPol(x, deg));
      blockvar[k] += xi * xi / finfo[i];
    }
    for (i = w2; i < w3; ++i) {
      xi           = (double)i/nj;
      x            = 3 * xi/2 - 1;
      xi           = cos(M_PI_2 * MeyerPol(x, deg));
      blockvar[k] += xi * xi / finfo[i];
    }
    blockvar[k]   /= nj;
    blockvar[k]    = log( blockvar[k] );
    blockcut[k]    = - log( nj * eps_cut );
  }
  
  // finest level case
  ++k;
  nj *= 2;
  w1  = ceil(nj/3.0);
  w2  = 2 * w1 + j % 2 -1;
  w3 = n2 - pow(2.0, j - 3) - 1;
  for (i = w1; i < w2; ++i) {
    xi           = (double)i/nj;
    x            = 3 * xi - 1;
    xi           = sin(M_PI_2 * MeyerPol(x, deg));
    blockvar[k] += xi * xi / finfo[i];
  }
  for (i = w2; i < w3; ++i) {
    xi           = (double)i/nj;
    x            = 3 * xi/2 - 1;
    xi           = cos(M_PI_2 * MeyerPol(x, deg));
    blockvar[k] += xi * xi / finfo[i];
  }
  for (i = w3; i < n2; ++i){
    blockvar[k] += 1.0 / finfo[i];
  }
  blockvar[k]   /= nj;
  blockvar[k]    = log( blockvar[k] );
  blockcut[k]    = - log( nj * eps_cut );
  
  k  = -1;
  for (j = j0; j < J; ++j) {
    ++k;
    if (blockvar[k] > blockcut[k]) {
      break;
    }
  }
  j1 = j - 1;
  
  return List::create(
    _["finfo"]       = finfo,
    _["blockVar"]    = blockvar,
    _["blockCutoff"] = blockcut,
    _["j0"]          = j0,
    _["finest"]      = j1
    );
}

// Function that computes the blur info for Direct case
List DirectChanInfo(int m, int n, NumericVector sigma, NumericVector alpha){

  int i, j, n2;
  
  double tmp, nsqrt, tmp2;
  
  NumericVector threshfft(m);
  NumericVector level(m);
  IntegerVector freq(m, NA_INTEGER);
  
  n2  = n/2 + 1;
  int J = log2((double)n);
  
  NumericMatrix finfo(n2, m);
  NumericMatrix fcut(n2, m);
  nsqrt = pow((double)n, 1.0/2);
  
  // Compute channel level thresholds
  for (j = 0; j < m; ++j) {
    tmp          = log(sigma[j]);
    threshfft[j] = tmp - log(pow((double)n, alpha[j]/2.0)) + 0.5 * log(fabs(log(nsqrt) - tmp));
  } 
  // Compute the cutoffs for each channel
  for (j = 0; j < m; ++j) {
    tmp  = 0.5 * alpha[j];
    tmp2 = threshfft[j];
    fcut(0, j) = R_NegInf;
    for (i = 1; i < n2; ++i) {
      fcut(i,j)  = tmp * log((double)i) + tmp2;
    }
    freq[j]  = n2 - 1;
  }
  // Return the info
  return List::create(
    _["decay"]       = finfo,
    _["cutoffs"]     = fcut,
    _["freqCutoffs"] = freq,
    _["finest"]      = J - 1
    );
}

// Function that computes the feasible levels from fft knowledge
List SmoothChanInfo(int m, int n, fftw_complex * g_fft, NumericVector sigma, NumericVector alpha){

  int i, j, n2, best;
  
  double tmp, nsqrt, tmp2;
  
  NumericVector threshfft(m);
  IntegerVector level(m);
  IntegerVector freq(m, NA_INTEGER);
  
  n2  = n/2 + 1;
  
  NumericMatrix finfo(n2, m);
  NumericMatrix fcut(n2, m);
  nsqrt = pow((double)n, 1.0/2);
  
  // Compute channel level thresholds
  for (j = 0; j < m; ++j) {
    tmp          = log(sigma[j]);
    threshfft[j] = tmp - log(pow((double)n, alpha[j]/2.0)) + 0.5 * log(fabs(log(nsqrt) - tmp));
  } 
  // Compute the cutoffs for each channel
  for (j = 0; j < m; ++j) {
    tmp  = 0.5 * alpha[j];
    tmp2 = threshfft[j];
    fcut(0,j) = R_NegInf;
    for (i = 1; i < n2; ++i) {
      finfo(i,j) = log(sqrt((g_fft[i + j * n2][0] * g_fft[i + j * n2][0]) + (g_fft[i + j * n2][1] * g_fft[i + j * n2][1])));
      fcut(i,j)  = tmp * log((double)i) + tmp2;
    }
  }
  
  for (j = 0; j < m; ++j) {
    for (i = 1; i < n2; ++i) {
      tmp = finfo(i, j) - 0.5 * alpha[j] * log((double)i);
      if (tmp < threshfft[j]) {
        freq[j]  = i;
        level[j] = floor(log2((double)i + 1.0)) - 1;
        break;
      }
    }
    // Threshold never met (direct convolution case)
    if (freq[j] == NA_INTEGER) {
      freq[j]  = n2 - 1;
      level[j] = floor(log2((double)freq[j]));
    }
  }
  // Find the best channel from fourier info
  best = 0;
  int current = freq[0];
  for (j = 1; j < m; ++j) {
    if (freq[j] > current) {
      current = freq[j];
      best    = j;
    } 
  }
  
  // Return the info
  return List::create(
    _["maxLevels"]   = level,
    _["decay"]       = finfo,
    _["cutoffs"]     = fcut,
    _["freqCutoffs"] = freq,
    _["bestChannel"] = best + 1,
    _["finest"]      = level[best]
  );
}

// Function that computes the feasible levels from fft knowledge
int FindBestChannel(int m, int n, fftw_complex * g_fft, NumericVector sigma, NumericVector alpha){

  int i, j, n2, best;
  
  double tmp, nsqrt;
  
  NumericVector threshfft(m);
  NumericVector level(m);
  NumericVector freq(m, -1.0);
  n2    = n/2 + 1;
  nsqrt = pow((double)n, 1.0/2);
  
  
  // Compute channel level thresholds
  for (j = 0; j < m; ++j) {
    tmp   = log(sigma[j]);
    threshfft[j] = tmp - log(pow((double)n, alpha[j]/2.0)) + 0.5 * log(fabs(log(nsqrt) - tmp));
  }
  // Compute the normalised Fourier decay levels.
  for (j = 0; j < m; ++j) {
    // Skip first level since Inf always > threshfft 
    for (i = 1; i < n2; ++i) {
      tmp  = log(sqrt((g_fft[i + j * n2][0] * g_fft[i + j * n2][0]) + (g_fft[i + j * n2][1] * g_fft[i + j * n2][1]))) - 0.5 * alpha[j] * log((double)i);
      if (tmp < threshfft[j]) {
        freq[j]  = i + 1;
        level[j] = floor(log2((double)i + 1)) - 1;
        break;
      }
    }
    // Threshold never met (direct convolution case)
    if (freq[j] == -1.0) {
      freq[j]  = n2 - 1;
      level[j] = log2((double)freq[j]);
    }
    
  }
  
  best = 1;
  int current = freq[0];
  for (j = 1; j < m; ++j) {
    if (freq[j] > current){
      current = freq[j];
      best    = j + 1;
    } 
  }
  
  return best;
}

// Compute theoretical Eta if it is not specified
double TheoreticalEta(NumericVector alpha, String resolution, int m, int n,
                      fftw_complex * g_multi_out, NumericVector sigma){
  
  double eta;

  int best_chan = 1;
  if (resolution == "smooth") {
    best_chan = FindBestChannel(m, n, g_multi_out, sigma, alpha);
    eta = 4 * sqrt(alpha[best_chan - 1]);
  } else {
    eta = 4 * sqrt(min(alpha));
  }
  
  return eta;
}

// Function that only computes the highest scale level
int HighestScale(int m, int n, fftw_complex * g_fft, NumericVector sigma,
                NumericVector alpha, String resolution = "smooth", int j0 = 3, int deg = 3){

  int i, j, n2, j1;
  double tmp;
  
  n2 = n/2 + 1;
  
  
  if (resolution == "smooth"){
    NumericVector threshfft(m);
    NumericVector level(m);
    IntegerVector freq(m, NA_INTEGER);
    int best = 1;
    double nsqrt = pow((double)n, 1.0/2); 
    // Compute channel level thresholds
    for (j = 0; j < m; ++j) {
      tmp   = log(sigma[j]);
      threshfft[j] = tmp - log(pow((double)n, alpha[j]/2.0)) + 0.5 * log(fabs(log(nsqrt) - tmp));
    }
    // Compute the normalised Fourier decay levels.
    for (j = 0; j < m; ++j) {
      // Skip first level since Inf always > threshfft
      for (i = 1; i < n2; ++i) {
        tmp  = log(sqrt((g_fft[i + j * n2][0] * g_fft[i + j * n2][0]) + (g_fft[i + j * n2][1] * g_fft[i + j * n2][1]))) - 0.5 * alpha[j] * log((double)i);
        if (tmp < threshfft[j]) {
          freq[j]  = i + 1;
          level[j] = floor(log2((double)i + 1)) - 1;
          break;
        }
      }
      // Threshold never met (direct convolution case)
      if (freq[j] == NA_INTEGER) {
        freq[j]  = n2 - 1;
        level[j] = log2((double)freq[j]);
      }
      
    }
    int current = freq[0];
    for (j = 1; j < m; ++j) {
      if (freq[j] > current) {
        current = freq[j];
        best    = j + 1;
      } 
    }
    j1 = level[best - 1];
  } else {
    int J = log2((double)n2);
    NumericVector eps(m);
    NumericVector finfo(n2);
    double x, xi;
    double eps_cut = 0;
    tmp     = 0;
    bool cond      = false;
    int w1, w2, w3, nj;
    nj = pow(2.0, j0 - 1);
    
    for (j = 0; j < m; ++j) {
      eps[j]   = pow((double)n, alpha[j])/(sigma[j] * sigma[j]);
      eps_cut += log(eps[j]);
    }
    
    // Compute the normalised Fourier decay levels.
    for (i = 1; i < n2; ++i) {
        for (j = 0; j < m; ++j)
          finfo[i] += ((g_fft[i + j * n2][0] * g_fft[i + j * n2][0]) + (g_fft[i + j * n2][1] * g_fft[i + j * n2][1])) * eps[j] * pow((double)i, 1 - alpha[j]);
    }
    
    nj = 1 << (j0 - 1); 
    for (j = j0; j < J; ++j) {
      nj *= 2;
      w1  = ceil(nj/3.0);
      w2  = 2 * w1 + j % 2 -1;
      w3  = w1 + nj;
      
      for (i = w1; i < w2; ++i) {
        xi   = (double)i/nj;
        x    = 3 * xi - 1;
        xi   = sin(M_PI_2 * MeyerPol(x, deg));
        tmp += xi * xi / finfo[i];
      }
      for (i = w2; i < w3; ++i) {
        xi   = (double)i/nj;
        x    = 3 * xi/2 - 1;
        xi   = cos(M_PI_2 * MeyerPol(x, deg));
        tmp += xi * xi / finfo[i];
      }
      if (tmp > 1/eps_cut) {
        cond = true;
        break;
      }
    }
    j1 = j - 1;
    // If necessary, check finest resolution.
    if (cond == false) {
      nj *= 2;
      w1  = ceil(nj/3.0);
      w2  = 2 * w1 + j % 2 -1;
      w3  = w1 + nj;
      w2  = 2 * w1 + j % 2 -1;
      w3 = n2 - pow(2.0, j - 3) - 1;
      for (i = w1; i < w2; ++i) {
        xi   = (double)i/nj;
        x    = 3 * xi - 1;
        xi   = sin(M_PI_2 * MeyerPol(x, deg));
        tmp += xi * xi / finfo[i];
      }
      for (i = w2; i < w3; ++i) {
        xi   = (double)i/nj;
        x    = 3 * xi/2 - 1;
        xi   = cos(M_PI_2 * MeyerPol(x, deg));
        tmp += xi * xi / finfo[i];
      }
      for (i = w3; i < n2; ++i){
        tmp += 1.0 / finfo[i];
      }
      if (tmp > 1/eps_cut) {
        ++j1;
      }
    }
  }
  
  
  return j1;
}

// Returns the Wavelet projection using the coefficients and scales as inputs.
// [[Rcpp::export]]
NumericVector multiProj(NumericVector beta, int j0 = 3, int j1 = NA_INTEGER, int deg = 3){
  
  int n, n2, nj, J, i, j, jmax, k, ks, kj, w1, w2, w3;
  double x, xi, cx, sx, p;

  n  = beta.size();
  n2 = n/2 + 1;
  J  = log2((double)n);
  
  if (j1 == NA_INTEGER)
    j1 = J - 1;
  
  jmax = j1;

   // FFTW specific definitions
  double *real_in, *real_out;
  fftw_complex *in, *out;
  fftw_plan real_p, back_p;

  real_in  = (double*)fftw_malloc(sizeof(double) * n);
  real_out = (double*)fftw_malloc(sizeof(double) * n);
  in       = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  out      = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2);

  real_p   = fftw_plan_dft_r2c_1d(n, real_in, out, FFTW_PATIENT);
  back_p   = fftw_plan_dft_c2r_1d(n, in, real_out, FFTW_PATIENT);
  
  NumericVector final_out(n);
  
  memset(real_in, 0, sizeof (double) * n);
  memset(out, 0, sizeof (fftw_complex) * n2);
  
  nj = 1 << j0;
  ks = floor((double)n/nj);
  cx = 1.0/sqrt((double)nj);
  w1 = ceil((double)nj/3);
  w2 = 2 * w1 + j0 % 2 - 1;
  k  = - ks;
  for (i = 0; i < nj; ++i) {
    k         += ks;
    real_in[k] = beta[i];
  }

  fftw_execute(real_p);
  
  memset(in, 0, sizeof (fftw_complex) * n);
  
  in[0][0] = cx * out[0][0];
  
  for (i = 1; i < w1; i++) {
    in[i][0]     = cx * out[i][0];
    in[i][1]     = cx * out[i][1];
    in[n - i][0] = cx * out[n2 - i][0];
    in[n - i][1] = - cx * out[n2 - i][1];
  }

  for (i = w1; i < w2; i++) {
    xi           = (double)i/nj;
    x            = 3 * xi - 1;
    xi           = cx * cos(M_PI_2 * MeyerPol(x, deg));
    in[i][0]     = xi * out[i][0];
    in[i][1]     = xi * out[i][1];
    in[n - i][0] = xi * out[n2 - i][0];
    in[n - i][1] = - xi * out[n2 - i][1];
  }

  memset(real_out, 0, sizeof (double) * n);
  
  fftw_execute(back_p);
  
  for (i = 0; i < n; ++i) {
    final_out[i] = real_out[i];
  }
  
  // If needed add the finest details
  if (j1 == J - 1) {
    
    memset(real_in, 0, sizeof (double) * n);
    memset(out, 0, sizeof (fftw_complex) * n2);
  
    nj = 1 << j1;
    ks = floor((double)n/nj);
    p  = 1.0/sqrt((double)nj);
    w1 = ceil((double)nj/3);
    w2 = 2 * w1 + j0 % 2 - 1;
    w3 = n2 - pow(2.0, j1 - 3) - 1;
    k  = - ks;
    kj = nj;
    for (i = 0; i < nj; ++i) {
      k         += ks;
      real_in[k] = beta[kj];
      kj++;
    }
    
    fftw_execute(real_p);
    
    memset(in, 0, sizeof (fftw_complex) * n);
    
    for (i = w1; i < w2; ++i) {
      xi           = (double)i/nj;
      x            = 3 * xi - 1;
      x            = p * sin(M_PI_2 * MeyerPol(x, deg));
      cx           = cos(M_PI * xi) * x;
      sx           = sin(- M_PI * xi) * x;
      in[i][0]     = cx * out[i][0] - sx * out[i][1];
      in[i][1]     = sx * out[i][0] + cx * out[i][1];
      in[n - i][0] = cx * out[n2 - i][0] + sx * out[n2 - i][1];
      in[n - i][1] = - sx * out[n2 - i][0] + cx * out[n2 - i][1];
    }
    
    for (i = w2; i < w3; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi/2 - 1;
      x              = p * cos(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(- M_PI * xi) * x;
      in[i][0]       = cx * out[i][0] - sx * out[i][1];
      in[i][1]       = sx * out[i][0] + cx * out[i][1];
      in[n - i][0]   = cx * out[n2 - i][0] + sx * out[n2 - i][1];
      in[n - i][1]   = - sx * out[n2 - i][0] + cx * out[n2 - i][1];
    }
    
    n2 -= 1;
    for (i = w3; i < n2; ++i) {
      in[i][0]       = - p * out[i][0];
      in[i][1]       = - p * out[i][1];
      in[n - i][0]   = - p * out[n2 - i][0];
      in[n - i][1]   = - p * out[n2 - i][1];
    }
    in[n2][0]   = - p * out[n2][0];
    in[n2][1]   = - p * out[n2][1];
    n2 += 1;
    
    memset(real_out, 0, sizeof (double) * n);
    
    fftw_execute(back_p);
    
    for (i = 0; i < n; ++i) {
      final_out[i] += real_out[i];
    }
  } else {
    jmax++;
  }
  
  memset(real_in, 0, sizeof (double) * n);
  memset(out, 0, sizeof (fftw_complex) * n2);
  
  
  for (j = j0; j < jmax; ++j) {
    nj = 1 << j;
    ks = floor((double)n/nj);
    p  = 1.0/sqrt((double)nj);
    w1 = ceil((double)nj/3);
    w2 = 2 * w1 + j0 % 2 - 1;
    w3 = w1 + nj;
    k  = - ks;
    kj = nj;
    for (i = 0; i < nj; ++i) {
      k         += ks;
      real_in[k] = beta[kj];
      kj++;
    }
    
    fftw_execute(real_p);
    
    memset(in, 0, sizeof (fftw_complex) * n);
    
    for (i = w1; i < w2; ++i) {
      xi           = (double)i/nj;
      x            = 3 * xi - 1;
      x            = p * sin(M_PI_2 * MeyerPol(x, deg));
      cx           = cos(M_PI * xi) * x;
      sx           = sin(- M_PI * xi) * x;
      in[i][0]     = cx * out[i][0] - sx * out[i][1];
      in[i][1]     = sx * out[i][0] + cx * out[i][1];
      in[n - i][0] = cx * out[n2 - i][0] + sx * out[n2 - i][1];
      in[n - i][1] = - sx * out[n2 - i][0] + cx * out[n2 - i][1];
    }
    
    for (i = w2; i < w3; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi/2 - 1;
      x              = p * cos(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(- M_PI * xi) * x;
      in[i][0]       = cx * out[i][0] - sx * out[i][1];
      in[i][1]       = sx * out[i][0] + cx * out[i][1];
      in[n - i][0]   = cx * out[n2 - i][0] + sx * out[n2 - i][1];
      in[n - i][1]   = - sx * out[n2 - i][0] + cx * out[n2 - i][1];
    }
    
    memset(real_out, 0, sizeof (double) * n);
    
    fftw_execute(back_p);
    
    for (i = 0; i < n; ++i) {
      final_out[i] += real_out[i];
    }
  }
  
  fftw_free(in);
  fftw_free(out);
  fftw_free(real_in);
  fftw_free(real_out);
  
  fftw_destroy_plan(real_p);
  fftw_destroy_plan(back_p);

  return final_out;
}

// [[Rcpp::export]]
NumericVector multiThresh(NumericMatrix signal, NumericMatrix G, NumericVector alpha = NumericVector::create(),
                          String resolution = "smooth", int j0 = 3, int j1 = NA_INTEGER, double eta = NA_REAL, int deg = 3) {
  
  int i, m, n, n2, j, J, nj, jmax, w1, w2, w3;
  double tmpd, x, xi;
  
    // FFTW specific definitions
  double       *g_m_real_in;
  fftw_complex *g_multi_out;
  fftw_plan     g_m_real_p;
  
  n  = signal.nrow();
  n2 = n/2 + 1;
  m  = signal.ncol();
  
  // Set default alpha to be 1 uniformly ("No dependence")
  if (alpha.size() == 0)
    alpha = rep(1.0, m);

  // Check if inputs agree
  if (!( m == alpha.size() && m == G.ncol() && n == G.nrow() ))
    stop("Dimension mismatch; signal, alpha and G");
  
  g_m_real_in  = (double*)fftw_malloc(sizeof(double) * n * m);
  g_multi_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2 * m);
  g_m_real_p   = fftw_plan_many_dft_r2c(1, &n, m, g_m_real_in, NULL, 1, n, g_multi_out, NULL, 1, n2, FFTW_ESTIMATE);
  
  for (j = 0; j < m; ++j) {
    for (i = 0; i < n; ++i) {
      g_m_real_in[i + j * n] = G(i,j);
    }
  }
  
  fftw_execute(g_m_real_p);
  
  NumericVector sigma = multiSigma(signal, deg);
  
  J = log2((double)n);
  
  // if j1 not specified, compute all possible levels
  if (j1 == NA_INTEGER)
    j1 = J - 1;
  
  // If eta not specified, use Theoretical value
  if (R_IsNA(eta))
    eta = TheoreticalEta(alpha, resolution, m, n, g_multi_out, sigma);
  
  // Power factors n^alpha/sigma^2
  NumericVector thrMat(n);
  
  // compute smallest required values
  nj = 1 << j0;
  w1 = ceil((double)nj/3.0);
  w3 = n - w1;
  
  for (i = w1; i < n2; ++i) {
    tmpd = 0;
    for (j = 0; j < m; ++j) {
       x     = pow((double)i, 1 - alpha[j]);
       xi    = pow((double)n, alpha[j])/(sigma[j] * sigma[j]);
       tmpd += x * ((g_multi_out[i + j * n2][0] * g_multi_out[i + j * n2][0]) + (g_multi_out[i + j * n2][1] * g_multi_out[i + j * n2][1])) * xi;
    }
    thrMat[i] = 1.0/tmpd;
  }
  for (i = n2; i < w3; ++i) {
    tmpd = 0;
    for (j = 0; j < m; ++j) {
       x     = pow((double)i, 1 - alpha[j]);
       xi    = pow((double)n, alpha[j])/pow((double)sigma[j], 2);
       tmpd += x * ((g_multi_out[n - i + j * n2][0] * g_multi_out[n - i + j * n2][0]) + (g_multi_out[n - i + j * n2][1] * g_multi_out[n - i + j * n2][1])) * xi;
    }    
    thrMat[i] = 1.0/tmpd;
  }
  
  NumericVector thr(j1 - j0 + 1);
  
  // Check if last level needed and alter
  if (j1 == J - 1) {
    jmax = j1 - 1;
  } else {
    jmax = j1;
  }
  nj = 1 << (jmax + 1);
  
  for (j = jmax; j >= j0; --j) {
    nj /= 2;
    w1 = ceil(nj/3.0);
    w2 = 2 * w1 + j % 2 -1;
    w3 = w1 + nj;
    
    for (i = w1; i < w2; ++i) {
      xi           = (double)i/nj;
      x            = 3 * xi - 1;
      xi           = sin(M_PI_2 * MeyerPol(x, deg));
      x            = 1.0/nj * xi * xi;
      thr[j - j0] += x * thrMat[i] + x * thrMat[n - i];
    }
    for (i = w2; i < w3; ++i) {
      xi           = (double)i/nj;
      x            = 3 * xi/2 - 1;
      xi           = cos(M_PI_2 * MeyerPol(x, deg));
      x            = 1.0/nj * xi * xi;
      thr[j - j0] += x * thrMat[i] + x * thrMat[n - i];
    }
    thr[j - j0] = sqrt(log((double)n) * eta * thr[j - j0]);
  }

  // Linear extrapolate the finest level if required
  if (j1 == J - 1) {
    tmpd         = thr[jmax - j0] - thr[jmax - j0 - 1];
    thr[j1 - j0] = thr[jmax - j0] + tmpd;
  }
  
  fftw_free(g_m_real_in);
  fftw_free(g_multi_out);
  
  fftw_destroy_plan(g_m_real_p);
  
  return thr;
}


// Helper function to compute thresholds fast (n, sigma, g_multi_out, alpha, j0, j1, eta, deg)
NumericVector MaxiThreshFFTW(int n, NumericVector sigma, fftw_complex * g_fft, NumericVector alpha,
                          int j0, int j1, double eta , int deg = 3) {
  
  int i, m, n2, j, J, nj, jmax, w1, w2, w3;
  double tmpd, x, xi;
  
  J = log2((double)n);
  n2 = n/2 + 1;
  m  = sigma.size();
  // Power factors n^alpha/sigma^2
  NumericVector thrMat(n);
  
  // compute smallest required values
  nj = 1 << j0;
  w1 = ceil((double)nj/3.0);
  w3 = n - w1;
  
  for (i = w1; i < n2; ++i) {
    tmpd = 0;
    for (j = 0; j < m; ++j) {
       x     = pow((double)i, 1 - alpha[j]);
       xi    = pow((double)n, alpha[j])/(sigma[j] * sigma[j]);
       tmpd += x * ((g_fft[i + j * n2][0] * g_fft[i + j * n2][0]) + (g_fft[i + j * n2][1] * g_fft[i + j * n2][1])) * xi;
    }
    thrMat[i] = 1.0/tmpd;
  }
  for (i = n2; i < w3; ++i) {
    tmpd = 0;
    for (j = 0; j < m; ++j) {
       x     = pow((double)i, 1 - alpha[j]);
       xi    = pow((double)n, alpha[j])/(sigma[j] * sigma[j]);
       tmpd += x * ((g_fft[n - i + j * n2][0] * g_fft[n - i + j * n2][0]) + (g_fft[n - i + j * n2][1] * g_fft[n - i + j * n2][1])) * xi;
    }    
    thrMat[i] = 1.0/tmpd;
  }
  
  NumericVector thr(j1 - j0 + 1);
  
  // Check if last level needed and alter
  if (j1 == J - 1) {
    jmax = j1 - 1;
  } else {
    jmax = j1;
  }
  nj = 1 << (jmax + 1);
  
  for (j = jmax; j >= j0; --j) {
    nj /= 2;
    w1 = ceil(nj/3.0);
    w2 = 2 * w1 + j % 2 -1;
    w3 = w1 + nj;
    
    for (i = w1; i < w2; ++i) {
      xi           = (double)i/nj;
      x            = 3 * xi - 1;
      xi           = sin(M_PI_2 * MeyerPol(x, deg));
      x            = 1.0/nj * xi * xi;
      thr[j - j0] += x * thrMat[i] + x * thrMat[n - i];
    }
    for (i = w2; i < w3; ++i) {
      xi           = (double)i/nj;
      x            = 3 * xi/2 - 1;
      xi           = cos(M_PI_2 * MeyerPol(x, deg));
      x            = 1.0/nj * xi * xi;
      thr[j - j0] += x * thrMat[i] + x * thrMat[n - i];
    }
    thr[j - j0] = sqrt(log((double)n) * eta * thr[j - j0]);
  }
  
  // Linear extrapolate the finest level if required
  if (j1 == J - 1) {
    tmpd         = thr[jmax - j0] - thr[jmax - j0 - 1];
    thr[j1 - j0] = thr[jmax - j0] + tmpd;
  }
  
  return thr;
}

// [[Rcpp::export]]
NumericVector multiEstimate(NumericMatrix signal, NumericMatrix G, 
                            NumericVector alpha = NumericVector::create(), 
                            String resolution = "smooth", String blur = "direct", NumericVector sigma = NumericVector::create(),
                            int j0 = 3, int j1 = NA_INTEGER, double eta = NA_REAL, 
                            NumericVector thresh = NumericVector::create(), String shrinkType = "hard", int deg = 3){

  // Regular definitions
  int i, j, J, n2, w1, w2, w3, n, jmax, m, nj;
  
  n  = signal.nrow();
  m  = signal.ncol();
  n2 = n/2 + 1;
  J  = log2((double)n);
  
  if (alpha.size() == 0)
    alpha = rep(1.0,m);
    
  
  if (!( m == alpha.size() && m == G.ncol() && n == G.nrow() ))
    stop("Dimension mismatch; sigma, alpha and G");
    
  double xi, x, p, cx, sx;
  
  // R objects
  NumericVector final_out(n);
  
  // FFTW specific definitions
  double *x_m_real_in, *g_m_real_in, *real_in, *real_out, *sig_real_out;
  fftw_complex *in, *out, *x_fft, *conj, *x_multi_out, *g_multi_out, *sig_in;
  fftw_plan x_m_real_p, g_m_real_p, real_p, backward_p, sigma_back_p;

  x_m_real_in  = (double*)fftw_malloc(sizeof(double) * n * m);
  sig_real_out = (double*)fftw_malloc(sizeof(double) * n * m);
  g_m_real_in  = (double*)fftw_malloc(sizeof(double) * n * m);
  real_in      = (double*)fftw_malloc(sizeof(double) * n);
  real_out     = (double*)fftw_malloc(sizeof(double) * n);
  in           = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  x_multi_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2 * m);
  sig_in       = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n * m);
  g_multi_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2 * m);
  out          = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2);
  // Complex FFTW Vars
  x_fft        = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  conj         = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);  
  // plans
  x_m_real_p   = fftw_plan_many_dft_r2c(1, &n, m, x_m_real_in, NULL, 1, n, x_multi_out, NULL, 1, n2, FFTW_ESTIMATE);
  g_m_real_p   = fftw_plan_many_dft_r2c(1, &n, m, g_m_real_in, NULL, 1, n, g_multi_out, NULL, 1, n2, FFTW_ESTIMATE);
  real_p       = fftw_plan_dft_r2c_1d(n, real_in, out, FFTW_PATIENT);
  // Need to determine appropriate c2r transform plan;
  backward_p   = fftw_plan_dft_c2r_1d(n, in, real_out,FFTW_PATIENT);
  
  // Need an appropriate c2r multi plan for scale estimation;
  sigma_back_p = fftw_plan_many_dft_c2r(1, &n, m, sig_in, NULL, 1, n, sig_real_out , NULL, 1, n, FFTW_ESTIMATE);

  // Initialise all to zero
  memset(x_m_real_in, 0, sizeof(double) * n * m);
  memset(sig_real_out, 0, sizeof(double) * n * m);
  memset(g_m_real_in, 0, sizeof(double) * n * m);
  memset(real_in, 0, sizeof(double) * n);
  memset(real_out, 0, sizeof(double) * n);
  memset(in, 0, sizeof(fftw_complex) * n);
  memset(x_multi_out, 0, sizeof(fftw_complex) * n2 * m);
  memset(sig_in, 0, sizeof(fftw_complex) * n * m);
  memset(g_multi_out, 0, sizeof(fftw_complex) * n2 * m);
  memset(out, 0, sizeof(fftw_complex) * n2);
  memset(x_fft, 0, sizeof(fftw_complex) * n);
  memset(conj, 0, sizeof(fftw_complex) * n);
  
  for ( j = 0; j < m; ++j ){
    for ( i = 0; i < n; ++i ){
      x_m_real_in[i + j * n] = signal(i,j);
    }
  }
  
  fftw_execute(x_m_real_p);
  
  // Compute fine levels for scale estimation
  j  = J - 1;
  nj = 1 << j;
  w1 = ceil(nj/3.);
  w2 = ceil(2.*nj/3.);
  w3 = n2 - pow(2.0, j - 3) - 1;
  p  = 1.0/pow((double)n, 1.0/2)/pow(2.0, j/2.0);
  
  // Input convolution matrix to compute FFT, also 
  for ( j = 0; j < m; ++j) {
    for ( i = 0; i < n; ++i) {
      g_m_real_in[i + j * n] = G(i,j);
    }
  }
  
  fftw_execute(g_m_real_p);
  
  j = sigma.size();
  if (j == 0) {
    sigma = multiSigma(signal,deg);
  } else {
    if (j != m) {
      stop("Dimension mismatch : Length of sigma must match the number of columns of Y and G");
    }
  }
  
  if (blur == "direct")
    j1 = J - 1;
  
  if (j1 == NA_INTEGER)
    j1 = HighestScale(m, n, g_multi_out, sigma, alpha, resolution, j0, deg);
  
  if (j1 < j0) {
    j0 = j1;
    Rf_warning("Warning: j1 < j0, highest scale reset to j1 = j0");
  }
  
  // Compute theoretical Eta if it is not specified
  if (R_IsNA(eta)) {
    eta = TheoreticalEta(alpha, resolution, m, n, g_multi_out, sigma);
  }
  // Build thresholds if needed
  if (thresh.size() == 0) 
    thresh = MaxiThreshFFTW(n, sigma, g_multi_out ,alpha ,j0 ,j1 ,eta , deg);
  
  // Compute the multichannel normalised Fourier coefficients
  mlwavedxfft(x_fft, m, n, x_multi_out, g_multi_out, sigma, alpha);
  
  // Coarse scale approximation (no thresholding)
  j = j0;
  
  x  = 1.0/n;
  nj = 1 << j;
  w1 = ceil((double)nj/3);
  w2 = 2 * w1 + j % 2 - 1;
  p  = pow(2.0, 1 - j0) * M_PI;
  cx = x/pow(2.0, j/2.0);
  
  in[0][0] = x * x_fft[0][0];
  
  for (i = 1; i < w1; i++) {
    in[i][0]     = x * x_fft[i][0];
    in[i][1]     = x * x_fft[i][1];
    in[n - i][0] = x * x_fft[n - i][0];
    in[n - i][1] = x * x_fft[n - i][1];
  }

  for (i = w1; i < w2; i++) {
    xi           = (double)i/nj;
    x            = 3 * xi - 1;
    xi           = pow(cos(M_PI_2 * MeyerPol(x, deg)), 2.0)/n;
    in[i][0]     = xi * x_fft[i][0];
    in[i][1]     = xi * x_fft[i][1];
    in[n - i][0] = xi * x_fft[n - i][0];
    in[n - i][1] = xi * x_fft[n - i][1];
  }

  fftw_execute(backward_p);
    
  // Start wavelet expansion and record coarse approximation
  for (i = 0; i < n; ++i)
    final_out[i] = real_out[i];
  
  // Check if final J-1 level needed then add it.
  if (j1 == J - 1) {
    jmax = j1;
    j    = J - 1;
    
    memset(in, 0, sizeof (fftw_complex) * n);
    memset(conj, 0, sizeof (fftw_complex) * n);
    
    nj = 1 << j;
    w1 = ceil((double)nj/3.0);
    w2 = 2 * w1 + j % 2 - 1;
    w3 = n2 - pow(2.0, j - 3) - 1;
    p  = 1.0/n/pow(2.0, j/2.0);
    
    for (i = w1; i < w2; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi - 1;
      x              = p * sin(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(- M_PI * xi) * x;
      in[i][0]       = cx * x_fft[i][0] - sx * x_fft[i][1];
      conj[i][0]     = cx;
      in[i][1]       = sx * x_fft[i][0] + cx * x_fft[i][1];
      conj[i][1]     = - sx;
      in[n - i][0]   = cx * x_fft[n - i][0] + sx * x_fft[n - i][1];
      conj[n - i][0] = cx;
      in[n - i][1]   = - sx * x_fft[n - i][0] + cx * x_fft[n - i][1];
      conj[n - i][1] = sx;
    }
    
    for (i = w2; i < w3; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi/2 - 1;
      x              = p * cos(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(- M_PI * xi) * x;
      in[i][0]       = cx * x_fft[i][0] - sx * x_fft[i][1];
      conj[i][0]     = cx;
      in[i][1]       = sx * x_fft[i][0] + cx * x_fft[i][1];
      conj[i][1]     = - sx;
      in[n - i][0]   = cx * x_fft[n - i][0] + sx * x_fft[n - i][1];
      conj[n - i][0] = cx;
      in[n - i][1]   = - sx * x_fft[n - i][0] + cx * x_fft[n - i][1];
      conj[n - i][1] = sx;
    }
    
    n2 -= 1;
    for (i = w3; i < n2; ++i) {
      in[i][0]       = - p * x_fft[i][0];
      in[i][1]       = - p * x_fft[i][1];
      conj[i][0]     = - p;
      in[n - i][0]   = - p * x_fft[n - i][0];
      in[n - i][1]   = - p * x_fft[n - i][1];
      conj[n - i][0] = - p;
    }
    in[n2][0]   = - p * x_fft[n2][0];
    in[n2][1]   = - p * x_fft[n2][1];
    conj[n2][0] = - p;
    n2 += 1;
    
    // Invert the product
    fftw_execute(backward_p);

    // Threshold the fft inversion
    ThresholdFFTW(real_out, real_in, n, thresh[j1 - j0], shrinkType);
    
    // fft forward
    fftw_execute(real_p);

    // Project back
    for (i = 0; i < n2; ++i) {
      in[i][0] = out[i][0] * conj[i][0] - out[i][1] * conj[i][1];
      in[i][1] = out[i][0] * conj[i][1] + out[i][1] * conj[i][0];
    }
    for (i = n2; i < n; ++i) {
      in[i][0] = out[n - i][0] * conj[i][0] + out[n - i][1] * conj[i][1];
      in[i][1] = out[n - i][0] * conj[i][1] - out[n - i][1] * conj[i][0];
    }

    fftw_execute(backward_p);

    for (i = 0; i < n; ++i)
      final_out[i] += real_out[i] * nj;
  } else {
    jmax = j1 + 1;
  }
  
  
  // Set factors first to simplify loop adjustments
  j   = j0 - 1;
  nj  = 1 << j;
  p   = 1.0/n/pow((double)nj, 1.0/2.0);
  
  for (j = j0; j < jmax; ++j) {
    
    // Perform product of x_fft and psij_fft (Keep Conj(wavj for later))
    memset(in, 0, sizeof (fftw_complex) * n);
    memset(conj, 0, sizeof (fftw_complex) * n);
    
    p  /= pow(2.0, 0.5);
    nj *= 2;
    i   = 1 << j;
    w1  = ceil(i/3.0);
    w2  = 2 * w1 + j % 2 - 1;
    w3  = w1 + nj;
  
    for (i = w1; i < w2; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi-1;
      x              = p * sin(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(-M_PI * xi) * x;
      in[i][0]       = cx * x_fft[i][0] - sx * x_fft[i][1];
      conj[i][0]     = cx;
      in[i][1]       = sx * x_fft[i][0] + cx * x_fft[i][1];
      conj[i][1]     = - sx;
      in[n - i][0]   = cx * x_fft[n - i][0] + sx * x_fft[n - i][1];
      conj[n - i][0] = cx;
      in[n - i][1]   = - sx * x_fft[n - i][0] + cx * x_fft[n - i][1];
      conj[n - i][1] = sx;
    }
    
    for (i = w2; i < w3; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi/2 - 1;
      x              = p * cos(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(- M_PI * xi) * x;
      in[i][0]       = cx * x_fft[i][0] - sx * x_fft[i][1];
      conj[i][0]     = cx;
      in[i][1]       = sx * x_fft[i][0] + cx * x_fft[i][1];
      conj[i][1]     = -sx;
      in[n - i][0]   = cx * x_fft[n - i][0] + sx * x_fft[n - i][1];
      conj[n - i][0] = cx;
      in[n - i][1]   = - sx * x_fft[n - i][0] + cx * x_fft[n - i][1];
      conj[n - i][1] = sx;
    }
    
    // Invert the product
    fftw_execute(backward_p);

    // Threshold the fft inversion
    ThresholdFFTW(real_out, real_in, n, thresh[j - j0], shrinkType);
    
    // fft forward
    fftw_execute(real_p);

    // Project back
    for (i = 0; i < n2; ++i) {
      in[i][0] = out[i][0] * conj[i][0] - out[i][1] * conj[i][1];
      in[i][1] = out[i][0] * conj[i][1] + out[i][1] * conj[i][0];
    }
    for (i = n2; i < n; ++i) {
      in[i][0] = out[n - i][0] * conj[i][0] + out[n - i][1] * conj[i][1];
      in[i][1] = out[n - i][0] * conj[i][1] - out[n - i][1] * conj[i][0];
    }
    fftw_execute(backward_p);
    for (i = 0; i < n; ++i)
      final_out[i] += real_out[i] * nj;
  }

  fftw_destroy_plan(real_p);
  fftw_destroy_plan(x_m_real_p);
  fftw_destroy_plan(g_m_real_p);
  fftw_destroy_plan(backward_p);
  fftw_destroy_plan(sigma_back_p);

  fftw_free(in);
  fftw_free(x_m_real_in);
  fftw_free(g_m_real_in);
  fftw_free(real_in);
  fftw_free(out);
  fftw_free(real_out);
  fftw_free(x_multi_out);
  fftw_free(g_multi_out);
  fftw_free(x_fft);
  fftw_free(conj);
  fftw_free(sig_real_out);
  fftw_free(sig_in);
  
  return final_out;
}


// [[Rcpp::export]]
List multiCoef(NumericMatrix signal, NumericMatrix G, NumericVector alpha = NumericVector::create(),
                          int j0 = 3, int j1 = NA_INTEGER, int deg = 3){

  // Regular definitions
  int i, j, J, k, bk, n, n2, w1, w2, w3, jmax, m, nj;
  
  n  = signal.nrow();
  m  = signal.ncol();
  n2 = n/2 + 1;
  J  = log2((double)n);
  NumericVector beta(n);  
  
  if (alpha.size() == 0)
    alpha = rep(1.0,m);
  
  double xi, x, c, s, p, bp, nf, cx, sx;
  
  // R objects
  NumericVector sigma(m);
  
  // FFTW specific definitions
  double *x_m_real_in, *g_m_real_in, *sig_real_out;
  fftw_complex *x_fft, *x_multi_out, *g_multi_out, *sig_in;
  fftw_plan x_m_real_p, g_m_real_p, sigma_back_p;

  x_m_real_in  = (double*)fftw_malloc(sizeof(double) * n * m);
  sig_real_out = (double*)fftw_malloc(sizeof(double) * n * m);
  g_m_real_in  = (double*)fftw_malloc(sizeof(double) * n * m);
  x_multi_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2 * m);
  sig_in       = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n * m);
  g_multi_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2 * m);
  // Complex FFTW Vars
  x_fft        = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  // plans
  x_m_real_p   = fftw_plan_many_dft_r2c(1, &n, m, x_m_real_in, NULL, 1, n, x_multi_out, NULL, 1, n2, FFTW_ESTIMATE);
  g_m_real_p   = fftw_plan_many_dft_r2c(1, &n, m, g_m_real_in, NULL, 1, n, g_multi_out, NULL, 1, n2, FFTW_ESTIMATE);
  // Need to determine appropriate c2r transform plan;
  
  // Need an appropriate c2r multi plan for scale estimation;
  sigma_back_p = fftw_plan_many_dft_c2r(1, &n, m, sig_in, NULL, 1, n, sig_real_out , NULL, 1, n, FFTW_ESTIMATE);
  
  memset(sig_real_out, 0, sizeof(double) * n * m);
  memset(g_m_real_in, 0, sizeof(double) * n * m);
  memset(x_multi_out, 0, sizeof(fftw_complex) * n2 * m);
  memset(sig_in, 0, sizeof(fftw_complex) * n * m);
  memset(g_multi_out, 0, sizeof(fftw_complex) * n2 * m);
  
  for (j = 0; j < m; ++j) {
    for (i = 0; i < n; ++i) {
      x_m_real_in[i + j * n] = signal(i,j);
    }
  }
  
  fftw_execute(x_m_real_p);
  
  // Compute fine levels for scale estimation
  j  = J - 1;
  nj = 1 << j;
  w1 = ceil(nj/3.0);
  w2 = ceil(2.0 * nj/3.0);
  w3 = n2 - pow(2.0, j - 3) - 1;
  p  = 1.0/pow((double)n, 1.0/2)/pow(2.0, j/2.0);
  
  // Input convolution matrix to compute FFT, also 
  for (j = 0; j < m; ++j) {
    g_m_real_in[j * n] = G(0, j);
    for (i = 1; i < w1; ++i) {
      g_m_real_in[i + j * n]     = G(i, j);
      g_m_real_in[n - i + j * n] = G(n - i,j);
    }
    for (i = w1; i < w2; ++i) {
      g_m_real_in[i + j * n]     = G(i, j);
      g_m_real_in[n - i + j * n] = G(n - i, j);
      xi                         = (double)i/nj;
      x                          = 3 * xi - 1;
      x                          = p * sin(M_PI_2 * MeyerPol(x, deg));
      c                          = cos(M_PI * xi) * x;
      s                          = sin(- M_PI * xi) * x;
      sig_in[i + j * n][0]       = c * x_multi_out[i + j * n2][0] - s * x_multi_out[i + j * n2][1];
      sig_in[i + j * n][1]       = s * x_multi_out[i + j * n2][0] + c * x_multi_out[i + j * n2][1];
      sig_in[n - i + j * n][0]   = sig_in[i + j * n][0];
      sig_in[n - i + j * n][1]   = - sig_in[i + j * n][1];
    }
    for (i = w2; i < w3; ++i) {
      g_m_real_in[i + j * n]     = G(i, j);
      g_m_real_in[n - i + j * n] = G(n - i, j);
      xi                         = (double)i/nj;
      x                          = 3 * xi/2 - 1;
      x                          = p * cos(M_PI_2 * MeyerPol(x, deg));
      c                          = cos(M_PI * xi) * x;
      s                          = sin(- M_PI * xi) * x;
      sig_in[i + j * n][0]       = c * x_multi_out[i + j * n2][0] - s * x_multi_out[i + j * n2][1];
      sig_in[i + j * n][1]       = s * x_multi_out[i + j * n2][0] + c * x_multi_out[i + j * n2][1];
      sig_in[n - i + j * n][0]   = sig_in[i + j * n][0];
      sig_in[n - i + j * n][1]   = - sig_in[i + j * n][1];
    }
    for (i = w3; i < n2; ++i) {
      g_m_real_in[i + j * n]     = G(i, j);
      g_m_real_in[n - i + j * n] = G(n - i, j);
      sig_in[i + j * n][0]       = - p * x_multi_out[i + j * n2][0];
      sig_in[i + j * n][1]       = - p * x_multi_out[i + j * n2][1];
      sig_in[n - i + j * n][0]   = sig_in[i + j * n][0];
      sig_in[n - i + j * n][1]   = - sig_in[i + j * n][1];
    }
  }
  
  fftw_execute(g_m_real_p);
  
  fftw_execute(sigma_back_p);
  
  // Estimate sigma and noise
  
  sigma = est_sigma_from_mat(sig_real_out, n, m);
  
  if (j1 == NA_INTEGER){
    j1 = J - 1;
  }
  
  if (j1 < j0) {
    j1 = j0;
  }
  
  // Compute the multichannel normalised Fourier coefficients
  mlwavedxfft(x_fft, m, n, x_multi_out, g_multi_out, sigma, alpha);

  // Coarse scale approximation (no thresholding)
  j = j0;
  
  x  = 1.0/n;
  nf = x / pow(2.0, j/2.0);
  nj = 1 << j;
  w1 = ceil((double)nj/3);
  w2 = 2 * w1 + j % 2 - 1;
  p  = pow(2.0,1 - j0) * M_PI;
  
  for (k = 0; k < nj; ++k)
    beta[k] = nf * x_fft[0][0];
  
  for (i = 1; i < w1; i++) {
    beta[0]   += nf * (x_fft[i][0] + x_fft[n - i][0]);
    for (k = 1; k < nj; ++k){
      c        = nf * cos(p * i * k);
      s        = nf * sin(p * i * k);
      beta[k] += c  * x_fft[i][0] - s * x_fft[i][1];
      c        = nf * cos(p * (n - i) * k);
      s        = nf * sin(p * (n - i) * k);
      beta[k] += c  * x_fft[n - i][0] - s * x_fft[n - i][1];
    }
  }

  for (i = w1; i < w2; i++) {
    xi = (double)i/nj;
    x  = 3 * xi - 1;
    xi = cos(M_PI_2 * MeyerPol(x, deg));
    x  = nf * xi;
    // Only real part survives for first beta coef
    beta[0] += x * x_fft[i][0] + x * x_fft[n - i][0];
    for (k = 1; k < nj; ++k) {
      c        = x * cos(p * i * k);
      s        = x * sin(p * i * k);
      beta[k] += c * x_fft[i][0] - s * x_fft[i][1];
      c        = x * cos(p * (n - i) * k);
      s        = x * sin(p * (n - i) * k);
      beta[k] += c * x_fft[n - i][0] - s * x_fft[n - i][1];
    }
  }
 
  // Check j1 values are feasible
  if (j1 >= J) {
    Rf_warning("j1 too large, set to maximum feasible  j1 = log2(n) - 1");
    j1 = J - 1;
  }
  if (j1 < j0 ) {
    Rf_warning("j1 too small (must be larger than j0), j1 set equal to j0");
    j1 = j0;
  }
  // Check if final J-1 level needed then add it.
  if (j1 == J - 1) {
    jmax = j1;
    j    = J - 1;
    
    nj = 1 << j;
    w1 = ceil((double)nj/3.0);
    w2 = 2 * w1 + j % 2 - 1;
    w3 = n2 - pow(2.0, j - 3) - 1;
        
    p  = 1.0/n/pow(2.0, j/2.0);
    bp = 4 * M_PI / n;
    
    for (i = w1; i < w2; ++i) {
      xi = (double)i/nj;
      x  = 3 * xi - 1;
      x  = p * sin(M_PI_2 * MeyerPol(x, deg));
      cx = x * cos(M_PI * xi);
      sx = x * sin( - M_PI * xi);
      c  = x_fft[i][0] * cx + x_fft[i][1] * sx;
      s  = cx * x_fft[i][1] - sx * x_fft[i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = x_fft[n - i][0] * cx - x_fft[n - i][1] * sx;
      s  = cx * x_fft[n - i][1] + sx * x_fft[n - i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
    
    for (i = w2; i < w3; ++i) {
      xi = (double)i/nj;
      x  = 3 * xi/2 - 1;
      x  = p * cos(M_PI_2 * MeyerPol(x, deg));
      cx = x * cos(M_PI * xi);
      sx = x * sin( - M_PI * xi); 
      c  = x_fft[i][0] * cx + x_fft[i][1] * sx;
      s  = cx * x_fft[i][1] - sx * x_fft[i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = x_fft[n - i][0] * cx - x_fft[n - i][1] * sx;
      s  = cx * x_fft[n - i][1] + sx * x_fft[n - i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
    
    n2 -= 1;
    
    for (i = w3; i < n2; ++i) {
      c  = - p * x_fft[i][0];
      s  = - p * x_fft[i][1];  
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = - p * x_fft[n - i][0];
      s  = - p * x_fft[n - i][1];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
    // Add the last term
    bk = nj;
    for (k = 0; k < nj; ++k) {
      c         = - p * x_fft[n2][0];
      s         = - p * x_fft[n2][1];
      beta[bk] += c * cos(bp * n2 * k) - s * sin(bp * n2 * k);
      bk++;
    }
    n2 += 1;
    
  } else {
    jmax = j1 + 1;
  }
  
  // Set factors first to simplify loop adjustments
  j   = j0 - 1;
  nj  = 1 << j;
  p   = 1.0/n/pow((double)nj, 0.5);
  bp  = M_PI * pow(2.0, 2 - j0);
  
  for (j = j0; j < jmax ; ++j) {
    
    p  /= pow(2.0, 0.5);
    bp /= 2;
    nj *= 2;
    w1  = ceil((double)nj/3.0);
    w2  = 2 * w1 + j % 2 - 1;
    w3  = w1 + nj;
    
    for (i = w1; i < w2; ++i) {
      xi = (double)i/nj;
      x  = 3 * xi - 1;
      x  = p * sin(M_PI_2 * MeyerPol(x, deg));
      cx = x * cos(M_PI * xi);
      sx = x * sin( - M_PI * xi);
      c  = x_fft[i][0] * cx + x_fft[i][1] * sx;
      s  = cx * x_fft[i][1] - sx * x_fft[i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = x_fft[n - i][0] * cx - x_fft[n - i][1] * sx;
      s  = cx * x_fft[n - i][1] + sx * x_fft[n - i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
    
    for (i = w2; i < w3; ++i) {
      xi = (double)i/nj;
      x  = 3 * xi/2 - 1;
      x  = p * cos(M_PI_2 * MeyerPol(x, deg));
      cx = x * cos(M_PI * xi);
      sx = x * sin( - M_PI * xi); 
      c  = x_fft[i][0] * cx + x_fft[i][1] * sx;
      s  = cx * x_fft[i][1] - sx * x_fft[i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = x_fft[n - i][0] * cx - x_fft[n - i][1] * sx;
      s  = cx * x_fft[n - i][1] + sx * x_fft[n - i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
  }
  
  fftw_destroy_plan(x_m_real_p);
  fftw_destroy_plan(g_m_real_p);
  fftw_destroy_plan(sigma_back_p);

  fftw_free(x_m_real_in);
  fftw_free(g_m_real_in);
  fftw_free(x_multi_out);
  fftw_free(g_multi_out);
  fftw_free(x_fft);
  fftw_free(sig_real_out);
  fftw_free(sig_in);

  List betaList = List::create(
    _["coef"] = beta,
    _["j0"]   = j0,
    _["deg"]  = deg
  );
  
  betaList.attr("class") = "waveletCoef";

  return betaList;
}

// [[Rcpp::export]]
List multiWaveD(NumericMatrix signal, NumericMatrix G, NumericVector alpha = NumericVector::create(),
                    String resolution = "smooth", String blur = "direct", int j0 = 3, int j1 = NA_INTEGER,
                    double eta = NA_REAL, NumericVector thresh = NumericVector::create(),
                    String shrinkType = "hard", int deg = 3){

  // Regular definitions
  int i, j, J, k, bk, n, n2, w1, w2, w3, jmax, m, nj;
  
  n  = signal.nrow();
  m  = signal.ncol();
  n2 = n/2 + 1;
  J  = log2((double)n);
  
  if (alpha.size() == 0) {
    alpha = rep(1.0, m);
  }
  
  if (!( m == alpha.size() && m == G.ncol())) {
    stop("Dimension mismatch; sigma, alpha and G");
  }
  
  double xi, x, c, s, p, bp, cx, sx;
  
  // R objects
  NumericVector final_out(n);
  NumericVector beta(n);
  NumericVector sigma(m);
  NumericMatrix noise(n,m);
  List channel;
  
  // FFTW specific definitions
  double *x_m_real_in, *g_m_real_in, *real_in, *real_out, *sig_real_out;
  fftw_complex *in, *out, *x_fft, *conj, *x_multi_out, *g_multi_out, *sig_in;
  fftw_plan x_m_real_p, g_m_real_p, real_p, backward_p, sigma_back_p;

  x_m_real_in  = (double*)fftw_malloc(sizeof(double) * n * m);
  sig_real_out = (double*)fftw_malloc(sizeof(double) * n * m);
  g_m_real_in  = (double*)fftw_malloc(sizeof(double) * n * m);
  real_in      = (double*)fftw_malloc(sizeof(double) * n);
  real_out     = (double*)fftw_malloc(sizeof(double) * n);
  in           = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  x_multi_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2 * m);
  sig_in       = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n * m);
  g_multi_out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2 * m);
  out          = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2);
  // Complex FFTW Vars
  x_fft      = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  conj       = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);  
  // plans
  x_m_real_p = fftw_plan_many_dft_r2c(1, &n, m, x_m_real_in, NULL, 1, n, x_multi_out, NULL, 1, n2, FFTW_ESTIMATE);
  g_m_real_p = fftw_plan_many_dft_r2c(1, &n, m, g_m_real_in, NULL, 1, n, g_multi_out, NULL, 1, n2, FFTW_ESTIMATE);
  real_p     = fftw_plan_dft_r2c_1d(n, real_in, out, FFTW_ESTIMATE);
  // Need to determine appropriate c2r transform plan;
  backward_p = fftw_plan_dft_c2r_1d(n, in, real_out,FFTW_ESTIMATE);
  
  // Need an appropriate c2r multi plan for scale estimation;
  sigma_back_p = fftw_plan_many_dft_c2r(1, &n, m, sig_in, NULL, 1, n, sig_real_out , NULL, 1, n, FFTW_ESTIMATE);
  
  memset(sig_real_out, 0, sizeof(double) * n * m);
  memset(g_m_real_in, 0, sizeof(double) * n * m);
  memset(real_in, 0, sizeof(double) * n);
  memset(real_out, 0, sizeof(double) * n);
  memset(in, 0, sizeof(fftw_complex) * n);
  memset(x_multi_out, 0, sizeof(fftw_complex) * n2 * m);
  memset(sig_in, 0, sizeof(fftw_complex) * n * m);
  memset(g_multi_out, 0, sizeof(fftw_complex) * n2 * m);
  memset(out, 0, sizeof(fftw_complex) * n2);
  memset(x_fft, 0, sizeof(fftw_complex) * n);
  memset(conj, 0, sizeof(fftw_complex) * n);
  
  
  for (j = 0; j < m; ++j) {
    for (i = 0; i < n; ++i) {
      x_m_real_in[i + j * n] = signal(i,j);
    }
  }
  
  fftw_execute(x_m_real_p);
  
  // Compute fine levels for scale estimation
  j  = J - 1;
  nj = 1 << j;
  w1 = ceil(nj/3.);
  w2 = ceil(2. * nj/3.);
  w3 = n2 - pow(2.0, j - 3) - 1;
  p  = 1.0/pow((double)n, 1.0/2)/pow(2.0, j/2.0);
  
  // Input convolution matrix to compute FFT, also 
  for (j = 0; j < m; ++j) {
    g_m_real_in[j * n] = G(0,j);
    for (i = 1; i < w1; ++i) {
      g_m_real_in[i + j * n]     = G(i,j);
      g_m_real_in[n - i + j * n] = G(n - i,j);
      sig_in[i + j * n][0]       = 0;
      sig_in[i + j * n][1]       = 0;
      sig_in[n - i + j * n][0]   = 0;
      sig_in[n - i + j * n][0]   = 0;
    }
    for (i = w1; i < w2; ++i) {
      g_m_real_in[i + j * n]     = G(i,j);
      g_m_real_in[n - i + j * n] = G(n - i,j);
      xi                         = (double)i/nj;
      x                          = 3 * xi - 1;
      x                          = p * sin(M_PI_2 * MeyerPol(x, deg));
      c                          = cos(M_PI * xi) * x;
      s                          = sin(- M_PI * xi) * x;
      sig_in[i + j * n][0]       = c * x_multi_out[i + j * n2][0] - s * x_multi_out[i + j * n2][1];
      sig_in[i + j * n][1]       = s * x_multi_out[i + j * n2][0] + c * x_multi_out[i + j * n2][1];
      sig_in[n - i + j * n][0]   = sig_in[i + j * n][0];
      sig_in[n - i + j * n][1]   = - sig_in[i + j * n][1];
    }
    for (i = w2; i < w3; ++i) {
      g_m_real_in[i + j * n]     = G(i,j);
      g_m_real_in[n - i + j * n] = G(n - i,j);
      xi                         = (double)i/nj;
      x                          = 3 * xi/2 - 1;
      x                          = p * cos(M_PI_2 * MeyerPol(x, deg));
      c                          = cos(M_PI * xi) * x;
      s                          = sin(- M_PI * xi) * x;
      sig_in[i + j * n][0]       = c * x_multi_out[i + j * n2][0] - s * x_multi_out[i + j * n2][1];
      sig_in[i + j * n][1]       = s * x_multi_out[i + j * n2][0] + c * x_multi_out[i + j * n2][1];
      sig_in[n - i + j * n][0]   = sig_in[i + j * n][0];
      sig_in[n - i + j * n][1]   = - sig_in[i + j * n][1];
    }
    for (i = w3; i < n2; ++i) {
      g_m_real_in[i + j * n]     = G(i,j);
      g_m_real_in[n - i + j * n] = G(n - i,j);
      sig_in[i + j * n][0]       = - p * x_multi_out[i + j * n2][0];
      sig_in[i + j * n][1]       = - p * x_multi_out[i + j * n2][1];
      sig_in[n - i + j * n][0]   = sig_in[i + j * n][0];
      sig_in[n - i + j * n][1]   = - sig_in[i + j * n][1];
    }
  }
  
  fftw_execute(g_m_real_p);
  
  fftw_execute(sigma_back_p);
  
  // Estimate sigma and noise
  
  noise = est_noise(sig_real_out, n, m);
  
  sigma = est_sigma(noise);
  // Compute feasible levels
  if (resolution == "smooth") {
    if (blur == "direct") {
      channel = DirectChanInfo(m, n, sigma, alpha);
    } else {
      channel = SmoothChanInfo(m, n,  g_multi_out, sigma, alpha);
    }
  } else {
    if (resolution == "block") {
      channel = BoxCarChanInfo(m, n, g_multi_out, sigma, alpha, j0, deg);
    } else {
      Rf_warning("resolution input type not recognised, assumed regular smooth resolution.");
      channel = SmoothChanInfo(m, n,  g_multi_out, sigma, alpha);
    }
  }
  
  // Check j1 value is feasible
  if (j1 == NA_INTEGER) {
    j1 = as<int>(channel["finest"]);
  }
  if (j1 >= J) {
    Rf_warning("j1 too large, set to maximum feasible  j1 = log2(n) - 1");
    j1 = J - 1;
  }
  if (j1 < j0) {
    Rf_warning("Estimated finest resolution level, j1, too small (must be at least j0), Set j1 = j0");
    j1 = j0;
  }
  // Compute theoretical Eta if it is not specified
  if (R_IsNA(eta)) {
    eta = TheoreticalEta(alpha, resolution, m, n, g_multi_out, sigma);
  }

  if (thresh.size() == 0)
    thresh = MaxiThreshFFTW(n, sigma, g_multi_out, alpha, j0, j1, eta, deg);
  
  NumericMatrix proj(n, j1 - j0 + 2);
  
  // Compute the multichannel normalised Fourier coefficients
  mlwavedxfft(x_fft, m, n, x_multi_out, g_multi_out, sigma, alpha);
  
  // Coarse scale approximation (no thresholding)
  j = j0;
  
  memset(in, 0, sizeof (fftw_complex) * n);
  
  x  = 1.0/n;
  nj = 1 << j;
  w1 = ceil((double)nj/3);
  w2 = 2 * w1 + j % 2 - 1;
  p  = pow(2.0, 1 - j0) * M_PI;
  cx = x/pow(2.0, j/2.0);
  
  in[0][0] = x * x_fft[0][0];
  
  // 
  for (k = 0; k < nj; ++k)
    beta[k] = cx * x_fft[0][0];
  
  for (i = 1; i < w1; i++) {
    in[i][0]     = x * x_fft[i][0];
    in[i][1]     = x * x_fft[i][1];
    in[n - i][0] = x * x_fft[n - i][0];
    in[n - i][1] = x * x_fft[n - i][1];
    beta[0]     += cx * (x_fft[i][0] + x_fft[n - i][0]);
    for (k = 1; k < nj; ++k) {
      c        = cx * cos(p * i * k);
      s        = cx * sin(p * i * k);
      beta[k] += c  * x_fft[i][0] - s * x_fft[i][1];
      c        = cx * cos(p * (n - i) * k);
      s        = cx * sin(p * (n - i) * k);
      beta[k] += c  * x_fft[n - i][0] - s * x_fft[n - i][1];
    }
  }

  for (i = w1; i < w2; i++) {
    xi           = (double)i/nj;
    x            = 3 * xi - 1;
    xi           = pow(cos(M_PI_2 * MeyerPol(x, deg)), 2.0)/n;
    in[i][0]     = xi * x_fft[i][0];
    in[i][1]     = xi * x_fft[i][1];
    in[n - i][0] = xi * x_fft[n - i][0];
    in[n - i][1] = xi * x_fft[n - i][1];
    xi           = cos(M_PI_2 * MeyerPol(x, deg));
    x            = cx * xi;
    // Only real part survives for first beta coef
    beta[0] += x * x_fft[i][0] + x * x_fft[n - i][0];
    for (k = 1; k < nj; ++k) {
      c        = x * cos(p * i * k);
      s        = x * sin(p * i * k);
      beta[k] += c * x_fft[i][0] - s * x_fft[i][1];
      c        = x * cos(p * (n - i) * k);
      s        = x * sin(p * (n - i) * k);
      beta[k] += c * x_fft[n - i][0] - s * x_fft[n - i][1];
    }
  }

  fftw_execute(backward_p);
    
  // Start wavelet expansion and record coarse approximation
  for (i = 0; i < n; ++i) {
    final_out[i] = real_out[i];
    proj(i,0)    = real_out[i];
  }
  
  // Check if final J-1 level needed then add it.
  if (j1 == J - 1) {
    jmax = j1;
    j    = J - 1;
    
    memset(in, 0, sizeof (fftw_complex) * n);
    memset(conj, 0, sizeof (fftw_complex) * n);
    
    nj = 1 << j;
    w1 = ceil((double)nj/3.0);
    w2 = 2 * w1 + j % 2 - 1;
    w3 = n2 - pow(2.0, j - 3) - 1;
    p  = 1.0/n/pow(2.0,j/2.0);
    bp = 4 * M_PI / n;
    
    for (i = w1; i < w2; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi - 1;
      x              = p * sin(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(- M_PI * xi) * x;
      in[i][0]       = cx * x_fft[i][0] - sx * x_fft[i][1];
      conj[i][0]     = cx;
      in[i][1]       = sx * x_fft[i][0] + cx * x_fft[i][1];
      conj[i][1]     = - sx;
      in[n - i][0]   = cx * x_fft[n - i][0] + sx * x_fft[n - i][1];
      conj[n - i][0] = cx;
      in[n - i][1]   = - sx * x_fft[n - i][0] + cx * x_fft[n - i][1];
      conj[n - i][1] = sx;
      c              = x_fft[i][0] * cx + x_fft[i][1] * sx;
      s              = cx * x_fft[i][1] - sx * x_fft[i][0];
      bk             = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = x_fft[n - i][0] * cx - x_fft[n - i][1] * sx;
      s  = cx * x_fft[n - i][1] + sx * x_fft[n - i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
    
    for (i = w2; i < w3; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi/2 - 1;
      x              = p * cos(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(- M_PI * xi) * x;
      in[i][0]       = cx * x_fft[i][0] - sx * x_fft[i][1];
      conj[i][0]     = cx;
      in[i][1]       = sx * x_fft[i][0] + cx * x_fft[i][1];
      conj[i][1]     = - sx;
      in[n - i][0]   = cx * x_fft[n - i][0] + sx * x_fft[n - i][1];
      conj[n - i][0] = cx;
      in[n - i][1]   = - sx * x_fft[n - i][0] + cx * x_fft[n - i][1];
      conj[n - i][1] = sx;
      c              = x_fft[i][0] * cx + x_fft[i][1] * sx;
      s              = cx * x_fft[i][1] - sx * x_fft[i][0];
      bk             = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = x_fft[n - i][0] * cx - x_fft[n - i][1] * sx;
      s  = cx * x_fft[n - i][1] + sx * x_fft[n - i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
    
    n2 -= 1;
    for (i = w3; i < n2; ++i) {
      in[i][0]       = - p * x_fft[i][0];
      in[i][1]       = - p * x_fft[i][1];
      conj[i][0]     = - p;
      in[n - i][0]   = - p * x_fft[n - i][0];
      in[n - i][1]   = - p * x_fft[n - i][1];
      conj[n - i][0] = - p;
      c  = - p * x_fft[i][0];
      s  = - p * x_fft[i][1];  
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = - p * x_fft[n - i][0];
      s  = - p * x_fft[n - i][1];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
    in[n2][0]   = - p * x_fft[n2][0];
    in[n2][1]   = - p * x_fft[n2][1];
    conj[n2][0] = - p;
    // Add the last term to the coefficients
    bk = nj;
    for (k = 0; k < nj; ++k) {
      c         = - p * x_fft[n2][0];
      s         = - p * x_fft[n2][1];
      beta[bk] += c * cos(bp * n2 * k) - s * sin(bp * n2 * k);
      bk++;
    }
    n2 += 1;
    
    // Invert the product
    fftw_execute(backward_p);

    // Threshold the fft inversion
    ThresholdFFTW(real_out, real_in, n, thresh[j1 - j0], shrinkType);
    
    // fft forward
    fftw_execute(real_p);

    // Project back
    for (i = 0; i < n2; ++i) {
      in[i][0] = out[i][0] * conj[i][0] - out[i][1] * conj[i][1];
      in[i][1] = out[i][0] * conj[i][1] + out[i][1] * conj[i][0];
    }
    for (i = n2; i < n; ++i) {
      in[i][0] = out[n - i][0] * conj[i][0] + out[n - i][1] * conj[i][1];
      in[i][1] = out[n - i][0] * conj[i][1] - out[n - i][1] * conj[i][0];
    }

    fftw_execute(backward_p);
    
    k = j1 - j0 + 1;
    for (i = 0; i < n; ++i) {
      x             = real_out[i] * nj;
      final_out[i] += x;
      proj(i,k)     = x;
    }    
  } else {
    jmax = j1 + 1;
  }
  
  // Set factors first to simplify loop adjustments
  j   = j0 - 1;
  nj  = 1 << j;
  p   = 1.0/n/pow((double)nj, 1.0/2.0);
  bp  = M_PI * pow(2.0, 2 - j0);
  
  for (j = j0; j < jmax; ++j) {
    
    // Perform product of x_fft and psij_fft (Keep Conj(wavj for later))
    memset(in, 0, sizeof (fftw_complex) * n);
    memset(conj, 0, sizeof (fftw_complex) * n);
    
    p  /= pow(2.0, 0.5);
    bp /= 2;
    nj *= 2;
    i   = 1 << j;
    w1  = ceil(i/3.0);
    w2  = 2 * w1 + j % 2 - 1;
    w3  = w1 + nj;
  
    for (i = w1; i < w2; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi - 1;
      x              = p * sin(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(-M_PI * xi) * x;
      in[i][0]       = cx * x_fft[i][0] - sx * x_fft[i][1];
      conj[i][0]     = cx;
      in[i][1]       = sx * x_fft[i][0] + cx * x_fft[i][1];
      conj[i][1]     = - sx;
      in[n - i][0]   = cx * x_fft[n - i][0] + sx * x_fft[n - i][1];
      conj[n - i][0] = cx;
      in[n - i][1]   = - sx * x_fft[n - i][0] + cx * x_fft[n - i][1];
      conj[n - i][1] = sx;
      c              = x_fft[i][0] * cx + x_fft[i][1] * sx;
      s              = cx * x_fft[i][1] - sx * x_fft[i][0];
      bk             = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = x_fft[n - i][0] * cx - x_fft[n - i][1] * sx;
      s  = cx * x_fft[n - i][1] + sx * x_fft[n - i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
    
    for (i = w2; i < w3; ++i) {
      xi             = (double)i/nj;
      x              = 3 * xi/2 - 1;
      x              = p * cos(M_PI_2 * MeyerPol(x, deg));
      cx             = cos(M_PI * xi) * x;
      sx             = sin(- M_PI * xi) * x;
      in[i][0]       = cx * x_fft[i][0] - sx * x_fft[i][1];
      conj[i][0]     = cx;
      in[i][1]       = sx * x_fft[i][0] + cx * x_fft[i][1];
      conj[i][1]     = -sx;
      in[n - i][0]   = cx * x_fft[n - i][0] + sx * x_fft[n - i][1];
      conj[n - i][0] = cx;
      in[n - i][1]   = - sx * x_fft[n - i][0] + cx * x_fft[n - i][1];
      conj[n - i][1] = sx;
      c              = x_fft[i][0] * cx + x_fft[i][1] * sx;
      s              = cx * x_fft[i][1] - sx * x_fft[i][0];
      bk             = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * i * k) - s * sin(bp * i * k);
        bk++;
      }
      c  = x_fft[n - i][0] * cx - x_fft[n - i][1] * sx;
      s  = cx * x_fft[n - i][1] + sx * x_fft[n - i][0];
      bk = nj;
      for (k = 0; k < nj; ++k) {
        beta[bk] += c * cos(bp * (n - i) * k) - s * sin(bp * (n - i) * k);
        bk++;
      }
    }
    
    // Invert the product
    fftw_execute(backward_p);

    // Threshold the fft inversion
    ThresholdFFTW(real_out, real_in, n, thresh[j - j0], shrinkType);
    
    // fft forward
    fftw_execute(real_p);

    // Project back
    for (i = 0; i < n2; ++i) {
      in[i][0] = out[i][0] * conj[i][0] - out[i][1] * conj[i][1];
      in[i][1] = out[i][0] * conj[i][1] + out[i][1] * conj[i][0];
    }
    for (i = n2; i < n; ++i) {
      in[i][0] = out[n - i][0] * conj[i][0] + out[n - i][1] * conj[i][1];
      in[i][1] = out[n - i][0] * conj[i][1] - out[n - i][1] * conj[i][0];
    }
    fftw_execute(backward_p);
    k = j - j0 + 1;
    for (i = 0; i < n; ++i) {
      x             = real_out[i] * nj;
      final_out[i] += x;
      proj(i,k)     = x;
    }
  }
  
  // Determine Thresholded wavelet coefficients
  List shrunk_coefs            = ThresholdCoef(beta, j0, j1, thresh, shrinkType);
  NumericVector beta_shrink    = shrunk_coefs["coef"];
  NumericVector percent_shrunk = shrunk_coefs["percent"];
  NumericVector level_max      = shrunk_coefs["max"];
  
  fftw_destroy_plan(real_p);
  fftw_destroy_plan(x_m_real_p);
  fftw_destroy_plan(g_m_real_p);
  fftw_destroy_plan(backward_p);
  fftw_destroy_plan(sigma_back_p);

  fftw_free(in);
  fftw_free(x_m_real_in);
  fftw_free(g_m_real_in);
  fftw_free(real_in);
  fftw_free(out);
  fftw_free(real_out);
  fftw_free(x_multi_out);
  fftw_free(g_multi_out);
  fftw_free(x_fft);
  fftw_free(conj);
  fftw_free(sig_real_out);
  fftw_free(sig_in);
  
  List betaList = List::create(
    _["coef"] = beta,
    _["j0"]   = j0,
    _["deg"]  = deg
    );
  betaList.attr("class") = "waveletCoef";
  
  List shrunkBetaList = List::create(
    _["coef"] = beta_shrink,
    _["j0"]   = j0,
    _["deg"]  = deg
    );
  shrunkBetaList.attr("class") = "waveletCoef";
  
  List multiWaveD = List::create(
    _["signal"]            = signal,
    _["G"]                 = G,
    _["j0"]                = j0,
    _["j1"]                = j1,
    _["resolutionMethod"]  = resolution,
    _["blurDetected"]      = blur,
    _["alpha"]             = alpha,
    _["sigmaEst"]          = sigma,
    _["blurInfo"]          = channel,
    _["eta"]               = eta,
    _["thresholds"]        = thresh,
    _["estimate"]          = final_out,
    _["coef"]              = betaList,
    _["shrinkCoef"]        = shrunkBetaList,
    _["percent"]           = percent_shrunk,
    _["levelMax"]          = level_max,
    _["shrinkType"]        = shrinkType,
    _["projections"]       = proj,
    _["noise"]             = noise,
    _["degree"]            = deg
  );
  
  multiWaveD.attr("class") = "mWaveD";
  
  return multiWaveD;
}

//[[Rcpp::export]]
bool directDetect(NumericMatrix x){
  
  int i, j, m, n, n2;
  
  bool checkFlag = true;
  
  // FFTW specific definitions
  double       *x_m_real;
  fftw_complex *x_m_out;
  fftw_plan     x_m_real_p;
  
  n  = x.nrow();
  m  = x.ncol();  
  n2 = n/2 + 1;

  x_m_real   = (double*)fftw_malloc(sizeof(double) * n * m);
  x_m_out    = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n2 * m);
  x_m_real_p = fftw_plan_many_dft_r2c(1, &n, m, x_m_real, NULL, 1, n, x_m_out , NULL, 1, n2, FFTW_ESTIMATE);
  
  for ( j = 0; j < m; ++j ){
    for ( i = 0; i < n; ++i ){
      x_m_real[i + j * n] = x(i,j);
    }
  }
  
  fftw_execute(x_m_real_p);
  
  for (j = 0; j < m; ++j) {
    for (i = 0; i < n2; ++i) {
      if ( fabs(x_m_out[i + j * n2][0] - 1.0) > DBL_EPSILON || fabs(x_m_out[i + j * n2][1]) > DBL_EPSILON) {
        checkFlag = false;
        //  break out of both loops
        j = m;
        break;
      }
    }
  }

  fftw_free(x_m_out);
  fftw_free(x_m_real);
  
  fftw_destroy_plan(x_m_real_p);
  
  return checkFlag;
}

//[[Rcpp::export]]
bool boxcarDetect(NumericMatrix x){
  
  bool result = true;
  
  NumericMatrix vals(2, x.ncol());
  
  for(int i = 0; i< x.ncol(); i++){
    NumericVector test = unique(x(_, i));
    if (test.size() != 2) {
      result = false;
      break;
    } else {
      vals(_, i) = test;
    }
  }
  
  if (is_true(all(vals(0, _) == 0)) && is_true(all(vals(1, _) > 0))) {
    result = true;
  } else {
    result = false;
  }
  
  return result;
}

//[[Rcpp::export]]
NumericMatrix directBlur(int n, int m) {
  
  NumericMatrix result(n, m);
  
  for (int i = 0; i < m; ++i) {
    result(0, i) = 1;
  }
  
  return result;
}

//[[Rcpp::export]]
double theoreticalEta(NumericVector alpha, String blur, ComplexMatrix g_fft, NumericVector sigma){
  
  double eta;

  int i, j, n2, m, n;
  int best = 1;
  m = g_fft.ncol();
  n = g_fft.nrow();
  
  // Check dimensions agree
  if (!( m == alpha.size() && m == sigma.size()))
    stop("Dimension mismatch; number of rows of g_fft, and lengths of alpha and sigma should agree");
  
  double tmp, nsqrt;
  
  NumericVector threshfft(m);
  NumericVector level(m);
  NumericVector freq(m, -1.0);
  n2    = n/2 + 1;
  nsqrt = pow((double)n, 1.0/2);
  
  
  // Compute channel level thresholds
  for (j = 0; j < m; ++j) {
    tmp   = log(sigma[j]);
    threshfft[j] = tmp - log(pow((double)n, alpha[j]/2.0)) + 0.5 * log(fabs(log(nsqrt) - tmp));
  }
  // Compute the normalised Fourier decay levels.
  for (j = 0; j < m; ++j) {
    // Skip first level since Inf always > threshfft 
    for (i = 1; i < n2; ++i) {
      tmp  = log(sqrt((g_fft(i,j).r * g_fft(i,j).r) + (g_fft(i,j).i * g_fft(i,j).i))) - 0.5 * alpha[j] * log((double)i);
      if (tmp < threshfft[j]) {
        freq[j]  = i + 1;
        level[j] = floor(log2((double)i + 1)) - 1;
        break;
      }
    }
    // Threshold never met (direct convolution case)
    if (freq[j] == -1.0) {
      freq[j]  = n2 - 1;
      level[j] = log2((double)freq[j]);
    }
    
  }
  
  best = 1;
  int current = freq[0];
  for (j = 1; j < m; ++j) {
    if (freq[j] > current) {
      current = freq[j];
      best    = j + 1;
    } 
  }
  eta = 4 * sqrt(alpha[best - 1]);
  
  return eta;
}

