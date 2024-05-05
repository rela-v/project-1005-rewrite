## Using minpack.lm
parStartVal <- list(A=3, phase=0, offset=0)
fitSinCurve <- function(xx, observed, parStart=parStartVal) {
  # print("fitSinCurve.R: Setting up for fitting.")
  # function to get predictions
  getPred <- function(parS, xx) {
    # Amplitude * sine((1/24 Hz frequency * pi )*(xx + phase)) + offset
    parS$A * sin((1/12*pi)*(xx+parS$phase))+parS$offset
  }
  # function to calculate errors
  residFun <- function(p, observed, xx) {
    observed - getPred(p, xx)
  }
  # print("fitSinCurve.R: Fitting sinusoid using LM algorithm.")
  # output of the LM algorithm's sinusoidal fit
  nls.out <- nls.lm(par=parStartVal, fn = residFun, observed = observed, xx = xx)  
  # get parameters of the fitted function
  apar <- nls.out$par
  
  # assign amplitude
  A0 <- apar$A
  # assign the sign (+/-) of the amplitude to 'asign'
  asign <- sign(A0)
  
  # Restrict A > 0
  A <- A0 * asign
  
  # Calculate rounded phase
  phase <- (apar$phase + ifelse(asign==1,0,12)) %% 24 
  
  offset <- apar$offset
  peak <- (12 * sign(A0) - 6 - phase) %%24
  if(peak > 18) peak = peak - 24
  # print("fitSinCurve.R: Calculating statistics.")
  SSE <- sum(nls.out$fvec^2)
  SST <- sum((observed - mean(observed))^2)
  R2 <- 1 - SSE/SST
  res <- list(A=A, phase=phase, offset=offset, peak=peak, R2=R2)
  # print("fitSinCurve.R: Protocol complete.")
  return(list(res, summary(nls.out)))
}
