##' Simulation-based circadian power calculation assuming Normal errors in the cosiner model.
##'
##' Calculate power of simulated data with independent Normal errors in the cosiner model.
##' @title Simulation-based circadian power calculation assuming Normal errors in the cosiner model.
##' @param B_MC Numer of Monte-Carlo simulations. Default is 10000.
##' @param n Sample size.
##' @param A Amplitude of the sinusoidal curve: \eqn{A * sin(2\pi/period * (phi + cts)) + M}.
##' @param sigma Standard deviation \eqn{\sigma} of each gene in the Normal error term \eqn{N(0,\Sigma)},
##' where \eqn{\Sigma} is a symmetric matrix of size GmoduleSize with diagonal elements being \eqn{\sigma^2} and off-diagonal
##' elements being \eqn{\sigma^2\rho}.
##' @param phi Phase shift of the sinusoidal curve. Default is 0.
##' @param period Period of the sinusoidal curve. Default is 24.
##' @param rho Correlation coefficient in the Normal error term \eqn{N(0,\Sigma)},
##' where \eqn{\Sigma} is a symmetric matrix of size GmoduleSize with diagonal elements being \eqn{\sigma^2} and off-diagonal
##' elements being \eqn{\sigma^2\rho}.
##' @param GmoduleSize Correlated genes are simulated in modules of size GmoduleSize. Only used when \eqn{\rho} is not 0.
##' @param cts Circadian times of samples.
##' @param alpha Type I error control. Default is 0.05.
##' @param M MESOR (Midline Statistic Of Rhythm, a rhythm-adjusted mean) of the sinusoidal curve. Default is 10.
##' @return A vector of empirical type I error and power.
##' @author Wei Zong, Zhiguang Huo
##' @export
##' @examples
##' B_MC = 10000
##' n = 100
##' cts = seq(0,24,length.out = n+1)[-1]
##' A = 1
##' sigma = 1
##' phi = 0
##' M = 10
##' simuPower_norm(B_MC, n, A, sigma, phi=0, period = 24, cts=cts, alpha = 0.05)

simuPower_norm <- function(B_MC=10000, n, A, sigma, phi = 0, period = 24, M = 10,
                           rho = 0, GmoduleSize = 50,
                           cts=NULL, alpha = 0.05){
  stopifnot(length(cts) == n)

  if(rho != 0 & B_MC%%GmoduleSize != 0){
    stop("The correlated gene module size GmoduleSize has to be a divisor of the umer of Monte-Carlo simulations B_MC.")
  }
  Gmodule = B_MC / GmoduleSize

  ##Type I error
  MC.data0.ls = lapply(1:Gmodule, function(b){
    set.seed(b)
    s1 = matrix(rho, GmoduleSize, GmoduleSize)
    s2 = diag(rep(1 - rho, GmoduleSize))
    s = (s1 + s2)*sigma^2
    error = MASS:: mvrnorm(n = n, mu=rep(0,GmoduleSize), s)
    curve = matrix(rep(0 + M,GmoduleSize),nrow = GmoduleSize, ncol = n, byrow = T)
    aMC.data = curve + t(error)
    return(aMC.data)
  })
  MC.data0 = do.call(rbind, MC.data0.ls)
  pvalues.LR0 = apply(MC.data0, 1, function(yy) LR_rhythmicity(cts,yy)$pvalue)
  F_MC_typeIerror = mean(pvalues.LR0<alpha)

  #Power
  MC.data.ls = lapply(1:Gmodule, function(b){
    set.seed(b)
    s1 = matrix(rho, GmoduleSize, GmoduleSize)
    s2 = diag(rep(1 - rho, GmoduleSize))
    s = (s1 + s2)*sigma^2
    error = MASS:: mvrnorm(n = n, mu=rep(0,GmoduleSize), s)
    curve = matrix(rep(A*sin(2*pi/period*(cts + phi)) + M,GmoduleSize),nrow = GmoduleSize, ncol = n, byrow = T)
    aMC.data = curve + t(error)
    return(aMC.data)
  })
  MC.data = do.call(rbind, MC.data.ls)
  pvalues.LR = apply(MC.data, 1, function(yy) LR_rhythmicity(cts,yy)$pvalue)
  F_MC_power = mean(pvalues.LR < alpha)

  return(c(F_MC_typeIerror=F_MC_typeIerror, F_MC_power=F_MC_power))
}

#internal functions from 'diffCircadian'
LR_rhythmicity <- function(tt,yy,period=24,method="LR",FN=TRUE){
  if(method=="Wald"){
    WaldTest(tt=tt, yy=yy, period = period, type=FN)
  }
  else if(method=="LR"){
    LRTest(tt=tt, yy=yy, period = period, type=FN)
  }
  else(("Please check your input! Method only supports 'Wald','LR','F' or 'Permutation'."))
}

LRTest <- function(tt,yy, period = 24,type=TRUE){
  fitCurveOut <- fitSinCurve(tt,yy,period=period)
  n <- length(yy)
  rss <- fitCurveOut$rss
  tss <- fitCurveOut$tss

  amp <- fitCurveOut$amp
  phase <- fitCurveOut$phase
  offset <- fitCurveOut$offset

  sigma02 <- 1/(n)*sum((yy-mean(yy))^2)
  sigmaA2 <- 1/(n)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)

  l0 <- -n/2*log(2*pi*sigma02)-1/(2*sigma02)*sum((yy-mean(yy))^2)
  l1 <- -n/2*log(2*pi*sigmaA2)-1/(2*sigmaA2)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)

  dfdiff <- (n-1)-(n-3)
  LR_stat <- -2*(l0-l1)

  if(type==FALSE){
    pvalue <- pchisq(LR_stat,dfdiff,lower.tail = F)
  }
  else if(type==TRUE){
    r <- 2
    k <- 3
    LR_stat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(LR_stat,df1 = r, df2 = n-k, lower.tail = F)
  }
  R2 <- 1-rss/tss
  res <- list(
    amp = amp,
    phase = phase,
    peakTime = (6 - phase) %% period,
    offset = offset,
    sigma02=sigma02, sigmaA2=sigmaA2,
    l0=l0,
    l1=l1,
    stat=LR_stat,
    pvalue=pvalue,R2=R2)
  return(res)
}

fitSinCurve <- function(tt, yy, period = 24, parStart = list(amp=3,phase=0, offset=0)){

  getPred <- function(parS, tt) {
    parS$amp * sin(2*pi/period * (tt + parS$phase)) + parS$offset
  }

  residFun <- function(p, yy, tt) yy - getPred(p,tt)

  nls.out <-minpack.lm::nls.lm(par=parStart, fn = residFun, yy = yy,	tt = tt)

  apar <- nls.out$par

  amp0 <- apar$amp
  asign <- sign(amp0)
  ## restrict amp > 0
  amp <- amp0 * asign

  phase0 <- apar$phase
  #phase <- (round(apar$phase) + ifelse(asign==1,0,12)) %% period
  phase <- (phase0 + ifelse(asign==1,0,period/2)) %% period
  offset <- apar$offset

  peak <- (period/2 * sign(amp0) - period/4 - phase) %%period
  if(peak > period/4*3) peak = peak - period

  A <- amp0 * cos(2*pi/period * phase0)
  B <- amp0 * sin(2*pi/period * phase0)

  rss <- sum(nls.out$fvec^2)
  tss <- sum((yy - mean(yy))^2)
  R2 <- 1 - rss/tss

  if(F){
    amp <- apar$amp
    phase <- apar$phase
    offset <- apar$offset
  }

  res <- list(amp=amp, phase=phase, offset=offset, peak=peak, A=A, B=B, tss=tss, rss=rss, R2=R2)
  res
}


