##' Simulation-based circadian power calculation assuming independent Normal errors with outliers in the cosiner model.
##'
##' Calculate power of simulated data assuming independent Normal errors with outliers in the cosiner model.
##' @title Simulation-based circadian power calculation assuming independent Normal errors with outliers in the cosiner model.
##' @param B_MC Numer of Monte-Carlo simulations. Default is 10000.
##' @param n Sample size.
##' @param A Amplitude of the sinusoidal curve: \eqn{A * sin(2\pi/period * (phi + cts)) + M}.
##' @param sigma \eqn{\sigma} in the independent Normal error term \eqn{N(0,\sigma)}.
##' @param outlier.perc The percentage of outlying samples. Default is 0.05.
##' @param phi Phase shift of the sinusoidal curve. Default is 0.
##' @param period Period of the sinusoidal curve. Default is 24.
##' @param M MESOR (Midline Statistic Of Rhythm, a rhythm-adjusted mean) of the sinusoidal curve. Default is 10.
##' @param cts Circadian times of samples.
##' @param alpha Type I error control. Default is 0.05.
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
##' outlier.perc = 0.05
##' simuPower_outlier(B_MC, n, A, sigma, outlier.perc, phi, period = 24, M=10, cts,
##' alpha = 0.05)
simuPower_outlier = function(B_MC=10000, n, A, sigma, outlier.perc=0.05, phi=0, period = 24, M=10, cts=NULL, alpha = 0.05){
  stopifnot(length(cts) == n)

  out.cts = sample(cts,round(outlier.perc*n))
  sin.cts = setdiff(cts, out.cts)
  cts_sort = c(sin.cts,out.cts)

  ##TypeI error
  MC.data0 = NULL
  for(b in 1:B_MC){
    set.seed(b)
    error = rnorm(length(sin.cts),0,sigma)
    yy.sin = 0 + M + error

    lb = min(c(yy.sin,M-A))
    ub = max(c(yy.sin,M+A))
    yy.out = runif(length(out.cts),lb,ub)

    yy = c(yy.sin, yy.out)
    MC.data0 = rbind(MC.data0, yy)
  }
  pvalues.LR0 = apply(MC.data0, 1, function(yy) LR_rhythmicity(cts_sort,yy)$pvalue)
  F_MC_typeIerror = mean(pvalues.LR0<alpha)

  ##Power
  MC.data = NULL
  for(b in 1:B_MC){
    set.seed(b)
    error = rnorm(length(sin.cts),0,sigma)
    yy.sin = A*sin(2*pi/period*(sin.cts + phi)) + M + error

    lb = min(c(yy.sin,M-A))
    ub = max(c(yy.sin,M+A))
    yy.out = runif(length(out.cts),lb,ub)

    yy = c(yy.sin, yy.out)
    MC.data = rbind(MC.data, yy)
  }
  pvalues.LR = apply(MC.data, 1, function(yy) LR_rhythmicity(cts_sort,yy)$pvalue)
  F_MC_power = mean(pvalues.LR < alpha)

  return(c(F_MC_typeIerror=F_MC_typeIerror, F_MC_power=F_MC_power))
}
