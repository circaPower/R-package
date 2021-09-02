##' Circadian times sampled from the density estimated from observed sampling times in real data.
##'
##' Generate circadian times at different sample sizes from kernel density estimated
##' from observed sampling times in real data.
##' @title Circadian Times Sampling from the Kernal Density estimated from Real Data
##' @param observed_cts Observed circadian times in real data. Lie between cts_lb and cts_ub.
##' @param cts_lb Lower bound of circadian times.
##' @param cts_ub Lower bound ofcircadian times.
##' @param n Number of time of death to generate.
##' @return A circadian times vector of length n
##' @author Wei Zong, Zhiguang Huo
##' @export
##' @examples
##' data(cts_Chen)
##' cts_kernel_generator(observed_cts=cts_Chen, cts_lb = 0, cts_ub = 24, n = 100)

cts_kernel_generator = function(observed_cts, cts_lb = 0, cts_ub = 24, n){
  if(min(observed_cts)<cts_lb | max(observed_cts)>cts_ub){
    stop("observed_cts should lie between cts_lb and cts_ub")
  }
  period = abs(cts_ub - cts_lb)
  dens = density(observed_cts, n=1000, from = cts_lb, to = cts_ub)
  sample_cts = sample(dens$x, n, prob = dens$y)
  sample_cts2 = ifelse(sample_cts>cts_ub,sample_cts-period,sample_cts)
  sample_cts3 = ifelse(sample_cts2<cts_lb,sample_cts2+period,sample_cts2)
  return(sample_cts3)
}
