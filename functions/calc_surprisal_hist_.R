require(assertthat)
require(crayon)
#' returns vector of suprisal values, takes hist.in (histogram), sample.v.in (vector, numeric) as inputs
#'
#'
calc_surprisal_hist_ <- function(hist.in, sample.v.in) {
  assert_that(class(hist.in) == "histogram")
  assert_that(class(sample.v.in) == "numeric")
  mid.v <- hist.in[["mids"]]
  surprise.v <- sample.v.in |> 
    vapply(FUN = function(x) {
      dist.v <- sqrt( (x-mid.v)^2 )
      min.dist.ind <- which(dist.v == min(dist.v))
      len.tmp <- length(min.dist.ind)
      if(len.tmp > 1) {
        cat("\n\t")
        cat(crayon::bgMagenta("resolving distance tie..."))
        min.dist.ind <- sample(x=min.dist.ind, size=1, prob = rep(1/len.tmp, len.tmp))
      }
      prob.ret <- hist.in$density[min.dist.ind] / sum(hist.in$density)
      return(-1*log(prob.ret))
    }, FUN.VALUE = numeric(1))
  return(surprise.v)
}