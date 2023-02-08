source("functions/sub_functions/compute_density_.R")
#'
#'
#'
est_bandwidth_ <- function(kern.v, ow.ppp, verbose=T) {
  if(verbose) cat("\n\nCALL FUNCT: est_bandwidth_ ...")
  return.list <- list()
  # for each item in kern.list
  for(i in seq_along(kern.v)) {
    # get estimation method
    kern.tmp <- kern.v[i]
    if(verbose) {
      cat( crayon::bgRed("\n\tkernel estimation method = ", kern.tmp) )
      cat( paste0("\t\ti = ", i)  )
    }
    
    ### est bandwidth
    if(verbose) cat("\n\t\test bandwidth...")
    t1 <- Sys.time()
    if(kern.tmp == "bw.ppl") {
      bw.tmp <- bw.ppl(ow.ppp, correction="translation", srange=c(1, 400), ns=32, shortcut=F)
    } else if (kern.tmp == "bw.diggle") {
      bw.tmp <- bw.diggle(ow.ppp, correction="translation", hmax=400, nr=512)
    } else {
      cat(crayon::bgMagenta("\n\nshouldn't happen\n\n"))
    }
    if(verbose) cat(paste0(" = ", bw.tmp, "\t\ttime : ", difftime(Sys.time(), t1, units="secs")))
    #
    # compute a set of densities using the computed bandwidth
    
    return.list[[kern.tmp]] <- bw.tmp
  }
  return(return.list)
}