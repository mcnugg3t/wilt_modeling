require(tidyverse)
require(crayon)
#' takes 2 column df - first one is counts, second is value. simply unfolds the counts
#'
#'
unfold_ <- function(df.2col, verbose=T) {
  if(verbose) {
    cat("\t")
    cat(crayon::bgGreen("FUNCT : unfold_"))
  }
  
  count.v <- df.2col[[1]]
  val.v <- df.2col[[2]]
  
  v.to.return <- c()
  
  n.iter <- length(count.v)
  inc.tmp <- n.iter %/% 10
  cat("\tunfolding...")
  for(i in seq_along(count.v)) {
    if(verbose & (i %% inc.tmp == 0)) cat(paste0(i %/% inc.tmp, "..")) 
    v.to.return <- c(
                    v.to.return, 
                    rep(val.v[i], count.v[i]))
  }
  
  return(v.to.return)
}