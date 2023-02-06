#'
#'
#'
detect_remove_vars_ <- function(var.names, patterns, verbose) {
  if(verbose) cat("\n\n\tDETECT VARS TO REMOVE...")
  to.return <- c()
  # loop over patterns
  for(i in seq_along(patterns)) {
    pat.tmp <- patterns[i]
    var.ind.tmp <- str_detect(var.names, pat.tmp)
    to.return <- c(to.return, var.names[var.ind.tmp])
  }
  
  return(to.return)
}