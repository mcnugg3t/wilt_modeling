require(crayon)
#'
#'
#'
create_interact_terms_ <- function(cont.var.v, verbose=T, DBG=F) {
  if(verbose) cat(crayon::bgGreen("\nFUNCT : create_interact_terms_"))
  
  interact.v <- c()
  n.end <- length(cont.var.v)
  
  for(i in seq_along(cont.var.v)) {
    if (i == n.end) next
    term1 <- cont.var.v[i]
    if(DBG) cat("\n\ti = ", i, " ... ", term1, "\n\t\tlooping ", (i+1), " to ", n.end)
    for(j in (i+1):n.end) {
      term2 <- cont.var.v[j]
      if(DBG) cat("\n\t\t\tj = ", j, " ... ", term2)
      interact.v <- c(interact.v, paste0(term1, " x ", term2))
    }
  }
  
  return(interact.v)
}