source("functions/sub_functions/compute_density_.R")
require(crayon)
#'
#'
#'
generate_densities_ <- function(ow.ppp, bw.list, verbose=T) {
  if(verbose) cat(crayon::bgGreen("\n\nFUNCT : generate_densities_ ...\n") )
  est.bw <- c(bw.list[[1]], bw.list[[2]])
  range.bw <- seq(from=est.bw[1], to=est.bw[2], length.out = 7)
  inc <- range.bw[2]- range.bw[1]; f1 <- range.bw[1]; l1 <- range.bw[7]
  range.bw <- c(f1-2*inc, f1-inc, range.bw, l1+inc, l1+2*inc)
  
  for(i in seq_along(range.bw)) { # loop over 
    bw.tmp <- range.bw[i]
    # compute a set of densities using the computed bandwidth
      if(verbose) cat( paste0("\n\n\ti = ", i, "\tbw val = ", bw.tmp) )
      dens.tmp <- compute_density_(
        ppp.in = ow.ppp, 
        bw.in = bw.tmp,
        verbose = T, 
        DBG = T
      )
      if(verbose) cat(paste0("\n\t\tdone computing, saving..."))
      saveRDS(dens.tmp, file=paste0("clean_data/density/dens_bw_", round(bw.tmp, 0), ".Rds"))
      rm(dens.tmp)
      gc()
    }
    
}
