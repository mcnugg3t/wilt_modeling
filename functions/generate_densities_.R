source("functions/sub_functions/compute_density_.R")
require(crayon)
#'
#'
#'
generate_densities_ <- function(ow.ppp, bw.list, adj.v, verbose=T) {
  if(verbose) cat(crayon::bgGreen("\n\nFUNCT : generate_densities_ ...\n") )
  # for each item in kern.list
  for(i in seq_along(bw.list)) {
    # extract estimation method and vector of adjustments from list
    kern.tmp <- names(bw.list)[i]
    bw.tmp <- bw.list[[i]]
    if(verbose) cat(paste0("\n\tkern.tmp = ", kern.tmp, "\tkernel bw.tmp = ", bw.tmp))
    # compute a set of densities using the computed bandwidth
    for(j in seq_along(adj.v)) {
      if(verbose) cat(paste0("\n\n\t\tj = ", j, "\tadjust val = ", adjust.v[j],"\tj = ", j))
      dens.tmp <- compute_density_(
        ppp.in = ow.ppp, 
        bw.in = bw.tmp, 
        adj.num = adjust.v[j],
        verbose = T, 
        DBG = T
      )
      if(verbose) cat(paste0("\n\t\tdone computing, saving..."))
      bw.adj <- bw.tmp*adj.v[j]
      saveRDS(dens.tmp, file=paste0("clean_data/density/", kern.tmp, "/dens_bw_", bw.adj, ".Rds"))
      rm(dens.tmp)
      gc()
    }
    
  }
}