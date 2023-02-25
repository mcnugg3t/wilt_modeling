source("functions/sub_functions/compute_density_.R")
require(crayon)
#'
#'
#'
generate_densities_ <- function(ow.ppp, bw.v, grd.int, verbose=T) {
  if(verbose) cat(crayon::bgGreen("\n\nFUNCT : generate_densities_ ...\n") )
  if(verbose) cat(paste0("\n\tbw.v : " , paste0(bw.v, collapse= " , ")))
  density.computed = list.files("clean_data/density/")
  for(i in seq_along(bw.v)) { # loop over 
    bw.tmp <- bw.v[i] # extract bw val
    
    if(grd.int == 30) {
      save.name.sml = paste0("g30_dens_bw_", round(bw.tmp, 0), ".Rds")
    } else {
      save.name.sml = paste0("dens_bw_", round(bw.tmp, 0), ".Rds")
    }
    
    # if bw already computed, next
    if(save.name.sml %in% density.computed) next
    save.name.tmp = paste0("clean_data/density/", save.name.sml)
    
    # compute a set of densities using the computed bandwidth
    if(verbose) cat( paste0("\n\n\ti = ", i, "\tbw val = ", bw.tmp) )
      dens.tmp <- compute_density_(
        ppp.in = ow.ppp, 
        bw.in = bw.tmp,
        verbose = T, 
        DBG = T
      )
    if(verbose) cat(paste0("\n\t\tdone computing, saving..."))
    saveRDS(dens.tmp, file= save.name.tmp)
    rm(dens.tmp)
    gc()
  }
}
