require(tidyverse)
require(spatstat)
require(crayon)
#'
#'
#'
simulation_envelopes_L_ <- function(dens.files, ow.ppp, verbose=T, DBG=F) {
  complete.files <- list.files("clean_data/envelope/")
  if(verbose) cat(crayon::bgGreen("\nFUNCT : simulation_envelopes_L_"))
  for(j in seq_along(dens.files)) {
    # setup
    fn.tmp <- dens.files[j]
    fp.tmp <- paste0("clean_data/sample/", fn.tmp)
    cat(paste0("\n\tj = ", j, "\tfl.tmp = ", fn.tmp))
    bw.tmp <- fn.tmp |> 
      str_remove("ppp_list_bw_") |> 
      str_remove(".Rds")
    save.sml <- paste0("env_", bw.tmp, ".Rds")
    if(save.sml %in% complete.files) next
    save.tmp <- paste0("clean_data/envelope/", save.sml)
    save.img <- paste0("clean_data/envelope/img_", bw.tmp, ".jpg")
    if(DBG) cat(paste0("\n\t\tsave.tmp = ", save.tmp, "\n\t\tsave.img = ", save.img, "\n\t"))
    #
    if(verbose) {
      cat(crayon::bgCyan("loading ppp list..."))
      cat("\n\t")
      cat(crayon::bgYellow("calc envelope...\n"))
    }
    ppp.list.tmp <- readRDS(fp.tmp)
    envelope.tmp <- envelope(
      ow.ppp, 
      fun=Lest, 
      funargs=list(rmax=700,correction="border"),
      simulate= ppp.list.tmp # list of point patterns
    )
    
    cat("\n\t")
    cat(crayon::bgWhite("saving plot..."))
    jpeg(file=save.img)
    envelope.tmp |> plot(main=paste0("bandwidth: ", bw.tmp))
    dev.off()
    cat("\n\t")
    cat(crayon::bgWhite("saving envelope..."))
    saveRDS(envelope.tmp, file=save.tmp)
  }
}