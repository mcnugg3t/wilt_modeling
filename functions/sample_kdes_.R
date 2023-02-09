require(spatstat)
require(tidyverse)
require(assertthat)
require(terra)
require(crayon)
#'
#'
#'
sample_kdes_ <- function(covar.rast, n.sim, verbose=T, DBG=F) {
  if(verbose) cat( crayon::bgGreen("\n\nCALL FUNCT: sample_kdes_") )
  
  # identify density files
  folds.tmp <- list.files("clean_data/density/")
  for(i in seq_along(folds.tmp)) {
    fld <- folds.tmp[i]
    if(verbose) cat(paste0("\n\n\nfolder : ", fld)) # debug print
    pth.tmp <- paste0("clean_data/density/", fld, "/")
    # identify kde files
    kdes.tmp <- list.files( pth.tmp )
    if(DBG) cat(paste0("\n\tfound : ", paste0(kdes.tmp, collapse="  ,  "))) # debug print
    # loop over them, for each...
    for(j in seq_along(kdes.tmp)) {
      # construct paths
      fl.pth <- paste0(pth.tmp, kdes.tmp[j])
      save.pth <- paste0("clean_data/sample/", fld, "/sample_",kdes.tmp[j])
      if(verbose) cat(paste0("\n\n\t\treading : ", fl.pth)) # debug print
      # read file -> df
      kde.df <- readRDS( fl.pth ) |> 
        as.data.frame()
      if(DBG) cat("\n\t\tscaling...")
      kde.df <- kde.df |> 
        mutate(val_scale = value/sum(value, na.rm=T))
      if(DBG) cat("\n\t\tchecking...")
      assert_that(length(kde.df[,3][kde.df[,3]<0]) == 0) # no negative
      assert_that((1-sum(kde.df$val_scale))^2 < 0.1e-6) # sums to 1
      # sample indices n.sim times
      if(DBG) cat("\n\t\tsampling indices...")
      ind.sample <- sample(
        1:nrow(kde.df), 
        size=n.sim, 
        replace=T, 
        prob=kde.df$val_scale)
      if(DBG) cat("\n\t\tgetting unique...")
      ind.unique <- ind.sample |> 
        as_tibble() |>
        rename(index = value) |> 
        group_by(index) |> 
        summarise(count = n()) |> 
        arrange(desc(count)) 
      ind.unique <- ind.unique |> 
        add_column(x = kde.df$x[ind.unique$index],
                   y = kde.df$y[ind.unique$index])
      # extract
      ind.extr <- ind.unique |> 
        select(x, y) |> 
        as.matrix()
      if(DBG) cat("\n\t\textracting...")
      extr.res <- terra::extract(covar.rast, ind.extr)
      if(DBG) cat("\n\t\tmerging...")
      extr.res.full <- cbind(ind.unique, extr.res)
      if(DBG) cat("\n\t\tsaving...")
      # save tibble under clean_data/sample/
      if(DBG) cat(paste0("\n\tsave path = ", save.pth) )
      saveRDS(extr.res.full, file=save.pth)
    }
  }
  
  
}