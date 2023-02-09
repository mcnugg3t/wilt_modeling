require(spatstat)
require(tidyverse)
require(assertthat)
require(terra)
require(crayon)
source("functions/unfold_.R")
#' takes covariate raster and nsim, loops over all KDEs and computes sampling dist -> saves
#' 
#'
sample_kdes_ <- function(covar.rast, join.dat, n.sim, n.breaks=100, verbose=T, DBG=F) {
  # setup
  if(verbose) cat( crayon::bgGreen("\n\nFUNCT : sample_kdes_") )
  folds.tmp <- list.files("clean_data/density/") # identify density files
  
  ##
  has.saved.hist <- F
  has.saved.tab <- F
  ##
  
  # small loop - over 2 folders
  for(i in seq_along(folds.tmp)) {
    fld <- folds.tmp[i]
    if(verbose) cat(paste0("\n\n\nfolder : ", fld)) # debug print
    pth.tmp <- paste0("clean_data/density/", fld, "/")
    kdes.tmp <- list.files( pth.tmp ) # identify kde files
    if(DBG) cat(paste0("\n\tfound : ", paste0(kdes.tmp, collapse="  ,  "))) # debug print
    
    ##
    ## BIG LOOP over KDEs
    for(j in seq_along(kdes.tmp)) {
      cat(paste0("\n\n\tj = ", j))
      # construct paths
      fn.tmp <- kdes.tmp[j]
      fl.pth <- paste0(pth.tmp, fn.tmp)
      bw.tmp <- fn.tmp |> 
        str_remove("dens_bw_") |> 
        str_remove(".Rds") |> 
        as.numeric() |> 
        round(2)
      
      save.pth <- paste0("clean_data/sample_dist/", fld, "/samp_dist_bw_", bw.tmp, ".Rds" )
      
      if(verbose) cat(paste0("\n\t\treading : ", fl.pth, "\n\t\tsave path : ", save.pth)) # debug print
      # read file -> df
      kde.df <- readRDS( fl.pth ) |> 
        as.data.frame()
      if(DBG) cat("\n\t\tscaling...")
      cutoff.tmp <- (0 - kde.df$value)^2 |> quantile(probs=c(0.6))
      if(DBG) cat(paste0("\n\t\t60th percentile cutoff = ", cutoff.tmp))
      kde.df <- kde.df |> 
        mutate(value = if_else( ((0-value)^2) < cutoff.tmp, 0, value)) |> # if distance between 0 and value is less than 60th percentile, set to 0 
        mutate(val_scale = value/sum(value, na.rm=T))
      if(DBG) cat("\n\t\tchecking...")
      val.scale <- kde.df |> select(val_scale)
      assert_that(length(val.scale[val.scale<0]) == 0) # no negative values
      assert_that((1-sum(val.scale))^2 < 0.1e-6) # sums to 1, or very close
      
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
      # add x and y coordinates
      ind.unique <- ind.unique |> 
        add_column(x = kde.df$x[ind.unique$index],
                   y = kde.df$y[ind.unique$index])
      # extract covariate values from only unique x and y coordinates
      ind.extr <- ind.unique |> 
        select(x, y) |> 
        as.matrix()
      if(DBG) cat("\n\t\textracting...")
      extr.res <- terra::extract(covar.rast, ind.extr)
      if(DBG) cat("\n\t\tmerging...")
      extr.res.full <- cbind(ind.unique, extr.res) |> 
        select(-ow_rast_10, -index, -x, -y, -`study area`)
      
      # join data
      extr.join <- extr.res.full |> 
        left_join(join.dat) |> 
        select(-soils_rast)
      
      rm(extr.res.full, extr.res, ind.extr, ind.unique, ind.sample)
      gc()
      
      ##
      ## ADD INTERACTIONS
      ##
      # for(z in seq_along(interact.v)) {
      #   # extract the var names - split on " x "
      #   # each column to be interacted - get vector, subtract by its minimum
      #}
      
      
      ##
      ## loop over each variable, and depending on the type, create its sampling distribution
      vars.v <- names(extr.join)
      sample.dist.list <- list()
      cat("\n\n\t\t\tlooping over variables...")
      
      for(k in seq_along(vars.v)) {
        var.tmp <- vars.v[k]
        cat(paste0("\n\t\t\t\tvar = ", var.tmp))
        if(var.tmp == "count") next
        # subset just count and the variable of interest
        dat.subs <- extr.join |> 
          select(count, !!sym(var.tmp) )
        # unfold table into vector
        val.v <- unfold_(dat.subs, verbose=T)
        
        # depending on the class of val.v
        class.tmp <- class(val.v)
        if(verbose) {
          cat("\n\t\t\t\t"); cat(crayon::bgMagenta("class = ", class.tmp))
        }
       
        if(class.tmp == "numeric") {
          cat("\t"); cat(crayon::bgRed("continuous branch"))
          sample.dist.tmp <- hist(val.v,
                                  breaks=n.breaks) # compute histogramn
        } else if(class.tmp == "character" | class.tmp == "factor" | class.tmp == "integer") {
          cat("\t"); cat(crayon::bgBlue("discrete branch"))
          sample.dist.tmp <- table(val.v) # compute table
        } else {
          cat(crayon::red("SHOULDN'T REACH"))
        }
        if(verbose) {
          cat("\n\t\t\t\t")
          cat(crayon::bgWhite("writing..."))
        }
        sample.dist.list[[var.tmp]] <- sample.dist.tmp
      } # end var loop
      if(DBG) cat("\n\t\tsaving...")
      # save list under clean_data/sample_dist/
      if(DBG) cat(paste0("\n\tsave path = ", save.pth) )
      saveRDS(sample.dist.list, file=save.pth)
      rm(sample.dist.list)
      gc()
    } ## end BIG LOOP
    
  } ## end small loop
}
