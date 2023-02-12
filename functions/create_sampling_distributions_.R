require(spatstat)
require(tidyverse)
require(assertthat)
require(terra)
require(crayon)
source("functions/unfold_.R")
source("functions/add_interact_.R")
source("functions/mask_density_.R")
#' takes covariate raster and nsim, loops over all KDEs, draws n.sim points with replacement, computes sampling distribution saves
#' also takes (optional) join.dat := single dataframe to join with the data extracted from covar.rast
#' and interact.v := vector of interactions to calculate
create_sampling_distributions_ <- function(covar.rast, join.dat, n.sim, 
                                           interact.v, rm.vars, interact,
                                           n.breaks=100, verbose=T, DBG=F) {
  # setup
  if(verbose) cat( crayon::bgGreen("\n\nFUNCT : sample_kdes_") )
  pth.tmp <- paste0("clean_data/density/")
  kdes.tmp <- list.files( pth.tmp ) # identify kde files
  
  ##
  ## LOOP OVER KDEs
  for(j in seq_along(kdes.tmp)) {
      cat(paste0("\n\n\tj = ", j)); t1 <- Sys.time()
      
      { # construct path
        fn.tmp <- kdes.tmp[j]
        fl.pth <- paste0(pth.tmp, fn.tmp)
        bw.tmp <- fn.tmp |> 
          str_remove("dens_bw_") |> 
          str_remove(".Rds") |> 
          as.numeric()
        save.pth <- paste0("clean_data/sample_dist/samp_dist_bw_", bw.tmp, ".Rds" )
      } # debug print
      {
        if(verbose) cat(paste0("\n\t\treading : ", fl.pth, "\n\t\tsave path : ", save.pth))
        kde.dens <- readRDS( fl.pth )
        if(verbose) cat("\n\t\tmasking...")
        kde.df <- mask_density_(
          dens.in = kde.dens, 
          sa.rast = covar.rast[["study area"]]) |>
          as.data.frame(xy=T) |> 
          mutate(value = value/sum(value, na.rm=T))
        assert_that(sum(kde.df$value<0) == 0) # no negative values
        assert_that((1-sum(kde.df$value))^2 < 1e-7) # sums to 1, or very close
      }
      
      
      { # sample indices n.sim times
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
          select(-soils_rast) ##**hardcoded
        
        rm(extr.res.full, extr.res, ind.extr, ind.unique, ind.sample)
        gc()
      }
      
      ##
      ## ADD INTERACTIONS
      if(interact) {
        if(verbose) cat(crayon::bgCyan("\n\nINTERACTING...\n\n"))
        extr.join <- add_interact_(extr.join, interact.v, T, F)
      }
      
      
      
      ##
      ## loop over each variable, and depending on the type, create its sampling distribution
      vars.v <- names(extr.join)
      sample.dist.list <- list()
      cat("\n\n\t\t\tlooping over variables...")
      
      for(k in seq_along(vars.v)) {
        var.tmp <- vars.v[k]
        cat(paste0("\n\t\t\t\tvar = ")); cat(crayon::underline(crayon::bold(var.tmp)))
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
          n.vals.tmp <- length(  unique(  val.v))
          cat(paste0("\n\n\nn.vals.tmp = ", n.vals.tmp, "\n\n\n"))
          n.breaks.tmp <- if_else(n.vals.tmp  < 100, n.vals.tmp, as.integer(100 ))
          sample.dist.tmp <- hist(val.v,
                                  breaks=n.breaks.tmp) # compute histogramn
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
      t2 <- Sys.time()
      if(verbose) cat(paste0("\n\ntime = ", difftime(t2, t1, units="secs")," s \tn.sim = ", n.sim))
    } ## end BIG LOOP
}
