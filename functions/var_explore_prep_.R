require(spatstat)
require(terra)
require(tidyverse)
require(assertthat)
source("functions/construct_ppp_.R")
#'
#'
#'
var_explore_prep_ <- function(rast.dat, pts.dat, kern.bw, interact.v, n.sim, verbose=T, DBG=T) {
  if(verbose) cat("\n\nVAR EXPLORE PREP (step 1 of 2)...")
  if(DBG) cat("\n\tsetup...")
  
  # construct ppp from study area and points
  ow.ppp <- construct_ppp_(rast.dat, pts.dat)
  
  ### calculate density
  if(verbose) cat("\n\tcalc density...")
  t1 <- Sys.time()
  ow.dens <- density(ow.ppp, sigma=kern.bw)
  if(verbose) cat(paste0("\ttime : ", difftime(Sys.time(), t1, units="secs")))
  
  ### manipulate density data  -> tibble, normalize, resample onto study area -> data.frame
  if(DBG) cat("\n\t\tconverting, normalizing...")
  
  ow.adj <- ow.dens |> 
    as_tibble() |> 
    mutate(val.norm = value/sum(value)) |> 
    select(-value)
  
  ow.adj.rast <- rast(ow.adj, crs=crs(study.area) ) |> 
    resample(study.area)
  
  ow.adj.tbl <- ow.adj.rast |> 
    as.data.frame(xy=T) |> 
    mutate(val.norm = abs(val.norm))
  
  #plot(ow.adj.rast[["val.norm"]])
  #plot(pts.dat, cex=0.5, col=rgb(red=0, green=0, blue=0, alpha=0.3), add=T)
  rm(ow.adj)
  gc()
  
  check.sum <- ow.adj.tbl |> 
    select(val.norm) |> 
    sum()
  
  if(DBG) cat(paste0("\n\tdensity sums to 1: " , assert_that(abs(check.sum-1) < 0.0001) ))
  
  if(verbose) cat("\n\tsampling indices ~ kernel density....")
  
  ### sample indices nsim times
  s.list <- list()
  for(i in 1:n.sim) {
    cat(paste0("\n\014index sampling: ", round(i/n.sim, 3)*100, " %"))
    s.tmp <- sample(
      x = 1:nrow(ow.adj.tbl), 
      size = rpois(n=1, lambda= ow.ppp$n), 
      replace=T,
      prob=ow.adj.tbl$val.norm)
    s.list[[i]] <- s.tmp
  }
  
  if(verbose) cat("\n\textracting covariates ~ sampled indices")
  # extract covariate values from data, compute mean and sd of each, store it
  # init empty tibble
  store.tbl <- tibble(
    var = character(),
    mu = numeric(),
    sig = numeric())
  
  source("functions/detect_remove_vars_.R")
  rmv.vars <- detect_remove_vars_(
    var.names = names(rast.dat),
    patterns = c("study area", "manage_rast", "wl2_cls", "soils_rast", "ow_rast"),
    verbose=F
  )
  
  # for each sample of indices...
  for(i in seq_along(s.list)) {
    cat(paste0("\n\014extracting covars: ", round(i/length(s.list), 3)*100, " %"))
    s.tmp <- s.list[[i]] # extract sample vector
    # get coords associated with indices, -> matrix
    crds.mat <- ow.adj.tbl |>  
      slice(s.tmp) |>  
      select(x, y) |> 
      as.matrix()
    
    # extract covar values from dat.10
    dat.subs <- terra::extract(rast.dat, crds.mat) |> 
      as_tibble() |> 
      select(-all_of(rmv.vars)) 
    
    ##
    ## add interactions
    ##
    interact.tmp <- interact.v[1]; interact.tmp
    terms.tmp <- str_split(interact.tmp, pattern=" x "); terms.tmp[[1]]
    vars.v.tmp <- terms.tmp[[1]]; vars.v.tmp
    dat.tst <- dat.subs |> 
      mutate(!!sym(interact.tmp) := !!sym(vars.v.tmp[1]) * !!sym(vars.v.tmp[2]) )
    
    for(j in seq_along(interact.v)) {
      interact.tmp <- interact.v[j]
      terms.tmp <- str_split(interact.tmp, pattern=" x ")
      vars.v.tmp <- terms.tmp[[1]]
      min.tmp.1 <- min(dat.subs[[vars.v.tmp[1]]])
      min.tmp.2 <- min(dat.subs[[vars.v.tmp[2]]])
      if(DBG) cat( paste0("\n\n\tcalc interaction: ", interact.tmp, "\n\t\t\tmin1 = ", min.tmp.1, "\tmin2 = ", min.tmp.2 ) )
      
      dat.subs <- dat.subs |> 
        mutate(
          var1.tmp = !!sym(vars.v.tmp[1]),
          var2.tmp = !!sym(vars.v.tmp[2]),
          !!sym(interact.tmp) := (var1.tmp - min.tmp.1 ) * (var2.tmp - min.tmp.2 )
          ) |> 
        select(-var1.tmp, -var2.tmp)
    }
    
    # calc mean and sd of each variable
    mean.v <- apply(dat.subs, MARGIN=2, FUN=mean, na.rm=T)
    sd.v <- apply(dat.subs, MARGIN=2, FUN=sd, na.rm=T )
    
    # store in tbl
    store.tbl <- store.tbl |> 
      add_row(
        var = names(mean.v),
        mu = mean.v,
        sig = sd.v
      )
  }  
  
  return(store.tbl)
}

# { # FOR DEBUGGING
#     dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
#     ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp")
#     band.v <- c(50, 100, 200, 300, 400, 500, 600)
#     i <- 1
#     band.tmp <- band.v[i]
#     rast.dat = dat.10 
#     pts.dat = ow.pts 
#     kern.band = band.tmp
#     interact.v = c("conv_ind x elev")
#     n.sim = 100
#     verbose = T 
#     DBG = T
#     j <- 1
# }

