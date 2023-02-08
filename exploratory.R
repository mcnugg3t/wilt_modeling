
{
  rm(list=ls())
  gc()
}

{
  library(spatstat)
  library(terra)
  library(tidyverse)
  library(assertthat)
} |> suppressPackageStartupMessages()

{
  dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
  dat.30 <- terra::rast("clean_data/joindat_30_1200.tif")
  ow.pts <- terra::vect("clean_data/ow_pts_clean.shp")
  source("functions/construct_ppp_.R")
  ow.ppp.10 <- construct_ppp_(dat.10, ow.pts, T)
  ow.ppp.30 <- construct_ppp_(dat.30, ow.pts, T)
  rm(dat.10, dat.30, ow.pts)
  gc()
  saveRDS(ow.ppp.10, file="clean_data/ow_ppp_10.Rds")
  saveRDS(ow.ppp.30, file="clean_data/ow_ppp_30.Rds")
}

# 1 - Besag's L-function : estimate + gradient
{
  ow.ppp <- readRDS("clean_data/ow_ppp.Rds")
  L.res <- Lest(ow.ppp, correction="translation")
  L.res |> plot()
  Lres.dat <- L.res |> 
    as.data.frame() |> 
    mutate(grad.L = (trans-lag(trans))/(r - lag(r)),
           grad.L.theo = (theo - lag(theo))/(r-lag(r))); Lres.dat
  Lres.dat |> 
    filter(r < 2000) |> 
    select(r, grad.L, grad.L.theo) |> 
    pivot_longer(cols=2:3, names_to = "type", values_to="val") |> 
    ggplot(aes(x=r, y=val, color=type)) + geom_smooth(span=0.10)
}

# generate densities across parameter grid
# move to its own function
{
  
  # setup
  rm(list=ls())
  gc()
  
  kern.v <- c("bw.ppl", "bw.diggle")
  ow.ppp.30 <- readRDS("clean_data/ow_ppp_30.Rds")
  
  source("functions/est_bandwidth_.R")
  bw.res <- est_bandwidth_(kern.v, ow.ppp.30, verbose=T)
  saveRDS(bw.res, file="clean_data/bw_res_30.Rds")
  
  ow.ppp.10 <- readRDS("clean_data/ow_ppp_10.Rds")
  bw.res <- readRDS("clean_data/bw_res_30.Rds")
  adjust.v <- c(0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4)
  
  source("functions/generate_densities_.R")
  generate_densities_(
    ow.ppp = ow.ppp.10,
    bw.list = bw.res,
    adj.v = adjust.v
    )
  
}

# from each density, create a list of point patterns
{
  
  source("functions/create_ppps_.R")
  
  #'
  #'
  #'
  create_ppps_ <- function(density.rast) {
    # resample density raster (if necessary) onto study grid
    # draw Poisson(500) indices
    # at each infection cell, draw infection count according to empirical dist
    # locations assigned by perturbing cell center by x + runif(-10, 10), y + runif(-10, 10)
    # convert to point pattern
  }
  
}

# simulation envelope for L under
# 1 - inhomogeneous poisson
# 2 - cox process
{
  tst.envelope <- envelope(
    ow.ppp, 
    Lest, 
    simulate#=  # list of point patterns
    )
}

# Fry plot locally + globally
{
  
}

#
# compute surprisal for each kde
#
{
  dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
  source("functions/sample_kdes_.R")
  sample_kdes_(
    covar.rast = dat.10, 
    n.sim=1e6, 
    verbose=T, 
    DBG=T)
}

# load data (1) all covariates in raster stack: 10 m x 10 m grid with 1200 m buffer 
#           (2) shapefile of all oak wilt points
{
  source("functions/var_explore_prep_.R")
  
  
  for(i in seq_along(band.v)) {
    band.tmp <- band.v[i]
    cat(crayon::bgRed("\nbandwidth = ", band.tmp))
    res <- var_explore_prep_(
      rast.dat = dat.10, 
      pts.dat = ow.pts, 
      kern.band = band.tmp,
      n.sim = 1000,
      verbose = T, 
      DBG = T
      )
    saveRDS(res, file=paste0("expl_data/res_", band.tmp, ".Rds") )
  }
}

{
  rm(list=ls())
  gc()
  band.v <- c(50, 100, 200, 300, 400, 500, 600)
  for(i in seq_along(band.v)) {
    band.tmp <- band.v[i]
    cat(paste0("\014\nBANDWIDTH = ", band.tmp) )
    dat.tmp <- readRDS(file=paste0("expl_data/res_", band.v[i], ".Rds") )
    plot(dat.tmp[[3]] + labs(title=paste0("bandwidth : ", band.tmp)))
    cat("\n\tenter for next...")
    usr.in <- readline(prompt="  ")
  }
}


  
# if you do any of the Neyman-Scott models with different kernels, can you recover anything like the observed figures?
# create a quadscheme?
