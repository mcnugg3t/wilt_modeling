
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

# load data (1) all covariates in raster stack: 10 m x 10 m grid with 1200 m buffer 
#           (2) shapefile of all oak wilt points
{
  source("functions/var_explore_.R")
  dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
  ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp")
  band.v <- c(50, 100, 200, 300, 400, 500, 600)
  for(i in seq_along(band.v)) {
    band.tmp <- band.v[i]
    cat(crayon::bgRed("\nbandwidth = ", band.tmp))
    res <- var_explore_(
      rast.dat = dat.10, 
      pts.dat = ow.pts, 
      kern.band = band.tmp,
      n.sim = 1000,
      verbose = T, 
      DBG = T)
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
