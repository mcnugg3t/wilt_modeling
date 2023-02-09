
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
  #ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp")
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
  ow.ppp <- readRDS("clean_data/ow_ppp_10.Rds")
  density.df <- readRDS("clean_data/density/bw.diggle/dens_adj_1.Rds") |> 
    as.data.frame() |> 
    mutate(val.old = value) |> 
    mutate(value = if_else( ((0-value)^2)<0.1e-7, 0, value ),
           value = value/sum(value, na.rm=T))
  density.df$value |> sum()
  density.df |> filter(value>0) |> nrow() # 3637 non-zero vals
  
  cln.dat <- terra::rast("clean_data/joindat_10_1200.tif")
  ow.df <- cln.dat[["ow_rast_10"]] |> 
    as.data.frame(xy=T) |> 
    filter(ow_rast_10 > 0)
  source("functions/create_ppps_.R")
  ppp.list <- create_ppps_(
    density.df,
    ow.df,
    n.sim=100,
    owin.in = ow.ppp$window
  )
  
  saveRDS(ppp.list, file="clean_data/ppp/ppp_test.Rds")
  #plot(density.rast)
  #plot(study.grid[["study area"]], add=T)
  
  
  
}

# simulation envelope for L under
# 1 - inhomogeneous poisson
# 2 - cox process
{
  ow.ppp <- readRDS("clean_data/ow_ppp_10.Rds")
  ppp.list <- readRDS("clean_data/ppp/ppp_test.Rds")
  tst.envelope <- envelope(
    ow.ppp, 
    fun=Lest, 
    funargs=list(correction="translation"),
    simulate= ppp.list # list of point patterns
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
<<<<<<< HEAD
dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp")
}

# construct spatstat objects
{ 
# study area bit-mask is one of the layers in the raster stack
study.area <- dat.10[["study area"]]
# convert study area raster object to dataframe with coordinates
sa.df <- study.area |> 
  as.data.frame(xy=T)
# construct spatstat owin object 
ext.sa <- ext(study.area)
ow.win.10 <- owin(
  xrange = c(ext.sa[1], ext.sa[2]),
  yrange = c(ext.sa[3], ext.sa[4]),
  mask = sa.df[,1:2]
  )
# construct spatstat ppp object
ow.crds <- crds(ow.pts)
ow.ppp <- ppp(
  x = ow.crds[,1], 
  y = ow.crds[,2], 
  window = ow.win.10) # 59 points outside of window + some duplicated
dummy.scheme <- ppp(
  x = sa.df$x, 
  y = sa.df$y)
ow.quad <- quadscheme(data=ow.ppp, dummy=dummy.scheme, method="grid")
}
#

# construct density and kppm objects
{
ow.dens <- density(ow.ppp, sigma=150)
#ow.kppm.thom <- kppm(ow.ppp, clusters="Thomas", rmax=1000)
#ow.kppm.cauch <- kppm(ow.ppp, clusters="Cauchy", rmax=1000)
#ow.kppm.vargam <- kppm(ow.ppp, clusters="VarGamma", rmax=1000)
}

# plot
{
  plot(ow.dens)
  plot(ow.dens^(1/2))
  plot(ow.dens^(2/3))
  #ow.kppm.cauch |> predict() |> plot()
  # plot((ow.dens)^(1/2))
  # plot(log.dens)
}

# kernel density -> tibble, normalize, resample onto study area -> data.frame
{
ow.adj <- ow.dens^(2/3) |> 
  as_tibble() |> 
  mutate(val.norm = value/sum(value)) |> 
  select(-value)
ow.adj.rast <- rast(ow.adj, crs=crs(study.area) ) |> 
  resample(study.area)
ow.adj.tbl <- ow.adj.rast |> 
  as.data.frame(xy=T)
rm(ow.adj)
gc()
}

# check density estimate
{
tst.sum <- values(ow.adj.rast)[!is.na(values(ow.adj.rast))] |> sum()
assert_that(abs(tst.sum-1) < 0.0001)
plot(ow.adj.rast[["val.norm"]])
plot(ow.pts, cex=0.5, col=rgb(red=0, green=0, blue=0, alpha=0.3),  add=T)
#rm(ow.adj.rast)
gc()
}

# sample indices nsim times
{
  n.sim <- 1000
  s.list <- list()
  for(i in 1:n.sim) {
    cat(paste0("\n\014", round(i/n.sim, 3)*100, " %"))
    s.tmp <- sample(
      x = 1:nrow(ow.adj.tbl), 
      size = rpois(n=1, lambda=nrow(ow.crds)), 
      replace=T,
      prob=ow.adj.tbl$val.norm
    )
    s.list[[i]] <- s.tmp
=======
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
>>>>>>> f948557d8a03caed183afa4c97d699aeab02fc06
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
