
rm(list=ls())
gc()

{
  library(spatstat)
  library(terra)
} |> suppressPackageStartupMessages()

dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp")

ow.crds <- crds(ow.pts)

study.area <- dat.10[["study area"]]

sa.df <- study.area |> 
  as.data.frame(xy=T)

ow.win.10 <- owin(
  xrange = c(ext.sa[1], ext.sa[2]),
  yrange = c(ext.sa[3], ext.sa[4]),
  mask = sa.df[,1:2]
  )

ow.ppp <- ppp(x = ow.crds[,1], y = ow.crds[,2], window = ow.win.10)

tst <- plot(ow.ppp)
