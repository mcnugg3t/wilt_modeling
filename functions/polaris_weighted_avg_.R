require(terra)
require(tidyverse)
#'
#' function assumes it will find files in folder, e.g. wd + raw_data/polaris/bd/bd_0_5.tif
#'
polaris_weighted_avg_ <- function(str.var, grd.temp, verbose=T, DBG=F) {
  if(verbose) cat(paste0("\n\nWEIGHTED AVERAGE...\t", str.var))
  #
  depths.v <- c(5, 10, 15, 30, 40)
  dn.v <- c("_0_5", "_5_15", "_15_30", "_30_60", "_60_100")
  #
  for(i in seq_along(depths.v)) {
    if(verbose) cat(paste0("\n\ti = ", i))
    fn.tmp <- paste0("raw_data/polaris/", str.var, "/", str.var, dn.v[i], ".tif")
    if(DBG) cat(paste0("\n\tfn = ", fn.tmp, "\n\t\tdepth = ", dn.v[i]))
    if(i == 1) {
      if(DBG) cat("\n\t\treading i = 1")
      return.rast <- terra::rast(fn.tmp) * depths.v[i]
    } else {
      if(DBG) cat()
      rast.tmp <- terra::rast(fn.tmp) * depths.v[i]
      return.rast <- return.rast + rast.tmp
    }
  }
  if(verbose) cat(paste0("\n\n count NA = ", sum(is.na(values(return.rast)))) )
  return.rast <- return.rast/100
  if(DBG) cat("\n\nwarping to grd.temp...")
  return.rast <- return.rast |> 
    project(crs(grd.temp)) |> 
    resample(grd.temp, method="bilinear")
  return(return.rast)
}


