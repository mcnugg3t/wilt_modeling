require(terra)
require(crayon)
#'
#'
map_wilt_30_ <- function(pred_tbl, DBG) {
  if(DBG) cat(crayon::bgGreen("\n\n\t\tmap_wilt_30_ ... "))
  # load raster data
  mod.dat.rast = terra::rast("mod_data/compare/wilt_split_rast.tif") # on the 30-m grid
  wilt.test = mod.dat.rast$wilt_test
  names(wilt.test) = "wilt"
  d = dim(mod.dat.rast)
  # form raster out of input pred_tbl and resample to loaded raster
  pred.rast = vect(pred_tbl, geom=c("x", "y"), crs=crs(mod.dat.rast)) |> 
    rasterize(mod.dat.rast, field="pred")
  names(pred.rast) = "pred"
  return.df = c(wilt.test, pred.rast) |> 
    as.data.frame(xy=T) |> 
    mutate(pred = if_else(is.na(pred), 0, pred))
  return(return.df)
}

#'
#'
map_wilt_10_ <- function(pred_tbl, DBG) {
  if(DBG) cat(crayon::bgGreen("\n\n\t\tmap_wilt_10_ ... "))
  # load raster data
  mod.dat.rast = terra::rast("mod_data/compare/wilt_split_rast_10.tif") # on the 10-m grid
  wilt.test = mod.dat.rast$wilt_test
  names(wilt.test) = "wilt"
  d = dim(mod.dat.rast)
  # form raster out of input pred_tbl and resample to loaded raster
  pred.rast = vect(pred_tbl, geom=c("x", "y"), crs=crs(mod.dat.rast)) |> 
    rasterize(mod.dat.rast, field="pred")
  names(pred.rast) = "pred"
  return.df = c(wilt.test, pred.rast) |> 
    as.data.frame(xy=T) |> 
    mutate(pred = if_else(is.na(pred), 0, pred))
  return(return.df)
}