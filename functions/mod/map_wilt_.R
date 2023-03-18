require(terra)
require(crayon)
#'
#'
map_wilt_ <- function(pred_tbl, DBG) {
  if(DBG) cat(crayon::bgGreen("\n\n\t\tmap_wilt_ ... "))
  # load raster data
  mod.dat.rast = terra::rast("mod_data/compare/wilt_split_rast.tif")
  wilt.test = mod.dat.rast$wilt_test; names(wilt.test) = "wilt"
  d = dim(mod.dat.rast)
  # form raster out of input pred_tbl and resample to loaded raster
  pred.rast = rast(pred_tbl, type="xyz", crs=crs(mod.dat.rast))
  pred.rast = terra::resample(pred.rast, mod.dat.rast)
  names(pred.rast) = "pred"
  
  return.df = c(wilt.test, pred.rast) |> 
    as.data.frame(xy=T) |> 
    mutate(pred = if_else(is.na(pred), 0, pred))
}