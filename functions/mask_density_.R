require(terra)
#' takes density object and raster with study area bitmask
#'
#'
mask_density_ <- function(dens.in, sa.rast) {
  # convert dens.in to raster using crs of sa.rast
  # dens.rast 
  t1 <- dens.in |> 
    as.data.frame() |> 
    as.matrix() |> 
    terra::rast(type="xyz", crs=crs(sa.rast)) |> 
    resample(sa.rast)
  t2 <- ifel(t1 < 1e-7, NA, t1)
  return(t2)
}

