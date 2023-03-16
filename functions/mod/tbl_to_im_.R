
#' Takes a tibble recording spatial mapping of (already scaled) environmental variables, vector of variable names, crs of spatial data
#'
#' returns list of spatstat im objects
#'
tbl_to_im_ <- function(tbl_dat, var.v, dat.crs, verbose, DBG=F) {
  require(terra)
  require(tidyverse)
  require(spatstat)
  im.return <- list()
  
  xmin = min(tbl_dat[["x"]], na.rm=T); xmax = max(tbl_dat[["x"]], na.rm=T)
  ymin = min(tbl_dat[["y"]], na.rm=T); ymax = max(tbl_dat[["y"]], na.rm=T)
  if(verbose) cat(paste0("\n\txmin = ", xmin, "\txmax = ", xmax, "\tymin = ", ymin, "\tymax = ", ymax))
  # loop over variables and for each...
  for(i in seq_along(var.v)) {
    var.tmp = var.v[i]
    if(verbose) {cat(paste0("\nVar : ", var.tmp, " (i = ", i, ")"))} # debug print
    dat.subs.mat <- tbl_dat |> # subset data
      dplyr::select(x, y, !!sym(var.tmp) ) |> 
      as.matrix() 
    
    rast.tmp = terra::rast(
      dat.subs.mat, 
      type="xyz", 
      crs=dat.crs, 
      extent= ext(c(xmin, xmax, ymin, ymax))) 
      
    dat.mat.tmp = rast.tmp |> 
      rev() |> 
      flip(direction="horizontal") |>  
      as.matrix(wide=T)
    
    #dat.mat.tmp |> image()
    
    im.tmp = im(
      mat = dat.mat.tmp, 
      xrange = c(xmin, xmax), 
      yrange = c(ymin, ymax), 
      unitname=c("meter", "meters"))
    
    #im.tmp |> plot()
    if(str_detect(var.tmp, pattern=" x ")) var.tmp = str_replace_all(var.tmp, pattern=" x ", replacement="_x_")
    im.return[[var.tmp]] = im.tmp
  }
  return(im.return)
}

# tbl_dat = dat.mod.30; var.v = vars.in; dat.crs = terra::crs('EPSG:3071'); verbose=T; DBG=T

