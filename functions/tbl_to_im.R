
#' Takes a tibble recording spatial mapping of (already scaled) environmental variables, vector of variable names, crs of spatial data
#'
#' returns list of spatstat im objects
#'
tbl_to_im <- function(tbl_dat, var.v, dat.crs, tmp.path, DEBUG=F) {
  require(stars)
  require(raster)
  require(tidyverse)
  require(spatstat)
  im.return <- list()
  for(var in var.v) { # for each variable
    if(DEBUG) {cat(paste0("\nVar : ", var))} # debug print
    dat.tmp <- tbl_dat %>% dplyr::select(x, y, !!var) # subset data
    if(DEBUG) {cat(paste0("\n\tconverting to stars..."))} # debug print
    img.tmp <- st_as_stars(dat.tmp) # convert to stars object
    st_crs(img.tmp) <- dat.crs # set CRS of stars object
    path.tmp <- paste0(tmp.path, var, ".tif")
    if(DEBUG) {cat(paste0("\n\t\t writing stars at...\t", path.tmp))}
    write_stars(img.tmp, dsn=path.tmp) # write raster
    if(DEBUG) {cat(paste0("\n\t\t\t reading as raster obj from...", path.tmp))}
    rast.tmp <- raster(path.tmp) # read
    mat.tmp <- rast.tmp %>% as.matrix() # convert to matrix
    if(DEBUG) {cat(paste0("\n\t\t\t\t transmat..."))} # debug print
    mat.conv <- transmat(mat.tmp, from="European", to="spatstat") # transmat to spatstat format
    im.return[[var]] <- im(mat.conv, xrange=c(639794, 652604), yrange=c(529893.8, 547893.8), unitname=c("meter", "meters"))
  }
  return(im.return)
}



