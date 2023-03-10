require(terra)
require(crayon)
#'
#'
#'
define_study_area_ <- function(grd.int, wilt.buffer.dist, oak.buffer.dist, verbose=T, DBG=T) {

  # if using 10 grid...
  if(grd.int == 10) {
    grd.temp <- rast("mid_data/10/grid/grd_template_10.tif")
    wl2.rast <- rast("mid_data/10/wiscland2/wl2_cls_10.tif")
    manage.rast <- rast("mid_data/10/manage/manage_rast_10.tif")
  # otherwise
  } else {
    grd.temp <- rast("mid_data/30/grid/grd_template.tif")
    wl2.rast <- rast("mid_data/30/wiscland2/wl2_cls_30.tif")
    manage.rast <- rast("mid_data/30/manage/manage_rast_30.tif")
  }
  ##
  ## mask grd.base in a series of steps 
  ##
  
  ##
  ## 1 : using OW buffer - e.g. include cells points within 1200 m of observed infection
  if(verbose) cat(crayon::bgBlue("\n\tconstructing wilt buffer = ", wilt.buffer.dist, " m..."))
  # read OW points
  ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp") 
  # construct buffer around each ow pt and aggregate to one vector object
  ow.pts.buff <- ow.pts |> 
    buffer(width=wilt.buffer.dist) |> 
    aggregate()
  if(verbose) cat("\tmasking...")
  mask.ow.buff.1 <- terra::mask(grd.temp, ow.pts.buff)
  if(DBG) cat( paste0("\n\t\tcurrent count cells with values = ", sum(!is.na(values(mask.ow.buff.1)))) )
  ##
  ## 2 : using oak buffer
  if(verbose) cat(crayon::bgBlue("\n\tconstructing oak buffer = ", oak.buffer.dist, " m...") )
  wl2.oak.subs <- ifel(wl2.rast == 4230, 4230, NA) # get oak cells
  oak.buff <- as.polygons(wl2.oak.subs) |>  
    buffer(width = oak.buffer.dist) |>  # buffer using param width
    aggregate() # aggregate
  if(verbose) cat("\tmasking...")
  mask.ow.buff.2 <- terra::mask(mask.ow.buff.1, oak.buff) # next mask
  if(DBG) cat( paste0("\n\t\tcurrent count cells with values = ", sum(!is.na(values(mask.ow.buff.2)))) )
  rm(mask.ow.buff.1)
  gc()
  
  ##
  ## 3  extract manage polygon data at mask.ow.buff pts
  if(verbose) cat(crayon::bgBlue("\n\textracting manage data..."))
  cnnf.mask <- ifel(
                  manage.rast %in% c("O", "T", "W"), 
                  NA,
                  manage.rast)
  mask.ow.buff.3 <- mask(mask.ow.buff.2, cnnf.mask)
  if(DBG) cat( paste0("\n\t\tcurrent count cells with values = ", sum(!is.na(values(mask.ow.buff.3)))) )
  if(grd.int == 10) {
    names(mask.ow.buff.3) <- "study area"
    terra::writeRaster(mask.ow.buff.3, filename="mid_data/10/study_area/sa_base.tif", overwrite=T)
  } else {
    names(mask.ow.buff.3) <- "study area"
    terra::writeRaster(mask.ow.buff.3, filename="mid_data/30/study_area/sa_base.tif", overwrite=T)
  }
  
}
# clean.dat <- terra::rast("clean_data/joindat_10_1200.tif")
# manage.tmp <- clean.dat[["manage_rast_10]]
# grd.int = 30; wilt.buffer.dist = 500; oak.buffer.dist=40; verbose=T; DBG=T
