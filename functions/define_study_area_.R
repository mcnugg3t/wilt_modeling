#'
#'
#'
define_study_area_ <- function(grd.int, buffer.dist, verbose=T, DBG=T) {
  # read OW points
  ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp") 
  # construct buffer and aggregate
  ow.pts.buff <- ow.pts |> 
    buffer(width=buffer.dist) |> 
    aggregate()
  if(grd.int == 10) {
    grd.temp <- rast("mid_data/10/grid/grd_template_10.tif")
    wl2.rast <- rast("mid_data/10/wiscland2/wl2_cls_10.tif")
    manage.rast <- rast("mid_data/10/manage/manage_rast_10.tif")
  } else {
    grd.temp <- rast("mid_data/30/grid/grd_template.tif")
    wl2.rast <- rast("mid_data/30/wiscland2/wl2_cls_30.tif")
    manage.rast <- rast("mid_data/30/manage/manage_rast_30.tif")
  }
  # mask grd.base by OW buffer
  mask.ow.buff <- terra::mask(grd.temp, ow.pts.buff)
  if(DBG) cat( paste0("\n\tinitial sum of values = ", sum(!is.na(values(mask.ow.buff)))) )
  # extract wl2 class values at mask.ow.buff pts
  extr.wl2 <- terra::extract(wl2.rast, crds(mask.ow.buff))
  # identify indices of points where something not oak was extracted - there are ~369k such indices out of the ~813k in extr.wl2
  extr.ind <- which(extr.wl2[,1] != 4230)
  # crop values where not oak
  values(mask.ow.buff)[!is.na(values(mask.ow.buff))][extr.ind] <- NA
  if(DBG) cat( paste0("\n\tpost-oak filt sum of values = ", sum(!is.na(values(mask.ow.buff)))) )
  rm(extr.wl2, extr.ind)
  gc()
  # extract manage polygon data at mask.ow.buff pts
  extr.manage <- terra::extract(manage.rast, crds(mask.ow.buff))
  # identify indices of points where something not USFS was extracted
  extr.ind <- which( extr.manage[,1] %in% c("O", "T", "W") ) 
  # crop values where not USFS
  values(mask.ow.buff)[!is.na(values(mask.ow.buff))][extr.ind] <- NA
  if(DBG) cat( paste0("\n\tpost-private filt sum of values = ", sum(!is.na(values(mask.ow.buff)))) )
  rm(extr.manage, extr.ind)
  gc()
  names(mask.ow.buff) <- "study area"
  values(mask.ow.buff)[!is.na(values(mask.ow.buff))] <- 1
  if(grd.int == 10) {
    terra::writeRaster(mask.ow.buff, filename="mid_data/10/study_area/sa_base.tif")
  } else {
    terra::writeRaster(mask.ow.buff, filename="mid_data/30/study_area/sa_base.tif")
  }
  
}

