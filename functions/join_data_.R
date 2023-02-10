require(crayon)
#'
#'
#'
join_data_ <- function(grd.int, verbose=T, DBG=T) {
  # SETUP
  if(verbose) cat("\nJOINING DATA...")
  fold.v <- c("groundwater", "bedrock", "manage", "polaris", "saga", "soils", "wilt", "wiscland2")
  sa.base <- rast( paste0("mid_data/", grd.int, "/study_area/sa_base.tif") )
  to.return <- sa.base
  
  # LOOP OVER FOLDERS
  for(i in seq_along(fold.v)) {
    fold.tmp <- fold.v[i]
    if(verbose) cat(paste0("\n\n\tfolder = ", fold.tmp, "\t(i = ", i, ")"))
    nest.fls.tmp <- list.files(path = paste0("mid_data/", grd.int, "/", fold.v[i])) # list files
    nest.fls.clean <- nest.fls.tmp[str_detect(nest.fls.tmp, pattern="(.tif$)")] # remove files that don't end in .tif
    if(DBG) {
      cat("\n\t\tnest.fls.tmp = ")
      cat(crayon::green(nest.fls.clean, collapse="  ,  "))
    }
    
    # loop over filenames that DO end in .tif - read each, resample it to the base grid (base.grd) appropriately
    for(j in seq_along(nest.fls.clean)) {
      nest.fl.tmp <- nest.fls.clean[j]
      var.tmp <- nest.fl.tmp |> str_remove(pattern=".tif")
      #if(var.tmp %in% var.skip) next
      if( var.tmp %in% c("wl2_cls_10", "wl2_cls_30") ) {
        method.tmp <- "near"
      } else {
        method.tmp <- "bilinear"
      }
      # assemble file path
      fl.read.tmp <- paste0( "mid_data/", grd.int, "/", fold.v[i], "/", nest.fls.clean[j] ) 
      if(DBG) cat(paste0("\n\n\t\t\tfl.read.tmp = ", fl.read.tmp, "\t\twhere var.tmp = ", var.tmp)) # dbg print
      # read file
      rast.tmp <- terra::rast(fl.read.tmp) |> 
        mask(sa.base)
      
      names(rast.tmp) <- var.tmp
      
      to.return <- c(to.return, rast.tmp)
      rm(rast.tmp)
      gc()
    }
  }
  return(to.return)
}

# test
#grd.temp.10 <- terra::rast("mid_data/grid/grd_template_10.tif")

#tst.fls <- list.files("mid_data/groundwater/")
#str_detect(tst.fls, pattern="(.tif$)")
# tif.tst1 <- terra::rast("mid_data/groundwater/gw_no_smooth.tif")
# str(values(tif.tst1))
# values(tif.tst1)[,1] |> class()

