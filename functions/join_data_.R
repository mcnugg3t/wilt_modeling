require(crayon)
#'
#'
#'
join_data_ <- function(base.grd, fold.skip, var.skip, verbose=T, DBG=T) {
  if(verbose) cat("\nJOINING DATA...")
  
  # generate folders - exclude .mxd
  mid.fold <- list.files(path="mid_data")
  skip <- which(str_detect(mid.fold, pattern="(.mxd)") )
  if(DBG) {
    cat("\n\tmid_data folders = ")
    cat(crayon::blue(mid.fold, collapse = "    "))
  }
  
  # loop over folders
  for(i in seq_along(mid.fold)) {
    if(i %in% skip | mid.fold[i] == "grid" | mid.fold[i] %in% fold.skip) next # skip grid folder
    if(verbose) cat(paste0("\n\n\tfolder = ", mid.fold[i], "\t(i = ", i, ")"))
    nest.fls.tmp <- list.files(path = paste0("mid_data/", mid.fold[i]))
    nest.fls.clean <- nest.fls.tmp[str_detect(nest.fls.tmp, pattern="(.tif$)")] # remove files that don't end in .tif
    if(DBG) {
      cat("\n\t\tnest.fls.tmp = ")
      cat(crayon::green(nest.fls.clean, collapse="  ,  "))
    }
    # loop over filenames that DO end in .tif - read each, resample it to the base grid (base.grd) appropriately
    for(j in seq_along(nest.fls.clean)) {
      nest.fl.tmp <- nest.fls.clean[j]
      var.tmp <- nest.fl.tmp |> str_remove(pattern=".tif")
      if(var.tmp %in% var.skip) next
      # project base.grd to foreign crs if necessary
      fl.read.tmp <- paste0( "mid_data/", mid.fold[i], "/", nest.fls.clean[j] ) # assemble file path
      if(DBG) cat(paste0("\n\t\t\tfl.read.tmp = ", fl.read.tmp, "\tvar.tmp = ", var.tmp, "\n\t\t\t\tmasking...")) # dbg print
      rast.tmp <- terra::rast(fl.read.tmp) # read file
      if(crs(rast.tmp) != crs(base.grd)) {
        cat(paste0("\n\ncrs of rast = ", paste0(crs(rast.tmp), collapse=" ")))
        base.grd.tmp <- base.grd |> 
          terra::project(crs(rast.tmp))
      } else {
        base.grd.tmp <- base.grd
      }
      if(DBG) cat(paste0("\n\t\t\t\textracting..."))
      extr.vals <- extract(rast.tmp , crds(base.grd.tmp) )
      write.inds <- which( !is.na(values(base.grd)) )
      base.grd[[var.tmp]] <- 0
      base.grd[[var.tmp]][write.inds] <- extr.vals
      rm(rast.tmp, extr.vals, write.inds)
      gc()
    }
  }
  return(base.grid)
}

# test
#grd.temp.10 <- terra::rast("mid_data/grid/grd_template_10.tif")

#tst.fls <- list.files("mid_data/groundwater/")
#str_detect(tst.fls, pattern="(.tif$)")
# tif.tst1 <- terra::rast("mid_data/groundwater/gw_no_smooth.tif")
# str(values(tif.tst1))
# values(tif.tst1)[,1] |> class()

