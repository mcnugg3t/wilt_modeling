require(crayon)
#'
#'
#'
join_data_ <- function(base.grd, verbose=T, DBG=T) {
  if(verbose) cat("\nJOINING DATA...")
  
  # generate folders - exclude .mxd
  mid.fold.raw <- list.files(path="mid_data")
  mid.fold <- str_remove_all(mid.fold.raw, pattern="(.mxd)")
  if(DBG) {
    cat("\n\tmid_data folders = ")
    cat(crayon::blue(mid.fold, collapse = "    "))
  }
  
  # loop over folders
  for(i in seq_along(mid.fold)) {
    if(mid.fold[i] == "grid") next # skip grid folder
    if(verbose) cat(paste0("\n\n\tfolder = ", mid.fold[i], "\t(i = ", i, ")"))
    nest.fls.tmp <- list.files(path = paste0("mid_data/", mid.fold[i]))
    nest.fls.clean <- nest.fls.tmp[str_detect(nest.fls.tmp, pattern="(.tif$)")] # remove files that don't end in .tif
    if(DBG) {
      cat("\n\t\tnest.fls.tmp = ")
      cat(crayon::green(nest.fls.clean, collapse="  ,  "))
    }
    # loop over filenames that DO end in .tif - read each, resample it to the base grid (base.grd) appropriately
    for(j in seq_along(nest.fls.clean)) {
      fl.read.tmp <- paste0( "mid_data/", mid.fold[i], "/", nest.fls.clean[j] )
      if(DBG) cat(paste0("\n\t\t\tfl.read.tmp = ", fl.read.tmp))
      # determine behavior based on type stored
      tif.tmp <- terra::rast(fl.read.tmp)
      type.tmp <- values(tif.tmp)[,1] |> class()
      if(DBG) cat(paste0("\n\t\t\t\tTYPE = ", toString(type.tmp)))
      # resample onto grid
    }
  }
  
}

# test
grd.temp.10 <- terra::rast("mid_data/grid/grd_template_10.tif")
join_data_(base.grd = grd.temp.10)
#tst.fls <- list.files("mid_data/groundwater/")
#str_detect(tst.fls, pattern="(.tif$)")
# tif.tst1 <- terra::rast("mid_data/groundwater/gw_no_smooth.tif")
# str(values(tif.tst1))
# values(tif.tst1)[,1] |> class()
tif.tst2 <- terra::rast("mid_data/manage/manage_rast.tif")
v.tmp <- values(tif.tst2)[,1]
c.tmp <- crds(tif.tst2)
tst.tbl <- tibble(
  x = c.tmp[,1],
  y = c.tmp[,1],
  val = v.tmp[!is.na(v.tmp)]
) |> arrange(desc(y), x  )

tst.tbl |> arrange(x, desc(y) )
