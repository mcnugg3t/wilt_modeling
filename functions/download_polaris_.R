#'
#'
#'
download_polaris_ <- function(write.wd, verbose=T, DBG=T) {
  if(verbose) cat("\nDOWNLOAD POLARIS...")
  #
  prop.v <- c("alpha", "bd", "clay", "hb", "ksat", "lambda", "n", "om", "ph", "sand", "silt", "theta_r", "theta_s")
  dep.v <- c("0_5", "5_15", "15_30", "30_60", "60_100", "100_200")
  #
  for(i in seq_along(prop.v)) {
    if(verbose) cat(paste0("\n\ti = ", i, ",\tproperty = ", prop.v[i]))
    dir.create("polaris/", prop.v[i])
    #
    url.base <- paste0("http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/", prop.v[i], "/p50/")
    if(DBG) cat(paste0("\n\turl base = ", url.base))
    #
    for(j in seq_along(dep.v)) {
      fn.tmp <- paste0("polaris/", prop.v[i], "/", prop.v[i], "_", dep.v[j], ".tif")
      url.tmp <- paste0(url.base, dep.v[j], "/lat4546_lon-89-88.tif")
      if(DBG) cat(paste0("\n\t\t\turl.tmp = ", url.tmp, "\n\t\t\t\tfn.tmp = ", fn.tmp))
      download.file(url = url.tmp, destfile=fn.tmp, mode="wb")
    }
  }
}

# download_polaris_()
# download.file(url = "http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/theta_s/p50/100_200/lat4546_lon-89-88.tif", destfile="tst.tif", mode="wb")