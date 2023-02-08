#' construct spatstat ppp object
#'
#'
construct_ppp_ <- function(rast.dat, pts.dat, verbose=T) {
  if(verbose) cat("\n\nCONSTRUCTING PPP...")
  study.area <- rast.dat[["study area"]]
  sa.df <- study.area |> 
    as.data.frame(xy=T)
  ext.sa <- ext(study.area)
  ow.win <- owin(
    xrange = c(ext.sa[1], ext.sa[2]),
    yrange = c(ext.sa[3], ext.sa[4]),
    mask = sa.df[,1:2])
  rm(sa.df)
  gc()
  # extract coordinates
  pts.crds <- crds(pts.dat)
  # construct spatstat ppp object
  if(verbose) {
    ow.ppp <- ppp(x = pts.crds[,1], y = pts.crds[,2], window = ow.win)
  } else {
    ow.ppp <- ppp(x = pts.crds[,1], y = pts.crds[,2], window = ow.win) |> 
      suppressWarnings()
  }
  return(ow.ppp)
}
