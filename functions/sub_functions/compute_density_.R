require(spatstat)
require(crayon)
#'
#'
#'
compute_density_ <- function(ppp.in, bw.in, verbose=T, DBG=T) {
  cat("\n\n\t")
  cat(crayon::bgGreen("SUB-FUNCT: compute_density_...") )
  if(verbose) cat(crayon::bgBlue("\tbw in = ", bw.in) )
  
  ### density
  if(verbose) cat("\n\tcalc density...")
  t1 <- Sys.time()
  dens.res <- density(ppp.in, sigma=bw.in, positive=T)
  if(verbose) cat(paste0("\ttime : ", difftime(Sys.time(), t1, units="secs")))
  
  return(dens.res)
}

# { # FOR DEBUGGING
#     dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
#     ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp")
#     band.v <- c(50, 100, 200, 300, 400, 500, 600)
#     i <- 1
#     band.tmp <- band.v[i]
#     rast.dat = dat.10 
#     pts.dat = ow.pts 
#     kern.band = band.tmp
#     interact.v = c("conv_ind x elev")
#     n.sim = 100
#     verbose = T 
#     DBG = T
#     j <- 1
# }
# difftime(t2, t1, units="secs")
