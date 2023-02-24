require(terra)
require(tidyverse)
require(crayon)
require(assertthat)
source("functions/sub_functions/create_one_ppp_list_.R")
source("functions/mask_density_.R")
#' top level function which calls create_one_ppp_list_
#' 
#'
create_ppps_ <- function(n.sim, verbose=T, DBG=F) {
  if(verbose) cat(crayon::bgGreen("\nFUNCT : create_ppps_"))
  
  # INIT
  {
    file.v <- list.files("clean_data/density/") # density files
    complete.v <- list.files("clean_data/sample/") # already calculated ppp lists
    cln.dat <- terra::rast("clean_data/joindat_10_1200.tif") # load full clean data
    ow.df <- cln.dat[["ow_rast_10"]] |>  # calc ow cell locations and counts
      as.data.frame(xy=T) |> 
      filter(ow_rast_10 > 0)
    sa.rast = cln.dat[["study area"]]
    ow.ppp = readRDS("clean_data/ow_ppp_10.Rds") # load ppp
    wilt.owin <- ow.ppp$window
    rm(cln.dat, ow.ppp)
    gc()
  }
  
  # LOOP
  for(j in seq_along(file.v)) { # over each density file
    #
    fl.tmp <- file.v[j]
    
    cat(paste0("\n\nj = ", j, "\n\tfl.tmp = ", fl.tmp))
    path.tmp <- paste0("clean_data/density/", fl.tmp)
    bw.num <- fl.tmp |> 
      str_remove("dens_bw_") |> 
      str_remove(".Rds") |> 
      as.numeric()
    save.sml <- paste0("ppp_list_bw_", bw.num, ".Rds" )
    if(save.sml %in% complete.v) next
    save.tmp <- paste0("clean_data/sample/", save.sml)
    if(verbose) cat(paste0("\n\tsave.tmp = ", save.tmp))
    #
    # read density
    if(verbose) cat(paste0("\n\treading density.df, masking..."))
    
    density.df <- readRDS(path.tmp) |> 
      mask_density_(sa.rast = sa.rast) |> 
      as.data.frame(xy=T) |> 
      mutate(value = if_else(
        condition = (is.na(value) | value<0),
        true=0,
        false=value)) |> 
      mutate(value = value/sum(value, na.rm=T))
    assert_that( (sum(density.df$value)-1)^2 < 1e-9 )
    assert_that(sum(is.na(density.df$value)) == 0)
    
    if(verbose) cat(crayon::bgMagenta("\tcreating ppp list..."))
    
    ppp.list <- create_one_ppp_list_(
      density.df,
      ow.df,
      n.sim = n.sim,
      owin.in = wilt.owin,
      n.ow.cells = nrow(ow.df)
    )
    
    cat("\n\t\tsaving...")
    saveRDS(ppp.list, file=save.tmp)
  }
}
  
