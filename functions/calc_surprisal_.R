require(tidyverse)
require(crayon)
source("functions/sub_functions/calc_surprisal_hist_.R")
source("functions/sub_functions/calc_surprisal_tab_.R")
#'
#'
#'
calc_surprisal_ <- function(dens.files, ow.pts.dat, bw.select.v, grd.int, verbose=T, DBG=F) {
  var.v <- names(ow.pts.dat)
  
  return.dat <- tibble(
    var = character(),
    bw = numeric(),
    surprisal = numeric()
  )
  
  for(i in seq_along(dens.files)) { # for each density file
    fn.tmp <- dens.files[i] # extr filename
    bw.tmp <- fn.tmp |>  # extr bandwidth from filename
      str_remove("samp_dist_bw_") |>
      str_remove(".Rds") |> 
      as.numeric()
    if(!(bw.tmp %in% bw.select.v)) next
    if(verbose) cat(paste0("\n\n\ti = ", i, "\tfn.tmp = ", fn.tmp, "   bw = ", bw.tmp)) # debug print
    
    fp.tmp <- paste0("clean_data/sample_dist/", grd.int, "/", fn.tmp) # assemble path
    sample.dist.tmp <- readRDS(fp.tmp) # read sampling dist file
    
    for(k in seq_along(var.v)) { # loop over variables - these need to match between ow.pts.dat and the sample.dist.tmp
      var.tmp <- var.v[k] # extract variable
      if(DBG) {cat("\n\t\t"); cat(crayon::bgBlue("k = ", k, "\t var = ", var.tmp))} # debug print
      sample.v <- ow.pts.dat[[var.tmp]]
      if(class(sample.v) == "factor") {sample.v <- as.numeric(sample.v) }
      dens.tmp <- sample.dist.tmp[[var.tmp]]
      #
      class.tmp <- class(dens.tmp)
      if(class.tmp == "histogram") {
        surp.v <- calc_surprisal_hist_(dens.tmp, sample.v)
      } else if (class.tmp == "table") {
        surp.v <- calc_surprisal_tab_(dens.tmp, sample.v)
      } else {
        cat(crayon::bgRed("\n\n\nSHOULDN'T REACH - calc_surprisal_.R\n\n\n"))
      }
      #
      n.dat <- length(surp.v)
      if(DBG) {cat("\n\t\t\t"); cat(crayon::bgWhite("adding to data..."))}
      return.dat <- return.dat |> 
        add_row(
          var = rep(var.tmp, n.dat),
          bw = rep(bw.tmp, n.dat),
          surprisal = surp.v
        )
    }
  } # end outer loop
  saveRDS(return.dat, file=paste0("clean_data/plot/information_plot_dat_", grd.int, ".Rds") )
}

