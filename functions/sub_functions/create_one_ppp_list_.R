require(terra)
require(tidyverse)
require(crayon)
#' this function generates list of n.sim random point patterns for a given density
#' 0) takes density object as df, ow.rast, n.sim 1) determines number of cells ~ poisson(500) and their locations ~ density
#'
create_one_ppp_list_ <- function(density.df, ow.df, n.sim, owin.in, n.ow.cells, verbose=T, DBG=F) {
  if(verbose) cat(crayon::bgGreen("\n\tFUNCT : create_ppps_"))
  
  if(verbose) cat( paste0("\n\tdrawing n.cells infected ~ Poisson(", n.ow.cells, ")") )
  
  n.pts = rpois(
    n=n.sim, 
    lambda=n.ow.cells)
  
  ind.mat <- matrix(
    data=0, 
    nrow=n.sim, 
    ncol=max(n.pts))
  
  if(verbose) cat(paste0("\n\tassigning cell totals to ", n.sim, " vectors with avg. length ", 
                         mean(n.pts), "  ~ density (takes time)\n\t\t1 -> 10  :  ") )
  
  prog.inc <- n.sim %/% 10 
  for(j in 1:n.sim) {
    if(j %% prog.inc == 0) cat( paste0( (j %/% prog.inc), ".." ) )
    samp.tmp <- sample(
      x = 1:nrow(density.df),
      size = n.pts[j], 
      replace = T, 
      prob = density.df$value)
    ind.mat[j,1:length(samp.tmp)] <- samp.tmp
  }
  
  ow.ct.tab <- ow.df$ow_rast_10 |> table()
  ow.ct.sum <- ow.ct.tab |> sum()
  ct.probs <- c(ow.ct.tab[1]/ow.ct.sum, ow.ct.tab[2]/ow.ct.sum, ow.ct.tab[3]/ow.ct.sum)
  
  
  ##
  ## BIG LOOP
  ##
  if(verbose) cat(paste0("\n\tcreating ppps (n=", n.sim  ,"), (takes time)\n\t\t2 -> 100  :  ") )
  ppp.list.return <- list()
  prog.inc <- n.sim %/% 50
  for(i in 1:n.sim) {
    if(i %% prog.inc == 0) cat( paste0( (i %/% prog.inc)*2, ".." ) )
    ind.tmp <- ind.mat[i,][ind.mat[i,] > 0] # get indices from matrix
    assert_that(length(ind.tmp) == n.pts[i]) # check
    
    ##
    ## create tibble with appropriate coordinates, add wilt counts within cells - empirical dist
    ind.tbl <- ind.tmp |>
      as_tibble() |>
      rename(index = value) |>
      add_column(
        x = density.df$x[ind.tmp],
        y = density.df$y[ind.tmp])
    ind.tbl <- ind.tbl |>
      add_column(
        count = sample(
          x=c(1,2,3),
          size=nrow(ind.tbl),
          replace=T,
          prob = ct.probs
        ))
    
    ##
    ## unfold ind.ext
    ind.ext <- ind.tbl
    while(sum(ind.ext$count > 1) > 0) {
      subs.rows <- ind.ext[ind.ext$count > 1,]
      ind.ext <- rbind(ind.ext, subs.rows)
      new.counts <- ind.ext$count[ind.ext$count > 1] - 1
      ind.ext$count[ind.ext$count > 1] <- new.counts
    }
    rm(ind.tbl)
    gc()
    # perturb each point
    coords.df <- ind.ext |> 
      select(x, y) |> 
      mutate(
        x = x+runif(n=nrow(ind.ext),min=-5,max=5),
        y = y+runif(n=nrow(ind.ext),min=-5,max=5)
      )
    
    # create ppp
    ppp.tmp <- ppp(
      x = coords.df[["x"]], 
      y = coords.df[["y"]], 
      window = owin.in)
    
    ppp.list.return[[i]] <- ppp.tmp
  } 
  ##
  ## END BIG LOOP
  return(ppp.list.return)
}
  
