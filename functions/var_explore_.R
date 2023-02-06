require(spatstat)
require(terra)
require(tidyverse)
require(assertthat)
#'
#'
#'
var_explore_ <- function(rast.dat, pts.dat, kern.band, n.sim, verbose=T, DBG=T) {
  if(verbose) cat("\n\nVAR EXPLORE...")
  if(DBG) cat("\n\tsetup...")
  ### construct spatstat objects
  study.area <- rast.dat[["study area"]]
  sa.df <- study.area |> 
    as.data.frame(xy=T)
  # construct spatstat owin object 
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
  ow.ppp <- ppp(x = pts.crds[,1], y = pts.crds[,2], window = ow.win) # 59 points outside of window + some duplicated

  if(verbose) cat("\n\tcalc density...")
  t1 <- Sys.time()
  ow.dens <- density(ow.ppp, sigma=kern.band)
  if(verbose) cat(paste0("\ttime : ", difftime(Sys.time(), t1, units="secs")))
  
  ### kernel density -> tibble, normalize, resample onto study area -> data.frame
  if(DBG) cat("\n\t\tconverting, normalizing...")
  ow.adj <- ow.dens |> 
      as_tibble() |> 
      mutate(val.norm = value/sum(value)) |> 
      select(-value)
  
  ow.adj.rast <- rast(ow.adj, crs=crs(study.area) ) |> 
      resample(study.area)
  
  ow.adj.tbl <- ow.adj.rast |> 
      as.data.frame(xy=T) |> 
      mutate(val.norm = abs(val.norm))
  
  plot(ow.adj.rast[["val.norm"]])
  plot(pts.dat, cex=0.5, col=rgb(red=0, green=0, blue=0, alpha=0.3), add=T)
  
  rm(ow.adj)
  gc()
  
  check.sum <- ow.adj.tbl |> 
    select(val.norm) |> 
    sum()
  
  if(DBG) cat(paste0("\n\tdensity sums to 1: " , assert_that(abs(check.sum-1) < 0.0001) ))
  
  if(verbose) cat("\n\tsampling indices ~ kernel density....")
  # sample indices nsim times

  s.list <- list()
  for(i in 1:n.sim) {
    cat(paste0("\n\014index sampling: ", round(i/n.sim, 3)*100, " %"))
    s.tmp <- sample(
        x = 1:nrow(ow.adj.tbl), 
        size = rpois(n=1, lambda=nrow(pts.crds)), 
        replace=T,
        prob=ow.adj.tbl$val.norm)
    s.list[[i]] <- s.tmp
  }
  
  if(verbose) cat("\n\textracting covariates ~ sampled indices")
  # extract covariate values from data, compute mean and sd of each, store it
  # init empty tibble
  store.tbl <- tibble(
      var = character(),
      mu = numeric(),
      sig = numeric())
  
  source("functions/detect_remove_vars_.R")
  rmv.vars <- detect_remove_vars_(
    var.names = names(rast.dat),
    patterns = c("study area", "manage_rast", "wl2_cls", "soils_rast", "ow_rast"),
    verbose=F
  )
  
  # for each sample of indices...
  for(i in seq_along(s.list)) {
      cat(paste0("\n\014extracting covars: ", round(i/length(s.list), 3)*100, " %"))
      s.tmp <- s.list[[i]] # extract sample vector
      # get coords associated with indices, -> matrix
      crds.mat <- ow.adj.tbl |>  
        slice(s.tmp) |>  
        select(x, y) |> 
        as.matrix()
      
      # extract covar values from dat.10
      dat.subs <- terra::extract(rast.dat, crds.mat) |> 
        as_tibble() |> 
        select(-all_of(rmv.vars)) |> 
        as.matrix()
      
      # calc mean and sd of each variable
      mean.v <- apply(dat.subs, MARGIN=2, FUN=mean, na.rm=T)
      sd.v <- apply(dat.subs, MARGIN=2, FUN=sd, na.rm=T )
      
      # store in tbl
      store.tbl <- store.tbl |> 
        add_row(
          var = names(mean.v),
          mu = mean.v,
          sig = sd.v
        )
  }
  
  # manipulate store.tbl :
  #   1. add index for re-pivot
  #   2. pivot longer
  #   3. group and compute center + scale params
  store.tbl.2 <- store.tbl |> 
      add_column(ind = 1:nrow(store.tbl)) |> 
      pivot_longer(cols=2:3, names_to=c("stat"), values_to="val") |> 
      group_by(var, stat) |> 
      mutate(stat.center = mean(val, na.rm=T),
             stat.scale = sd(val, na.rm=T),
             val.scale = (val-stat.center)/stat.scale)
  
  # store scaling parameters
  scale.tbl <- store.tbl.2 |> 
      group_by(var, stat) |> 
      summarise(center = mean(stat.center),
                scale = mean(stat.scale))
    
  plot.tbl <- store.tbl.2 |> 
      select(-val, -stat.center, -stat.scale) |>
      pivot_wider(id_cols=c(var, ind), names_from=stat, values_from=val.scale)
    
  # calc observed covar mean and sd
  obs.dat <- terra::extract(
      x = rast.dat, 
      y = pts.crds)
    
  obs.samp <- obs.dat |> 
      as_tibble() |> 
      select(-all_of(rmv.vars))
    
  obs.mu <- apply(obs.samp, MARGIN=2, FUN=mean, na.rm=T)
  obs.sd <- apply(obs.samp, MARGIN=2, FUN=sd, na.rm=T)
    
  obs.tbl <- tibble(
      var = names(obs.samp),
      mu = obs.mu,
      sig = obs.sd) |> 
    pivot_longer(cols = 2:3, names_to = "stat", values_to="val") |> 
    left_join(scale.tbl, by=c("var", "stat")) |> 
    mutate(val.scale = (val - center)/scale ) |> 
    select(-val, -center, -scale) |> 
    pivot_wider(id_cols=1, names_from="stat", values_from=val.scale)
    # we want to scale the observed values onto the sampling dist - each variable has its own
    
  p1 <- ggplot() +
    geom_hex(data = plot.tbl, aes(x=mu, y=sig)) +
    geom_point(data=obs.tbl, aes(x=mu, y=sig), size=8, shape=10, color="red") +
    facet_wrap(~var) +
    theme(strip.text.x = element_text(size=18))
  
  return(list("sample.dat" = plot.tbl,
              "obs.dat" = obs.tbl,
              "plot" = p1))
}