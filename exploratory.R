
{
  rm(list=ls())
  gc()
}

{
  library(spatstat)
  library(terra)
  library(tidyverse)
  library(assertthat)
} |> suppressPackageStartupMessages()

# load data (1) all covariates in raster stack: 10 m x 10 m grid with 1200 m buffer 
#           (2) shapefile of all oak wilt points
{
  #dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
  dat.tmp <- terra::rast("clean_data/joindat_30_1200.tif")
  ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp")
}

# construct spatstat objects
{ 
  #study.area <- dat.tmp[["study area"]]
  study.area <- dat.tmp[["study area"]]
  # convert study area raster object to dataframe with coordinates
  sa.df <- study.area |> 
    as.data.frame(xy=T)
  # construct spatstat owin object 
  ext.sa <- ext(study.area)
  ow.win <- owin(
    xrange = c(ext.sa[1], ext.sa[2]),
    yrange = c(ext.sa[3], ext.sa[4]),
    mask = sa.df[,1:2]
  )
  # extract coordinates
  ow.crds <- crds(ow.pts)
  # construct spatstat ppp object
  ow.ppp <- ppp(x = ow.crds[,1], y = ow.crds[,2], window = ow.win) # 59 points outside of window + some duplicated
}
#
# nrow(sa.df) # 676.5k rows

# construct density object
{
  ow.dens <- density(ow.ppp, sigma=150)
}

{
  x.v <- seq(-300, 300, by=1)
  y.v <- dnorm(x.v, mean=0, sd=150)
  plot(x.v, y.v)
}

# plot
plot(ow.dens)
plot(ow.dens^(2/3))
plot(ow.dens^(1/2))

# kernel density -> tibble, normalize, resample onto study area -> data.frame
{
ow.adj <- ow.dens^(2/3) |> 
  as_tibble() |> 
  mutate(val.norm = value/sum(value)) |> 
  select(-value)
ow.adj.rast <- rast(ow.adj, crs=crs(study.area) ) |> 
  resample(study.area)
ow.adj.tbl <- ow.adj.rast |> 
  as.data.frame(xy=T)
rm(ow.adj)
gc()
}

# check density estimate
{
tst.sum <- values(ow.adj.rast)[!is.na(values(ow.adj.rast))] |> sum()
assert_that(abs(tst.sum-1) < 0.0001)
plot(ow.adj.rast[["val.norm"]])
plot(ow.pts, cex=0.5, col=rgb(red=0, green=0, blue=0, alpha=0.3),  add=T)
#rm(ow.adj.rast)
gc()
}

# sample indices nsim times
{
  n.sim <- 1000
  s.list <- list()
  for(i in 1:n.sim) {
    cat(paste0("\n\014", round(i/n.sim, 3)*100, " %"))
    s.tmp <- sample(
      x = 1:nrow(ow.adj.tbl), 
      size = rpois(n=1, lambda=nrow(ow.crds)), 
      replace=T,
      prob=ow.adj.tbl$val.norm
    )
    s.list[[i]] <- s.tmp
  }
}

# extract covariate values from data, compute mean and sd of each, store it
{
  # init empty tibble
  store.tbl <- tibble(
    var = character(),
    mu = numeric(),
    sig = numeric()
  )
  # for each sample of indices
  for(i in seq_along(s.list)) {
    cat(paste0("\n\014extracting covariates : ", round(i/length(s.list), 3)*100, " %"))
    s.tmp <- s.list[[i]] # extract sample vector
    crds.mat <- ow.adj.tbl |> # get coords associated 
      slice(s.tmp) |>  
      select(x, y) |> 
      as.matrix()
  
  # extract values from dat.10
  dat.subs <- terra::extract(dat.tmp, crds.mat) |> 
    as_tibble() |> 
    select(-`study area`, -manage_rast_30, -wl2_cls_30) |> 
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
}

# manipulate store.tbl :
#   1. add index for re-pivot
#   2. pivot longer
#   3. group and compute center + scale params
{
  tbl.tmp <- store.tbl |>
    filter(var != "ow_rast_30") 
  store.tbl.2 <- tbl.tmp |> 
    add_column(ind = 1:nrow(tbl.tmp)) |> 
    pivot_longer(cols=2:3, names_to=c("stat"), values_to="val") |> 
    group_by(var, stat) |> 
    mutate(stat.center = mean(val, na.rm=T),
           stat.scale = sd(val, na.rm=T),
           val.scale = (val-stat.center)/stat.scale
           )
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
    x = dat.tmp, 
    y = ow.crds
  )

  obs.samp <- obs.dat |> 
    as_tibble() |> 
    select(-`study area`, -manage_rast_30, -wl2_cls_30, -ow_rast_30)

  obs.mu <- apply(obs.samp, MARGIN=2, FUN=mean, na.rm=T)
  obs.sd <- apply(obs.samp, MARGIN=2, FUN=sd, na.rm=T)

  obs.tbl <- tibble(
    var = names(obs.samp),
    mu = obs.mu,
    sig = obs.sd
    ) |> 
    pivot_longer(cols = 2:3, names_to = "stat", values_to="val") |> 
    left_join(scale.tbl, by=c("var", "stat")) |> 
    mutate(val.scale = (val - center)/scale ) |> 
    select(-val, -center, -scale) |> 
    pivot_wider(id_cols=1, names_from="stat", values_from=val.scale)
  # we want to scale the observed values onto the sampling dist - each variable has its own

  ggplot() +
    geom_hex(data = plot.tbl, aes(x=mu, y=sig)) +
    geom_point(data=obs.tbl, aes(x=mu, y=sig), size=8, shape=10, color="red") +
    facet_wrap(~var) +
    theme(strip.text.x = element_text(size=18))
}
  
# if you do any of the Neyman-Scott models with different kernels, can you recover anything like the observed figures?
# create a quadscheme?
