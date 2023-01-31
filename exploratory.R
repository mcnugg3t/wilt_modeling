
rm(list=ls())
gc()

{
  library(spatstat)
  library(terra)
  library(tidyverse)
} |> suppressPackageStartupMessages()
# load data (1) all covariates in raster stack: 10 m x 10 m grid with 1200 m buffer 
#           (2) shapefile of all oak wilt points
dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp")
# study area bit-mask is one of the layers in the raster stack
study.area <- dat.10[["study area"]]
# convert terra::raster object to dataframe with coordinates
sa.df <- study.area |> 
  as.data.frame(xy=T)

# construct spatstat owin object 
ext.sa <- ext(study.area)
ow.win.10 <- owin(
  xrange = c(ext.sa[1], ext.sa[2]),
  yrange = c(ext.sa[3], ext.sa[4]),
  mask = sa.df[,1:2]
  )

# construct spatstat ppp object
ow.crds <- crds(ow.pts)
ow.ppp <- ppp(x = ow.crds[,1], y = ow.crds[,2], window = ow.win.10) # 59 points outside of window + some duplicated

ow.ppm <- ppm(ow.ppp, )

ow.dens <- density(ow.ppp, sigma=300)
ow.thin.dens <- density(rthin(ow.ppp, P = 0.2), sigma=300)
log.dens <- log(ow.dens)
plot(ow.dens)
plot(ow.dens^(2/3))
plot((ow.dens))
plot(log.dens)
{
plot(ow.thin.dens)
plot(sqrt(ow.thin.dens))
plot(log(ow.thin.dens))
}

{
ow.adj <- ow.dens^(2/3) |> 
  as_tibble() |> 
  mutate(val.norm = value/sum(value)) |> 
  select(-value)
ow.adj.rast <- rast(ow.adj, crs=crs(study.area) ) |> 
  resample(study.area)
ow.adj <- ow.adj.rast |> 
  as.data.frame(xy=T)
gc()
}

values(ow.adj.rast)[!is.na(values(ow.adj.rast))] |> sum()

plot(ow.adj.rast[["val.norm"]])
plot(ow.pts, add=T)
rm(ow.adj.rast)
gc()

# sample indices
n.samp <- 1000
s.list <- list()
for(i in 1:n.samp) {
  s.tmp <- sample(
    x = 1:nrow(ow.adj), 
    size = rpois(n=1, lambda=nrow(ow.crds)), 
    replace=T,
    prob=ow.adj$val.norm
  )
  s.list[[i]] <- s.tmp
}

#tst |> ggplot(aes(x=value)) + geom_density()

store.tbl <- tibble(
  var = character(),
  mu = numeric(),
  sig = numeric()
)

# for each sample of indices
for(i in seq_along(s.list)) {
  s.tmp <- s.list[[i]] # extract sample vector
  crds.mat <- ow.adj |> # get coords associated 
    slice(s.tmp) |>  
    select(x, y) |> 
    as.matrix()
  
  # extract values from dat.10
  dat.subs <- terra::extract(dat.10, crds.mat) |> 
    as_tibble() |> 
    select(-`study area`, -manage_rast_10, -wl2_cls_10)
  
  # mean and sd of each variable
  mean.v <- apply(dat.subs, MARGIN=2, FUN=mean, na.rm=T)
  sd.v <- apply(dat.subs, MARGIN=2, FUN=sd, na.rm=T )
  # store in vector / matrix
  store.tbl <- store.tbl |> 
    add_row(
      var = names(mean.v),
      mu = mean.v,
      sig = sd.v
    )
}

store.tbl.2 <- store.tbl |>
  add_column(ind = 1:nrow(store.tbl)) |> 
  pivot_longer(cols=2:3, names_to=c("stat"), values_to="val") |> 
  group_by(var, stat) |> 
  mutate(stat.center = mean(val, na.rm=T),
         stat.scale = sd(val, na.rm=T),
         val.scale = (val-stat.center)/stat.scale
         )

scale.tbl <- store.tbl.2 |> 
  group_by(var, stat) |> 
  summarise(center = mean(stat.center),
            scale = mean(stat.scale))

plot.tbl <- store.tbl.2 |> 
  select(-val, -stat.center, -stat.scale) |>
  pivot_wider(id_cols=c(var, ind), names_from=stat, values_from=val.scale)

# now we want the mean and sd observed
obs.dat <- terra::extract(
  x = dat.10, 
  y = ow.crds
  )

obs.samp <- obs.dat |> 
  as_tibble() |> 
  select(-`study area`, -manage_rast_10, -wl2_cls_10)

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

store.tbl |> 
  ggplot(aes(x=mu, y=sig)) +
  geom_hex() +
  geom_point(data=obs.tbl, aes(x=mu, y=sig)) +
  facet_wrap(~var)

# if you do any of the Neyman-Scott models with different kernels, can you recover anything like the observed figures?

# create a quadscheme
