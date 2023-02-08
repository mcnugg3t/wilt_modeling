############ FEB 9 ##########

{
  rm(list=ls())
  gc()
}

{
  library(spatstat)
  library(tidyverse)
  library(assertthat)
  library(terra)
  library(crayon)
} |> suppressPackageStartupMessages()

list.files("clean_data/sample/bw.diggle/")

# e.g. for bd
{ # setup
  tst.sample <- readRDS("clean_data/sample/bw.diggle/sample_dens_adj_1.Rds") |> 
    select(count, bd)
  n.expand <- tst.sample |> filter(count > 1) |> nrow() # number of rows with count >1
  count.v <- tst.sample$count[1:n.expand] # subset relevant counts
  val.v <- tst.sample$bd[1:n.expand] # subset relevant values
  val.minus.v <- tst.sample$bd[-(1:n.expand)] # subset values not needing to be repeated
}

{
  t1 <- Sys.time()
  for(i in seq_along(count.v)) {
    cat(paste0("\014", round(i/length(count.v), 3)*100, " %"))
    ct.tmp <- count.v[i]
    val.tmp <- val.v[i]
    tst.sample <- tst.sample |> 
      add_row(
        count = rep(1, ct.tmp),
        !!sym("bd") := rep(val.tmp, ct.tmp)
        )
  }
  t2 <- Sys.time()
  cat(crayon::bgRed("\n\n\ntime = ", difftime(t2, t1, units="secs")) )
}

{
  t1 <- Sys.time()
  covar.v <- c()
  for(i in seq_along(count.v)) {
    cat(paste0("\014", round(i/length(count.v), 3)*100, " %"))
    ct.tmp <- count.v[i]
    val.tmp <- val.v[i]
    covar.v <- c(covar.v, rep(val.tmp, ct.tmp))
  }
  covar.v <- c(covar.v, val.minus.v)
  t2 <- Sys.time()
  cat(crayon::bgRed("\n\n\ntime = ", difftime(t2, t1, units="secs")) )
}

###
###
###

# read kde
tst.dens <- readRDS("clean_data/density/bw.ppl/dens_adj_0.6.Rds")
# check kde
tst.dens.df <- tst.dens |> 
  as.data.frame()
n.neg <- tst.dens.df |> 
  filter(value < 0) |> 
  nrow()
# check no negative values
assert_that(n.neg == 0)

dens.df.scaled <- tst.dens.df |> 
  mutate(val_scale = value/sum(value, na.rm=T))
# check valid prob dist - sums to 1 after scaling
assert_that((1-sum(dens.df.scaled$val_scale))^2 < 0.1e-6) 
rm(tst.dens.df, tst.dens)
gc()

# sample from kde n.sim times
n.sim <- 1e5
ind.sample <- sample(1:nrow(dens.df.scaled), size=n.sim, replace=T, prob=dens.df.scaled$val_scale)
ind.unique <- ind.sample |> 
  as_tibble() |>
  rename(index = value) |> 
  group_by(index) |> 
  summarise(count = n()) |> 
  arrange(desc(count)) |> 
  add_column(x = dens.df.scaled$x[ind.unique$index],
             y = dens.df.scaled$y[ind.unique$index])

# save here?

ind.extr <- ind.unique |> 
  select(x, y) |> 
  as.matrix()

dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")

extr.res <- terra::extract(dat.10, ind.extr)
extr.res.full <- cbind(ind.unique, extr.res)
rm(extr.res)
gc()
# compute sampling distribution of a given covariate
# say, bd
var.tmp <- "bd"
count.v.tmp <- extr.res.full[["count"]]
var.v.tmp <- extr.res.full[[var.tmp]]
covar.v <- c()
for(i in seq_along(count.v.tmp) ) {
  cat(paste0("\014\n ", round(i/length(count.v.tmp), 3)*100, " %"))
  covar.v <- c(covar.v, rep(var.v.tmp[i], count.v.tmp[i]) )
}

var.hist <- covar.v |> 
  hist(breaks=500)


# what we have : some observed values
# what we want : the suprisal of each value (for a given bw) 
# then calculated across various bw values
# eventually : merging both 

# for a given covariate value, and a given bandwidth, we want to be able to go look up its surprisal

{
  tst.dat.clean <- tst.dens.df |> 
    mutate(val_corr = if_else(value < 0.0000001, 0, value) ) |>
    mutate(val_corr = val_corr / sum(val_corr))
  tst.dat.clean |> filter(val_corr > 0) |> nrow()
}
  mutate(val_corr = val_corr / sum(val_corr, na.rm=T))


############ FEB 7 ##########
{
  x.v <- runif(100000, 0, 1)
  y.v <- runif(100000, 0, 40)
  int.v <- x.v*y.v
  dat.viz <- tibble(
    x = x.v,
    y = y.v,
    interact = int.v
  )
  dat.viz |> 
    ggplot(aes(x=x, y=y, color=interact ) ) + 
    geom_point()
}
############ FEB 6 ##########
{
  rm(list=ls())
  gc()
  setwd("D:/Backed Up/Desktop/wilt_modeling")
}

{
  library(terra)
  library(tidyverse)
  library(assertthat)
} |> suppressPackageStartupMessages()

# join soil information
{
  dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
  soil.rast <- dat.10[["soils_rast"]]
  soil.rast.df <- soil.rast |> 
    as.data.frame(xy=T) |> 
    rename(soil.id = soils_rast)
  rm(dat.10, soil.rast)
  gc()
  soil.dat <- terra::vect("raw_data/cnnf_soil/CNNF_Soil_Layer_2.shp")
  soil.df <- soil.dat |> 
    as.data.frame() |> 
    select(-LTA, -LT, -LTP, -ROADS_ANAL, -BIOMASS_HA, -SHAPE_AREA, -SHAPE_LEN) |> 
    rename(soil.id = OBJECTID)
  soil.join.dat <- left_join(soil.rast.df, 
                             soil.df, 
                             by="soil.id")
  rm(soil.dat, soil.id, soil.rast.df)
  gc()
}

# analyze soil data
{
  plot.dat <- soil.join.dat |> 
    mutate(
      soil.name = as.factor(SOIL_NAME),
      surface.texture = as.factor(SURFACE_TE),
      drainage.class = as.factor(DRAINAGE_C),
      compaction = as.factor(COMPACTION),
      erosion = as.factor(EROSION_DI),
      hydric = as.factor(HYDRIC_SOI),
      equipment = as.factor(EQUIPMENT_)
      ) |> 
    select(-SOIL_NAME, -SURFACE_TE, -DRAINAGE_C, -EROSION_DI, -HYDRIC_SOI, -EQUIPMENT_)
  
  soil.join.dat |> 
    select()
  
  smry.dat <- plot.dat |> 
    group_by(soil.name) |> 
    summarise(rel_area = n()/nrow(plot.dat)*100) |> 
    arrange(desc(rel_area))
  
  smry.dat |> 
    ggplot() +
    geom_col(aes(x=soil.name, y=rel_area, fill=soil.name)) +
    theme(axis.text.x = element_text(angle=90))
}

