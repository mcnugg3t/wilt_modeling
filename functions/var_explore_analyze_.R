require(spatstat)
require(terra)
require(tidyverse)
require(assertthat)
#'
#'
#'
var_explore_analyze_ <- function(rast.dat, pts.dat, kern.band, n.sim, verbose=T, DBG=T) {
  
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
  pts.crds <- crds(pts.dat)
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