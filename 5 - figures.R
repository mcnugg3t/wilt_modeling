
{
  rm(list=ls())
  gc()
  library(tidyverse)
  library(assertthat)
  library(plotly)
  library(terra)
}

{ ########### MENGES & LOUCKS 1984 #############
  x.v = seq(0, 500, by=1)
  y.v = 0.628*log(x.v, base=10) - 0.656
  plot.tbl <- tibble(
    x = x.v,
    y = y.v
  ) |> 
    mutate(grad = (lag(y)-y)/(lag(x)-x)  ) |>
    mutate(
      y = if_else(
        condition= is.infinite(y)|is.na(y),
        true=0,
        false=y),
      grad = if_else(condition= is.infinite(grad)|is.na(grad),
                  true=0,
                  false=grad),
      grad = grad/max(grad)
      ) |>
    filter(y > 0 & y <= 1) |> 
    mutate(grad = grad / (max(grad)*0.5) ) |> 
    pivot_longer(cols = 2:3, 
                 names_to = "type",
                 values_to = "val")
  
  plot.tbl |> 
    ggplot(aes(x=x,y=val, color=type)) + 
    geom_line(linewidth=2) + 
    ylim(0, 1) + 
    theme_minimal() + 
    labs(title="Cumulative density of new infections vs distance from old infections (d)", x="d", y="Cumulative Density") +
    theme(title=element_text(size=30), axis.text = element_text(size=24))
}

{ # plot KDE consistent with clustering pattern
  rm(list=ls())
  gc()
  library(spatstat)
  library(plotly)
  dens.dat.plt <- readRDS("clean_data/density/dens_bw___.Rds") |> 
    as.data.frame() 
  max.dens = max(dens.dat.plt$value, na.rm=T)
  inc.tmp = max.dens/10
  fig <- plot_ly(
    data = dens.dat.plt,
    type = "surface",
    contours = list(
      z = list(show=T, start=0, end=max.dens, size=inc.tmp)
    ),
    x = ~x, 
    y = ~y, 
    z = ~value)  |> add_surface()
  fig
}

{ ######### RATE OF INFECTION VS OAK CLS PROB #############
  rm.vars <- readRDS("clean_data/rm_vars.Rds")
  dat.cln <- terra::rast("clean_data/joindat_10_1200.tif") |> 
    subset(subset=rm.vars, negate=T)
  names(dat.cln)
  
  dat.plt <- terra::rast("clean_data/joindat_10_1200.tif") |> 
    subset(subset=c("ow_rast_10", "wl2_cls_10", "wl2_oakprob_10")) |> 
    as.data.frame() |> 
    filter(wl2_cls_10 == 4230) |> 
    mutate(tile10 = ntile(wl2_oakprob_10, 10)) |> 
    group_by(tile10) |> 
    summarise(rate_infect = sum(ow_rast_10, na.rm=T)/n())
  
  dat.plt |> 
    mutate(rate_infect = rate_infect * 100) |> 
    ggplot(aes(x=tile10, y=rate_infect)) +
    geom_col() +
    labs(title = "Oak wilt infection rate versus oak density decile (Wiscland2)", 
         x="Decile: oak classification probability (WiscLand2)",
         y="Infection rate (per km2)") +
    theme(
      plot.title = element_text(size=24),
      axis.text = element_text(size=18),
      axis.title = element_text(size=25))
}

{
  # read base data
  dat.plt <- terra::rast("clean_data/joindat_10_1200.tif") |> 
    subset(subset=c("alpha", "hb", "ksat", "lambda", "n", "om", "theta_r", "theta_s"), negate=T) |> 
    as.data.frame(xy=T)  
  
  { # base layer = sparse sample, where height = elev - 
    base.layer <- dat.plt |> 
      mutate(z = elev - bedrock_10/100) |>
      select(x, y, z) |> 
      slice_sample(prop=0.1) |> 
      add_column(color.assign=0) #|> 
      #add_column(alph = 0.0)
      
    # surface layer - yellow where sand is high, brown where sand is low
    soil.layer <- dat.plt |> 
      mutate(z = elev) |> 
      select(x, y, z, sand) |> 
      mutate(color.assign = ntile(sand, 4)) |> 
      select(-sand) |> 
      #add_column(alph = 0.1) |> 
      slice_sample(prop=0.2)
      
    comb.dat <- rbind(base.layer, soil.layer)
  }
  
  {
    p1 <- plot_ly(data=comb.dat, x=~x, y=~y, z=~z, color=~color.assign, 
                  colors = c("black", "brown", "orange", "goldenrod", "yellow") ) |> 
      add_markers() |> 
      layout()
    p1
  }
  
}
