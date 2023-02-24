
##
##### INIT #### 
{
  rm(list=ls())
  gc()
  library(spatstat)
  library(terra)
  library(tidyverse)
  library(assertthat)
} |> suppressPackageStartupMessages()


##
{ ##### 0 - construct spatstat ppp objects and save #### 
  dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
  dat.30 <- terra::rast("clean_data/joindat_30_1200.tif")
  ow.pts <- terra::vect("clean_data/ow_pts_clean.shp")
  #ow.pts <- terra::vect("mid_data/wilt/ow_pts_comb.shp")
  source("functions/construct_ppp_.R")
  ow.ppp.10 <- construct_ppp_(dat.10, ow.pts, T)
  ow.ppp.30 <- construct_ppp_(dat.30, ow.pts, T)
  rm(dat.10, dat.30, ow.pts)
  gc()
  saveRDS(ow.ppp.10, file="clean_data/ow_ppp_10.Rds")
  saveRDS(ow.ppp.30, file="clean_data/ow_ppp_30.Rds")
}
##
{ ##### 1 - Besag's L-function : estimate + gradient #### 
  ow.ppp <- readRDS("clean_data/ow_ppp_10.Rds")
  L.res <- Lest(ow.ppp, correction="translation")
  L.res |> plot(main = "Empirical L-function (translation correction) compared to theoretical")
  
  Lres.dat <- L.res |> 
    as.data.frame() |> 
    mutate(grad.L = (trans-lag(trans))/(r - lag(r)),
           grad.L.theo = (theo - lag(theo))/(r-lag(r))); Lres.dat
  
  # Lres.dat |> 
  #   filter(r < 2000) |> 
  #   select(r, grad.L, grad.L.theo) |> 
  #   pivot_longer(cols=2:3, names_to = "type", values_to="val")
  
  Lres.dat |> 
    filter(r < 2000) |> 
    select(r, trans, theo) |> 
    pivot_longer(cols=2:3, names_to="type", values_to="val") |> 
    ggplot(aes(x=r, y=val, color=type)) +
    geom_smooth()
  
  Lres.dat |> 
    filter(r < 2000) |> 
    select(r, grad.L) |> 
    ggplot(aes(x = r, y = grad.L)) + 
      geom_point() +
      geom_smooth(span=0.2) +
      xlim(0, 700) +
      labs(title="Gradient empirical L-function, smoothed with span=0.20", x="r", y="L(r)") +
      theme(plot.title = element_text(size=20))
}
##
{ ##### 2 - generate densities across bw values #### 
  #
  # estimate bandwidth using 30-m x 30-m quadrature scheme
  { 
    rm(list=ls())
    gc()
    kern.v <- c("bw.ppl", "bw.diggle")
    ow.ppp.30 <- readRDS("clean_data/ow_ppp_30.Rds")
    source("functions/est_bandwidth_.R")
    bw.res <- est_bandwidth_(kern.v, ow.ppp.30, verbose=T)
    saveRDS(bw.res, file="clean_data/bw_res_30.Rds")
  }
  
  # create bandwidth vector
  { 
    bw.v <- readRDS("clean_data/bw_res_30.Rds")
    range.bw <- seq(from=bw.v[1], to=bw.v[2], length.out = 7)
    inc <- range.bw[2]- range.bw[1]
    f1 <- range.bw[1]
    l1 <- range.bw[7]
    range.bw <- c(f1-2*inc, f1-inc, range.bw, l1+inc, l1+2*inc)
  }
  
  #
  # generate densities
  { 
    rm(list=ls())
    gc()
    ow.ppp.10 <- readRDS("clean_data/ow_ppp_10.Rds")
    source("functions/generate_densities_.R")
    # wide search
    generate_densities_(
      ow.ppp = ow.ppp.10,
      bw.v = range.bw
    )
    # fine search
    generate_densities_(
      ow.ppp = ow.ppp.10,
      bw.v = c(65, 70, 80, 85, 90, 95, 105, 110)
      )
  }
}
##
{ ##### 3 - from each density, create a list of point patterns ###### 
  #
  { # prep
    rm(list=ls())
    gc()
    source("functions/create_ppps_.R")
  }
  create_ppps_(n.sim=100, verbose=T, DBG=T)
}
##
{ ##### 4 - simulation envelope for L under varying bw values ########### 
  { # setup
    rm(list=ls())
    gc()
    dens.files <- list.files("clean_data/sample/")
    ow.ppp <- readRDS("clean_data/ow_ppp_10.Rds")
    source("functions/simulation_envelopes_L_.R")
  }
  simulation_envelopes_L_(dens.files, ow.ppp, verbose=T, DBG=T)
}
##
{ ####### 5 - Fry plot locally + globally ########
  
}
##
{ ####### 6 - COVARIATE SAMPLING DISTRIBUTIONS ########
  
  { # prep
    rm(list=ls())
    gc()
    
    # variables to remove
    rm.vars <- c("alpha", "hb", "ksat", "lambda", "n", "om", "theta_r", "theta_s")
    saveRDS(rm.vars, file="clean_data/rm_vars.Rds")
    
    # all possible interactions between continuous variables
    continuous.variables <- c("gw_10", "bd", "clay", "ph", "sand", "silt", 
                              "aspect", "channel_dist_s5", "channel_dist_s7", 
                              "conv_ind", "elev", "hillshade", "ls", "plan_curv",
                              "prof_curv", "rsp", "slope", "topo_wet", "total_catch", 
                              "val_depth", "wl2_oakprob_10")
    source("functions/create_interact_terms_.R")
    interact.v <- create_interact_terms_(continuous.variables, verbose=T, DBG=F)
    saveRDS(interact.v, file="clean_data/interact_v.Rds")
    
    soils.df <- readRDS("clean_data/soils_df.Rds")
    dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
    dat.in <- terra::subset(dat.10, subset=rm.vars, negate=T)
    rm(dat.10)
    gc()
    source("functions/create_sampling_distributions_.R")
  }
  create_sampling_distributions_(covar.rast = dat.in,
                                   join.dat = soils.df,
                                   n.sim=5e4,
                                   interact = T,
                                   interact.v = interact.v,
                                   rm.vars = rm.vars,
                                   verbose=T, 
                                   DBG=T)
}
##
{ ####### 7 - SURPRISAL DISTRIBUTIONS ########
  
  { # setup for surprisal
    rm(list=ls())
    gc()
    
    soils.dat.df <- readRDS("clean_data/soils_df.Rds")
    rm.vars <- readRDS("clean_data/rm_vars.Rds")
    interact.v <- readRDS("clean_data/interact_v.Rds")
    
    ow.pts.dat <- terra::rast("clean_data/joindat_10_1200.tif") |> 
      as.data.frame() |>
      filter(ow_rast_10 > 0) |> 
      left_join(soils.dat.df) |> 
      select(-`study area`, -soils_rast, -ow_rast_10) |> 
      select(-all_of(rm.vars))
    
    source("functions/add_interact_.R")
    interact.dat <- add_interact_(ow.pts.dat, interact.v)
    var.v <- names(interact.dat); var.v
    rm(ow.pts.dat)
    gc()
    
    return.dat <- tibble(
      var = character(),
      bw = numeric(),
      surprisal = numeric()
    )
  }
  
  { ## CALC SURPRISAL
    dens.files <- list.files("clean_data/sample_dist")
    bw.select.v <- c(70, 74, 80, 85, 90, 95, 99, 105, 110)
    source("functions/calc_surprisal_.R")
    calc_surprisal_(dens.files, interact.dat, bw.select.v, verbose=T, DBG=F)    
  }
  
  { ## PLOT SUPRISAL
    plot.dat <- readRDS("clean_data/plot/information_plot_dat.Rds") |> 
      mutate(surprisal = if_else(is.infinite(surprisal), max(surprisal)+runif(n=1, min=0, max=4), surprisal))
    
    tst.v = plot.dat |> 
      filter(var == "ph x hillshade") |> 
      pull(surprisal)
    
    # for a given variable, what are median and SD of surprisal across all bw values
    p1.dat <- plot.dat |> 
      group_by(var) |> 
      summarise(med.ic = median(surprisal, na.rm=T),
                sd.ic = sd(surprisal, na.rm=T)) |> 
      arrange(desc(med.ic)) |> 
      mutate(rank = row_number())
    
    # remove all interactions except top 10
    is.interact <- p1.dat$var |> str_detect(pattern=" x ") |> as.numeric(); is.interact
    subset.rows = as.numeric(!is.interact); subset.rows
    subset.rows[1:10] = 1; 
    
    
    p2.dat = p1.dat[which(subset.rows == 1), ] 
    
    p2.dat |> 
      ggplot(aes(x=sd.ic, y=med.ic, color=var)) +
        geom_label(aes(label=var), size=5) +
        labs(title = "Variable Information Content: Median vs SD", x="Standard deviation of information content", y="Median information content") +
        theme(plot.title=element_text(size=30),
              axis.title = element_text(size=20))
    
    p1.dat |> 
      filter(rank <= 70) |> 
      ggplot(aes(x=med.ic, y=sd.ic, size=rank, color=var)) + geom_point() + geom_label(aes(x=med.ic, y=sd.ic, label=var))
    
    p2.dat <- plot.dat |> 
      group_by(var, bw) |> 
      summarise(med.ic = median(surprisal)) |> 
      group_by(var) |> 
      summarise(max.med.ic = max(med.ic)) |> 
      arrange(desc(max.med.ic)) |> 
      mutate(index = row_number())
    p2.dat |> ggplot(aes(x=index, y=max.med.ic)) + geom_col()
    source("functions/plot_surprisal_.R")
    plot_surprisal_(plot.dat)
  }
}

