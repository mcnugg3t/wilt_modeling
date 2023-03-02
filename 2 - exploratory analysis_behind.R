
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
  dat.10 <- terra::rast("clean_data/joindat_10_800_50.tif")
  dat.30 <- terra::rast("clean_data/joindat_30_800_50.tif")
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
{ ##### 1 - Inhomogeneous L-function : estimate + gradient #### 
  rm(list=ls())
  gc()
  ow.ppp <- readRDS("clean_data/ow_ppp_10.Rds")
  dens.tst <- readRDS("clean_data/density/dens_bw_174.Rds")
  L.inhom <- Linhom(ow.ppp, dens.tst, correction="translation")
  L.inhom |> as.data.frame() |> as_tibble() |>
    mutate(theo = theo - r,
           trans = trans - r) |> 
    pivot_longer(cols=2:3, names_to="type", values_to="val") |>  
    ggplot(aes(x=r, y=val, color=type)) + geom_point() + xlim(0, 500) + ylim(-500, 100)
}
##
{ ##### 1 - Besag's L-function : estimate + gradient #### 
  rm(list=ls())
  gc()
  ow.ppp <- readRDS("clean_data/ow_ppp_10.Rds")
  #L.res <- Lest(ow.ppp, correction="translation")
  L.bootstrap <- lohboot(ow.ppp, Lest, correction="translation")
  plot(L.bootstrap)
  saveRDS(L.bootstrap, file="clean_data/exploratory/L_bootstrap.Rds")
  
  
  L.res |> plot(main = "Empirical L-function (translation correction) compared to theoretical")
  
  Lres.dat <- L.bootstrap |> 
    as.data.frame() |> 
    mutate(grad.L = (trans-lag(trans))/(r - lag(r)),
           grad.L.theo = (theo - lag(theo))/(r-lag(r)),
           grad.lo = (lo-lag(lo))/(r-lag(r)),
           grad.hi = (hi-lag(hi))/(r-lag(r))); Lres.dat
  
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
    select(r, grad.L, grad.lo, grad.hi) |> 
    pivot_longer(cols=2:4, names_to="type", values_to="val") |> 
    ggplot(aes(x = r, y = val, color=type)) + 
      geom_point() +
      geom_smooth(span=0.2) +
      xlim(0, 700) +
      ylim(0, 5) +
      labs(title="Gradient empirical L-function, smoothed with span=0.20", x="r", y="L(r)") +
      theme(plot.title = element_text(size=20))
}

##
{ ##### 2 - inhomogeneous L-function simulation envelope
  ow.ppp.10 <- readRDS("clean_data/ow_ppp_10.Rds")
  lam <- readRDS("clean_data/density/dens_bw_80.Rds")
  E <- envelope(ow.ppp.10, Linhom, sigma=124, correction=c("bord.modif", "translation"),
                simulate=expression(rpoispp(lam)),
                use.theory=T, nsim=3, global=T)
  plot(E, . - r ~ r)
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
    bw.res <- est_bandwidth_(kern.v, ow.ppp.30, verbose=T) # estimate as inhomogeneous poisson, estimate as cox
    saveRDS(bw.res, file="clean_data/bw_res_30.Rds")
    bw.res
  }
  
  # create bandwidth vector
  { 
    bw.v <- readRDS("clean_data/bw_res_30.Rds") |> unlist()
    range.bw <- seq(from=bw.v[1], to=bw.v[2], length.out = 7)
    inc <- range.bw[2]- range.bw[1]
    f1 <- range.bw[1]
    l1 <- range.bw[7]
    range.bw <- c(f1-inc, range.bw, l1+inc); range.bw
  }
  
  #
  # generate densities
  { 
    { # 10 m grid
      ow.ppp.10 <- readRDS("clean_data/ow_ppp_10.Rds")
      source("functions/generate_densities_.R")
      # wide search - grid 10
      generate_densities_(
        ow.ppp = ow.ppp.10,
        bw.v = range.bw,
        grd.int = 10
      )
      # fine search - grid 10
      generate_densities_(
        ow.ppp = ow.ppp.10,
        bw.v = c(65, 70, 80, 85, 90, 95, 105, 110),
        grd.int = 10
      )
    }
    
    { # 30 m grid
      ow.ppp.30 <- readRDS("clean_data/ow_ppp_30.Rds")
      # wide search - grid 30
      generate_densities_(
        ow.ppp = ow.ppp.30,
        bw.v = range.bw,
        grd.int = 30
      )
      # fine search - grid 30
      generate_densities_(
        ow.ppp = ow.ppp.30,
        bw.v = c(65, 70, 80, 85, 90, 95, 105, 110),
        grd.int = 30
      )
    }
  }
}

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
    
    continuous.vars.30 <- c("gw_30", "bd", "clay", "ph", "sand", "silt", 
                            "aspect", "channel_net_s5", "channel_net_s7", 
                            "conv_ind", "elev", "hillshade", "ls", "plan_curv",
                            "prof_curv", "rsp", "slope", "topo_wetness", "total_catch", 
                            "val_dep", "wl2_oakprob_30")
    source("functions/create_interact_terms_.R")
    interact.v <- create_interact_terms_(continuous.variables, verbose=T, DBG=F)
    interact.v.30 <- create_interact_terms_(continuous.vars.30, verbose=T, DBG=F)
    saveRDS(interact.v, file="clean_data/interact_v.Rds")
    saveRDS(interact.v.30, file="clean_data/interact_v_30.Rds")
    
    # load data and subset
    soils.df <- readRDS("clean_data/soils_df.Rds")
    dat.10 <- terra::rast("clean_data/joindat_10_800_50.tif")
    dat.in <- terra::subset(dat.10, subset=rm.vars, negate=T)
    
    dat.30 <- terra::rast("clean_data/joindat_30_800_50.tif")
    dat.in.30 <- terra::subset(dat.30, subset=rm.vars, negate=T)
    
    rm(dat.10, dat.30)
    gc()
    
    source("functions/create_sampling_distributions_.R")
  }
  
  { # 10 - m grid create sampling distributions with n.sim points
    create_sampling_distributions_(covar.rast = dat.in,
                                     join.dat = soils.df,
                                     n.sim=2.5e4,
                                     interact = T,
                                     inter.v = interact.v,
                                     rm.vars = rm.vars,
                                     grd.int = 10,
                                     verbose=T, 
                                     DBG=T)
    
    create_sampling_distributions_(covar.rast = dat.in.30,
                                   join.dat = soils.df,
                                   n.sim=2.5e4,
                                   interact = T,
                                   inter.v = interact.v.30,
                                   rm.vars = rm.vars,
                                   grd.int = 30,
                                   verbose=T, 
                                   DBG=T)
  }
}
##
{ ####### 7 - SURPRISAL DISTRIBUTIONS ########
  
  { # setup for surprisal
    rm(list=ls())
    gc()
    
    soils.dat.df <- readRDS("clean_data/soils_df.Rds")
    rm.vars <- readRDS("clean_data/rm_vars.Rds")
    interact.v <- readRDS("clean_data/interact_v.Rds")
    interact.v.30 <- readRDS("clean_data/interact_v_30.Rds")
    
    ow.pts.dat <- terra::rast("clean_data/joindat_10_800_50.tif") |> 
      as.data.frame() |>
      filter(ow_rast_10 > 0) |> 
      left_join(soils.dat.df) |> 
      select(-`study area`, -soils_rast, -ow_rast_10) |> 
      select(-all_of(rm.vars))
    
    ow.pts.dat.30 <- terra::rast("clean_data/joindat_30_800_50.tif") |> 
      as.data.frame() |>
      filter(ow_rast_30 > 0) |> 
      left_join(soils.dat.df) |> 
      select(-`study area`, -soils_rast, -ow_rast_30) |> 
      select(-all_of(rm.vars))
    
    source("functions/add_interact_.R")
    interact.dat <- add_interact_(ow.pts.dat, interact.v)
    interact.dat.30 <- add_interact_(ow.pts.dat.30, interact.v.30)
    
    rm(ow.pts.dat, ow.pts.dat.30)
    gc()
    
  }
  
  { ## CALC SURPRISAL
    dens.files <- list.files("clean_data/sample_dist/10/")
    dens.files.30 <- list.files("clean_data/sample_dist/30/")
    bw.select.v <- c(23, 49, 74, 99, 124, 149, 174, 199, 224, 249, 274) #**add here**
    source("functions/calc_surprisal_.R")
    calc_surprisal_(
      dens.files, 
      interact.dat, 
      bw.select.v, grd.int = 10,
      verbose=T, DBG=F)
    calc_surprisal_(
      dens.files.30, 
      interact.dat.30, 
      bw.select.v, grd.int = 30,
      verbose=T, DBG=F)
  }
  
  { ## PLOT SUPRISAL
    {
      plot.dat <- readRDS("clean_data/plot/information_plot_dat_10.Rds") 
      max.finite.surp <- max(plot.dat$surprisal[!is.infinite(plot.dat$surprisal)], na.rm=T); max.finite.surp
      plot.dat <- plot.dat |> 
        mutate(surprisal = if_else(
                              condition = is.infinite(surprisal), 
                              true = max.finite.surp + runif(n=1, min=0, max=4),
                              false = surprisal))
    }
    {
      plot.dat.30 <- readRDS("clean_data/plot/information_plot_dat_30.Rds") 
      max.finite.surp.30 <- max(plot.dat.30$surprisal[!is.infinite(plot.dat.30$surprisal)], na.rm=T); max.finite.surp.30
      plot.dat.30 <- plot.dat.30 |> 
        mutate(surprisal = if_else(
          condition = is.infinite(surprisal), 
          true = max.finite.surp.30 + runif(n=1, min=0, max=4),
          false = surprisal))
    }
    
    {
      # remove all interactions except top 10
      
      # filter to all vars and top 10 interaction terms (using median ic)
      filter_topN_interacts_ <- function(dat.in, n.top) {
        var.unique <- dat.in$var |> # unique variables
          unique()
        int.ind <- var.unique |>   # indices of those variables with interact _ x _
          str_detect(" x ")
        interact.vars <- var.unique[int.ind]; interact.vars # subset interact terms
        noninteract.vars <- var.unique[!int.ind]; noninteract.vars # and non-interact terms
        top.n.interact <- dat.in |> # get var names for top 10 interact terms
          filter(var %in% interact.vars) |> 
          group_by(var) |> 
          summarise(med_ic = median(surprisal, na.rm=T)) |> 
          arrange(desc(med_ic)) |> 
          select(var) |> 
          slice(1:n.top) |> 
          pull(var)
        return.dat <- dat.in |> # filter - either a plain covariate or one of the top 10 interaction terms
          filter(var %in% noninteract.vars | var %in% top.n.interact) |> 
          filter(!(var %in% c("SLOPE_CLAS", "SOIL_NAME", "EROSION_DI")))
        return(return.dat)
      }
      p1.dat <- filter_topN_interacts_(plot.dat, 20)
      p2.dat <- filter_topN_interacts_(plot.dat.30, 20)
      
      median.ic.val <- median(p1.dat$surprisal, na.rm=T)
      p1.dat |> 
        group_by(as.factor(bw), var) |> 
        ggplot(aes(y=surprisal, group=bw, fill=bw)) +
          geom_boxplot() + 
          facet_wrap(~var) + 
          geom_hline(yintercept=median.ic.val, color="red")
      
      median.ic.val <- median(p2.dat$surprisal, na.rm=T)
      p2.dat |> 
        group_by(as.factor(bw), var) |> 
        ggplot(aes(y=surprisal, group=bw, fill=bw)) +
          geom_boxplot() + 
          facet_wrap(~var) + 
          geom_hline(yintercept=median.ic.val, color="red")
      
      vis_1var_ <- function(in.dat, in.var) {
        return <- in.dat |> 
          filter(var == in.var) |> 
          group_by(as.factor(bw)) |> 
          ggplot(aes(y=surprisal, group=bw, fill=bw)) +
          geom_boxplot() + geom_hline(yintercept=median.ic.val, color="red")
      }
      
      vis_1var_(p1.dat, "slope") |> plot()
      vis_1var_(p1.dat, "channel_dist_s7") |> plot()
      vis_1var_(p1.dat, "gw_10") |> plot()
      
      vis_1var_(p2.dat, "clay") |> plot() # yes
      vis_1var_(p2.dat, "gw_30") |> plot() # yes/maybe
      vis_1var_(p2.dat, "rsp") |> plot() # yes
      vis_1var_(p2.dat, "val_dep") |> plot() # yes/maybe - more at higher BW
      vis_1var_(p2.dat, "wl2_oakprob_30") |> plot() # yes
    }
    
    
    p1.dat.smry <- p1.dat |> 
      group_by(var) |> 
      summarise(sd.ic = sd(surprisal, na.rm=T),
                med.ic = median(surprisal, na.rm=T))
    
    p2.dat.smry <- p2.dat |> 
      group_by(var) |> 
      summarise(sd.ic = sd(surprisal, na.rm=T),
                med.ic = median(surprisal, na.rm=T))
    
    p1.dat.smry |> 
      ggplot(aes(x=sd.ic, y=med.ic, color=var)) +
        geom_label(aes(label=var), size=5) +
        labs(title = "Variable Information Content: Median vs SD", x="Standard deviation of information content", y="Median information content") +
        theme(plot.title=element_text(size=30),
              axis.title = element_text(size=20))
    
    p2.dat.smry |> 
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





##
# { ##### 3 - from each density, create a list of point patterns ###### 
#   #
#   { # prep
#     rm(list=ls())
#     gc()
#     source("functions/create_ppps_.R")
#   }
#   create_ppps_(n.sim=100, verbose=T, DBG=T)
# }
##
# { ##### 4 - simulation envelope for L under varying bw values ########### 
#   { # setup
#     rm(list=ls())
#     gc()
#     dens.files <- list.files("clean_data/sample/")
#     ow.ppp <- readRDS("clean_data/ow_ppp_10.Rds")
#     source("functions/simulation_envelopes_L_.R")
#   }
#   simulation_envelopes_L_(dens.files, ow.ppp, verbose=T, DBG=T)
# }
##
