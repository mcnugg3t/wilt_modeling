#--------------- Setup--------------
{
  rm(list=ls())
  gc()
  library(tidyverse)
  library(assertthat)
  library(plotly)
  library(terra)
}
#--------------- Fig 1--------------

{
  dat.rast = terra::rast("clean_data/joindat_10_800_50.tif")
  one.lyr = dat.rast$`study area`
  #sa.poly = one.lyr |> terra::as.polygons() |> terra::
  sa.poly
  terra::writeVector(sa.poly, file="plt_data/sa_poly.shp")
}

#--------------- Fig 2 - Loh boot --------------
{
  L.res = readRDS("clean_data/exploratory/L_bootstrap.Rds")
  Lres.dat <- L.res |> 
    as.data.frame() |> 
    mutate(grad.L = (trans-lag(trans))/(r - lag(r)),
           grad.L.theo = (theo - lag(theo))/(r-lag(r)),
           grad.lo = (lo-lag(lo))/(r-lag(r)),
           grad.hi = (hi-lag(hi))/(r-lag(r))) 
  Lres.dat  
  
  Lres.dat |> 
    filter(r < 1000) |> 
    ggplot() + 
    geom_ribbon(aes(x=r, ymin=lo, ymax=hi), fill="gray") +
    geom_line(aes(x=r, y=trans), color="black", linewidth=1.5) +
    geom_line(aes(x=r, y=theo), color="red", linetype="dashed") +
    labs(x="r", y="L(r)") +
    theme(axis.title = element_text(size=32), axis.text = element_text(size=24))
  
  Lres.dat |> 
    filter(r < 1000) |> 
    ggplot() +
    geom_point( aes(x=r, y=grad.L), alpha=0.1 ) +
    geom_smooth( aes(x=r, y=grad.L), method="loess", span=0.35, linewidth = 1.5 ) +
    geom_line(aes(x=r, y=grad.L.theo), linewidth = 1.2, color="red", linetype="dashed") +
    labs(x="r", y="L'(r)") +
    ylim(0, 5) +
    theme(axis.title = element_text(size=32), axis.text = element_text(size=24))
  library(spatstat)
  ow.ppp = readRDS("clean_data/ow_ppp_10.Rds")
  fryplot(ow.ppp, dmax=1000, axes=T)
  
}
#--------------- Fig 3--------------
{
  compare.dat = readRDS("plt_data/model_compare.Rds")
  
  cpm.mods = unique(compare.dat$mod_name)[
    compare.dat$mod_name |> unique() |> stringr::str_detect("kppm")
  ]; cpm.mods
  
  compare.dat |> 
    filter(mod_name %in% cpm.mods) |> 
    filter(mod_name != "kppm_30_ha2_thom_pq") |> 
    pivot_wider(id_cols=1:2, names_from=measure, values_from=val) |>
    ggplot(aes(x=area_surveyed, y=recall, color=mod_name)) + 
    geom_line(linewidth=1.2) + 
    scale_color_brewer(type="div", palette="RdBu") +
    labs(x="Proportion of study area", y="Recall") +
    theme(axis.title = element_text(size=32),
          axis.text = element_text(size=24),
          legend.title = element_text(size=32),
          legend.text = element_text(size=24))
  
  compare.dat$mod_name |> unique()
  
  compare.res = compare.dat |> 
    mutate(Model =   if_else(mod_name == "gam_10_h0_logit", "GAM H0 10-m",
                     if_else(mod_name == "gam_10_ha1_logit", "GAM HA1 10-m",
                     if_else(mod_name == "gam_10_ha2_logit", "GAM HA2 10-m",
                     if_else(mod_name == "gam_30_h0_logit", "GAM H0 30-m",
                     if_else(mod_name == "gam_30_ha1_logit", "GAM HA1 30-m",
                     if_else(mod_name == "gam_30_ha2_logit", "GAM HA2 30-m",
                     if_else(mod_name == "kppm_30_cauch_pq", "CPM H0 30-m",
                     if_else(mod_name == "kppm_30_ha2_cauch_wc", "CPM HA1 30-m", "REMOVED"))))))))) |> 
                     #if_else(mod_name == "kppm_30_ha2_thom_wc", "CPM HA1 30-m", "REMOVED")))))) |> 
    select(-mod_name) |> 
    relocate(Model) |> 
    filter(Model != "REMOVED") |> 
    select(-scale_thresh)
  
  names.tmp = compare.res$Model |> unique(); names.tmp
  
  compare.res |> str()
  
  # FIG 3 C
  compare.res |> # compare kppm models with ha2 models
    filter(Model %in% c("GAM HA2 30-m", "CPM HA1 30-m", "CPM H0 30-m")) |> 
    pivot_wider(id_cols=1:2, names_from=measure, values_from=val) |>
    ggplot(aes(x=area_surveyed, y=recall, color=Model)) + 
    geom_line(linewidth=1.2) + 
    scale_color_brewer(type="div", palette="RdBu") +
    labs(x="Proportion of study area", y="Recall") +
    theme(axis.title = element_text(size=32),
          axis.text = element_text(size=24),
          legend.title = element_text(size=32),
          legend.text = element_text(size=24))
  
  gam.names = names.tmp[str_detect(names.tmp, pattern="GAM")]; gam.names
  
  # FIG 3A
  compare.res |> # compare GAMs to one another
    filter(Model %in% c(gam.names, "CPM HA1 30-m") ) |> 
    filter(str_detect(Model, "10-m", negate=T)) |> 
    mutate(Model = str_remove(Model, " 30-m")) |> 
    #filter(!(Model %in% c("GAM HA1 10-m", "GAM HA1 30-m") )) |> 
    #filter((mod_name %in% c("gam_10_h0_logit", "gam_30_h0_logit", "gam_10_ha1_logit", "gam_10_ha2_logit", "gam_30_ha2_logit", "gam_30_ha1_logit"))) |> 
    #filter((mod_name %in% c("gam_10_h0_logit", "gam_30_h0_logit", "gam_30_ha1_logit", "gam_10_ha1_logit", "gam_10_ha2_logit"))) |> 
    pivot_wider(id_cols=1:2, names_from=measure, values_from=val) |>
    ggplot(aes(x=area_surveyed, y=recall, color=Model)) + 
    geom_line(linewidth=1.2) + 
    scale_colour_viridis_d(option = "A") +
    #scale_color_brewer(type="qual", palette="RdBu") +
    labs(x="Proportion of study area", y="Recall") +
    theme(axis.title = element_text(size=32),
          axis.text = element_text(size=24),
          legend.title = element_text(size=32),
          legend.text = element_text(size=24))
  
  compare.res |> # compare GAMs to one another
    filter(Model %in% gam.names) |> 
    #filter(mod_name %in% gam.names) |> 
    #filter((mod_name %in% c("gam_10_h0_logit", "gam_30_h0_logit", "gam_30_ha1_logit", "gam_10_ha1_logit", "gam_10_ha2_logit"))) |> 
    #filter(!(mod_name %in% c("gam_10_ha2_logit", "gam_10_ha2_pois", "gam_30_ha2_pois"))) |> 
    pivot_wider(id_cols=1:2, names_from=measure, values_from=val) |>
    ggplot(aes(x=area_surveyed, y=recall, color=Model)) + 
    geom_line(linewidth=1.2) + 
    scale_colour_viridis_d(option = "A") +
    #scale_color_brewer(type="qual", palette="RdBu") +
    labs(x="Proportion of study area", y="Recall") +
    theme(axis.title = element_text(size=32),
          axis.text = element_text(size=24),
          legend.title = element_text(size=32),
          legend.text = element_text(size=24)) + xlim(0, 0.03) + ylim(0, 0.175)
}

#--------------- Fig 4--------------
#--------------- Fig 5--------------
{
  {
    rm( list=ls() )
    gc()
    library(mgcv)
    library(tidyverse)
    library(mgcViz)
  } |> suppressPackageStartupMessages()
  
  {
    mod.10 = readRDS("E:/tmp/wilt_models/gam_10_h0_logit.Rds")
    dat.10 = readRDS("mod_data/model_in/dat_10_weighted.Rds")
    dev.resid.10 = mgcv::residuals.gam(mod.10, type="deviance")
    pear.resid.10 = mgcv::residuals.gam(mod.10, type="pearson")
    dat.10 = dat.10 |> 
      add_column(resid_deviance = dev.resid.10,
                 resid_pearson = pear.resid.10,
                 fitted = mod.10$fitted.values) |> 
      mutate(resid_deviance_contrast = log(log(log(log(log(log(resid_deviance + 1)+1)+1)+1)+1)+1) )
  }
  
  g1 = getViz(mod.10)
  check(g1,
        type="deviance"
        )
  dat.10 |> filter(resid_deviance > 1) |> ggplot(aes(x=fitted, y=resid_deviance)) + geom_hex() +
    theme(
      axis.text = element_text(size=24),
      axis.title = element_text(size=32),
      legend.title = element_text(size=24),
      legend.text = element_text(size=18))
  
  
  names.tmp = dat.10 |> names()
  
  dat.subs = dat.10 |> slice_sample(prop=0.1) |> filter(resid_deviance <= 1)
  
  ggplot() +
    geom_smooth(data=dat.10, aes(x=fitted, y=resid_deviance)) + 
    geom_point(data=dat.subs, aes(x=fitted, y=resid_deviance), alpha=0.1)
  
  for(i in 4:16) {
    var.tmp = names.tmp[i]
    p1 = dat.10 |> 
      select(all_of(c(var.tmp, "resid_deviance") ) ) |> 
      ggplot(aes_string("x"=var.tmp, "y"="resid_deviance")) + 
      geom_point(data= dat.subs, aes_string("x"=var.tmp, "y"="resid_deviance"), alpha=0.1) +
      geom_smooth(method=)
    p1 |> plot()
    cat("\nany key for next")
    prpt = readline(prompt="...")
  }
  
  wilt.pts = dat.10 |> 
    filter(ow_rast_10 > 0)
  
  {
    legend = c("red")
    p1 = dat.10 |> 
      filter(resid_deviance < 2) |> 
      ggplot() + 
      geom_tile(data=dat.10, aes(x=x,y=y, fill=resid_deviance_contrast )) + scale_fill_viridis_c(option="D") +
      geom_point(data=wilt.pts, aes(x=x,y=y, color=legend), alpha=0.1) + scale_color_manual(label="Infection sites", values=legend) +
      theme(axis.text = element_text(size=18),
            axis.title = element_text(size=24),
            legend.title = element_text(size=24),
            legend.text = element_text(size=18))
    p1 |> plot()
  }
  
  
  {
    mod.30 = readRDS("E:/tmp/wilt_models/gam_30_h0_logit.Rds")
    dat.30 = readRDS("mod_data/model_in/dat_30_weighted.Rds")
    dev.resid.30 = mgcv::residuals.gam(mod.30, type="deviance")
    pear.resid.30 = mgcv::residuals.gam(mod.30, type="pearson")
    dat.30 = dat.30 |> 
      add_column(resid_deviance = dev.resid.30,
                 resid_pearson = pear.resid.30,
                 fitted = mod.30$fitted.values)
  }
  
  g2 = getViz(mod.30)
  check(g2,
        type="deviance"
  )
  
  dat.30 |> 
    filter(resid_deviance > 1) |> 
    ggplot(aes(x=fitted, y=resid_deviance)) + geom_hex() +
    theme(axis.text = element_text(size=24),
          axis.title = element_text(size=32),
          legend.title = element_text(size=24),
          legend.text = element_text(size=18))
  
  names.tmp = dat.30 |> names(); names.tmp
  
  dat.subs = dat.30 |> filter(resid_deviance <= 2) |> slice_sample(prop=0.05)
  
  for(i in 4:12) {
    var.tmp = names.tmp[i]
    p.dat = dat.30 |> 
      select(all_of(c(var.tmp, "resid_deviance") ) ) |> 
      filter(resid_deviance <= 1) 
    p1 = ggplot() + 
      geom_point(data = dat.subs, aes_string("x"=var.tmp, "y"="resid_deviance"), alpha=0.1) +
      geom_smooth(data = p.dat, aes_string("x"=var.tmp, "y"="resid_deviance"))
    p1 |> plot()
    cat("\nany key for next")
    prpt = readline(prompt="...")
  }
  
  wilt.pts = dat.30 |> 
    filter(ow_rast_30 > 0)
  
  dat.30$resid_deviance |> hist()
  
  dat.mod.30 = dat.30 |> 
    mutate(resid_deviance_contrast = log(log(log(resid_deviance + 1)+1)+1) )
  {
    legend = c("red")
   p2 = ggplot() + 
      geom_tile(data=dat.mod.30, aes(x=x,y=y,fill=resid_deviance_contrast )) + 
      geom_point(data=wilt.pts, aes(x=x,y=y, color=cols), alpha=0.1) + 
      scale_fill_viridis_c(option="D") +
      scale_color_manual(label="Infection sites",values=legend) + 
      theme(axis.text = element_text(size=18),
            axis.title = element_text(size=24),
            legend.title = element_text(size=24),
            legend.text = element_text(size=18))
   p2 |> plot()
    }
  
}
#---------------- Fig 6 ---------
{
  mod.10 = readRDS("E:/tmp/wilt_models/gam_10_h0_logit.Rds")
  mod.10 |> summary()
  
  g1 = getViz(mod.10)
  o <- plot( sm(g1, 2) )
  o
  o + l_fitLine(colour = "red") + 
    l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
    l_points(shape = 19, size = 1, alpha = 0.1)
  
  plot(sm(g1, 1)) + l_fitRaster() + l_fitContour() + l_points()
}
#--------------- Table data --------------

{ ### table 2 ###
  library(terra)
  dat.10 = readRDS("mod_data/model_in/dat_10.Rds")
  dat.10.rast = terra::rast("clean_data/joindat_10_800_50.tif")
  names(dat.10.rast)
  topo.df = dat.10.rast$topo_wet |> as.data.frame()
  topo.df |> summarise(min = min(topo_wet, na.rm=T),
                       max = max(topo_wet, na.rm=T))
  
  sand.df = dat.10.rast$sand |> as.data.frame()
    sand.df |> summarise(min = min(sand, na.rm=T), 
                         max = max(sand, na.rm=T))
  
  bd.df = dat.10.rast$bd |> as.data.frame()
  bd.df |> summarise(min = min(bd, na.rm=T),
                     max = max(bd, na.rm=T))
  
  ph.df = dat.10.rast$ph |> as.data.frame()
  ph.df |> summarise(min = min(ph, na.rm=T),
                     max = max(ph, na.rm=T))
  
  dat.10 |> names()
  
  dat.10 |> select(elev) |> summarise(min = min(elev, na.rm=T),
                                      max = max(elev, na.rm=T))
  dat.10 |> select(plan_curv) |> summarise(min = min(plan_curv, na.rm=T),
                                      max = max(plan_curv, na.rm=T))
  dat.10 |> pull(aspect) |> hist()
  dat.10 |> pull(slope) |> hist()
  dat.10 |> select(channel_dist_s5) |> summarise(min = min(channel_dist_s5, na.rm=T),
                                               max = max(channel_dist_s5, na.rm=T))
  dat.10 |> pull(slope) |> hist()
  dat.30 = readRDS("mod_data/model_in/dat_30.Rds")
  dat.30 |> names()
}

#--------------- old, misc ------------------
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
#-------------- reply to reviewers --------------
dat.30 = terra::rast("clean_data/joindat_30_800_50.tif")
ow.pts = terra::vect("clean_data/ow_pts_clean.shp")
{
  dat.30$sand |> plot(main="Soil Sand (%)")
  points(ow.pts, cex=0.5, col="red", alpha=0.2)
}

