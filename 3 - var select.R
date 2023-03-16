{ ###### PKG #####
  rm(list=ls())
  gc()
  library(spatstat)
  library(terra)
  library(tidyverse)
  library(assertthat)
  library(ggpubr)
} |> suppressPackageStartupMessages()

{ ####### TEST TP SMOOTH ########
  #
  rm(list=ls())
  gc()
  source("functions/mod/pear_corr_.R")
  source("functions/mod/var_vs_response_.R")
  source("functions/mod/var2_vs_response_.R")
  #
  dat.30 = terra::rast("clean_data/joindat_30_800_50.tif") |> 
    as.data.frame(xy=T) 
  names(dat.30)
  #
  dat.mod = dat.30 |> 
    select(-c(3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 17, 18, 19, 22, 29, 30, 33, 34, 36)) |> 
    relocate(ow_rast_30)
  names(dat.mod)
  #
  dat.mod |> pear_corr_()
  
  p1 = var_vs_response_(dat.mod, "ow_rast_30", "hillshade"); plot(p1)
  (var2_vs_response_(dat.mod, "ow_rast_30", "aspect", "hillshade"))
  
  
  library(mgcv)
  {
    t1 <- Sys.time()
    # tst.gam <- mgcv::bam(ow_rast_30 ~ 
    #                        s(x, y) + te(wl2_oakprob_30, slope, elev, aspect, bs=c("ps", "ps", "ps", "cc")),
    #                      data=dat.30, family=poisson(link="log"), method="fREML", nthreads=12)
    tst.gam <- mgcv::bam(ow_rast_30 ~
                           s(x, y) + 
                           te(slope, hillshade, aspect, bs=c("ps", "ps", "cc")) +
                           te(topo_wetness, total_catch, conv_ind, bs=c("ps", "ps", "ps")),
                         data=dat.30, family=poisson(link="log"), method="fREML", nthreads=12)
    t2 <- Sys.time()
    cat(paste0("\ntime = ", difftime(t2, t1, units="mins"), " min"))
    cat(paste0("\nAIC = ", AIC(tst.gam)))
    #saveRDS(tst.gam, file="model_data/models/gam_30_freml.Rds")
    tst.gam |> summary()
  }
  
  
  { ####### visualize ########
    library(mgcViz)
    {
      b = getViz(tst.gam)
      plot( sm(b, 1) )
      plotRGL( sm(b, 2), fix=c("wl2_oakprob_30" = 50, "slope"= 0.05), residuals=F )
    }
  }
}

{ ####### Load data, add interactions, save ########
  { ###### init #######
    source("functions/add_interact_.R")
    interact.v.10 <- readRDS("clean_data/interact_v.Rds")
    interact.v.30 <- readRDS("clean_data/interact_v_30.Rds")
    vars.filt.10 <- readRDS("mod_data/full_dat/view_vars_10.Rds")
    vars.filt.30 <- readRDS("mod_data/full_dat/view_vars_30.Rds")
  }
  
  { ####### 10-m add interactions and filter to vars, save #######
    dat.intr.10 <- terra::rast("clean_data/joindat_10_800_50.tif") |> 
      as.data.frame(xy=T) |>
      mutate(wl2_oakprob_10 = if_else(wl2_cls_10 == 4230, true=wl2_oakprob_10, false=0)) |> 
      add_interact_(
        interact.v = interact.v.10, 
        verbose=T, DBG=T) |> 
      select( all_of(c("ow_rast_10", "x", "y", vars.filt.10)) ) |> 
      drop_na()
    saveRDS(dat.intr.10, file="mod_data/full_dat/dat_interact_10.Rds")
    rm(dat.intr.10)
    gc()
  }
  
  { ####### 30-m add interactions and filter to vars, save #######
    dat.intr.30 <- terra::rast("clean_data/joindat_30_800_50.tif") |> 
      as.data.frame(xy=T) |> 
      mutate(wl2_oakprob_30 = if_else(wl2_cls_30 == 4230, true=wl2_oakprob_30, false=0)) |>
      add_interact_(
        interact.v = interact.v.30, 
        verbose=T, DBG=T) |> 
      select( all_of(c("ow_rast_30", "x", "y", vars.filt.30)) ) |> 
      drop_na()
    saveRDS(dat.intr.30, file="mod_data/full_dat/dat_interact_30.Rds")
    rm(dat.intr.30)
    gc()
  }
}

{ ####### Identify covarying predictors and plot response vs candidate covars ########
    dat.intr.30 <- readRDS("mod_data/full_dat/dat_interact_30.Rds")
    dat.intr.10 <- readRDS("mod_data/full_dat/dat_interact_10.Rds")
    source("functions/mod/pear_corr_.R")
    
    dat.intr.10 |> 
      select(-x, -y, -ow_rast_10) |> 
      pear_corr_()
    dat.intr.30 |> 
      select(-x, -y, -ow_rast_30) |> 
      pear_corr_()
    #vars.filt.10 <- readRDS("mod_data/full_dat/view_vars_10.Rds")
    
    compare.list.30 <- list(
      c("conv_ind x prof_curv", "plan_curv x prof_curv", "ph x prof_curv", "bd x prof_curv", "prof_curv"), # 1
      c("bd x conv_ind", "conv_ind", "conv_ind x plan_curv", "bd x plan_curv", "plan_curv"), # 2
      c("ph x rsp", "rsp", "channel_net_s5", "channel_net_s7"), # 3
      c("conv_ind x hillshade", "hillshade x plan_curv", "ph x hillshade", "bd x hillshade", "hillshade"), # 4
      c("bd x ph", "ph x conv_ind", "ph x plan_curv"), # 5
      c("ph x topo_wetness", "topo_wetness"), # 6
      c("sand x conv_ind", "sand x plan_curv", "ph x sand", "bd x sand", "sand"), # 7
      c("slope", "topo_wetness"), # 8
      c("bd x ph", "ph x sand"), # 9
      c("ph x sand", "clay"), # 10
      c("topo_wetness", "bd x conv_ind", "conv_ind x hillshade"), # 11
      c("elev", "channel_net_s5") # 12
    )
    
    compare.list.10 <- list(
      c("bd x hillshade", "hillshade"),
      c("prof_curv", "conv_ind x prof_curv", "plan_curv x prof_curv", "bd x prof_curv", "ph x prof_curv"),
      c("ph x rsp", "rsp", "channel_dist_s5", "channel_dist_s7", "val_depth"),
      c("ph x topo_wet", "topo_wet"),
      c("ph x sand", "sand x prof_curv", "bd x sand", "sand", "sand x conv_ind", "sand x plan_curv"),
      c("sand x conv_ind", "sand x plan_curv"),
      c("bd x conv_ind", "bd x plan_curv"),
      c("bd x ph", "ph x conv_ind", "ph x plan_curv"),
      c("ph x wl2_oakprob_10", "wl2_oakprob_10")
    )
    
    # for each comparison, plot covar quintiles (left) and deciles (right) vs disease rate 
    for(i in seq_along(compare.list.10)) {
      vars.extr <- compare.list.10[[i]]
      cat(paste0("\ni = ", i, "\tvar.tmp = ", paste0(vars.extr, collapse=" , ")))
      dir.create(paste0("expl_data/var_compare/10/", i))
      for(j in seq_along(vars.extr)) {
          var.tmp <- vars.extr[j]
          plot.dat <- dat.intr.10 |> 
            select( "ow_rast_10", all_of(var.tmp) ) |> 
            drop_na() |> 
            mutate(quintile = ntile(!!sym(var.tmp), n=5),
                   decile = ntile(!!sym(var.tmp), n=10))
          left.plot <- plot.dat |> 
            group_by(quintile) |> 
            summarise(ow_rate = sum(ow_rast_10)/n() ) |> 
            ggplot(aes(x=quintile, y=ow_rate, fill=ow_rate)) + geom_col()
          right.plot <- plot.dat |> 
            group_by(decile) |> 
            summarise(ow_rate = sum(ow_rast_10)/n() ) |> 
            ggplot(aes(x=decile, y=ow_rate, fill=ow_rate)) + geom_col()
          fig <- ggarrange(left.plot, right.plot, labels=c(paste0(var.tmp, " quintile"), paste0(var.tmp, " decile")))
          jpeg(filename=paste0("expl_data/var_compare/10/", i, "/", j, ".jpg"))
          plot(fig)
          dev.off()
      }
    }
}

{ ###### Select predictors to keep #######
  source('functions/mod/pear_corr_.R')
  source("functions/mod/condition_index_.R")
  source("functions/mod/var_vs_response_.R")
  
    { ####### Select predictors for 30-m data, check cond_ind after removal ######
      names(dat.intr.30)
      var.subs.30 <- dat.intr.30 |>
        select(ow_rast_30,
               wl2_oakprob_30, elev, aspect, # removed val_dep, gw_30, slope, clay,
               `bd x prof_curv`, # 1
               `bd x conv_ind`, # 2
               `channel_net_s5`, #`channel_net_s7`, # 3
               `conv_ind x hillshade`, # 4
               #`bd x ph`, # 5 - not linear
               `topo_wetness`, # 6
               `ph x sand`  # 7
        )
      var.subs.30 |> pear_corr_()
      
      (var.subs.30 |> 
        select(-ow_rast_30) |> 
        condition_index_())
      
      mod.dat.30 <- var.subs.30 |> 
        add_column(y = dat.intr.30[["y"]], .after="ow_rast_30") |> 
        add_column(x = dat.intr.30[["x"]], .after="ow_rast_30")
      
      saveRDS(mod.dat.30, file="mod_data/model_in/dat_30.Rds")
    }
    
    { ####### Select predictors for 10-m data, check cond_ind after removal, vis all ######
      names(dat.intr.10)
      var.subs.10 <- dat.intr.10 |>  #**important*
        select(ow_rast_10,
               elev, aspect, slope, channel_dist_s5, wl2_oakprob_10, plan_curv, val_depth, # removed sand, rsp, channel_dist_s7, ph x rsp, val_depth
               `bd x hillshade`, # 1
               `bd x prof_curv`, # 2
               # 3 is channel_dist_s5
               `ph x topo_wet`, # 4
               `ph x sand`, # 5
               #`sand x conv_ind`, #**non-linear* 6 
               `bd x conv_ind`, # 7
               `bd x ph`, #**nonlin* # 8 
              )
      var.subs.10 |> pear_corr_()
      
      (var.subs.10 |> 
        select(-ow_rast_10) |> 
        condition_index_())
      
      { #### plot comparison ###
        p1 = var_vs_response_(var.subs.10, "ow_rast_10", "ph x topo_wet")
        p1
        p2 = var_vs_response_(var.subs.10, "ow_rast_10", "slope")
        p2
      }
      
      mod.dat.10 <- var.subs.10 |> 
        add_column(y = dat.intr.10[["y"]], .after="ow_rast_10") |> 
        add_column(x = dat.intr.10[["x"]], .after="ow_rast_10")
      saveRDS(mod.dat.10, file="mod_data/model_in/dat_10.Rds")
  }
}
  
{ ###### FIT GAM on 10-m data ########
      rm(list=ls())
      gc()
      
      dat.10 <- readRDS("mod_data/model_in/dat_10.Rds")
      dat.10 |> str()
      
      
      { ## Replace spaces in variable names with underscores
        nm.tmp = names(dat.10)
        und.ind = str_detect(nm.tmp, pattern="_") |> 
          which()
        names(dat.10)[und.ind] <- nm.tmp[und.ind] |> 
          str_replace_all(pattern=" ", replacement="_")
        names(dat.10)
        dat.10 |> str()
      }
      
      { ####### scale covariates ########
        for(i in 4:ncol(dat.10)) {
          dat.10[,i] <- scale(dat.10[,i])
        }
      }
      
      library(mgcv)
      {
        t1 <- Sys.time()
        mod.tst <- bam(
          ow_rast_10 ~ s(x, y) + s(wl2_oakprob_10, bs="cs"),
                      data = dat.10,
                      family=poisson(link="log"), method="fREML", nthreads=6,
                      #coef=coef.last
                       )
        coef.last <- mod.tst$coefficients
        t2 <- Sys.time()
        cat(paste0("\ndifftime = ", difftime(t2, t1, units="mins"), " min"))
        summary(mod.tst)
        plot(mod.tst)
      }
}


{ ###### 30-m GAM #########
  rm(list=ls())
  gc()
  library(mgcv)
  dat.mod.30 <- readRDS("mod_data/model_in/dat_30.Rds")
  
  {
    dat.mod.30 |> str()
    ind.tmp = stringr::str_detect(names(dat.mod.30), pattern=" x ") |> 
      which(); ind.tmp
    names(dat.mod.30)[ind.tmp] <- names(dat.mod.30)[ind.tmp] |> 
      stringr::str_replace_all(pattern=" x ", replacement="_x_")
    str(dat.mod.30)
  }
  
  for(i in seq_len(ncol(dat.mod.30))) {
    if(i >= 4) {
      dat.mod.30[,i] <- scale(dat.mod.30[,i])
    }
  }
  dat.mod.30 |> names()
  {
    t1 <- Sys.time()
    tst.gam <- mgcv::bam(ow_rast_30 ~ 
                           s(x, y) + 
                           s(wl2_oakprob_30, bs="cs") + 
                           # s(elev, aspect, bs=c("ps", "cc") ) + 
                           s(elev, bs="cs") +
                           s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs"),
                           #s(topo_wetness, bs="cs") + s(conv_ind_x_hillshade, bs="cs") + #,
                           # bd_x_prof_curv + bd_x_conv_ind + conv_ind_x_hillshade + ph_x_sand,
                           # s(bd_x_prof_curv, bs="cs") + s(bd_x_conv_ind, bs="cs") + 
                           # s(conv_ind_x_hillshade, bs="cs") + s(ph_x_sand, bs="cs"),
                       data=dat.mod.30, family=poisson(link="log"), method="fREML", nthreads=6)
    t2 <- Sys.time()
    cat(paste0("\ntime = ", difftime(t2, t1, units="mins"), " min"))
    cat(paste0("\nAIC = ", AIC(tst.gam)))
    #saveRDS(tst.gam, file="model_data/models/gam_30_freml.Rds")
  }
  tst.gam |> summary()
  tst.gam |> plot()
  library(mgcViz)
  b = getViz(tst.gam)
  plot( sm(b, 1) )
  plot( sm(b, 2) )
  plot( sm(b, 3) )
  plot( sm(b, 4) )
}

{ ###### 10-m GAM #########
  rm(list=ls())
  gc()
  library(mgcv)
  
  { ###### replace variable names with spaces with underscores ######
    dat.mod.10 <- readRDS("model_data/var12_10_moddat.Rds")
    dat.mod.10 |> str()
    ind.tmp = stringr::str_detect(names(dat.mod.10), pattern=" x ") |> 
      which(); ind.tmp
    names(dat.mod.10)[ind.tmp] <- names(dat.mod.10)[ind.tmp] |> 
      stringr::str_replace_all(pattern=" x ", replacement="_x_")
    dat.mod.10 |> names()
  }
  # scale
  for(i in seq_len(ncol(dat.mod.10))) {
    if(i >= 2) {
      dat.mod.10[,i] <- scale(dat.mod.10[,i])
    }
  }
  #dat.mod.10 |> str()
  
  { ###### fit gam ######
    t1 <- Sys.time()
    tst.gam <- mgcv::bam(ow_rast_10 ~ 
                           s(x, y) + s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(aspect, bs="cc") +
                           # te(elev, aspect, bs=c("ps", "cc")) +
                           s(ph_x_wl2_oakprob_10, bs="cs") + s(bd_x_hillshade, bs="cs") + s(bd_x_plan_curv, bs="cs"),
                           #ph_x_topo_wet + ph_x_rsp + ph_x_wl2_oakprob_10 + bd_x_sand + bd_x_ph + bd_x_hillshade + bd_x_plan_curv + conv_ind_x_prof_curv,
                         data=dat.mod.10, family=poisson(link="log"), method="fREML", nthreads=6)
    t2 <- Sys.time()
    cat(paste0("\ntime = ", difftime(t2, t1, units="mins"), " min"))
    cat(paste0("\nAIC = ", AIC(tst.gam)))
  }
  { ###### visualize ######
    tst.gam |> summary()
    #tst.gam |> plot()
    library(mgcViz)
    b = getViz(tst.gam)
    plotRGL( sm(b, 1) )
    plot( sm(b, 2) )
    plotRGL( sm(b, 3), residuals=T )
    plot( sm(b, 4) )
  }
}

{ ######### Get CPM data ready ###########
  rm(list=ls())
  gc()
  
  { ####### create quadschemes ########
    ow.ppp.10 = readRDS("clean_data/ow_ppp_10.Rds")
    ow.ppp.30 = readRDS("clean_data/ow_ppp_30.Rds")
    quads.10 = quadscheme(data=ow.ppp.10, method="grid", eps=30)
    quads.10
    saveRDS(quads.10, file="mod_data/kppm/quad_10_30.Rds")
    
    quads.30 = quadscheme(data=ow.ppp.30, method="grid", eps=30)
    quads.30
    saveRDS(quads.30, file="mod_data/kppm/quad_30_30.Rds")
    
  }
  
  { ####### create im lists ########
    source("functions/mod/tbl_to_im_.R")
    
    dat.mod.30 <- readRDS("mod_data/model_in/dat_30.Rds")
    vars.in = names(dat.mod.30)[4:ncol(dat.mod.30)]; vars.in
    
    im.list.30 = tbl_to_im_(
      dat.mod.30, 
      vars.in, 
      terra::crs('EPSG:3071'),
      verbose=T, DBG=T)
    saveRDS(im.list.30, file="mod_data/kppm/im_list_30.Rds")
    
    dat.mod.10 <- readRDS("mod_data/model_in/dat_10.Rds")
    vars.in = names(dat.mod.10)[4:ncol(dat.mod.10)]; vars.in
    im.list.10 = tbl_to_im_(
      dat.mod.10,
      vars.in, 
      terra::crs('EPSG:3071'),
      verbose=T, DBG=T)
    saveRDS(im.list.10, file="mod_data/kppm/im_list_10.Rds")
  }
}

{ #######  Fit 30-m KPPM #########
  quads.30 = readRDS("mod_data/kppm/quad_30_30.Rds")
  im.list.30 = readRDS("mod_data/kppm/im_list_30.Rds"); names(im.list.30)
  {
    kppm.30.thom = kppm(
      X = quads.30, 
      trend = ~.,
      clusters="Thomas",
      penalised = T,
      data=im.list.30,
      subset=ow.ppp.30$window,
      method = "palm",
      improve.type = "quasi")
    
    kppm.30.mat = kppm(
      X = quads.30, 
      trend = ~.,
      clusters="MatClust",
      penalised = T,
      data=im.list.30,
      subset=ow.ppp.30$window,
      method = "palm",
      improve.type = "quasi")
    
    kppm.30.cauch = kppm(
      X = quads.30, 
      trend = ~.,
      clusters="Cauchy",
      penalised = T,
      data=im.list.30,
      subset=ow.ppp.30$window,
      method = "palm",
      improve.type = "quasi")
    
    kppm.30.vargam = kppm(
      X = quads.30, 
      trend = ~.,
      clusters="VarGamma",
      penalised = T,
      data=im.list.30,
      subset=ow.ppp.30$window,
      method = "palm",
      improve.type = "quasi")
    
    kppm.30
    kppm.30 |> plot()
  }
}

{ ####### Fit 10-m KPPM #########
  
  im.list.10 = readRDS("mod_data/kppm/im_list_10.Rds"); names(im.list.10)
  {
    kppm.10.full = kppm(
      X = ow.ppp.10, 
      trend = ~.,
      clusters="Thomas",
      data=im.list.10,
      rmax=1000,
      method = "palm",
      improve.type = "quasi")
    
    kppm.10.subs = kppm(
      X = ow.ppp.10, 
      trend = ~ elev + channel_dist_s5 + ,
      clusters="Thomas",
      data=im.list.10,
      rmax=1000,
      method = "palm",
      improve.type = "quasi")
  }
  {  
    kppm.10
    kppm.10 |> summary()
    kppm.10 |> plot()
  }
}
