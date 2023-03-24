{ ###### PKG #####
  rm(list=ls())
  gc()
  #library(spatstat)
  library(terra)
  library(tidyverse)
  library(assertthat)
  library(ggpubr)
} |> suppressPackageStartupMessages()

#-------------- setup, select vars ---------------------

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
    source("functions/mod/condition_index_.R")
    
    dat.intr.10 |> 
      select(-x, -y, -ow_rast_10) |> 
      pear_corr_()
    dat.intr.30 |> 
      select(-x, -y, -ow_rast_30) |> 
      pear_corr_()
    
    (dat.intr.10 |> 
      condition_index_())
    (dat.intr.30 |> 
      condition_index_())
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

#-------------- GAM ---------------------

{ ###### GAM 10-m ########
  { ##### INIT #####
    {
      rm(list=ls())
      gc()
      library(tidyverse)
      library(assertthat)
      library(mgcv)
    } |> suppressPackageStartupMessages()
    { ##### data setup #####
      save.pth = "E:/tmp/wilt_models/"
      dat.10 = readRDS("mod_data/compare/wilt_split_tbl_10.Rds")
      dat.names = names(dat.10); dat.names
      
      { ###### correct var names #####
        for(i in seq_along(dat.names)) {
          nam.tmp = dat.names[i]
          if(str_detect(nam.tmp, pattern=" x ")) {
            nam.repl = str_replace(nam.tmp, pattern=" x ", replacement="_x_")
            names(dat.10)[i] = nam.repl
          }
        }
        names(dat.10)
      }
      
      { ###### scale vars #####
        for(i in 4:16) {
          dat.10[,i] = scale(dat.10[,i])
        }
      }
      
      { ##### add logit response & weights #####
        n.no = sum(dat.10$ow_rast_10 == 0); n.no
        n.yes = sum(dat.10$ow_rast_10 > 0); n.yes
        dat.10 = dat.10 |> 
          mutate(logit_resp = if_else(ow_rast_10 > 0, 1, 0),
                 logit_resp = as.integer(logit_resp),
                 weights_raw = if_else(ow_rast_10 == 0, 1/n.no, 1/n.yes),
                 logit_weights_raw = if_else(ow_rast_10 > 0, ow_rast_10*weights_raw, weights_raw),
                 weights = weights_raw/mean(weights_raw),
                 logit_weights = logit_weights_raw/mean(logit_weights_raw) 
          )
        names(dat.10)
      }
      dat.10$logit_resp |> class()
      saveRDS(dat.10, file="mod_data/model_in/dat_10_weighted.Rds")
    }
  } 
  
  { ####### GAMs for hypothesis testing ########
    dat.10 = readRDS("mod_data/model_in/dat_10_weighted.Rds")
    dat.10 |> str()
    save.pth = "E:/tmp/wilt_models/"
    { ###### H0 - spatial + hosts + env ###### 
      { ##### H0 poisson #####
        t1 = Sys.time()
        gam.10.h0.pois = mgcv::bam(ow_rast_10 ~ s(x, y) + s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + 
                                     s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc") + s(channel_dist_s5, bs="cs") + s(plan_curv, bs="cs"),
                                   data=dat.10, 
                                   #weights=dat.10$weights, 
                                   family=poisson(link="log"), 
                                   method="fREML", 
                                   nthreads=6)
        t2 = Sys.time()
        cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
        saveRDS(gam.10.h0.pois, file= paste0(save.pth, "gam_10_h0_pois.Rds") )
        (summary(gam.10.h0.pois))
        rm(gam.10.h0.pois)
      }
      
      { ##### H0 logit #####
        t1 = Sys.time()
        gam.10.h0.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + 
                                      s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc") + s(channel_dist_s5, bs="cs") + s(plan_curv, bs="cs"),
                                    data = dat.10, 
                                    #weights= dat.10$logit_weights,
                                    family=binomial(link="logit"), 
                                    method="fREML", 
                                    nthreads=6)
        t2 = Sys.time()
        cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
        saveRDS(gam.10.h0.logit, file=paste0(save.pth, "gam_10_h0_logit.Rds") )
        (summary(gam.10.h0.logit))
        rm(gam.10.h0.logit)
      }
    }
    
    { ######### HA1 - spatial + hosts ##########
      # ha1 - poisson
      t1 = Sys.time()
      gam.10.ha1.pois = mgcv::bam(ow_rast_10 ~ s(x, y) + s(wl2_oakprob_10, bs="cs"),
                                  data=dat.10, 
                                  #weights=dat.10$weights,
                                  family=poisson(link="log"), 
                                  method="fREML", 
                                  nthreads=6)
      t2 = Sys.time()
      cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
      saveRDS(gam.10.ha1.pois, file= paste0(save.pth, "gam_10_ha1_pois.Rds") )
      rm(gam.10.ha1.pois)
      
      # ha1 - logit
      t1 = Sys.time()
      gam.10.ha1.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_10, bs="cs"),
                                   data=dat.10,
                                   #weights=dat.10$logit_weights,
                                   family=binomial(link="logit"), 
                                   method="fREML",
                                   nthreads=6)
      t2 = Sys.time()
      cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
      saveRDS(gam.10.ha1.logit, file= paste0(save.pth, "gam_10_ha1_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(gam.10.ha1.logit)
    }
    
    { ######### HA2 - hosts + env no spatial ##########
      # ha2 - poisson
      t1 = Sys.time()
      gam.10.ha2.pois = mgcv::bam(ow_rast_10 ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + 
                                    s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc") + s(channel_dist_s5, bs="cs") + s(plan_curv, bs="cs"),
                                  data = dat.10,
                                  #weights = dat.10$weights ,
                                  family=poisson(link="log"), 
                                  method="fREML", 
                                  nthreads=6)
      t2 = Sys.time()
      cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
      saveRDS(gam.10.ha2.pois, file= paste0(save.pth, "gam_10_ha2_pois.Rds") )
      rm(gam.10.ha2.pois)
      
      # ha2 - logit
      {
        t1 = Sys.time()
        gam.10.ha2.logit = mgcv::bam(logit_resp ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + 
                                                s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc") + s(channel_dist_s5, bs="cs") + s(plan_curv, bs="cs"),
                                   data=dat.10,
                                   #weights= dat.10$logit_weights,
                                   family=binomial(link="logit"), 
                                   method="fREML", 
                                   nthreads=6)
        t2 = Sys.time()
        cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
        saveRDS(gam.10.ha2.logit, file= paste0(save.pth, "gam_10_ha2_logit.Rds") ); #gam.10.ha2.logit = readRDS(paste0(save.pth, "gam_10_ha2_logit.Rds"))
        cat(crayon::bgRed("\n\nDONE\n\n")) # gam.10.ha2.logit |> summary()
        rm(gam.10.ha2.logit)
      }
    }
  }
  
  { ####### Variable importance GAMs #######
    dat.10 = readRDS("mod_data/model_in/dat_10_weighted.Rds")
    save.pth = "E:/tmp/wilt_models/var_importance/"
    #dat.10 |> str()
    { ##### null / saturated ######
      #
      # SATURATED MODEL
      #
      # sat.mod = readRDS("E:/tmp/wilt_models/gam_10_h0_logit.Rds")
      # sat.mod |> summary()
      # rm(sat.mod)
      # gam.10.h0.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc") + s(channel_dist_s5, bs="cs") + s(plan_curv, bs="cs"),
      #                             data = dat.10, 
      #                             family=binomial(link="logit"), 
      #                             method="fREML", 
      #                             nthreads=6)
      
      # sat model
      gam.10.sat = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + s(bd_x_hillshade, bs="cs"),
                             data=dat.10,
                             family=binomial(link="logit"), 
                             method="fREML", 
                             nthreads=6)
      saveRDS(gam.10.sat, file=paste0(save.pth, "gam_10_sat.Rds"))
      rm(gam.10.sat)
      cat(crayon::bgRed("\n\ndone -1 of 5\n\n"))
      
      
      # null model
      mean.resp = mean(dat.10$logit_resp); mean.resp;
      gam.10.null = mgcv::bam(logit_resp ~ rep(mean.resp, nrow(dat.10) ),
                              data=dat.10,
                              family=binomial(link="logit"), 
                              method="fREML", 
                              nthreads=6)
      saveRDS(gam.10.null, file=paste0(save.pth, "gam_10_null.Rds"))
      gam.10.null |> summary()
      AIC(gam.10.null)
      rm(gam.10.null)
      cat(crayon::bgRed("\n\ndone 0 of 5\n\n"))
      
      # remove x,y
      gam.10.rmv.xy = mgcv::bam(logit_resp ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + s(bd_x_hillshade, bs="cs"),
                                data=dat.10,
                                family=binomial(link="logit"), 
                                method="fREML", 
                                nthreads=6)
      saveRDS(gam.10.rmv.xy, file=paste0(save.pth, "gam_10_rmv_xy.Rds"))
      rm(gam.10.rmv.xy)
      cat(crayon::bgRed("\n\ndone 1 of 5\n\n"))
      
      # remove wl2 oak class prob
      gam.10.rmv.wl2 = mgcv::bam(logit_resp ~ s(x, y) + s(elev, bs="cs") + s(val_depth, bs="cs") + s(bd_x_hillshade, bs="cs"),
                                 data=dat.10,
                                 family=binomial(link="logit"), 
                                 method="fREML", 
                                 nthreads=6)
      saveRDS(gam.10.rmv.wl2, file=paste0(save.pth, "gam_10_rmv_wl2.Rds"))
      rm(gam.10.rmv.wl2)
      cat(crayon::bgRed("\n\ndone 2 of 5\n\n"))
      
      # remove elev
      gam.10.rmv.elev = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_10, bs="cs") + s(val_depth, bs="cs") + s(bd_x_hillshade, bs="cs"),
                                  data=dat.10,
                                  family=binomial(link="logit"), 
                                  method="fREML", 
                                  nthreads=6)
      saveRDS(gam.10.rmv.elev, file=paste0(save.pth, "gam_10_rmv_elev.Rds"))
      rm(gam.10.rmv.elev)
      cat(crayon::bgRed("\n\ndone 3 of 5\n\n"))
      
      # remove val_dep
      gam.10.rmv.valdep = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(bd_x_hillshade, bs="cs"),
                                    data=dat.10,
                                    family=binomial(link="logit"), 
                                    method="fREML", 
                                    nthreads=6)
      saveRDS(gam.10.rmv.valdep, file=paste0(save.pth, "gam_10_rmv_valdep.Rds"))
      rm(gam.10.rmv.valdep)
      cat(crayon::bgRed("\n\ndone 4 of 5\n\n"))
      
      # remove bd x hillshade
      gam.10.rmv.bdxhillshade = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs"),
                                          data=dat.10,
                                          family=binomial(link="logit"), 
                                          method="fREML", 
                                          nthreads=6)
      saveRDS(gam.10.rmv.bdxhillshade, file=paste0(save.pth, "gam_10_rmv_bdxhillshade.Rds"))
      rm(gam.10.rmv.bdxhillshade)
      cat(crayon::bgRed("\n\ndone 5 of 5\n\n"))
    }
    { ##### ha2#########
      
      #sat.mod = readRDS("E:/tmp/wilt_models/gam_10_ha2_logit.Rds")
      #sat.mod |> summary()
      
      # sat mod
      ha2.10.sat = mgcv::bam(logit_resp ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + s(bd_x_hillshade, bs="cs") + s(channel_dist_s5, bs="cs") + s(aspect, bs="cc"),
                             data=dat.10,
                             family=binomial(link="logit"), 
                             method="fREML", 
                             nthreads=6)
      saveRDS(ha2.10.sat, file=paste0(save.pth, "ha2_10_sat.Rds"))
      rm(ha2.10.sat)
      
      # null model
      # 
      mean.resp = mean(dat.10$logit_resp); mean.resp; 
      ha2.10.null = mgcv::bam(logit_resp ~ rep(mean.resp, nrow(dat.10)),
                                               data=dat.10,
                                               family=binomial(link="logit"), 
                                               method="fREML", 
                                               nthreads=6)
      saveRDS(ha2.10.null, file=paste0(save.pth, "ha2_10_null.Rds"))
      rm(ha2.10.null)
      
      # rm wl2 model
      #
      ha2.10.rm.wl2 = mgcv::bam(logit_resp ~ s(elev, bs="cs") + s(val_depth, bs="cs") + s(bd_x_hillshade, bs="cs") + s(channel_dist_s5, bs="cs") + s(aspect, bs="cc"),
                             data=dat.10,
                             family=binomial(link="logit"), 
                             method="fREML", 
                             nthreads=6)
      saveRDS(ha2.10.rm.wl2, file=paste0(save.pth, "ha2_10_rm_wl2.Rds"))
      rm(ha2.10.rm.wl2)
      
      # rm elev
      #
      ha2.10.rm.elev = mgcv::bam(logit_resp ~ s(wl2_oakprob_10, bs="cs") + s(val_depth, bs="cs") + s(bd_x_hillshade, bs="cs") + s(channel_dist_s5, bs="cs") + s(aspect, bs="cc"),
                             data=dat.10,
                             family=binomial(link="logit"), 
                             method="fREML", 
                             nthreads=6)
      saveRDS(ha2.10.rm.elev, file=paste0(save.pth, "ha2_10_rm_elev.Rds"))
      rm(ha2.10.rm.elev)
      
      # rm val_dep
      #
      ha2.10.rm.valdep = mgcv::bam(logit_resp ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(bd_x_hillshade, bs="cs") + s(channel_dist_s5, bs="cs") + s(aspect, bs="cc"),
                                 data=dat.10,
                                 family=binomial(link="logit"), 
                                 method="fREML", 
                                 nthreads=6)
      saveRDS(ha2.10.rm.valdep, file=paste0(save.pth, "ha2_10_rm_valdep.Rds"))
      rm(ha2.10.rm.valdep)
      
      # rm bd_x_hillshade
      #
      ha2.10.rm.bdxhs = mgcv::bam(logit_resp ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") +  s(val_depth, bs="cs") + s(channel_dist_s5, bs="cs") + s(aspect, bs="cc"),
                                 data=dat.10,
                                 family=binomial(link="logit"), 
                                 method="fREML", 
                                 nthreads=6)
      saveRDS(ha2.10.rm.bdxhs, file=paste0(save.pth, "ha2_10_rm_bdxhs.Rds"))
      rm(ha2.10.rm.bdxhs)
      
      # rm chan net s5
      #
      ha2.10.rm.cn = mgcv::bam(logit_resp ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc"),
                                 data=dat.10,
                                 family=binomial(link="logit"), 
                                 method="fREML", 
                                 nthreads=6)
      saveRDS(ha2.10.rm.cn, file=paste0(save.pth, "ha2_10_rm_cn.Rds"))
      rm(ha2.10.rm.cn)
      
      # rm aspect
      #
      ha2.10.rm.asp = mgcv::bam(logit_resp ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + s(channel_dist_s5, bs="cs") + s(bd_x_hillshade, bs="cs"),
                               data=dat.10,
                               family=binomial(link="logit"), 
                               method="fREML", 
                               nthreads=6)
      saveRDS(ha2.10.rm.asp, file=paste0(save.pth, "ha2_10_rm_asp.Rds"))
      rm(ha2.10.rm.asp)
    }
  }
  
  { ##### Prec-Rec GAMs (10-m) #####
    #{ # setup
    #   dat.10 = readRDS("mod_data/model_in/dat_10_weighted.Rds")
    #   save.pth = "E:/tmp/wilt_models/train_test/"
    #   n.no = sum(dat.10$wilt_train == 0); n.no
    #   n.yes = sum(dat.10$wilt_train > 0); n.yes
    #   dat.10 = dat.10 |> 
    #     mutate(logit_train = if_else(wilt_train > 0, 1, 0),
    #            logit_train = as.integer(logit_resp),
    #            weights_raw = if_else(wilt_train == 0, 1/n.no, 1/n.yes),
    #            logit_weights_raw = if_else(wilt_train > 0, wilt_train*weights_raw, weights_raw),
    #            weights = weights_raw/mean(weights_raw),
    #            logit_weights = logit_weights_raw/mean(logit_weights_raw) 
    #     )
    #   saveRDS(dat.10, file="mod_data/model_in/dat_10_weighted.Rds")
    # }
    dat.10 = readRDS("mod_data/model_in/dat_10_weighted.Rds")
    save.pth = "E:/tmp/wilt_models/train_test/"
    { ###### H0 - spatial + hosts + env ###### 
      { ##### H0 poisson #####
        t1 = Sys.time()
        train.gam.10.h0.pois = mgcv::bam(wilt_train ~ s(x, y) + s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + 
                                           s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc") + s(channel_dist_s5, bs="cs") + s(plan_curv, bs="cs"),
                                         data=dat.10, 
                                         #weights=dat.10$weights, 
                                         family=poisson(link="log"), 
                                         method="fREML", 
                                         nthreads=6)
        t2 = Sys.time()
        cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
        (summary(train.gam.10.h0.pois))
        saveRDS(train.gam.10.h0.pois, file= paste0(save.pth, "train_gam_10_h0_pois.Rds") )
        rm(train.gam.10.h0.pois)
      }
      { ##### H0 - logit ######
        t1 = Sys.time()
        train.gam.10.h0.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + 
                                            s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc") + s(channel_dist_s5, bs="cs") + s(plan_curv, bs="cs"),
                                          data=dat.10,
                                          #weights = dat.10$logit_weights,
                                          family=binomial(link="logit"), 
                                          method="fREML", 
                                          nthreads=6)
        t2 = Sys.time()
        cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
        saveRDS(train.gam.10.h0.logit, file= paste0(save.pth, "train_gam_10_h0_logit.Rds") )
        cat(crayon::bgRed("\n\nDONE\n\n"))
        rm(train.gam.10.h0.logit)
      }
    }
    
    { ######### HA1 - spatial + hosts ##########
      # ha1 - poisson
      t1 = Sys.time()
      train.gam.10.ha1.pois = mgcv::bam(wilt_train ~ s(x, y) + s(wl2_oakprob_10, bs="cs"),
                                        data=dat.10,
                                        #weights=dat.10$weights,
                                        family=poisson(link="log"), 
                                        method="fREML", 
                                        nthreads=6)
      t2 = Sys.time()
      cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
      saveRDS(train.gam.10.ha1.pois, file= paste0(save.pth, "train_gam_10_ha1_pois.Rds") )
      rm(train.gam.10.ha1.pois)
      
      # ha1 - logit
      t1 = Sys.time()
      train.gam.10.ha1.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_10, bs="cs"),
                                         data=dat.10, 
                                         #weights=dat.10$logit_weights,
                                         family=binomial(link="logit"), 
                                         method="fREML",
                                         nthreads=6)
      t2 = Sys.time()
      cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
      saveRDS(train.gam.10.ha1.logit, file= paste0(save.pth, "train_gam_10_ha1_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(train.gam.10.ha1.logit)
    }
    
    { ######### HA2 - hosts + env no spatial ##########
      # ha2 - poisson
      t1 = Sys.time()
      train.gam.10.ha2.pois = mgcv::bam(wilt_train ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + 
                                          s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc") + s(channel_dist_s5, bs="cs") + s(plan_curv, bs="cs"),
                                        data=dat.10, 
                                        #weights=dat.10$weights,
                                        family=poisson(link="log"), 
                                        method="fREML", 
                                        nthreads=6)
      t2 = Sys.time()
      cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
      saveRDS(train.gam.10.ha2.pois, file= paste0(save.pth, "train_gam_10_ha2_pois.Rds") )
      rm(train.gam.10.ha2.pois)
      
      # ha2 - logit
      t1 = Sys.time()
      train.gam.10.ha2.logit = mgcv::bam(logit_resp ~ s(wl2_oakprob_10, bs="cs") + s(elev, bs="cs") + s(val_depth, bs="cs") + 
                                           s(bd_x_hillshade, bs="cs") + s(aspect, bs="cc") + s(channel_dist_s5, bs="cs") + s(plan_curv, bs="cs"),
                                         data=dat.10,
                                         #weights=dat.10$logit_weights,
                                         family=binomial(link="logit"), 
                                         method="fREML", 
                                         nthreads=6)
      t2 = Sys.time()
      cat(paste0("\n\ntime = ", difftime(t2, t1, units="mins"), " min"))
      saveRDS(train.gam.10.ha2.logit, file= paste0(save.pth, "train_gam_10_ha2_logit.Rds") )
      rm(train.gam.10.ha2.logit)
    }
    cat(crayon::bgRed("\n\n\n\nDONE\n\n\n\n"))
  }
}

{ ###### GAM 30-m #########
  
  { ########## setup #########
    {
      rm(list=ls())
      gc()
      library(mgcv)
      library(terra)
      library(tidyverse)
    } |> suppressPackageStartupMessages()
    
    dat.30 = readRDS("mod_data/compare/wilt_split_tbl.Rds")
    names(dat.30)
    
    { ###### replace spaces in variable names with underscores ########
      dat.30 |> str()
      ind.tmp = stringr::str_detect(names(dat.30), pattern=" x ") |> 
        which(); ind.tmp
      names(dat.30)[ind.tmp] <- names(dat.30)[ind.tmp] |> 
        stringr::str_replace_all(pattern=" x ", replacement="_x_")
      dat.30 |> names()
    }
    
    { ######## scale vars - not spatial position #########
      for(i in 4:12) {
        dat.30[,i] <- scale(dat.30[,i])
      }
    }
    { ##### add logit response & weights #####
      n.no = sum(dat.30$ow_rast_30 == 0); n.no
      n.yes = sum(dat.30$ow_rast_30 > 0); n.yes
      dat.30 = dat.30 |> 
        mutate(logit_resp = if_else(ow_rast_30 > 0, 1, 0),
               logit_resp = as.integer(logit_resp),
               weights_raw = if_else(ow_rast_30 == 0, 1/n.no, 1/n.yes),
               logit_weights_raw = if_else(ow_rast_30 > 0, ow_rast_30*weights_raw, weights_raw),
               weights = weights_raw/mean(weights_raw),
               logit_weights = logit_weights_raw/mean(logit_weights_raw) 
        )
    }
    saveRDS(dat.30, file="mod_data/model_in/dat_30_weighted.Rds")
  }
  
  { ######### fit hypothesis GAMs ########### 
    dat.30 = readRDS("mod_data/model_in/dat_30_weighted.Rds")
    names(dat.30)
    save.pth = "E:/tmp/wilt_models/"

    { ###### H0 - spatial + hosts + env ###### 
      {######## h0 - poisson
        t1 = Sys.time()
        gam.30.h0.pois = mgcv::bam(ow_rast_30 ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(aspect, bs="cc") +
                                                s(bd_x_prof_curv, bs="cs") + s(bd_x_conv_ind, bs="cs") + s(channel_net_s5, bs="cs") + s(conv_ind_x_hillshade, bs="cs") +
                                                s(topo_wetness, bs="cs") + s(ph_x_sand, bs="cs"),
                                   data=dat.30,
                                   #weights=dat.30$weights,
                                   family=poisson(link="log"), 
                                   method="fREML", nthreads=6)
        t2 = Sys.time()
        cat(paste0("\ntime = ", difftime(t2, t1, units="mins"), " min"))
        saveRDS(gam.30.h0.pois, file= paste0(save.pth, "gam_30_h0_pois.Rds") ) # gam.30.h0.pois = readRDS(paste0(save.pth, "gam_30_h0_pois.Rds"))
      }
      gam.30.h0.pois |> summary()
      rm(gam.30.h0.pois)
      
      { #### h0 - logit
        gam.30.h0.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(aspect, bs="cc") +
                                                  s(bd_x_prof_curv, bs="cs") + s(bd_x_conv_ind, bs="cs") + s(channel_net_s5, bs="cs") + s(conv_ind_x_hillshade, bs="cs") +
                                                  s(topo_wetness, bs="cs") + s(ph_x_sand, bs="cs"),
                                    data=dat.30, 
                                    #weights=dat.30$logit_weights, 
                                    family=binomial(link="logit"), 
                                    method="fREML", 
                                    nthreads=6)
        saveRDS(gam.30.h0.logit, file= paste0(save.pth, "gam_30_h0_logit.Rds") ) # gam.30.h0.logit = readRDS(paste0(save.pth, "gam_30_h0_logit.Rds")) 
        cat(crayon::bgRed("\n\nDONE\n\n")) # gam.30.h0.logit |> summary()
        rm(gam.30.h0.logit)
      }
    }
    
    { ######### HA1 - spatial + hosts ##########
      # ha1 - poisson
      gam.30.ha1.pois = mgcv::bam(ow_rast_30 ~ s(x, y) + s(wl2_oakprob_30, bs="cs"),
                                 data=dat.30, 
                                 #weights=dat.30$weights,
                                 family=poisson(link="log"), 
                                 method="fREML", 
                                 nthreads=6)
      saveRDS(gam.30.ha1.pois, file= paste0(save.pth, "gam_30_ha1_pois.Rds") )
      rm(gam.30.ha1.pois)
      
      # ha1 - logit
      gam.30.ha1.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs"),
                                  data=dat.30, 
                                  #weights=dat.30$logit_weights,
                                  family=binomial(link="logit"), 
                                  method="fREML",
                                  nthreads=6)
      saveRDS(gam.30.ha1.logit, file= paste0(save.pth, "gam_30_ha1_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(gam.30.ha1.logit)
    }
    
    { ######### HA2 - hosts + env no spatial ##########
      # ha2 - poisson
      gam.30.ha2.pois = mgcv::bam(ow_rast_30 ~ s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(aspect, bs="cc") +
                                                s(bd_x_prof_curv, bs="cs") + s(bd_x_conv_ind, bs="cs") + s(channel_net_s5, bs="cs") + s(conv_ind_x_hillshade, bs="cs") +
                                                s(topo_wetness, bs="cs") + s(ph_x_sand, bs="cs"),
                                  data=dat.30, 
                                  #weights=dat.30$weights,
                                  family=poisson(link="log"), 
                                  method="fREML", 
                                  nthreads=6)
      saveRDS(gam.30.ha2.pois, file= paste0(save.pth, "gam_30_ha2_pois.Rds") ) #gam.30.ha2.pois= readRDS(paste0(save.pth, "gam_30_ha2_pois.Rds"))
      
      # ha2 - logit
      gam.30.ha2.logit = mgcv::bam(logit_resp ~ s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(aspect, bs="cc") +
                                                  s(bd_x_prof_curv, bs="cs") + s(bd_x_conv_ind, bs="cs") + s(channel_net_s5, bs="cs") + s(conv_ind_x_hillshade, bs="cs") +
                                                  s(topo_wetness, bs="cs") + s(ph_x_sand, bs="cs"),
                                   data=dat.30,
                                   #weights=dat.30$logit_weights,
                                   family=binomial(link="logit"), 
                                   method="fREML", 
                                   nthreads=6)
      saveRDS(gam.30.ha2.logit, file= paste0(save.pth, "gam_30_ha2_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(gam.30.ha2.pois, gam.30.ha2.logit)
    }
  }
  
  { ####### fit var importance GAMs #######
    save.pth = "E:/tmp/wilt_models/var_importance/"
    #
    # SATURATED MODEL
    #
    #sat.mod = readRDS("E:/tmp/wilt_models/gam_30_h0_logit.Rds")
    #sat.mod |> summary()
    # gam.30.h0.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(aspect, bs="cc") +
    #                               s(bd_x_prof_curv, bs="cs") + s(bd_x_conv_ind, bs="cs") + s(channel_net_s5, bs="cs") + s(conv_ind_x_hillshade, bs="cs") +
    #                               s(topo_wetness, bs="cs") + s(ph_x_sand, bs="cs"),
    #                             data=dat.30, 
    #                             #weights=dat.30$logit_weights, 
    #                             family=binomial(link="logit"), 
    #                             method="fREML", 
    #                             nthreads=6)
    
    # sat
    gam.30.sat = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(channel_net_s5, bs="cs"),
                           data=dat.30,
                           family=binomial(link="logit"), 
                           method="fREML", 
                           nthreads=6)
    saveRDS(gam.30.sat, file = paste0(save.pth, "gam_30_sat.Rds"))
    rm(gam.30.sat)
    
    # null
    mean.resp = mean(dat.30$logit_resp)
    gam.30.null = mgcv::bam(logit_resp ~ rep(mean.resp, nrow(dat.30)),
                            data=dat.30,
                            family=binomial(link="logit"), 
                            method="fREML", 
                            nthreads=6)
    saveRDS(gam.30.null, file = paste0(save.pth, "gam_30_null.Rds"))
    rm(gam.30.null)
    
    # remove x,y
    gam.30.rmv.xy = mgcv::bam(logit_resp ~ s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(channel_net_s5, bs="cs"),
                               data=dat.30,
                               family=binomial(link="logit"), 
                               method="fREML", 
                               nthreads=6)
    saveRDS(gam.30.rmv.xy, file=paste0(save.pth, "gam_30_rmv_xy.Rds"))
    rm(gam.30.rmv.xy)
    cat(crayon::bgRed("\n\ndone 1 of 5\n\n"))
    
    # remove wl2 oak class prob
    gam.30.rmv.wl2 = mgcv::bam(logit_resp ~ s(x, y) + s(elev, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(channel_net_s5, bs="cs"),
                                data=dat.30,
                                family=binomial(link="logit"), 
                                method="fREML", 
                                nthreads=6)
    saveRDS(gam.30.rmv.wl2, file=paste0(save.pth, "gam_30_rmv_wl2.Rds"))
    rm(gam.30.rmv.wl2)
    cat(crayon::bgRed("\n\ndone 2 of 5n\n"))
    
    # remove elev
    gam.30.rmv.elev = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(channel_net_s5, bs="cs"),
                               data=dat.30,
                               family=binomial(link="logit"), 
                               method="fREML", 
                               nthreads=6)
    saveRDS(gam.30.rmv.elev, file=paste0(save.pth, "gam_30_rmv_elev.Rds"))
    rm(gam.30.rmv.elev)
    cat(crayon::bgRed("\n\ndone 3 of 5n\n"))
    
    # remove bd x prof_curv
    gam.30.rmv.bdxprofcurv = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(channel_net_s5, bs="cs"),
                                data=dat.30,
                                family=binomial(link="logit"), 
                                method="fREML", 
                                nthreads=6)
    saveRDS(gam.30.rmv.bdxprofcurv, file=paste0(save.pth, "gam_30_rmv_bdxprofcurv.Rds"))
    rm(gam.30.rmv.bdxprofcurv)
    cat(crayon::bgRed("\n\ndone 4 of 5n\n"))
    
    # remove channel_net_s5
    gam.30.rmv.chan5 = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(bd_x_prof_curv, bs="cs"),
                               data=dat.30,
                               family=binomial(link="logit"), 
                               method="fREML", 
                               nthreads=6)
    saveRDS(gam.30.rmv.chan5, file=paste0(save.pth, "gam_30_rmv_chan5.Rds"))
    rm(gam.30.rmv.chan5)
    cat(crayon::bgRed("\n\ndone 5 of 5 n\n"))
  }
  
  { ##### fit PREC-REC GAMs ######
    { ###### prep data ########
      dat.30 = readRDS("mod_data/model_in/dat_30_weighted.Rds")
    #   { ##### add logit response & weights #####
    #     n.no = sum(dat.30$wilt_train == 0); n.no
    #     n.yes = sum(dat.30$wilt_train > 0); n.yes
    #     dat.30 = dat.30 |> 
    #       mutate(logit_train = if_else(wilt_train > 0, 1, 0),
    #              logit_train = as.integer(logit_train),
    #              weights_raw = if_else(wilt_train == 0, 1/n.no, 1/n.yes),
    #              logit_weights_raw = if_else(wilt_train > 0, wilt_train*weights_raw, weights_raw),
    #              weights = weights_raw/mean(weights_raw),
    #              logit_weights = logit_weights_raw/mean(logit_weights_raw) 
    #       )
    #   }
    # }
    
    save.pth = "E:/tmp/wilt_models/train_test/"
    
    { ###### H0 - spatial + hosts + env ###### 
      {# h0 - poisson
      train.gam.30.h0.pois = mgcv::bam(wilt_train ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(channel_net_s5, bs="cs"),
                                 data=dat.30,
                                 #weights=dat.30$weights,
                                 family=poisson(link="log"), 
                                 method="fREML", nthreads=6)
      
      saveRDS(train.gam.30.h0.pois, file= paste0(save.pth, "train_gam_30_h0_pois.Rds") )
      rm(train.gam.30.h0.pois)
      }
      
      { ###### H0 logit ########
        train.gam.30.h0.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(channel_net_s5, bs="cs"),
                                    data=dat.30,
                                    #weights=dat.30$logit_weights,
                                    family=binomial(link="logit"), 
                                    method="fREML", 
                                    nthreads=6)
        
        saveRDS(train.gam.30.h0.logit, file= paste0(save.pth, "train_gam_30_h0_logit.Rds") )
        cat(crayon::bgRed("\n\nDONE\n\n"))
        rm(train.gam.30.h0.logit)
      }
    }
    
    { ######### HA1 - spatial + hosts ##########
      # ha1 - poisson
      train.gam.30.ha1.pois = mgcv::bam(wilt_train ~ s(x, y) + s(wl2_oakprob_30, bs="cs"),
                                        data=dat.30, 
                                        #weights=dat.30$weights,
                                        family=poisson(link="log"),
                                        method="fREML", 
                                        nthreads=6)
      saveRDS(train.gam.30.ha1.pois, file= paste0(save.pth, "train_gam_30_ha1_pois.Rds") )
      
      # ha1 - logit
      train.gam.30.ha1.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs"),
                                   data=dat.30, 
                                   #weights=dat.30$logit_weights,
                                   family=binomial(link="logit"), 
                                   method="fREML",
                                   nthreads=6)
      saveRDS(train.gam.30.ha1.logit, file= paste0(save.pth, "train_gam_30_ha1_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(train.gam.30.ha1.pois, train.gam.30.ha1.logit)
    }
    
    { ######### HA2 - hosts + env no spatial ##########
      # ha2 - poisson
      train.gam.30.ha2.pois = mgcv::bam(wilt_train ~ s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(channel_net_s5, bs="cs") + s(ph_x_sand, bs="cs"),
                                        data=dat.30,
                                        #weights=dat.30$weights,
                                        family=poisson(link="log"), 
                                        method="fREML", 
                                        nthreads=6)
      saveRDS(train.gam.30.ha2.pois, file= paste0(save.pth, "train_gam_30_ha2_pois.Rds") )
      
      # ha2 - logit
      train.gam.30.ha2.logit = mgcv::bam(logit_resp ~ s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(channel_net_s5, bs="cs") + s(ph_x_sand, bs="cs"),
                                        data=dat.30,
                                        #weights=dat.30$weights,
                                        family=binomial(link="logit"), 
                                        method="fREML", 
                                        nthreads=6)
      saveRDS(train.gam.30.ha2.logit, file= paste0(save.pth, "train_gam_30_ha2_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(train.gam.30.ha2.pois, train.gam.30.ha2.logit)
    }
    
  }
  }
}

#-------------- CPM ---------------------
{ ##### setup #####
  rm(list=ls())
  gc()
  library(tidyverse)
  library(terra)
  library(spatstat)
  library(assertthat)
}

{ ######### Prep CPM data ###########
  
  { ####### Create quadschemes and save ########
    ow.ppp.10 = readRDS("clean_data/ow_ppp_10.Rds")
    ow.ppp.30 = readRDS("clean_data/ow_ppp_30.Rds")
    train.ppp.10 = readRDS("clean_data/train_ppp_10.Rds")
    train.ppp.30 = readRDS("clean_data/train_ppp_30.Rds")
    
    #quads.10 = quadscheme(data=ow.ppp.10, method="grid", eps=10)
    quads.30 = quadscheme(data=ow.ppp.30, method="grid", eps=30)
    #quadtrain.10 = quadscheme(data=train.ppp.10, method="grid", eps=10)
    quadtrain.30 = quadscheme(data=train.ppp.10, method="grid", eps=30)
    
    #saveRDS(quads.10, file="mod_data/kppm/quad_10_10.Rds")
    saveRDS(quads.30, file="mod_data/kppm/quad_30_30.Rds")
    #saveRDS(quadtrain.10, file="mod_data/kppm/quadtrain_10_10.Rds")
    saveRDS(quadtrain.30, file="mod_data/kppm/quadtrain_30_30.Rds")
    
    rm(list=ls())
    gc()
  }
  
  { ####### create im lists ########
    source("functions/mod/tbl_to_im_.R")
    
    dat.mod.30 <- readRDS("mod_data/model_in/dat_30.Rds")
    dat.mod.30.scale = dat.mod.30
    for(i in 4:ncol(dat.mod.30.scale)) {
      dat.mod.30.scale[,i] = scale(dat.mod.30[,i])
    }
    saveRDS(dat.mod.30.scale, "mod_data/model_in/dat_30_scaled.Rds")
    
    vars.in = names(dat.mod.30.scale)[4:ncol(dat.mod.30.scale)]; vars.in
    
    im.list.30 = tbl_to_im_(
      dat.mod.30.scale, 
      vars.in, 
      terra::crs('EPSG:3071'),
      verbose=T, DBG=T)
    saveRDS(im.list.30, file="mod_data/kppm/im_list_30.Rds")
    rm(im.list.30, dat.mod.30, dat.mod.30.scale)
    
    
    dat.mod.10 <- readRDS("mod_data/model_in/dat_10.Rds")
    dat.mod.10.scale = dat.mod.10
    for(i in 4:ncol(dat.mod.10.scale)) {
      dat.mod.10.scale[,i] = scale(dat.mod.10[,i])
    }
    saveRDS(dat.mod.10.scale, file="mod_data/model_in/dat_10_scaled.Rds")
    vars.in = names(dat.mod.10.scale)[4:ncol(dat.mod.10.scale)]; vars.in
    
    im.list.10 = tbl_to_im_(
      dat.mod.10.scale,
      vars.in, 
      terra::crs('EPSG:3071'),
      verbose=T, DBG=T)
    saveRDS(im.list.10, file="mod_data/kppm/im_list_10.Rds")
    rm(dat.mod.10, dat.mod.10.scale, im.list.10)
  }
  
  { ##### check im #######
    im.list.10 = readRDS("mod_data/kppm/im_list_10.Rds")
    im.list.10 |> names()
    im.list.10[[1]] |> plot()
    im.list.10[[13]] |> plot()
    
    im.list.10 = readRDS("mod_data/kppm/im_list_10.Rds")
    
    im.list.30 = readRDS("mod_data/kppm/im_list_30.Rds")
    im.list.30 |> names()
    im.list.30[[2]] |> plot()
    rm(im.list.10, im.list.30)
  }
}

{ #######  Fit CPMs on 30-m data #########
  { ###### Setup #####
    quads.30 = readRDS("mod_data/kppm/quad_30_30.Rds")
    quadtrain.30 = readRDS("mod_data/kppm/quadtrain_30_30.Rds")
    im.list.30 = readRDS("mod_data/kppm/im_list_30.Rds"); names(im.list.30)
  }
  
  { ####### Hypothesis models ########
    save.pth = "E:/tmp/wilt_models/"
    
    kppm.30.thom.pq = kppm(
      X = quads.30, 
      #trend = ~.,
      trend = ~ wl2_oakprob_30 + elev + bd_x_prof_curv + channel_net_s5,
      clusters="Thomas",
      penalised = T,
      data= im.list.30,
      method = "palm",
      improve.type = "quasi")
    #summary(kppm.30.thom.pq)
    saveRDS(kppm.30.thom.pq, file=paste0(save.pth, "kppm_30_thom_pq.Rds"))
    rm(kppm.30.thom.pq)
    
    kppm.30.thom.wclik = kppm(
      X = quads.30, 
      #trend = ~.,
      trend = ~ elev + bd_x_prof_curv + channel_net_s5,
      clusters="Thomas",
      penalised = T,
      #use.gam = T,
      data= im.list.30,
      method = "clik2",
      improve.type = "wclik1")
    summary(kppm.30.thom.wclik)
    saveRDS(kppm.30.thom.wclik, file=paste0(save.pth, "kppm_30_thom_wclik.Rds"))
    rm(kppm.30.thom.wclik)
    
    kppm.30.mat.pq = kppm(
      X = quads.30, 
      #trend= ~.,
      trend = ~ elev + bd_x_prof_curv + channel_net_s5,
      clusters="MatClust",
      penalised = T,
      data=im.list.30,
      method = "palm",
      improve.type = "quasi")
    summary(kppm.30.mat.pq)
    saveRDS(kppm.30.mat.pq, file=paste0(save.pth, "kppm_30_mat_pq.Rds"))
    rm(kppm.30.mat.pq)
    
    kppm.30.mat.wclik1 = kppm(
      X = quads.30, 
      #trend= ~.,
      trend = ~ elev + bd_x_prof_curv + channel_net_s5,
      clusters="MatClust",
      penalised = T,
      data=im.list.30,
      method = "clik2",
      improve.type = "wclik1")
    summary(kppm.30.mat.wclik1)
    saveRDS(kppm.30.mat.wclik1, file=paste0(save.pth, "kppm_30_mat_wclik1.Rds"))
    rm(kppm.30.mat.wclik1)
    
    kppm.30.cauch.pq = kppm(
      X = quads.30, 
      #trend = ~.,
      trend = ~ elev + bd_x_prof_curv + channel_net_s5,
      clusters="Cauchy",
      penalised = T,
      data=im.list.30,
      method = "palm",
      improve.type = "quasi")
    summary(kppm.30.cauch.pq)
    saveRDS(kppm.30.cauch.pq, file=paste0(save.pth, "kppm_30_cauch_pq.Rds"))
    rm(kppm.30.cauch.pq)
    
    kppm.30.cauch.wclik1 = kppm(
      X = quads.30, 
      #trend = ~.,
      trend = ~ elev + bd_x_prof_curv + channel_net_s5,
      clusters="Cauchy",
      penalised = T,
      data=im.list.30,
      method = "clik2",
      improve.type = "wclik1")
    summary(kppm.30.cauch.wclik1)
    saveRDS(kppm.30.cauch.wclik1, file=paste0(save.pth, "kppm_30_cauch_wclik1.Rds"))
    rm(kppm.30.cauch.wclik1)
    
    # kppm.30.vargam = kppm(
    #   X = quads.30, 
    #   trend = ~.,
    #   clusters="VarGamma",
    #   penalised = T,
    #   data=im.list.30,
    #   #subset=ow.ppp.30$window,
    #   method = "palm",
    #   improve.type = "quasi")
    # summary(kppm.30.vargam)
    # saveRDS(kppm.30.vargam, file=paste0(save.pth, "kppm_30_vargam.Rds"))
    
    { ####### alternative #######
      save.pth = "E:/tmp/wilt_models/"
      kppm.30.thom.pq.ha1 = kppm(
        X = quads.30, 
        #trend = ~.,
        trend = ~ wl2_oakprob_30,
        clusters="Thomas",
        penalised = T,
        data= im.list.30,
        method = "palm",
        improve.type = "quasi")
      saveRDS(kppm.30.thom.pq.ha1, file=paste0(save.pth, "kppm_30_thom_pq_ha1.Rds"))
      
      kppm.30.thom.wclik1.ha1 = kppm(
        X = quads.30, 
        #trend = ~.,
        trend = ~ wl2_oakprob_30,
        clusters="Thomas",
        penalised = T,
        data= im.list.30,
        method = "clik2",
        improve.type = "wclik1")
      saveRDS(kppm.30.thom.wclik1.ha1, file=paste0(save.pth, "kppm_30_thom_wclik1_ha1.Rds"))
      rm(kppm.30.thom.pq.ha1, kppm.30.thom.wclik1.ha1)
      cat(crayon::bgRed("\n\n1 of 3"))
      
      
      kppm.30.mat.pq.ha1 = kppm(
        X = quads.30, 
        #trend= ~.,
        trend = ~ wl2_oakprob_30,
        clusters="MatClust",
        penalised = T,
        data=im.list.30,
        method = "palm",
        improve.type = "quasi")
      saveRDS(kppm.30.mat.pq.ha1, file=paste0(save.pth, "kppm_30_mat_pq_ha1.Rds"))
      
      kppm.30.mat.wclik1.ha1 = kppm(
        X = quads.30, 
        #trend= ~.,
        trend = ~ wl2_oakprob_30,
        clusters="MatClust",
        penalised = T,
        data=im.list.30,
        method = "clik2",
        improve.type = "wclik1")
      saveRDS(kppm.30.mat.wclik1.ha1, file=paste0(save.pth, "kppm_30_mat_wclik1_ha1.Rds"))
      rm(kppm.30.mat.pq.ha1, kppm.30.mat.wclik1.ha1)
      
      kppm.30.cauch.pq.ha1 = kppm(
        X = quads.30, 
        #trend = ~.,
        trend = ~ wl2_oakprob_30,
        clusters="Cauchy",
        penalised = T,
        data=im.list.30,
        method = "palm",
        improve.type = "quasi")
      saveRDS(kppm.30.cauch.pq.ha1, file=paste0(save.pth, "kppm_30_cauch_pq_ha1.Rds"))
      
      kppm.30.cauch.wclik1.ha1 = kppm(
        X = quads.30, 
        #trend = ~.,
        trend = ~ wl2_oakprob_30,
        clusters="Cauchy",
        penalised = T,
        data=im.list.30,
        method = "clik2",
        improve.type = "wclik1")
      saveRDS(kppm.30.cauch.wclik1.ha1, file=paste0(save.pth, "kppm_30_cauch_wclik1_ha1.Rds"))
      
      rm( kppm.30.cauch.pq.ha1, kppm.30.cauch.wclik1.ha1 )
      
      cat(crayon::bgRed("\n\n3 of 3"))
    }
  }
  
  { ########## Prec + Rec models #########
    save.pth = "E:/tmp/wilt_models/train_test/"
    
    { #### HA2 ####
      tkppm.ha2.pq.thom = kppm(
        X = quadtrain.30, 
        trend = ~ wl2_oakprob_30,
        clusters="Thomas",
        penalised = T,
        data= im.list.30,
        method = "palm",
        improve.type = "quasi")
      #summary(train.kppm.30.thom)
      saveRDS(tkppm.ha2.pq.thom, file=paste0(save.pth, "tkppm_ha2_pq_thom.Rds"))
      rm(tkppm.ha2.pq.thom)
      
      tkppm.ha2.wc.thom = kppm(
        X = quadtrain.30, 
        trend = ~ wl2_oakprob_30,
        clusters="Thomas",
        penalised = T,
        data= im.list.30,
        method = "clik2",
        improve.type = "wclik1")
      #summary(train.kppm.30.thom)
      saveRDS(tkppm.ha2.wc.thom, file=paste0(save.pth, "tkppm_ha2_wc_thom.Rds"))
      rm(tkppm.ha2.wc.thom)
      
      tkppm.ha2.pq.cauch = kppm(
        X = quadtrain.30, 
        trend = ~ wl2_oakprob_30,
        clusters="Cauchy",
        penalised = T,
        data= im.list.30,
        method = "palm",
        improve.type = "quasi")
      #summary(train.kppm.30.thom)
      saveRDS(tkppm.ha2.pq.cauch, file=paste0(save.pth, "tkppm_ha2_pq_cauch.Rds"))
      rm(tkppm.ha2.pq.cauch)
      
      tkppm.ha2.wc.cauch = kppm(
        X = quadtrain.30, 
        trend = ~ wl2_oakprob_30,
        clusters="Cauchy",
        penalised = T,
        data= im.list.30,
        method = "clik2",
        improve.type = "wclik1")
      #summary(train.kppm.30.thom)
      saveRDS(tkppm.ha2.wc.cauch, file=paste0(save.pth, "tkppm_ha2_wc_cauch.Rds"))
      rm(tkppm.ha2.wc.cauch)
    }
    train.kppm.30.pq.thom = kppm(
      X = quadtrain.30, 
      trend = ~ wl2_oakprob_30 + elev + bd_x_prof_curv + channel_net_s5,
      clusters="Thomas",
      penalised = T,
      data= im.list.30,
      method = "palm",
      improve.type = "quasi")
    #summary(train.kppm.30.thom)
    saveRDS(train.kppm.30.pq.thom, file=paste0(save.pth, "train_kppm_30_pq_thom.Rds"))
    rm(train.kppm.30.pq.thom)
    
    train.kppm.30.wclik1.thom = kppm(
      X = quadtrain.30, 
      trend = ~ wl2_oakprob_30 + elev + bd_x_prof_curv + channel_net_s5,
      clusters="Thomas",
      penalised = T,
      data= im.list.30,
      method = "clik2",
      improve.type = "wclik1")
    #summary(train.kppm.30.thom)
    saveRDS(train.kppm.30.wclik1.thom, file=paste0(save.pth, "train_kppm_30_wclik1_thom.Rds"))
    rm(train.kppm.30.wclik1.thom)
    
    train.kppm.30.pq.mat = kppm(
      X = quadtrain.30, 
      trend = ~ elev + bd_x_prof_curv + channel_net_s5,
      clusters="MatClust",
      penalised = T,
      data=im.list.30,
      method = "palm",
      improve.type = "quasi")
    #summary(train.kppm.30.mat)
    saveRDS(train.kppm.30.pq.mat, file=paste0(save.pth, "train_kppm_30_pq_mat.Rds"))
    rm(train.kppm.30.pq.mat)
    
    train.kppm.30.wclik1.mat = kppm(
      X = quadtrain.30, 
      trend = ~ elev + bd_x_prof_curv + channel_net_s5,
      clusters="MatClust",
      penalised = T,
      data=im.list.30,
      method = "clik2",
      improve.type = "wclik1")
    #summary(train.kppm.30.mat)
    saveRDS(train.kppm.30.wclik1.mat, file=paste0(save.pth, "train_kppm_30_wclik1_mat.Rds"))
    rm(train.kppm.30.wclik1.mat)
    
    train.kppm.30.pq.cauch = kppm(
      X = quadtrain.30, 
      trend = ~ elev + bd_x_prof_curv + channel_net_s5,
      clusters="Cauchy",
      penalised = T,
      data=im.list.30,
      method = "palm",
      improve.type = "quasi")
    #summary(train.kppm.30.cauch)
    saveRDS(train.kppm.30.pq.cauch, file=paste0(save.pth, "train_kppm_30_pq_cauch.Rds"))
    rm( train.kppm.30.pq.cauch )
    
    train.kppm.30.wclik1.cauch = kppm(
      X = quadtrain.30, 
      trend = ~ elev + bd_x_prof_curv + channel_net_s5,
      clusters="Cauchy",
      penalised = T,
      data=im.list.30,
      method = "clik2",
      improve.type = "wclik1")
    #summary(train.kppm.30.cauch)
    saveRDS(train.kppm.30.wclik1.cauch, file=paste0(save.pth, "train_kppm_30_wclik1_cauch.Rds"))
    rm( train.kppm.30.wclik1.cauch )
    # train.kppm.30.cauch.gam = kppm(
    #   X = quadtrain.30, 
    #   trend = ~ elev + aspect + bd_x_prof_curv + channel_net_s5,
    #   clusters="Cauchy",
    #   data=im.list.30,
    #   use.gam=T,
    #   method = "clik2",
    #   improve.type = "wclik1")
    # train.kppm.30.cauch.gam |> predict() |> plot()
    # saveRDS(train.kppm.30.cauch.gam, file=paste0(save.pth, "train_kppm_30_cauch_gam.Rds"))
  }
  
  rm(quads.30, quadtrain.30, im.list.30)
}

{ ####### Fit 10-m KPPM #########
  { ##### setup #####
    {
      rm(list=ls())
      gc()
    }
    quads.10 = readRDS("mod_data/kppm/quad_10_10.Rds")
    quadtrain.10 = readRDS("mod_data/kppm/quadtrain_10_10.Rds")
    im.list.10 = readRDS("mod_data/kppm/im_list_10.Rds"); names(im.list.10)
    save.pth = "E:/tmp/wilt_models/"
  }
  
  { ####### hypothesis models ##########
    kppm.10.thom = kppm(
      X = quads.10, 
      #trend = ~.,
      trend = ~ elev + channel_dist_s5 + val_depth + bd_x_hillshade + bd_x_hillshade + bd_x_prof_curv,
      clusters="Thomas",
      penalised = T,
      data= im.list.10,
      method = "palm",
      improve.type = "quasi")
    #summary(kppm.10.thom)
    saveRDS(kppm.10.thom, file=paste0(save.pth, "kppm_10_thom.Rds"))
    rm( kppm.10.thom )
    
    kppm.10.mat = kppm(
      X = quads.10, 
      #trend= ~.,
      trend = ~ elev + channel_dist_s5 + val_depth + bd_x_hillshade + bd_x_hillshade + bd_x_prof_curv,
      clusters="MatClust",
      penalised = T,
      data=im.list.10,
      method = "palm",
      improve.type = "quasi")
    #summary(kppm.10.mat)
    saveRDS(kppm.10.mat, file=paste0(save.pth, "kppm_10_mat.Rds"))
    rm( kppm.10.mat )
    
    kppm.10.cauch = kppm(
      X = quads.10, 
      #trend = ~.,
      trend = ~ elev + channel_dist_s5 + val_depth + bd_x_hillshade + bd_x_hillshade + bd_x_prof_curv,
      clusters="Cauchy",
      penalised = T,
      data=im.list.10,
      method = "palm",
      improve.type = "quasi")
    #summary(kppm.10.cauch)
    saveRDS(kppm.10.cauch, file=paste0(save.pth, "kppm_10_cauch.Rds"))
    rm( kppm.10.cauch )
    cat(crayon::bgRed("\n\n DONE \n\n"))
    
    { #### ha1 ####
      kppm.10.thom.ha1 = kppm(
        X = quads.10, 
        #trend = ~.,
        trend = ~ wl2_oakprob_10,
        clusters="Thomas",
        penalised = T,
        data= im.list.10,
        method = "palm",
        improve.type = "quasi")
      #summary(kppm.10.thom)
      saveRDS(kppm.10.thom.ha1, file=paste0(save.pth, "kppm_10_thom_ha1.Rds"))
      rm( kppm.10.thom.ha1 )
      cat(crayon::bgRed("\n\n1 of 3"))
      
      kppm.10.mat.ha1 = kppm(
        X = quads.10, 
        #trend= ~.,
        trend = ~ wl2_oakprob_10,
        clusters="MatClust",
        penalised = T,
        data=im.list.10,
        method = "palm",
        improve.type = "quasi")
      #summary(kppm.10.mat)
      saveRDS(kppm.10.mat.ha1, file=paste0(save.pth, "kppm_10_mat_ha1.Rds"))
      rm( kppm.10.mat.ha1 )
      cat(crayon::bgRed("\n\n2 of 3"))
      
      kppm.10.cauch.ha1 = kppm(
        X = quads.10, 
        #trend = ~.,
        trend = ~ wl2_oakprob_10,
        clusters="Cauchy",
        penalised = T,
        data=im.list.10,
        method = "palm",
        improve.type = "quasi")
      #summary(kppm.10.cauch)
      saveRDS(kppm.10.cauch.ha1, file=paste0(save.pth, "kppm_10_cauch_ha1.Rds"))
      rm( kppm.10.cauch.ha1 )
      cat(crayon::bgRed("\n\n3 of 3"))
    }
  }
  
  {
    save.pth = "E:/tmp/wilt_models/train_test/"
    # kppm.10.full = kppm(
    #   X = quadtrain.10,
    #   trend = ~.,
    #   clusters="Thomas",
    #   data=im.list.10,
    #   rmax=1000,
    #   method = "palm",
    #   improve.type = "quasi")
    # summary(kppm.10.full)
    
    # kppm.10.subs = kppm(
    #   X = quadtrain.10,
    #   trend = ~ elev + channel_dist_s5 + ,
    #   clusters="Thomas",
    #   data=im.list.10,
    #   rmax=1000,
    #   method = "palm",
    #   improve.type = "quasi")
    
    train.kppm.10.thom = kppm(
      X = quadtrain.10, 
      trend = ~ elev + channel_dist_s5 + val_depth + bd_x_hillshade + bd_x_hillshade + bd_x_prof_curv,
      clusters="Thomas",
      penalised = T,
      data= im.list.10,
      method = "palm",
      improve.type = "quasi")
    #summary(train.kppm.10.thom)
    saveRDS(train.kppm.10.thom, file=paste0(save.pth, "train_kppm_10_thom.Rds"))
    rm(train.kppm.10.thom)
    
    train.kppm.10.mat = kppm(
      X = quadtrain.10, 
      trend = ~ elev + channel_dist_s5 + val_depth + bd_x_hillshade + bd_x_hillshade + bd_x_prof_curv,
      clusters="MatClust",
      penalised = T,
      data=im.list.10,
      method = "palm",
      improve.type = "quasi")
    #summary(train.kppm.10.mat)
    saveRDS(train.kppm.10.mat, file=paste0(save.pth, "train_kppm_10_mat.Rds"))
    rm(train.kppm.10.mat)
    
    train.kppm.10.cauch = kppm(
      X = quadtrain.10, 
      trend = ~ elev + channel_dist_s5 + val_depth + bd_x_hillshade + bd_x_hillshade + bd_x_prof_curv,
      clusters="Cauchy",
      penalised = T,
      data=im.list.10,
      method = "palm",
      improve.type = "quasi")
    #summary(train.kppm.10.cauch)
    saveRDS(train.kppm.10.cauch, file=paste0(save.pth, "train_kppm_10_cauch.Rds"))
    rm(train.kppm.10.cauch)
  }
}

#-------------- compare ---------------------

{ ####### CPM ######
  { #### setup ###
    rm(list=ls())
    gc()
    library(tidyverse)
    library(spatstat)
    { ###### compare cpm function #######
      compare_cpm <- function(cpm_sat, cpm_subs, DBG=F) {
        delta.dof <- as.integer(attr(logLik(cpm_sat), "df")) - as.integer(attr(logLik(cpm_subs), "df")) # delta degrees of freedom in sat (large) vs subs (smaller)
        if(DBG) {cat(paste0("\n\n\ndelta.dof = ", as.integer(delta.dof)) )}
        chi.stat <- -2*(as.numeric(logLik(cpm_sat)) - as.numeric(logLik(cpm_subs)))
        if(DBG) {cat(paste0("\nchi-stat = ", chi.stat))}
        p.val <- pchisq(chi.stat, df=delta.dof, lower.tail=F)
        if(DBG) cat(paste0("\np = ", p.val))
        return(p.val)
      }
    }
    mod.pth = "E:/tmp/wilt_models/"
  }
  
  { #### hypothesis ####
    { ##### 10-m AIC & Chi-sq ####
      mod.pth = "E:/tmp/wilt_models/"
      
      { ##### cauch #####
        kppm.10.cauch.sat = readRDS( paste0(mod.pth, "kppm_10_cauch.Rds" )  ) 
        kppm.10.cauch.ha1 = readRDS(paste0(mod.pth, "kppm_10_cauch_ha1.Rds"))
        extractAIC(kppm.10.cauch.sat) # 2305630
        extractAIC(kppm.10.cauch.ha1) # 2284975 ##2305630 - 2284975
        compare_cpm(cpm_sat = kppm.10.cauch.sat, cpm_subs = kppm.10.cauch.ha1, DBG=T) # 0
        rm(kppm.10.cauch.sat, kppm.10.cauch.ha1)
        gc()
      }
      
      { ##### thomas ####
        kppm.10.thom.sat = readRDS( paste0(mod.pth, "kppm_10_thom.Rds" )  )
        kppm.10.thom.ha1 = readRDS(paste0(mod.pth, "kppm_10_thom_ha1.Rds"))
        extractAIC(kppm.10.thom.sat)
        extractAIC(kppm.10.thom.ha1)
        compare_cpm(cpm_ha = kppm.10.thom.ha1, cpm_h0=kppm.10.thom.sat)
        rm(kppm.10.thom.sat, kppm.10.thom.ha1)
        gc()
      }
      
      { ##### matern ####
        kppm.10.mat.sat = readRDS( paste0(mod.pth, "kppm_10_mat.Rds" )  ) # 'log Lik.' -1152230 (df=6) $clfit$value -1199838
        kppm.10.mat.ha1 = readRDS(paste0(mod.pth, "kppm_10_mat_ha1.Rds")) # 'log Lik.' -1142463 (df=2) $clfit$value [1] -1181988
        extractAIC(kppm.10.mat.sat) # 2304472
        extractAIC(kppm.10.mat.ha1) # 2284930
        compare_cpm(cpm_subs = kppm.10.mat.ha1, cpm_sat = kppm.10.mat.sat) # 0
        rm(kppm.10.mat.sat, kppm.10.mat.ha1)
        gc()
      }
      
    }
      
    { ##### 30-m AIC & Chi-sq ####
      mod.pth = "E:/tmp/wilt_models/"
      { #### cauch ###
        kppm.30.cauch.sat = readRDS( paste0(mod.pth, "kppm_30_cauch.Rds" )  ) # 'log Lik.' -1199838 (df=5)
        kppm.30.cauch.ha1 = readRDS( paste0(mod.pth, "kppm_30_cauch_ha1.Rds") ) # 'log Lik.' -1181988 (df=2)
        (extractAIC(kppm.30.cauch.sat)) # 2399686
        (extractAIC(kppm.30.cauch.ha1)) # 2363981
        (compare_cpm(cpm_subs = kppm.30.cauch.ha1, cpm_sat=kppm.30.cauch.sat, DBG=T))
      }
      
      { #### thom ###
        kppm.30.thom.sat = readRDS( paste0(mod.pth, "kppm_30_thom.Rds" )  )
        kppm.30.mat.ha1 = readRDS( paste0(mod.pth, "") )
        (extractAIC(kppm.30.thom.sat))
        (extractAIC(kppm.30.thom.ha1))
        (compare_cpm(cpm_ha = kppm.30.thom.ha1, cpm_h0=kppm.30.thom.sat))
      }
      
      { ### mat ####
        kppm.30.mat.sat = readRDS( paste0(mod.pth, "kppm_30_mat.Rds" )  )
        kppm.30.mat.ha1 = readRDS( paste0(mod.pth, "k") )
        (extractAIC(kppm.30.mat.sat))
        (extractAIC(kppm.30.mat.ha1))
        (compare_cpm(cpm_ha = kppm.30.mat.ha1, cpm_h0=kppm.30.mat.sat))
      }
    
    }
  }
  
  { ##### var importance ####
    { #### 10-m ####
      
    }
    
    { #### 30-m ####
      
    }
  }
  
}

{ ###### GAMs ######
  { ##### setup ######
    rm(list=ls())
    gc()
    library(tidyverse)
    library(mgcv)
    load_mod_ <- function(str.in) {
      mod.pth = "E:/tmp/wilt_models/"
      return(readRDS(paste0(mod.pth, str.in, ".Rds")))
    }
  }
  
  { ###### 10-m ######
    {#### load models for AIC comparison ###### 
      hypothesis.models = list.files("E:/tmp/wilt_models/") 
      
      gam.10.h0.logit = load_mod_("gam_10_h0_logit")
      gam.10.h0.pois = load_mod_("gam_10_h0_pois")
      
      gam.10.ha1.logit = load_mod_("gam_10_ha1_logit")
      gam.10.ha1.pois = load_mod_("gam_10_ha1_pois")
      
      gam.10.ha2.logit = load_mod_("gam_10_ha2_logit")
      gam.10.ha2.pois = load_mod_("gam_10_ha2_pois")
    }
    
    { ####### summaries #######
      gam.10.h0.logit |> summary() # wl2, elev, val_dep, bd x hillshade - possible aspect / plan_curv
      gam.10.ha2.logit |> summary() # wl2, elev, val_dep, bd x hillshade, channel_dist_s5, aspect - possible plan_curv
    }
    
    { ######## AIC comparison ##########
      AIC(gam.10.h0.pois) # 8953.762
      AIC(gam.10.ha1.pois) # 8989.706
      AIC(gam.10.ha2.pois) # 9075.332
      
      AIC(gam.10.h0.logit) # 8308.182
      AIC(gam.10.ha1.logit) # 8337.225
      AIC(gam.10.ha2.logit) # 8416.19
    }
    
    { ######## hypothesis test comparison ###########
      anova.gam(gam.10.h0.logit, gam.10.ha1.logit, test="Chisq")
      #anova.gam(gam.10.ha1.logit, gam.10.h0.logit, test="Chisq")
      
      anova.gam(gam.10.h0.logit, gam.10.ha2.logit, test="Chisq")
      #anova.gam(gam.10.ha2.logit, gam.10.h0.logit, test="Chisq")
      
      anova.gam(gam.10.h0.pois, gam.10.ha1.pois, test="Chisq")
      anova.gam(gam.10.h0.pois, gam.10.ha2.pois, test="Chisq")
    }
    
    { ######## variable importance #########
      { ##### setup #####
        dev_calc_ = function(mod_subs, mod_sat, mod_null) {
          dev.ratio = (deviance(mod_subs) - deviance(mod_sat)) / (deviance(mod_null))
        }
        mod.pth = "E:/tmp/wilt_models/var_importance/"
        { #### h0 ####
          { ##### load h0 10-m logit gams #######
            mod.null = readRDS( paste0(mod.pth, "gam_10_null.Rds") )
            mod.sat = readRDS( paste0(mod.pth, "gam_10_sat.Rds") )
            mod.rmv.xy = readRDS( paste0(mod.pth, "gam_10_rmv_xy.Rds") )
            mod.rmv.wl2 = readRDS( paste0(mod.pth, "gam_10_rmv_wl2.Rds") )
            mod.rmv.valdep = readRDS( paste0(mod.pth, "gam_10_rmv_valdep.Rds") )
            mod.rmv.elev = readRDS( paste0(mod.pth, "gam_10_rmv_elev.Rds") )
            mod.rmv.bdxhs = readRDS( paste0(mod.pth, "gam_10_rmv_bdxhillshade.Rds") )
          }
          { ####### calc, tbl, plot ########
            var.import = c(
              "xy" = dev_calc_(mod_subs = mod.rmv.xy, mod_sat = mod.sat, mod_null = mod.null),
              "wiscland2" = dev_calc_(mod_subs = mod.rmv.wl2, mod_sat = mod.sat, mod_null = mod.null),
              "val depth" = dev_calc_(mod_subs = mod.rmv.valdep, mod_sat = mod.sat, mod_null = mod.null),
              "elev" = dev_calc_(mod_subs = mod.rmv.elev, mod_sat = mod.sat, mod_null = mod.null),
              "bd x hillshade" = dev_calc_(mod_subs = mod.rmv.bdxhs, mod_sat = mod.sat, mod_null = mod.null)
            )
            var.import.tbl = tibble(
              var = names(var.import), 
              val = as.numeric(var.import)) |> 
              mutate(
                relative_importance = val/max(val),
                var = as.factor(var)
              )
            
            saveRDS(var.import.tbl, file="plt_data/var_import_tbl_10.Rds")
            
            var.import.tbl = readRDS("plt_data/var_import_tbl_10.Rds") |> 
              mutate(prop_import = val/sum(val))
            var.import.tbl |> 
              ggplot(aes(x=var, y=relative_importance, fill=var)) + geom_col()
            
            gam.10.logit.h0 = readRDS("E:/tmp/wilt_models/gam_10_h0_logit.Rds")
            gam.10.logit.h0 |> summary()
          }
        }
        
        { #### ha2 ####
          { ###### load ha2 10-m logit gams ######
            mod.ha.null = readRDS( paste0(mod.pth, "ha2_10_null.Rds") )
            mod.ha.sat = readRDS( paste0(mod.pth, "ha2_10_sat.Rds") )
            mod.ha.rmv.wl2 = readRDS( paste0(mod.pth, "ha2_10_rm_wl2.Rds") )
            mod.ha.rmv.valdep = readRDS( paste0(mod.pth, "ha2_10_rm_valdep.Rds") )
            mod.ha.rmv.elev = readRDS( paste0(mod.pth, "ha2_10_rm_elev.Rds") )
            mod.ha.rmv.cn = readRDS(paste0(mod.pth, "ha2_10_rm_cn.Rds"))
            mod.ha.rmv.bdxhs = readRDS( paste0(mod.pth, "ha2_10_rm_bdxhs.Rds") )
            mod.ha.rmv.asp = readRDS( paste0(mod.pth, "ha2_10_rm_asp.Rds") )
          }
          
          {
            var.import = c(
              "wl2" = dev_calc_(mod_subs = mod.ha.rmv.wl2, mod_sat = mod.ha.sat, mod_null = mod.ha.null),
              "valdep" = dev_calc_(mod_subs = mod.ha.rmv.valdep, mod_sat = mod.ha.sat, mod_null = mod.ha.null),
              "elev" = dev_calc_(mod_subs = mod.ha.rmv.elev, mod_sat = mod.ha.sat, mod_null = mod.ha.null),
              "chan_net" = dev_calc_(mod_subs = mod.ha.rmv.cn, mod_sat = mod.ha.sat, mod_null = mod.ha.null),
              "bd x hillshade" = dev_calc_(mod_subs = mod.ha.rmv.bdxhs, mod_sat = mod.ha.sat, mod_null = mod.ha.null),
              "aspect" = dev_calc_(mod_subs = mod.ha.rmv.asp, mod_sat = mod.ha.sat, mod_null = mod.ha.null)
            )
            
            var.import.tbl = tibble(
              var = names(var.import), 
              val = as.numeric(var.import)) |> 
              mutate(
                relative_importance = val/sum(val),
                var = as.factor(var)
              )
            saveRDS(var.import.tbl, file="plt_data/var_import_ha2_10.Rds")
          }
        }
        
      }
      
    
  }
  
  { ###### 30-m ######
    { ##### load models ###### 
      gam.30.h0.logit = load_mod_("gam_30_h0_logit")
      gam.30.h0.pois = load_mod_("gam_30_h0_pois")
      
      gam.30.ha1.logit = load_mod_("gam_30_ha1_logit")
      gam.30.ha1.pois = load_mod_("gam_30_ha1_pois")
      
      gam.30.ha2.logit = load_mod_("gam_30_ha2_logit")
      gam.30.ha2.pois = load_mod_("gam_30_ha2_pois")
    }
    
    { ###### summaries ######
      gam.30.h0.logit |> summary() # wl2, elev, bd x prof_curv, channel_net_s5
      gam.30.ha2.logit |> summary() # wl2, elev, bd x prof_curv, channel_net_s5, ph x sand
    }
    
    { ##### AIC ######
      AIC(gam.30.h0.logit) # 5906.83
      AIC(gam.30.h0.pois) # 6428.835
      
      AIC(gam.30.ha1.logit) # 5936.182
      AIC(gam.30.ha1.pois) # 6473.728
      
      AIC(gam.30.ha2.logit) # 6045.435
      AIC(gam.30.ha2.pois) # 6579.864
    }
    
    { ######## hypothesis test comparison ###########
      anova.gam(gam.30.h0.logit, gam.30.ha1.logit, test="Chisq")
      anova.gam(gam.30.h0.logit, gam.30.ha2.logit, test="Chisq")
      
      anova.gam(gam.30.h0.pois, gam.30.ha1.pois, test="Chisq")
      anova.gam(gam.30.h0.pois, gam.30.ha2.pois, test="Chisq")
    }
    
    { ######## variable importance #########
      dev_calc_ = function(mod_subs, mod_sat, mod_null) {
        dev.ratio = (deviance(mod_subs) - deviance(mod_sat)) / (deviance(mod_null))
      }
      
      mod.pth = "E:/tmp/wilt_models/var_importance/"
      
      mod.null = readRDS(paste0(mod.pth, "gam_30_null.Rds"))
      mod.sat = readRDS( paste0(mod.pth, "gam_30_sat.Rds") )
      
      mod.rmv.xy = readRDS(paste0(mod.pth, "gam_30_rmv_xy.Rds"))
      mod.rmv.wl2 = readRDS(paste0(mod.pth, "gam_30_rmv_wl2.Rds"))
      mod.rmv.elev = readRDS(paste0(mod.pth, "gam_30_rmv_elev.Rds"))
      mod.rmv.chan5 = readRDS(paste0(mod.pth, "gam_30_rmv_chan5.Rds"))
      mod.rmv.bdxprof = readRDS(paste0(mod.pth, "gam_30_rmv_bdxprofcurv.Rds"))
      
      var.import = c(
        "x, y" = dev_calc_(mod_subs=mod.rmv.xy, mod_sat=mod.sat, mod_null=mod.null),
        "wiscland2_oakprob" = dev_calc_(mod_subs=mod.rmv.wl2, mod_sat=mod.sat, mod_null=mod.null),
        "elevation" = dev_calc_(mod_subs=mod.rmv.elev, mod_sat=mod.sat, mod_null=mod.null),
        "channel_network_distance_s5" = dev_calc_(mod_subs=mod.rmv.chan5, mod_sat=mod.sat, mod_null=mod.null),
        "bd x prof_curv" = dev_calc_(mod_subs=mod.rmv.bdxprof, mod_sat=mod.sat, mod_null=mod.null))
      var.import.tbl = tibble(var = names(var.import), val=as.numeric(var.import)) |> 
        mutate(relative_importance = val/max(val),
               var = as.factor(var))
      saveRDS(var.import.tbl, file="plt_data/var_import_tbl_30.Rds")
      
      
      var.import.tbl = readRDS("plt_data/var_import_tbl_30.Rds") |> 
        mutate(rel_import = val / sum(val))
      var.import.tbl |>  #**plot*
        ggplot(aes(x=var, y=relative_importance, fill=var)) + geom_col()
      
      gam.30.logit.h0 = readRDS("E:/tmp/wilt_models/gam_30_h0_logit.Rds")
      
      gam.30.logit.h0 |> summary()
    }
    
  }
  
  
  
  
  
  {
    
  }
  
  }
}

{ ######## Precision + Recall on train/test split ########
  { #### init ####
    rm(list=ls())
    gc()
    library(terra)
    library(spatstat)
    library(tidyverse)
  } |> suppressPackageStartupMessages()
  
  { #### add train-test partition to dat.30 ####
    { ###### 30-m #######
      # covariate data
      dat.30 = readRDS("mod_data/model_in/dat_30.Rds")
      dat.rast = terra::rast("clean_data/joindat_30_800_50.tif")
      ow.pts = terra::vect("clean_data/ow_pts_clean.shp")
      nrow(ow.pts)
      test.ind = which(ow.pts$year >= 2018); length(test.ind)
      train.pts = ow.pts[-test.ind, ]
      test.pts = ow.pts[test.ind, ]
      
      train.rast = rasterize(train.pts, dat.rast, field=1, sum=T); names(train.rast) = "wilt_train"
      test.rast = rasterize(test.pts, dat.rast, field=1, sum=T); names(test.rast) = "wilt_test"
      comb.rast = c(train.rast, test.rast)
      names(comb.rast)
      terra::writeRaster(comb.rast, filename="mod_data/compare/wilt_split_rast.tif", overwrite=T)
      rm(dat.rast, ow.pts, train.pts, test.pts, train.rast, test.rast)
      comb.tbl = comb.rast |> 
        as.data.frame(xy=T)
      dat.30 = left_join(dat.30, comb.tbl, by=c("x", "y"))  |> 
        mutate(wilt_train = if_else(is.na(wilt_train), 0, wilt_train),
               wilt_test = if_else(is.na(wilt_test), 0, wilt_test) )
      saveRDS(dat.30, file="mod_data/compare/wilt_split_tbl.Rds")
      sum(dat.30$wilt_train)
      sum(dat.30$wilt_test)
    }
    
    { ###### 10-m #######
      # covariate data
      dat.10 = readRDS("mod_data/model_in/dat_10.Rds")
      dat.rast = terra::rast("clean_data/joindat_10_800_50.tif")
      ow.pts = terra::vect("mid_data/wilt/ow_pts_comb.shp")
      nrow(ow.pts)
      test.ind = which(ow.pts$year >= 2018); length(test.ind)
      train.pts = ow.pts[-test.ind, ]
      test.pts = ow.pts[test.ind, ]
      
      train.rast = rasterize(train.pts, dat.rast, field=1, sum=T); names(train.rast) = "wilt_train"
      test.rast = rasterize(test.pts, dat.rast, field=1, sum=T); names(test.rast) = "wilt_test"
      comb.rast = c(train.rast, test.rast)
      names(comb.rast)
      terra::writeRaster(comb.rast, filename="mod_data/compare/wilt_split_rast_10.tif", overwrite=T)
      rm(dat.rast, ow.pts, train.pts, test.pts, train.rast, test.rast)
      comb.tbl = comb.rast |> 
        as.data.frame(xy=T)
      dat.10 = left_join(dat.10, comb.tbl, by=c("x", "y"))  |> 
        mutate(wilt_train = if_else(is.na(wilt_train), 0, wilt_train),
               wilt_test = if_else(is.na(wilt_test), 0, wilt_test) )
      saveRDS(dat.10, file="mod_data/compare/wilt_split_tbl_10.Rds")
      sum(dat.10$wilt_train)
      sum(dat.10$wilt_test)
    }
    
  }

  {##### assemble data #####
    # models for comparison
    mod.pth = "E:/tmp/wilt_models/train_test/"
    
    gam_30_h0_pois = readRDS( paste0(mod.pth, "train_gam_30_h0_pois.Rds") )
    gam_30_h0_logit = readRDS( paste0(mod.pth, "train_gam_30_h0_logit.Rds") )
    gam_30_ha1_pois = readRDS( paste0(mod.pth, "train_gam_30_ha1_pois.Rds") )
    gam_30_ha1_logit = readRDS( paste0(mod.pth, "train_gam_30_ha1_logit.Rds") )
    gam_30_ha2_pois = readRDS( paste0(mod.pth, "train_gam_30_ha2_pois.Rds") )
    gam_30_ha2_logit = readRDS( paste0(mod.pth, "train_gam_30_ha2_logit.Rds") )
    
    gam_10_h0_pois = readRDS( paste0(mod.pth, "train_gam_10_h0_pois.Rds") )
    gam_10_h0_logit = readRDS( paste0(mod.pth, "train_gam_10_h0_logit.Rds"))
    gam_10_ha1_pois = readRDS( paste0(mod.pth, "train_gam_10_ha1_pois.Rds") )
    gam_10_ha1_logit = readRDS( paste0(mod.pth, "train_gam_10_ha1_logit.Rds") )
    gam_10_ha2_pois = readRDS( paste0(mod.pth, "train_gam_10_ha2_pois.Rds") )
    gam_10_ha2_logit = readRDS( paste0(mod.pth, "train_gam_10_ha2_logit.Rds") )
    
    kppm_30_cauch_pq = readRDS( paste0(mod.pth, "train_kppm_30_pq_cauch.Rds") )
    kppm_30_cauch_wc = readRDS( paste0(mod.pth, "train_kppm_30_wclik1_cauch.Rds") )
    
    kppm_30_mat_pq = readRDS( paste0(mod.pth, "train_kppm_30_pq_mat.Rds") )
    kppm_30_mat_wc = readRDS( paste0(mod.pth, "train_kppm_30_wclik1_mat.Rds") )
    
    kppm_30_thom_pq = readRDS( paste0(mod.pth, "train_kppm_30_pq_thom.Rds") )
    kppm_30_thom_wc = readRDS( paste0(mod.pth, "train_kppm_30_wclik1_thom.Rds") )
    
    kppm_ha_thom_pq = readRDS( paste0(mod.pth, "tkppm_ha2_pq_thom.Rds"))
    kppm_ha_thom_wc = readRDS( paste0(mod.pth, "tkppm_ha2_wc_thom.Rds"))
    kppm_ha_cauch_pq = readRDS( paste0(mod.pth, "tkppm_ha2_pq_cauch.Rds"))
    kppm_ha_cauch_wc = readRDS( paste0(mod.pth, "tkppm_ha2_wc_cauch.Rds"))
    
    dat.30 = readRDS("mod_data/compare/wilt_split_tbl.Rds")
    quads.30 = readRDS("mod_data/kppm/quad_30_30.Rds")
    dat.10 = readRDS("mod_data/compare/wilt_split_tbl_10.Rds")
    
    mod.list = list(
      "gam_30_h0_logit" = gam_30_h0_logit,
      "gam_30_h0_pois" = gam_30_h0_pois,
      "gam_30_ha1_logit" = gam_30_ha1_logit,
      "gam_30_ha1_pois" = gam_30_ha1_pois,
      "gam_30_ha2_logit" = gam_30_ha2_logit,
      "gam_30_ha2_pois" = gam_30_ha2_pois,
      "gam_10_h0_logit" = gam_10_h0_logit,
      "gam_10_h0_pois" = gam_10_h0_pois,
      "gam_10_ha1_logit" = gam_10_ha1_logit,
      "gam_10_ha1_pois" = gam_10_ha1_pois,
      "gam_10_ha2_logit" = gam_10_ha2_logit,
      "gam_10_ha2_pois" = gam_10_ha2_pois,
      "kppm_30_cauch_pq" = kppm_30_cauch_pq,
      "kppm_30_cauch_wc" = kppm_30_cauch_wc,
      "kppm_30_thom_pq" = kppm_30_thom_pq,
      "kppm_30_thom_wc" = kppm_30_thom_wc,
      "kppm_30_mat_pq" = kppm_30_mat_pq,
      "kppm_30_mat_wc" = kppm_30_mat_wc,
      "kppm_30_ha2_thom_pq" = kppm_ha_thom_pq,
      "kppm_30_ha2_thom_wc" = kppm_ha_thom_wc,
      "kppm_30_ha2_cauch_pq" = kppm_ha_cauch_pq,
      "kppm_30_ha2_cauch_wc" = kppm_ha_cauch_wc
      #,
      #"kppm_10_cauch" = kppm_10_cauch,
      #"kppm_10_thom" = kppm_10_thom,
      #"kppm_10_mat" = kppm_10_mat
    )
    
    names(mod.list)
    mod.v = c( rep(1,12), rep(2, 10) ) # 1 for GAM, 2 for CPM
    grid.v = c( 
                rep(30, 6), rep(10, 6), 
                rep(30, 10)#, rep(10, 3) 
                ) # 30 for 30-m grid, 10 for 10-m grid
    
    model.in.data = list(
      "GAM data 30" = dat.30,
      "GAM data 10" = dat.10,
      "quads 30" = quads.30
    )
  }
  
  { ###### run comparison #####
    source("functions/mod/compare_models_tbl_.R")
    compare.res = compare_models_tbl_(mod.list, mod.v, grid.v, model.in.data, verbose=T, DBG=T)
    compare.res |> str()
    cat(crayon::bgRed("\n\nDONE\n\n"))
    saveRDS(compare.res, file="plt_data/model_compare.Rds")
    rm(compare.res)
  }
  
  compare.res = readRDS("plt_data/model_compare.Rds")
  compare.res |> 
    pivot_wider(id_cols=1:2, names_from=measure, values_from=val) |> 
    ggplot(aes(x=area_surveyed, y=recall, color=mod_name)) + geom_line(linewidth=1.2)
  
  compare.res$mod_name |> unique()
  
  compare.res |> 
    filter(mod_name %in% c("gam_30_ha2_logit", #"gam_30_h0_pois", 
                           #"gam_10_ha2_logit", #"gam_10_h0_pois", 
                           #"gam_30_ha2_logit", #"gam_30_ha2_pois",
                           #"gam_10_ha2_logit"#, "gam_10_ha2_pois"
                           "kppm_30_ha2_cauch_wc", "kppm_30_ha2_thom_wc",
                           "kppm_30_thom_wc", "kppm_30_mat_pq"#,
                           #"kppm_30_cauch_pq", "kppm_30_cauch_wc", "kppm_30_thom_pq", 
                           #"kppm_30_mat_wc" 
                           ) ) |> 
    pivot_wider(id_cols=1:2, names_from=measure, values_from=val) |>
    ggplot(aes(x=area_surveyed, y=recall, color=mod_name)) + geom_line(linewidth=1.2) + scale_color_brewer(type="div", palette="RdBu") #+
    #xlim(0, 0.2)
  
}


{
  
  
}


#---------- end ---------------

