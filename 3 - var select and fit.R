{ ###### PKG #####
  rm(list=ls())
  gc()
  library(spatstat)
  library(terra)
  library(tidyverse)
  library(assertthat)
  library(ggpubr)
} |> suppressPackageStartupMessages()

#-------------- setup ---------------------

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

#-------------- GAM ---------------------

{ ###### Fit GAMs on 10-m data ########
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

{ ###### Fit GAMs on 30-m data #########
  
  { ########## setup #########
    {
      rm(list=ls())
      gc()
      library(mgcv)
      library(terra)
      library(tidyverse)
    } |> suppressPackageStartupMessages()
    
    #dat.mod.30 <- readRDS("mod_data/model_in/dat_30.Rds")
    dat.mod.30 = readRDS("mod_data/compare/wilt_split_tbl.Rds")
    names(dat.mod.30)
    
    { ###### replace spaces in variable names with underscores ########
      dat.mod.30 |> str()
      ind.tmp = stringr::str_detect(names(dat.mod.30), pattern=" x ") |> 
        which(); ind.tmp
      names(dat.mod.30)[ind.tmp] <- names(dat.mod.30)[ind.tmp] |> 
        stringr::str_replace_all(pattern=" x ", replacement="_x_")
      str(dat.mod.30)
    }
    { ######## scale vars - not spatial position #########
      for(i in seq_len(ncol(dat.mod.30))) {
        if(i >= 4) {
          dat.mod.30[,i] <- scale(dat.mod.30[,i])
        }
      }
    }
    dat.mod.30 |> names()
    
    dat.mod.30 = dat.mod.30 |> 
      mutate(logit_resp = if_else(ow_rast_30>0, 1, 0)) |> 
      mutate(logit_train = if_else(wilt_train > 0, 1, 0)) |> 
      mutate(wilt_train = if_else(wilt_train < 0, 0, 1))
  }
  
  { ##########  fit test GAM ############
    {
      { ##### POISSON ####
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
      }
      
      { ###### LOGIT #######
        t1 <- Sys.time()
        tst.gam = mgcv::bam(logit_resp ~ 
                              s(x, y) + 
                              s(wl2_oakprob_30, bs="cs") + 
                              s(elev, bs="cs") +
                              s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs") +
                              conv_ind_x_hillshade,
                            # s(elev, aspect, bs=c("ps", "cc") ) + 
                            #te(elev, aspect, bs=c("ps", "cc") ) +
                            #s(topo_wetness, bs="cs") + s(conv_ind_x_hillshade, bs="cs") + #,
                            # bd_x_prof_curv + bd_x_conv_ind + conv_ind_x_hillshade + ph_x_sand,
                            # s(bd_x_prof_curv, bs="cs") + s(bd_x_conv_ind, bs="cs") + 
                            # s(conv_ind_x_hillshade, bs="cs") + s(ph_x_sand, bs="cs"),
                            data=dat.mod.30, family=binomial(link="logit"), method="fREML", nthreads=6)
        t2 <- Sys.time()
        cat(paste0("\ntime = ", difftime(t2, t1, units="mins"), " min"))
        cat(paste0("\nAIC = ", AIC(tst.gam)))
      }
      
      
    }
    #saveRDS(tst.gam, file="model_data/models/gam_30_freml.Rds")
    { ####### diagnostics ########
      tst.gam |> summary()
      tst.gam |> plot()
      library(mgcViz)
      b = getViz(tst.gam)
      plot( sm(b, 1) )
      plot( sm(b, 2) )
      plot( sm(b, 3) )
      plot( sm(b, 4) )
    }
  }
  
  { ######### fit hypothesis GAMs ########### 
    names(dat.mod.30)
    save.pth = "E:/tmp/wilt_models/"

    { ###### H0 - spatial + hosts + env ###### 
      # h0 - poisson
      gam.30.h0.pois = mgcv::bam(ow_rast_30 ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") +
                                   s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs"),
                                 data=dat.mod.30, family=poisson(link="log"), method="fREML", nthreads=6)
      saveRDS(gam.30.h0.pois, file= paste0(save.pth, "gam_30_h0_pois.Rds") )
      
      # h0 - logit
      gam.30.h0.logit = mgcv::bam(logit_resp ~ s(x, y) + 
                                    s(wl2_oakprob_30, bs="cs") + 
                                    s(elev, bs="cs") +
                                    s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs") +
                                    conv_ind_x_hillshade,
                                  data=dat.mod.30, family=binomial(link="logit"), method="fREML", nthreads=6)
      saveRDS(gam.30.h0.logit, file= paste0(save.pth, "gam_30_h0_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(gam.30.h0.pois, gam.30.h0.logit)
    }
    
    { ######### HA1 - spatial + hosts ##########
      # ha1 - poisson
      gam.30.ha1.pois = mgcv::bam(ow_rast_30 ~ s(x, y) + s(wl2_oakprob_30, bs="cs"),
                                 data=dat.mod.30, family=poisson(link="log"), method="fREML", nthreads=6)
      saveRDS(gam.30.ha1.pois, file= paste0(save.pth, "gam_30_ha1_pois.Rds") )
      
      # ha1 - logit
      gam.30.ha1.logit = mgcv::bam(logit_resp ~ s(x, y) + s(wl2_oakprob_30, bs="cs"),
                                  data=dat.mod.30, 
                                  family=binomial(link="logit"), 
                                  method="fREML",nthreads=6)
      AIC(gam.30.ha1.logit)
      saveRDS(gam.30.ha1.logit, file= paste0(save.pth, "gam_30_ha1_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(gam.30.ha1.pois, gam.30.ha1.logit)
    }
    
    { ######### HA2 - hosts + env no spatial ##########
      # ha2 - poisson
      gam.30.ha2.pois = mgcv::bam(ow_rast_30 ~ s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(aspect, bs="cc") + 
                                    s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(ph_x_sand, bs="cs"),
                                  data=dat.mod.30, family=poisson(link="log"), method="fREML", nthreads=6)
      gam.30.ha2.pois |> summary()
      saveRDS(gam.30.ha2.pois, file= paste0(save.pth, "gam_30_ha2_pois.Rds") )
      
      # ha2 - logit
      gam.30.ha2.logit = mgcv::bam(logit_resp ~ s(wl2_oakprob_30, bs="cs") + s(elev, bs="cc") + s(aspect, bs="cc") + 
                                     s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(ph_x_sand, bs="cs") +
                                     s(conv_ind_x_hillshade, bs="cs"),
                                   data=dat.mod.30, family=binomial(link="logit"), method="fREML", nthreads=6)
      saveRDS(gam.30.ha2.logit, file= paste0(save.pth, "gam_30_ha1_logit.Rds") )
      AIC(gam.30.ha2.logit)
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(gam.30.ha2.pois, gam.30.ha2.logit)
    }
  }
  
  { ##### fit PREC-REC GAMs ######
    names(dat.mod.30)
    save.pth = "E:/tmp/wilt_models/train_test/"
    
    { ###### H0 - spatial + hosts + env ###### 
      # h0 - poisson
      train.gam.30.h0.pois = mgcv::bam(wilt_train ~ s(x, y) + s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") +
                                   s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs"),
                                 data=dat.mod.30, family=poisson(link="log"), method="fREML", nthreads=6)
      saveRDS(train.gam.30.h0.pois, file= paste0(save.pth, "train_gam_30_h0_pois.Rds") )
      
      # h0 - logit
      train.gam.30.h0.logit = mgcv::bam(logit_train ~ s(x, y) + 
                                    s(wl2_oakprob_30, bs="cs") + 
                                    s(elev, bs="cs") +
                                    s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs") +
                                    conv_ind_x_hillshade,
                                  data=dat.mod.30, family=binomial(link="logit"), method="fREML", nthreads=6)
      saveRDS(train.gam.30.h0.logit, file= paste0(save.pth, "train_gam_30_h0_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(train.gam.30.h0.pois, train.gam.30.h0.logit)
    }
    
    { ######### HA1 - spatial + hosts ##########
      # ha1 - poisson
      train.gam.30.ha1.pois = mgcv::bam(wilt_train ~ s(x, y) + s(wl2_oakprob_30, bs="cs"),
                                  data=dat.mod.30, family=poisson(link="log"), method="fREML", nthreads=6)
      saveRDS(train.gam.30.ha1.pois, file= paste0(save.pth, "train_gam_30_ha1_pois.Rds") )
      
      # ha1 - logit
      train.gam.30.ha1.logit = mgcv::bam(logit_train ~ s(x, y) + s(wl2_oakprob_30, bs="cs"),
                                   data=dat.mod.30, 
                                   family=binomial(link="logit"), 
                                   method="fREML",nthreads=6)
      saveRDS(train.gam.30.ha1.logit, file= paste0(save.pth, "train_gam_30_ha1_logit.Rds") )
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(train.gam.30.ha1.pois, train.gam.30.ha1.logit)
    }
    
    { ######### HA2 - hosts + env no spatial ##########
      # ha2 - poisson
      train.gam.30.ha2.pois = mgcv::bam(wilt_train ~ s(wl2_oakprob_30, bs="cs") + s(elev, bs="cs") + s(aspect, bs="cc") + 
                                    s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(ph_x_sand, bs="cs"),
                                  data=dat.mod.30, family=poisson(link="log"), method="fREML", nthreads=6)
      #train.gam.30.ha2.pois |> summary()
      cat(crayon::bgRed("\nAIC = ", AIC(train.gam.30.ha2.pois)))
      saveRDS(train.gam.30.ha2.pois, file= paste0(save.pth, "train_gam_30_ha2_pois.Rds") )
      
      # ha2 - logit
      train.gam.30.ha2.logit = mgcv::bam(logit_train ~ s(wl2_oakprob_30, bs="cs") + s(elev, bs="cc") + s(aspect, bs="cc") + 
                                     s(channel_net_s5, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(ph_x_sand, bs="cs") +
                                     s(conv_ind_x_hillshade, bs="cs"),
                                   data=dat.mod.30, family=binomial(link="logit"), method="fREML", nthreads=6)
      saveRDS(train.gam.30.ha2.logit, file= paste0(save.pth, "train_gam_30_ha1_logit.Rds") )
      AIC(train.gam.30.ha2.logit)
      cat(crayon::bgRed("\n\nDONE\n\n"))
      rm(train.gam.30.ha2.pois, train.gam.30.ha2.logit)
    }
  }
}

#-------------- CPM ---------------------

{ ######### Prep CPM data ###########
  
  {
    rm(list=ls())
    gc()
    library(tidyverse)
    library(terra)
    library(spatstat)
    library(assertthat)
  }
  
  { ####### Create quadschemes and save ########
    ow.ppp.10 = readRDS("clean_data/ow_ppp_10.Rds")
    ow.ppp.30 = readRDS("clean_data/ow_ppp_30.Rds")
    train.ppp.10 = readRDS("clean_data/train_ppp_10.Rds")
    train.ppp.30 = readRDS("clean_data/train_ppp_30.Rds")
    
    quads.10 = quadscheme(data=ow.ppp.10, method="grid", eps=30)
    quads.30 = quadscheme(data=ow.ppp.30, method="grid", eps=30)
    quadtrain.10 = quadscheme(data=ow.ppp.10, method="grid", eps=30)
    quadtrain.30 = quadscheme(data=ow.ppp.10, method="grid", eps=30)
    
    saveRDS(quads.10, file="mod_data/kppm/quad_10_30.Rds")
    saveRDS(quads.30, file="mod_data/kppm/quad_30_30.Rds")
    saveRDS(quadtrain.10, file="mod_data/kppm/quad_10_30.Rds")
    saveRDS(quadtrain.30, file="mod_data/kppm/quad_30_30.Rds")
    
    rm(list=ls())
    gc()
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
    rm(dat.mod.30, im.list.30, dat.mod.10, im.list.10)
  }
}

{ #######  Fit CPMs on 30-m data #########
  library(spatstat)
  { ###### Setup #####
    quads.30 = readRDS("mod_data/kppm/quad_30_30.Rds")
    quadtrain.30 = readRDS("mod_data/kppm/quad_30_30.Rds")
    im.list.30 = readRDS("mod_data/kppm/im_list_30.Rds"); names(im.list.30)
    save.pth = "E:/tmp/wilt_models/"
  }
  
  { ####### Hypothesis models ########
    kppm.30.thom = kppm(
      X = quads.30, 
      trend = ~ elev + aspect + bd_x_prof_curv + channel_net_s5,
      clusters="Thomas",
      penalised = T,
      data= im.list.30,
      method = "palm",
      improve.type = "quasi")
    summary(kppm.30.thom)
    saveRDS(kppm.30.thom, file=paste0(save.pth, "kppm_30_thom.Rds"))
    
    kppm.30.mat = kppm(
      X = quads.30, 
      trend = ~ elev + aspect + bd_x_prof_curv + channel_net_s5,
      clusters="MatClust",
      penalised = T,
      data=im.list.30,
      method = "palm",
      improve.type = "quasi")
    summary(kppm.30.mat)
    saveRDS(kppm.30.mat, file=paste0(save.pth, "kppm_30_mat.Rds"))
    
    kppm.30.cauch = kppm(
      X = quads.30, 
      trend = ~ elev + aspect + bd_x_prof_curv + channel_net_s5,
      clusters="Cauchy",
      penalised = T,
      data=im.list.30,
      method = "palm",
      improve.type = "quasi")
    summary(kppm.30.cauch)
    saveRDS(kppm.30.cauch, file=paste0(save.pth, "kppm_30_cauch.Rds"))
    
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
    
  }
  
  { ########## Prec + Rec models #########
    save.pth = "E:/tmp/wilt_models/train_test/"
    
    train.kppm.30.thom = kppm(
      X = quadtrain.30, 
      trend = ~ elev + aspect + bd_x_prof_curv + channel_net_s5,
      clusters="Thomas",
      penalised = T,
      data= im.list.30,
      method = "palm",
      improve.type = "quasi")
    summary(train.kppm.30.thom)
    saveRDS(train.kppm.30.thom, file=paste0(save.pth, "train_kppm_30_thom.Rds"))
    
    train.kppm.30.mat = kppm(
      X = quadtrain.30, 
      trend = ~ elev + aspect + bd_x_prof_curv + channel_net_s5,
      clusters="MatClust",
      penalised = T,
      data=im.list.30,
      method = "palm",
      improve.type = "quasi")
    summary(train.kppm.30.mat)
    saveRDS(train.kppm.30.mat, file=paste0(save.pth, "train_kppm_30_mat.Rds"))
    
    train.kppm.30.cauch = kppm(
      X = quadtrain.30, 
      trend = ~ elev + aspect + bd_x_prof_curv + channel_net_s5,
      clusters="Cauchy",
      penalised = T,
      data=im.list.30,
      method = "palm",
      improve.type = "quasi")
    summary(train.kppm.30.cauch)
    saveRDS(train.kppm.30.cauch, file=paste0(save.pth, "train_kppm_30_cauch.Rds"))
    
    train.kppm.30.cauch.gam = kppm(
      X = quadtrain.30, 
      trend = ~ elev + aspect + bd_x_prof_curv + channel_net_s5,
      clusters="Cauchy",
      data=im.list.30,
      use.gam=T,
      method = "clik2",
      improve.type = "wclik1")
    train.kppm.30.cauch.gam |> predict() |> plot()
    saveRDS(train.kppm.30.cauch.gam, file=paste0(save.pth, "train_kppm_30_cauch_gam.Rds"))
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

#-------------- compare ---------------------
{
  {
    rm(list=ls())
    gc()
    library(terra)
    library(spatstat)
    library(tidyverse)
  } |> suppressPackageStartupMessages()
  
  {##### assemble data #####
    
    # models for comparison
    mod.pth = paste0("E:/tmp/wilt_models/train_test/")
    gam_30_h0_logit = readRDS( paste0(mod.pth, "train_gam_30_h0_logit.Rds") )
    gam_30_h0_pois = readRDS( paste0(mod.pth, "train_gam_30_h0_pois.Rds") )
    
    kppm_30_cauch = readRDS( paste0(mod.pth, "train_kppm_30_cauch.Rds") )
    kppm_30_mat = readRDS( paste0(mod.pth, "train_kppm_30_mat.Rds") )
    kppm_30_thom = readRDS( paste0(mod.pth, "train_kppm_30_thom.Rds") )
    
    
    { #### add train-test partition to dat.30 ####
      # covariate data
      dat.30 = readRDS("mod_data/model_in/dat_30.Rds")
      
      dat.rast = terra::rast("clean_data/joindat_30_800_50.tif")
      ow.pts = terra::vect("mid_data/wilt/ow_pts_comb.shp")
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
      #sum(dat.30$wilt_train)
      #sum(dat.30$wilt_test)
    }
    
    dat.30 = readRDS("mod_data/compare/wilt_split_tbl.Rds")
    quads.30 = readRDS("mod_data/kppm/quad_30_30.Rds")
    
    mod.list = list(
      "gam_30_h0_logit" = gam_30_h0_logit,
      "gam_30_h0_pois" = gam_30_h0_pois,
      "kppm_30_cauch" = kppm_30_cauch,
      "kppm_30_thom" = kppm_30_thom,
      "kppm_30_mat" = kppm_30_mat
    )
    
    names(mod.list)
    mod.v = c(1,1,2,2,2) # 1 for GAM, 2 for CPM
    
    model.in.data = list(
      "GAM data" = dat.30,
      "quads" = quads.30
    )
    
  }
  
  {
    source("functions/mod/compare_models_tbl_.R")
    tst = compare_models_tbl_(mod.list, mod.v, model.in.data, verbose=T, DBG=T)
    tst |> str()
  }
  
  tst |> 
    pivot_wider(id_cols=1:2, names_from=measure, values_from=val) |> 
    ggplot(aes(x=area_surveyed, y=recall, color=mod_name)) + geom_point()
}


{
  
  
}