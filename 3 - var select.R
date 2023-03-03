{ ###### INIT #####
  rm(list=ls())
  gc()
  library(spatstat)
  library(terra)
  library(tidyverse)
  library(assertthat)
  library(ggpubr)
} |> suppressPackageStartupMessages()

{ ####### LOAD DATA, SUBSET ########
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
    saveRDS(dat.inter.30, file="mod_data/full_dat/dat_interact_30.Rds")
    rm(dat.intr.30)
    gc()
  }
  
  { ####### PLOT 10-m covar vs. disease rate relationships ########
    dat.intr.10 <- readRDS("mod_data/full_dat/dat_interact_10.Rds")
    source("functions/mod/pear_corr_.R")
    corr.plot.10 <- dat.intr.10 |> 
      select(-x, -y, -ow_rast_10) |> 
      pear_corr_()
    plot(corr.plot.10)
    vars.filt.10 <- readRDS("mod_data/full_dat/view_vars_10.Rds")
    
    compare.list <- list(
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
    for(i in seq_along(compare.list)) {
      vars.extr <- compare.list[[i]]
      cat(paste0("\ni = ", i, "\tvar.tmp = ", paste0(vars.extr, collapse=" , ")))
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
          jpeg(filename=paste0("data_img/var_compare/10/", i, "/", j, ".jpg"))
          plot(fig)
          dev.off()
      }
    }
    
    { ####### check cond_ind after removal, vis all ######
      names(dat.intr.10)
      var.subs <- dat.intr.10 |>  #**important*
        select(-x, -y,
               -hillshade, # 1
               -prof_curv, -`conv_ind x prof_curv`, -`plan_curv x prof_curv`, -`ph x prof_curv`, # 2
               -`ph x rsp`, -`rsp`, -channel_dist_s5, # 3
               -topo_wet,# 4
               -`ph x sand`, -`bd x sand`, -sand, -`sand x conv_ind`, -`sand x plan_curv`,# 5
               -`sand x plan_curv`,# 6
               -`bd x plan_curv`,# 7
               -`ph x conv_ind`, -`ph x plan_curv`,# 8
               -`ph x wl2_oakprob_10`,# 9
               -gw_10, -`ph x hillshade`, -`bd x prof_curv`, -`sand x prof_curv`, -slope# extra
                )
      var.subs |> pear_corr_() |> plot()
      source("functions/mod/condition_index_.R")
      ci.val <- var.subs |> 
        select(-ow_rast_10) |> 
        condition_index_(); ci.val
        
      for(i in seq_len(ncol(var.subs))) {
        if(i == 1) next
        var.tmp <- colnames(var.subs)[i]
        cat(paste0("\ni = ", i, "\tvar.tmp = ", var.tmp))
        plot.dat <- var.subs |> 
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
        plot(fig)
        prpt = readline(prompt="\nany key for next...")
      }
      save.dat <- var.subs |> 
        add_column(y = dat.intr.10[["y"]], .after="ow_rast_10") |> 
        add_column(x = dat.intr.10[["x"]], .after="ow_rast_10")
      saveRDS(save.dat, file="mod_data/model_in/dat_10.Rds")
  }
  
    {
      rm(list=ls())
      gc()
      dat.10 <- readRDS("mod_data/model_in/dat_10.Rds")
      library(mgcv)
      names(dat.10)
      names(dat.10)[c(5:7, 11:12)] <- names(dat.10)[c(5:7, 11:12)] |> str_replace_all(pattern=" ", replacement="_")
      names(dat.10)
      mod.tst <- bam(
        ow_rast_10 ~ s(x, y) + aspect + bd_x_conv_ind + bd_x_hillshade + bd_x_ph + 
                        channel_dist_s7 + clay + elev + hillshade_x_prof_curv + ph_x_topo_wet + 
                        plan_curv + val_depth + wl2_oakprob_10,
                    data = dat.10,
                    family=poisson(link="log"), method="fREML", nthreads=6
                     )
      summary(mod.tst)
    }
}
  
  { ####### PLOT 30-m covar vs. disease relationships ########
    dat.intr.30 <- readRDS("mod_data/full_dat/dat_interact_30.Rds")
    vars.filt.30 <- readRDS("mod_data/full_dat/view_vars_30.Rds")
  }
  
  {
    require(corrplot)
    require(RColorBrewer)
    source("functions/mod/pear_corr_.R")
  }
  # 10
  pear_corr_(dat.subs.10 |> 
               select(-ow_rast_10, -x, -y, 
                      -gw_10, # bc elev
                      -clay,  # bc
                      -`sand x plan_curv`, #
                      -`sand x conv_ind`, #
                      -`bd x prof_curv`, #
                      -`sand x prof_curv`, # 
                      -`ph x sand`,
                      -`hillshade x prof_curv`,
                      -`rsp`,
                      -`ph x plan_curv`,
                      -`ph x prof_curv`,
                      -`ph x conv_ind`,
                      -`bd x conv_ind`,
                      -`plan_curv x prof_curv`)) #
  dat.keep.10 <- dat.subs.10 |> 
    select(-gw_10, # bc elev
           -clay,  # bc
           -`sand x plan_curv`, #
           -`sand x conv_ind`, #
           -`bd x prof_curv`, #
           -`sand x prof_curv`, # 
           -`ph x sand`,
           -`hillshade x prof_curv`,
           -`rsp`,
           -`ph x plan_curv`,
           -`ph x prof_curv`,
           -`ph x conv_ind`,
           -`bd x conv_ind`,
           -`plan_curv x prof_curv`)
  saveRDS(dat.keep.10, file="model_data/var12_10_moddat.Rds")
  
  # 30
  pear_corr_(dat.subs.30 |> 
               select(-ow_rast_30, -x, -y,
                      -`ph x hillshade`, 
                      -`bd x sand`,
                      -`ph x prof_curv`,
                      -`ph x plan_curv`,
                      -`ph x conv_ind`,
                      -`bd x conv_ind`,
                      -rsp)
  )
  dat.keep.30 <- dat.subs.30 |> 
    select(-`ph x hillshade`, 
           -`bd x sand`,
           -`ph x prof_curv`,
           -`ph x plan_curv`,
           -`ph x conv_ind`,
           -`bd x conv_ind`,
           -rsp)
  saveRDS(dat.keep.30, file="model_data/var8_30_moddat.Rds")
}

{ ###### condition index ######
  
  source("functions/mod/condition_index_.R")
  
  # moddat.10 <- readRDS("model_data/var12_10_moddat.Rds") |> 
  #   select(-ow_rast_10, -x, -y) # only predictors, no x-y
  # moddat.30 <- readRDS("model_data/var8_30_moddat.Rds") |>
  #   select(-ow_rast_30, -x, -y) # only predictors, no x-y
  # cor.10 <- cor(moddat.10)
  # cor.30 <- cor(moddat.30)
  # eig.10 <- cor.10 |> eigen()
  # eig.30 <- cor.30 |> eigen()
  # 
  # cond.ind.10 <- sqrt(max(eig.10$values)/min(eig.10$values))
  # cond.ind.30 <- sqrt(max(eig.30$values)/min(eig.30$values))
  # 
  # cat(paste0("\ncondition index for 10-m moddat: ", cond.ind.10))
  # cat(paste0("\ncondition index for 30-m moddat: ", cond.ind.30))
}


{ ###### 30-m GAM #########
  rm(list=ls())
  gc()
  library(mgcv)
  dat.mod.30 <- readRDS("model_data/var8_30_moddat.Rds")
  dat.mod.30 |> str()
  names(dat.mod.30)[8:11] <- c("bd_x_prof_curv", "bd_x_plan_curv", "bd_x_hillshade", "ph_x_sand")
  for(i in seq_len(ncol(dat.mod.30))) {
    if(i >= 2) {
      dat.mod.30[,i] <- scale(dat.mod.30[,i])
    }
  }
  dat.mod.30 |> str()
  {
    t1 <- Sys.time()
    tst.gam <- mgcv::bam(ow_rast_30 ~ 
                           s(x, y) + aspect + s(val_dep, bs="cs") + s(wl2_oakprob_30, bs="cs") + s(bd_x_prof_curv, bs="cs") + s(bd_x_plan_curv, bs="cs") + 
                           s(bd_x_hillshade, bs="cs") + s(ph_x_sand, bs="cs"),
                       data=dat.mod.30, family=poisson(link="log"), method="fREML", nthreads=6)
    t2 <- Sys.time()
    cat(paste0("\ntime = ", difftime(t2, t1, units="mins"), " min"))
    cat(paste0("\nAIC = ", AIC(tst.gam)))
    #saveRDS(tst.gam, file="model_data/models/gam_30_freml.Rds")
    tst.gam |> summary()
  }
  tst.gam |> plot()
}

{ ###### 10-m GAM #########
  rm(list=ls())
  gc()
  library(mgcv)
  dat.mod.10 <- readRDS("model_data/var12_10_moddat.Rds")
  dat.mod.10 |> str()
  names(dat.mod.10)[7:14] <- c("ph_x_topo_wet", "ph_x_rsp", "ph_x_wl2_oakprob_10", "bd_x_sand",
                               "bd_x_ph", "bd_x_hillshade", "bd_x_plan_curv", "conv_ind_x_prof_curv")
  for(i in seq_len(ncol(dat.mod.10))) {
    if(i >= 2) {
      dat.mod.10[,i] <- scale(dat.mod.10[,i])
    }
  }
  dat.mod.10 |> str()
  {
    t1 <- Sys.time()
    tst.gam <- mgcv::bam(ow_rast_10 ~ 
                           s(x, y) + aspect + s(elev, bs="cs") + s(ph_x_wl2_oakprob_10, bs="cs") + s(bd_x_hillshade, bs="cs") +
                           s(bd_x_hillshade, bs="cs") + s(bd_x_plan_curv, bs="cs"),
                           #s(bd_x_prof_curv, bs="cs") + s(bd_x_plan_curv, bs="cs") + 
                           #s(bd_x_hillshade, bs="cs") + ,
                         data=dat.mod.10, family=poisson(link="log"), method="fREML", nthreads=6)
    t2 <- Sys.time()
    cat(paste0("\ntime = ", difftime(t2, t1, units="mins"), " min"))
    cat(paste0("\nAIC = ", AIC(tst.gam)))
    saveRDS(tst.gam, file="model_data/models/gam_10_freml.Rds")
    tst.gam |> summary()
  }
  tst.gam |> plot()
}



{ ############ visualize model results
  tst.mod <- readRDS("model_data/models/gam_30_freml.Rds") # load model
  tst.mod$fitted.values # fitted
  # load input data
  # load ow pts
  # plot
  
}