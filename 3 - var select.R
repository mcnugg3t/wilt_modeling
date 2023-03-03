{ ###### INIT #####
  rm(list=ls())
  gc()
  library(spatstat)
  library(terra)
  library(tidyverse)
  library(assertthat)
} |> suppressPackageStartupMessages()

{ ####### LOAD DATA, SUBSET ########
  
  source("functions/add_interact_.R")
  interact.v.10 <- readRDS("clean_data/interact_v.Rds")
  interact.v.30 <- readRDS("clean_data/interact_v_30.Rds")
  # 10
  dat.intr.10 <- terra::rast("clean_data/joindat_10_800_50.tif") |> 
    as.data.frame(xy=T) |> 
    add_interact_(
      interact.v = interact.v.10, 
      verbose=T, DBG=T) 
  dat.intr.10 |> str()
  
  dat.subs.10 <- dat.intr.10 |> 
   mutate(wl2_oakprob_10 = if_else(wl2_cls_10 == 4230, true=wl2_oakprob_10, false=0)) |> 
   select(
        c("ow_rast_10", "x", "y", "aspect", "rsp", "wl2_oakprob_10", "clay", "gw_10", "elev",
          "ph x sand", "ph x conv_ind",  "ph x topo_wet", "ph x plan_curv", "ph x prof_curv", "ph x rsp", "ph x wl2_oakprob_10",
          "sand x prof_curv", "sand x plan_curv",  "sand x conv_ind", "sand x conv_ind",
          "bd x sand", "bd x ph", "bd x hillshade", "bd x prof_curv", "bd x conv_ind", "bd x plan_curv",
          "conv_ind x prof_curv",
          "plan_curv x prof_curv", 
          "hillshade x prof_curv")
    )
  rm(dat.intr.10)
  gc()
  # 30
  dat.intr.30 <- terra::rast("clean_data/joindat_30_800_50.tif") |> 
    as.data.frame(xy=T) |> 
    add_interact_(
      interact.v = interact.v.30, 
      verbose=T, DBG=T) 
  dat.intr.30 |> str()
  dat.subs.30 <- dat.intr.30 |> 
    mutate(wl2_oakprob_30 = if_else(wl2_cls_30 == 4230, true=wl2_oakprob_30, false=0)) |> 
    select("ow_rast_30", "x", "y", "aspect", "clay", "wl2_oakprob_30", "val_dep", "rsp",
           "bd x prof_curv", "bd x plan_curv", "bd x hillshade", "bd x sand", "bd x conv_ind", 
           "ph x sand", "ph x prof_curv", "ph x hillshade", "ph x plan_curv", "ph x conv_ind")
  
  rm(dat.intr.30)
  gc()
  
  dat.subs.10 <- dat.subs.10 |> drop_na()
  dat.subs.30 <- dat.subs.30 |> drop_na()
  
  {
    require(corrplot)
    require(RColorBrewer)
    pear_corr_ <- function(dat.in) {
      M <- cor(dat.in)
      p1 <- corrplot(M, type="upper", order="hclust",
                     col=brewer.pal(n=8, name="RdYlBu"))
    }
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
  moddat.10 <- readRDS("model_data/var12_10_moddat.Rds") |> 
    select(-ow_rast_10, -x, -y) # only predictors, no x-y
  moddat.30 <- readRDS("model_data/var8_30_moddat.Rds") |>
    select(-ow_rast_30, -x, -y) # only predictors, no x-y
  cor.10 <- cor(moddat.10)
  cor.30 <- cor(moddat.30)
  eig.10 <- cor.10 |> eigen()
  eig.30 <- cor.30 |> eigen()
  
  cond.ind.10 <- sqrt(max(eig.10$values)/min(eig.10$values))
  cond.ind.30 <- sqrt(max(eig.30$values)/min(eig.30$values))
  
  cat(paste0("\ncondition index for 10-m moddat: ", cond.ind.10))
  cat(paste0("\ncondition index for 30-m moddat: ", cond.ind.30))
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