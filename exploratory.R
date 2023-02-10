
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
##### 0 - construct spatstat ppp objects and save #### 
{ 
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
##### 1 - Besag's L-function : estimate + gradient #### 
{
  ow.ppp <- readRDS("clean_data/ow_ppp.Rds")
  L.res <- Lest(ow.ppp, correction="translation")
  L.res |> plot()
  Lres.dat <- L.res |> 
    as.data.frame() |> 
    mutate(grad.L = (trans-lag(trans))/(r - lag(r)),
           grad.L.theo = (theo - lag(theo))/(r-lag(r))); Lres.dat
  Lres.dat |> 
    filter(r < 2000) |> 
    select(r, grad.L, grad.L.theo) |> 
    pivot_longer(cols=2:3, names_to = "type", values_to="val") |> 
    ggplot(aes(x=r, y=val, color=type)) + geom_smooth(span=0.10)
}
##
##### 2 - generate densities across parameter grid #### 
{
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
  #
  # generate densities
  { 
    rm(list=ls())
    gc()
    ow.ppp.10 <- readRDS("clean_data/ow_ppp_10.Rds")
    bw.res <- readRDS("clean_data/bw_res_30.Rds")
    source("functions/generate_densities_.R")
    generate_densities_(
      ow.ppp = ow.ppp.10,
      bw.list = bw.res
      )
  }
}

##
## 3 - from each density, create a list of point patterns
{
  #
  # prep
  {
    rm(list=ls())
    gc()
    cln.dat <- terra::rast("clean_data/joindat_10_1200.tif") # load full clean data
    ow.df <- cln.dat[["ow_rast_10"]] |>  # calc ow cell locations and counts
      as.data.frame(xy=T) |> 
      filter(ow_rast_10 > 0)
    rm(cln.dat)
    gc()
    ow.ppp <- readRDS("clean_data/ow_ppp_10.Rds") # load ppp
    wilt.owin <- ow.ppp$window
    file.v <- list.files("clean_data/density/")
  }
  
  #
  # big loop
  {
  source("functions/create_ppps_.R")
      for(j in seq_along(file.v)) { # over each density file
        # prep
        fl.tmp <- file.v[j]
        cat(paste0("\n\n\nj = ", j, "\n\tfl.tmp = ", fl.tmp))
        path.tmp <- paste0("clean_data/density/", fl.tmp)
        bw.num <- fl.tmp |> 
          str_remove("dens_bw_") |> 
          str_remove(".Rds") |> 
          as.numeric()
        save.tmp <- paste0("clean_data/sample/ppp_list_bw_", bw.num, ".Rds" )
        cat(paste0("\n\tsave.tmp = ", save.tmp))
        
        # read density
        cat(paste0("\n\treading density.df, na's -> 0..."))
        density.tmp <- readRDS(path.tmp) |> 
          as.data.frame() |> 
          mutate(val.old = value) |> 
          mutate(value = if_else(is.na(value), 0, value))
        cutoff.tmp <- (0-density.tmp$value)^2 |> quantile(probs=c(0.6))
        density.df <- density.tmp |> 
          mutate(value = if_else( ((0-value)^2) < cutoff.tmp, 0, value), # if distance between 0 and value is < 60th_percentile, set to 0
                 value = value/sum(value, na.rm=T))
        n.dens = density.df |> filter(value > 0) |> nrow()
        cat(paste0("\n\tcutoff = ", cutoff.tmp, " (", n.dens, ")"))
        cat(crayon::bgMagenta("\tcreating ppp list..."))
        ppp.list <- create_ppps_(
          density.df,
          ow.df,
          n.sim=100,
          owin.in = wilt.owin
        )
        cat("\n\t\tsaving...")
        saveRDS(ppp.list, file=save.tmp)
      }
    }
}

# simulation envelope for L under:
# 1 - inhomogeneous poisson
# 2 - cox process
{
  { # setup
    rm(list=ls())
    gc()
    ow.ppp <- readRDS("clean_data/ow_ppp_10.Rds")
    files.tmp <- list.files("clean_data/sample/")
  }
  # loop over point patterns
  for(j in seq_along(files.tmp)) {
    # setup
    fn.tmp <- files.tmp[j]
    fp.tmp <- paste0("clean_data/sample/", fn.tmp)
    cat(paste0("\n\n\tj = ", j, "\tfl.tmp = ", fn.tmp))
    bw.tmp <- fn.tmp |> 
        str_remove("ppp_list_bw_") |> 
        str_remove(".Rds")
    save.tmp <- paste0("clean_data/envelope_", bw.tmp, ".Rds")
    save.img <- paste0("clean_data/envelope/img_", bw.tmp, ".jpg")
    cat(paste0("\n\t\tsave.tmp = ", save.tmp, "\n\t\tsave.img = ", save.img, "\n\t"))
    #
    cat(crayon::bgCyan("loading ppp list..."))
    ppp.list.tmp <- readRDS(fp.tmp)
    cat("\n\t")
    cat(crayon::bgYellow("calc envelope...\n"))
    envelope.tmp <- envelope(
        ow.ppp, 
        fun=Lest, 
        funargs=list(rmax=700,correction="border"),
        simulate= ppp.list.tmp # list of point patterns
    )
    cat("\n\t")
    cat(crayon::bgWhite("saving plot..."))
    jpeg(file=save.img)
    envelope.tmp |> plot(main=paste0("bandwidth: ", bw.tmp))
    dev.off()
    cat("\n\t")
    cat(crayon::bgWhite("saving envelope..."))
    saveRDS(envelope.tmp, file=save.tmp)
    }
}

# Fry plot locally + globally (?)
{
  
}

#
# compute sampling distribution for each kde - both continuous and discrete variables
#
{
  # prep
  {
    rm(list=ls())
    gc()
    soils.df <- readRDS("clean_data/soils_df.Rds")
    dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
  }
  
  #
  {
    source("functions/sample_kdes_.R")
    ##
    ##
    ##
    rm.vars <- c()
    ##
    ##
    ##
    sample_kdes_(
      covar.rast = dat.10,
      join.dat = soils.df,
      n.sim=1e5, 
      verbose=T, 
      DBG=T)
  }
}

##
## Calculate surprisal distributions
{
  
  { # setup
    rm(list=ls())
    gc()
    soils.dat.df <- readRDS("clean_data/soils_df.Rds")

    ow.pts.dat <- terra::rast("clean_data/joindat_10_1200.tif") |> 
      as.data.frame() |>
      filter(ow_rast_10 > 0) |> 
      left_join(soils.dat.df) |> 
      select(-`study area`, -soils_rast, -ow_rast_10) |> 
      select(-rm.vars)
    var.v <- names(ow.pts.dat); var.v
    source("functions/calc_surprisal_hist_.R")
    source("functions/calc_surprisal_tab_.R")
    return.dat <- tibble(
      var = character(),
      kernel = character(),
      bw = numeric(),
      surprisal = numeric()
    )
    fold.v <- list.files("clean_data/sample_dist")
  }
  
  
  
  { ## RUN BIG LOOP
    for(i in seq_along(fold.v)) { # loop over folders
        fold.tmp <- fold.v[i]
        fold.path.tmp <- paste0("clean_data/sample_dist/", fold.tmp, "/")
        cat(paste0("\n\ni = ", i, "\tfolder = ", fold.tmp))
        dens.files <- list.files(fold.path.tmp)
        
        for(j in seq_along(dens.files)) {
          fn.tmp <- dens.files[j]
          bw.tmp <- fn.tmp |> 
            str_remove("samp_dist_bw_") |>
            str_remove(".Rds") |> 
            as.numeric()
          cat(paste0("\n\n\tj = ", j, "\tfn.tmp = ", fn.tmp, "   bw = ", bw.tmp))
          fp.tmp <- paste0(fold.path.tmp, fn.tmp)
          
          # read sampling dist file
          sample.dist.tmp <- readRDS(fp.tmp)
          #
          
          for(k in seq_along(var.v)) { # inner loop
            var.tmp <- var.v[k]
            cat("\n\t\t")
            cat(crayon::bgBlue("k = ", k, "\t var = ", var.tmp))
            sample.v <- ow.pts.dat[[var.tmp]]
            dens.tmp <- sample.dist.tmp[[var.tmp]]
            #
            class.tmp <- class(dens.tmp)
            if(class.tmp == "histogram") {
              surp.v <- calc_surprisal_hist_(dens.tmp, sample.v)
            } else if (class.tmp == "table") {
              surp.v <- calc_surprisal_tab_(dens.tmp, sample.v)
            }
            #
            n.dat <- length(surp.v)
            cat("\n\t\t")
            cat(crayon::bgWhite("adding to data..."))
            return.dat <- return.dat |> 
              add_row(
                var = rep(var.tmp, n.dat),
                kernel = rep(fold.tmp, n.dat),
                bw = rep(bw.tmp, n.dat),
                surprisal = surp.v
              )
        }
      }
    } # end outer loop
  }
  
  {
      var.v <- return.dat$var |> unique(); var.v
      
      # infinite values
      return.dat |> 
        group_by(var, bw) |> 
        summarise(count_inf = sum(is.infinite(surprisal))) |> 
        arrange(desc(count_inf))
      
      # max non-infinite value
      max.non.inf <- return.dat |> 
        filter(!is.infinite(surprisal)) |> 
        select(surprisal) |> 
        min(); max.non.inf
      
      replace.inf <- max.non.inf-4
      
      plot.dat <- return.dat |> 
        mutate(surprisal = if_else(is.infinite(surprisal), replace.inf, surprisal))
      
      
      for(i in seq_along(var.v)) {
        
        { # prep for plot
          var.tmp <- var.v[i]
          cat("\nViewing : ")
          cat(crayon::bgBlue(var.tmp))
          p1.dat <- plot.dat |> 
            filter(var == var.tmp) |> 
            mutate(kernel = as.factor(kernel),
                 bw = as.factor(bw)) |> 
            group_by(kernel, bw)
          
          surp.v <- p1.dat$surprisal
          #mean.surprisal <- mean(-1*surp.v, na.rm=T); mean.surprisal
          #median.surprisal = median(log(-1*surp.v), na.rm=T); median.surprisal
          
          plot.dat <- p1.dat |> 
            ungroup() |> 
            mutate(
              bw.num = as.numeric(levels(bw))[bw],
              plt.x = bw.num + runif( n=nrow(p1.dat), min=-5, max=5 ),
              plt.y = log(-1*surprisal) )
          
          central.dat <- plot.dat |> 
            group_by(bw.num) |> 
            summarise(mean = mean(plt.y),
                      median = median(plt.y)) |>
            mutate(skew = (mean-median)^2) |> 
            ungroup() |> 
            pivot_longer(cols=2:3, names_to="centrality", values_to = "val")
      }
      
      # want to plot the mean and median surprisal for each
      
      
      {
        clrs <- c("mean" = "violet", "median" = "deeppink")
        plot.tmp <- plot.dat |> 
          ggplot(aes(x=log(plt.x, base=10), y=log(plt.y) )) + 
          theme(
            plot.title = element_text(size=24),
            axis.text.x = element_text(size=14, angle=0),
            panel.border = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill="black")) + 
          geom_hex(bins=30) +
          geom_point(data=central.dat, aes(x=log(bw.num, base=10), y=log(val), size=skew, color=centrality)) + 
          labs(title = var.tmp, x = "log(bandwidth)", y="log(-surprisal)" ) +
          scale_color_manual(values = clrs) +
          scale_fill_viridis_c(option="D", name="frequency")
        plot(plot.tmp)
      }
      
      
      cat("\nPress any key for next...")
      prpt <- readline(prompt=" ")
    }
    
  }
  
  
}

