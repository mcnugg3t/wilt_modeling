
{
  rm(list=ls())
  gc()
  library(tidyverse)
  library(assertthat)
  library(plotly)
  library(terra)
}

{ ########## READ DENSITY #########
  rm(list=ls())
  gc()
  dens.dat <- readRDS("clean_data/density/dens_bw_99.Rds")
  tst1 <- dens.dat$v # 3111 rows, 2878 columns
  tst2 <- dens.dat$v |> as.vector()
  assert_that( 
    all(
      tst2[1:3111][!is.na(tst2[1:3111])] == tst1[,1][!is.na(tst1[,1])]) ) # check that as.vector unfolds column-wise := all non-NA values are same between first 3111 values as.vector and first column of matrix
  # so it unfolds column-wise, with the result that coordinates should repeat the first x value, for each value in v.val 
  dens.tbl <- tibble(
    x = numeric(),
    y = numeric()
  )
  # for each column
  for(i in seq_len(ncol(tst1))) {
    dens.tbl <- dens.tbl |> 
      add_row(
        x = rep(dens.dat$xcol[i], 3111),
        y = dens.dat$yrow
      )
  }
  dens.tbl <- dens.tbl |> 
    add_column(val = tst2)
  
  fig <- plot_ly(
    data = dens.tbl[!is.na(dens.tbl$val),], 
    type="scatter3d",
    mode="markers") |> 
    add_markers(x = ~ x,
                y = ~ y,
                z = ~ val,
                color = ~val)
  fig
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
