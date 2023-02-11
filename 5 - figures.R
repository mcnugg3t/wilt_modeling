
{
  rm(list=ls())
  gc()
  library(tidyverse)
  library(assertthat)
  library(plotly)
  library(terra)
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
