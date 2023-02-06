
{
  rm(list=ls())
  gc()
  setwd("D:/Backed Up/Desktop/wilt_modeling")
}

{
  library(terra)
  library(tidyverse)
  library(assertthat)
} |> suppressPackageStartupMessages()

# join soil information
{
  dat.10 <- terra::rast("clean_data/joindat_10_1200.tif")
  soil.rast <- dat.10[["soils_rast"]]
  soil.rast.df <- soil.rast |> 
    as.data.frame(xy=T) |> 
    rename(soil.id = soils_rast)
  rm(dat.10, soil.rast)
  gc()
  soil.dat <- terra::vect("raw_data/cnnf_soil/CNNF_Soil_Layer_2.shp")
  soil.df <- soil.dat |> 
    as.data.frame() |> 
    select(-LTA, -LT, -LTP, -ROADS_ANAL, -BIOMASS_HA, -SHAPE_AREA, -SHAPE_LEN) |> 
    rename(soil.id = OBJECTID)
  soil.join.dat <- left_join(soil.rast.df, 
                             soil.df, 
                             by="soil.id")
  rm(soil.dat, soil.id, soil.rast.df)
  gc()
}

# analyze soil data
{
  plot.dat <- soil.join.dat |> 
    mutate(
      soil.name = as.factor(SOIL_NAME),
      surface.texture = as.factor(SURFACE_TE),
      drainage.class = as.factor(DRAINAGE_C),
      compaction = as.factor(COMPACTION),
      erosion = as.factor(EROSION_DI),
      hydric = as.factor(HYDRIC_SOI),
      equipment = as.factor(EQUIPMENT_)
      ) |> 
    select(-SOIL_NAME, -SURFACE_TE, -DRAINAGE_C, -EROSION_DI, -HYDRIC_SOI, -EQUIPMENT_)
  
  soil.join.dat |> 
    select()
  
  smry.dat <- plot.dat |> 
    group_by(soil.name) |> 
    summarise(rel_area = n()/nrow(plot.dat)*100) |> 
    arrange(desc(rel_area))
  
  smry.dat |> 
    ggplot() +
    geom_col(aes(x=soil.name, y=rel_area, fill=soil.name)) +
    theme(axis.text.x = element_text(angle=90))
}

