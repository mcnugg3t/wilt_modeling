
######## INIT ##########
rm(list=ls())
gc()
#
# SET PATH TO THIS FOLDER IF NOT USING RPROJECT
#
this.wd <- "D:/Backed Up/Desktop/wilt_modeling/"

{ # load packages
library(tidyverse)
library(assertthat)
library(sf)
library(terra)
} |> suppressPackageStartupMessages()

{ # construct intermediate paths
  main.wd <- this.wd
  dat.wd <- paste0(this.wd, "raw_data/") # raw data
  mid.wd <- paste0(this.wd, "mid_data/")
}
#
######## Oak Wilt Points ##########
set.wd(dat.wd)

{ # clean, merge, and save
ow.pts.2021 <- sf::st_read("wilt_points/2021/2021_oakwilt_review_sites.shp") |> 
    select(2:6, 49) |> 
    filter(Active == "yes") |> 
    add_column(year = 2021) |> 
    select(-Active)

ow.pts.2004.2020 <- sf::st_read("wilt_points/2004-2021/All_Oak_Wilt_waypoints_2004-present.shp") 

pts.clean <- ow.pts.2004.2020 |> 
  select(2:6, 32, 35) |> 
  mutate(dtime = lubridate::as_datetime(ltime),
         year = lubridate::isoyear(dtime)) |> 
  select(1:5, 10)

assert_that( st_crs(ow.pts.2021) == st_crs(pts.clean) )

tst <- sf::st_coordinates(pts.clean) |> as_tibble()

pts.comb <- rbind(ow.pts.2021, pts.clean) |> 
  st_zm(drop=T, what="ZM")

st_write(pts.comb, dsn=paste0(mid.wd, "ow_pts_comb.shp"))
}
# rasterize counts onto grid templates

{
grd.temp.10 <- terra::rast("mid_data/grd_template_10.tif")
ow.pts.comb <- terra::vect("mid_data/ow_pts_comb.shp")
ow.rast.10 <- terra::rasterizeGeom(ow.pts.comb, grd.temp.10, fun="max")
#plot(ow.rast.10)
max(values(ow.rast.10), na.rm=T)
terra::writeRaster(ow.rast.10, filename="mid_data/ow_rast_10.tif", overwrite=T)
  
grd.temp.30 <- terra::rast("mid_data/grd_template.tif")
ow.rast.30 <- terra::rasterizeGeom(ow.pts.comb, grd.temp.30, fun="max")
max(values(ow.rast.30), na.rm=T)
terra::writeRaster(ow.rast.30, filename="mid_data/ow_rast_30.tif", overwrite=T)  
}



#
######## WiscLand2 Data ##########

## WiscLand2 level 3 classification
#
{ # read raw WL2.3 data, bounding box, project bbox to WL2.3 crs 
setwd(dat.wd) # move into raw data directory
wl2.class <- terra::rast("wiscland2_lvl_3/wiscland2_level3.tif")
sa.bbox <-terra::vect("full_area/elev_pull_bbox.shp") |>  # Read bbox shapefile defining full study area 
  project(crs(wl2.class))
assert_that(crs(sa.bbox) == crs(wl2.class))
}

{ # Crop Wiscland2 data to study area bbox, save
wl2.crop <- terra::crop(wl2.class, sa.bbox, snap="near") 
terra::writeRaster(wl2.crop, filename=paste0(mid.wd, "wl2_crop.tif"))
}

{ # define grid template using WiscLand2 level 3 grid - background value = 0, save
grd.template <- wl2.crop
values(grd.template) <- 0
terra::writeRaster(grd.template, filename=paste0(mid.wd, "grd_template.tif"))
}

{
wl2.cls.prob <- terra::rast("wiscland2_lvl_3/wiscland2_level3_4230_conf.tif") |> 
    project(crs(grd.template)) |> 
    resample(grd.template, method="near")
assert_that(crs(wl2.cls.prob) == crs(grd.template))
terra::writeRaster(wl2.cls.prob, filename=paste0(mid.wd, "wl2_cls_prob.tif"))
}


#
######## POLARIS DATA ##########
setwd(this.wd)

source("functions/download_polaris_.R")
download_polaris_()

source("functions/weighted_average_.R")
property.v <- c("alpha", "bd", "clay", "hb", "ksat", "lambda", "n", "om", "ph", "sand", "silt", "theta_r", "theta_s")
for(p in seq_along(property.v)) {
  p.tmp <- property.v[p]
  cat(paste0("\np = ", p, "\tprop = ", p.tmp))
  rast.tmp <- weighted_average_(this.wd, p.tmp)
  terra::writeRaster(rast.tmp, paste0("mid_data/polaris_wa_100_", p.tmp, ".tif"))
}

#
######## CNNF Manage Areas ##########
rm(list=ls())
gc()
{
grd.template <- terra::rast("mid_data/grid/grd_template.tif")
cnnf.manage <- terra::vect("raw_data/manage_areas/CNNF_Management_Areas.shp")
# create and save tibble of ma code descriptions
macode.tbl <- tibble(
  ind = 1:24,
  MA = unique(cnnf.manage$MA),
  MA_descrip = unique(cnnf.manage$MA_DESCRIP)
)
saveRDS(macode.tbl, file="mid_data/manage/macode_tbl.Rds")
manage.rast <- terra::rasterize(cnnf.manage, grd.template, field="MA")
plot(manage.rast)
values(manage.rast)[,1]
terra::writeRaster(manage.rast, filename="mid_data/manage/manage_rast.tif")
}

#
######## GW Depth ##########
rm(list=ls())
gc()
{
gw.raw <- terra::rast("raw_data/groundwater/GW.tif")
sa.bbox <- terra::vect("raw_data/full_area/elev_pull_bbox.shp") 
sa.bbox.crop <- sa.bbox |> 
  terra::project(crs(gw.raw))
gw.crop <- terra::crop(gw.raw, sa.bbox.crop) |> 
  terra::project(crs(sa.bbox))
terra::writeRaster(gw.crop, filename="mid_data/gw_no_smooth.tif")
grd.template <- terra::rast("mid_data/grd_template.tif")
gw.smooth.cubic <- terra::resample(gw.crop, grd.template, method="cubic")
gw.smooth.bilin <- terra::resample(gw.crop, grd.template, method="bilinear")
plot(gw.smooth.cubic)
plot(gw.smooth.bilin)
terra::writeRaster(gw.smooth.cubic, filename="mid_data/gw_smooth_cubic.tif")
terra::writeRaster(gw.smooth.bilin, filename="mid_data/gw_smooth_bilin.tif")
}
#
######## SAGA DATA ##########
rm(list=ls())
gc()
setwd(this.wd)
# create and save 10m grid template, 30m DEM for large-scale description
{
dem10 <- terra::rast("raw_data/Wilt_DEM_10_final.tif")
grd.template <- terra::rast("mid_data/grd_template.tif")
grd.template10 <- dem10 |> terra::project(crs(grd.template))
assert_that(crs(grd.template10) == crs(grd.template))
terra::writeRaster(grd.template10, filename="mid_data/grd_template_10.tif")
dem30 <- terra::resample(dem10, grd.template, method="average", threads=T)
terra::writeRaster(dem30, filename="raw_data/Wilt_DEM_30.tif")
}

#
######## CREATE EXPLORE DATA ##########
rm(list=ls())
gc()

# want a dataframe capturing the full shebang - e.g. at 10 m and 30 m
# will then crop the dataframe based on something external determining SA

# 10 m
{}
grd.temp.10 <- terra::rast("mid_data/grid/grd_template_10.tif")

# loop over mid_data - skip grid
# for each .tif file that's not in exclude list

source("functions/join_data_.R")