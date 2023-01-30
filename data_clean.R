
######## INIT ##########
rm(list=ls())
gc()
{ # load packages
library(tidyverse)
library(assertthat)
library(sf)
library(terra)
} |> suppressPackageStartupMessages()
#
######## WiscLand2 Data & Grid Templates ##########
rm(list=ls())
gc()
## WiscLand2 level 3 classification
#
# crop WL2 data to the study area extent - slightly different for 10 m and 30 m grids
# Crop Wiscland2 data to study area bbox, save
{
  { 
    wl2.class <- terra::rast("raw_data/wiscland2_lvl_3/wiscland2_level3.tif")
    crs(wl2.class)
    sa.bbox <- terra::vect("raw_data/full_area/elev_pull_bbox.shp") |>  # Read bbox shapefile defining full study area 
      project(crs(wl2.class))
    assert_that(crs(sa.bbox) == crs(wl2.class))
    wl2.crop <- terra::crop(wl2.class, sa.bbox, snap="near") 
    terra::writeRaster(wl2.crop, filename="mid_data/30/wiscland2/wl2_cls_30.tif", overwrite=T)
  }
  # create 30 m grid template & save
  { 
    grd.template <- wl2.crop
    values(grd.template) <- 0
    terra::writeRaster(grd.template, filename="mid_data/30/grid/grd_template.tif", overwrite=T)
  }
  # resample cls prob onto 30 m grid
  {
    wl2.cls.prob <- terra::rast("raw_data/wiscland2_lvl_3/wiscland2_level3_4230_conf.tif") |>
      project(crs(grd.template)) |> 
      resample(grd.template, method="bilinear")
    assert_that(crs(wl2.cls.prob) == crs(grd.template))
    terra::writeRaster(wl2.cls.prob, filename="mid_data/30/wiscland2/wl2_oakprob_30.tif", overwrite=T)
  }
  # create 10 m grid template & save
  { 
    grd.template.10 <- rast("raw_data/Wilt_DEM_10_final.tif")
    values(grd.template.10) <- 0
    terra::writeRaster(grd.template.10, filename="mid_data/10/grid/grd_template_10.tif", overwrite=T)
  }
  # resample wl2 and class prob onto 30 m grid
  {
    wl2.class.10 <- wl2.class |> 
      project(crs(grd.template.10)) |> 
      resample(grd.template.10, method="near")
    assert_that(crs(wl2.class.10) == crs(grd.template.10))
    writeRaster(wl2.class.10, filename="mid_data/10/wiscland2/wl2_cls_10.tif")
    wl2.classprob.10 <- wl2.cls.prob |> 
      project(crs(grd.template.10)) |> 
      resample(grd.template.10, method="near")
    assert_that(crs(wl2.class.10) == crs(grd.template.10))
    writeRaster(wl2.classprob.10, filename="mid_data/10/wiscland2/wl2_oakprob_10.tif")
  }
}
#
######## Oak Wilt Points ##########
# clean, merge, and save
{
  rm(list=ls())
  gc()
  
  { 
  ow.pts.2021 <- sf::st_read("raw_data/wilt_points/2021/2021_oakwilt_review_sites.shp") |> 
      select(2:6, 49) |> 
      filter(Active == "yes") |> 
      add_column(year = 2021) |> 
      select(-Active)
  ow.pts.2004.2020 <- sf::st_read("raw_data/wilt_points/2004-2021/All_Oak_Wilt_waypoints_2004-present.shp") 
  pts.clean <- ow.pts.2004.2020 |> 
    select(2:6, 32, 35) |> 
    mutate(dtime = lubridate::as_datetime(ltime),
           year = lubridate::isoyear(dtime)) |> 
    select(1:5, 10)
  assert_that( st_crs(ow.pts.2021) == st_crs(pts.clean) )
  pts.comb <- rbind(ow.pts.2021, pts.clean) |> 
    st_zm(drop=T, what="ZM")
  st_write(pts.comb, dsn="mid_data/wilt/ow_pts_comb.shp", overwrite=T)
  }
  
  rm(list=ls())
  gc()
  
  { # rasterize counts onto grid templates
  grd.temp.10 <- terra::rast("mid_data/10/grid/grd_template_10.tif")
  ow.pts.comb <- terra::vect("mid_data/wilt/ow_pts_comb.shp") |> 
    project(crs(grd.temp.10))
  ow.rast.10 <- terra::rasterizeGeom(ow.pts.comb, grd.temp.10, fun="max")
  #plot(ow.rast.10)
  #max(values(ow.rast.10), na.rm=T)
  terra::writeRaster(ow.rast.10, filename="mid_data/10/wilt/ow_rast_10.tif", overwrite=T)
    
  grd.temp.30 <- terra::rast("mid_data/30/grid/grd_template.tif")
  ow.pts.comb <- terra::vect("mid_data/wilt/ow_pts_comb.shp") |> 
    project(crs(grd.temp.30))
  ow.rast.30 <- terra::rasterizeGeom(ow.pts.comb, grd.temp.30, fun="max")
  max(values(ow.rast.30), na.rm=T)
  terra::writeRaster(ow.rast.30, filename="mid_data/30/wilt/ow_rast_30.tif", overwrite=T)  
  }
  
}
#
######## POLARIS DATA ##########

rm(list=ls())
gc()
source("functions/download_polaris_.R")
download_polaris_()

source("functions/polaris_weighted_avg_.R")
property.v <- c("alpha", "bd", "clay", "hb", "ksat", "lambda", "n", "om", "ph", "sand", "silt", "theta_r", "theta_s")

for(grd in c(10, 30)) {
  if(grd == 10) {
    grd.tmp <- rast("mid_data/10/grid/grd_template_10.tif")
  } else {
    grd.tmp <- rast("mid_data/30/grid/grd_template.tif")
  }
  for(p in seq_along(property.v)) {
    p.tmp <- property.v[p]
    cat(paste0("\np = ", p, "\tprop = ", p.tmp))
    rast.tmp <- polaris_weighted_avg_(str.var = p.tmp, grd.temp = grd.tmp, verbose=T, DBG=T)
    terra::writeRaster(rast.tmp, paste0("mid_data/", toString(grd),"/polaris/", p.tmp, ".tif"), overwrite=T)
  }
}
  
#
######## CNNF Manage Areas ##########

{
  rm(list=ls())
  gc()
  
  { # 30
    cnnf.manage <- terra::vect("raw_data/manage_areas/CNNF_Management_Areas.shp")
    grd.template <- terra::rast("mid_data/30/grid/grd_template.tif")
    # create and save tibble of ma code descriptions
    macode.pre.tbl <- tibble(
      MA = unique(cnnf.manage$MA),
      MA_descrip = unique(cnnf.manage$MA_DESCRIP)
      )
    manage.rast.30 <- cnnf.manage |> 
      project(crs(grd.template)) |> 
      terra::rasterize(grd.template, field="MA")
    terra::writeRaster(manage.rast.30, filename="mid_data/30/manage/manage_rast_30.tif", overwrite=T)
  }
  
  
  { # 10
    grd.template.10 <- terra::rast("mid_data/10/grid/grd_template_10.tif")
    manage.rast.10 <- cnnf.manage |> 
      project(crs(grd.template.10)) |> 
      rasterize(grd.template.10, field="MA")
    terra::writeRaster(manage.rast.10, filename="mid_data/10/manage/manage_rast_10.tif", overwrite=T)
  }
  
  {# construct ma code tbl
  v.tmp <- values(manage.rast.30)[,1]
  c.tmp <- crds(manage.rast.30)
  ma.tbl <- tibble(
    x = c.tmp[,1],
    y = c.tmp[,2],
    val = v.tmp[!is.na(v.tmp)]
  ) |> arrange(desc(y), x  )
  u.vals <- unique(ma.tbl$val)
  ma.code.correct <- tibble(
      ind = u.vals, 
      MA = levels(manage.rast.30)[[1]][,"MA"][u.vals+1]
    ) |> 
    left_join(macode.pre.tbl)
  saveRDS(ma.code.correct, file="mid_data/manage/macode_tbl.Rds")
  }
}
#
######## GW Depth ##########

{
  rm(list=ls())
  gc()
  { # 30
  grd.temp.30 <- rast("mid_data/30/grid/grd_template.tif")
  gw.raw <- terra::rast("raw_data/groundwater/GW.tif")
  ext.to.crop <- grd.temp.30 |> 
    project(crs(gw.raw)) |> 
    ext()
  gw.crop <- gw.raw |> crop(ext.to.crop)
  # resample to grid
  gw.30 <- gw.crop |> 
    project(crs(grd.temp.30)) |> 
    resample(grd.temp.30, method="bilinear")
  writeRaster(gw.30, filename="mid_data/30/groundwater/gw_30.tif")
  }
  
  { # 10
    grd.temp.10 <- rast("mid_data/10/grid/grd_template_10.tif")
    gw.10 <- gw.crop |> 
      project(crs(grd.temp.10)) |> 
      resample(grd.temp.10, method="bilinear")
    writeRaster(gw.10, filename="mid_data/10/groundwater/gw_10.tif")
  }
}

#
######## SOIL LAYER ##########
{
  rm(list=ls())
  gc()
  { # 30
    cnnf.soil <- vect("raw_data/cnnf_soil/CNNF_Soil_Layer_2.shp")
    grd.temp.30 <- rast("mid_data/30/grid/grd_template.tif")
    ext.to.crop <- grd.temp.30 |> 
      project(crs(cnnf.soil)) |> 
      ext()
    cnnf.crop <- cnnf.soil |> 
      crop(ext.to.crop) |> 
      project(crs(grd.temp.30))
    cnnf.rast <- rasterize(cnnf.crop, grd.temp.30, field="OBJECTID")
    writeRaster(cnnf.rast, filename="mid_data/30/soils/soils_rast.tif", overwrite=T)
  }
  rm(grd.temp.30, ext.to.crop, cnnf.crop, cnnf.rast)
  gc()
  { # 10
    grd.temp.10 <- rast("mid_data/10/grid/grd_template_10.tif")
    ext.to.crop <- grd.temp.10 |> 
      project(crs(cnnf.soil)) |> 
      ext()
    cnnf.crop <- cnnf.soil |> 
      crop(ext.to.crop) |> 
      project(crs(grd.temp.10))
    cnnf.rast <- rasterize(cnnf.crop, grd.temp.10, field="OBJECTID")
    writeRaster(cnnf.rast, filename="mid_data/10/soils/soils_rast.tif", overwrite=T)
  }
  
}
#
######## CREATE EXPLORE DATA ##########

{
  rm(list=ls())
  gc()
  source("functions/define_study_area_.R")
  {
    define_study_area_(grd.int = 10, buffer.dist = 1200, verbose=T, DBG=T)
    define_study_area_(grd.int = 30, buffer.dist = 1200, verbose=T, DBG=T)
  }
  
  {
    source("functions/join_data_.R")
    joindat.10 <- join_data_(grd.int=10)
    writeRaster(joindat.10, filename="clean_data/joindat_10_1200.tif", overwrite=T)
    cat(crayon::bgRed("\n\nDONE JOINING 10m x 10m data"))
    rm(joindat.10)
    gc()
    joindat.30 <- join_data_(grd.int=30)
    writeRaster(joindat.30, filename="clean_data/joindat_30_1200.tif", overwrite=T)
    cat(crayon::bgRed("\n\nDONE JOINING 30m x 30m data"))
  }
  
  rm(list=ls())
  gc()
  
  {
    joindat.10 <- rast("clean_data/joindat_10_1200.tif")
    for(i in seq_along(names(joindat.10)) ) {
        cat(paste0("\nplotting : ", names(joindat.10)[i]))
        plot(joindat.10[[i]])
        cat("\npress enter for next...")
        in1 <- readline("")
    }
    rm(joindat.10)
    gc()
  }
  
  {
    joindat.30 <- rast("clean_data/joindat_30_1200.tif")
    for(i in seq_along(names(joindat.30)) ) {
      cat(paste0("\nplotting : ", names(joindat.30)[i]))
      plot(joindat.30[[i]])
      cat("\npress enter for next...")
      in1 <- readline("")
    }
    rm(joindat.30)
    gc()
  }
  
}
# define study area

  





tst.rast <- rast("mid_data/wiscland2/wl2_cls_prob.tif")
crs(tst.rast) == crs(sa.base)

tst <- tst.rast |> resample(sa.base) |> mask(sa.base)

{
  i = 4
  j = 1
  verbose = T
  DBG = T
  base.grd = sa.base
  fold.skip= c("grid", "saga_30", "manage", "study_area")
  var.skip = c("gw_no_smooth", "gw_smooth_cubic", "ow_rast_30", "wilt_buff_rast")
}
