

######## INIT ##########

# set path to folder
this.wd <- "D:/Backed Up/Desktop/wilt_modeling/"
# set raw data WD
dat.wd <- paste0(this.wd, "raw_data/")
# load packages
{
library(tidyverse)
library(assertthat)
library(sf)
library(stars)
} |> suppressPackageStartupMessages()

#
######## WiscLand2 Data ##########

##
## WiscLand2 level 3 classification
##

# Set working directory
setwd(dat.wd) 
# Read shapefile defining full study area
full.area <- st_read("full_area/elev_pull_bbox.shp") 
# Check crs
assert_that(st_crs(full.area) == st_crs(3071)) 
# Read WiscLand2 data and reproject to project default CRS (3071)
wl2.class <- stars::read_stars("wiscland2_lvl_3/wiscland2_level3.tif", proxy = F) 
{
st_crs(wl2.class) <- st_crs(3071) # CRS is NAD83 HARN, but strings are non-standard causing errors later
} |> suppressWarnings()
  # Crop Wiscland2 data to the full study area
wl2.crop <- st_crop(wl2.class, full.area) 
# rename attribute
setNames(wl2.crop, c("wl2_3_cls"))
# check CRS
assert_that(st_crs(wl2.crop) == st_crs(3071))
# write cropped data for later use
stars::write_stars(obj = wl2.crop, dsn=paste0(dat.wd, "wiscland2_lvl_3/wl2_cls_crop.tif")) # write cropped classification data
# define grid template using WiscLand2 level 3 grid - background value = 0
grd.template <- wl2.crop %>% 
  mutate(val=0) %>% 
  select(val)
# check CRS
assert_that(st_crs(grd.template) == st_crs(3071))
# write grid template
stars::write_stars(obj=grd.template, dsn=paste0(dat.wd, "study_area/grd_template.tif"))
# remove variables no longer needed
rm(wl2.class); rm(wl2.crop); gc()

##
##
##
## ** BOOKMARK **
##
##
##

##
## WiscLand2 Oak classification probability
##

# Load raw data for Oak classification probability
wl2.clsprob <- stars::read_stars("wiscland2_lvl_3/wiscland2_level3_4230_conf.tif", proxy=F) # read classification probability data -- same issue as comment above - CRS wrongly specified in file
# Same as previously - CRS is NAD83 HARN WTM, but strings are nonstandard and cause issues later
{
st_crs(wl2.clsprob) <- st_crs(3071) # fixed by assignment
} |> suppressWarnings()
# warp because WiscLand2 classification probability is defined on a **VERY SLIGHTLY** different grid than the WL2 data itself
wl2.clsprob.crop <- stars::st_warp(src=wl2.clsprob, dest=grd.template, use_gdal=T, no_data_value=-100, method="bilinear") 
# rename attribute
setNames(wl2.clsprob.crop, c("wl2_cls_prob")) 
# check CRS
assert_that(st_crs(wl2.clsprob.crop) == st_crs(grd.template))
# write
write_stars(obj=wl2.clsprob.crop, dsn=paste0(dat.wd, "wiscland2_lvl_3/wl2_prob_crop.tif")) 
# remove variables no longer needed and collect garbage
rm(wl2.clsprob); rm(wl2.clsprob.crop); gc()

#
######## POLARIS ##########







