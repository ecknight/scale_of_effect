library(tidyverse)
library(raster)
library(rgdal)
library(sp)
library(rgeos)
library(sf)

#EVERYTHING IS IN 10TM

#1. Get extent----
setwd("/Volumes/ECK001/GIS/Projects/Scale")
extent.shp <- read_sf("avie_dep_extent_sa.shp")
extent.r <- raster("wetlpr_30m_sa.tif") %>% 
  reclassify(cbind(-Inf, Inf, 1)) %>% 
  crop(extent.shp) %>% 
  mask(extent.shp)
plot(extent.r)
writeRaster(extent.r, "avie_dep_extent_sa.tif", overwrite=TRUE)

#2. Clip all layers of interest by extent----

#reproject lcc layer
lcc <- raster("lcc_sa.tif")
lcc.10tm <- projectRaster(lcc, crs=crs(extent.r))
writeRaster(lcc.10tm, "lcc_10tm_sa.tif")

interest <- c("avie_sp1.tif",
              "cultivation_sa.tif",
              "dep_excessmoist.tif",
              "dep_moisture.tif",
              "dep_stype1.tif",
              "dep_vegcomp.tif",
              "fire_sa.tif",
              "harvest_sa.tif",
              "industrial_sa.tif",
              "roadlines_sa.tif",
              "seismiclines_sa.tif",
              "wetlpr_30m_sa.tif",
              "lcc_10tm_sa.tif",
              "dep_nutrient.tif",
              "wellsmerge_sa.tif")

for(i in 15:length(interest)){
  r <- raster(interest[i])
  r.resamp <- r %>% 
    resample(extent.r) %>% 
    crop(extent.r) %>% 
    mask(extent.r)
  writeRaster(r.resamp, paste0("1SameExtent/", interest[i]), overwrite=TRUE)
  plot(r.resamp)
}

#3. Reclassify to actual layers I want----
setwd("/Volumes/ECK001/GIS/Projects/Scale/1SameExtent/")

#3a. AVIE----
avi <- raster("avie_sp1.tif")

#3ai. Pine
rclmat <- matrix(c(0, 6.5, 0, 
                   6.5, 7.5, 1,
                   7.5, 9.5, 0,
                   9.5, 11.5, 1,
                   11.5, 13, 0, 
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
pine <- reclassify(avi,rclmat)
plot(pine)
names(pine) <- "pine"
writeRaster(pine, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/pine.tif", overwrite=TRUE)

#3b. LCC----
lcc <- raster("lcc_10tm_sa.tif")

#3bi. Coniferous
rclmat <- matrix(c(0, 2.5, 1, 
                   2.5, 18, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
conifer <- reclassify(lcc, rclmat)
plot(conifer)
names(conifer) <- "conifer"
writeRaster(conifer, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/conifer.tif", overwrite=TRUE)

#3bii. Deciduous
rclmat <- matrix(c(0, 4.5, 0, 
                   4.5, 5.5, 1,
                   5.5, 18, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
decid <- reclassify(lcc, rclmat)
plot(decid)
names(decid) <- "decid"
writeRaster(decid, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/decid.tif", overwrite=TRUE)

#3biii. Mixedwood
rclmat <- matrix(c(0, 5.5, 0, 
                   5.5, 6.5, 1,
                   6.5, 18, 0,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
mixed <- reclassify(lcc, rclmat)
plot(mixed)
names(mixed) <- "mixed"
writeRaster(mixed, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/mixed.tif", overwrite=TRUE)

#3biv. Water
rclmat <- matrix(c(0, 17.5, 0, 
                   17.5, 18, 1,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
water <- reclassify(lcc, rclmat)
plot(water)
names(water) <- "water"
writeRaster(water, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/water.tif")

#3c. HFI----

#3ci. Agriculture
cult <- raster("cultivation_sa.tif")
rclmat <- matrix(c(1, 5, 1, 
                   NA, NA, 0), 
                 ncol=3, byrow=TRUE) 
ag <- cult %>% 
  resample(extent.r) %>% 
  crop(extent.r) %>% 
  mask(extent.r) %>% 
  reclassify(rclmat) %>% 
  mask(extent.r)
plot(ag)
plot(extent.shp, add=TRUE, col=NA)
names(ag) <- "ag"
writeRaster(ag, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/ag.tif")

#3cii. Industrial
industrial <- raster("industrial_sa.tif")
rclmat <- matrix(c(1, 16, 1, 
                   NA, NA, 0), 
                 ncol=3, byrow=TRUE) 
industry <- industrial %>% 
  resample(extent.r) %>% 
  crop(extent.r) %>% 
  mask(extent.r) %>% 
  reclassify(rclmat) %>% 
  mask(extent.r)
plot(industry)
plot(extent.shp, add=TRUE, col=NA)
names(industry) <- "industry"
writeRaster(industry, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/industry.tif")

#3ciii. Seismic
seismiclines <- raster("seismiclines_sa.tif")
rclmat <- matrix(c(1, 2.5, 1,
                   2.5, 3.5, 0,
                   3.5, 4, 1,
                   NA, NA, 0), 
                 ncol=3, byrow=TRUE) 
seismic <- seismiclines %>% 
  resample(extent.r) %>% 
  crop(extent.r) %>% 
  mask(extent.r) %>% 
  reclassify(rclmat) %>% 
  mask(extent.r)
plot(seismic)
plot(extent.shp, add=TRUE, col=NA)
names(seismic) <- "seismic"
writeRaster(seismic, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/seismic.tif")

#3civ. All roads
roadlines <- raster("roadlines_sa.tif")
rclmat <- matrix(c(1, 15.5, 1,
                   15.5, 16.5, 0,
                   16.5, 19.5, 1,
                   19.5, 23, 0,
                   NA, NA, 0), 
                 ncol=3, byrow=TRUE) 
roads <- roadlines %>% 
  resample(extent.r) %>% 
  crop(extent.r) %>% 
  mask(extent.r) %>% 
  reclassify(rclmat) %>% 
  mask(extent.r)
plot(roads)
plot(extent.shp, add=TRUE, col=NA)
names(roads) <- "roads"
writeRaster(roads, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/roads.tif")

#3cv. Gravel roads
roadlines <- raster("roadlines_sa.tif")
rclmat <- matrix(c(1, 6.5, 0,
                   6.5, 7.5, 1,
                   7.5, 10.5, 0,
                   10.5, 11.5, 1,
                   11.5, 23, 0,
                   NA, NA, 0), 
                 ncol=3, byrow=TRUE) 
gravel <- roadlines %>% 
  resample(extent.r) %>% 
  crop(extent.r) %>% 
  mask(extent.r) %>% 
  reclassify(rclmat) %>% 
  mask(extent.r)
plot(gravel)
plot(extent.shp, add=TRUE, col=NA)
names(gravel) <- "gravel"
writeRaster(gravel, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/gravel.tif")

#3d. DEP----

#3di. Moisture
moistureclass <- raster("dep_moisture.tif")
rclmat <- matrix(c(9.5, 10, 1, 
                   8.5, 9.5, 2,
                   7.5, 8.5, 3,
                   5.5, 6.5, 4,
                   1, 1.5, 5,
                   2.5, 3.5, 6,
                   4.5, 5.5, 7,
                   1.5, 2.5, 8,
                   6.5, 7.5, 9,
                   3.5, 4.5, NA,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
moisture <- reclassify(moistureclass, rclmat) %>% 
  calc(fun=function(x){x/9})
plot(moisture)
plot(extent.shp, add=TRUE, col=NA)
names(moisture) <- "moisture"
writeRaster(moisture, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/moisture.tif", overwrite=TRUE)

#3dii. Nutrients
nutrientclass <- raster("dep_nutrient.tif")
rclmat <- matrix(c(1.5, 2.5, 1, 
                   4.5, 5.5, 2,
                   1, 1.5, 3,
                   2.5, 3.5, 4,
                   5.5, 6.5, 5,
                   3.5, 4.5, NA,
                   NA, NA, NA), 
                 ncol=3, byrow=TRUE) 
nutrient <- reclassify(nutrientclass, rclmat) %>% 
  calc(fun=function(x){x/5})
plot(nutrient)
plot(extent.shp, add=TRUE, col=NA)
names(nutrient) <- "nutrient"
writeRaster(nutrient, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/nutrient.tif", overwrite=TRUE)

#3e. Wetland----
wetland <- raster("wetlpr_30m_sa.tif")
names(wetland) <- "wetland"
writeRaster(wetland, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/wetland.tif", overwrite=TRUE)

#3f. Fire----
fireage <- raster("fire_sa.tif")
fire <- fireage %>% 
  calc(fun=function(x){1/(2016-x)})
plot(fire)
plot(extent.shp, add=TRUE, col=NA)
names(fire) <- "fire"
writeRaster(fire, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/fire.tif", overwrite=TRUE)

#3g. Clearcut----
harvestage <- raster("harvest_sa.tif")
harvest <- harvestage %>% 
  calc(fun=function(x){1/(2016-x)})
plot(harvest)
plot(extent.shp, add=TRUE, col=NA)
names(harvest) <- "harvest"
writeRaster(harvest, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/harvest.tif", overwrite=TRUE)

#3h. Wells----

#3cvi. Wells
wellsmerge <- raster("wellsmerge_sa.tif")
wells <- wellsmerge %>% 
  calc(fun=function(x){1/(2016-x)})
plot(wells)
plot(extent.shp, add=TRUE, col=NA)
names(wells) <- "wells"
writeRaster(wells, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/wells.tif", overwrite=TRUE)

#4. Read layers in----
setwd("/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/")
files <- list.files(pattern="*.tif")

layers <- data.frame()
for(i in 1:length(files)){
  name <- str_sub(files[i], -100, -5)
  layer <- raster(files[i])
  names(layer) <- name
  assign(name, layer)
  rm(layer)
  layers <- rbind(layers, data.frame(layer=name))
}

#Check if they stack
stack <- stack(ag, conifer, decid, fire, gravel, harvest, industry, mixed, moisture, nutrient, pine, roads, seismic, water, wells, wetland)
plot(stack) #They do! hooray!

#5. Calculate moving windows----
setwd("/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/")
files <- list.files(pattern="*.tif")
files <- list.files(pattern="wetland.tif")
radii <- c(100, 200, 400, 800, 1600, 3200, 6400, 12800)
loop <- expand.grid(files=files, radius=radii)

for(i in 1:nrow(loop)){
  name.i <- str_sub(loop$files[i], -100, -5)
  layer.i <- raster(as.character(loop$files[i]))
  radius.i <- loop$radius[i]
  if(name.i %in% c("moisture", "nutrient")){
    layer.focal <- focal(layer.i, focalWeight(layer.i, d=radius.i, type='circle'), na.rm=TRUE)
  }
  else{
    layer.focal <- focal(layer.i, focalWeight(layer.i, d=radius.i, type='circle'))
  }
  names(layer.focal) <- paste0(name.i,"-",radius.i)
  writeRaster(layer.focal, paste0("/Volumes/ECK001/GIS/Projects/Scale/3MovingWindow/",name.i,"-",radius.i,".tif"), overwrite=TRUE)
  print(paste0("Completed raster ", name.i, "-", radius.i, " - ", i, " of ", nrow(loop), " rasters"))
}

#6. Calculate extent for 12800 ag-----
setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/LAPRModel/Analysis/2019/TIFs")
files <- "nutrient.tif"
radii <- c(12800)
loop <- expand.grid(files=files, radius=radii)

for(i in 1:nrow(loop)){
  name.i <- str_sub(loop$files[i], -100, -5)
  layer.i <- raster(as.character(loop$files[i]))
  radius.i <- loop$radius[i]
  layer.focal <- focal(layer.i, focalWeight(layer.i, d=radius.i, type='circle'), na.rm=TRUE)
  names(layer.focal) <- paste0(name.i,"-",radius.i)
  writeRaster(layer.focal, paste0(name.i,"-",radius.i,".tif"), overwrite=TRUE)
  print(paste0("Completed raster ", name.i, "-", radius.i, " - ", i, " of ", nrow(loop), " rasters"))
}

#7. Fix fire, harvest, wells, wetland so that NAs are 0s within the extent layer----
extent.r <- raster("/Volumes/ECK001/GIS/Projects/Scale/avie_dep_extent_sa.tif")
setwd("/Volumes/ECK001/GIS/Projects/Scale/1SameExtent/")

#Fire
fireage <- raster("fire_sa.tif")
fire <- fireage %>% 
  calc(fun=function(x){1/(2016-x)}) %>% 
  calc(fun=function(x){ifelse(is.na(x), 0, x)}) %>% 
  crop(extent.r) %>% 
  mask(extent.r)
plot(fire)
plot(extent.shp, add=TRUE, col=NA)
names(fire) <- "fire"
writeRaster(fire, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/fire.tif", overwrite=TRUE)

#Clearcut
harvestage <- raster("harvest_sa.tif")
harvest <- harvestage %>% 
  calc(fun=function(x){1/(2016-x)}) %>% 
  calc(fun=function(x){ifelse(is.na(x), 0, x)}) %>% 
  crop(extent.r) %>% 
  mask(extent.r)
plot(harvest)
plot(extent.shp, add=TRUE, col=NA)
names(harvest) <- "harvest"
writeRaster(harvest, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/harvest.tif", overwrite=TRUE)

#Wells
wellsmerge <- raster("wellsmerge_sa.tif")
wells <- wellsmerge %>% 
  calc(fun=function(x){1/(2016-x)}) %>% 
  calc(fun=function(x){ifelse(is.na(x), 0, x)}) %>% 
  crop(extent.r) %>% 
  mask(extent.r)
plot(wells)
plot(extent.shp, add=TRUE, col=NA)
names(wells) <- "wells"
writeRaster(wells, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/wells.tif", overwrite=TRUE)

#Wetland
setwd("/Volumes/ECK001/GIS/Projects/Scale/1SameExtent")
wetlandpr <- raster("wetlpr_30m_sa.tif")
wetland <- wetlandpr %>% 
  calc(fun=function(x){ifelse(is.na(x), 0, x)}) %>% 
  crop(extent.r) %>% 
  mask(extent.r)
plot(wetland)
plot(extent.shp, add=TRUE, col=NA)
names(wetland) <- "wetland"
writeRaster(wetland, "/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/wetland.tif", overwrite=TRUE)

#8. Rerun moving windows----
setwd("/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/")
files <- list.files(pattern="*.tif")
#radii <- c(100, 200, 400, 800, 1600, 3200)
radii <- c(6400, 12800)
#loop <- expand.grid(files=files, radius=radii) %>% 
#  dplyr::filter(files %in% c("fire.tif", "harvest.tif", "wells.tif", "wetland.tif"))
loop <- expand.grid(files=files, radius=radii) %>% 
  dplyr::filter(files %in% c("wetland.tif"))

for(i in 1:nrow(loop)){
  name.i <- str_sub(loop$files[i], -100, -5)
  layer.i <- raster(as.character(loop$files[i]))
  radius.i <- loop$radius[i]
  layer.focal <- focal(layer.i, focalWeight(layer.i, d=radius.i, type='circle'))
  names(layer.focal) <- paste0(name.i,"-",radius.i)
  writeRaster(layer.focal, paste0("/Volumes/ECK001/GIS/Projects/Scale/3MovingWindow/",name.i,"-",radius.i,".tif"), overwrite=TRUE)
  print(paste0("Completed raster ", name.i, "-", radius.i, " - ", i, " of ", nrow(loop), " rasters"))
}

#9. Clip Arc-processed layers by ag extent----
setwd("/Volumes/ECK001/GIS/Projects/Scale/3MovingWindow")
files <- data.frame(file = list.files(pattern="*6400_arc.tif")) %>% 
  separate(file, into=c("layer", "extent", "extension"), remove=FALSE, sep="_") %>% 
  mutate(file = as.character(file)) %>% 
  dplyr::filter(layer!="ag") 

files <- data.frame(file = list.files(pattern="*12800_arc2.tif")) %>% 
  separate(file, into=c("layer", "extent", "extension"), remove=FALSE, sep="_") %>% 
  mutate(file = as.character(file)) %>% 
  dplyr::filter(layer!="ag") 

ag.6400 <- raster("ag-6400.tif")
ag.12800 <- raster("ag-12800.tif")

for(i in 1:nrow(files)){
  rast.i <- raster(files$file[i])
  extent.i <- files$extent[i]
  if(extent.i=="6400"){
    rast.clip.i <- rast.i %>% 
      crop(ag.6400) %>% 
      mask(ag.6400)
    newname.i <- paste0(files$layer[i], "-", files$extent[i], ".tif")
    writeRaster(rast.clip.i, newname.i, overwrite=TRUE)
  }
  else
  {
    rast.clip.i <- rast.i %>% 
      crop(ag.12800) %>% 
      mask(ag.12800)
    newname.i <- paste0(files$layer[i], "-", files$extent[i], ".tif")
    writeRaster(rast.clip.i, newname.i, overwrite=TRUE)
    
  }
  print(paste0("Completed raster ", files$file[i]))
}

#10. Clip moisture & nutrients by step 2 layers----
setwd("/Volumes/ECK001/GIS/Projects/Scale/3MovingWindow")
files <- data.frame(file = list.files(pattern="*_toclip.tif")) %>% 
  separate(file, into=c("filename", "extension"), remove=FALSE, sep="_") %>% 
  separate(filename, into=c("layer", "extent"), remove=FALSE) %>% 
  mutate(file = as.character(file))

moisture <- raster("/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/moisture.tif")
nutrient <- raster("/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/nutrient.tif")

for(i in 1:nrow(files)){
  rast.i <- raster(files$file[i])
  extent.i <- files$extent[i]
  if(extent.i=="moisture"){
    rast.clip.i <- rast.i %>% 
      crop(moisture) %>% 
      mask(moisture)
    newname.i <- paste0(files$layer[i], "-", files$extent[i], ".tif")
    writeRaster(rast.clip.i, newname.i)
  }
  else
  {
    rast.clip.i <- rast.i %>% 
      crop(nutrient) %>% 
      mask(nutrient)
    newname.i <- paste0(files$layer[i], "-", files$extent[i], ".tif")
    writeRaster(rast.clip.i, newname.i)
    
  }
  print(paste0("Completed raster ", files$file[i]))
}
