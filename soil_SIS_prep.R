####################################
## Prepare Ireland soil SIS data 
## 
## This script takes about 12 hours to run on WG's work laptop.
## 
## author: Willson Gaul wgaul@hotmail.com 
## created: 25 Feb 2020
## last modified: 23 July 2021
####################################

sis <- readOGR(dsn = "./data/SOIL_SISNationalSoils_shp/Data/SOIL_SISNationalSoils_Shp/", 
               layer = "SOIL_SISNationalSoils")
sis <- spTransform(sis, CRS("+init=epsg:29903"))

soil <- readOGR(dsn = "./data/Soils/Data/", 
                layer = "soils_ie_f1")
soil <- spTransform(soil, CRS("+init=epsg:29903"))

# make 10km, 1km, and 100m square template rasters
irish_hec_raster <- raster(xmn = -60000, xmx = 450000, ymn = -70000, ymx = 550000, 
                           crs = CRS("+init=epsg:29903"), vals = 1)
res(irish_hec_raster) <- 10000
irish_hec_shp <- st_as_sf(rasterToPolygons(irish_hec_raster), 
                          crs = CRS(irish_hec_raster))

irish_1km_raster <- raster(xmn = 10000, xmx = 380000, ymn = -30000, ymx = 500000, 
                           crs = CRS("+init=epsg:29903"), vals = 1)
res(irish_1km_raster) <- 1000
df_1km <- data.frame(rasterToPoints(irish_1km_raster, spatial = F))
colnames(df_1km) <- c("eastings", "northings")


### Soil Drainage calculation --------------------------------------------------
### Calculate percent of grid squares covered by drainage code
# get unique soil drainage codes
sis$DRAINAGE <- as.character(sis$DRAINAGE)
sl_dr <- unique(sis$DRAINAGE)
names(sl_dr) <- sl_dr
sl_dr <- as.list(sl_dr)
sl_dr <- sl_dr[!is.na(sl_dr)]

## calculate percent cover of Drainage classes for 10km squares
for(i in 1:length(sl_dr)) { # loop through soil associations
  lyr <- names(sl_dr)[i]
  sl <- sis[sis$DRAINAGE == lyr & !is.na(sis$DRAINAGE), ]
  # calculate percent of grid cell covered
  sl_r <- rasterize(sl, irish_hec_raster, getCover = TRUE) 
  sl_dr[[i]] <- sl_r
}

soil_drainage_10km_brick <- brick(sl_dr)

## calculate percent cover of Drainage classes for 1km squares
for(i in 1:length(sl_dr)) { # loop through soil associations
  lyr <- names(sl_dr)[i]
  sl <- sis[!is.na(sis$DRAINAGE) & sis$DRAINAGE == lyr, ]
  # calculate percent of grid cell covered
  sl_r <- rasterize(sl, irish_1km_raster, getCover = TRUE) 
  sl_dr[[i]] <- sl_r
}

soil_drainage_1km_brick <- brick(sl_dr)

## save raster bricks
writeRaster(soil_drainage_10km_brick, filename = "soil_drainage_10km_brick.grd", 
            overwrite = TRUE, bylayer = TRUE, 
            suffix = names(soil_drainage_10km_brick))
saveRDS(soil_drainage_10km_brick, file = "soil_drainage_10km_brick.rds")
writeRaster(soil_drainage_1km_brick, filename = "soil_drainage_1km_brick.grd", 
            overwrite = TRUE, bylayer = TRUE, 
            suffix = names(soil_drainage_1km_brick))
saveRDS(soil_drainage_1km_brick, file = "soil_drainage_1km_brick.rds")

rm(sl_dr, sl_r, sl, sis)
### end soil drainage ----------------------------------------------------------


### Soil Association calculation ------------------------------------------------
### Calculate percent of grid squares covered by each soil Association
# get unique soil Association codes
# sl_ass <- unique(sis$Associatio)
# names(sl_ass) <- sl_ass
# sl_ass <- as.list(sl_ass)
# 
# ## calculate percent cover of Associations for 10km squares
# for(i in 1:length(sl_ass)) { # loop through soil associations
#   lyr <- names(sl_ass)[i]
#   sl <- sis[sis$Associatio == lyr, ]
#   # calculate percent of grid cell covered
#   sl_r <- rasterize(sl, irish_hec_raster, getCover = TRUE)
#   sl_ass[[i]] <- sl_r
# }
# # brick the rasters that give the percent of cell covered by each Association
# soil_10km_brick <- brick(sl_ass)
# 
# ## calculate percent cover of Associations for 1km squares
# for(i in 1:length(sl_ass)) { # loop through soil associations
#   lyr <- names(sl_ass)[i]
#   sl <- sis[sis$Associatio == lyr, ]
#   # calculate percent of grid cell covered
#   sl_r <- rasterize(sl, irish_1km_raster, getCover = TRUE)
#   sl_ass[[i]] <- sl_r
# }
# # brick the rasters that give the percent of cell covered by each Association
# soil_1km_brick <- brick(sl_ass)
### end soil association -------------------------------------------------------


### IFS Soil Code ------------------------------------------------------------
### Calculate percent of grid squares covered by IFS Soil
# get unique codes
soil$IFS_SOIL <- as.character(soil$IFS_SOIL)
sl_ifs <- unique(soil$IFS_SOIL)
names(sl_ifs) <- sl_ifs
sl_ifs <- as.list(sl_ifs)
sl_ifs <- sl_ifs[!is.na(sl_ifs)]

## calculate percent cover of Drainage classes for 10km squares
for(i in 1:length(sl_ifs)) { # loop through soil associations
  lyr <- names(sl_ifs)[i]
  sl <- soil[soil$IFS_SOIL == lyr & !is.na(soil$IFS_SOIL), ]
  # calculate percent of grid cell covered
  sl_r <- rasterize(sl, irish_hec_raster, getCover = TRUE) 
  sl_ifs[[i]] <- sl_r
}

soil_IFS_10km_brick <- brick(sl_ifs)

## calculate percent cover of Drainage classes for 1km squares
for(i in 1:length(sl_ifs)) { # loop through soil codes
  lyr <- names(sl_ifs)[i]
  sl <- soil[soil$IFS_SOIL == lyr & !is.na(soil$IFS_SOIL), ]
  # calculate percent of grid cell covered
  sl_r <- rasterize(sl, irish_1km_raster, getCover = TRUE) 
  sl_ifs[[i]] <- sl_r
}

soil_IFS_1km_brick <- brick(sl_ifs)

## save raster bricks
writeRaster(soil_IFS_10km_brick, filename = "soil_IFS_10km_brick.grd", 
            overwrite = TRUE, bylayer = TRUE, 
            suffix = names(soil_IFS_10km_brick))
saveRDS(soil_IFS_10km_brick, file = "soil_IFS_10km_brick.rds")
writeRaster(soil_IFS_1km_brick, filename = "soil_IFS_1km_brick.grd", 
            overwrite = TRUE, bylayer = TRUE,
            suffix = names(soil_IFS_1km_brick))
saveRDS(soil_IFS_1km_brick, file = "soil_IFS_1km_brick.rds")

rm(soil_IFS_10km_brick, soil_IFS_1km_brick, sl_ifs, sl_r, soil)
### end IFS soil code ----------------------------------------------------------
