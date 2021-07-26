#####################################
## Interpolate E-OBS maximum summer temp (98th quantile) to Irish hectad scale
## 
## inputs:  * eobs.RData - data prepared by Jon Yearsley and shared in the 
##            Grassland dropbox folder, summer 2017
##          * ireland_coastline shapefile (I got this from Jon Yearsley but I
##                think it might be a naturalEarth file originally)
##          * IE_10km_hecs.shp shapefile of Irish hectads created in QGIS by wg
## outputs: krg_summer_tx_predict - a spatialGridDataFrame holding 
##          interpolated values of 98th quantile of maximum summer temperature 
##          over the years 1995-2016 at the hectad scale.  
##          This includes many ocean blocks where 
##          predictions are likely poor.  This should be masked before being
##          used for analysis.  This is the object to use for most things.
##          
##          krg_summer_tx_rast - a rasterized version of krg_summer_tx_predict  
##          
##          krig_summer_tx_map - a rasterBrick version of krg_summer_tx_predict,
##          useful for easy plotting. krg_summer_tx_predict is probably the 
##          object to use for analysis.
## 
##          year_summer_tx - a list containing an element for each year 
##          1995-2016 and each element is a list of length 3 
##          holding the krig predictions as a spatial grid, as a raster, and
##          as a rasterBrick (for easy plotting)
##          
## NOTES: This script first finds the 98th quantile of maximum daily summer
## temperatures for each year, then takes the mean over all years to give 
## mean summer max temp.  It then interpolates that to the hectad scale.
## 
## author: Willson Gaul
## created: 13 Dec 2017
## last modified: 3 March 2020 - change how template raster is produced
#############################################

load("./data/eobs.RData")
print_plots <- F

## ----------------------- load coastline for masking -------------------------
# Load Ireland coastline
ir <- readOGR(dsn='./data/', layer='ireland_coastline')
ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))

### ----------------- prepare hectad raster -----------------------------------
# make 10km square template raster
irish_hec_raster <- raster(xmn = -60000, xmx = 450000, ymn = -70000, ymx = 550000,
                           crs = CRS("+init=epsg:29903"), vals = 1)
res(irish_hec_raster) <- 10000

irish_1km_raster <- raster(xmn = 10000, xmx = 380000, ymn = -30000, ymx = 500000, 
                           crs = CRS("+init=epsg:29903"), vals = 1)
res(irish_1km_raster) <- 1000
### --------------------- end hectad raster -----------------------------------

## calculate 98th quantile of maximum daily temp for each year -----------------
# tx is maximum daily temperature (deg C). 
eobs_data$location <- paste(eobs_data$latitude, eobs_data$longitude, sep = "_")
eobs_data$year <- lubridate::year(eobs_data$date)
eobs_data$month <- lubridate::month(eobs_data$date)
tx_summary <-  select(eobs_data, longitude, latitude, 
                             location, month, year, tx) %>%
  filter(month < 10 & month > 3) %>%
  group_by(location, year) %>%
  mutate(annual_tx = stats::quantile(tx, probs = 0.98)) %>%
  select(-tx, -month) %>%
  distinct()

# calculate mean max summer temp over all years
tx_summary <- ungroup(tx_summary) %>%
  group_by(location) %>%
  mutate(mean_tx = mean(annual_tx, na.rm = T))

mean_tx <- select(tx_summary, -c(year, annual_tx)) %>%
  distinct()

## end calculating annual min winter temp -------------------------------------

## make both yearly and average yearly precip into spatial dfs -----------------
## yearly 
# make yearly summaries into spatial data frames for interpolation
spat_tx <- tx_summary[, 3:ncol(tx_summary)]
coordinates(spat_tx) <- as.matrix(tx_summary[, 1:2])
proj4string(spat_tx) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north -----
spat_tx <- spTransform(spat_tx, CRS("+init=epsg:29903"))
dimnames(spat_tx@coords)[[2]][which(
  dimnames(spat_tx@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_tx@coords)[[2]][which(
  dimnames(spat_tx@coords)[[2]] == "latitude")] <- "northings"

## mean annual
# make multi-year temp summaries into spatial data frames for interpolation
spat_mean_tx <- mean_tx[, 3:ncol(mean_tx)]
coordinates(spat_mean_tx) <- as.matrix(mean_tx[, 1:2])
proj4string(spat_mean_tx) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_mean_tx <- spTransform(spat_mean_tx, 
                            CRS("+init=epsg:29903"))
dimnames(spat_mean_tx@coords)[[2]][which(
  dimnames(spat_mean_tx@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_mean_tx@coords)[[2]][which(
  dimnames(spat_mean_tx@coords)[[2]] == "latitude")] <- "northings"

## ----------------------- interpolate by krigging ----------------------------
# Interpolate data in spat_tx and spat_mean_tx to the raster 
# grid irish_hec_raster
irish_spat_grid <- as(irish_hec_raster, 'SpatialGrid')
irish_spat_grid_1km <- as(irish_1km_raster, 'SpatialGrid')

## mean min (98th quantile) winter temp over all years
# create empirical variogram for mean min winter temp
gs_mean_tx <- gstat(id = "mean_tx", 
                    formula = mean_tx ~ 1, 
                    data = spat_mean_tx) # make gstat object
v_mean_tx <- variogram(gs_mean_tx)
#vgm(as.character(vgm()[, 1])) # show vgm options
f_var_mean_tx <- fit.variogram(v_mean_tx, 
                               vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", 
                                     "Ste"))) 
f_var_mean_tx
if(print_plots) {
  plot(variogramLine(f_var_mean_tx, max(v_mean_tx$dist)), type = 'l')
  points(v_mean_tx[, 2:3], pch = 20, col = 'red')
}

# use variogram in kriging interpolation
krg_mean_tx <- gstat(formula = mean_tx~1, 
                     data = spat_mean_tx, model = f_var_mean_tx)

krg_mean_tx_predict <- predict(krg_mean_tx, irish_spat_grid)
names(krg_mean_tx_predict) <- c("mean_tx", "variance")
krg_mean_tx_predict_1km <- predict(krg_mean_tx, irish_spat_grid_1km)
names(krg_mean_tx_predict_1km) <- c("mean_tx", "variance")

# make rasters
krg_mean_tx_rast <- raster::raster(krg_mean_tx_predict)
krg_mean_tx_rast_1km <- raster::raster(krg_mean_tx_predict_1km)

# masking with shapefile, but I think the shapefile cuts off some hectads that
# do have data (some sp. datasets have 1014 hectads)
if(print_plots) {
  plot(raster::mask(krg_mean_tx_rast, ir_TM75))
}
# krig_mean_tx_map <- raster::brick(krg_mean_tx_predict)
# krig_mean_tx_map <- raster::mask(krig_mean_tx_map, ir_TM75) 
# names(krig_mean_tx_map) <- c("mean_tx", "variance")
# if(print_plots) {
#   plot(krig_mean_tx_map, main = "Average 98th quantile of summer tx")
# }
### save outputs ---------------------------------------------------------------
saveRDS(krg_mean_tx_rast, file = "summer_tx_hectad.rds")
saveRDS(krg_mean_tx_rast_1km, file = "summer_tx_1km.rds")
