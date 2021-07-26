###################################
## Interpolate E-OBS total annual precipitation (1995-2016) to Irish hectad scale
## 
## inputs:  * eobs.RData - data prepared by Jon Yearsley and shared in the 
##            Grassland dropbox folder, summer 2017
##          * ireland_coastline shapefile (I got this from Jon Yearsley but I
##                think it might be a naturalEarth file originally)
##          * IE_10km_hecs.shp shapefile of Irish hectads created in QGIS by wg
## outputs: krg_mean_rr_predict - a spatialGridDataFrame holding interpolated
##          values of mean annual precipitation over the years 1995-2016 
##          (excluding 2010-2012) at the hectad scale.  
##          This includes many ocean blocks where 
##          predictions are likely poor.  This should be masked before being
##          used for analysis.  This is the object to use for most things.
##          
##          krg_mean_rr_rast - a rasterized version of krg_mean_rr_predict.  
##          
##          krig_mean_rr_map - a rasterBrick version of krg_mean_rr_predict, 
##          useful for easy plotting.  krg_mean_rr_predict is probably the 
##          object to use for analysis.
## 
##          year_tot_rr - a list containing an element for each year 1995-2016
##          (excluding 2010-2012), and each element is a list of length 3 
##          holding the krig predictions as a spatial grid, as a raster, and
##          as a rasterBrick (for easy plotting)
##          
##              
## NOTES: This script first summarizes daily precipitation to produce annual
## precipitation amounts, then interpolates that to the hectad scale. 
##
## TODO:  - linear interpolation
##        - compare krig and linear results
# alternatviely, try approxfun() or the R equivalent of Matlab "interp" "interp2" for 2d.
## 
## author: Willson Gaul
## created: 25 Sep 2017
## last modified: 3 March 2020 - change how template raster is produced
#################################### 

load("./data/eobs.RData")

print_plots <- F
do_linear <- F # do linear interpolation also (at bottom script)?

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

### calculate annual precipitation --------------------------------------------
# rr is daily precipitation in mm 
# There are some NA values in rr.  Diagnostic code at the end of this script
# looks at these.
# It seems all the NA values are from 2010-2012, and affect many hectads 
# (the majority of hectads in the Republic).  So, as an initial way of dealing
# with this, I propose just not using 2010-2012 in calculating average annual
# precipitation.

## calculate average annual precipitation from 1995-2015 (excluding 2010-2012)
eobs_data$location <- paste(eobs_data$latitude, eobs_data$longitude, sep = "_")
eobs_data$year <- lubridate::year(eobs_data$date)
precip_summary <-  select(eobs_data, longitude, latitude, 
                          location, year, rr) %>%
  group_by(location, year) %>%
  mutate(annual_rr_mm = sum(rr)) %>%
  select(-rr) %>%
  distinct()

# remove sums for the years that have NAs
drop_yrs <- unique(eobs_data$year[which(is.na(eobs_data$rr))])
precip_summary$annual_rr_mm[which(precip_summary$year %in% drop_yrs)] <- NA

## calculate mean annual precipitation
precip_summary <- ungroup(precip_summary) %>%
  group_by(location) %>%
  mutate(mean_annual_rr_mm = mean(annual_rr_mm, na.rm = T))

mean_annual_rr <- select(precip_summary, -c(year, annual_rr_mm)) %>%
  distinct()

## end calculating annual precipitation --------------------------------------

## make both yearly and average yearly precip into spatial dfs -----------------
## ------------------ convert eobs lat/long to OSI east/north ----------------
## yearly
# make precipitation summaries into spatial data frames for interpolation
spat_precip <- precip_summary[, 3:ncol(precip_summary)]
coordinates(spat_precip) <- as.matrix(precip_summary[, 1:2])
# set proj4string based on WGS84 here: http://spatialreference.org/ref/epsg/4326/
proj4string(spat_precip) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_precip <- spTransform(spat_precip, CRS("+init=epsg:29903"))
dimnames(spat_precip@coords)[[2]][which(
  dimnames(spat_precip@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_precip@coords)[[2]][which(
  dimnames(spat_precip@coords)[[2]] == "latitude")] <- "northings"

## mean annual
# make precipitation summaries into spatial data frames for interpolation
spat_mean_rr <- mean_annual_rr[, 3:ncol(mean_annual_rr)]
coordinates(spat_mean_rr) <- as.matrix(mean_annual_rr[, 1:2])
# set proj4string based on WGS84 here: http://spatialreference.org/ref/epsg/4326/
proj4string(spat_mean_rr) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_mean_rr <- spTransform(spat_mean_rr, CRS("+init=epsg:29903"))
dimnames(spat_mean_rr@coords)[[2]][which(
  dimnames(spat_mean_rr@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_mean_rr@coords)[[2]][which(
  dimnames(spat_mean_rr@coords)[[2]] == "latitude")] <- "northings"
## -------------------------- end coordinate conversion ------------------------

## ----------------------- interpolate by krigging ----------------------------
# I want to interpolate data in spat_eobs to the raster grid irish_hec_raster
# following this: http://rspatial.org/analysis/rst/4-interpolation.html
irish_spat_grid <- as(irish_hec_raster, 'SpatialGrid')
irish_1km_grid <- as(irish_1km_raster, 'SpatialGrid')

## mean annual precipitation
# create empirical variogram for mean annual precipitaion
gs_mean_rr <- gstat(id = "mean_yrl_rr", formula = mean_annual_rr_mm~1, 
                 data = spat_mean_rr) # make gstat object
v_mean_rr <- variogram(gs_mean_rr, width = 1000)
#vgm(as.character(vgm()[, 1])) # show vgm options
f_var_mean_rr <- fit.variogram(v_mean_rr, 
                               vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", 
                                "Ste"))) 
f_var_mean_rr
if(print_plots) {
  plot(variogramLine(f_var_mean_rr, max(v_mean_rr$dist)), type = 'l', 
       ylim = c(0, max(v_mean_rr[, 3])))
  points(v_mean_rr[, 2:3], pch = 20, col = 'red')
}

# use variogram in kriging interpolation
krg_mean_rr <- gstat(formula = mean_annual_rr_mm~1, 
                     data = spat_mean_rr, model = f_var_mean_rr)

krg_mean_rr_predict <- predict(krg_mean_rr, irish_spat_grid)
names(krg_mean_rr_predict) <- c("mean_annual_rr", "variance")
krg_mean_rr_predict_1km <- predict(krg_mean_rr, irish_1km_grid)
names(krg_mean_rr_predict_1km) <- c("mean_annual_rr", "variance")

## !!! ---- RESULT -------- !!!
## krg_mean_rr_predict has the results for mean annual precipitation from 
# 1995-2016 (excluding 2010-2012).  
# make rasters
krg_mean_rr_rast <- raster::raster(krg_mean_rr_predict)
krg_mean_rr_rast_1km <- raster::raster(krg_mean_rr_predict_1km)

# masking with shapefile, but I think the shapefile cuts off some hectads that
# do have data (some sp. datasets have 1014 hectads)
if(print_plots) {
  plot(raster::mask(krg_mean_rr_rast, ir_TM75))
  plot(raster::mask(krg_mean_rr_rast_1km, ir_TM75))
}

# krig_mean_rr_map <- raster::brick(krg_mean_rr_predict)
# # krig_mean_rr_map <- raster::mask(krig_mean_rr_map, ir_TM75) 
# names(krig_mean_rr_map) <- c("mean_annual_rr", "variance")
# if(print_plots) {
#   plot(krig_mean_rr_map, main = "Average annual rr")
# }

### Check NAs and errors ------------------------------------------------------
# these plots only work after creating spat_eobs based on 'nas' rather than
# based on eobs data.  Don't need to view every time, diagnostic only
view_precip_nas <- F
if(view_precip_nas == T) {
  nas <- eobs_data[which(is.na(eobs_data$rr)), ]
  nas$location <- paste0(nas$latitude, nas$longitude)
  nas$year <- lubridate::year(nas$date)
  
  for (i in 1:length(unique(nas$location))) {
    print(unique(nas$location)[i])
    print(table(nas$year[nas$location == unique(nas$location)[i]]))
  }
  ## ------------------ convert eobs lat/long to OSI east/north ----------------
  ## yearly
  # make precipitation summaries into spatial data frames for interpolation
  spat_eobs <- nas[, 3:ncol(nas)] 
  coordinates(spat_eobs) <- as.matrix(nas[, 1:2])
  # set proj4string based on WGS84 here: http://spatialreference.org/ref/epsg/4326/
  proj4string(spat_eobs) <- CRS(
    "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  # convert lat/long to OSI east/north
  spat_eobs <- spTransform(spat_eobs, CRS("+init=epsg:29903"))
  dimnames(spat_eobs@coords)[[2]][which(
    dimnames(spat_eobs@coords)[[2]] == "longitude")] <- "eastings"
  dimnames(spat_eobs@coords)[[2]][which(
    dimnames(spat_eobs@coords)[[2]] == "latitude")] <- "northings"
  
  plot(ir_TM75)
  plot(spat_eobs[, ], add = T) # plot hectads with NA values
  # look at individual hectads based on locations in 'nas'
  plot(spat_eobs[which(spat_eobs$location == "51.875-8.125"), ], col = "red", 
       add = T)
  plot(spat_eobs[which(spat_eobs$location == "51.875-8.375"), ], col = "red", 
       add = T)
  plot(spat_eobs[which(spat_eobs$location == "53.875-9.375"), ], col = "red", 
       add = T)
}
## end diagnosing NAs

### save outputs ---------------------------------------------------------------
saveRDS(krg_mean_rr_rast, file = "annual_precip_hectad.rds")
saveRDS(krg_mean_rr_rast_1km, file = "annual_precip_1km.rds")
