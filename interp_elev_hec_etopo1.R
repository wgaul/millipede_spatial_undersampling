#################################
## Interpolate ETOPO1 elevation to IE hectads
## 
## inputs:  - the 'ETOPO1_Ice_g_gmt4.grd' dataset that I downloaded 8 May 2019
##          - 'ireland_coastline' shapefile
##          - IE_10km_hecs.shp shapefile
## outputs: - "elevation_hec_ETOPO1.RData" an .RData file with a raster holding 
##            elevation values at the hectad scale
##
## author: W. Gaul
## created: 24 Nov 2017
## last modified: 10 May 2019 changed location of grid 
##      25 Oct 2019 - changed directory to millipede map of ignorance project
##      4 March 2020 - make a 1km raster also
################################

print_plots <- F

## Load Ireland coastline
ir <- readOGR(dsn='./data/', layer='ireland_coastline')
ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))
rm(ir)

## ----------------- prepare hectad raster -----------------------------------
# make 10km square template raster
irish_hec_raster <- raster(xmn = -60000, xmx = 450000, ymn = -70000, ymx = 550000, 
                           crs = CRS("+init=epsg:29903"), vals = 1)
res(irish_hec_raster) <- 10000

# make 1km template raster
irish_1km_raster <- raster(xmn = 10000, xmx = 380000, ymn = -30000, ymx = 500000, 
                           crs = CRS("+init=epsg:29903"), vals = 1)
res(irish_1km_raster) <- 1000

## load elevation data
d  <- raster('./data/ETOPO1_Ice_g_gmt4.grd')
proj4string(d) <- CRS("+init=epsg:4326")
## end load data --------------------------------------------------------------

# crop to IE extent
orig_elev <- crop(d, extent(c(-11, -5, 51, 55.5)))
orig_elev <- projectRaster(orig_elev, crs = CRS("+init=epsg:29903"))

if(print_plots) {
  plot(orig_elev)
  plot(ir_TM75, add = T)
  plot(mask(orig_elev, ir_TM75))
}

# interpolate
# I want to interpolate data in orig_elev to the raster grid irish_hec_raster

# create empirical variogram for elevation
elev_spatDF <- as(orig_elev, "SpatialPointsDataFrame")
names(elev_spatDF@data) <- "elevation"
# gs_elev <- gstat(id = "elevation", formula = elevation~1, 
#                  data = elev_spatDF) # make gstat object
# v <- variogram(gs_elev, width = 1000)

v <- variogram(elevation~1, elev_spatDF)
f_var <- fit.variogram(v, vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", 
                                "Ste")))

f_var
if(print_plots){
  plot(variogramLine(f_var, max(v$dist)), type = 'l', 
       ylim = c(0, max(v[, 3])))
  points(v[, 2:3], pch = 20, col = 'red')
}
# use variogram in kriging interpolation
# use only 24 nearest neighbors to save memory
krg <- gstat(NULL, "elevation", elevation~1, elev_spatDF, model = f_var, 
             nmax = 24) 
elev_hec <- interpolate(irish_hec_raster, krg)
elev_1km <- interpolate(irish_1km_raster, krg)

if(print_plots){
  par(mfrow = c(1, 2))
  plot(mask(orig_elev, ir_TM75))
  plot(mask(elev_hec, ir_TM75))
  par(mfrow = c(1, 1))
}

## save results ---------------------------------------------------------------
# this is the prefered format for single objects
saveRDS(elev_hec, file = "elevation_hec_ETOPO1.rds") 
saveRDS(elev_1km, file = "elevation_1km_ETOPO1.rds")