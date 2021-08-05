#############################
## prepare data for millipede maps of ignorance project
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 25 Oct 2019
## last modified: 29 July 2021
##############################

### load millipede data
# NBDC training and test data
mill <- read_tsv("./data/GBIF_27_Feb_2021_download/occurrence.txt")
# remove unwanted columns
mill <- mill[ , -c(2, 5:11, 13:15, 18:26, 28:36, 38:48, 51:54, 67, 72:79, 82:86,
                   95, 104:105, 107:108, 140:166, 212:215)]

hec_names <- read_csv("./data/Irish_land_hectads.csv")
options("scipen"=100, "digits"=4) # make sure R doesn't use scientific notation
hec_names$en <- paste(hec_names$eastings, hec_names$northings, sep = "_")
# put hectad names onto mill data frame (for use calculating spatial evenness 
# of training data)
mill_en <- lat_lon_to_osi(
  mill, lat = "decimalLatitude", lon = "decimalLongitude", 
  orig_crs = "+init=epsg:4326", 
  precision = 10000)[, c(ncol(mill)+1, ncol(mill)+2)]
# add 5000 meters to make coordinates match hec_names (which seems to mark the
# middle of hectads)
mill_en <- mill_en + 5000
mill_en$en <- paste(mill_en$decimalLongitude.1, 
                    mill_en$decimalLatitude.1, 
                     sep = "_")
mill_en <- left_join(mill_en, hec_names, by = "en")
options("scipen"=0, "digits"=7) # restore defaults

# add hectad column to mill
mill$hectad <- mill_en$hectad
rm(mill_en)

# Load Ireland coastline
ir <- readOGR(dsn='./data/', layer='ireland_coastline')
ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))

### prepare and load environmental data ---------------------------------------
if(!all(file.exists("annual_precip_hectad.rds") & 
        file.exists("annual_precip_1km.rds") &
        file.exists("summer_tx_hectad.rds") & 
        file.exists("summer_tx_1km.rds") & 
        file.exists("winter_tn_hectad.rds") &
        file.exists("winter_tn_1km.rds") & 
        file.exists("elevation_hec_ETOPO1.rds") & 
        file.exists("elevation_1km_ETOPO1.rds") & 
        file.exists("corine_label_1_hectad.RData") & 
        file.exists("corine_label_1_1km.RData") & 
        file.exists("corine_label_2_hectad.RData") & 
        file.exists("corine_label_2_1km.RData") &
        file.exists("soil_IFS_10km_brick.rds") & 
        file.exists("soil_IFS_1km_brick.rds") & 
        file.exists("soil_drainage_10km_brick.rds") & 
        file.exists("soil_drainage_1km_brick.rds"))) {
  source("./eobs_annual_precipitation_hectad.R")
  source("./eobs_max_summer_temp_hec.R")
  source("./eobs_min_winter_temp_hec.R")
  source("./interp_elev_hec_etopo1.R")
  source("./prep_corine.R")
  source("./soil_SIS_prep.R")
}
annual_precip_hectad <- read_rds("annual_precip_hectad.rds")
# load("annual_precip_hectad.RData")
rr_rast_1k <- read_rds("annual_precip_1km.rds")
# load("summer_tx_hectad.RData")
summer_tx_hectad <- read_rds("summer_tx_hectad.rds")
mean_tx_rast_1k <- read_rds("summer_tx_1km.rds")
# load("winter_tn_hectad.RData")
winter_tn_hectad <- read_rds("winter_tn_hectad.rds")
mean_tn_rast_1k <- read_rds("winter_tn_1km.rds")
elev <- read_rds("elevation_hec_ETOPO1.rds")
elev_1k <- read_rds("elevation_1km_ETOPO1.rds")
load("corine_label_1_hectad.RData")
load("corine_label_1_1km.RData")
load("corine_label_2_hectad.RData")
load("corine_label_2_1km.RData")
soil_drainage_10km_brick <- read_rds("soil_drainage_10km_brick.rds")
soil_drainage_1km_brick <- read_rds("soil_drainage_1km_brick.rds")
soil_IFS_10km_brick <- read_rds("soil_IFS_10km_brick.rds")
soil_IFS_1km_brick <- read_rds("soil_IFS_1km_brick.rds")

# remove Swamp layer (no data there as far as I can tell)
soil_IFS_10km_brick <- subset(soil_IFS_10km_brick,
                              which(names(soil_IFS_10km_brick) != "Swamp"))
soil_IFS_1km_brick <- subset(soil_IFS_1km_brick,
                              which(names(soil_IFS_1km_brick) != "Swamp"))

# pred_brick_10k <- brick(list(
#   mean_tn = resample(krg_mean_tn_rast, krg_mean_rr_rast), 
#   mean_tx = resample(krg_mean_tx_rast, krg_mean_rr_rast), 
#   mean_rr = krg_mean_rr_rast, 
#   artificial_surfaces = resample(artificial_surfaces_l1_rast, krg_mean_rr_rast),
#   forest_seminatural_l1 = resample(forest_seminatural_l1_rast, 
#                                    krg_mean_rr_rast),
#   wetlands_l1 = resample(wetlands_l1_rast, krg_mean_rr_rast), 
#   pasture_l2 = resample(pasture_l2_rast, krg_mean_rr_rast), 
#   arable_l2 = resample(arable_land_l2_rast, krg_mean_rr_rast), 
#   elev = resample(elev, krg_mean_rr_rast), 
#   soil_drainage = resample(soil_drainage_10km_brick, krg_mean_rr_rast), 
#   soil_IFS = resample(soil_IFS_10km_brick, krg_mean_rr_rast)))
# mask pred brick by one of the CORINE layers to get only Irish land cells
# pred_brick_10k <- mask(pred_brick_10k, pred_brick_10k$artificial_surfaces)
# # scale and centre environmental predictors over study extent
# pred_brick_10k <- scale(pred_brick_10k)

pred_brick_1k <- brick(list(
  mean_tn = mean_tn_rast_1k, 
  mean_tx = mean_tx_rast_1k, 
  mean_rr = rr_rast_1k, 
  artificial_surfaces = artificial_surfaces_l1_rast_1km, 
  forest_seminatural_l1 = forest_seminatural_l1_rast_1km, 
  wetlands_l1 = wetlands_l1_rast_1km, 
  pasture_l2 = pasture_l2_rast_1km,
  arable_l2 = arable_land_l2_rast_1km,
  elev = elev_1k)) #soil_drainage = soil_drainage_1km_brick,soil_IFS = soil_IFS_1km_brick
# pred_brick_1k <- mask(pred_brick_1k, ir_TM75) 
# pred_brick_1k <- scale(pred_brick_1k) # shouldn't matter for RF
rm(ir_TM75, ir)

if(analysis_resolution == 10000) pred_brick <- pred_brick_10k else
  if(analysis_resolution == 1000) pred_brick <- pred_brick_1k else
    stop("Analysis resolution does not match any of the resolutions at which we have prepared predictor varialbes.  Analysis resolution must be either 10000 or 1000.")

# make hec_names spatial 
hec_names_spat <- SpatialPointsDataFrame(
  coords = hec_names[, c("eastings", "northings")], 
  data = hec_names, proj4string = CRS("+init=epsg:29903"))
# make sure millipede data is in same projection as predictor data
hec_names_spat <- spTransform(hec_names_spat, 
                              raster::projection(pred_brick))
### end prepare predictor variables -------------------------------------------

### prepare millipede data ----------------------------------------------------
# remove records with no date
mill <- mill[!is.na(mill$eventDate), ]
# keep only records with coordinate uncertainty <= 1km
mill <- mill[!is.na(mill$coordinateUncertaintyInMeters) & 
               mill$coordinateUncertaintyInMeters <= 1000, ]
# keep only records from 1970 to present
mill <- mill[mill$year >= 1970, ]
mill <- mill[mill$taxonRank == "SPECIES", ]
# About a third of records are at monthly temporal precision, but that will 
# be ok I think.  
warning("Make day_of_year into month_of_year instead?")

## make checlist ID variable
mill$checklist_ID <- paste0(mill$recordedBy, mill$eventDate, mill$locality, 
                            mill$decimalLatitude, mill$decimalLongitude)
## calculate list lengths for each checklist
mill$list_length <- NA
for(i in 1:length(unique(mill$checklist_ID))) {
  id <- unique(mill$checklist_ID)[i]
  ll <- length(which(mill$checklist_ID == id))
  mill$list_length[mill$checklist_ID == id] <- ll
}

## make Julian day and year variables
mill$day_of_year <- yday(mill$eventDate)
mill$sin_doy <- sin((2*pi*mill$day_of_year) / 365)
mill$cos_doy <- cos((2*pi*mill$day_of_year) / 365)
# mill$temp_resolution <- mill$EndDate - mill$StartDate # can't do this w/GBIF data

# join predictor variables to millipede data
# make millipede data spatial 
mill_spat <- SpatialPointsDataFrame(coords = mill[, c("decimalLongitude", 
                                                      "decimalLatitude")], 
                                    data = mill, 
                                    proj4string = CRS("+init=epsg:4326"))
# make sure millipede data is in same projection as predictor data
mill_spat <- spTransform(mill_spat, raster::projection(pred_brick))

pred_df <- data.frame(raster::extract(pred_brick, mill_spat, method = "simple",
                                      df = TRUE, cellnumbers = TRUE))
mill <- cbind(mill, pred_df)
### end prepare millipede data ------------------------------------------------

## scale years, day of year, and spatial coordinate dimensions
# eastings and northings should be scaled in the same way so space isn't more
# important in one direction
east_mean <- mean(mill$decimalLongitude)
north_mean <- mean(mill$decimalLatitude)
mill$eastings_csc <- mill$decimalLongitude - east_mean
mill$northings_csc <- mill$decimalLatitude - north_mean
spat_sd_eastings <- sd(mill$eastings_csc)
spat_sd_northings <- sd(mill$northings_csc)
mill$eastings_csc <- mill$eastings_csc/spat_sd_eastings
mill$northings_csc <- mill$northings_csc/spat_sd_northings

mill$year_csc <- mill$year - mean(mill$year)
mill$year_csc <- mill$year_csc/sd(mill$year_csc)
mill$doy_csc <- mill$day_of_year - 182.5 # center day of year
mill$doy_csc <- mill$doy_csc/sd(mill$doy_csc)

# # I think 10 km in space is about as big as a year and as big as about 8 julian
# # days when all variables are scaled:
# max(mill$year) - min(mill$year)
# (max(mill$northings) - min(mill$northings))/10000
# 365/43

