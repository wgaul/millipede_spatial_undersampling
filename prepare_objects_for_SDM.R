#############################
## Prepare observation data for millipede SDMs
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.  It
## relies on the outputs from 'prepare_data.R'.  This puts the millipede data 
## into a format for use on sonit by fit_SDM.R. This script also makes the
## spatial block cross-validation folds using spatial packages that are not
## available on sonic.  This script saves the required objects as .rds files
## which can be loaded on sonic for fitting SDMs
## 
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 24 Jan 2020
## last modified: 27 July 2020
##############################
library(blockCV)
library(sf)

# spread mill data frame to have a column for each species
# mill_fewer_vars <- select(mill, -DataSetHash, -SurveyHash, -RecorderHash, 
#                           -CommonName, -TaxonGroupName, -StartMonth, -EndMonth, 
#                           -Genus, -Species, -TaxonName, 
#                           -(Name_matched_TNRS:Accepted_name_lsid_TNRS), 
#                           -ID, -cells) 
mill_fewer_vars <- select(mill, occurrenceID, recordedBy, eventDate, 
                          year, month, day, decimalLatitude, decimalLongitude,
                          coordinateUncertaintyInMeters, species, checklist_ID,
                          list_length, day_of_year, sin_doy, cos_doy, 
                          mean_tn:doy_csc)
# rename species column to match name used in later code
colnames(mill_fewer_vars)[colnames(mill_fewer_vars) == "species"] <- 
  "Genus_species"

mill_wide <- mill_fewer_vars %>% 
  mutate(observed = 1) %>%  
  spread(Genus_species, value = observed, fill = 0) %>%
  select(-occurrenceID) %>%
  group_by(checklist_ID) %>%
  summarise_at(vars("Adenomeris gibbosa":"Thalassisobates littoralis"), sum) %>%
  left_join(., mill_fewer_vars[, colnames(mill_fewer_vars) %nin% 
                                 c("occurrenceID", "recordedBy", 
                                   "Genus_species")], 
            by = "checklist_ID") %>%
  unique()

## these were the columns removed in the left_join function when using NBDC
## data.  Keeping them here for reference. 27 July 2021.
# c("occurrenceID", "StartDate", 
#   "EndDate",
#   "SiteName", "GridReference", 
#   "Latitude", "Longitude", 
#   "Genus_species")

mill_wide <- SpatialPointsDataFrame(
  coords = mill_wide[, c("eastings", "northings")], 
  data = mill_wide, proj4string = CRS("+init=epsg:29903"))

if(analysis_resolution == 1000) {
  # get only checklists with spatial resolution of 1km or less
  mill_wide <- mill_wide[mill_wide$Precision <= 1000, ]
}

if(make_spatial_blocks) {
  ### Make spatial block CV folds -------------------------------------------
  # Find best block size using only a subset of all data (for memory efficiency)
  # spatialAutoRange(pred_brick, sampleNumber = 1000, 
  #                  nCores = 3, showPlots = TRUE, plotVariograms = FALSE)
  
  ## make spatial blocks for CV testing
  # Make multiple different CV splits
  fold_assignments <- list()
  for(i in 1:n_cv_trials) {
    for(j in 1:length(cv_block_sizes)) {
      if(cv_block_sizes[j] != "random") {
        fold_assignments[[length(fold_assignments) + 1]] <- spatialBlock(
          mill_spat, 
          theRange = as.numeric(as.character(cv_block_sizes[j])), 
          k = n_folds, 
          selection = "random", 
          iteration = 5, 
          showBlocks = TRUE, 
          xOffset = runif(n = 1, min = 0, max = 1), 
          yOffset = runif(n = 1, min = 0, max = 1),
          rasterLayer = pred_brick$pasture_l2, 
          biomod2Format = FALSE)
      }
      if(cv_block_sizes[j] == "random") {
        fold_assignments[[length(fold_assignments) + 1]] <- list(
          blocks = sample(1:n_folds, size = nrow(mill_wide), replace = T),
          range = "random")
      }
    }
  }
  
  
  # make many smaller blocks for spatial undersampling of test (and training) 
  # data get many different block configurations by making a df, where each 
  # column is a different spatial block configuration and the rows are the 
  # block assignemnts for each millipede record.
  
  checklists_spat <- mill_wide[ , "checklist_ID"] # make df of checklists
  block_subsamp <- checklists_spat
  
  for(i in 1:n_subsamp_block_draws) {
    b_subsamp <- spatialBlock(checklists_spat, 
                              theRange = block_range_spat_undersamp,
                              k = n_folds, 
                              selection = "random", 
                              iteration = 5, 
                              # rows = 20, cols = 20, 
                              xOffset = runif(n = 1, min = 0, max = 1), 
                              yOffset = runif(n = 1, min = 0, max = 1),
                              showBlocks = FALSE, 
                              rasterLayer = pred_brick$pasture_l2, 
                              biomod2Format = FALSE)
    # add spatial subsampling grid cell ID
    block_subsamp <- st_join(
      st_as_sf(block_subsamp), 
      st_as_sf(b_subsamp$blocks[, names(
        b_subsamp$blocks) == "layer"]))
  }
  block_subsamp <- data.frame(block_subsamp)
  # remove geometry column
  block_subsamp <- block_subsamp[, -grepl(".*geomet.*", 
                                                  colnames(block_subsamp))]
} else {
  block_subsamp <- readRDS("block_subsamp.rds")
  fold_assignments <- readRDS("fold_assignments.rds")
}
### end make spatial blocks ---------------------------------------------------

# make new data with standardized recording effort
if(analysis_resolution == 10000) {
  newdata_temp <- cbind(hec_names,  
                        data.frame(raster::extract(pred_brick, hec_names_spat, 
                                                   df = TRUE, 
                                                   method = "simple", 
                                                   cellnumbers = TRUE)))
} else
  if(analysis_resolution == 1000) {
    # Load Ireland coastline
    ir <- readOGR(dsn='./data/', layer='ireland_coastline')
    ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))
    
    newdata_temp <- as.data.frame(mask(pred_brick, ir_TM75))
    newdata_temp <- cbind(newdata_temp, 
                          data.frame(coordinates(mask(pred_brick, ir_TM75))))
    newdata_temp <- newdata_temp[complete.cases(newdata_temp), ]
    colnames(newdata_temp)[colnames(newdata_temp) == "x"] <- "eastings"
    colnames(newdata_temp)[colnames(newdata_temp) == "y"] <- "northings"
    rm(ir, ir_TM75)
  } else stop("Analysis resolution must be either 10000 or 1000 in order to created 'newdata' for predictions.")

newdata <- data.frame()
for(i in 1:18) { # get predictions for every 20 days throught year
  dt <- newdata_temp
  dt$day_of_year <- i*20
  newdata <- bind_rows(newdata, dt)
}
newdata$sin_doy <- sin((2*pi*newdata$day_of_year) / 365)
newdata$cos_doy <- cos((2*pi*newdata$day_of_year) / 365)
newdata$list_length <- 2 # median list length: median(mill_wide_df$list_length)

newdata <- SpatialPointsDataFrame(
  coords = newdata[, c("eastings", "northings")], 
  data = newdata, 
  proj4string = CRS("+init=epsg:29903"))

rm(newdata_temp)
### save objects for use on sonic ---------------------------------------------
saveRDS(mill_wide, "mill_wide.rds")
if(make_spatial_blocks) {
  try(saveRDS(block_subsamp, "block_subsamp.rds"))
  try(saveRDS(fold_assignments, "fold_assignments.rds"))
}
saveRDS(newdata, "newdata.rds")

try(rm(agricultural_l1_rast, agricultural_l1_rast_1km, arable_land_l2_rast, 
       arable_land_l2_rast_1km, artificial_non_ag_vegetated_l2_rast, 
       artificial_non_ag_vegetated_l2_rast_1km, artificial_surfaces_l1_rast, 
       artificial_surfaces_l1_rast_1km, clc_l1_props_1km, clc_l1_props_hecs, 
       elev, elev_1k, forest_l2_rast, forest_l2_rast_1km, 
       forest_seminatural_l1_rast, forest_seminatural_l1_rast_1km, 
       heterogeneous_ag_l2_rast, heterogeneous_ag_l2_rast_1km, 
       industrial_commercial_transport_l2_rast, 
       industrial_commercial_transport_l2_rast_1km, inland_wetlands_l2_rast, 
       inland_wetlands_l2_rast_1km, krg_mean_rr_predict, krg_mean_rr_rast, 
       krg_mean_tn_predict, krg_mean_tn_rast, krg_mean_tx_predict, 
       krg_mean_tx_rast, maritime_wetlands_l2_rast, 
       maritime_wetlands_l2_rast_1km, mean_tn_rast_1k, mean_tx_rast_1k, 
       mill_fewer_vars, mill_wide_df, mill_spat, mine_dump_construction_l2_rast,
       mine_dump_construction_l2_rast_1km, open_space_no_veg_l2_rast, 
       open_space_no_veg_l2_rast_1km, pasture_l2_rast, pasture_l2_rast_1km, 
       permanent_crops_l2_rast, permanent_crops_l2_rast_1km, rr_rast_1k, 
       scrub_herbaceous_l2_rast, scrub_herbaceous_l2_rast_1km, 
       soil_drainage_10km_brick, soil_drainage_1km_brick, soil_IFS_10km_brick, 
       soil_IFS_1km_brick, urban_fabric_l2_rast, urban_fabric_l2_rast_1km, 
       wetlands_l1_rast, wetlands_l1_rast_1km, krig_mean_rr_map, 
       krig_mean_tn_map, krig_mean_tx_map))