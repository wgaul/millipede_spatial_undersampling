################################
## This script organises the analysis for the millipede maps of ignorance
## 
## TODO:
##  - put all environmental data in this project directory
##  - add soil data
##  - add geology data
##  
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 25 Oct 2019
## last modified: 21 Sep 2021
#################################

rm(list = ls())
dbg <- F
calc_1k_distances <- F # run distances for 1km grid (might take a long time)
seed <- 23012020  # 23 Jan 2020


run_rf <- F
make_spatial_blocks <- F # takes a few minutes. Set to T for final run
get_partial_dependence <- F # calculate partial dependence (time consuming)
run_evals <- F

analysis_resolution <- 1000 # analysis resolution (10000 or 1000 m grid squares)
n_folds <- 3 # number of cross-validation folds to use
n_cv_trials <- 33 # number of different cross-validation fold layouts to use
cv_block_sizes <- c("random") # sizes of CV spatial blocks (in meters), or "random" for random (not spatial block) CV
n_subsamp_block_draws <- 3000 # number of spatial subsampling block configurations to make
block_range_spat_undersamp <- 30000 # spatial undersampling grid block size (m)

if(make_spatial_blocks == T) set.seed(seed) # only set on 1st run to creat spatial blocks

library(wgutil)
library(Hmisc)
library(rgdal)
library(gdalUtils)
library(ggplot2)
library(fields)
library(gstat)
library(raster)
library(parallel)
library(sf) # can't install sf on sonic as of 8 Jan 2020
library(fasterize)
library(rgdal)
library(dismo)
#library(rgeos)
library(lubridate)
library(tidyverse)
library(randomForest)

source("functions_maps_of_ignorance.R")

n_cores <- 2

# select species to fit models to
sp_to_fit <- list("Macrosternodesmus palicola", "Boreoiulus tenuis",
                  "Ommatoiulus sabulosus", "Blaniulus guttulatus",
                  "Glomeris marginata", "Cylindroiulus punctatus")
# sp_to_fit <- "Cylindroiulus punctatus"
names(sp_to_fit) <- sp_to_fit

# define environmental predictors for each species
sp_predictors <- list(
  "Macrosternodesmus palicola" = c("mean_rr", "artificial_surfaces", 
                                   "arable_l2"),
  "Boreoiulus tenuis" = c("mean_rr", "elev", "arable_l2", 
                          "artificial_surfaces"),
  "Ommatoiulus sabulosus" = c("mean_rr", "elev", "artificial_surfaces", 
                              "wetlands_l1"), 
  "Blaniulus guttulatus" = c("mean_tn", "mean_rr", "elev", 
                             "artificial_surfaces", "arable_l2", 
                             "forest_seminatural_l1"),
  "Cylindroiulus latestriatus" = c("mean_tn", "mean_rr", "elev", 
                                   "artificial_surfaces", "arable_l2", 
                                   "pasture_l2", "forest_seminatural_l1"),
  "Glomeris marginata" = c("mean_tn", "mean_rr", "elev", "artificial_surfaces", 
                           "arable_l2", "forest_seminatural_l1"), 
  "Cylindroiulus punctatus" = c("mean_tn", "mean_rr", "elev", 
                                "artificial_surfaces", "arable_l2", 
                                "pasture_l2", "forest_seminatural_l1"))

source("prepare_data.R")
source("prepare_objects_for_SDM.R")
# mod_names <- c("month_ll_rf", "spat_ll_rf")
# mod_names <- c("env_ll_rf", "env_spat_ll_rf")
mod_names <- c("month_ll_rf", "spat_ll_rf", "env_ll_rf", "env_spat_ll_rf")
mods_for_pd_plots <- c("env_spat_ll_rf")

if(run_rf) source("fit_rf.R")

## evaluate models
# Specify the model(s) using the string that was used as the beginning of the 
# file names that hold the fitted model objects.  e.g. here use "day_ll_rf" 
# for the file day_ll_rf_noSubSamp_fits_Julus_scaninavius.rds
if(run_evals) {
  for(mod_name in mod_names) {
  evals <- data.frame() # data frame to hold evaluation results
  source("evaluate_models.R")
  }
}

## make plots
source("plots_main.R")
source("plots_supplementary.R")


## print numbers for manuscript ------------------------------------------------
summary(mill$year) # years
nrow(mill)
table(mill$coordinateUncertaintyInMeters)
# number of checklists with 1 km resolution or better
length(unique(mill$checklist_ID[mill$coordinateUncertaintyInMeters <= 1000]))

#*#*# not used right now
# number of checklists with 10 km resolution or better
# length(unique(mill$checklist_ID[mill$coordinateUncertaintyInMeters <= 10000]))


## end print numbers for manuscript -------------------------------------------



make_sampling_plots <- F # map sampling coverage in env. and geographic space
if(make_sampling_plots) source("sampling_coverage_maps.R")



save.image("millipede_maps_sonic.RData")
