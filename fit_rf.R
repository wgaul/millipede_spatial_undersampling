#############################
## Fit random forest SDMs for millipedes
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.  It
## relies on the outputs from 'prepare_objects_for_SDM.R'
## 
## TODO:  - make all testing be on spatially sub-sampled data (and use the
##          same spatially sub-sampled dataset to test all models)
##        - add list length model
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 23 Jan 2020
## last modified: 21 June 2020
##############################

on_sonic <- F

if(!dir.exists("./saved_objects")) dir.create("./saved_objects")

# load objects that are needed to fit SDMs
mill_wide <- readRDS("mill_wide.rds")
newdata <- readRDS("newdata.rds")
fold_assignments <- readRDS("fold_assignments.rds")
block_subsamp <- readRDS("block_subsamp.rds")

# keep only fold assignments for cross-validation sizes to run
kp <- sapply(fold_assignments, FUN = function(x, cv_size) {
  if(x$range %in% cv_size) {TRUE} else
    FALSE
}, cv_size = cv_block_sizes)
fold_assignments <- fold_assignments[kp] 

### fit random forest ---------------------------------------------------------
fit_rf <- function(test_fold, sp_name, sp_df, pred_names, newdata, 
                   sp_df_original, mtry = NULL, block_cv_range, pred_brick, 
                   block_range_spat_undersamp) {
  # ARGS: test_fold - integer giving the CV fold to be used as test data
  #       sp_name - character string giving the name of the column with 
  #                 detection/non-detection data to use
  #       sp_df - data frame to be used for model fitting (possibly 
  #               undersampled), with species observations and a column called 
  #               "folds" giving the cross-validation folds each record is in
  #       pred_names - character vector giving the variables to use as 
  #       predictors
  #       newdata - new data for predictions with standardized recording effort
  #       sp_df_original - the original (not undersampled) data
  if(is_tibble(sp_df)) {
    sp_df <- data.frame(sp_df)
  }
  sp_name <- gsub(" ", ".", sp_name)
  colnames(sp_df) <- gsub(" ", ".", colnames(sp_df))
 
  # label which observations are in the test fold
  sp_df$test_fold <- sp_df$folds == test_fold
  sp_df_original$test_fold <- sp_df_original$folds == test_fold

  if(is.null(mtry)) {
    mtry <- floor(sqrt(length(pred_names)))
  }
  
  # use 1000 trees which is hopefully high enough to get stable variable
  # importance if I want it (see manual linked in help documentation)
  mod <- tryCatch({
    randomForest(
      x = sp_df[sp_df$folds != test_fold, which(colnames(sp_df) %in% 
                                                  pred_names)],
      y = factor(sp_df[sp_df$folds != test_fold, which(colnames(sp_df) == 
                                                         sp_name)]),
      ntree = 1000, 
      mtry = mtry,   
      nodesize = 1, 
      replace = TRUE, classwt = NULL, 
      importance = TRUE, 
      keep.forest = TRUE)}, error = function(x) NA)

  # make predictions for model testing (predict to ALL folds, and predict to 
  # ALL locations for which there is data, not just the spatially 
  # undersampled locations)
  # For some reason, predict will fail if any columns have NAs, even if those
  # are not predictor columns.  Because in some rare cases a fold is not 
  # assigned to some cells (b/c they are on the boundary of a fold), there are
  # occasionally NAs in the "folds" column, so that column should not be passed
  # to predict.
  f_pred <- tryCatch({
    predict(mod, newdata = sp_df_original[, which(
      colnames(sp_df_original) != "folds")], 
            type = "prob")[, "1"]}, 
    error = function(x) NA)
  # select columns to keep in df of predictions
  preds <- sp_df_original[ , c("checklist_ID", "eastings", "northings", 
                               "hectad", "folds", "test_fold", sp_name)]
  preds$pred <- f_pred # add predictions to df

  # make dataframe of predictions with standardized recording effort
  # predict to ALL sites (both in training and test sets)
  newdata$pred <- tryCatch({
    predict(mod, newdata = newdata[, which(colnames(newdata) != "folds")], 
            type = "prob")[, "1"]}, 
    error = function(x) NA)
  
  # drop predictor variables from predictions dataframe 
  if(analysis_resolution == 10000) {
    newdata <- select(newdata, hectad, eastings, northings, 
                      day_of_year, pred)
  } else if(analysis_resolution == 1000) {
    newdata <- select(newdata, eastings, northings, day_of_year, pred)
    newdata$en <- paste0(round(newdata$eastings), "_", 
                         round(newdata$northings))
    # summarise standardized predictions over the entire year
    newdata <- group_by(newdata, en) %>%
      summarise(mean_pred = mean(pred), eastings = mean(eastings), 
                northings = mean(northings))
  } else stop("Analysis resolution must be either 10000 or 1000")


  
  ## calculate class balance (proportion of checklists with a detection)
  prop_dets <- tryCatch({
    length(which(sp_df[sp_df$folds != test_fold, 
                                  which(colnames(sp_df) == sp_name)] != 0)) / 
    nrow(sp_df)}, error = function (x) NA)
  
  ## calculate Simpson's evenness for training and test datasets
  # Only calculate this for fold # 1.  The same dataset is used in multiple 
  # folds, so the value will be identical for all folds using that dataset.
  if(test_fold == 1) {
    if(analysis_resolution == 10000){
      # calculate at hectad scale
      # calculate number of checklists per hectad
      table_nobs <- data.frame(table(sp_df$hectad))
      colnames(table_nobs) <- c("hectad", "nrec")
      # get all hectads (use newdata df to get names of all hectads)
      all_hectads <- newdata[newdata$day_of_year == newdata$day_of_year[1], ]
      all_hectads <- left_join(all_hectads, table_nobs, by = "hectad")
      all_hectads$nrec[is.na(all_hectads$nrec)] <- 0
      simps_train_hec <- simpson_even(as.numeric(all_hectads$nrec))
    } else simps_train_hec <- NA
   
    ## calculate at spatial subsampling block scale (e.g. 30k X 30km)
    lst_spat <- SpatialPointsDataFrame(
      coords = sp_df[, c("eastings", "northings")], 
      data = sp_df, proj4string = CRS("+init=epsg:29903"))
    lst_spat <- lst_spat[ , "checklist_ID"] # make df of checklists
    blks <- spatialBlock(lst_spat, 
                         theRange = block_range_spat_undersamp,
                         k = n_folds, selection = "random", iteration = 5, 
                         xOffset = runif(n = 1, min = 0, max = 1), 
                         yOffset = runif(n = 1, min = 0, max = 1),
                         showBlocks = FALSE, 
                         rasterLayer = pred_brick$pasture_l2, 
                         biomod2Format = FALSE)
    # add spatial subsampling grid cell ID to each hectad
    all_hectads <- hec_names
    all_hectads <- SpatialPointsDataFrame(
      coords = all_hectads[, c("eastings", "northings")], 
      data = all_hectads, proj4string = CRS("+init=epsg:29903"))
    all_hectads <- st_join(st_as_sf(all_hectads), 
                           st_as_sf(blks$blocks[, names(
                             blks$blocks) == "layer"]))
    all_hectads <- data.frame(all_hectads)
    colnames(all_hectads)[which(
      colnames(all_hectads) == "layer")] <- "subsamp_blocks"
    # remove geometry column
    all_hectads <- all_hectads[, -which(grepl(".*geomet.*", 
                                              colnames(all_hectads)))]
    table_nobs <- data.frame(table(sp_df$hectad)) # count nrec per hectad
    colnames(table_nobs) <- c("hectad", "nrec")
    # add number of recs
    all_hectads <- left_join(all_hectads, table_nobs, by = "hectad")
    all_hectads$nrec[is.na(all_hectads$nrec)] <- 0
    # sum number of recs per subsample block
    all_hectads <- group_by(all_hectads, subsamp_blocks) %>%
      summarise(nrec = sum(nrec))
    simps_train_subsampBlock <- simpson_even(as.numeric(all_hectads$nrec))
    rm(all_hectads, blks, lst_spat)
  } else {
    simps_train_hec <- NA
    simps_train_subsampBlock <- NA
  }
  
  
  # return fitted model, and predictions for this model
  tryCatch(list(
    m = mod, preds = preds, standardized_preds = newdata, 
    train_sites = unique(newdata$hectad[newdata$folds != test_fold]), 
    test_sites = unique(newdata$hectad[newdata$folds == test_fold]), 
    block_cv_range = block_cv_range, 
    n_detections_test_fold = try(sum(sp_df[sp_df$folds == test_fold, 
                                           sp_name])),
    simpson_training_hectad = simps_train_hec, 
    simpson_training_subsampBlock = simps_train_subsampBlock, 
    proportion_detections = prop_dets), 
    error = function(x) "No list exported from fit_rf.")
}


call_fit_rf <- function(fold_assignments, sp_name, test_fold, sp_df, 
                        pred_names, spatial.under.sample, block_subsamp, 
                        mtry = NULL, pred_brick, block_range_spat_undersamp, 
                        newdata, ...) {
  # Function to call fit_rf on each species in a list of species dfs
  # First spatially sub-sample non-detection records to even absence records
  # and improve class balance (presuming sp. is rare)  
  # Then fit RF.
  #
  # ARGS: fold_assignments - object of class SpatialBlock with the 
  #         cross-validation fold that each observation (row) of sp_df is
  #         assigned to
  #       sp_name - character string giving the name of the column with species
  #           detection/non-detection data
  #       test_fold - integer givin the CV fold to be used for testing
  #       sp_df - data frame with species observations and predictor variables
  #       pred_names - character vector with names of columns to be used as 
  #           predictor variables
  #       spatial.under.sample - T/F indicating whether to perform spatial
  #           under sampling (Robinson et al. 2018).  If TRUE, the spatial
  #           undersampling grid cells must be provided as a column named
  #           "spat_subsamp_cell"

  # add CV fold assignments to sp_df
  if(class(fold_assignments) == "SpatialBlock") {
    sp_df <- st_join(
      st_as_sf(sp_df), 
      st_as_sf(fold_assignments$blocks[, names(fold_assignments$blocks) ==
                                         "folds"]))
  }
  if(class(fold_assignments) == "list") { # add random CV blocks
    sp_df$folds <- fold_assignments$blocks
  }
  
  # convert to df from spatial
  if(is_tibble(sp_df) | "sf" %in% class(sp_df) | 
     "SpatialPointsDataFrame" %in% class(sp_df)) {
    sp_df <- data.frame(sp_df)
  }
  colnames(sp_df) <- gsub(" ", ".", colnames(sp_df))
  sp_name <- gsub(" ", ".", sp_name)
 
  # convert species record counts to p/a
  sp_df[colnames(sp_df) == sp_name] <- pa(sp_df[colnames(sp_df) == sp_name])
  sp_df_original <- sp_df
  
  if(spatial.under.sample) {
    # add a column of subsampling block assignments by chosing randomly from
    # the many allocations created in "prepare_objects_for_SDM.R"
    sp_df$spat_subsamp_cell <- block_subsamp[, sample(2:(ncol(block_subsamp)-1), 
                                                         size = 1)]
    # spatially sub-sample absence checklists to 1 per cell
    # separate presence and absence checklists.  Keep all presence checklists.
    presences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 1, ]
    absences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 0, ]
    cell_abs_tab <- table(absences$spat_subsamp_cell)
    keep_ab_rows <- c()
    for(ri in 1:length(unique(absences$spat_subsamp_cell))) {
      cell <- unique(absences$spat_subsamp_cell)[ri]
      keep_ab_rows <- c(keep_ab_rows, 
                        sample(which(absences$spat_subsamp_cell == cell), 
                               size = 1))
    }
    absences <- absences[keep_ab_rows, ] 
    
    # combine spatially sub-sampled non-detection data with all detection data
    sp_df <- bind_rows(absences, presences)
  }

  # if using block CV, assign newdata grid cells to CV fold
  # if using random CV, the folds pertain to checklists, not locations, so all
  # newdata locations are essentially in test set
  if(class(fold_assignments) == "SpatialBlock") {
    newdata <- st_join(
      st_as_sf(newdata), 
      st_as_sf(fold_assignments$blocks[, names(fold_assignments$blocks) == 
                                         "folds"]))
  }
  newdata <- data.frame(newdata)
  
  # fit RF
  lapply(1:n_folds, FUN = fit_rf, sp_name = sp_name, 
         sp_df = sp_df, pred_names = pred_names, newdata = newdata, 
         sp_df_original = sp_df_original, mtry = mtry, 
         block_cv_range = fold_assignments$range, pred_brick = pred_brick,
         block_range_spat_undersamp = block_range_spat_undersamp)
}


#### Fit Random Forest -------------------------------------------------------
if("day_ll_rf" %in% mod_names) {
  ### fit Day of Year + List Length models ------------------------------------
  # train with raw data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    day_ll_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf, 
      sp_df = mill_wide, 
      pred_names = c("sin_doy", "cos_doy", "list_length"), 
      block_subsamp = block_subsamp, 
      block_range_spat_undersamp = block_range_spat_undersamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = FALSE, 
      mtry = 1, 
      mc.cores = n_cores)
    try(print(pryr::object_size(day_ll_rf_fits)))
    try(saveRDS(day_ll_rf_fits, paste0("./saved_objects/", 
                                       "day_ll_rf_noSubSamp_fits_", 
                                       gsub(" ", "_", sp_name), 
                                       analysis_resolution, ".rds")))
    rm(day_ll_rf_fits)
  }
  
  # train with spatially subsampled data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    day_ll_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf, 
      sp_df = mill_wide, 
      pred_names = c("sin_doy", "cos_doy", "list_length"), 
      block_subsamp = block_subsamp, 
      block_range_spat_undersamp = block_range_spat_undersamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = TRUE, 
      mtry = 1,
      mc.cores = n_cores)
    try(print(pryr::object_size(day_ll_rf_fits)))
    try(saveRDS(day_ll_rf_fits, paste0("./saved_objects/", 
                                       "day_ll_rf_SubSamp_fits_", 
                                       gsub(" ", "_", sp_name),
                                       analysis_resolution, ".rds")))
    rm(day_ll_rf_fits)
  }
} ### end Day of Year + List Length ------------------------------------------

if("spat_ll_rf" %in% mod_names) {
  ### fit Spatial + List Length + DOY models -----------------------------------
  # train with raw data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    spat_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf, 
      sp_df = mill_wide, 
      pred_names = c("eastings", "northings", 
                     "sin_doy", "cos_doy", "list_length"),
      block_subsamp = block_subsamp, 
      block_range_spat_undersamp = block_range_spat_undersamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = FALSE, 
      mtry = 2, 
      mc.cores = n_cores)
    try(print(pryr::object_size(spat_rf_fits)))
    try(saveRDS(spat_rf_fits, paste0("./saved_objects/", 
                                     "spat_ll_rf_noSubSamp_fits_", 
                                     gsub(" ", "_", sp_name),
                                     analysis_resolution, ".rds")))
    rm(spat_rf_fits)
  }
  
  # train with spatially undersampled data 
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    spat_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf, 
      sp_df = mill_wide, 
      pred_names = c("eastings", "northings", 
                     "sin_doy", "cos_doy", "list_length"),
      block_subsamp = block_subsamp, 
      block_range_spat_undersamp = block_range_spat_undersamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = TRUE, 
      mtry = 2, 
      mc.cores = n_cores)
    try(print(pryr::object_size(spat_rf_fits)))
    try(saveRDS(spat_rf_fits, paste0("./saved_objects/", 
                                     "spat_ll_rf_SubSamp_fits_", 
                                     gsub(" ", "_", sp_name),
                                     analysis_resolution,
                                     ".rds")))
    rm(spat_rf_fits)
  }
} ### end spatial + List Length + DOY models ---------------------------------


if("env_ll_rf" %in% mod_names) {
  ### fit environmental + LL + DOY model --------------------------------------
  # train with raw data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    env_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf, 
      sp_df = mill_wide, 
      pred_names = c(sp_predictors[[which(names(sp_predictors) == sp_name)]],
                     "sin_doy", "cos_doy", "list_length"),
      block_subsamp = block_subsamp, 
      block_range_spat_undersamp = block_range_spat_undersamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = FALSE, 
      mc.cores = n_cores)
    try(print(pryr::object_size(env_rf_fits)))
    try(saveRDS(env_rf_fits, paste0("./saved_objects/", 
                                    "env_ll_rf_noSubSamp_fits_", 
                                    gsub(" ", "_", sp_name),
                                    analysis_resolution, ".rds")))
    rm(env_rf_fits)
  }
  
  # train with spatially subsampled data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    env_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf, 
      sp_df = mill_wide, 
      pred_names = c(sp_predictors[[which(names(sp_predictors) == sp_name)]],
                     "sin_doy", "cos_doy", "list_length"),
      block_subsamp = block_subsamp, 
      block_range_spat_undersamp = block_range_spat_undersamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = TRUE, 
      mc.cores = n_cores)
    try(print(pryr::object_size(env_rf_fits)))
    try(saveRDS(env_rf_fits, paste0("./saved_objects/", 
                                    "env_ll_rf_SubSamp_fits_", 
                                    gsub(" ", "_", sp_name),
                                    analysis_resolution, ".rds")))
    rm(env_rf_fits)
  }
} ### end environmental + LL + DOY model --------------------------------------

if("env_spat_ll_rf" %in% mod_names) {
  ### fit environmental + lat + long + LL + DOY model -------------------------
  # train with raw data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    env_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf, 
      sp_df = mill_wide, 
      pred_names = c(sp_predictors[[which(names(sp_predictors) == sp_name)]],
                     "eastings", "northings", 
                     "sin_doy", "cos_doy", "list_length"),
      block_subsamp = block_subsamp, 
      block_range_spat_undersamp = block_range_spat_undersamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = FALSE, 
      mc.cores = n_cores)
    try(print(pryr::object_size(env_rf_fits)))
    try(saveRDS(env_rf_fits, paste0("./saved_objects/", 
                                    "env_spat_ll_rf_noSubSamp_fits_", 
                                    gsub(" ", "_", sp_name),
                                    analysis_resolution, ".rds")))
    rm(env_rf_fits)
  }
  
  # train with spatially subsampled data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    env_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf, 
      sp_df = mill_wide, 
      pred_names = c(sp_predictors[[which(names(sp_predictors) == sp_name)]],
                     "eastings", "northings", 
                     "sin_doy", "cos_doy", "list_length"),
      block_subsamp = block_subsamp, 
      block_range_spat_undersamp = block_range_spat_undersamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = TRUE, 
      mc.cores = n_cores)
    try(print(pryr::object_size(env_rf_fits)))
    try(saveRDS(env_rf_fits, paste0("./saved_objects/", 
                                    "env_spat_ll_rf_SubSamp_fits_", 
                                    gsub(" ", "_", sp_name),
                                    analysis_resolution, ".rds")))
    rm(env_rf_fits)
  }
} ### end environmental + Lat + Lon + LL + DOY model --------------------------
### end fit random forest ----------------------------------------------------


### Get elements of the fitted models for results -----------------------------
for(mod_name in mod_names) {
  for(i in 1:length(sp_to_fit)) {
    ### trained with raw data
    fits <- readRDS(paste0("./saved_objects/", mod_name, "_noSubSamp_fits_", 
                           gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                           ".rds"))
    sp_name <- names(sp_to_fit)[i]
    
    ## Get predictions to actual data
    predictions <- lapply(fits, FUN = function(x) {
      lapply(x[sapply(x, is.list)], FUN = function(x) {
        preds <- tryCatch(x$preds, error = function(x) NA)
        preds$cv <- as.character(x$block_cv_range)
        return(preds)})})
    predictions <- bind_rows(lapply(
      predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
    # save results
    try(saveRDS(predictions, 
                paste0("./saved_objects/", "actual_predictions_", 
                       mod_name, "_noSubSamp_", 
                       gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                       ".rds")))
    
    ## Get predictions with standardized survey effort
    mill_predictions <- lapply(fits, FUN = function(x) {
      lapply(x[sapply(x, is.list)], FUN = function(x) {
        preds <- tryCatch(x$standardized_preds, error = function(x) NA)
        preds$cv <- as.character(x$block_cv_range)
        return(preds)})})
    mill_predictions <- bind_rows(lapply(
      mill_predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
    # save results
    try(saveRDS(mill_predictions, 
                paste0("./saved_objects/", "standard_predictions_", 
                       mod_name, "_noSubSamp_", 
                       gsub(" ", "_", sp_name), analysis_resolution, 
                       ".rds")))
    
    ## get variable importance (averaged over 3 folds) 
    var_imp <- lapply(fits, FUN = function(x) {
      bind_rows(lapply(x[sapply(x, is.list)], FUN = function(y) {
        df <- data.frame(y$m$importance)
        df$variable <- rownames(df)
        df[, c("variable", "MeanDecreaseGini")]
        df$species <- sp_name
        df$model <- mod_name
        df$train_dat <- "raw"
        df$cv <- as.character(y$block_cv_range)
        return(df)}))}
    )
    var_imp <- bind_rows(var_imp)
    # save results
    try(saveRDS(var_imp, 
                paste0("./saved_objects/", 
                       "var_importance_", mod_name, "_noSubSamp_", 
                       gsub(" ", "_", sp_name), analysis_resolution, 
                       ".rds")))
    
    if(get_partial_dependence & mod_name %in% mods_for_pd_plots) {
      ## get partial dependence
      # have to do this in for loops because when using lapply it seems R can't
      # interpret a variable holding the name of the variable to get plot for
      # in the partialPlot call.
      pd_list <- list()
      # loop through versions of this model for this species
      for(fi in 1:length(fits)) { 
        fits_pd <- list()
        for(xi in 1:length(fits[[fi]])) { # loop through all cv folds
          mod <- fits[[fi]][[xi]]$m
          vars <- as.character(rownames(importance(mod))) # get variable names
          # get partial dependence for each variable
          pl <- list()
          for(vb in 1:length(vars)) {
            pd <- data.frame(partialPlot(mod, 
                                         pred.data = data.frame(mill_wide), 
                                         x.var = vars[vb], 
                                         which.class = "1", plot = FALSE))
            pd$variable <- vars[vb]
            pd$species <- sp_name
            pd$model <- mod_name
            pd$train_data <- "raw"
            pd$cv <- as.character(fits[[fi]][[xi]]$block_cv_range)
            pl[[vb]] <- pd
          }
          pl <- bind_rows(pl)
         
          # add partial dependence for day of year
          ## calculate partial dependence for that same mean_rr variable by hand
          # This is my best guess about how to do this based on eq. 53 from that 
          # Friedman (2001) article
          
          # make day of year values to predict to
          doy_vals <- seq(from = 1, to = 365, length.out = 50)
          # create a vector to hold the mean prediction for each doy value
          dependence <- c() 
          
          for(dv in 1:length(doy_vals)) {
            # copy the training data to use for getting predictions
            ndat <- data.frame(mill_wide) 
            # fix day of year values to the value we are testing
            ndat$day_of_year <- doy_vals[dv] 
            # calculate sin and cos of day of year
            ndat$sin_doy <- sin((2*pi*ndat$day_of_year) / 365) 
            ndat$cos_doy <- cos((2*pi*ndat$day_of_year) / 365)
            # get predicted values for each case in the training data but 
            # with the day of year (and cos_doy and sin_doy) values fixed to 
            # the value that we are testing
            probs <- as.numeric(as.character(predict(mod, newdata = ndat, 
                                                     type = "prob")[, "1"]))
            
            # make any probabilities of zero be slightly positive to prevent 
            # problems with taking the log
            probs[probs == 0] <- 0.0001
            
            # get the partial dependence following eq 10.52 of Elements of 
            # Stat. Learning
            dependence[dv] <- mean(log(probs)) - 
              ((mean(log(probs)) + mean(log(1-probs)))/2)
          }
          
          # add dependence for DOY to df with other dependences
          pd_doy <- data.frame(
            x = doy_vals, y = dependence, 
            variable = "day_of_year", species = sp_name, 
            model = mod_name, train_data = "raw", 
            cv = as.character(fits[[fi]][[xi]]$block_cv_range))
          pl <- bind_rows(pl, pd_doy)
          
          # put partial dependence df in list
          fits_pd[[xi]] <- pl
        }
        pd_list[[fi]] <- bind_rows(fits_pd)
      }
      pd_list <- bind_rows(pd_list)
      # save results
      try(saveRDS(pd_list, 
                  paste0("./saved_objects/", 
                         "partial_dependence_", mod_name, "_noSubSamp_", 
                         gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                         ".rds")))
    } # end get partial dependence 
    rm(fits)
    
    ######### trained with spatial subsampling #############################
    fits <- readRDS(paste0("./saved_objects/", mod_name, "_SubSamp_fits_", 
                           gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                           ".rds"))
    
    ## Get predictions to actual data
    predictions <- lapply(fits, FUN = function(x) {
      lapply(x[sapply(x, is.list)], FUN = function(x) {
        preds <- tryCatch(x$preds, error = function(x) NA)
        preds$cv <- as.character(x$block_cv_range)
        return(preds)})})
    predictions <- bind_rows(lapply(
      predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
    # save results
    try(saveRDS(predictions, 
                paste0("./saved_objects/", "actual_predictions_", 
                       mod_name, "_SubSamp_", 
                       gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                       ".rds")))
    
    # Get predictions with standardized survey effort
    mill_predictions <- lapply(fits, FUN = function(x) {
      lapply(x[sapply(x, is.list)], FUN = function(x) {
        preds <- tryCatch(x$standardized_preds, error = function(x) NA)
        preds$cv <- as.character(x$block_cv_range)
        return(preds)})})
    mill_predictions <- bind_rows(lapply(
      mill_predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
    # save results
    try(saveRDS(mill_predictions, 
                paste0("./saved_objects/", 
                       "standard_predictions_", mod_name, "_SubSamp_", 
                       gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                       ".rds")))
    
    # get variable importance (averaged over 5 folds) 
    var_imp <- lapply(fits, FUN = function(x) {
      bind_rows(lapply(x[sapply(x, is.list)], FUN = function(y) {
        df <- data.frame(y$m$importance)
        df$variable <- rownames(df)
        df[, c("variable", "MeanDecreaseGini")]
        df$species <- sp_name
        df$model <- mod_name
        df$train_dat <- "spat_subsamp"
        df$cv <- as.character(y$block_cv_range)
        return(df)}))}
    )
    var_imp <- bind_rows(var_imp)
    # save results
    try(saveRDS(var_imp, 
                paste0("./saved_objects/", 
                       "var_importance_", mod_name, "_SubSamp_", 
                       gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                       ".rds")))
    
    if(get_partial_dependence & mod_name %in% mods_for_pd_plots) {
      ## get partial dependence
      # have to do this in for loops because when using lapply it seems R can't
      # interpret a variable holding the name of the variable to get plot for
      # in the partialPlot call.
      pd_list <- list()
      # loop through versions of this model for this species
      for(fi in 1:length(fits)) { 
        fits_pd <- list()
        for(xi in 1:length(fits[[fi]])) { # loop through all cv folds
          mod <- fits[[fi]][[xi]]$m
          vars <- as.character(rownames(importance(mod))) # get variable names
          # get partial dependence for each variable
          pl <- list()
          for(vb in 1:length(vars)) {
            pd <- data.frame(partialPlot(mod, 
                                         pred.data = data.frame(mill_wide), 
                                         x.var = vars[vb], 
                                         which.class = "1", plot = FALSE))
            pd$variable <- vars[vb]
            pd$species <- sp_name
            pd$model <- mod_name
            pd$train_data <- "spat_subsamp"
            pd$cv <- as.character(fits[[fi]][[xi]]$block_cv_range)
            pl[[vb]] <- pd
          }
          pl <- bind_rows(pl)
          
          # add partial dependence for day of year
          ## calculate partial dependence for that same mean_rr variable by hand
          # This is my best guess about how to do this based on eq. 53 from that
          # Friedman (2001) article
          
          # make day of year values to predict to
          doy_vals <- seq(from = 1, to = 365, length.out = 50)
          # create a vector to hold the mean prediction for each doy value
          dependence <- c() 
          
          for(dv in 1:length(doy_vals)) {
            # copy the training data to use for getting predictions
            ndat <- data.frame(mill_wide) 
            # fix day of year values to the value we are testing
            ndat$day_of_year <- doy_vals[dv] 
            # calculate sin and cos of day of year
            ndat$sin_doy <- sin((2*pi*ndat$day_of_year) / 365) 
            ndat$cos_doy <- cos((2*pi*ndat$day_of_year) / 365)
            # get predicted values for each case in the training data but 
            # with the day of year (and cos_doy and sin_doy) values fixed 
            # to the value that we are testing
            probs <- as.numeric(as.character(predict(mod, newdata = ndat, 
                                                     type = "prob")[, "1"]))
            
            # make any probabilities of zero be slightly positive to prevent 
            # problems with taking the log
            probs[probs == 0] <- 0.0001
            
            # get the partial dependence following eq 10.52 of Elements of 
            # Stat. Learning
            dependence[dv] <- mean(log(probs)) - 
              ((mean(log(probs)) + mean(log(1-probs)))/2)
          }
          
          # add dependence for DOY to df with other dependences
          pd_doy <- data.frame(
            x = doy_vals, y = dependence, 
            variable = "day_of_year", species = sp_name, 
            model = mod_name, train_data = "spat_subsamp", 
            cv = as.character(fits[[fi]][[xi]]$block_cv_range))
          pl <- bind_rows(pl, pd_doy)
          
          fits_pd[[xi]] <- pl
        }
        pd_list[[fi]] <- bind_rows(fits_pd)
      }
      pd_list <- bind_rows(pd_list)
      # save results
      try(saveRDS(pd_list, 
                  paste0("./saved_objects/", 
                         "partial_dependence_", mod_name, "_SubSamp_", 
                         gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                         ".rds")))
    } # end get partial dependence plots
    try(rm(fits, pd_list, var_imp, mill_predictions))
  }
}

print(Sys.time())

if(on_sonic) quit(save = "no")
