#############################
## Fit random forest SDMs for millipedes using randomly undersampled data
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.  It
## relies on the outputs from 'prepare_objects_for_SDM.R'
## 
## This script was added in April 2022 in response to reviewer comments
## to test random undersampling, in order to determine whether spatial 
## undersampling provides a benefit above and beyond a non-spatial undersampling.
## 
## TODO:  
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 11 April 2022
## last modified: 14 April 2022
##############################

if(!dir.exists("./saved_objects")) dir.create("./saved_objects")

# load objects that are needed to fit SDMs
mill_wide <- readRDS("mill_wide.rds")
newdata <- readRDS("newdata.rds")
fold_assignments <- readRDS("fold_assignments.rds")

### get mean proportion of checklists that are detection lists for each species 
## in the spatially undersampled data
# I am hard-coding this right now to save time while working on reviews.  
# I calculated these on 14 April 2022 from the results used for the initial
# submission to the journal.
# It would be possible (and better) to generate this in the script by getting 
# the information from the fitted models objects.
class_balance_df <- data.frame(
  species = c("Macrosternodesmus palicola", "Boreoiulus tenuis",
              "Ommatoiulus sabulosus", "Blaniulus guttulatus",
              "Glomeris marginata", "Cylindroiulus punctatus"), 
  mean_proportion_detections = c(0.219, 0.246, 0.260, 0.359, 0.486, 0.557))


# keep only fold assignments for cross-validation sizes to run
kp <- sapply(fold_assignments, FUN = function(x, cv_size) {
  if(x$range %in% cv_size) {TRUE} else
    FALSE
}, cv_size = cv_block_sizes)
fold_assignments <- fold_assignments[kp] 

### fit random forest ---------------------------------------------------------
fit_rf_random_undersample <- function(test_fold, sp_name, sp_df, pred_names, 
                                      newdata, sp_df_original, mtry = NULL, 
                                      block_cv_range, pred_brick) {
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
  preds <- sp_df_original[ , c("checklist_ID", "eastings", 
                               "northings", "folds", "hectad", 
                               "test_fold", sp_name)]
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
                      month, pred)
  } else if(analysis_resolution == 1000) {
    newdata <- select(newdata, eastings, northings, 
                      month, pred)
    newdata$en <- paste0(round(newdata$eastings), "_", 
                         round(newdata$northings))
    # summarise standardized predictions over the entire year
    newdata <- group_by(newdata, en) %>%
      summarise(mean_pred = mean(pred), 
                eastings = mean(eastings), 
                northings = mean(northings))
  } else stop("Analysis resolution must be either 10000 or 1000")
  
 
  ## calculate class balance (proportion of checklists with a detection) for
  ## data IN THIS TRAINING FOLD
  prop_dets <- tryCatch({
    length(which(sp_df[sp_df$folds != test_fold, 
                       which(colnames(sp_df) == sp_name)] != 0)) / 
      nrow(sp_df[sp_df$folds != test_fold, ])}, error = function (x) NA)
  
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
      all_hectads <- newdata[newdata$month == newdata$month[1], ]
      all_hectads <- left_join(all_hectads, table_nobs, by = "hectad")
      all_hectads$nrec[is.na(all_hectads$nrec)] <- 0
      simps_train_hec <- simpson_even(as.numeric(all_hectads$nrec))
    } else simps_train_hec <- NA
  }
  
  # return fitted model, and predictions for this model
  tryCatch(list(
    m = mod, preds = preds, standardized_preds = newdata, 
    train_sites = tryCatch(unique(newdata$hectad[newdata$folds != test_fold]),
                           error = function(x) NA),
    test_sites = tryCatch(unique(newdata$hectad[newdata$folds == test_fold]), 
                          error = function(x) NA),
    n_detections_test_fold = try(sum(sp_df[sp_df$folds == test_fold, 
                                           sp_name])),
    simpson_training_hectad = tryCatch(simps_train_hec, error = function(x) NA), 
    proportion_detections = prop_dets, 
    data_for_this_model = sp_df[, colnames(sp_df) %in% 
                                  c("checklist_ID", "folds", "test_fold")]), 
    error = function(x) "No list exported from fit_rf_random_undersample")
}

call_fit_rf_random_undersample <- function(fold_assignments, sp_name, test_fold, 
                                           sp_df, pred_names, 
                                           class_balance_df,
                                           random.under.sample,  mtry = NULL, 
                                           pred_brick,  newdata, ...) {
  # Function to call fit_rf_random_undersample on each species in a list of 
  # species dfs
  # First randomly sub-sample non-detection records to improve class balance 
  # Then fit RF.
  #
  # ARGS: fold_assignments - object of class SpatialBlock with the 
  #         cross-validation fold that each observation (row) of sp_df is
  #         assigned to
  #       sp_name - character string giving the name of the column with species
  #           detection/non-detection data
  #       test_fold - integer giving the CV fold to be used for testing
  #       sp_df - data frame with species observations and predictor variables
  #       pred_names - character vector with names of columns to be used as 
  #           predictor variables
  #       class_balance_df - a data frame with the mean proportion of checklists 
  #           that are detection checklists in the spatially undersampled data 
  #           for each species
  #       random.under.sample - T/F indicating whether to perform random
  #           under sampling.

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

  if(random.under.sample) {
    # sub-sample absence checklists to match the average class balance from
    # checklists in the spatially undersampled datasets for this species
    # separate presence and absence checklists.  Keep all presence checklists.
    presences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 1, ]
    absences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 0, ]

    class_balance_df$species <- gsub(" ", ".", class_balance_df$species)
    
    # what proportion of data are detections in spatially undersampled data for
    # this species?
    this_prop_dets <- class_balance_df$mean_proportion_detections[which(
      class_balance_df$species == sp_name)]
    # calculate how many non-detections are needed to match that proportion 
    # of detections
    n_nondets_needed <- (nrow(presences) - 
      (nrow(presences) * this_prop_dets)) / this_prop_dets

    absences <- absences[sample(1:nrow(absences), size = n_nondets_needed), ] 
    
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
  
  newdata <- data.frame(newdata) # must be df (nob-spatial) for fit_rf()
 
  # fit RF
  lapply(1:n_folds, FUN = fit_rf_random_undersample, sp_name = sp_name, 
         sp_df = sp_df, pred_names = pred_names, newdata = newdata, 
         sp_df_original = sp_df_original, mtry = mtry, pred_brick = pred_brick)
}


#### Fit Random Forest -------------------------------------------------------
if("month_ll_rf" %in% mod_names) {
  ### fit Month + List Length models ------------------------------------
  # train with randomly undersampled data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    month_ll_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf_random_undersample, 
      sp_df = mill_wide, 
      class_balance_df = class_balance_df,
      random.under.sample = TRUE, 
      pred_names = c("sin_month", "cos_month", "list_length"), 
      block_subsamp = block_subsamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = FALSE, 
      mtry = 1, 
      mc.cores = n_cores)
    try(print(pryr::object_size(month_ll_rf_fits)))
    try(saveRDS(month_ll_rf_fits, paste0("./saved_objects/", 
                                         "month_ll_rf_randomSubSamp_fits_", 
                                         gsub(" ", "_", sp_name), 
                                         analysis_resolution, ".rds")))
    rm(month_ll_rf_fits)
  }
} 
### end month of Year + List Length ------------------------------------------


if("spat_ll_rf" %in% mod_names) {
  ### fit Spatial + List Length + month models -----------------------------
  # train with randomly undersampled data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    spat_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf_random_undersample, 
      sp_df = mill_wide, 
      class_balance_df = class_balance_df,
      random.under.sample = TRUE, 
      pred_names = c("eastings", "northings", 
                     "sin_month", "cos_month", "list_length"), 
      block_subsamp = block_subsamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = FALSE, 
      mtry = 1, 
      mc.cores = n_cores)
    try(print(pryr::object_size(spat_rf_fits)))
    try(saveRDS(spat_rf_fits, paste0("./saved_objects/", 
                                         "spat_ll_rf_randomSubSamp_fits_", 
                                         gsub(" ", "_", sp_name), 
                                         analysis_resolution, ".rds")))
    rm(spat_rf_fits)
  }
} 
### end Spatial + List Length + month ------------------------------------------

if("env_ll_rf" %in% mod_names) {
  ### fit environmental + LL + month models ------------------------------------
  # train with randomly undersampled data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    env_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf_random_undersample, 
      sp_df = mill_wide, 
      class_balance_df = class_balance_df,
      random.under.sample = TRUE, 
      pred_names = c(sp_predictors[[which(names(sp_predictors) == sp_name)]],
                     "sin_month", "cos_month", "list_length"), 
      block_subsamp = block_subsamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = FALSE, 
      mtry = 1, 
      mc.cores = n_cores)
    try(print(pryr::object_size(env_rf_fits)))
    try(saveRDS(env_rf_fits, paste0("./saved_objects/", 
                                         "env_ll_rf_randomSubSamp_fits_", 
                                         gsub(" ", "_", sp_name), 
                                         analysis_resolution, ".rds")))
    rm(env_rf_fits)
  }
} 
### end environmental + LL + month ------------------------------------------


if("env_spat_ll_rf" %in% mod_names) {
  ### fit environmental + lat + long + LL + month models -----------------------
  # train with randomly undersampled data
  for(i in 1:length(sp_to_fit)) {
    sp_name <- names(sp_to_fit)[i]
    env_rf_fits <- mclapply(
      fold_assignments, 
      sp_name = sp_name, 
      FUN = call_fit_rf_random_undersample, 
      sp_df = mill_wide, 
      class_balance_df = class_balance_df,
      random.under.sample = TRUE, 
      pred_names = c(sp_predictors[[which(names(sp_predictors) == sp_name)]],
                     "eastings", "northings", 
                     "sin_month", "cos_month", "list_length"), 
      block_subsamp = block_subsamp, 
      pred_brick = pred_brick, newdata = newdata,
      spatial.under.sample = FALSE, 
      mtry = 1, 
      mc.cores = n_cores)
    try(print(pryr::object_size(env_rf_fits)))
    try(saveRDS(env_rf_fits, paste0("./saved_objects/", 
                                         "env_spat_ll_rf_randomSubSamp_fits_", 
                                         gsub(" ", "_", sp_name), 
                                         analysis_resolution, ".rds")))
    rm(env_rf_fits)
  }
} 
### end environmental + lat + long + LL + month -------------------------------
### end fit random forest ----------------------------------------------------

### Get elements of the fitted models for results -----------------------------
# TODO: Here 12 April
for(mod_name in mod_names) {
  for(i in 1:length(sp_to_fit)) {
    ######### trained with random undersampling #############################
    fits <- readRDS(paste0("./saved_objects/", mod_name, 
                           "_randomSubSamp_fits_", 
                           gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                           ".rds"))
    
    ## Get predictions to actual data
    predictions <- lapply(fits, FUN = function(x) {
      lapply(x[sapply(x, is.list)], FUN = function(x) {
        preds <- tryCatch(x$preds, error = function(x) NA)
        return(preds)})})
    predictions <- bind_rows(lapply(
      predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
    # save results
    try(saveRDS(predictions, 
                paste0("./saved_objects/", "actual_predictions_", 
                       mod_name, "_randomSubSamp_", 
                       gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                       ".rds")))
    
    # Get predictions with standardized survey effort
    mill_predictions <- lapply(fits, FUN = function(x) {
      lapply(x[sapply(x, is.list)], FUN = function(x) {
        preds <- tryCatch(x$standardized_preds, error = function(x) NA)
        return(preds)})})
    mill_predictions <- bind_rows(lapply(
      mill_predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
    # save results
    try(saveRDS(mill_predictions, 
                paste0("./saved_objects/", 
                       "standard_predictions_", mod_name, "_randomSubSamp_", 
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
        df$train_dat <- "random_subsamp"
        return(df)}))}
    )
    var_imp <- bind_rows(var_imp)
    # save results
    try(saveRDS(var_imp, 
                paste0("./saved_objects/", 
                       "var_importance_", mod_name, "_randomSubSamp_", 
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
            pd$train_data <- "random_subsamp"
            pl[[vb]] <- pd
          }
          pl <- bind_rows(pl)
          
          # add partial dependence for month
          ## calculate partial dependence for that same mean_rr variable by hand
          # This is my best guess about how to do this based on eq. 53 from that
          # Friedman (2001) article
          
          # make month values to predict to
          month_vals <- seq(from = 1, to = 12)
          # create a vector to hold the mean prediction for each month value
          dependence <- c() 
          
          for(dv in 1:length(month_vals)) {
            # copy the training data to use for getting predictions
            ndat <- data.frame(mill_wide) 
            # fix month values to the value we are testing
            ndat$month <- month_vals[dv] 
            # calculate sin and cos of month
            ndat$sin_month <- sin((2*pi*ndat$month) / 12) 
            ndat$cos_month <- cos((2*pi*ndat$month) / 12)
            # get predicted values for each case in the training data but 
            # with the month (and cos_month and sin_month) values fixed 
            # to the value that we are testing
            probs <- as.numeric(as.character(predict(mod, newdata = ndat, 
                                                     type = "prob")[, "1"]))
            
            # make any probabilities of zero be slightly positive 
            # and probabililties of one be slightly less 
            # to prevent problems with taking the log
            probs[probs == 0] <- 0.0001
            probs[probs == 1] <- 0.9999
            
            # get the partial dependence following eq 10.52 of Elements of 
            # Stat. Learning
            dependence[dv] <- mean(log(probs)) - 
              ((mean(log(probs)) + mean(log(1-probs)))/2)
            if(is.infinite(dependence[dv])) {
              warning(paste0("Infinite values in partial dependence calculations.  fi = ", fi, " , xi = ", xi))
            }
          }
          
          # add dependence for month to df with other dependences
          pd_month <- data.frame(
            x = month_vals, y = dependence, 
            variable = "month", species = sp_name, 
            model = mod_name, train_data = "spat_subsamp", 
            cv = as.character(fits[[fi]][[xi]]$block_cv_range))
          pl <- bind_rows(pl, pd_month)
          
          fits_pd[[xi]] <- pl
        }
        pd_list[[fi]] <- bind_rows(fits_pd)
      }
      pd_list <- bind_rows(pd_list)
      # save results
      try(saveRDS(pd_list, 
                  paste0("./saved_objects/", 
                         "partial_dependence_", mod_name, "_randomSubSamp_", 
                         gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                         ".rds")))
    } # end get partial dependence plots
    try(rm(fits, pd_list, var_imp, mill_predictions))
  }
}

try(rm(class_balance_df))
print(Sys.time())
