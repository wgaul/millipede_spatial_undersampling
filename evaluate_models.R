#############################
## Evaluate SDMs for millipedes
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
## 
## TODO:  - test on original data
##        - test on spatially sub-sampled data
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 12 May 2020
## last modified: 16 June 2020
##############################
library(pROC)
library(psych)

## Get spatially subsampled test points ---------------------------------------
test_points_ss <- sapply(sp_to_fit, FUN = function(x) NULL)
names(test_points_ss) <- gsub(" ", ".", names(test_points_ss))
mill_wide_df <- data.frame(mill_wide)
for(i in 1:length(test_points_ss)) {
  mill_wide_df[, colnames(mill_wide_df) == names(test_points_ss)[i]] <- pa(
    mill_wide_df[, colnames(mill_wide_df) == names(test_points_ss)[i]])
  
  test_points_ss[[i]] <- list()
  # make a large number of different spatially subsampled test points datasets 
  # to choose from so each test (i.e. each model in a CV fold) is tested with a 
  # different set.  Use the subsample block designations created earlier.
  for(j in 1:n_subsamp_block_draws) {
    # add a column of subsampling block assignments by chosing randomly from
    # the many allocations created in "prepare_objects_for_SDM.R"
    mill_wide_df$spat_subsamp_cell <- block_subsamp[, sample(
      2:(ncol(block_subsamp)-1), size = 1)]

    # separate presence and absence checklists.  Keep all presence checklists.
    presences <- mill_wide_df[mill_wide_df[colnames(mill_wide_df) == 
                                             names(test_points_ss)[i]] > 0, ]
    absences <- mill_wide_df[mill_wide_df[colnames(mill_wide_df) == 
                                            names(test_points_ss)[i]] == 0, ]
    
    # spatially sub-sample absence checklists to 2 per cell
    # Ideally, testing would use data that is spatially subsampled without 
    # separating detections and non-detections.  However, given the very small 
    # number of detections we have, in order to get enough detections in the 
    # test datasets I am opting to use the same "spatial undersampling" that 
    # I used for model training - i.e. subsampling non-detections but keeping 
    # all detections

    keep_ab_rows <- c()
    for(ri in 1:length(unique(absences$spat_subsamp_cell))) {
      cell <- unique(absences$spat_subsamp_cell)[ri]
      keep_ab_rows <- c(keep_ab_rows, 
                        sample(which(absences$spat_subsamp_cell == cell), 
                               size = 2))
    }
    absences <- absences[keep_ab_rows, ] 
    # combine spatially sub-sampled non-detection data with all detection data
    # and store this test dataset
    test_points_ss[[i]][[j]] <- bind_rows(absences, presences)
  }
}

# measure spatial evennes using Simpson's evenness for test datasets
evenness_test <- lapply(test_points_ss, FUN = function(x) {
  sapply(x, FUN = function(y) {
    simpson_even(as.numeric(table(y$hectad)))
  })
})
## end get spatially subsampled test points ----------------------------------

## evaluate models ------------------------------------------------------------
# trained with raw data
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0("./saved_objects/", mod_name, "_noSubSamp_fits_", 
                         gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                         ".rds"))
  sp_name <- names(sp_to_fit)[i]
  
  ## test against original (raw) data 
  # AUC --------------------------------------------------
  aucs <- lapply(fits, FUN = function(x, sp_name) {
    auc_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        as.numeric(pROC::roc(
          response = y$preds[y$preds$test_fold == T, 
                             colnames(y$preds) == gsub(" ", ".", sp_name)], 
          predictor = y$preds$pred[y$preds$test_fold == T])$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  simps_trains_hec <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$simpson_training_hectad}, error = function(x) NA)
    })})
  simps_trains_blk <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$simpson_training_subsampBlock}, error = function(x) NA)
    })})
  
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", metric = "AUC", 
                   value = unlist(aucs), 
                   block_cv_range = unlist(block_ranges), 
                   simpson_training_hectad = unlist(simps_trains_hec), 
                   simpson_training_subsampBlock = unlist(simps_trains_blk), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges, simps_trains_hec, simps_trains_blk) # end AUC ----------
  
  # Cohen's Kappa --------------------------------------
  kappa_calc <- function(x, resp, pred) {
    ## Function to calculate Cohen's Kappa using the best threshold
    # make a 2 column df with observed responses (0 or 1) and a T/F whether the 
    # predicted value is above threshold x.  These will be used to compute
    # Cohen's kappa for agreement between categorical response (0, 1) and 
    # categorical predictor (0, 1)
    # ARGS: x - cutoff to test
    #       resp - vector of observed responses (0 or 1)
    #       pred - numeric vector of continuous predictions
    vals <- data.frame(resp = factor(resp[!is.na(resp)]), 
                       pred = factor(as.numeric(pred[!is.na(pred)] > x)))
    psych::cohen.kappa(vals)$kappa
  }
  
  kappas <- seq(from = 0, to = 1, by = 0.1) # Thresholds to try
  names(kappas) <- kappas
  
  # get kappa for each fold
  kp <- lapply(fits, function(x, sp_name, kappas) {
    kp_iter <- lapply(x, FUN = function(y, sp_name, kappas) {
      tryCatch({
        k_res <- sapply(kappas, kappa_calc, 
                        resp = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == sp_name], 
                        pred = y$preds$pred[y$preds$test_fold == T])
        # return kappa and threshold that maximised kappa
        k_res[k_res == max(k_res)][1]}, 
        error = function(x) NA)
    }, sp_name = sp_name, kappas = kappas)
  }, sp_name = gsub(" ", ".", sp_to_fit[[i]]), kappas = kappas)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", metric = "Kappa", 
                   value = unlist(kp), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end Kappa ----------------------------------
  
  # sensitivity --------------------------------------------------
  # use threshold that maximised Kappa for each model
  sens <- mapply(FUN = function(x, thresh, sp_name) {
    sens_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        sensitivity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", metric = "sensitivity", 
                   value = unlist(sens), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end sensitivity --------------------------------------
  
  # specificity --------------------------------------------------
  # use threshold that maximised Kappa for each model
  specif <- mapply(FUN = function(x, thresh, sp_name) {
    specif_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        specificity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})

  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", metric = "specificity", 
                   value = unlist(specif), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end specificity --------------------------------------
  
  # Brier score --------------------------------------------------------------
  brier <- lapply(fits, FUN = function(x, sp_name) {
    br_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        pr <- y$preds$pred[y$preds$test_fold == T] # predicted values
        # observed values
        ob <- y$preds[y$preds$test_fold == T, 
                      colnames(y$preds) == gsub(" ", ".", sp_name)]
        # squared error
        sq_er <- (pr - ob)^2
        mean(sq_er) # return mean squared error (Brier score)
      }, 
      error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", metric = "Brier", 
                   value = unlist(brier), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  # end Brier score -----------------------------------------------------------
  ## end test on raw data -----------------------------------------------------
  
  ## test against spatially subsampled data -----------------------------------
  aucs <- lapply(fits, FUN = function(x, sp_name, test_points_ss) {
    auc_iter <- lapply(x, FUN = function(y, sp_name, test_points_ss) {
      # select test dataset to use for this test
      test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                               size = 1)]]
      # subset to spatially  subsampled test data points that are in the test 
      # fold for this model
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" or y$preds$test_fold == T are in different 
      # CV folds b/c cv block boundary is within the hectad and checklists are 
      # at finer spatial resolution, or because CV used randomly chosen 
      # checklists, regardless of what hectad the were in.  E.g. with random 
      # CV, if one checklist from a hectad is in the test fold and another is 
      # not, the hectad will be in "y$preds$hectad[y$preds$test_fold == T]" and
      # so the training checklist will be included by the previous line of code,
      # but will now be dropped by this following line of code.)
      test_data <- test_data[test_data$test_fold == T, ] 
      n_dets_nondets <- tryCatch({
        data.frame(table(test_data[, colnames(test_data) == gsub(" ", ".", 
                                                                 sp_name)]))}, 
                                 error = function(x) NA)
      auc <- tryCatch({
        as.numeric(pROC::roc(response = test_data[, colnames(test_data) == 
                                                    gsub(" ", ".", sp_name)], 
                             predictor = test_data$pred)$auc)}, 
        error = function(x) NA)
      list(auc = auc, n_dets_nondets = n_dets_nondets)
    }, sp_name = sp_name, test_points_ss = test_points_ss)
  }, sp_name = sp_to_fit[[i]], 
  test_points_ss = test_points_ss[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  simps_trains_hec <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$simpson_training_hectad}, error = function(x) NA)
    })})
  simps_trains_blk <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$simpson_training_subsampBlock}, error = function(x) NA)
    })})
  
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "spat_subsamp", metric = "AUC", 
                   value = unlist(lapply(aucs, FUN = function(x) {
                     lapply(x, FUN = function (y) {y$auc})})), 
                   block_cv_range = unlist(block_ranges), 
                   simpson_training_hectad = unlist(simps_trains_hec), 
                   simpson_training_subsampBlock = unlist(simps_trains_blk), 
                   proportion_detections = unlist(prop_dets), 
                   n_dets_in_test = unlist(lapply(aucs, FUN = function(x) {
                     lapply(x, FUN = function (y) {
                       if("1" %in% as.character(y$n_dets_nondets$Var1)) {
                         y$n_dets_nondets$Freq[y$n_dets_nondets$Var1 == "1"]
                       } else 0
                     })
                   })), 
                   n_nondets_in_test = unlist(lapply(aucs, FUN = function(x) {
                     lapply(x, FUN = function (y) {
                       if("0" %in% as.character(y$n_dets_nondets$Var1)) {
                         y$n_dets_nondets$Freq[y$n_dets_nondets$Var1 == "0"]
                       } else 0})
                     })))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges, simps_trains_hec, simps_trains_blk) # end AUC ----------
  
  # Kappa ------------------------------------
  # get kappa for each fold
  kp <- lapply(fits, function(x, sp_name, kappas, test_points_ss) {
    kp_iter <- lapply(x, FUN = function(y, sp_name, kappas, test_points_ss) {
      # select test dataset to use for this test
      test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                          size = 1)]]
      # subset to spatially  subsampled test data points that are in the test 
      # fold for this model
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" or y$preds$test_fold == T are actually in 
      # different CV folds.
      test_data <- test_data[test_data$test_fold == T, ]
      tryCatch({
        k_res <- sapply(kappas, kappa_calc, 
                        resp = test_data[, colnames(test_data) == sp_name], 
                        pred = test_data$pred)
        # return kappa and threshold that maximised kappa
        k_res[k_res == max(k_res)][1]}, 
        error = function(x) NA)
    }, sp_name = sp_name, kappas = kappas, test_points_ss = test_points_ss)
  }, sp_name = gsub(" ", ".", sp_to_fit[[i]]), kappas = kappas, 
  test_points_ss = test_points_ss[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "spat_subsamp", metric = "Kappa", 
                   value = unlist(kp), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end Kappa ----------------------------------
  
  # sensitivity --------------------------------------------------
  sens <- mapply(FUN = function(x, thresh, sp_name, test_points_ss) {
    sens_iter <- mapply(FUN = function(y, thresh, sp_name, test_points_ss) {
      # select test dataset to use for this test 
      test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                          size = 1)]]
      # subset to spatially  subsampled test data points that are in the test 
      # fold for this model
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      # join predictions to test data
      test_data <- data.frame(left_join(test_data, y$preds))
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" or y$preds$test_fold == T are actually in 
      # different CV folds.
      test_data <- test_data[test_data$test_fold == T, ]
      tryCatch({
        sensitivity(threshold = thresh, 
                    response = test_data[, colnames(test_data) == 
                                           gsub(" ", ".", sp_name)], 
                    predictions = test_data$pred)}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name, 
                                  test_points_ss = test_points_ss))
  }, fits, kp, 
  MoreArgs = list(sp_name = sp_to_fit[[i]], 
                  test_points_ss = test_points_ss[[i]]), 
  SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "spat_subsamp",  
                   metric = "sensitivity", 
                   value = unlist(sens), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end sensitivity --------------------------------------
  
  # specificity --------------------------------------------------
  specif <- mapply(FUN = function(x, thresh, sp_name, test_points_ss) {
    specif_iter <- mapply(FUN = function(y, thresh, sp_name, test_points_ss) {
      # select test dataset to use for this test
      test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                          size = 1)]]
      # subset to spatially  subsampled test data points that are in the test 
      # fold for this model
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" or y$preds$test_fold == T are actually in 
      # different CV folds.
      test_data <- test_data[test_data$test_fold == T, ]
      tryCatch({
        specificity(threshold = thresh, 
                    response = test_data[, colnames(test_data) == 
                                           gsub(" ", ".", sp_name)], 
                    predictions = test_data$pred)}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name, 
                                  test_points_ss = test_points_ss))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]], 
                               test_points_ss = test_points_ss[[i]]),
  SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "spat_subsamp", 
                   metric = "specificity", 
                   value = unlist(specif), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end specificity --------------------------------------
  
  # Brier score --------------------------------------------------------------
  brier <- lapply(fits, FUN = function(x, sp_name, test_points_ss) {
    br_iter <- lapply(x, FUN = function(y, sp_name, test_points_ss) {
      tryCatch({
        # select test dataset to use for this test
        test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                            size = 1)]]
        # subset to spatially  subsampled test data points that are in the test 
        # fold for this model
        test_data <- test_data[test_data$hectad %in% y$test_sites |
                                 test_data$hectad %in% 
                                 y$preds$hectad[y$preds$test_fold == T], ]
        test_data <- left_join(test_data, y$preds) # join predictions to test data
        # drop checklists not in the test fold (sometimes checlists in hectads 
        # that are in "y$test_sites" or y$preds$test_fold == T are in different 
        # CV folds b/c cv block boundary is within the hectad and checklists are 
        # at finer spatial resolution, or because CV used randomly chosen 
        # checklists, regardless of what hectad the were in. 
        test_data <- test_data[test_data$test_fold == T, ] 
        pr <- test_data$pred # predicted values
        # observed values
        ob <- test_data[, colnames(test_data) == gsub(" ", ".", sp_name)]
        sq_er <- (pr - ob)^2 # squared error
        mean(sq_er) # return mean squared error (Brier score)
      }, 
      error = function(x) NA)
    }, sp_name = sp_name, test_points_ss = test_points_ss)
  }, sp_name = sp_to_fit[[i]], test_points_ss = test_points_ss[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "spat_subsamp", metric = "Brier", 
                   value = unlist(brier), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  # end Brier score -----------------------------------------------------------
  ## end test with spatially subsampled data ----------------------------------
  rm(fits)
}

# trained with spatial subsampling
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0("./saved_objects/", mod_name, "_SubSamp_fits_", 
                         gsub(" ", "_", sp_to_fit[[i]]), analysis_resolution, 
                         ".rds"))
  sp_name <- names(sp_to_fit)[i]
  
  ## test against raw data -----------------------------------------------
  aucs <- lapply(fits, FUN = function(x, sp_name) {
    auc_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        as.numeric(pROC::roc(
          response = y$preds[y$preds$test_fold == T, 
                             colnames(y$preds) == gsub(" ", ".", sp_name)], 
          predictor = y$preds$pred[y$preds$test_fold == T])$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", metric = "AUC", 
                   value = unlist(aucs), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges, n_det_nondet) # end AUC ---------------------------------
  
  # Kappa ----------------------------------------------
  # get Kappa for each fold
  kp <- lapply(fits, function(x, sp_name, kappas) {
    kp_iter <- lapply(x, FUN = function(y, sp_name, kappas) {
      tryCatch({
        k_res <- sapply(kappas, kappa_calc, 
                        resp = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == sp_name], 
                        pred = y$preds$pred[y$preds$test_fold == T])
        # return kappa and threshold that maximised kappa
        k_res[k_res == max(k_res)][1]}, 
        error = function(x) NA)
    }, sp_name = sp_name, kappas = kappas)
  }, sp_name = gsub(" ", ".", sp_to_fit[[i]]), kappas = kappas)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", metric = "Kappa", 
                   value = unlist(kp), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end Kappa ----------------------------------
  
  # sensitivity --------------------------------------------------
  sens <- mapply(FUN = function(x, thresh, sp_name) {
    sens_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        sensitivity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", 
                   metric = "sensitivity", 
                   value = unlist(sens), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end sensitivity --------------------------------------
  
  # specificity --------------------------------------------------
  specif <- mapply(FUN = function(x, thresh, sp_name) {
    specif_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        specificity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", 
                   metric = "specificity", 
                   value = unlist(specif), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end specificity --------------------------------------
  
  # Brier score --------------------------------------------------------------
  brier <- lapply(fits, FUN = function(x, sp_name) {
    br_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        pr <- y$preds$pred[y$preds$test_fold == T] # predicted values
        # observed values
        ob <- y$preds[y$preds$test_fold == T, 
                      colnames(y$preds) == gsub(" ", ".", sp_name)]
        # squared error
        sq_er <- (pr - ob)^2
        mean(sq_er) # return mean squared error (Brier score)
      }, 
      error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", metric = "Brier", 
                   value = unlist(brier), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges)
  # end Brier score -----------------------------------------------------------
  ## end test against raw data -------------------------------------------
  
  ## test against spatially subsampled data -----------------------------------
  aucs <- lapply(fits, FUN = function(x, sp_name, test_points_ss) {
    auc_iter <- lapply(x, FUN = function(y, sp_name, test_points_ss) {
      # select test dataset to use for this test
      test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                          size = 1)]]
      # subset to spatially  subsampled test data points that are in the test 
      # fold for this model
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" or y$preds$test_fold == T are in different 
      # CV folds b/c cv block boundary is within the hectad and checklists are 
      # at finer spatial resolution, or because CV used randomly chosen 
      # checklists, regardless of what hectad the were in.  E.g. with random 
      # CV, if one checklist from a hectad is in the test fold and another is 
      # not, the hectad will be in "y$preds$hectad[y$preds$test_fold == T]" and
      # so the training checklist will be included by the previous line of code,
      # but will now be dropped by this following line of code.)
      test_data <- test_data[test_data$test_fold == T, ] 
      n_dets_nondets <- tryCatch({
        data.frame(table(test_data[, colnames(test_data) == gsub(" ", ".", 
                                                                 sp_name)]))}, 
        error = function(x) NA)
      auc <- tryCatch({
        as.numeric(pROC::roc(response = test_data[, colnames(test_data) == 
                                                    gsub(" ", ".", sp_name)], 
                             predictor = test_data$pred)$auc)}, 
        error = function(x) NA)
      list(auc = auc, n_dets_nondets = n_dets_nondets)
    }, sp_name = sp_name, test_points_ss = test_points_ss)
  }, sp_name = sp_to_fit[[i]], 
  test_points_ss = test_points_ss[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  simps_trains_hec <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$simpson_training_hectad}, error = function(x) NA)
    })})
  simps_trains_blk <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$simpson_training_subsampBlock}, error = function(x) NA)
    })})
  
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", metric = "AUC", 
                   value = unlist(lapply(aucs, FUN = function(x) {
                     lapply(x, FUN = function (y) {y$auc})})), 
                   block_cv_range = unlist(block_ranges), 
                   simpson_training_hectad = unlist(simps_trains_hec), 
                   simpson_training_subsampBlock = unlist(simps_trains_blk), 
                   proportion_detections = unlist(prop_dets), 
                   n_dets_in_test = unlist(lapply(aucs, FUN = function(x) {
                     lapply(x, FUN = function (y) {
                       if("1" %in% as.character(y$n_dets_nondets$Var1)) {
                         y$n_dets_nondets$Freq[y$n_dets_nondets$Var1 == "1"]
                       } else 0
                     })
                   })), 
                   n_nondets_in_test = unlist(lapply(aucs, FUN = function(x) {
                     lapply(x, FUN = function (y) {
                       if("0" %in% as.character(y$n_dets_nondets$Var1)) {
                         y$n_dets_nondets$Freq[y$n_dets_nondets$Var1 == "0"]
                       } else 0})
                   })))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges, simps_trains_hec, simps_trains_blk) 
  
  # Kappa ------------------------------------
  # get kappa for each fold
  kp <- lapply(fits, function(x, sp_name, kappas, test_points_ss) {
    kp_iter <- lapply(x, FUN = function(y, sp_name, kappas, test_points_ss) {
      # select test dataset to use for this test
      test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                          size = 1)]]
      # subset to spatially  subsampled test data points that are in the test 
      # fold for this model
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" or y$preds$test_fold == T are actually in 
      # different CV folds.
      test_data <- test_data[test_data$test_fold == T, ]
      tryCatch({
        k_res <- sapply(kappas, kappa_calc, 
                        resp = test_data[, colnames(test_data) == sp_name], 
                        pred = test_data$pred)
        # return kappa and threshold that maximised kappa
        k_res[k_res == max(k_res)][1]}, 
        error = function(x) NA)
    }, sp_name = sp_name, kappas = kappas, test_points_ss = test_points_ss)
  }, sp_name = gsub(" ", ".", sp_to_fit[[i]]), kappas = kappas, 
  test_points_ss = test_points_ss[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", metric = "Kappa", 
                   value = unlist(kp), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end Kappa ----------------------------------
  
  # sensitivity --------------------------------------------------
  sens <- mapply(FUN = function(x, thresh, sp_name, test_points_ss) {
    sens_iter <- mapply(FUN = function(y, thresh, sp_name, test_points_ss) {
      # select test dataset to use for this test
      test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                          size = 1)]]
      # subset to spatially  subsampled test data points that are in the test 
      # fold for this model
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" or y$preds$test_fold == T are actually in 
      # different CV folds.
      test_data <- test_data[test_data$test_fold == T, ]
      tryCatch({
        sensitivity(threshold = thresh, 
                    response = test_data[, colnames(test_data) == 
                                           gsub(" ", ".", sp_name)], 
                    predictions = test_data$pred)}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name, 
                                  test_points_ss = test_points_ss))
  }, fits, kp, 
  MoreArgs = list(sp_name = sp_to_fit[[i]], 
                  test_points_ss = test_points_ss[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", 
                   metric = "sensitivity", 
                   value = unlist(sens), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end sensitivity --------------------------------------
  
  # specificity --------------------------------------------------
  specif <- mapply(FUN = function(x, thresh, sp_name, test_points_ss) {
    specif_iter <- mapply(FUN = function(y, thresh, sp_name, test_points_ss) {
      # select test dataset to use for this test
      test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                          size = 1)]]
      # subset to spatially  subsampled test data points that are in the test 
      # fold for this model
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" or y$preds$test_fold == T are actually in 
      # different CV folds.
      test_data <- test_data[test_data$test_fold == T, ]
      tryCatch({
        specificity(threshold = thresh, 
                    response = test_data[, colnames(test_data) == 
                                           gsub(" ", ".", sp_name)], 
                    predictions = test_data$pred)}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name, 
                                  test_points_ss = test_points_ss))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]], 
                               test_points_ss = test_points_ss[[i]]),
  SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", 
                   metric = "specificity", 
                   value = unlist(specif), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end specificity --------------------------------------
  
  # Brier score --------------------------------------------------------------
  brier <- lapply(fits, FUN = function(x, sp_name, test_points_ss) {
    br_iter <- lapply(x, FUN = function(y, sp_name, test_points_ss) {
      tryCatch({
        # select test dataset to use for this test
        test_data <- test_points_ss[[sample(1:length(test_points_ss), 
                                            size = 1)]]
        # subset to spatially  subsampled test data points that are in the test 
        # fold for this model
        test_data <- test_data[test_data$hectad %in% y$test_sites |
                                 test_data$hectad %in% 
                                 y$preds$hectad[y$preds$test_fold == T], ]
        test_data <- left_join(test_data, y$preds) # join predictions to test data
        # drop checklists not in the test fold (sometimes checlists in hectads 
        # that are in "y$test_sites" or y$preds$test_fold == T are in different 
        # CV folds b/c cv block boundary is within the hectad and checklists are 
        # at finer spatial resolution, or because CV used randomly chosen 
        # checklists, regardless of what hectad the were in. 
        test_data <- test_data[test_data$test_fold == T, ] 
        pr <- test_data$pred # predicted values
        # observed values
        ob <- test_data[, colnames(test_data) == gsub(" ", ".", sp_name)]
        sq_er <- (pr - ob)^2 # squared error
        mean(sq_er) # return mean squared error (Brier score)
      }, 
      error = function(x) NA)
    }, sp_name = sp_name, test_points_ss = test_points_ss)
  }, sp_name = sp_to_fit[[i]], test_points_ss = test_points_ss[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  prop_dets <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$proportion_detections}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", metric = "Brier", 
                   value = unlist(brier), 
                   block_cv_range = unlist(block_ranges), 
                   proportion_detections = unlist(prop_dets))
  evals <- bind_rows(evals, ev)
  # end Brier score -----------------------------------------------------------
  ## end test on spatially subsampled data ------------------------------------
  rm(fits)
}
### end evaluation ---------------------------------------------------

print("Writing out evals file.")

write_csv(evals, paste0("./saved_objects/evals_", mod_name, ".csv"))

rm(test_points_ss, evals)
