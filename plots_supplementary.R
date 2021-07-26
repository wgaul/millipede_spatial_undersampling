#############################
## Plots for millipedes SDMs - Supplementary materials
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 10 June 2020
## last modified: 9 Dec 2020
##############################
t_size <- 20

library(GGally)

## predictor variable correlations ------------------------------------------
# make sure predictor variables are in the order I specify in the colnames 
# argument to ggpairs
pred_vals <- data.frame(eastings = mill_wide$eastings, 
                        northings = mill_wide$northings, 
                        day_of_year = mill_wide$day_of_year, 
                        list_length = mill_wide$list_length, 
                        mean_tn = mill_wide$mean_tn, 
                        mean_rr = mill_wide$mean_rr, 
                        elev = mill_wide$elev, 
                        artificial = mill_wide$artificial_surfaces, 
                        arable = mill_wide$arable_l2, 
                        wetlands = mill_wide$wetlands_l1, 
                        forest = mill_wide$forest_seminatural_l1, 
                        pasture = mill_wide$pasture_l2)
pred_cor_plot <- ggpairs(
  data = pred_vals, 
  # title = "Predictor variable values\non millipede checklists", 
  axisLabels = "none", 
  columnLabels = c("eastings", "northings", "day of\nyear", "checklist\nlength",
                   "min.\ntemp.", "precip.", "elevation", 
                   "artificial\nsurfaces", "arable\nland", 
                   "wetlands", "forest\nand\nsemi-\nnatural", "pasture"))
## end predictor variable correlations --------------------------------------
plot_preds <- mask(pred_brick_1k, ir_TM75) # copy predictor variables
# select only predictors I used in models
plot_preds <- subset(plot_preds, c("mean_tn", "mean_rr", "elev", 
                                   "artificial_surfaces", "arable_l2", 
                                   "wetlands_l1", 
                                   "forest_seminatural_l1", "pasture_l2"))
# rename layers for readability
names(plot_preds) <- c("minimum\ntemperature", "precipitation", "elevation", 
                       "artificial\nsurfaces", "arable land", "wetlands", 
                       "forest seminatural", "pasture")
## predictor variable maps --------------------------------------------------

## end predictor variable maps ------------------------------------------------

## checklist length plots -----------------------------------------------------
hist(mill_wide$list_length, breaks = 15)
# calculated median and mean list length in each hectad
ll_df <- group_by(data.frame(mill_wide), hectad) %>%
  summarise(median_ll = median(list_length),
            mean_ll = mean(list_length), 
            n_lists = n()) %>% 
  left_join(., hec_names)

list_length_map <- ggplot(data = ll_df, 
                          aes(x = eastings, y = northings, fill = median_ll)) + 
  geom_tile() + 
  scale_fill_continuous(name = "Median\nchecklist\nlength") + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

# ggplot(data = ll_df, aes(x = eastings, y = northings, fill = n_lists)) +
#   geom_tile() +
#   scale_fill_continuous(name = "Number\nof\nchecklists") +
#   theme_bw()
## end checklist length plots -------------------------------------------------

## spatial evenness of training and test datasets ----------------------------
# spat_evenness_boxplot_bothCV <- ggplot(
#   data = evals[!is.na(evals$simpson_training_subsampBlock), ], 
#   aes(x = factor(block_cv_range), y = simpson_training_subsampBlock)) + 
#   geom_boxplot() + 
#   facet_wrap(~train_data) + 
#   ggtitle("Spatial Evenness of training data\ncalculated at subsample block scale\nincluding cells with zero observations")
# spat_evenness_boxplot_bothCV

spat_evenness_boxplot <- ggplot(
  data = evals[!is.na(evals$simpson_training_subsampBlock) & 
                 evals$block_cv_range == "random" & 
                 evals$test_data == "spat_subsamp" & 
                 evals$metric == "AUC", ], 
  aes(x = factor(train_data, 
                 levels = c("raw", "spat_subsamp"), 
                 labels = c("raw", "spatially\nunder-sampled")), 
      y = simpson_training_subsampBlock)) + 
  geom_boxplot() + 
  facet_wrap(~species) + 
  xlab("Training data") + ylab("Simpson's evenness") +
  # ggtitle("Spatial Evenness of training data\ncalculated at subsample block scale\nincluding cells with zero observations") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
spat_evenness_boxplot
# number of datasets used in boxplot = 
# number of species * number of models * number of model runs (33)
4*33
### end spatial evenness of datasets plot --------------------------------------

## plot results using random CV -----------------------------------------------
# how many detections and non-detections in test folds?
ndet_df <- pivot_longer(evals[evals$test_data == "spat_subsamp", ], 
                        n_dets_in_test:n_nondets_in_test, 
                        names_to = "detection", values_to = "count")
ggplot(data = ndet_df, 
       aes (x = factor(as.character(block_cv_range)), y = count, 
            color = detection)) + 
  geom_boxplot() + 
  ylab("Number in test dataset") + 
  ggtitle("When testing on spatially subsampled data") + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))


# class balance - proportion of observations that are detections in training data
class_balance_boxplot <- ggplot(
  data = evals[!is.na(evals$proportion_detections), ], 
  aes(x = factor(train_data, 
                 levels = c("raw", "spat_subsamp"), 
                 labels = c("raw", "spatially\nunder-sampled")), 
      y = proportion_detections)) + 
  geom_boxplot() + 
  facet_wrap(~species) + 
  ylab("Proportion of checklists with a detection\nin training CV folds") + 
  xlab("Training data") +
  geom_abline(slope = 0, intercept = 0.5, linetype = "dashed") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))




## look at effect of block CV range size
ggplot(data = evals[evals$metric == "AUC" & !is.na(evals$block_cv_range), ], 
       aes(x = factor(block_cv_range), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year\n(DOY) +\nList Length\n", 
                        "\nEnvironment +\nDOY +\nList Length\n", 
                        "\nLat + Lon +\nDOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(train_data, 
                               levels = c("raw", "spat_subsamp"), 
                               labels = c("raw", "spatially\nundersampled"))) + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8)



### plot AUC with all combinations of training and test data for random CV ----
ggplot(data = evals[evals$metric == "AUC" & 
                      as.character(evals$block_cv_range) == "random", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "spat_ll_rf","env_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year\n(DOY) +\nList Length\n", 
                        "\nLat + Lon +\nDOY +\nList Length\n",
                        "\nEnvironment +\nDOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("AUC\n(Cross-Validated)") + 
  ggtitle(paste0("Random CV\nmodel resolution: ", analysis_resolution)) + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))
### end AUC ----------------------------------------------------------------


## plot all performance metrics for best model - block CV ---------------------
evals_median_blockCV <- filter(evals, model == "env_spat_ll_rf" & 
                              block_cv_range == "30000" & 
                              test_data == "spat_subsamp") %>%
  select(species, train_data, metric, value) %>%
  group_by(species, metric, train_data) %>% 
  summarise(median = median(value))


ggplot(data = evals_median_blockCV, 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = median, 
           group = factor(species))) + 
  geom_point(aes(color = factor(species))) + 
  geom_line(aes(color = factor(species))) + 
  facet_wrap(~factor(metric, 
                     levels = c("AUC", "sensitivity", "specificity", "Kappa", 
                                "Brier"), ordered = T)) + 
  xlab("Training Data") + 
  ylab("value") + 
  ggtitle(paste0("30 km block CV\nmodel resolution: ", analysis_resolution)) + 
  scale_color_viridis_d(name = "Species", #option = "magma", 
                        begin = 0, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))
## end plot performance for best model - block CV -----------------------------



### plot variable importance --------------------------------------------------
### plot variable importance only for best model -----------------------
# read in variable importance results
vimp <- list.files("./saved_objects/")
vimp <- vimp[grepl("var_import.*", vimp)]
vimp <- vimp[grepl(paste0(".*", analysis_resolution, ".rds"), vimp)]
vimp <- vimp[grepl(".*env_spat_ll.*", vimp)]
vimp <- lapply(vimp, function(x) readRDS(paste0("./saved_objects/", x)))
# average the variable importance from each CV fold
vimp <- lapply(vimp, FUN = function(x) {
  group_by(x, variable, cv) %>%
    summarise(MeanDecreaseGini = mean(MeanDecreaseGini), 
              species = unique(species), model = unique(model), 
              train_dat = unique(train_dat))
})
vimp <- bind_rows(vimp)

# make environmental variable names pretty
vimp$variable <- gsub("_", " ", vimp$variable)
vimp$variable <- gsub("l2", "", vimp$variable)
vimp$variable <- gsub("l1", "", vimp$variable)
vimp$variable <- gsub("mean rr", "precipitation", vimp$variable)
vimp$variable <- gsub("mean tn", "minimum temperature", vimp$variable)
vimp$variable <- gsub("arable", "arable land", vimp$variable)
vimp$variable <- gsub(" doy", "(day of year)", vimp$variable)

vimp_plots_best <- lapply(sp_to_fit, FUN = function(x, v_df) {
  dat <- v_df[v_df$cv == "random" & v_df$species == x, ]
  dat <- dat[order(dat$MeanDecreaseGini, decreasing = FALSE), ]
  ggplot(data = dat, 
         aes(x = factor(variable, levels = dat$variable, 
                        labels = dat$variable, ordered = T), 
             y = MeanDecreaseGini)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    ylab("Mean decrease Gini index") + xlab("") + 
    ggtitle(dat$species[1]) +
    # ggtitle(paste0(dat$species[1], "\n", dat$model[1], "\n", 
    #                "CV: ", dat$cv[1], 
    #                " analysis resolution: ", analysis_resolution)) + 
    facet_wrap(~factor(train_dat, 
                       levels = c("raw", "spat_subsamp"), 
                       labels = c("raw", "spatially\nunder-\nsampled")), 
               scales = "free") + 
    theme_bw() + 
    theme(text = element_text(size = 0.7*t_size))
}, v_df = vimp)

# for(i in 1:length(vimp_plots_best)) print(vimp_plots_best[i])
multiplot(plotlist = vimp_plots_best, cols = 2) # variable importance for best model


## graph variable importance by CV strategy for best model and both raw and 
## subsampled training data
# CV strategy does not change variable importance conlcusions
# read in variable importance results
vimp <- list.files("./saved_objects/")
vimp <- vimp[grepl("var_import.*", vimp)]
vimp <- vimp[grepl(paste0(".*", analysis_resolution, ".rds"), vimp)]
vimp <- vimp[grepl(".*env_spat_ll.*", vimp)]
vimp <- lapply(vimp, function(x) readRDS(paste0("./saved_objects/", x)))
# average the variable importance from each CV fold
vimp <- lapply(vimp, FUN = function(x) {
  group_by(x, variable, cv) %>%
    summarise(MeanDecreaseGini = mean(MeanDecreaseGini), 
              species = unique(species), model = unique(model), 
              train_dat = unique(train_dat))
})
vimp_plots <- lapply(vimp, FUN = function(x) {
  ggplot(data = x, 
         aes(x = variable, y = MeanDecreaseGini)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    ggtitle(paste0(x$species[1], "\n", x$model[1], ", trained with: ", 
                   x$train_dat[1], "\n")) + 
    facet_wrap(~cv)
})

for(i in 1:length(vimp_plots)) print(vimp_plots[i])
### end plot variable importance ----------------------------------------------


### plot partial dependence ---------------------------------------------------
# plot pd from raw and spatially under-sampled training data
## make pd plots using random CV, spatially under-sampled training data, and
## the most complex model
# read in variable importance results
vimp <- list.files("./saved_objects/")
vimp <- vimp[grepl("var_import.*", vimp)]
vimp <- vimp[grepl(paste0(".*", analysis_resolution, ".rds"), vimp)]
vimp <- vimp[grepl(".*env_spat_ll.*", vimp)]
vimp <- lapply(vimp, function(x) readRDS(paste0("./saved_objects/", x)))
# average the variable importance from each CV fold
vimp <- lapply(vimp, FUN = function(x) {
  group_by(x, variable, cv) %>%
    summarise(MeanDecreaseGini = mean(MeanDecreaseGini), 
              species = unique(species), model = unique(model), 
              train_dat = unique(train_dat))
})
vimp <- bind_rows(vimp)

# read in partial dependence files
pd <- list.files("./saved_objects/")
pd <- pd[grepl("partial_depen.*", pd)]
pd <- pd[grepl(paste0(".*", analysis_resolution, ".rds"), pd)]
pd <- pd[grepl(paste0(".*", mods_for_pd_plots, ".*", collapse = "|"), pd)]
names(pd) <- pd
pd <- lapply(pd, function(x) readRDS(paste0("./saved_objects/", x)))
# average the dependence for each variable from each CV fold
pd <- lapply(pd, FUN = function(x) {
  group_by(x, x, variable, cv) %>%
    summarise(y = mean(y), 
              species = unique(species), model = unique(model), 
              train_data = unique(train_data)) %>%
    filter(variable != "cos_doy" & variable != "sin_doy")
})
pd <- bind_rows(pd)
# make eastings and northings be in km instead of m
pd$x[pd$variable == "eastings" | pd$variable == "northings"] <- 
  pd$x[pd$variable == "eastings" | pd$variable == "northings"] / 1000

# make plots
pd_raw_plots <- lapply(sp_to_fit, FUN = function(x, dat, vimp) {
  vi <- vimp[vimp$species == x, ] # get variable importance for this sp.
  vi <- vi[order(vi$MeanDecreaseGini, decreasing = T), ]
  vi <- vi[-which(grepl(".*doy", vi$variable))[2], ]
  vi$variable <- gsub(".*doy", "day_of_year", vi$variable)
  # make better variable names
  vi$variable <- gsub("arabl.*", "arable\nland", vi$variable)
  vi$variable <- gsub("artific.*", "artificial\nsurfaces", vi$variable)
  vi$variable <- gsub("elev", "elevation", vi$variable)
  vi$variable <- gsub("fores.*", "forest and\nsemi-natural\nland", vi$variable)
  vi$variable <- gsub("list_l.*", "checklist\nlength", vi$variable)
  vi$variable <- gsub("mean_rr", "annual\nprecipitation", vi$variable)
  vi$variable <- gsub("mean_tn", "annual\nminimum\ntemperature", vi$variable)
  vi$variable <- gsub("pastu.*", "pasture", vi$variable)
  vi$variable <- gsub("wetla.*", "wetlands", vi$variable)
  vi$variable <- gsub("day_o.*", "day of\nyear", vi$variable)
  
  # make variable column a factor
  vi$variable <- factor(vi$variable, levels = vi$variable, 
                        ordered = T)
  
  pdat <- dat[dat$species == x, ] # get pd data for this sp.
  # make better variable names
  pdat$variable <- gsub("arabl.*", "arable\nland", pdat$variable)
  pdat$variable <- gsub("artific.*", "artificial\nsurfaces", pdat$variable)
  pdat$variable <- gsub("elev", "elevation", pdat$variable)
  pdat$variable <- gsub("fores.*", "forest and\nsemi-natural\nland", 
                        pdat$variable)
  pdat$variable <- gsub("list_l.*", "checklist\nlength", pdat$variable)
  pdat$variable <- gsub("mean_rr", "annual\nprecipitation", pdat$variable)
  pdat$variable <- gsub("mean_tn", "annual\nminimum\ntemperature", 
                        pdat$variable)
  pdat$variable <- gsub("pastu.*", "pasture", pdat$variable)
  pdat$variable <- gsub("wetla.*", "wetlands", pdat$variable)
  pdat$variable <- gsub("day_o.*", "day of\nyear", pdat$variable)
  # make variable column a factor
  pdat$variable <- factor(pdat$variable, levels = levels(vi$variable))
  

  ggplot(data = pdat,
         aes(x = x, y = y)) +
    geom_point() +
    geom_line() +
    facet_wrap(~variable, scales = "free_x") +
    ylab("Partial dependence") +
    xlab("") +
    ggtitle(pdat$species[1]) + 
    theme_bw() + 
    theme(text = element_text(size = t_size), 
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
  
}, dat = pd[pd$train_data == "raw" & pd$cv == "random", ],
vimp = vimp[vimp$train_dat == "spat_subsamp" & vimp$cv == "random", ])
### end plot partial dependence ---------------------------------------------

### plot predictions with standardized list length ----------------------------
# load standardized predictions
stpred <- list.files("./saved_objects/")
stpred <- stpred[grepl("standard_pre.*", stpred)]
stpred <- stpred[grepl(paste0(".*", analysis_resolution, ".rds"), stpred)]
stpred <- stpred[grepl(paste0(".*", mods_for_pd_plots, ".*", collapse = "|"), 
                       stpred)]
names(stpred) <- gsub("standard.*tions_", "", stpred)
names(stpred) <- gsub(".rds", "", names(stpred))

# # plot average of predictions from all 5 folds (so 4 predictions will be to
# # training data, one prediction to test data in each grid cell)
# prediction_plots <- mapply(FUN = function(x, nm, mill) {
#   dat <- readRDS(paste0("./saved_objects/", x)) # load object
#   # map only predictions from random CV
#   dat <- dat[dat$cv == "random", ] # use only random CV results
#   sp <- gsub(".*Samp_", "", nm) # get species name
#   sp <- gsub("1000.*", "", sp)
#   sp <- gsub("_", " ", sp)
#   mod <- gsub("_SubSamp.*|_noSub.*", "", nm) # get model name
#   train_data <- gsub("^.*_rf_", "", nm) # get training data type name
#   train_data <- gsub("_.*_.*$", "", train_data)
#   
#   # get average predictions for each grid cell (averaging over all folds)
#   dat <- group_by(dat, en) %>%
#     summarise(mean_prediction = mean(mean_pred, na.rm = T), 
#               eastings = mean(eastings), northings = mean(northings))
#   ggplot() + 
#     geom_tile(data = dat, 
#               aes(x = eastings, y = northings, fill = mean_prediction)) + 
#     ggtitle(sp) + 
#     scale_fill_continuous(name = "") + 
#     theme_bw()
# }, stpred, names(stpred), MoreArgs = list(mill = mill), SIMPLIFY = FALSE)
# 
# multiplot(plotlist = prediction_plots[c(5, 2, 6, 11, 8, 12)], 
#           layout = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = T))
# multiplot(plotlist = prediction_plots[c(1, 4, 3, 7, 10, 9)], 
#           layout = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = T))
# for(i in 1:length(prediction_plots)) print(prediction_plots[i])

## using code from main text
sp_map_list <- list()
pred_se_cor_plots <- list()
for(i in 1:length(sp_to_fit)) {
  sn <- sp_to_fit[[i]]
  # load standardized predictions
  sp_preds <- list(
    raw = readRDS(
      paste0("./saved_objects/standard_predictions_env_spat_ll_rf_noSubSamp_", 
             gsub(" ", "_", sn), "1000.rds")), 
    spat_subsamp = readRDS(
      paste0("./saved_objects/standard_predictions_env_spat_ll_rf_SubSamp_", 
             gsub(" ", "_", sn), "1000.rds")))
  
  sp_se_dat <- lapply(sp_preds, function(dat) {
    # map only predictions from random CV
    dat <- dat[dat$cv == "random", ] # use only random CV results
    imod <- "env_spat_ll_rf"
    
    # get average predictions for each grid cell (averaging over all folds)
    dat <- group_by(dat, en) %>%
      summarise(mean_prediction = mean(mean_pred, na.rm = T), 
                se = se(mean_pred), 
                eastings = mean(eastings), northings = mean(northings))
    dat
  })

  # plot correlation of predction and se
  sp_se_dat$raw$training_data <- "raw"
  sp_se_dat$spat_subsamp$training_data <- "spat_subsamp"
  cor_dat <- bind_rows(sp_se_dat)
  pred_se_cor_plots[[length(pred_se_cor_plots) + 1]] <- ggplot(
    data = cor_dat, aes(x = mean_prediction, y = se)) + 
    geom_point() + 
    facet_wrap(~factor(training_data, levels = c("raw", "spat_subsamp"), 
                       labels = c("raw", "spatially\nundersampled"))) + 
    xlab("Mean prediction") + ylab("Standard error") + 
    ggtitle(sn) + 
    theme_bw() + 
    theme(text = element_text(size = 0.55*t_size))
  rm(cor_dat)
  
  maps_prediction <- lapply(
    sp_se_dat, function(dat, annot, ir_TM75) {
      # map mean predictions
      ggplot() + 
        geom_sf(data = st_as_sf(ir_TM75), fill = NA) + 
        geom_tile(data = dat, 
                  aes(x = eastings, y = northings, fill = mean_prediction)) +
        ylab("") + xlab("Longitude") + 
        guides(fill = guide_colorbar(title = "", 
                                     barwidth = unit(0.4 * t_size, 
                                                     "points"))) +
        theme_bw() + 
        theme(text = element_text(size = 0.55*t_size), 
              axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), 
              plot.margin = unit(c(-0.25, -0.4, 0, -0.5), "lines"))
    }, annot = annot, ir_TM75 = ir_TM75)
  
  maps_ciWidth <- lapply(sp_se_dat, function(dat, annot, ir_TM75) {
    # map se of mean prediction
    ggplot() +
      geom_sf(data = st_as_sf(ir_TM75), fill = NA) +
      geom_tile(data = dat,
                aes(x = eastings, y = northings, fill = se)) +
      ylab("") + xlab("Longitude") +
      scale_fill_gradient(low = "white", high = "red") +  
      guides(fill = guide_colorbar(title = "",
                                   barwidth = unit(0.4 * t_size, "points"))) +
      theme_bw() +
      theme(text = element_text(size = 0.55*t_size),
            axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
            plot.margin = unit(c(0, -0.6, -0.25, -0.5), "lines"))
  }, annot = annot, ir_TM75 = ir_TM75)
  
  ## add raw observations for this sp
  # get species observations from mill_wide
  sobs <- data.frame(mill_wide)[, names(mill_wide) == sn]
  
  map_observed <- ggplot() + 
    geom_sf(data = st_as_sf(ir_TM75), fill = NA) + 
    geom_sf(
      data = st_as_sf(mill_wide[sobs == 0, ]), 
      color = "dark grey", size = 0.02*t_size) + 
    geom_sf(data = st_as_sf(mill_wide[sobs > 0, ]), 
            color = "dark orange", size = 0.02*t_size) + 
    geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) + 
    geom_text(data = annot[c(2, 4), ], aes(x = x1, y = y1, label = label), 
              size = 0.12*t_size) + 
    geom_segment(data = annot[3, ], aes(x = x1, xend = x2, y = y1, yend = y2), 
                 arrow = arrow(length = unit(0.1, "npc"))) + 
    ylab("Latitude") + xlab("Longitude") + 
    ggtitle("(a)", subtitle = gsub(" ", "\n", sn)) +
    theme_bw() + 
    theme(text = element_text(size = 0.55*t_size), 
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), 
          plot.margin = unit(c(-0.25, -0.4, -0.5, -0.5), "lines"))
  
  # combine maps in one ordered list
  sp_maps <- list(observed = map_observed, 
                  raw_pred = maps_prediction$raw, 
                  spat_subsamp_pred = maps_prediction$spat_subsamp, 
                  raw_ciWidth = maps_ciWidth$raw, 
                  spat_subsamp_ciWidth = maps_ciWidth$spat_subsamp)
  # add titles
  sp_maps[[1]] <- sp_maps[[1]] + ggtitle("(a)") + ylab("Latitude")
  sp_maps[[2]] <- sp_maps[[2]] + ggtitle("(b)")
  sp_maps[[3]] <- sp_maps[[3]] + ggtitle("(c)")
  sp_maps[[4]] <- sp_maps[[4]] + ggtitle("(d)") + 
    ylab("Latitude") + xlab("Longitude")
  sp_maps[[5]] <- sp_maps[[5]] + ggtitle("(e)") + xlab("Longitude")
  
  # put maps for this species into sp_map_list
  sp_map_list[[length(sp_map_list) + 1]] <- sp_maps
}
### end plot standardized predictions -----------------------------------------


### plot DOY demonstration to show transformed variables -----------------------
smod <- readRDS("./saved_objects/env_spat_ll_rf_SubSamp_fits_Macrosternodesmus_palicola1000.rds")
doy_dat <- data.frame(newdata)
# add doy 1 to the data, and add more days 
doy1 <- doy_dat[doy_dat$day_of_year == 20, ]
doy1$day_of_year <- 1 
doy1$cos_doy <- cos((2*pi*doy1$day_of_year) / 365)
doy1$sin_doy <- sin((2*pi*doy1$day_of_year) / 365)
doy_copy <- doy_dat
doy_copy$day_of_year <- doy_copy$day_of_year - 10
doy_copy$cos_doy <- cos((2*pi*doy_copy$day_of_year) / 365)
doy_copy$sin_doy <- sin((2*pi*doy_copy$day_of_year) / 365)

doy_dat <- bind_rows(doy_dat, doy1, doy_copy) # add additional days
# get standardized predictions to each grid cell on each day
doy_dat$prediction <- predict(smod[[1]][[1]]$m, newdata=doy_dat, 
                              type = "prob")[, "1"]

doy_dat <- group_by(doy_dat, day_of_year, sin_doy, cos_doy) %>%
  summarise(mean_pred = mean(prediction))

doy_plot <- ggplot(data = doy_dat, 
                   aes(x = sin_doy, y = cos_doy, color = day_of_year, 
                       size = mean_pred)) + 
  geom_point() + 
  xlab("Dsin") + ylab("Dcos") +
  scale_color_continuous(name = "Day of\nyear") + 
  scale_size_continuous("Mean predicted\nprobability\nof being\nrecorded") + 
  theme_bw() + 
  theme(text = element_text(size = t_size))
doy_plot
rm(doy1, smod)
### end plot DOY demonstration ------------------------------------------------

### re-make atlas map for B. tenuis -------------------------------------------
bten_data <- mill_wide[mill_wide$`Boreoiulus tenuis` > 0, ]

# make hec_names spatial 
hec_names_spat <- SpatialPointsDataFrame(
  coords = hec_names[, c("eastings", "northings")], 
  data = hec_names, proj4string = CRS("+init=epsg:29903"))
# make sure millipede data is in same projection as predictor data
hec_names_spat <- spTransform(hec_names_spat, 
                              raster::projection(pred_brick))

## indicate B. tenuis presence record
hec_names_spat$B_tenuis <- NA
hec_names_spat$B_tenuis[hec_names_spat$hectad %in% mill_wide$hectad] <- 0
hec_names_spat$B_tenuis[hec_names_spat$hectad %in% bten_data$hectad] <- 1

Bten_atlas_map <- ggplot() + 
  geom_sf(data = st_as_sf(ir_TM75), fill = NA) + 
  geom_sf(data = st_as_sf(hec_names_spat[!is.na(hec_names_spat$B_tenuis) & 
                                           hec_names_spat$B_tenuis == 0, ]), 
          color = "grey", size = 0.04*t_size) +
  geom_sf(data = st_as_sf(hec_names_spat[!is.na(hec_names_spat$B_tenuis) & 
                                           hec_names_spat$B_tenuis == 1, ]), 
          color = "black", size = 0.04*t_size) +
  geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) + 
  geom_text(data = annot[c(2, 4), ], aes(x = x1, y = y1, label = label)) + 
  geom_segment(data = annot[3, ], aes(x = x1, xend = x2, y = y1, yend = y2), 
               arrow = arrow(length = unit(0.1, "npc"))) + 
  ylab("Latitude") + xlab("Longitude") + 
  ggtitle("(a)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
Bten_atlas_map
### end re-make atlas map ----------------------------------------------------


#### save plots to files -------------------------------------------------------
ggsave("FigS1.jpg", class_balance_boxplot + 
         theme(text = element_text(size = t_size*0.6)), 
       width = 16, height = 16, units = "cm", device = "jpg")
ggsave("FigS2.jpg", spat_evenness_boxplot + 
         theme(text = element_text(size = t_size*0.6)), 
       width = 16, height = 16, units = "cm", device = "jpg")
ggsave("FigS3.jpg", plot(plot_preds), 
       width = 15, height = 15, units = "cm", device = "jpg")
ggsave("FigS4.jpg", doy_plot + 
         theme(text = element_text(size = t_size)), 
       width = 15, height = 15, units = "cm", device = "jpg")
ggsave("FigS5.jpg", multiplot(plotlist = vimp_plots_best[1:3], cols = 1), 
       width = 15, height = 20, units = "cm", device = "jpg")
ggsave("FigS6.jpg", multiplot(plotlist = vimp_plots_best[4:6], cols = 1), 
       width = 15, height = 22, units = "cm", device = "jpg")


ggsave("FigS18.jpg", pred_se_cor_plots[[1]] + pred_se_cor_plots[[2]] + 
         pred_se_cor_plots[[3]] + pred_se_cor_plots[[4]] + 
         pred_se_cor_plots[[5]] + pred_se_cor_plots[[6]] + 
         plot_layout(ncol = 2), 
       width = 16, height = 16, units = "cm", device = "jpg")
ggsave("FigS19.jpg", pred_cor_plot, width = 20, height = 20, 
       units = "cm", device = "jpg")
ggsave("FigS20.jpg", list_length_map, width = 17, height = 15, 
       units = "cm", device = "jpg")




ggsave("FigSPD1.jpg", pd_raw_plots[["Macrosternodesmus palicola"]] + 
       pd_plots[["Macrosternodesmus palicola"]] + 
         plot_annotation(tag_levels = "a") & 
         theme(text = element_text(size = 0.7*t_size)), 
       width = 20, height = 14, units = "cm", device = "jpg")
ggsave("FigSPD2.jpg", pd_raw_plots[["Boreoiulus tenuis"]] + 
         pd_plots[["Boreoiulus tenuis"]] + 
         plot_annotation(tag_levels = "a") & 
         theme(text = element_text(size = 0.7*t_size)), 
       width = 20, height = 14, units = "cm", device = "jpg")
ggsave("FigSPD3.jpg", pd_raw_plots[["Ommatoiulus sabulosus"]] + 
         pd_plots[["Ommatoiulus sabulosus"]] + 
         plot_annotation(tag_levels = "a") & 
         theme(text = element_text(size = 0.7*t_size)), 
       width = 20, height = 14, units = "cm", device = "jpg")
ggsave("FigSPD4.jpg", pd_raw_plots[["Blaniulus guttulatus"]] + 
         pd_plots[["Blaniulus guttulatus"]] + 
         plot_annotation(tag_levels = "a") & 
         theme(text = element_text(size = 0.6*t_size), 
               axis.text.x = element_text(angle = 45, 
                                          hjust = 1, vjust = 1)), 
       width = 20, height = 18, units = "cm", device = "jpg")
ggsave("FigSPD5.jpg", pd_raw_plots[["Glomeris marginata"]] + 
         pd_plots[["Glomeris marginata"]] + 
         plot_annotation(tag_levels = "a") & 
         theme(text = element_text(size = 0.6*t_size), 
               axis.text.x = element_text(angle = 45, 
                                          hjust = 1, vjust = 1)), 
       width = 20, height = 18, units = "cm", device = "jpg")
ggsave("FigSPD6.jpg", pd_raw_plots[["Cylindroiulus punctatus"]] + 
         pd_plots[["Cylindroiulus punctatus"]] + 
         plot_annotation(tag_levels = "a") & 
         theme(text = element_text(size = 0.6*t_size), 
               axis.text.x = element_text(angle = 45, 
                                          hjust = 1, vjust = 1)), 
       width = 20, height = 18, units = "cm", device = "jpg")

ggsave("FigSP1.jpg", sp_map_list[[1]][[1]] + sp_map_list[[1]][[2]] + 
         sp_map_list[[1]][[3]] + 
         plot_spacer() + sp_map_list[[1]][[4]] + sp_map_list[[1]][[5]], 
       width = 20, height = 20*0.6, units = "cm", device = "jpg")
ggsave("FigSP2.jpg", sp_map_list[[2]][[1]] + sp_map_list[[2]][[2]] + 
         sp_map_list[[2]][[3]] + 
         plot_spacer() + sp_map_list[[2]][[4]] + sp_map_list[[2]][[5]], 
       width = 20, height = 20*0.6, units = "cm", device = "jpg")
ggsave("FigSP3.jpg", sp_map_list[[3]][[1]] + sp_map_list[[3]][[2]] + 
         sp_map_list[[3]][[3]] + 
         plot_spacer() + sp_map_list[[3]][[4]] + sp_map_list[[3]][[5]], 
       width = 20, height = 20*0.6, units = "cm", device = "jpg")
ggsave("FigSP4.jpg", sp_map_list[[4]][[1]] + sp_map_list[[4]][[2]] + 
         sp_map_list[[4]][[3]] + 
         plot_spacer() + sp_map_list[[4]][[4]] + sp_map_list[[4]][[5]], 
       width = 20, height = 20*0.6, units = "cm", device = "jpg")
ggsave("FigSP5.jpg", sp_map_list[[5]][[1]] + sp_map_list[[5]][[2]] + 
         sp_map_list[[5]][[3]] + 
         plot_spacer() + sp_map_list[[5]][[4]] + sp_map_list[[5]][[5]], 
       width = 20, height = 20*0.6, units = "cm", device = "jpg")
ggsave("FigSP6.jpg", sp_map_list[[6]][[1]] + sp_map_list[[6]][[2]] + 
         sp_map_list[[6]][[3]] + 
         plot_spacer() + sp_map_list[[6]][[4]] + sp_map_list[[6]][[5]], 
       width = 20, height = 20*0.6, units = "cm", device = "jpg")

ggsave("Fig_Btenuis_atlasSDM.jpg", Bten_atlas_map +
         (sp_map_list[[2]][[3]] + ggtitle("(b)")) + plot_spacer() + 
         (sp_map_list[[2]][[5]] + ggtitle("(c)")), 
       width = 14, height = 14, units = "cm", device = "jpg")

# ggsave("FigSP7.jpg", pred_se_cor_plots[[1]], 
#        width = 15, height = 10, units = "cm", device = "jpg")
# ggsave("FigSP8.jpg", pred_se_cor_plots[[2]], 
#        width = 15, height = 10, units = "cm", device = "jpg")
# ggsave("FigSP9.jpg", pred_se_cor_plots[[3]], 
#        width = 15, height = 10, units = "cm", device = "jpg")
# ggsave("FigSP10.jpg", pred_se_cor_plots[[4]], 
#        width = 15, height = 10, units = "cm", device = "jpg")
# ggsave("FigSP11.jpg", pred_se_cor_plots[[5]], 
#        width = 15, height = 10, units = "cm", device = "jpg")
# ggsave("FigSP12.jpg", pred_se_cor_plots[[6]], 
#        width = 15, height = 10, units = "cm", device = "jpg")

