##########################
## Functions for maps of ignorance
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 25 Oct 2019
## last modified: 9 Jan 2020
###########################

euc_dist <- function(x, y) {
  # compute the euclidean distance between points x and y 
  # ARGS: x, y - vectors of matching lengths giving the coordinates of each
  #         point x and y in any number of dimensions
  # VALUE:  The euclidean distance between x and y
  if(length(x) != length(y)) stop("x and y must be vectors of the same length.")
  sides <- x - y # get distance between points in each dimension
  sqrt(sum(sides^2)) # return euclidean distance
}

dist_to_nearest_record <- function(df, query_points, coords, parallel = FALSE, 
                                   ncores = NULL, ...) {
  # calculate the euclidean distance to the nearest record in df from the 
  # points in query_points.
  # All dimensions should be scaled before being sent into this function.
  # ARGS: df - data frame of observations containing columns for coordinates
  #       query_points - data frame giving the points from which distances to 
  #           nearest records are desired.  This should have two columns for 
  #           coordinates with the same names as those in df and coords
  #       coords - character vector of any length giving the names of the columns
  #           containing the coordinates of each record
  # VALUE: a data frame with all the points in query_points and the distances 
  #         from those points to the nearest records in df
  list2env(list(...), envir = environment())
  # Check that dimensions have been scaled by making sure SDs are within 2
  # (all SDs should be close to 1 after scaling, so differences between SDs should
  # be small, though SDs might not be exactly 1 if some sites are removed after
  # scaling variables)
  df_sds <- lapply(df[, coords], function(x) sd(x, na.rm = T))
  df_sds <- as.matrix(dist(df_sds))
  if(any(df_sds > 2)) warning("Did you scale spatial dimensions of records before using dist_to_nearest_record()?")
  q_sds <- lapply(query_points[, coords], function(x) sd(x, na.rm = T))
  q_sds <- as.matrix(dist(q_sds))
  if(any(q_sds > 2)) warning("Did you scale dimensions of query points before using dist_to_nearest_record()?.")
  
  dists <- query_points # make df to hold resulting minimum distances
  dists$dist_to_nearest_rec <- NA
  
  # subset columns to only coords and coerce to data frame (in case it is a tibble)
  df <- data.frame(df[, coords]) 
  query_points <- data.frame(query_points[, coords])
  # coerce all columns to numeric
  for (i in 1:ncol(df)) {df[, i] <- as.numeric(df[, i])}
  for (i in 1:ncol(query_points)) {
    query_points[, i] <- as.numeric(query_points[, i])}
  
  if(parallel) {
    cl <- makeCluster(ncores)
    for (i in 1:nrow(query_points)) {
      # find smallest distance between query point i and a record in df
      ds <- parApply(cl = cl, X = df, MARGIN = 1, FUN = euc_dist, 
                     y = query_points[i, ], chunk.size = chunk.size)
      if(any(!is.na(ds))) {
        dists$dist_to_nearest_rec[i] <- min(ds, na.rm = TRUE)
      } else dists$dist_to_nearest_rec[i] <- NA
    }
    stopCluster(cl)
  }
  if(!parallel) {
    for (i in 1:nrow(query_points)) {
      # find smallest distance between query point i and a record in df
      dists$dist_to_nearest_rec[i] <- min(apply(df, MARGIN = 1, FUN = euc_dist, 
                                                y = query_points[i, ]))
    }
  }
  
  # make a relative distance colulmn also so most distant grid query point has distance = 1
  dists$dist_to_nearest_rec_relative <- dists$dist_to_nearest_rec/max(
    dists$dist_to_nearest_rec, na.rm = T)
  dists # return df with distance from each query point to nearest record
}