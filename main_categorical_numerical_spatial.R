# Packages
################################################################################
library(spatstat)
library(spatstat.data)
library(sf)
library(dplyr)
library(MASS)
library(ddalpha)


# install.packages("devtools")
# devtools::install_github(repo = "hannahblo/ddandrda", ref = "hb_04.11")
library(ddandrda)



# small running example
################################################################################
# Step 1: get smaller subset of the point pattern
gor <- rescale(gorillas, 1000, unitname = "km")
gex <- lapply(gorillas.extra, rescale, s = 1000, unitname = "km")

# construct artificial example for demonstration
# consider the ppp in a small quadrat
quadrats_np <- quadrats(gex$vegetation, 10)
ass_quadrats <- tileindex(gor$x, gor$y, quadrats_np)
gor_reduced <- gor[c(which(ass_quadrats == "33"),
                     which(ass_quadrats == "34"),
                     which(ass_quadrats == "43"),
                     which(ass_quadrats == "44"))]

new_xrange <- quadrats_np$image$xcol[
  apply(quadrats_np$image$v, 2, function(r) any(r %in% c(31, 32, 41, 42)))]
new_yrange <- quadrats_np$image$yrow[
  apply(quadrats_np$image$v, 1, function(r) any(r %in% c(31, 32, 41, 42)))]
new_owin <- as.owin(list(xrange = c(581.5,582.7),
                         yrange = c(676.3, 677.5)))

vegetation_np <- tess(image = gex$vegetation)
ass_vegetation_reduced <- tileindex(gor_reduced$x, gor_reduced$y, vegetation_np)
gor_transition <- gor_reduced[which(ass_vegetation_reduced == "Transition")]
gor_disturbed <- gor_reduced[which(ass_vegetation_reduced == "Disturbed")]
gor_grassland <- gor_reduced[which(ass_vegetation_reduced == "Grassland")]
gor_primary <- gor_reduced[which(ass_vegetation_reduced == "Primary")]

set.seed(325)
random_primary <- sample(seq(1, 85), 20)
ex_ppp <- gor_primary[random_primary[c(1,2,3,4,5,6,8,10,11,12,14,15, 17)]]
# ex_ppp <- superimpose(ex_ppp, gor_grassland[4])
# ex_ppp <- superimpose(ex_ppp, gor_disturbed[3])
ex_ppp <- unmark(ex_ppp)

# check duplications run
ex_ppp <- superimpose(ex_ppp, c(x = ex_ppp$x[1], y = ex_ppp$y[1]))

# plot(unmark(gor_primary[random_primary]), pch = 1, ces = 0.5)
# points(gor_grassland[4], pch = 19, cex = 1, col = "darkblue")
# points(gor_disturbed[3], pch = 17, cex = 1, col = "darkgreen")
# plot(ex_ppp)
# plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")),
#      las = 2)
# points(ex_ppp, pch = 16)


plot(gex$vegetation, clipwin = new_owin, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")),
     las = 2,
     main = NULL)
points(gor_primary[random_primary[c(1,2,3,4,5,6,8,10,11,12,14,15, 17)]],
       pch = 16)
points(gor_grassland[4],
       pch = 15)
points(gor_disturbed[3],
       pch = 17)

plot(ex_ppp, clipwin = new_owin, main = NULL)

plot(unmark(gor_primary[random_primary[c(1,2,3,4,5,6,8,10,11,12,14,15, 17)]]),
     clipwin = new_owin, pch = 16,
     main = NULL)
points(gor_grassland[4],
       pch = 15)
points(gor_disturbed[3],
       pch = 17)

# step 2 data preparation
plot(gex$elevation)
mm <- transmat(gex$elevation, from = "spatstat", to = "European")
grid_numeric <- c(t(apply(mm, 2, rev)))

plot(gex$vegetation)
mm <- transmat(gex$vegetation, from = "spatstat", to = "European")
grid_nominal <- c(t(apply(mm, 2, rev)))

grid_spatial <- expand.grid(gex$vegetation$xcol, gex$vegetation$yrow)
observed <- as.data.frame(matrix(c(ex_ppp$x, ex_ppp$y), ncol = 2))
observed_count <- observed %>% group_by_all() %>% count
observed <- unname(as.matrix(observed_count[, c(1,2)]))
empirical_prob <- observed_count[, 3]

grid_observed_x <- unlist(lapply(observed[, 1], FUN = function(x) {which.min(abs(gex$vegetation$xcol - x))}))
grid_observed_y <- unlist(lapply(observed[, 2], FUN = function(y) {which.min(abs(gex$vegetation$yrow - y))}))
index_grid_observed <- (grid_observed_y - 1) * length(gex$vegetation$xcol) + grid_observed_x



# Step 3: compute ufg based on reduced ppp, full grid all (spatial, ordinal, nominal) component
start_time <- Sys.time()
erg_2 <- ddandrda::compute_ufg_grid_cns(observed = observed,
                          empirical_prob = as.vector(empirical_prob),
                          observed_in_grid = index_grid_observed,
                          grid_spatial = grid_spatial,
                          grid_numeric = grid_numeric,
                          grid_nominal = grid_nominal,
                          lower_bound_ufg = 2,
                          upper_bound_ufg = 4)
Sys.time() - start_time # Time difference of 32.77225 secs
# saveRDS(erg_2, file = "erg_2_redppp_fullgrid_34.RDS")

gex$depth <- gex$elevation
# get matrix with inpu
matrix_erg <- transmat(matrix(erg_2$depth, nrow = length(gex$vegetation$xcol)),
                       from = "European", to = "Cartesian")
matrix_erg <- transmat(matrix_erg, from = "European", to = "spatstat")
gex$depth$v <- matrix_erg
plot(gex$depth)
points(ex_ppp)

par(mfrow = c(1,3))
plot(gex$elevation)
points(ex_ppp)
plot(gex$depth, clipwin = new_owin, main = NULL, las = 2)
# points(ex_ppp)
plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")))
points(ex_ppp)


# Step 4: compute ufg based on reduced ppp, full grid with spatial and nominal component
# but without elevation/numeric component
grid_numeric <- rep(1, length(grid_numeric))
# reduced version but full grid
start_time <- Sys.time()
erg_3 <- ddandrda::compute_ufg_grid_cns(observed,
                          empirical_prob = empirical_prob,
                          observed_in_grid = index_grid_observed,
                          grid_spatial = grid_spatial,
                          grid_numeric = grid_numeric,
                          grid_nominal = grid_nominal,
                          lower_bound_ufg = 3,
                          upper_bound_ufg = 4)
Sys.time() - start_time
# saveRDS(erg_3, file = "erg_3_redppp_nonumeric_fullgrid_34")

gex$depth <- gex$elevation
# get matrix with input
matrix_erg <- transmat(matrix(erg_3$depth, nrow = length(gex$vegetation$xcol)),
                       from = "European", to = "Cartesian")
matrix_erg <- transmat(matrix_erg, from = "European", to = "spatstat")
gex$depth$v <- matrix_erg
plot(gex$depth, clipwin = new_owin, main = NULL, las = 2)
points(ex_ppp)

par(mfrow = c(1,2))
plot(gex$depth, clipwin = new_owin, main = NULL, las = 2)
points(ex_ppp)
plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")))
points(ex_ppp)



# Step 5: compute ufg based on reduced ppp, full grid, but only spatial component
grid_numeric <- rep(1, length(grid_numeric))
grid_nominal <- rep("a", length(grid_nominal))
# reduced version but full grid
start_time <- Sys.time()
erg_4 <- ddandrda::compute_ufg_grid_cns(observed,
                          empirical_prob = empirical_prob,
                          observed_in_grid = index_grid_observed,
                          grid_spatial = grid_spatial,
                          grid_numeric = grid_numeric,
                          grid_nominal = grid_nominal,
                          lower_bound_ufg = 3,
                          upper_bound_ufg = 3)
Sys.time() - start_time
# saveRDS(erg_4, file = "erg_4_redppp_nonumeric_nocateg_fullgrid_34.RDS")

gex$depth <- gex$elevation
# get matrix with input
matrix_erg <- transmat(matrix(erg_4$depth, nrow = length(gex$vegetation$xcol)),
                       from = "European", to = "Cartesian")
matrix_erg <- transmat(matrix_erg, from = "European", to = "spatstat")
gex$depth$v <- matrix_erg

par(mfrow = c(1,1))
plot(gex$depth)
points(ex_ppp)

# compare to ddalpha
simpl_d <- depth.simplicial(grid_spatial, observed, exact = TRUE)

gex$depth_sd <- gex$elevation
matrix_erg <- transmat(matrix(simpl_d, nrow = length(gex$vegetation$xcol)),
                       from = "European", to = "Cartesian")
matrix_erg <- transmat(matrix_erg, from = "European", to = "spatstat")
gex$depth_sd$v <- matrix_erg

par(mfrow = c(1,2))
plot(gex$depth)
plot(gex$depth_sd)


# entire gorillas example year 2006
################################################################################

# data subset
gor <- rescale(gorillas, 1000, unitname = "km")
gor <- subset(gor, date < "2007-01-01",  drop = TRUE)
gor <- unmark(gor) # two duplicates
gex <- lapply(gorillas.extra, rescale, s = 1000, unitname = "km")

# data preparation
plot(gex$elevation)
mm <- transmat(gex$elevation, from = "spatstat", to = "European")
grid_numeric <- c(t(apply(mm, 2, rev)))

plot(gex$vegetation)
mm <- transmat(gex$vegetation, from = "spatstat", to = "European")
grid_nominal <- c(t(apply(mm, 2, rev)))

grid_spatial <- expand.grid(gex$vegetation$xcol, gex$vegetation$yrow)
observed <- as.data.frame(matrix(c(gor$x, gor$y), ncol = 2))
observed_count <- observed %>% group_by_all() %>% count
observed <- unname(as.matrix(observed_count[, c(1,2)]))
empirical_prob <- as.vector(observed_count[, 3])
empirical_prob <- pull(empirical_prob, n)

grid_observed_x <- unlist(lapply(observed[, 1], FUN = function(x) {which.min(abs(gex$vegetation$xcol - x))}))
grid_observed_y <- unlist(lapply(observed[, 2], FUN = function(y) {which.min(abs(gex$vegetation$yrow - y))}))
index_grid_observed <- (grid_observed_y - 1) * length(gex$vegetation$xcol) + grid_observed_x


# data information
plot(gor, pch = 16, ces = 0.5, main = NULL)

vegetation_np <- tess(image = gex$vegetation)
quadratcount_vegetation <- quadratcount(gor,tess = vegetation_np)

summary(grid_numeric[index_grid_observed])

# computation
start_time <- Sys.time()
erg_4 <- ddandrda::compute_ufg_grid_cns(observed,
                          empirical_prob = empirical_prob,
                          observed_in_grid = index_grid_observed,
                          grid_spatial = grid_spatial,
                          grid_numeric = grid_numeric,
                          grid_nominal = grid_nominal,
                          lower_bound_ufg = 2,
                          upper_bound_ufg = 4)
Sys.time() - start_time
saveRDS(erg_4, file = "erg_4_121ppp_fullgrid_234.RDS")
# erg_4 <- readRDS("erg_4_121ppp_fullgrid_234.RDS")


# Summaries
summary(erg_4$depth)
summary(erg_4$depth[index_grid_observed])
max(erg_4$depth)
max(erg_4$depth[index_grid_observed])



# Plotting
# The depth function from different perspectives

# solely the depth function
gex$depth <- gex$elevation
matrix_erg <- transmat(matrix(erg_4$depth, nrow = length(gex$vegetation$xcol)),
                       from = "European", to = "Cartesian")
matrix_erg <- transmat(matrix_erg, from = "European", to = "spatstat")
gex$depth$v <- matrix_erg

index_NA <- which(is.na(gex$elevation$v))
gex$depth$v[index_NA] <- NA
par(mfrow = c(1,1))
par(bg = 'white')
plot(gex$depth,  main = NULL, las = 2)
# points(gor)


# the depth function in comparison to elevation and vegetation component
par(mfrow = c(1,3))
plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")),
     las = 2)
points(gor)
plot(gex$depth,main = NULL, las = 2)
points(gor)
plot(gex$elevation)
points(gor)




# further detail on data set --> help to analyse the result
# Analysis of Vegetation
par(mfrow = c(1,1))
quadratcount(gor, tess = vegetation_np)
plot(quadratcount(gor, tess = vegetation_np), las = 2)
vegetation_assign <- tileindex(gor$x, gor$y, vegetation_np)
vegetation_assign
plot(vegetation_assign)


gor_transition <- gor[which(vegetation_assign == "Transition")]
gor_secondary <- gor[which(vegetation_assign == "Secondary")]
gor_primary <- gor[which(vegetation_assign == "Primary")]
gor_grassland <- gor[which(vegetation_assign == "Grassland")]
gor_colonising <- gor[which(vegetation_assign == "Colonising")]
gor_disturbed <- gor[which(vegetation_assign == "Disturbed")]

# Plotting the ppp splitted by the vegetation
par(mfrow = c(2,3))
plot(gor_transition, cols = "black", main = "Transition")
plot(gor_secondary, cols = "black", main = "Secondary")
plot(gor_primary, cols = "black", main = "Primary")
plot(gor_grassland, cols = "black", main = "Grassland")
plot(gor_colonising, cols = "black", main = "Colonising")
plot(gor_disturbed, cols = "black", main = "Disturbed")

# Plotting vegetation with the ppp splitted by the vegetation
par(mfrow = c(2,3))
plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")),
     ribbon = FALSE,
     main = "Transition")
points(gor_transition, col = "black")
plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")),
     ribbon = FALSE,
     main = "Secondary")
points(gor_secondary, col = "black")
plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")),
     ribbon = FALSE,
     main = "Primary")
points(gor_primary, col = "black")
plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")),
     ribbon = FALSE,
     main = "Grassland")
points(gor_grassland, col = "black")
plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")),
     ribbon = FALSE,
     main = "Colonising")
points(gor_colonising, col = "black")
plot(gex$vegetation, col = rev(c("gold", "darkgoldenrod1", "lightpink", "lightcoral", "mediumorchid1", "lightslateblue")),
     ribbon = FALSE,
     main = "Disturbed")
points(gor_disturbed, col = "black")


# Analysis of Elevation
par(mfrow = c(1,1))
gex$tess_elevation <- gex$elevation


divison_elevaion <- function(numeric_observed,
                             probs_quantile = c(0, 0.25, 0.5, 0.75, 1),
                             gex) {
  quantile_elev_obs <- quantile(numeric_observed, probs = probs_quantile)

  gex$tess_elevation <- gex$elevation

  above_max <- (gex$elevation$v >= max(grid_numeric[index_grid_observed]))
  below_min <- (gex$elevation$v <= min(grid_numeric[index_grid_observed]))
  gex$tess_elevation$v[above_max] <- "Above Max"
  gex$tess_elevation$v[below_min] <- "Below Min"

  for (i in seq(1, length(probs_quantile) - 1)) {
    between_i_i1 <-  (gex$elevation$v >= quantile_elev_obs[i]) &
      (gex$elevation$v < quantile_elev_obs[i + 1])
    gex$tess_elevation$v[between_i_i1] <- paste0("between ", quantile_elev_obs[i],
                                                 " and ", quantile_elev_obs[i + 1])
  }

  return(gex)
}




max(grid_numeric[index_grid_observed])
max(grid_numeric, na.rm = TRUE)
length(is.na(gex$elevation$v))
above_max <- which(gex$elevation$v > max(grid_numeric[index_grid_observed]))
below_min <- which(gex$elevation$v < min(grid_numeric[index_grid_observed]))
between <- setdiff(seq(1, length(is.na(gex$elevation$v))), c(above_max, below_min))
gex$tess_elevation$v[above_max] <- "Above Max"
gex$tess_elevation$v[below_min] <- "Below Min"
gex$tess_elevation$v[between] <- "Between Min and Max"
gex$tess_elevation$v[is.na(gex$elevation$v)] <- NA
# gex$tess_elevation$v[] <- "Above observed maximum"

tess(image = gex$tess_elevation)
plot(tess(image = gex$tess_elevation),
     las = 2,
     main = "Areas with elevation below, above or inbetween the observed elevation values")



test <- divison_elevaion(numeric_observed = grid_numeric[index_grid_observed],
                         probs_quantile = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                         gex = gex)

quantile_elev_obs <- quantile(grid_numeric[index_grid_observed], probs = c(0, 0.25, 0.5, 0.75, 1))

above_max <- (gex$elevation$v >= max(grid_numeric[index_grid_observed]))
below_min <- (gex$elevation$v <= min(grid_numeric[index_grid_observed]))
between_min_1st <-  (gex$elevation$v >= min(grid_numeric[index_grid_observed])) &
  (gex$elevation$v < quantile_elev_obs[2])
between_1st_med <-  (gex$elevation$v >= quantile_elev_obs[2]) &
  (gex$elevation$v < quantile_elev_obs[3])
between_med_3rd <-  (gex$elevation$v >= quantile_elev_obs[3]) &
  (gex$elevation$v < quantile_elev_obs[4])
between_3rd_max <-  (gex$elevation$v >= quantile_elev_obs[4]) &
  (gex$elevation$v < quantile_elev_obs[5])
gex$tess_elevation$v[below_min] <- 1#  "Below Min"
gex$tess_elevation$v[between_min_1st] <- 2# "Between Min and 1st Qu."
gex$tess_elevation$v[between_1st_med] <- 3 # "Between 1st Qu. and Median"
gex$tess_elevation$v[between_med_3rd] <- 4# "Between Medain and 3rd Qu. "
gex$tess_elevation$v[between_3rd_max] <- 5# "Between 3rd Qu. and Max"
gex$tess_elevation$v[above_max] <- 6# "Above Max"

tess(image = test$tess_elevation)
plot(tess(image = gex$tess_elevation),
     ribbon = FALSE,
     las = 2,
     main = "Tesslation based on discretised elevation", col = rainbow(6, start = 0, end = 0.8, alpha = 1))
legend("right", title = NULL,
       c( "Below Min", "Between Min and 1st Qu.", "Between 1st Qu. and Median",
          "Between Medain and 3rd Qu. ",  "Between 3rd Qu. and Max", "Above Max"),
       fill = rainbow(6, start = 0, end = 0.8), horiz = FALSE)



# Analysis of number of points within fine grid
# Defining regular grid and counting number of observations within each pixel
par(mfrow = c(1,1))
grid_tess <- quadrats(gor, nx = 30, ny = 30)
hist((as.numeric(quadratcount(gor, tess = grid_tess))))
summary( as.factor(as.numeric(quadratcount(gor, tess = grid_tess))))

plot(quadratcount(gor, tess = grid_tess), main = "Regular grid: Counting number of observation within each pixel")
plot(intensity(quadratcount(gor, tess = grid_tess), image = TRUE),
     main = "Regular grid: Counting number of observation within each pixel",
     ribbon = FALSE,
     col = rainbow(5, start = 0.5, end = 1, alpha = 1))
legend("right", title = "Count",
       c( "0", "1.", "2", "3", "4"),
       fill = rainbow(5, start = 0.5, end = 1, alpha = 1), horiz = FALSE)

plot(tess(image = grid_tess))

points(gor, col = "black")

plot(grid_tess)

grid_tess$tiles <- quadratcount(gor, tess = grid_tess)

gex$tess_observed <- gex$elevation
gex$tess_observed$v <- quadratcount(gor, tess = grid_tess)
plot(gex$tess_observed)


# Comparing to projected centralities
maxufg_index <- which(erg_4$depth == max(erg_4$depth))
grid_spatial[maxufg_index, ] # 582.7589 677.1507
grid_numeric[maxufg_index] # 1805
grid_nominal[maxufg_index] # primary


summary(as.factor(grid_nominal[index_grid_observed]))
# 1  2  3  4  5  6
# 24  1  7 76  5  6
# --> compare to above boxplot to see what each number means (4 is primary)

summary(grid_numeric[index_grid_observed])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1340    1637    1792    1781    1954    2053


simpldd <- ddalpha::depth.simplicial(grid_spatial, grid_spatial[index_grid_observed, ])
maxsimpldd_index <- which(simpldd == max(simpldd))
grid_spatial[maxsimpldd_index, ] # 582.6054 677.0585
grid_numeric[maxsimpldd_index] # 1900
grid_nominal[maxsimpldd_index] # primary


# concerning the difference in the numeric component of ufg depth and simplicial depth
length(which(grid_numeric[index_grid_observed] >= grid_numeric[maxsimpldd_index]))/length(index_grid_observed) # 0.3445378
length(which(grid_numeric[index_grid_observed] < grid_numeric[maxsimpldd_index]))/length(index_grid_observed) # 0.6554622

length(which(grid_numeric[index_grid_observed] >= grid_numeric[maxufg_index]))/length(index_grid_observed) #  0.4789916
length(which(grid_numeric[index_grid_observed] < grid_numeric[maxufg_index]))/length(index_grid_observed) # 0.5210084

maxnumeric_index <- which(grid_numeric== median(grid_numeric[index_grid_observed]))
grid_spatial[maxnumeric_index, ]
grid_nominal[maxnumeric_index]
grid_spatial[maxnumeric_index, ][which(grid_nominal[maxnumeric_index] == 4), ]

marks_ppp_compare <- c(rep("Elevation median", length(maxnumeric_index)), "Max simplicial depth", "Max ufg depth")
ppp_compare <- ppp(x = c(grid_spatial[maxnumeric_index, 1], grid_spatial[maxsimpldd_index, 1], grid_spatial[maxufg_index, 1]),
                   y = c(grid_spatial[maxnumeric_index, 2], grid_spatial[maxsimpldd_index, 2], grid_spatial[maxufg_index, 2]),
                   marks= marks_ppp_compare,
                   window = gor$window)

plot(ppp_compare,
     main = "Different Centrality Approaches",
     las =2)

# data information
gor_medians <- superimpose(gor, list(x = c(grid_spatial[maxsimpldd_index, 1], grid_spatial[maxufg_index, 1]),
                                     y = c(grid_spatial[maxsimpldd_index, 2], grid_spatial[maxufg_index, 2])))
marks(gor_medians) <- c(rep("Observation", gor$n), "Simplicial-Median", "Ufg-Median")
plot(gor_medians,
     cols = c("black", "red", "blue"),
     pch = c(1, 19, 19),
     ces = 0.5, main = NULL,
     leg.side = "right")


# # Gegenbsp starshaped und quasiconcavity --> old
# ################################################################################
# library(MASS)
# set.seed(837)
# p_1 <- MASS::mvrnorm(n = 1000, mu = c(0,0), Sigma = diag(c(0.005, 0.005)))
# p_2 <- MASS::mvrnorm(n = 1000, mu = c(0.5,1), Sigma = diag(c(0.005, 0.005)))
# p_3 <- MASS::mvrnorm(n = 1000, mu = c(1,0), Sigma = diag(c(0.005, 0.005)))
# # p_4 <- MASS::mvrnorm(n = 10, mu = c(1,1),  Sigma = diag(c(0.01, 0.01)))
#
# grid_axis <- seq(-0.5, 1.5, 0.01)
#
# observed <- rbind(rbind(p_1, p_2), p_3)
# grid_spatial <- expand.grid(grid_axis, grid_axis)
# grid_observed_x <- unlist(lapply(observed[, 1], FUN = function(x) {which.min(abs(grid_axis - x))}))
# grid_observed_y <- unlist(lapply(observed[, 2], FUN = function(y) {which.min(abs(grid_axis - y))}))
# index_grid_observed <- (grid_observed_y - 1) * length(grid_axis) + grid_observed_x
#
#
# erg_star <- depth.simplicial(grid_spatial, observed, exact = TRUE)
#
# matrix_erg <- transmat(matrix(erg_star, nrow = length(grid_axis)),
#                        from = "European", to = "Cartesian")
# matrix_erg <- transmat(matrix_erg, from = "European", to = "spatstat")
# im_starsh <- im(matrix_erg, xrange = c(-1,2), yrange = c(-1, 2))
# par(mfrow = c(1,1))
# plot(im_starsh) # heat.colors(10) oder  topo.colors(10)
# # points(observed, pch = 20)

