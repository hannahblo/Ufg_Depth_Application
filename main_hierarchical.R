################################################################################
# Setup
################################################################################
# Allbus Analysis
library(foreign)
library(plyr)
library(ggplot2)
library(ggsci)
library(scales)
library(ddandrda)

# download the Allbus 2021 data
# ( https://search.gesis.org/research_data/ZA5280 - last accessed Nov 2024)


################################################################################
#  needed additional functions:
################################################################################

# This function generates from hierarchical nominal data that are given in the
# format like the ISCO-08 format of the Allbus data set (ZA5280_v2-0-1.sav) the
# corresponding formal context
compute_hierarchical_scaling_vec <- function(values){

  data_values <- rep("",length(values))
  for (k in (1:length(data_values))) {
    # Error handling/ threatment of leading zeros that may be omitted otherwise
    if (nchar(values[k]) == 3) {data_values[k] <- paste("0",
                                                      as.character(values[k]),
                                                      sep = "") }
    if (nchar(values[k]) == 4) {data_values[k] <- as.character(values[k])}
    if (nchar(values[k]) <= 2) {print("ERROR"); return(NULL)}
  }
  m <- length(data_values)
  values <- rep("",m*4)
  t <- 1
  for (l in c(1,2,3,4)) {
    for (k in (1:m)) {
      values[t] <- substr(as.character(data_values[k]), 1, l)
      t <- t + 1
    }
  }
  unique_values <- unique(values)
  X <- NULL
  for (value in unique_values) {
    temp <- rep(0,m)
    L <- nchar(value)
    for (k in (1:m)) {
      temp[k] <- substr(as.character(data_values[k]), 1, L) == value
    }
    X <- cbind(X,temp)
  }
  colnames(X) <- unique_values
  return(X)}



################################################################################
# Computation of the ufg-depth
################################################################################



# read the Allbus 2021 data ( https://search.gesis.org/research_data/ZA5280 )
# For this please change working directory accordingly
dat <- read.spss("ZA5280_v2-0-1.sav", use.value.labels = FALSE)

# select only persons, for which the ISCO-08 status (occupational status
# according to ISCO-08) is available (dat$isco08>0)
indexs <- which(dat$isco08 > 0)
# data points coded with 4 digits from 0-9
# Note that in the original data leading zeros are omitted, the function
# 'compute_hierarchical_scaling_vec' cares for that
x <- dat$isco08[indexs]
# different weights are used for east and west because of an oversampling of
# east germany
weights <- dat$wghthew[indexs]
context <- compute_hierarchical_scaling_vec(x)
rownames(context) <- x
# make a weighted data table (first column contains the names of the
# hierarchical categories, y_weighted are the weighted counts)
weighted_repr <- ddandrda::get_weighted_representation(x = cbind(x,context),
                                                     y = weights)
# Construct the weighted context without ties
weighted_context <- weighted_repr$x[,-1]
rownames(weighted_context) <- weighted_repr$x[,1]
# extract the weighted counts
y_weighted <- weighted_repr$y_weighted
# x_weighted contains the categories for the weighted context
x_weighted <- weighted_repr$x_weighted[,1]



# compute contributions of 1-element and 2-element ufg premises, as well as the
# final vector that contains ufg-depths
ufg_1 <- ddandrda::ufg_1_depth_hierarchical(weighted_context,weights = y_weighted)
ufg_2 <- ddandrda::ufg_2_depth_hierarchical(weighted_context,weights = y_weighted)
# Compute ufg-depth with C_1=C_2=1
ufg   <- ufg_1$depths + ufg_2$depths








# data points with highest depths:
i1 <- which(ufg_1$depth == max(ufg_1$depths))
x_weighted[i1]
# [1] 4110
#deepest point w.r.t. ufg_1: 4110: General Office Clerks

i2 <- which(ufg_2$depths == max(ufg_2$depths))
x_weighted[i2]
# [1] 3341 3343 3342 3344
# deepest points w.r.t. ufg_2
# (all subcategories of category 334: Administrative and Specialized Secretaries)
#
# 3341: Office Supervisors
# 3342: Legal Secretaries
# 3343: Administrative and Executive Secretaries
# 3344: Medical Secretaries


i <- which(ufg == max(ufg))
x_weighted[i]
# [1] 3221
#deepest point w.r.t. ufg depth: 3221: Nursing Associate Professionals

max(ufg)
# maximal depth
# [1] 0.9266822


# Further analysis
length(unique(ufg))
# 285 unique depth values



length(unique(ddandrda::compute_quasiconcave_hull(ufg,weighted_context)))
# [1] 3 contour sets sets of D^{qc}



# modus on the finest hierarchical level
(table(x))[which.max(table(x))]

# 4110
#  112 ( 112 cases with category 4110)
# categorial modus: 4110: General Office Clerks (of course identical to the
# ufg-1 mode of line 181)

# Top down approach: Iterative mode starting from the most coarse level followed
# by successively looking at more fine-grained levels
j <- seq_len(nrow(context))
for (k in (1:4)) {
  o <- order(colSums(context[j,]), decreasing = TRUE)
  j <- which(context[,o[k]] == 1)
}
colnames(context)[o[4]]
# Median according to the top down approach
# [1] "3343" Administrative and Executive Secretaries


# smallest ufg-depth:

which(ufg == min(ufg))
# [1] 281
x_weighted[281]
# [1] 6210
# datapoint with smallest ufg-depth: 6210: Forestry and Related Workers

# computation of the (corresponding intents of the) upper level sets
depth_values = sort(unique(ufg), decreasing = TRUE)

for (threshold in depth_values) {
  indexs <- which(ufg >= threshold)
  extent <- rep(0,length(x_weighted))
  extent[indexs] <- 1
  intent <- ddandrda::calculate_psi(extent, weighted_context)
  print(colnames(weighted_context)[which(intent == 1)])
  if (all(intent == 0)) {break}
}

#'
#
# 
# generalized Tukey depth:
 D_tukey <- ddandrda::compute_tukeys_depth(context, context,
                                         row_weights=weights)

table(D_tukey)

#D_tukey
# 0.709984634088231 0.746869964580187 
#             1913               787 






################################################################################
# Histograms:
################################################################################

# function for logarithmic scaling
f <- function(x){log(1 + 0.1*x)}

# Make Boxplot of data values between x_min and x-max
x_min <- 3000
x_max <- 5000
indexs <- which(x >= x_min & x <= x_max)
df <- data.frame(x = x[indexs])
lwd1 <- 1
lwd2 <- 0.8
cols = (pal_bmj(palette = c("default"), alpha = 1))(n = 6)
CC = 1/1056
pdf("allbusplot1.pdf")
ggplot(df, aes(x = x))  +
  theme(axis.text.x = element_text(size = 10, angle = 90)) +
  scale_x_continuous(breaks = seq(3000, 4999, by = 100)) +
  geom_vline(alpha = 0.5, xintercept = 3343, col = cols[4], lwd = lwd2) +
  geom_vline(xintercept = 3221,alpha = 0.5,col = cols[5], lwd = lwd2) +
  geom_vline(xintercept = 4110, col = cols[6], lwd = lwd2, alpha = 0.5) +
  geom_histogram(aes(y = after_stat(f(count))),
                 breaks = seq(x_min,x_max,1000), col = "black", lwd = lwd1,
                 alpha = 0.3) +
  labs(y = "counts (log scale: y=ln(1+0.1*counts))",
       x = "ISCO-08 classification of occupation") +
  geom_histogram(aes(y = after_stat(f(count))), breaks = seq(x_min,x_max,100),
                 alpha = 0.2, col = cols[1], lwd = lwd1) +
  geom_histogram(aes(y = after_stat(f(count))),
                 breaks = seq(x_min, x_max, 10), col = cols[2], alpha = 0.5,
                 lwd = lwd1) +
  geom_histogram(aes(y = after_stat(f(count))), breaks = seq(x_min, x_max, 1),
                 col = cols[3], alpha = 1, lwd = lwd1)
dev.off()

# further graphic

x_min <- 3200
x_max <- 3400
indexs <- which(x >= x_min & x <= x_max)
df <- data.frame(x = x[indexs])
lwd1 <- 1
lwd2 <- 0.8
cols <- (pal_bmj(palette = c("default"), alpha = 1))(n = 6)

pdf("allbusplot2.pdf")


ggplot(df, aes(x = x))  +
  theme(axis.text.x = element_text(size = 10, angle = 90)) +
  scale_x_continuous(breaks = seq(3000, 4999, by = 100)) +
  geom_vline(alpha = 0.5, xintercept = 3343, col = cols[4], lwd = lwd2) +
  geom_vline(xintercept = 3221, alpha = 0.5, col = cols[5], lwd = lwd2)  +
  labs(y = "counts (log sacle: y=ln(1+0.1*counts))",
       x = "ISCO-08 classification of occupation")  +
  geom_histogram(aes(y = after_stat(f(count))),
                 breaks = seq(x_min, x_max,100), alpha = 0.2, col = cols[1],
                 lwd = lwd1) +
  geom_histogram(aes(y = after_stat(f(count))),
                 breaks = seq(x_min, x_max, 10), col = cols[2], alpha = 0.5,
                 lwd = lwd1) +
  geom_histogram(aes(y = after_stat(f(count))),
                 breaks = seq(x_min, x_max, 1), col = cols[3], alpha = 1,
                 lwd = lwd1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 16))
dev.off()


