####### Script Information ########################
# Brandon P.M. Edwards
# Project/Observer Effects on Detectability
# 1-simulate-data-example.R
# Created January 2023
# Last Updated January 2023

####### Import Libraries and External Files #######

library(bSims)
library(doParallel)
library(foreach)
library(dirmult)
library(detect)

# This is the (fairly) generic function to simulate point counts
source("functions/sim-pc.R")

####### Set Constants #############################

set.seed(seed = 12121216,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")


n_cores <- 3 # bSims can simulate datasets in parallel, so modify these as needed.
n_obs <- 50 # How many total observations in a data set?
n_sim <- 1 # How many datasets would you like to simulate?

phi <- 0.4 
tau <- 100

#' The density doesn't really matter here because we're just interested in 
#' recovering values of phi and/or tau
den <- 10 

#' n_protocols refers to number of binning methods. I didn't make it as flexible 
#' as I should have because this number controls both the number of time AND
#' distance methods
n_protocols <- 4

# Number of time bins and distance bins in each protocol/binning method
n_time_bins <- c(3, 4, 10, 7)
n_dist_bins <- c(3, 4, 7, 11)

# Maximum times and distances (columns) for each 
max_times <- matrix(data = c(3, 4, 5, NA, NA, NA, NA, NA, NA, NA, 
                             2, 4, 5, 6, NA, NA, NA, NA, NA, NA,
                             1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                             2, 3, 4, 5, 6, 7, 8, NA, NA, NA),
                    nrow = 4, byrow = TRUE)
max_dist <- matrix(data = c(50, 100, 400, NA, NA, NA, NA, NA, NA, NA, NA,
                            25, 50, 100, 150, NA, NA, NA, NA, NA, NA, NA, 
                            25, 50, 75, 100, 125, 150, 400, NA, NA, NA, NA, 
                            10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150),
                   nrow = 4, byrow = TRUE)

####### Simulate Datasets #########################

#' sim_pc() is the function wrapper I wrote to simulate a harmonized dataset.
#' In a nutshell, it will iterate over all binning methods (after randomly allocating
#' the total n_obs to each binning method) and simulate one dataset per method,
#' then combine them at the end. This function therefore returns a list of
#' dataframes:
#' * rem: the counts dataset under the removal sampling methods
#' * dis: the counts dataset under the distance sampling methods
#' * dist_design: the design matrix that will be fed into cmulti()
#' * time_design: the design matrix that will be fed into cmulti()
#' WARNING this function takes a while to run!! 
sim_data <- sim_pc(n_obs = n_obs,
                   phi = phi,
                   tau = tau,
                   den = den,
                   n_protocols = n_protocols,
                   n_time_bins = n_time_bins,
                   n_dist_bins = n_dist_bins,
                   max_times = max_times,
                   max_dist = max_dist,
                   n_cores = n_cores)

# Extracting the components listed above
rem_data <- sim_data$rem
rem_design <- sim_data$time_design
dis_data <- sim_data$dis
dis_design <- sim_data$dist_design

removal_model <- cmulti(as.matrix(rem_data[, 3:ncol(rem_data)]) | as.matrix(rem_design) ~ 1, type = "rem")
message(paste0("Modelled Phi: ", exp(coef(removal_model)[1]))) # should be around set value of phi

distance_model <- cmulti(as.matrix(dis_data[, 3:ncol(dis_data)]) | as.matrix(dis_design) ~ 1, type = "dis")
message(paste0("Modelled Tau: ", exp(coef(distance_model)[1]))) # should be around set value of tau

