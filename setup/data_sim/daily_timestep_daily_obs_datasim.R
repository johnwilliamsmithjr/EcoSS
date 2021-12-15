## daily_timestep_daily_obs_datasim.R
## generate data for daily timestep with daily observations
## load required libraries
library(zoo)
library(coda)
library(imputeTS)
library(Rcpp)
library(truncnorm)
library(tmvtnorm)

## read in driver data
drivers_ev = read.csv('./data/temp_data.csv', header = T)
drivers_ev$mint = na_interpolation(drivers_ev$mint)
drivers_ev$maxt <- na_interpolation(drivers_ev$maxt)
drivers_ev$rad <- na_interpolation(drivers_ev$rad)
drivers_ev$rad <- drivers_ev$rad * .0864

## Reads in initial conditions
drivers_ev <- do.call(rbind, replicate(2, coredata(drivers_ev), simplify = FALSE))
drivers_ev$projectday <- as.numeric(rownames(drivers_ev))

nday <- 730
drivers_ev$ca <- c(rep(396.52, 365), rep(398.65, 365))

drivers_ev <- drivers_ev[1:nday,]

## set parameter values, checks for user changes

if (!('p_e' %in% ls())){
  p_e = c(1/500, 0.33, 0.15, 0.27, 1/(2*365), 1/(25*365), 1/(2*365), .1, 1/(250* 365), 0.1725, 3)  
  cat('No changes detected in process parameters, using default parameter set\n')
}

if (!('init_var' %in% ls())){
  init_var <- c(1,4,.5,.5,4)  
  cat('No changes detected in IC precisions, using default precision set\n')
}

if (!('init_mean' %in% ls())){
  init_mean <- c(400,7500,700,100,8500)
  cat('No changes detected in IC means, using default initial means\n')
}

if (!('var_add' %in% ls())){
  var_add <- c(2, 9, (2)^2, .25, 16)  
  cat('No changes detected in process precisions, using default precision set\n')
}

if (!('var_obs' %in% ls())){
  var_obs <- c(4, 16, 1^2, 1^2, 25)
  cat('No changes detected in observation precisions, using default precision set\n')
}

prec_add <- 1 / var_add
prec_obs <- 1 / var_obs
prec_init <- 1 / init_var

## transform p1 and p8 for reparameterization
tmp1 <- p_e[1] * p_e[8]
tmp2 <- p_e[1] - tmp1
p_e[1] <- tmp2
p_e[8] <- tmp1

drivers_ev <- do.call(rbind, replicate(6, coredata(drivers_ev), simplify = FALSE))
drivers_ev$projectday <- as.numeric(rownames(drivers_ev))
nmonths <- 24*3
nday <- 730*3
drivers_ev$ca <- c(rep(396.52, 365), rep(398.65, 365), rep(400.83, 365),
                   rep(404.24, 365), rep(406.55, 365), rep(408.52, 365))
drivers_ev <- drivers_ev[1:nday,]

## Sources R scripts to generate data
source('./setup/data_sim/multiyearDALEC.R')