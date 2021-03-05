## daily_timestep_yearly_obs_setup.R
## generate data and set up MCMC

library(zoo)
library(coda)
library(imputeTS)
library(Rcpp)
library(truncnorm)
library(tmvtnorm)
## read in driver data
drivers_ev = read.csv('./temp_data.csv', header = T)
drivers_ev$mint = na.interpolation(drivers_ev$mint)
drivers_ev$maxt <- na.interpolation(drivers_ev$maxt)
drivers_ev$rad <- na.interpolation(drivers_ev$rad)
drivers_ev$rad <- drivers_ev$rad * .0864
## Reads in initial conditions
drivers_ev <- do.call(rbind, replicate(2, coredata(drivers_ev), simplify = FALSE))
drivers_ev$projectday <- as.numeric(rownames(drivers_ev))

nday <- 730
drivers_ev$ca <- c(rep(396.52, 365), rep(398.65, 365))

drivers_ev <- drivers_ev[1:nday,]

## set parameter values
p_e = c(1/500, 0.33, 0.15, 0.27, 1/(2*365), 1/(25*365), 1/(2*365), .1, 1/(250* 365), 0.1725, 3)

## transform p1 and p8
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

init_var <- c(1,4,.5,.5,4)
init_mean <- c(100,9200,100,20,11000)
var_add <- c(2, 9, (2)^2, .25, 16)
var_obs <- c(4, 16, 1^2, 1^2, 25)
prec_add <- 1 / var_add
prec_obs <- 1 / var_obs
prec_init <- 1 / init_var
## Sources R scripts to generate data
source('./multiyearDALEC.R')
source('./MCMC_funs.R')

library(Rcpp)
library(truncnorm)
library(tmvtnorm)
## Load matrixStats library for colMeans2
library(matrixStats)
## Load laGP for GP initializations of latent states
library(laGP)
library(lhs)
library(HDInterval)

## Compile C++ scripts from src folder 
sourceCpp('./src/lmvnd_ll.cpp')
sourceCpp('./src/likelihood_SS.cpp')
sourceCpp('./src/lgamma.cpp')
sourceCpp('./src/Cpred.cpp')
## Load required R functions
source('./cpp_ss_likelihood.R')
source('./get_latent_sample.R')
source('./latent_initialization.R')

## ------------ Sampler Test Run -----------##
# Set carbon stock variables from data generation scripts
Cf <- Cf_e
Cw <- Cw_e
Cr <- Cr_e
Clit <- Clit_e
Csom <- Csom_e
# Set driver variables
mint <- drivers_ev$mint
maxt <- drivers_ev$maxt
lat <- lat_e 
Nit <- Nit_e 
LMA <- LMA_e
Ca <- drivers_ev$ca 
yearday <- drivers_ev$yearday
rad <- drivers_ev$rad
# Set chain length and burn in for MCMC
chain_length <- 10000
burn <- 2000
#num_gen <- 1
options(warn = 2)
## ----- Preparing MCMC ------ ##
# Compute LAI estimate 
LAI <- pmax(.1, Cf.obs/LMA_e)
# Column bind matrices and build initial Cpred matrix
C_obs <- cbind(Cf.obs, Cw.obs, Cr.obs, Clit.obs, Csom.obs)
missing_data_ind = seq(from = 1, to = nday, by = 365)
missing_data_ind = c(missing_data_ind, nday)
C_obs[-missing_data_ind,] <- NA
C <- cbind(Cf, Cw, Cr, Clit, Csom)

# Initalize matrix for p draws to be stored
pars <- matrix(NA, nrow = chain_length, ncol = 16)

#G <- array(NA, dim = c(chain_length, 730, num_gen))
G <- matrix(NA, nrow = chain_length, ncol = nday)
G[1,] <- G_obs

phi_g <- 1 / .025
tau_g <- 1 / .04
# create matrix for acceptance for p_i 


# Initialize arrays for C_samples
C_samples <- array(NA, dim = c(nday, 5, chain_length))
C.0 <- array(NA, dim = c(chain_length, 5))
C.0[1,] <- init_mean
C_samples[1,,1] <- init_mean
# Initialize first entry to observed data
C_samples[, , 1] <- latent_initialization(C_obs, missing_data_ind, 1, nday, type = c('pwl', 'pwl', 'pwl', 'pwl', 'pwl'))

## function to update carbon stocks for particle filter
stock_updates <- function(stocks, p, G, drivers, prec, i){
  C <- stocks
  C_update <- rep(NA, dim(stocks)[2])
  ## average temperature for a day
  tbar <- .5*(drivers$mint[i] + drivers$maxt[i])
  ## update cf
  C_update[1] <- (1 - p[5])*C[i-1,1] + G[i]*(1-p[2])*p[3] + rnorm(1,0,1/sqrt(prec[1]))
  ## update cw
  C_update[2] <- (1 - p[6])*C[i-1,2] + G[i]*(1 - p[2])*(1 - p[3])*(1 - p[4]) + rnorm(1,0,1/sqrt(prec[2]))
  ## update cr
  C_update[3] <- (1 - p[7])*C[i-1,3] + G[i]*(1 - p[2])*(1 - p[3])*p[4] + rnorm(1,0,1/sqrt(prec[3])) 
  ## update clit
  C_update[4] <- (1 - .5*exp(p[10]*tbar)*(p[1] + p[8]))*C[i-1,4] + p[5]*C[i-1,1] + p[7]*C[i-1,3] + rnorm(1,0,1/sqrt(prec[4]))
  ## update csom
  C_update[5] <- (1 - .5*exp(p[10]*tbar)*p[9])*C[i-1,5] + .5*exp(p[10]*tbar)*p[1]*C[i-1,4] + p[6]*C[i-1,2] + rnorm(1,0,1/sqrt(prec[5]))
  return(C_update)
}

particle_update <- function(y, prec, numpart, init_mean, init_prec, updatefun, obs_prec, p, ...){
  x <- array(NA, dim = c(numpart, dim(y)[1], dim(y)[2]))
  for (k in 1:dim(x)[3]){
    x[,1,k] <- rnorm(numpart, mean = init_mean[k], sd = 1 / sqrt(1 / init_prec[k]))
    w <- dnorm(x[,1,k], y[1,k], sd = 1/sqrt(init_prec[k]))
    x[,1,k] <- sample(x[,1,k], numpart, prob = w, replace = TRUE) 
  }
  for (i in 2:dim(y)[1]){
    ## update states
    for (j in 1:numpart){
      x[j,i,] <- updatefun(stocks = x[j,,], i = i, prec = prec, p = p, ...)
      
    }
    if (all(!(is.na(y[i,])))){
      for (k in 1:dim(x)[3]){
        w <- dnorm(x[,i,k], mean = y[i,k], sd = 1 / sqrt(obs_prec[k]))
        if (sum(w) == 0) w <- rep(1, 100)
        x[,i,k] <- sample(x[,i,k], numpart, prob = w, replace = TRUE) 
      }
    }
  }
  return(x)
}
plower <- c(1e-6, .2, .01, .01, 1e-4, 1e-6, 1e-6, 1e-5, 1e-6, .05, 2, 0, 0, 0, 0, 0)
pupper <- c(.01, .7, .5, .5, .1, .01, .01, .1, .01, .2, 20, NA, NA, NA, NA, NA)

pars <- matrix(NA, nrow = 1, ncol = 16)
pars[1,] <- c(.5*(pupper[1:11] + plower[1:11]), rep(1,5))
pars[1,1:11] <- c(.002, .35, .13, .3, .0015, .00012, .0015, .00025, .000012, .17, 3.5)

pft <- particle_update(y = C_obs, prec = prec_add, numpart = 100,
                       init_mean = init_mean, init_prec = prec_init, updatefun = stock_updates,
                       G = G_e, p = pars[1,1:11], drivers = drivers_ev, obs_prec = prec_obs)
C_samples[,1,1] <- colMeans2(pft[,,1])
C_samples[,2,1] <- colMeans2(pft[,,2])
C_samples[,3,1] <- colMeans2(pft[,,3])
C_samples[,4,1] <- colMeans2(pft[,,4])
C_samples[,5,1] <- colMeans2(pft[,,5])
# Set vector of upper and lower bounds for accept reject sampling
acceptance <- matrix(0, nrow = chain_length, ncol = 16)
fluxes_obs <- list(Lr_obs = Lr_obs, D_e_obs = D_e_obs, Rh1_obs = Rh1_obs, Rh2_obs = Rh2_obs, 
                   Ar_obs = Ar_obs, Af_obs = Af_obs, Lf_obs = Lf_obs, Ra_obs = Ra_obs, Lw_obs = Lw_obs)

sd <- rep(0, 16)
sd[11] <- .05
sd[1] <- .01
sd[2] <- .002
sd[3] <- .002
sd[4] <- .001
sd[5] <- .03
sd[6] <- .3 / 400
sd[7] <- .016
sd[8] <- .03
sd[9] <- .01
sd[10] <- .008
sd[11] = .009
sd[12:16] <- .17
colnames(pars) <- c(paste0('p', 1:11), paste0('precadd', 1:5))
old <- cpp_ss_likelihood_ex(Cobs = C_obs, C = C_samples[,,1], var_obs = var_obs, pars = pars, init_mean = init_mean,
                            init_sd = init_var, N = nday, G = G[1,], G_e = G_obs, drivers = drivers_ev, 
                            LAI = pmax(C_samples[,1,1]/LMA_e), fluxes = fluxes_obs)
accept <- rep(0, 16)
cov.block <- NULL
thresh2 <- matrix(0, ncol = 16, nrow = chain_length)
p_star_mat <- matrix(0, ncol = 16, nrow = chain_length)
update = TRUE
block <- c(2,3,4,11)

block <- vector('list', 4)
block[[1]] <- c(11)
block[[2]] <- c(1,8)
block[[3]] <- c(9,10)
block[[4]] <- c(2,3,4)
burn <- 2000