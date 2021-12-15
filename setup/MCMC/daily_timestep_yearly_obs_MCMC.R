## daily_timestep_yearly_obs_MCMC.R
## sets up necessary arrays, initializations to perform MCMC

source('./setup/MCMC/MCMC_funs.R')
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
source('./setup/MCMC/cpp_ss_likelihood.R')
source('./setup/MCMC/get_latent_sample.R')
source('./setup/MCMC/latent_initialization.R')

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
if (!('chain_length' %in% ls())){
  chain_length <- 10000
  cat('No changes detected in chain length, using default of 10000\n')
}
if (!('burn' %in% ls())){
  burn <- 2000
  cat('No changes detected in burnin period, using default of 2000\n')
}
options(warn = 2)
## ----- Preparing MCMC ------ ##
# Compute LAI estimate 
LAI <- pmax(.1, Cf.obs/LMA_e)
# Column bind matrices and build initial Cpred matrix
C_obs <- cbind(Cf.obs, Cw.obs, Cr.obs, Clit.obs, Csom.obs)

if (!('observation_index' %in% ls())){
  missing_data_ind = seq(from = 1, to = nday, by = 365)
  cat('No changes detected in observation index, using yearly observation data\n')
} else{
  missing_data_ind = observation_index
}
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


# Initialize arrays for C_samples
C_samples <- array(NA, dim = c(nday, 5, chain_length))
C.0 <- array(NA, dim = c(chain_length, 5))
C.0[1,] <- init_mean
C_samples[1,,1] <- init_mean
# Initialize first entry to observed data

if (!('stock_inits' %in% ls())){
  if (!('init_type' %in% ls())){
    C_samples[, , 1] <- latent_initialization(C_obs, missing_data_ind, 1, nday, type = c('pwl', 'pwl', 'pwl', 'pwl', 'pwl'))
    cat('Initializing latent states using default, piece wise linear interpolation\n')
  } else if (init_type == 'pwl'){
    C_samples[, , 1] <- latent_initialization(C_obs, missing_data_ind, 1, nday, type = c('pwl', 'pwl', 'pwl', 'pwl', 'pwl'))
    cat('Initializing latent states using piece wise linear interpolation\n')
  } else if (init_type == 'gp'){
    C_samples[, , 1] <- latent_initialization(C_obs, missing_data_ind, 1, nday, type = c('gp', 'gp', 'gp', 'gp', 'gp'))
    cat('Initializing latent states using Gaussian Process\n')
  } else{
    cat('Initialization type not recognized, please use one of pwl or gp for latent states\n')
  }
} else{
  C_samples[,,1] <- stock_inits
}


## lower and upper bounds of parameters
plower <- c(1e-6, .2, .01, .01, 1e-4, 1e-6, 1e-6, 1e-5, 1e-6, .05, 1, 0, 0, 0, 0, 0)
pupper <- c(.01, .7, .5, .5, .1, .01, .01, .1, .01, .2, 20, NA, NA, NA, NA, NA)

pars <- matrix(NA, nrow = 1, ncol = 16)

if (!('par_inits' %in% ls())){
  pars[1,] <- c(.5*(pupper[1:11] + plower[1:11]), rep(1,5))
  cat('No changes detected to parameter initializations, using default\n')
} else{
  if (length(par_inits) != 16){
    cat('Initializations for parameters must be length 16, please correct this\n')
  }
  pars[1,] <- par_inits
}


# Set vector of upper and lower bounds for accept reject sampling
acceptance <- matrix(0, nrow = chain_length, ncol = 16)
fluxes_obs <- list(Lr_obs = Lr_obs, D_e_obs = D_e_obs, Rh1_obs = Rh1_obs, Rh2_obs = Rh2_obs, 
                   Ar_obs = Ar_obs, Af_obs = Af_obs, Lf_obs = Lf_obs, Ra_obs = Ra_obs, Lw_obs = Lw_obs)

## checks for user tuned standard deviations for MH
if (!('sd' %in% ls())){
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
  cat('No changes detected in proposal standard deviations, using default values\n')
}

## set colnames of parameters
colnames(pars) <- c(paste0('p', 1:11), paste0('precadd', 1:5))

## compute initial likelihood for MCMC
old <- cpp_ss_likelihood_ex(Cobs = C_obs, C = C_samples[,,1], var_obs = var_obs, pars = pars, init_mean = init_mean,
                            init_sd = init_var, N = nday, G = G[1,], G_e = G_obs, drivers = drivers_ev, 
                            LAI = pmax(C_samples[,1,1]/LMA_e), fluxes = fluxes_obs)
## initialize diagnostic vectors for acceptance, block covariance, threshold of samples, proposed values, 
## and whether or not to give updates
accept <- rep(0, 16)
cov.block <- NULL
thresh2 <- matrix(0, ncol = 16, nrow = chain_length)
p_star_mat <- matrix(0, ncol = 16, nrow = chain_length)

if (!('update' %in% ls())){
  update = TRUE
  cat('No option detected for MCMC update, using default (update = TRUE)\n')
}

## check if block is given as user argument, otherwise use default block

if (!('block' %in% ls())){
  block <- vector('list', 4)
  block[[1]] <- c(11)
  block[[2]] <- c(1,8)
  block[[3]] <- c(9,10)
  block[[4]] <- c(2,3,4) 
  cat('No changes detected in block updates, using default blocks\n')
}