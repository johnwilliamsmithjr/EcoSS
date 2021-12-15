## file: monthly_timestep_monthly_obs_MCMC.R
## sets up MCMC to be run for the monthly latent state model
## with monthly observation resolution

## source MCMC functions
source('./setup/MCMC/MCMC_funs.R')

## Compile C++ scripts from src folder 
sourceCpp('./src/lmvnd_ll.cpp')
sourceCpp('./src/likelihood_SS.cpp')
sourceCpp('./src/lgamma.cpp')
sourceCpp('./src/Cpred_rp2.cpp')
## Load required R functions
source('./setup/MCMC/cpp_ss_likelihood.R')
source('./setup/MCMC/get_latent_sample.R')
source('./setup/MCMC/latent_initialization.R')

## likelihood
ss_monthly_ll <- function(C_obs, C, init_mean, init_var, prec_obs, pars, fluxes, G, ndays, drivers, init_ev, r = 1, ...){
  native_step = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6))
  p = pars[1:11]
  prec_add = pars[12:16]
  ## Initial Condition Likelihood 
  LL <- sum(dnorm(x = C_obs[1,], mean = init_mean, sd = sqrt(init_var), log = TRUE))
  ## Build Cpred
  LAI_ = Cf_approx / LMA_e
  preds <- getDALECpredictions2(p_e = p, C = C, LAI = LAI_, maxt = drivers$maxt[1:nday],
                                mint = drivers$mint[1:nday], G_obs = fluxes$G, ind = native_step-1, init_mean = init_mean, native_step = native_step,
                                ...)
  na_ind <- which(is.na(C_obs[,1]))
  Cpred <- preds$Cslow
  Cpred2 <- preds$Cfast
  
  ## Process Model Likelihood
  for (i in 1:dim(C_obs)[2]){
    LL <- LL + sum(dnorm(x = C[,i], mean = Cpred[,i], sd = 1 / sqrt(prec_add[i]), log = TRUE))
  }
  ## Observation Model Likelihood
  obs_ind <- which(!is.na(C_obs[,1]))
  
  if (1 %in% obs_ind) obs_ind <- obs_ind[-which(obs_ind == 1)]
  maxt <- drivers$maxt
  mint <- drivers$mint
  for (i in obs_ind){
    LL <- LL + sum(dnorm(x = C_obs[i,], mean = C[i,], sd = sqrt(1 / prec_obs), log = TRUE))
  }
  
  maxt <- maxt[1 : max(cumsum(ndays))]
  mint <- mint[1 : max(cumsum(ndays))]
  ## Flux Likelihood
  LL = LL + Lr_ll(p = p, C = Cpred2, fluxes$Lr_obs)
  #print(LL)
  LL = LL + De_ll(p = p, C = Cpred2, maxt = maxt, mint = mint,
                  D_e_obs = fluxes$D_e_obs)
  
  LL = LL + Rh1_ll(p = p, C = Cpred2, Rh1_obs = fluxes$Rh1_obs, 
                   maxt = maxt, mint = mint)
  
  LL = LL + Rh2_ll(p = p, C = Cpred2, Rh2_obs = fluxes$Rh2_obs, 
                   maxt = maxt, mint = mint)
  
  LL = LL + Ar_ll(p = p, G = G, Ar_obs = fluxes$Ar_obs)
  
  LL = LL + Af_ll(p = p, G = G, Af_obs = fluxes$Af_obs)
  
  LL = LL + Lf_ll(p = p, C = Cpred2, Lf_obs = fluxes$Lf_obs)
  
  LL = LL + Ra_ll(p = p, G = G, Ra_obs = fluxes$Ra_obs)
  
  LL = LL + Lw_ll(p = p, Lw_obs = fluxes$Lw_obs, C = Cpred2)
  
  LL = LL + sum(dnorm(preds$G, G_obs, sd = .001, log = TRUE))
  
  ## Priors
  prior <- sum(log(1 / prec_add))
  return(r*LL + prior)
}

## ensure that the proper r is being used for the run

if (!('r' %in% ls())){
  r <- 1
  cat('No r value supplied for datacloning, using default (r=1)\n')
}

## define getDALECpredictions2, a function to extract predictions on the fast and slow timescale to be fed into the MCMC routine and likelihood
getDALECpredictions2 <- function(p_e, C, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, N, G_obs, ind, init_mean, native_step, preds = NULL){
  ## get predictions from Cpred_m.cpp
  preds <- Cpred_m(p_e, C, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, N, G_obs, ind, init_mean)
  ## store fast timescale and slow timescale predictions
  Cfast <- cbind(preds$Cfpred, preds$Cwpred, preds$Crpred, preds$Clitpred, preds$Csompred)
  Cslow <- cbind(preds$Cfpred[native_step], preds$Cwpred[native_step], preds$Crpred[native_step], preds$Clitpred[native_step], preds$Csompred[native_step])
  ## store gpp
  G <- preds$G
  ## store tlapm coefficients for native timestep
  A <- cbind(preds$A_cf, preds$A_cw, preds$A_cr, preds$A_clit, preds$A_csom)
  b <- cbind(preds$b_cf, preds$b_cw, preds$b_cr, preds$b_clit, preds$b_csom)
  #ind = ind + 1
  ## initialize tpalm coefficients for new timestep
  newA <- matrix(NA, ncol = 5, nrow = length(native_step))
  newB <- matrix(NA, ncol = 5, nrow = length(native_step))
  ## create coefficients
  for (i in 1:length(native_step)){
    for (j in 1:5){
      if (i == 1){
        newA[i,j] <- prod(A[1:native_step[i], j])
        bs <- 0
        for (t in 1:(native_step[i]-1)){
          bs <- bs + prod(A[(t+1):native_step[i], j]) * b[t,j]
        }
        newB[i,j] <- bs + b[native_step[i], j]
      } else{
        newA[i,j] <- prod(A[(native_step[i-1]+1) : native_step[i], j])
        bs <- 0
        for (t in (native_step[i-1]+1):(native_step[i]-1)){
          bs <- bs + prod(A[(t+1):native_step[i], j]) * b[t,j]
          newB[i,j] <- bs + b[native_step[i], j]
        }
      }
    }
  }
  ## output
  return(list(A = newA, b = newB, Cslow = Cslow, Cfast = Cfast, G = G, oldA = A, oldB = b))
}

## create carbon stock observation dataframe
C_obs <- cbind(Cf.obs, Cw.obs, Cr.obs, Clit.obs, Csom.obs)
C_obs <- C_obs[1:nmonths,]

## create flux list
fluxes_obs <- list(Lr_obs = Lr_obs, D_e_obs = D_e_obs, Rh1_obs = Rh1_obs, Rh2_obs = Rh2_obs, 
                   Ar_obs = Ar_obs, Af_obs = Af_obs, Lf_obs = Lf_obs, Ra_obs = Ra_obs, Lw_obs = Lw_obs, G = G_obs)

## store true latent states
Cf <- Cf_e
Cw <- Cw_e
Cr <- Cr_e
Clit <- Clit_e
Csom <- Csom_e

# Set driver variables
mint <- drivers_ev$mint[1:nday] 
maxt <- drivers_ev$maxt[1:nday] 
lat <- lat_e 
Nit <- Nit_e 
LMA <- LMA_e
Ca <- drivers_ev$ca[1:nday] 
yearday <- drivers_ev$yearday[1:nday] 
rad <- drivers_ev$rad[1:nday]

# Set chain length and burn in for MCMC

if (!('chain_length' %in% ls())){
  chain_length <- 10000
  cat('No changes detected in chain length, using default of 10000\n')
}
if (!('burn' %in% ls())){
  burn <- 2000
  cat('No changes detected in burnin period, using default of 2000\n')
}

## set warnings as errors
options(warn = 2)
## ----- Preparing MCMC ------ ##
# Compute LAI estimate to be fed into p_11 optimization routine
LAI <- pmax(.1, Cf.obs/LMA_e)
# Column bind matrices and build initial Cpred matrix
C_obs <- cbind(Cf.obs, Cw.obs, Cr.obs, Clit.obs, Csom.obs)

## reset observation index
missing_data_ind = seq(from = 1, to = nmonths, by = 1)
#missing_data_ind = c(missing_data_ind, nmonths)
ndays <- c(31, diff(cumsum(ndays)[missing_data_ind]))
C_obs[-missing_data_ind,] <- NA
C <- cbind(Cf, Cw, Cr, Clit, Csom)

# Initalize matrix for p draws to be stored
pars <- matrix(NA, nrow = chain_length, ncol = 16)
# Set vector of upper and lower bounds for accept reject sampling
plower <- c(1e-6, .2, .01, .01, 1e-4, 1e-6, 1e-6, 1e-5, 1e-6, .05, 1, 0, 0, 0, 0, 0)
pupper <- c(.01, .7, .5, .5, .1, .01, .01, .1, .01, .2, 20, NA, NA, NA, NA, NA)

## transform p1 and p8 (model reparameterization from standard DALECev)
plower[1] <- plower[1] + plower[8]
pupper[1] <- pupper[1] + pupper[8]
plower[8] <- 1e-6
pupper[8] <- 1
#G <- array(NA, dim = c(chain_length, nday, num_gen))
G <- matrix(NA, nrow = chain_length, ncol = nday)
G[1,] <- G_obs

## set g observation and error variance
phi_g <- 1 / (1e-6)
tau_g <- 1 / (.05)^2

# Initialize arrays for C_samples and initial conditions
C_samples <- array(NA, dim = c(nmonths, 5, chain_length))
C.0 <- array(NA, dim = c(chain_length, 5))
C.0[1,] <- init_mean
C_samples[1,,1] <- init_mean
# Initialize first entry to observed data

## latent state initializations

if (!('stock_inits' %in% ls())){
  if (!('init_type' %in% ls())){
    C_samples[, , 1] <- latent_initialization(C_obs, missing_data_ind, 1, nmonths, type = c('pwl', 'pwl', 'pwl', 'pwl', 'pwl'))
    cat('Initializing latent states using default, piece wise linear interpolation\n')
  } else if (init_type == 'pwl'){
    C_samples[, , 1] <- latent_initialization(C_obs, missing_data_ind, 1, nmonths, type = c('pwl', 'pwl', 'pwl', 'pwl', 'pwl'))
    cat('Initializing latent states using piece wise linear interpolation\n')
  } else if (init_type == 'gp'){
    C_samples[, , 1] <- latent_initialization(C_obs, missing_data_ind, 1, nmonths, type = c('gp', 'gp', 'gp', 'gp', 'gp'))
    cat('Initializing latent states using Gaussian Process\n')
  } else{
    cat('Initialization type not recognized, please use one of pwl or gp for latent states\n')
  }
} else{
  C_samples[,,1] <- stock_inits
}

## Parameter initializations

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

## initialize proposal stepsizes

if (!('sd' %in% ls())){
  sd <- rep(0, 16)
  sd[11] <- .02
  sd[1] <- .01
  sd[2] <- .002*5
  sd[3] <- .002*5
  sd[4] <- .001*5
  sd[5] <- .03
  sd[6] <- .3
  sd[7] <- .008
  sd[8] <- .03
  sd[9] <- .01
  sd[10] <- .008
  sd[11] = .009
  sd[12:16] <- .17
  cat('No changes detected in proposal standard deviations, using default values\n')
}

## set column names
colnames(pars) <- c(paste0('p', 1:11), paste0('precadd', 1:5))


preds <- getDALECpredictions2(p_e = pars[1:11], C = C_samples[,,1], LAI = Cf_approx / LMA_e, maxt = maxt, mint = mint,
                              Ca = Ca, lat_e = lat_e, yearday = yearday, Nit_e = Nit_e, rad = rad, N = nday, G_obs = G_e, ind = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6))-1, 
                              init_mean = C.0[1,], native_step = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6)))


old <- ss_monthly_ll(C_obs = C_obs, C = C_samples[,,1], init_mean = init_mean, init_var = init_var, prec_obs = prec_obs, pars = pars, fluxes = fluxes_obs,
                     G = G[1,], drivers = drivers_ev, init_ev = init_ev, ndays = ndays, Ca = Ca, lat_e = lat_e, yearday = yearday, rad = rad, N = nday, Nit_e = Nit_e, r = r)

## initialize acceptance matrix and other diagnostic tests
acceptance <- matrix(0, nrow = chain_length, ncol = 16)
accept <- rep(0, 16)
cov.block <- NULL
thresh2 <- matrix(0, ncol = 16, nrow = chain_length)
p_star_mat <- matrix(0, ncol = 16, nrow = chain_length)

if (!('update' %in% ls())){
  update = TRUE
  cat('No option detected for MCMC update, using default (update = TRUE)\n')
}

if (!('block' %in% ls())){
  block <- vector('list', 2)
  block[[1]] <- c(2,3,4)
  block[[2]] <- c(9, 10)
  cat('No changes detected in block updates, using default blocks\n')
}

preds <- getDALECpredictions2(p_e = p_e, C = C, LAI = Cf_approx / LMA_e, maxt = maxt, mint = mint,
                              Ca = Ca, lat_e = lat_e, yearday = yearday, Nit_e = Nit_e, rad = rad, N = 2190, G_obs = G_e, 
                              ind = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6))-1, 
                              init_mean = init_mean, native_step = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6)))

test <- cbind(Cf_approx, Cr_approx, Cw_approx, Clit_approx, Csom_approx)
Cpredi <- preds$Cslow

## specify transformation
tf <- rep('log', 16)