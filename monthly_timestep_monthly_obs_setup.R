## Program: monthly_timestep_monthly_obs.R
## Author: John W Smith Jr
## Date: 10/19/2018
# Load c++ functions through Rcpp
library(Rcpp)
library(coda)
library(truncnorm)
library(tmvtnorm)
library(matrixStats)
library(laGP)
## load lhs library for initial process parameter guesses
library(lhs)
## load HDInterval to obtain credible intervals
library(HDInterval)

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
r <- 1

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
## set working directory to proper directory and initialize variance
init_var <- c(1,4,.5,.5,4)

## source testMCMC.R to obtain mcmc update and likelihood functions
source('./MCMC_funs.R')

## Compile C++ scripts from src folder 
sourceCpp('./src/lmvnd_ll.cpp')
sourceCpp('./src/likelihood_SS.cpp')
sourceCpp('./src/lgamma.cpp')
sourceCpp('./src/Cpred_rp2.cpp')
## Load required R functions
source('./cpp_ss_likelihood.R')
source('./get_latent_sample.R')
source('./latent_initialization.R')

## load imputeTS to impute missing driver data
library(imputeTS)

## set initial conditions
init_var <- c(1,4,.5,.5,4)
init_mean <- c(400,7500,700,100,8500)

## set observation and additive (process) variance for simulation
var_add <- c(40, 300, 50, 1, 50)
var_obs <- c(10, 10, 1, 1^2, 10)
prec_add <- 1 / var_add
prec_obs <- 1 / var_obs
prec_init <- 1 / init_var

## ---------------- MODEL BLOCK -----------------------

## set new timestep
ndays <- rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6)

## set observation index
missing_data_ind = c(seq(from = 1, to = 24*3, by = 1))

## aggregated canopy model for GPP calculations
GPP = function(lai, p11, maxt, mint, Ca, lat, yearday, Nit, rad){
  ## Creates GPP function in order to compute the GPP for a day based on values input from the drivers
  
  ## Sets parameters
  psid = -2
  rtot = 1
  trange = .5*(maxt - mint)
  a = c(p11, 
        0.0156935, 
        4.22273,
        208.868,
        0.0453194, 
        0.37836,
        7.19298,
        0.011136,
        2.1001, 
        0.789798)
  gs = abs((psid))^a[10] / (a[6]*rtot + trange)
  pp = lai*Nit/gs*a[1]*exp(a[8]*maxt)
  qq = a[3] - a[4]
  ci = .5*(Ca + qq - pp + ((Ca+qq-pp)^2 - 
                             4*(Ca*qq - pp*a[3]))^(.5) )
  e0 = (a[7]*lai**2 / (lai**2 + a[9]) )
  dec = -23.4*cos((360*(yearday + 10)/365)*pi/180)*pi/180
  mult = tan(lat)*tan(dec)
  
  ## Changes dayl variable depending on value of mult
  if (mult >= 1){
    dayl = 24
  } else if(mult <= -1){
    dayl = 0
  } else{
    dayl = 24*acos(-mult)/pi
  }
  cps = e0*rad*gs*(Ca - ci) / (e0*rad + gs*(Ca - ci))
  gpp = cps*(a[2]*dayl + a[5])
  return(gpp)
}
## Sets global constants
nday = 730
psid = -2
rtot = 1

## read in driver data
drivers_ev = read.csv('./temp_data.csv', header = T)
drivers_ev$mint = na.interpolation(drivers_ev$mint)

## use missing data interpolation to get a full driver data set
drivers_ev$maxt <- na.interpolation(drivers_ev$maxt)
drivers_ev$rad <- na.interpolation(drivers_ev$rad)
drivers_ev$rad <- drivers_ev$rad * .0864

## Reads in initial conditions
init_ev = read.csv('./EV1_initial.csv', header = F)

## load zoo library to help with creating drivers dataframe
library(zoo)
drivers_ev <- do.call(rbind, replicate(6, coredata(drivers_ev), simplify = FALSE))
drivers_ev$projectday <- as.numeric(rownames(drivers_ev))
nmonths <- 24*3
nday <- 730*3
drivers_ev$ca <- c(rep(396.52, 365), rep(398.65, 365), rep(400.83, 365),
                   rep(404.24, 365), rep(406.55, 365), rep(408.52, 365))
drivers_ev <- drivers_ev[1:nday,]

## Sets values for p and initializes relevant variables
p_e = c(1/500, 0.33, 0.15, 0.27, 1/(2*365), 1/(25*365), 1/(2*365), .1, 1/(250* 365), 0.1725, 3.000000000)

## initialize storage for fluxes and stocks
Cf_e = rep(NA, nmonths)
Cr_e = rep(NA, nmonths)
Cw_e = rep(NA, nmonths)
Clit_e = rep(NA, nmonths)
Csom_e = rep(NA, nmonths)
G_e = rep(NA, nday)
Ra_e = rep(NA, nday)
Af_e = rep(NA, nday)
Aw_e = rep(NA, nday)
Ar_e = rep(NA, nday)
Lf_e = rep(NA, nday)
Lw_e = rep(NA, nday)
Lr_e = rep(NA, nday)
Rh1_e = rep(NA, nday)
Rh2_e = rep(NA, nday)
D_e = rep(NA, nday)
NEE_e = rep(NA, nday)
LAI = rep(NA, nday)
init_ev[12,] <- init_mean[1]
Cf_e[1] = init_ev[12,]
init_ev[13,] <- init_mean[3]
Cr_e[1] = init_ev[13,]
init_ev[14,] <- init_mean[2]
Cw_e[1] = init_ev[14,]
init_ev[15,] <- init_mean[4]
Clit_e[1] = init_ev[15,]
init_ev[16,] <- init_mean[5]
Csom_e[1] = init_ev[16,]
Clit_multiplier <- rep(NA, nday)
Cf_approx <- rep(NA, nday)
Cr_approx <- rep(NA, nday)
Cw_approx <- rep(NA, nday)
Clit_approx <- rep(NA, nday)
Csom_approx <- rep(NA, nday)
Clit_som_multiplier <- rep(NA, nday)
Csom_multiplier <- rep(NA, nday)
## Scales lat, sets Nit and LMA
lat_e = init_ev[17,]
lat_e = lat_e * pi / 180
Nit_e = init_ev[18,]
LMA_e = init_ev[19,]
ind <- 1

## Iterates through forloop to compute carbon fluxes
for (i in 1:nday){
  ## j is a dummy variable that was needed to match the results from the original F90 model, which did not use indices in the loop
  ## and instead wrote to the .csv file at every iteration in the loop
  if (i == 1){
    j = 1
  } else{
    j = i-1
  }
  ## Computes LAI, passes LAI, p11, and initial conditions and drivers into GPP function to calculate G_e
  if (i == 1){
    ## compute fast timestep fluxes and stocks
    Cf_approx[i] <- Cf_e[1]
    Cr_approx[i] <- Cr_e[1]
    Cw_approx[i] <- Cw_e[1]
    Clit_approx[i] <- Clit_e[1]
    Csom_approx[i] <- Csom_e[1]
    LAI[i] = max(.1, Cf_approx[i]/LMA_e)
    G_e[i] = GPP(lai = LAI[i], p11 = p_e[11], maxt = drivers_ev$maxt[i], mint = drivers_ev$mint[i], Ca = drivers_ev$ca[i], 
                 lat = lat_e, yearday = drivers_ev$yearday[i], Nit = Nit_e, rad = drivers_ev$rad[i])
    Trate = .5*exp(p_e[10]*.5*(drivers_ev$maxt[i] + drivers_ev$mint[i]))
    Ra_e[i] = p_e[2]*G_e[i]
    Af_e[i] = (G_e[i] - Ra_e[i])*p_e[3]
    Ar_e[i] = (G_e[i] - Ra_e[i] - Af_e[i])*p_e[4]
    Aw_e[i] = G_e[i] - Ra_e[i] - Af_e[i] - Ar_e[i]
    Lf_e[i] = p_e[5]*Cf_approx[j]
    Lw_e[i] = p_e[6]*Cw_approx[j]
    Lr_e[i] = p_e[7]*Cr_approx[j]
    Rh1_e[i] = p_e[1]*p_e[8]*Clit_approx[j]*Trate
    Rh2_e[i] = p_e[9]*Csom_approx[j]*Trate
    D_e[i] = (p_e[1]-p_e[1]*p_e[8])*Clit_approx[j]*Trate
    NEE_e[i] <- Ra_e[i] + Rh1_e[i] + Rh2_e[i] - G_e[i]
    Cf_approx[i] <- Cf_approx[i] * (1 - p_e[5]) + G_e[i]*(1 - p_e[2])*p_e[3]
    Cr_approx[i] <- Cr_approx[i] * (1 - p_e[7]) + G_e[i] * (1 - p_e[2]) * (1 - p_e[3]) * p_e[4]
    Cw_approx[i] <- Cw_approx[i] * (1 - p_e[6]) + G_e[i] * (1 - p_e[2]) * (1 - p_e[3]) * (1 - p_e[4])
    Clit_approx[i] <- Clit_approx[i] * (1 - .5*exp(p_e[10] * .5 * (drivers_ev$mint[i] + drivers_ev$maxt[i]))*(p_e[1])) +
      p_e[5] * Cf_e[1] + p_e[7] * Cr_e[1]
    Csom_approx[i] <- (1 - .5*exp(p_e[10] * .5 * (drivers_ev$mint[i] + drivers_ev$maxt[i]))*p_e[9]) * Csom_approx[i] +
      .5*exp(.5 * p_e[10] * (drivers_ev$mint[i] + drivers_ev$maxt[i])) * (p_e[1]-p_e[1]*p_e[8]) * Clit_e[1] + p_e[6] * Cw_e[1]
    
  } else{
    ## compute fast timestep fluxes and stocks
    LAI[i] = max(.1, Cf_approx[i-1]/LMA_e)
    G_e[i] = GPP(lai = LAI[i], p11 = p_e[11], maxt = drivers_ev$maxt[i], mint = drivers_ev$mint[i], Ca = drivers_ev$ca[i], 
                 lat = lat_e, yearday = drivers_ev$yearday[i], Nit = Nit_e, rad = drivers_ev$rad[i])
    Trate = .5*exp(p_e[10]*.5*(drivers_ev$maxt[i] + drivers_ev$mint[i]))
    Ra_e[i] = p_e[2]*G_e[i]
    Af_e[i] = (G_e[i] - Ra_e[i])*p_e[3]
    Ar_e[i] = (G_e[i] - Ra_e[i] - Af_e[i])*p_e[4]
    Aw_e[i] = G_e[i] - Ra_e[i] - Af_e[i] - Ar_e[i]
    Lf_e[i] = p_e[5]*Cf_approx[j]
    Lw_e[i] = p_e[6]*Cw_approx[j]
    Lr_e[i] = p_e[7]*Cr_approx[j]
    Rh1_e[i] = p_e[1]*p_e[8]*Clit_approx[j]*Trate
    Rh2_e[i] = p_e[9]*Csom_approx[j]*Trate
    D_e[i] = (p_e[1]-p_e[1]*p_e[8])*Clit_approx[j]*Trate
    NEE_e[i] <- Ra_e[i] + Rh1_e[i] + Rh2_e[i] - G_e[i]
    if ((i-1) %in% cumsum(ndays)){
      tmp <- c(Cf_e[ind], Cr_e[ind], Cw_e[ind], Clit_e[ind], Csom_e[ind])
      tmp2 <- c(Cf_approx[i-1], Cr_approx[i-1], Cw_approx[i-1], Clit_approx[i-1], Csom_approx[i-1])
      Cf_approx[i-1] <- tmp[1]
      Cr_approx[i-1] <- tmp[2]
      Cw_approx[i-1] <- tmp[3]
      Clit_approx[i-1] <- tmp[4]
      Csom_approx[i-1] <- tmp[5]
    }
    Cf_approx[i] <- Cf_approx[i-1] * (1 - p_e[5]) + G_e[i]*(1 - p_e[2])*p_e[3]
    Cr_approx[i] <- Cr_approx[i-1] * (1 - p_e[7]) + G_e[i] * (1 - p_e[2]) * (1 - p_e[3]) * p_e[4]
    Cw_approx[i] <- Cw_approx[i-1] * (1 - p_e[6]) + G_e[i] * (1 - p_e[2]) * (1 - p_e[3]) * (1 - p_e[4])
    Clit_approx[i] <- Clit_approx[i-1] * (1 - .5*exp(p_e[10] * .5 * (drivers_ev$mint[i] + drivers_ev$maxt[i]))*(p_e[1])) +
      p_e[5] * Cf_approx[i-1] + p_e[7] * Cr_approx[i-1]
    Csom_approx[i] <- (1 - .5*exp(p_e[10] * .5 * (drivers_ev$mint[i] + drivers_ev$maxt[i]))*p_e[9]) * Csom_approx[i-1] +
      .5*exp(.5 * p_e[10] * (drivers_ev$mint[i] + drivers_ev$maxt[i])) * (p_e[1]-p_e[1]*p_e[8]) * Clit_approx[i-1] + p_e[6] * Cw_approx[i-1]
    
    
    if (i %in% cumsum(ndays)){
      ## compute slow timestep stocks
      ind <- which(cumsum(ndays) == i)
      Cf_e[ind] <- Cf_approx[i] + rnorm(1, 0, sqrt(var_add[1]))
      Cw_e[ind] <- Cw_approx[i] + rnorm(1, 0, sqrt(var_add[2]))
      Cr_e[ind] <- Cr_approx[i] + rnorm(1, 0, sqrt(var_add[3]))
      Clit_e[ind] <- Clit_approx[i] + rnorm(1, 0, sqrt(var_add[4]))
      Csom_e[ind] <- Csom_approx[i] + rnorm(1, 0, sqrt(var_add[5]))
    }
    
  }
  if ((i-1) %in% cumsum(ndays)){
    ## compute slow timestep stocks
    Cf_approx[i-1] <- tmp2[1]
    Cr_approx[i-1] <- tmp2[2]
    Cw_approx[i-1] <- tmp2[3]
    Clit_approx[i-1] <- tmp2[4]
    Csom_approx[i-1] <- tmp2[5]
  }
}
## ----- SYNTHESIZE DATA ---------------

## Adds observational noise
Cf.obs = Cf_e + rnorm(nmonths, 0, sqrt(var_obs[1]))
Cw.obs = Cw_e + rnorm(nmonths, 0, sqrt(var_obs[2]))
Cr.obs = Cr_e + rnorm(nmonths, 0, sqrt(var_obs[3]))
Clit.obs= Clit_e + rnorm(nmonths, 0, sqrt(var_obs[4]))
Csom.obs = Csom_e + rnorm(nmonths, 0, sqrt(var_obs[5]))
G_obs = G_e + rnorm(nday, 0, sd = .001)
Lr_obs <- Lr_e + rnorm(nday, 0, .001)

## add observational noise to d_e flux and create d_e likelihood function
D_e_obs <- D_e + rnorm(730, 0, .001)
.5*exp(p_e[10]*.5*(drivers_ev$maxt[i] + drivers_ev$mint[i]))
De_ll <- function(p, C, maxt, mint, D_e_obs){
  dnorm(D_e_obs[1], (p[1]-p[1]*p[8])*C[1,4]*.5*exp(p[10]*.5*(maxt[1] + mint[1])), sd = .001, log = TRUE)
  + sum(dnorm(D_e_obs[2:2190], (p[1]-p[1]*p[8])*C[(1:2189),4]*.5*exp(p[10]*.5*(maxt[2:2190] + mint[2:2190])), sd = .001, log = TRUE))
}
## add observational noise to ra flux and create likelihood function
Ra_obs <- Ra_e + rnorm(nday, 0, .001)

Ra_ll <- function(p, G, Ra_obs){
  sum(dnorm(p[2]*G, Ra_obs, .001, log = TRUE) )
}

## add observational noise to rh1 flux and create likelihood function
Rh1_obs <- Rh1_e + rnorm(nday, 0, .0001)

Rh1_ll <- function(p, C, Rh1_obs, maxt, mint){
  sum(dnorm(p[1]*p[8]*C[(1:2189),4]*.5*exp(p[10]*.5*(maxt[2:2190] + mint[2:2190])), Rh1_obs[2:2190], .0001, log = TRUE) ) +
    dnorm(p[1]*p[8]*C[(1),4]*.5*exp(p[10]*.5*(maxt[1] + mint[1])), log = TRUE)
}

## add observational noise to ar flux and create likelihood function
Ar_obs <- Ar_e + rnorm(nday, 0, .001)
Ar_ll <- function(p, G, Ar_obs){
  sum(dnorm(p[4]*(1-p[3])*(1-p[2])*G, Ar_obs, .001, log = TRUE) )
}

## add observational noise to rh2 flux and create likelihood function
Rh2_obs <- Rh2_e + rnorm(nday, 0, .0005)
Rh2_ll <- function(p, C, Rh2_obs, maxt, mint){
  sum(dnorm(p[9]*C[(1:2189),5]*.5*exp(p[10]*.5*(maxt[2:2190] + mint[2:2190])), Rh2_obs[2:2190], .0005, log = TRUE) ) +
    dnorm(p[9]*C[(1),5]*.5*exp(p[10]*.5*(maxt[1] + mint[1])), Rh2_obs[1], .0005, log = TRUE)
}

## add observational noise to af flux and create likelihood function
Af_obs <- Af_e + rnorm(nday, 0, .01)

Af_ll <- function(p, Af_obs, G){
  sum(dnorm(G*(p[3])*(1-p[2]), mean = Af_obs, sd = .01, log = TRUE))
}

## add observational noise to lf flux and create likelihood function
Lf_obs <- Lf_e + rnorm(nday, 0, .005)
Lf_ll <- function(p, C, Lf_obs){
  sum(dnorm(Lf_obs[2:2190], p[5]*C[(1:2189),1], sd = .005, log = TRUE)) +
    dnorm(Lf_obs[1], p[5]*C[(1),1], sd = .05, log = TRUE)
}

## add observational noise to nee flux
NEE_obs <- NEE_e + rnorm(nday, 0, .5)

## add observational noise to lw flux and create likelihood function
Lw_obs <- Lw_e + rnorm(nday, 0, .005)
Lw_ll <- function(p, Lw_obs, C){
  sum(dnorm(Lw_obs[2:2190], p[6]*C[(1:2189),2], .005, log = TRUE)) +
    dnorm(Lw_obs[1], p[6]*C[1, 2], .005, log = TRUE)
}

## add observational noise to lr flux and create likelihood function
Lr_ll <- function(p, C, Lr_obs){
  Lr_ll <- sum(dnorm(Lr_obs[2:2190], p[7]*C[(1:2189),3], sd = .001, log = TRUE)) +
    dnorm(Lr_obs[1], p[7]*C[1, 3], sd = .001, log = TRUE)
  return(Lr_ll)
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
chain_length <- 10000
burn <- 2000

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
C_samples[, , 1] <- latent_initialization(C_obs, missing_data_ind, 1, nmonths, type = c('pwl', 'pwl', 'pwl', 'pwl', 'pwl'))
pars <- matrix(NA, nrow = 1, ncol = 16)
pars[1,] <- c(.5*(pupper[1:11] + plower[1:11]), rep(1,5))
pars[1, 12:16] <- rep(1, 5)
## initialize acceptance matrix
acceptance <- matrix(0, nrow = chain_length, ncol = 16)

## initialize proposal stepsize
sd <- rep(.3, 11)
sd[11] <- .02
sd[1] <- .01
sd[2] <- .002
sd[3] <- .002
sd[4] <- .001
sd[5] <- .03
sd[6] <- .3
sd[7] <- .008
sd[8] <- .03
sd[9] <- .01
sd[10] <- .008
sd[11] = .009
sd[12:16] <- .17
colnames(pars) <- c(paste0('p', 1:11), paste0('precadd', 1:5))
preds <- getDALECpredictions2(p_e = pars[1:11], C = C_samples[,,1], LAI = Cf_approx / LMA_e, maxt = maxt, mint = mint,
                              Ca = Ca, lat_e = lat_e, yearday = yearday, Nit_e = Nit_e, rad = rad, N = nday, G_obs = G_e, ind = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6))-1, 
                              init_mean = C.0[1,], native_step = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6)))


old <- ss_monthly_ll(C_obs = C_obs, C = C_samples[,,1], init_mean = init_mean, init_var = init_var, prec_obs = prec_obs, pars = pars, fluxes = fluxes_obs,
                     G = G[1,], drivers = drivers_ev, init_ev = init_ev, ndays = ndays, Ca = Ca, lat_e = lat_e, yearday = yearday, rad = rad, N = nday, Nit_e = Nit_e, r = r)

trueLL <- ss_monthly_ll(C_obs = C_obs, C = C, init_mean = init_mean, init_var = init_var, prec_obs = prec_obs, pars = c(p_e, prec_add), fluxes = fluxes_obs,
                        G = G_e, drivers = drivers_ev, init_ev = init_ev, ndays = ndays, Ca = Ca, lat_e = lat_e, yearday = yearday, rad = rad, N = nday, Nit_e = Nit_e, r = r)

old <- old

accept <- rep(0, 16)
cov.block <- NULL
thresh2 <- matrix(0, ncol = 16, nrow = chain_length)
p_star_mat <- matrix(0, ncol = 16, nrow = chain_length)
update = TRUE
block <- vector('list', 2)
block[[1]] <- c(2,3,4)
block[[2]] <- c(9, 10)
#block[[3]] <- c(8, 10)

preds <- getDALECpredictions2(p_e = p_e, C = C, LAI = Cf_approx / LMA_e, maxt = maxt, mint = mint,
                              Ca = Ca, lat_e = lat_e, yearday = yearday, Nit_e = Nit_e, rad = rad, N = 2190, G_obs = G_e, 
                              ind = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6))-1, 
                              init_mean = init_mean, native_step = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6)))

test <- cbind(Cf_approx, Cr_approx, Cw_approx, Clit_approx, Csom_approx)
Cpredi <- preds$Cslow

## Set up MCMC 
sd <- sd
sd[c(2,3,4)] <- sd[c(2,3,4)] * 5
sd <- sd 
tf <- rep('log', 16)
old <- old * r
