## Program: monthly_timestep_yearly_obs_datasim.R
## Date: 10/19/2018

## load libraries
library(Rcpp, verbose = FALSE)
library(coda, verbose = FALSE)
library(truncnorm, verbose = FALSE)
library(tmvtnorm, verbose = FALSE)
library(matrixStats, verbose = FALSE)
library(laGP, verbose = FALSE)
library(imputeTS, verbose = FALSE)
## load lhs library for initial process parameter guesses
library(lhs, verbose = FALSE)
## load HDInterval to obtain credible intervals
library(HDInterval, verbose = FALSE)

if (!('p_e' %in% ls())){
  p_e = c(1/500, 0.33, 0.15, 0.27, 1/(2*365), 1/(25*365), 1/(2*365), .1, 1/(250* 365), 0.1725, 3.000000000)
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
  var_add <- c(40, 300, 50, 1, 50)
  cat('No changes detected in process precisions, using default precision set\n')
}

if (!('var_obs' %in% ls())){
  var_obs <- c(10, 10, 1, 1^2, 10)
  cat('No changes detected in observation precisions, using default precision set\n')
}

## set observation and additive (process) variance for simulation
#var_add <- c(40, 300, 50, 1, 50)
#var_obs <- c(10, 10, 1, 1^2, 10)
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
drivers_ev = read.csv('./data/temp_data.csv', header = T)
drivers_ev$mint = na_interpolation(drivers_ev$mint)

## use missing data interpolation to get a full driver data set
drivers_ev$maxt <- na_interpolation(drivers_ev$maxt)
drivers_ev$rad <- na_interpolation(drivers_ev$rad)
drivers_ev$rad <- drivers_ev$rad * .0864

## Reads in initial conditions
init_ev = read.csv('./data/EV1_initial.csv', header = F)

## load zoo library to help with creating drivers dataframe
library(zoo, quietly = TRUE, verbose = FALSE)
drivers_ev <- do.call(rbind, replicate(6, coredata(drivers_ev), simplify = FALSE))
drivers_ev$projectday <- as.numeric(rownames(drivers_ev))
nmonths <- 24*3
nday <- 730*3
drivers_ev$ca <- c(rep(396.52, 365), rep(398.65, 365), rep(400.83, 365),
                   rep(404.24, 365), rep(406.55, 365), rep(408.52, 365))
drivers_ev <- drivers_ev[1:nday,]

## initializes relevant variables

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
