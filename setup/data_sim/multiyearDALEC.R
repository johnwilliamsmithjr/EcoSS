## MultiYearDALEC.R
## Replicate dataframe 10 times and run DALEC

## Description: Simulates data using the Reflex Evergreen model. Adds both process and observational noise into the data.

## ---------------- MODEL BLOCK -----------------------
## set.seed(777) 

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

psid = -2
rtot = 1

## Reads in initial conditions
if ('EV1_initial.csv' %in% dir('./data/')){
  init_ev = read.csv('./data/EV1_initial.csv', header = F)
} else{
  cat('Please add file EV1_initial.csv to the directory. \n')
}

## quick check on projectday

drivers_ev$projectday <- as.numeric(rownames(drivers_ev))

## allocates space for fluxes
Cf_e = rep(NA, nday)
Cr_e = rep(NA, nday)
Cw_e = rep(NA, nday)
Clit_e = rep(NA, nday)
Csom_e = rep(NA, nday)
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
Cf_e[1] = init_mean[1]
Cr_e[1] = init_mean[3]
Cw_e[1] = init_mean[2]
Clit_e[1] = init_mean[4]
Csom_e[1] = init_mean[5]

## Scales lat, sets Nit and LMA
lat_e = init_ev[17,]
lat_e = lat_e * pi / 180
Nit_e = init_ev[18,]
LMA_e = init_ev[19,]

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
  LAI[i] = max(.1, Cf_e[j]/LMA_e)
  G_e[i] = GPP(lai = LAI[i], p11 = p_e[11], maxt = drivers_ev$maxt[i], mint = drivers_ev$mint[i], Ca = drivers_ev$ca[i], 
               lat = lat_e, yearday = drivers_ev$yearday[i], Nit = Nit_e, rad = drivers_ev$rad[i])
  ## Calculates Trate
  Trate = .5*exp(p_e[10]*.5*(drivers_ev$maxt[i] + drivers_ev$mint[i]))
  ## Calculates values needed to compute carbon pools
  Ra_e[i] = p_e[2]*G_e[i]
  Af_e[i] = (G_e[i] - Ra_e[i])*p_e[3]
  Ar_e[i] = (G_e[i] - Ra_e[i] - Af_e[i])*p_e[4]
  Aw_e[i] = G_e[i] - Ra_e[i] - Af_e[i] - Ar_e[i]
  Lf_e[i] = p_e[5]*Cf_e[j]
  Lw_e[i] = p_e[6]*Cw_e[j]
  Lr_e[i] = p_e[7]*Cr_e[j]
  Rh1_e[i] = p_e[8]*Clit_e[j]*Trate
  Rh2_e[i] = p_e[9]*Csom_e[j]*Trate
  D_e[i] = p_e[1]*Clit_e[j]*Trate
  ## Compute Carbon Pools, add process error noise
  Cf_e[i] = Cf_e[j] + Af_e[i] - Lf_e[i] + rnorm(1, 0, sqrt(var_add[1]))
  Cw_e[i] = Cw_e[j] + Aw_e[i] - Lw_e[i] + rnorm(1, 0, sqrt(var_add[2]))
  Cr_e[i] = Cr_e[j] + Ar_e[i] - Lr_e[i] + rnorm(1, 0, sqrt(var_add[3]))
  Clit_e[i] = Clit_e[j] + Lf_e[i] + Lr_e[i] - Rh1_e[i] - D_e[i] + rnorm(1, 0, sqrt(var_add[4]))
  Csom_e[i] = Csom_e[j] + D_e[i] - Rh2_e[i] + Lw_e[i] + rnorm(1, 0, sqrt(var_add[5]))
  NEE_e[i] = Ra_e[i] + Rh1_e[i] + Rh2_e[i] - G_e[i]
}
## NEE, Lf_e, Lw_e, 
## Stocks annually
## 
## -------------------------------------

## ----- SYNTHESIZE DATA ---------------

## Adds observational noise
Cf.obs = Cf_e + rnorm(nday, 0, sqrt(var_obs[1]))
Cw.obs = Cw_e + rnorm(nday, 0, sqrt(var_obs[2]))
Cr.obs = Cr_e + rnorm(nday, 0, sqrt(var_obs[3]))
Clit.obs= Clit_e + rnorm(nday, 0, sqrt(var_obs[4]))
Csom.obs = Csom_e + rnorm(nday, 0, sqrt(var_obs[5]))
G_obs = G_e + rnorm(nday, 0, sd = .2)
Lr_obs <- Lr_e + rnorm(nday, 0, .001)

D_e_obs <- D_e + rnorm(nday, 0, .01)
De_ll <- function(p, C, maxt, mint, D_e_obs){
  sum(dnorm(D_e_obs, p[1]*C[,4]*.5*exp(p[10]*.5*(maxt + mint)), sd = .01, log = TRUE))
}

Ra_obs <- Ra_e + rnorm(nday, 0, .08)

Ra_ll <- function(p, G, Ra_obs){
  sum(dnorm(p[2]*G, Ra_obs, .08, log = TRUE) )
}
Rh1_obs <- Rh1_e + rnorm(nday, 0, .2)

Rh1_ll <- function(p, C, Rh1_obs, maxt, mint){
  sum(dnorm(p[8]*C[,4]*.5*exp(p[10]*.5*(maxt + mint)), Rh1_obs, .2, log = TRUE) )
}

Ar_obs <- Ar_e + rnorm(nday, 0, .0005)
Ar_ll <- function(p, G, Ar_obs){
  sum(dnorm(p[4]*(1-p[3])*(1-p[2])*G, Ar_obs, .005, log = TRUE) )
}

#Lr_ll <- function(p, C, maxt, mint){
#  sum(dnorm(D_e_obs, p[1]*C[,4]*.5*exp(p[10]*.5*(maxt + mint)), sd = .1, log = TRUE))
#}

Rh2_obs <- Rh2_e + rnorm(nday, 0, .01)
Rh2_ll <- function(p, C, Rh2_obs, maxt, mint){
  sum(dnorm(p[9]*C[,5]*.5*exp(p[10]*.5*(maxt + mint)), Rh2_obs, .01, log = TRUE) )
}

Af_obs <- Af_e + rnorm(nday, 0, .1)

Af_ll <- function(p, Af_obs, G){
  sum(dnorm(G*(p[3])*(1-p[2]), mean = Af_obs, sd = .1, log = TRUE))
}

Lf_obs <- Lf_e + rnorm(nday, 0, .2)
Lf_ll <- function(p, C, Lf_obs){
  sum(dnorm(Lf_obs, p[5]*C[,1], sd = .2, log = TRUE))
}

Lw_obs <- Lw_e + rnorm(nday, 0, .005)
Lw_ll <- function(p, Lw_obs, C){
  sum(dnorm(Lw_obs, p[6]*C[,2], .005, log = TRUE))
}
