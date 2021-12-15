## updateMCMC function
## Author: John W Smith
## Function to update the state of the parameters in a markov chain monte carlo

## Inputs:
## pars: Argument of class matrix (or forced to be class matrix) containing previous draws of the MCMC parameters
## burn: Number of initial iterations to be "thrown out" while the algorithm burns in
## lh_fun: Function being evaluated to compute the model likelihood
## blocks: parameters to be sampled jointly 
## i: current iteration
## upper, lower: respective bounds for parameters
## transform: transformation to be applied
## adapt.length: length after burn in to wait to find the covariance structure between parameters being sampled jointly
## adapt.stop: length at which to stop computing the empirical covariance
## sd.proposal: standard deviation of proposal distributions
## lh_old: likelihood of previous iterations
## acceptance: matrix of acceptances for parameters
## thin: Option to thin out samples to help with memory and autocorrelation
## block.cov: covariance of parameters being sampled jointly
## check.accept: point at which to start adapting standard deviations for proposals

updateMCMC <- function(pars, parnames = colnames(pars), burn, lh_fun, block = NULL, i, upper, lower, transform = rep('none', length(parnames)),
         adapt.length = 500, adapt.stop = 50000, sd.proposal, lh_old, acceptance, thin = 4, block.cov = NULL, check.accept = 3000, r = 1,
         cov.adjust = rep(1, length(block)), method = "gibbs",  ...){
  ## Ensure that par is of class matrix
  if (!(is.null(dim(pars)))){
    new_par <- pars[i,]
  } else{
    new_par <- pars
  }
  if (is.null(parnames)) parnames <- names(pars)
  if (is.null(parnames)) parnames <- rep(0, length(upper))
  new_accept <- matrix(0, ncol = length(parnames), nrow = thin)
  ## Initializes storage for log threshold values and proposed points
  thresh <- rep(0,length(parnames))
  proposal <- rep(0,length(parnames))
  for (k in 1:thin){  
    ## Repeats loop for thin interations
    if (i < burn + adapt.length){
      ## Checks to see if adaptive sampling is to performed
      for (j in 1:length(parnames)){
        ## Loops over MCMC parameters
        if (i >= check.accept & i%%500 == 0){
          ## Tunes proposal standard deviations if wanted
          if (mean(acceptance[(i-499):i,j]) < .2){
            sd.proposal[j] <- sd.proposal[j] / 1.1
          } 
          if (mean(acceptance[(i-499):i,j]) > .4){
            sd.proposal[j] <- sd.proposal[j] * 1.1
          } 
        }
        ## Samples for log transformed variables
        if (transform[j] == 'log'){
          ## check if upper or lower bounds are NA
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            ## propose point
            p_star <- new_par
            tmp <- rtruncnorm(1, a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])
            ## compute proposal correction term
            prop.correction <- log(dtruncnorm(log(new_par[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm(log(p_star[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(new_par[j]), sd = sd.proposal[j]))
          } else{
            ## if neither upper or lower are NA
            ## propose point
            p_star <- new_par
            tmp <- rnorm(1, mean = log(p_star[j]), sd = sd.proposal[j])
            prop.correction <- 0
          }
          ## Transform proposed point back
          p_star[j] <- exp(tmp)
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          ## acceptance step
          #print(j)
          #print(lh_new)
          #print(prop.correction)
          #print(lh_old)
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
        } 
        if (transform[j] == 'none'){
          ## check if upper or lower bounds are NA
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            ## propose point
            p_star <- new_par
            tmp <- rtruncnorm(1, a = plower[j], b = pupper[j], mean = p_star[j], sd = sd.proposal[j])
            ## compute proposal correction term
            prop.correction <- log(dtruncnorm(new_par[j], a = plower[j], b = pupper[j], mean = p_star[j], sd = sd.proposal[j])) - 
              log(dtruncnorm(p_star[j], a = plower[j], b = pupper[j], mean = new_par[j], sd = sd.proposal[j]))
          } else{
            ## if neither upper or lower are NA
            ## propose point
            p_star <- new_par
            tmp <- rnorm(1, mean = p_star[j], sd = sd.proposal[j])
            prop.correction <- 0
          }
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          ## acceptance step
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
        } 
        ## store threshold value and proposed point
        thresh[j] <- threshold
        proposal[j] <- exp(tmp)
      }
    }
    ## check to see if current iteration is large enough to adapt
    if (i >= burn + adapt.length){
      if (i %% 500 == 0 & i < adapt.stop){
          ## create space for block.cov 
          block.cov <- vector('list', length(block))
        for (b in 1:length(block)){
          ## iterate through list of parameters to be block sampled and compute empirical correlation structures
          block.cov[[b]] <- cov(matrix(pars[burn:i, (block[[b]]) ], ncol = length(block[[b]])))
          block.cov[[b]] <- block.cov[[b]] / cov.adjust[b]
          
        }
      }
      for (j in (1:length(parnames))[-unlist(block)]){
        if (transform[j] == 'log'){
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            p_star <- new_par
            tmp <- rtruncnorm(1, a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])
            prop.correction <- log(dtruncnorm(log(new_par[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm(log(p_star[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(new_par[j]), sd = sd.proposal[j]))
          } else{
            p_star <- new_par
            tmp <- rnorm(1, mean = log(p_star[j]), sd = sd[j])
            prop.correction <- 0
          }
          p_star[j] <- exp(tmp)
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
        } 
        if (transform[j] == 'none'){
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            p_star <- new_par
            tmp <- rtruncnorm(1, a = plower[j], b = pupper[j], mean = p_star[j], sd = sd.proposal[j])
            prop.correction <- log(dtruncnorm(new_par[j], a = plower[j], b = pupper[j], mean = p_star[j], sd = sd.proposal[j])) - 
              log(dtruncnorm(p_star[j], a = plower[j], b = pupper[j], mean = new_par[j], sd = sd.proposal[j]))
          } else{
            p_star <- new_par
            tmp <- rnorm(1, mean = p_star[j], sd = sd.proposal[j])
            prop.correction <- 0
          }
          p_star[j] <- tmp
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
        } 
      }
      for (b in 1:length(block)){
        blocks <- block[[b]]
        p_star <- new_par
      if (all(is.na(pupper[blocks]))){
        if (method == 'gibbs'){
          start.value = pars[i, blocks]
          burn.in.samples = 500
          tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], sigma = block.cov[[b]], algorithm = method, start.value = start.value, 
                          burn.in.samples = burn.in.samples)
        } else{
            tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], sigma = block.cov[[b]], algorithm = method)
        }
      } else{
        if (method == 'gibbs'){
          start.value = pars[i, blocks]
          burn.in.samples = 500
          tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]], algorithm = method, start.value = start.value, 
                          burn.in.samples = burn.in.samples)
        } else{
          tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]], algorithm = method)
        }
      }
        #print(tmp)
      p_star[blocks] <- tmp
      lh_new <- lh_fun(pars = p_star, ...)
      if (any(is.na(pupper[blocks]))){
        prop.correction <- log(dtmvnorm(new_par[blocks], lower = plower[blocks], mean = p_star[blocks], sigma = block.cov[[b]])) - 
          log(dtmvnorm(x = p_star[blocks], mean = new_par[blocks], lower = plower[blocks], sigma = block.cov[[b]]))
      } else{
        if (is.nan(log(dtmvnorm(new_par[blocks], lower = plower[blocks], upper = pupper[blocks], mean = p_star[blocks], sigma = block.cov[[b]])) - 
                   log(dtmvnorm(x = p_star[blocks], mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]])))){
          prop.correction <- 0
        } else{
          prop.correction <- log(dtmvnorm(new_par[blocks], lower = plower[blocks], upper = pupper[blocks], mean = p_star[blocks], sigma = block.cov[[b]])) - 
            log(dtmvnorm(x = p_star[blocks], mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]])) 
        }
      }
      #print(prop.correction)
      threshold <- r*(lh_new - lh_old) + prop.correction
      if (log(runif(1,0,1)) < threshold){
        lh_old <- lh_new
        new_accept[k,blocks] <- 1
        new_par <- p_star
      }
      }
    }
  }
  pars <- rbind(pars, new_par)
  acceptance <- rbind(acceptance, colMeans(new_accept))
  return(list(pars = pars, acceptance = acceptance, lh_old = lh_old, cov = block.cov, sd = sd.proposal, thresh = thresh, p_star = proposal))
}

## cpp_ss_likeihood_ex function
## Author: John W Smith Jr

## Inputs:
## Cobs: matrix of observed carbon stocks. If no data is observed at a given timestep, values should be set to NA
## C: matrix of latent carbon stocks. Dimension should match that of Cobs. No entries should be missing
## var_obs: vector of observational variances
## pars: vector of process parameters, expected in the following order: 1:11 - process parameters. 12:16 - process precisions
## init_mean: vector of initial means for the carbon stocks
## init_var: vector of initial standard deviations for the carbon stocks
## N: number of days in simulation
## G: latent GPP - currently defunct and not used.
## G_e: observed GPP
## drivers: a list with required drivers for the DALEC model
## LAI: Leaf area index, must be the same row dimension as C and Cobs
## fluxes: a list with necessary fluxes. cpp_ss_likelihood uses every flux in the DALEC model.
cpp_ss_likelihood_ex = function(Cobs, C, var_obs, pars, init_mean, init_sd, N, G, G_e, drivers, LAI, fluxes, ind){
  ## set up parameters in form that C++ function wants
  Sigma <- var_obs
  sd <- 1 / pars[12:16]
  p <- pars[1:11]
  ## Get driver data from driver list
  Ca = drivers$ca[1:N]
  maxt = drivers$maxt[1:N]
  mint = drivers$mint[1:N]
  rad = drivers$rad[1:N]
  yearday = drivers$yearday[1:N]
  ## compute indices for missing data
  if (any(is.na(Cobs[,1]))){
    ind <- c(1:(dim(Cobs)[1]))[-which(is.na(Cobs[,1]))]
  } else{
    ind <- c(1:(dim(Cobs)[1]))
  }
  ind <- ind - 1
  missing_data_ind <- ind + 1
  ## build cpred
  buildCpredi <- Cpred(p_ = p, C_ = C, LAI_ = LAI, max_t_ = maxt, min_t_ = mint, Ca_ = Ca, lat_e_ = lat_e, yearday_ = yearday,
                       Nit_e_ = Nit_e, rad_ = rad, N = dim(Cobs)[1], G_ = G)
  Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                  buildCpredi$Clitpred, buildCpredi$Csompred)
  G_pred <- buildCpredi$G
  LL = 0
  for (i in 1:5){
    LL = LL + lgammad(sd[i], 10, 1)
  }
  ## Compute likelihood
  LL = LL + SSLL(Cobs_ <- Cobs, C_ <- C, Sigma_ <- Sigma, Cpred_ <- Cpredi, sd_ <- sd, init_mean_ <- init_mean, 
                 init_sd_ <- init_sd, N <- N, G_ <- G, G_e, rep(.025, length(G_e)), g_tau_ <- rep(.2^2, length(G_e)),
                 Gpred_ <- G_pred, ind = ind)
  LL = LL + Lr_ll(p = p, C = C, fluxes$Lr_obs)
  LL = LL + De_ll(p = p, C = C, maxt = maxt, mint = mint,
                  D_e_obs = fluxes$D_e_obs)
  LL = LL + Rh1_ll(p = p, C = C, Rh1_obs = fluxes$Rh1_obs, 
                   maxt = maxt, mint = mint)
  LL = LL + Rh2_ll(p = p, C = C, Rh2_obs = fluxes$Rh2_obs, 
                   maxt = maxt, mint = mint)
  LL = LL + Ar_ll(p = p, G = G, Ar_obs = fluxes$Ar_obs)
  LL = LL + Af_ll(p = p, G = G, Af_obs = fluxes$Af_obs)
  LL = LL + Lf_ll(p = p, C = C, Lf_obs = fluxes$Lf_obs)
  LL = LL + Ra_ll(p = p, G = G, Ra_obs = fluxes$Ra_obs)
  LL = LL + Lw_ll(p = p, Lw_obs = fluxes$Lw_obs, C = C)
  return(LL)
}

cpp_ss_likelihood_ex_bu = function(Cobs, C, var_obs, pars, init_mean, init_sd, N, G, G_e, drivers, LAI, fluxes){
  Sigma <- var_obs
  sd <- 1 / pars[12:16]
  p <- pars[1:11]
  Ca = drivers$ca[1:N]
  maxt = drivers$maxt[1:N]
  mint = drivers$mint[1:N]
  rad = drivers$rad[1:N]
  yearday = drivers$yearday[1:N]
  if (any(is.na(Cobs[,1]))){
    ind <- c(1:(dim(Cobs)[1]))[-which(is.na(Cobs[,1]))]
  } else{
    ind <- c(1:(dim(Cobs)[1]))
  }
  ind <- ind - 1
  missing_data_ind <- ind + 1
  buildCpredi <- Cpred(p_ = p, C_ = C, LAI_ = LAI, max_t_ = maxt, min_t_ = mint, Ca_ = Ca, lat_e_ = lat_e, yearday_ = yearday,
                       Nit_e_ = Nit_e, rad_ = rad, N = dim(Cobs)[1], G_ = G)
  Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                  buildCpredi$Clitpred, buildCpredi$Csompred)
  G_pred <- buildCpredi$G
  LL = 0
  for (i in 1:5){
    LL = LL + lgammad(sd[i], 10, 1)
  }
  #print(LL)
  LL = LL + SSLL(Cobs_ <- Cobs, C_ <- C, Sigma_ <- Sigma, Cpred_ <- Cpredi, sd_ <- sd, init_mean_ <- init_mean, 
                 init_sd_ <- init_sd, N <- N, G_ <- G, G_e, rep(.025, length(G_e)), g_tau_ <- rep(.2^2, length(G_e)),
                 Gpred_ <- G_pred, ind = ind)
  #print(LL)
  #LL = LL + dmvnorm(G, G_e, sigma = diag(.025, length(G_e)), log = T)
  LL = LL + Lr_ll(p = p, C = C[missing_data_ind,], fluxes$Lr_obs[missing_data_ind])
  #print(LL)
  LL = LL + De_ll(p = p, C = C[missing_data_ind,], maxt = maxt[missing_data_ind], mint = mint[missing_data_ind],
                  D_e_obs = fluxes$D_e_obs[missing_data_ind])
  #print(LL)
  LL = LL + Rh1_ll(p = p, C = C[missing_data_ind,], Rh1_obs = fluxes$Rh1_obs[missing_data_ind], 
                   maxt = maxt[missing_data_ind], mint = mint[missing_data_ind])
  #print(LL)
  LL = LL + Rh2_ll(p = p, C = C[missing_data_ind,], Rh2_obs = fluxes$Rh2_obs[missing_data_ind], 
                   maxt = maxt[missing_data_ind], mint = mint[missing_data_ind])
  #print(LL)
  LL = LL + Ar_ll(p = p, G = G[missing_data_ind], Ar_obs = fluxes$Ar_obs[missing_data_ind])
  #print(LL)
  LL = LL + Af_ll(p = p, G = G[missing_data_ind], Af_obs = fluxes$Af_obs[missing_data_ind])
  #print(LL)
  LL = LL + Lf_ll(p = p, C = C[missing_data_ind,], Lf_obs = fluxes$Lf_obs[missing_data_ind])
  #print(LL)
  LL = LL + Ra_ll(p = p, G = G[missing_data_ind], Ra_obs = fluxes$Ra_obs[missing_data_ind])
  #print(LL)
  
  #return(list(LL = LL, gp = G_pred, cp = Cpredi))
  return(LL)
}

cpp_ss_likelihood_g4 = function(Cobs, C, var_obs, pars, init_mean, init_sd, N, G, G_e, drivers, LAI, fluxes){
  Sigma <- var_obs
  sd <- 1 / pars[12:16]
  p <- pars[1:11]
  Ca = drivers$ca[1:N]
  maxt = drivers$maxt[1:N]
  mint = drivers$mint[1:N]
  rad = drivers$rad[1:N]
  yearday = drivers$yearday[1:N]
  if (any(is.na(Cobs[,1]))){
    ind <- c(1:(dim(Cobs)[1]))[-which(is.na(Cobs[,1]))]
  } else{
    ind <- c(1:(dim(Cobs)[1]))
  }
  ind <- ind - 1
  missing_data_ind <- ind + 1
  buildCpredi <- Cpred(p_ = p, C_ = C, LAI_ = LAI, max_t_ = maxt, min_t_ = mint, Ca_ = Ca, lat_e_ = lat_e, yearday_ = yearday,
                       Nit_e_ = Nit_e, rad_ = rad, N = dim(Cobs)[1], G_ = G)
  Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                  buildCpredi$Clitpred, buildCpredi$Csompred)
  G_pred <- buildCpredi$G
  LL = 0
  for (i in 1:5){
    LL = LL + lgammad(sd[i], 10, 1)
  }
  #print(LL)
  LL = LL + SSLL(Cobs_ <- Cobs, C_ <- C, Sigma_ <- Sigma, Cpred_ <- Cpredi, sd_ <- sd, init_mean_ <- init_mean, 
                 init_sd_ <- init_sd, N <- N, G_ <- G, G_e, rep(.025, length(G_e)), g_tau_ <- rep(.2^2, length(G_e)),
                 Gpred_ <- G_pred, ind = ind)
  #print(LL)
  #LL = LL + dmvnorm(G, G_e, sigma = diag(.025, length(G_e)), log = T)
  LL = LL + Lr_ll(p = p, C = C[missing_data_ind,], fluxes$Lr_obs[missing_data_ind])
  #print(LL)
  LL = LL + De_ll(p = p, C = C[missing_data_ind,], maxt = maxt[missing_data_ind], mint = mint[missing_data_ind],
                  D_e_obs = fluxes$D_e_obs[missing_data_ind])
  #print(LL)
  LL = LL + Rh1_ll(p = p, C = C[missing_data_ind,], Rh1_obs = fluxes$Rh1_obs[missing_data_ind], 
                   maxt = maxt[missing_data_ind], mint = mint[missing_data_ind])
  #print(LL)
  LL = LL + Rh2_ll(p = p, C = C[missing_data_ind,], Rh2_obs = fluxes$Rh2_obs[missing_data_ind], 
                   maxt = maxt[missing_data_ind], mint = mint[missing_data_ind])
  #print(LL)
  LL = LL + Ar_ll(p = p, G = G[missing_data_ind], Ar_obs = fluxes$Ar_obs[missing_data_ind])
  #print(LL)
  LL = LL + Af_ll(p = p, G = G[missing_data_ind], Af_obs = fluxes$Af_obs[missing_data_ind])
  #print(LL)
  LL = LL + Lf_ll(p = p, C = C[missing_data_ind,], Lf_obs = fluxes$Lf_obs[missing_data_ind])
  #print(LL)
  LL = LL + Ra_ll(p = p, G = G[missing_data_ind], Ra_obs = fluxes$Ra_obs[missing_data_ind])
  #print(LL)
  LL = LL + Lw_ll(p = p, C = )
  #return(list(LL = LL, gp = G_pred, cp = Cpredi))
  return(LL)
}

cpp_ss_likelihood_nf = function(Cobs, C, var_obs, pars, init_mean, init_sd, N, G, G_e, drivers, LAI, fluxes = NULL){
  Sigma <- var_obs
  sd <- 1 / pars[12:16]
  p <- pars[1:11]
  Ca = drivers$ca[1:N]
  maxt = drivers$maxt[1:N]
  mint = drivers$mint[1:N]
  rad = drivers$rad[1:N]
  yearday = drivers$yearday[1:N]
  if (any(is.na(Cobs[,1]))){
    ind <- c(1:(dim(Cobs)[1]))[-which(is.na(Cobs[,1]))]
  } else{
    ind <- c(1:(dim(Cobs)[1]))
  }
  ind <- ind - 1
  missing_data_ind <- ind + 1
  buildCpredi <- Cpred(p_ = p, C_ = C, LAI_ = LAI, max_t_ = maxt, min_t_ = mint, Ca_ = Ca, lat_e_ = lat_e, yearday_ = yearday,
                       Nit_e_ = Nit_e, rad_ = rad, N = dim(Cobs)[1], G_ = G)
  Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                  buildCpredi$Clitpred, buildCpredi$Csompred)
  G_pred <- buildCpredi$G
  LL = 0
  for (i in 1:5){
    LL = LL + lgammad(sd[i], 10, 1)
  }
  #print(LL)
  LL = LL + SSLL(Cobs_ <- Cobs, C_ <- C, Sigma_ <- Sigma, Cpred_ <- Cpredi, sd_ <- sd, init_mean_ <- init_mean, 
                 init_sd_ <- init_sd, N <- N, G_ <- G, G_e, rep(.025, length(G_e)), g_tau_ <- rep(.2^2, length(G_e)),
                 Gpred_ <- G_pred, ind = ind)
  #return(list(LL = LL, gp = G_pred, cp = Cpredi))
  return(LL)
}


updateLatentStates <- function(Afun, bfun, pars, iter, Samples, init_mean, init_prec,
                               data, proc_prec, obs_prec, time.length, ...){
  if (!(is.array(Samples))) stop('Samples must be of type "array"')
  jmax <- dim(Samples)[1]
  imax <- dim(Samples)[2]
  for (i in sample(1:5,5)){
    j <- 1
    type <- 'first'
    print(j)
    print(i)
    t1 <- Sys.time()
    atp <- Afun(pars, i, j+1, Samples, iter, ...)
    btp <- bfun(pars, i, j+1, Samples, iter, ...)
    Samples[j, i, iter] <- get_latent_sample(atp = atp, btp = btp, initial_prec = init_prec[i], init_mean = init_mean[i], 
                                             add_prec = proc_prec[i], data = data[1,i], type = type, xtp = Samples[j+1, i, iter-1])
    t5 <- Sys.time()
    print(t5 - t1)
    #print(j)
    #print(atp)
    #print(btp)
    #print(Samples[j+1, i, iter-1])
    for (j in 2:time.length){
      print(j)
      print(i)
      type <- 'middle'
      t1 <- Sys.time()
      at <- Afun(pars, i, j, Samples, iter, ...)
      atp <- Afun(pars, i, j+1, Samples, iter, ...)
      #t2 <- Sys.time()
      bt <- bfun(pars, i, j, Samples, iter, ...)
      btp <- bfun(pars, i, j+1, Samples, iter, ...)
      #t3 <- Sys.time()
      #print(t2 - t1)
      #print(t3 - t2)

      #t4 <- Sys.time()
      Samples[j, i, iter] <- get_latent_sample(at = at, atp = atp, bt = bt, btp = btp,
                                                 add_prec = proc_prec[i], obs_prec = obs_prec[i],
                                                 xtp1 = Samples[j+1, i, iter-1], xtm1 = Samples[j-1, i, iter-1],
                                                 data = data[j-1, i], type = 'middle')
      t5 <- Sys.time()
      print(t5 - t1)
    }
    type <- 'last'
    at <- Afun(pars, i, jmax, Samples, iter, ...)
    bt <- bfun(pars, i, jmax, Samples, iter, ...)
    print(jmax)
    #profvis({
    Samples[jmax, i, iter] <- get_latent_sample(at = at, bt = bt, obs_prec = obs_prec[i], 
                                                    add_prec = proc_prec[i], xtm1 = Samples[jmax-1, i, iter-1],
                                                    data = data[jmax-1, i], type = 'last')
    #})
  }
  return(Samples[,,iter])
}

updateLatentStates2 <- function(Afunc, bfunc, pars, iter, Samples, init_mean, init_prec,
                               data, proc_prec, obs_prec, time.length, ...){
  jmax <- dim(Samples)[1]
  imax <- dim(Samples)[2]
  A <- Afunc(pars, Samples, ...)
  b <- bfunc(pars, Samples, ...)
  newsamp <- matrix(0, nrow = 731, ncol = 5)
  for (i in 1:imax){
    j <- 1
    type <- 'first'

    t1 <- Sys.time()
    atp <- A[j, i]
    btp <- b[j, i]
    datapt <- data[1,i]
    xtp1 = Samples[j+1, i]
    
    newsamp[j,i] <- get_latent_sample(atp = atp, btp = btp, initial_prec = init_prec[i], init_mean = init_mean[i], 
                                      add_prec = proc_prec[i], data = datapt, type = type, xtp1 = xtp1)
    t5 <- Sys.time()
    sum <- 0
    for (j in 2:time.length){
      type <- 'middle'

      at <- A[j-1, i]
      atp <- A[j, i]
      
      bt <- b[j-1, i]
      btp <- b[j, i]
      
      datapt <- data[j-1,i]
      xtp1 = Samples[j+1, i]
      xtm1 = Samples[j-1, i]
      
      newsamp[j,i] <- get_latent_sample(at = at, atp = atp, bt = bt, btp = btp,
                                        add_prec = proc_prec[i], obs_prec = obs_prec[i],
                                        xtp1 = xtp1, xtm1 = xtm1,
                                        data = datapt, type = 'middle')
    }
    #print(sum)
    type <- 'last'
    at <- A[jmax-1, i]
    bt <- b[jmax-1, i]
    #print(jmax)
    j <- jmax
    datapt <- data[j-1,i]
    xtm1 = Samples[j-1, i]
  
    newsamp[j,i] <- get_latent_sample(at = at, bt = bt, obs_prec = obs_prec[i], 
                                      add_prec = proc_prec[i], xtm1 = xtm1,
                                      data = datapt, type = 'last')
  }
  return(newsamp)
}

Afun2 <- function(pars, Samples, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, N, G, ...){
  if(dim(Samples)[1] == 731){
    Samples = Samples[2:731,]
  }
  buildA <- Cpred(p_ = pars[1,1:11], C_ = Samples, LAI_ = LAI, max_t_ = maxt, min_t_ = mint,
                  Ca_ = Ca, lat_e_ = lat_e, yearday_ = yearday, Nit_e_ = Nit_e, rad_ = rad, N = N, G_ = G, ...)
  A <- cbind(buildA$A_cf, buildA$A_cw, buildA$A_cr, buildA$A_clit, buildA$A_csom)
  return(A)
}

bfun2 <- function(pars, Samples, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, N, G, ...){
  if(dim(Samples)[1] == 731){
    Samples = Samples[2:731,]
  }
  buildB <- Cpred(p_ = pars[1,1:11], C_ = Samples, LAI_ = LAI, max_t_ = maxt, min_t_ = mint,
                  Ca_ = Ca, lat_e_ = lat_e, yearday_ = yearday, Nit_e_ = Nit_e, rad_ = rad, N = N, G_ = G, ...)
  b <- cbind(buildB$b_cf, buildB$b_cw, buildB$b_cr, buildB$b_clit, buildB$b_csom)
  return(b)
}

Afun <- function(pars, i, j, Samples, iter, mint, maxt, ...){
  j <- j-1
  if (i == 1){
    index <- which(colnames(pars) == 'p5')
    p5 <- pars[, index]
    at = 1 - p5
  }
  if (i == 2){
    index <- which(colnames(pars) == 'p6')
    p6 <- pars[, index]
    at = 1 - p6
  }
  if (i == 3){
    index <- which(colnames(pars) == 'p7')
    p7 <- pars[, index]
    at = 1 - p7
  }
  if (i == 4){
    index1 <- which(colnames(pars) == 'p10')
    index2 <- which(colnames(pars) == 'p1')
    index3 <- which(colnames(pars) == 'p8')
    p10 <- pars[, index1]
    p1 <- pars[, index2]
    p8 <- pars[, index3]
    at <- (1 - .5*exp(.5*p10*(mint[j] + maxt[j]))*(p8 + p1))
  }
  if (i == 5){
    index1 <- which(colnames(pars) == 'p9')
    index2 <- which(colnames(pars) == 'p10')
    p9 <- pars[, index1]
    p10 <- pars[, index2]
    at <- 1 - .5*p9*exp(.5*p10*(maxt[j] + mint[j]))
  }
  return(at)
}

bfun <- function(pars, i, j, Samples, iter, mint, maxt, G){
  j <- j-1
  if (i == 1){
    index1 <- which(colnames(pars) == 'p2')
    index2 <- which(colnames(pars) == 'p3')
    p2 <- pars[, index1]
    p3 <- pars[, index2]
    bt <- G[iter-1, j]*(1-p2)*p3
  }
  if (i == 2){
    index1 <- which(colnames(pars) == 'p2')
    index2 <- which(colnames(pars) == 'p3')
    index3 <- which(colnames(pars) == 'p4')
    p2 <- pars[, index1]
    p3 <- pars[, index2]
    p4 <- pars[, index3]
    bt <- G[iter-1, j]*(1-p2)*(1-p3)*(1-p4)
  }
  if (i == 3){
    index1 <- which(colnames(pars) == 'p2')
    index2 <- which(colnames(pars) == 'p3')
    index3 <- which(colnames(pars) == 'p4')
    p2 <- pars[, index1]
    p3 <- pars[, index2]
    p4 <- pars[, index3]
    bt <- G[iter-1, j]*(1-p2)*(1-p3)*(p4)
  }
  if (i == 4){
    index1 <- which(colnames(pars) == 'p5')
    index2 <- which(colnames(pars) == 'p7')
    p5 <- pars[, index1]
    p7 <- pars[, index2]
    bt <- p5*Samples[j, 1] + p7*Samples[j, 3]
  }
  if (i == 5){
    index1 <- which(colnames(pars) == 'p6')
    index2 <- which(colnames(pars) == 'p1')
    index3 <- which(colnames(pars) == 'p10')
    p6 <- pars[, index1]
    p1 <- pars[, index2]
    p10 <- pars[, index3]
    bt <- p6*Samples[j, 2] + .5*p1*exp(.5*p10*(mint[j] + maxt[j]))*Samples[j,4]
  }
  return(bt)
}

updateLatentStates3 <- function(A = NULL, Afunc = NULL, b = NULL, bfunc = NULL, pars, Samples, IC, init_mean, 
                                initial_prec, add_prec, obs_prec, data_obs, ...)
{
  ## Function checks to make sure that it has a way to access A values and b values,
  ## otherwise the function is stopped
  if (is.null(A) & is.null(Afunc)){
    stop('One of arguments A or Afunc cannot be null')
  }
  if (is.null(b) & is.null(bfunc)){
    stop('One of arguments b or bfunc cannot be null')
  }
  ## If A is not provided or b is not provided, values are computed using given functions
  if (is.null(A)){
    A <- Afunc(pars, Samples, ...)
  }
  if (is.null(b)){
    b <- bfunc(pars, Samples, ...)
  }
  ## Determines loop indices and allocates space for new samples
  I <- dim(Samples)[1]
  J <- dim(Samples)[2]
  newsamp <- matrix(0, nrow = I, ncol = J)
  newIC <- matrix(0, nrow = 1, ncol = J)
  ## Sample new initial conditions
  for (j in 1:J){
      newIC[1,j] <- get_latent_sample(atp = A[1,j], btp = b[1,j], initial_prec = initial_prec[j], 
                    add_prec = add_prec[j], xtp1 = Samples[1,j], data = data_obs[1,j], 
                    type = 'first', init_mean = init_mean[j])
  }
  ## Sample new start points
  for (j in 1:J){
    newsamp[1,j] <- get_latent_sample(at = A[1,j], atp = A[2,j], bt = b[1,j], btp = b[2,j],
                    add_prec = add_prec[j], obs_prec = obs_prec[j], xtp1 = Samples[2,j],
                    xtm1 = newIC[1,j], data = data_obs[1,j], type = 'middle')
  }
  ## Sample interior latent states
  for (i in 2:(I-1))
    {
      for (j in 1:J)
        {
            newsamp[i,j] <- get_latent_sample(at = A[i,j], atp = A[i+1,j], bt = b[i,j], btp = b[i+1,j],
                            add_prec = add_prec[j], obs_prec = obs_prec[j], xtm1 = newsamp[i-1,j],
                            xtp1 = Samples[i+1,j], data = data_obs[i,j], type = 'middle')
        }
    }
  ## Sample terminal latent states
  for (j in 1:J){
    newsamp[I,j] <- get_latent_sample(at = A[I,j], bt = b[I, j], obs_prec = obs_prec[j], add_prec = add_prec[j],
                    xtm1 = newsamp[I-1,j], data = data_obs[I, j], type = 'last')
  }
  ## Return new initial conditions and new samples
  return(list(newsamp = newsamp, newIC = newIC))
}

updateMCMC2 <- function(pars, parnames = colnames(pars), burn, lh_fun, block = NULL, i, upper, lower, transform = rep('none', length(parnames)),
                        adapt.length = 500, adapt.stop = 50000, sd.proposal, lh_old, acceptance, thin = 4, block.cov = NULL, check.accept = 3000, r = 1,
                        cov.adjust = rep(1, length(block)), method = "gibbs",  ...){
  ## Ensure that par is of class matrix
  if (!(is.null(dim(pars)))){
    new_par <- pars[i,]
  } else{
    new_par <- pars
  }
  if (is.null(parnames)) parnames <- names(pars)
  if (is.null(parnames)) parnames <- rep(0, length(upper))
  new_accept <- matrix(0, ncol = length(parnames), nrow = thin)
  ## Initializes storage for log threshold values and proposed points
  thresh <- rep(0,length(parnames))
  proposal <- rep(0,length(parnames))
  for (k in 1:thin){  
    ## Repeats loop for thin interations
    if (i < burn + adapt.length){
      ## Checks to see if adaptive sampling is to performed
      for (j in 1:length(parnames)){
        ## Loops over MCMC parameters
        if (i >= check.accept & i%%500 == 0){
          ## Tunes proposal standard deviations if wanted
          if (mean(acceptance[(i-499):i,j]) < .2){
            sd.proposal[j] <- sd.proposal[j] / 1.1
          } 
          if (mean(acceptance[(i-499):i,j]) > .4){
            sd.proposal[j] <- sd.proposal[j] * 1.1
          } 
        }
        ## Samples for log transformed variables
        if (transform[j] == 'log'){
          ## check if upper or lower bounds are NA
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            ## propose point
            p_star <- new_par
            tmp <- rtruncnorm(1, a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])
            ## compute proposal correction term
            prop.correction <- log(dtruncnorm(log(new_par[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm(log(p_star[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(new_par[j]), sd = sd.proposal[j]))
          } else{
            ## if neither upper or lower are NA
            ## propose point
            p_star <- new_par
            tmp <- rnorm(1, mean = log(p_star[j]), sd = sd.proposal[j])
            prop.correction <- 0
          }
          ## Transform proposed point back
          p_star[j] <- exp(tmp)
          lh_new <- lh_fun(pars = p_star, ...)
          
          threshold <- r*(lh_new - lh_old) + prop.correction
          
          ## acceptance step
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
          ## store threshold value and proposed point
          thresh[j] <- threshold
          proposal[j] <- exp(tmp)
        }
        if (transform[j] == 'none'){
          ## check if upper or lower bounds are NA
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            ## propose point
            p_star <- new_par
            tmp <- rtruncnorm(1, a = plower[j], b = pupper[j], mean = p_star[j], sd = sd.proposal[j])
            ## compute proposal correction term
            prop.correction <- log(dtruncnorm(new_par[j], a = plower[j], b = pupper[j], mean = log(p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm(p_star[j], a = plower[j], b = pupper[j], mean = new_par[j], sd = sd.proposal[j]))
          } else{
            ## if neither upper or lower are NA
            ## propose point
            p_star <- new_par
            tmp <- rnorm(1, mean = p_star[j], sd = sd.proposal[j])
            prop.correction <- 0
          }
          ## Transform proposed point back
          p_star[j] <- tmp
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          ## acceptance step
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
          ## store threshold value and proposed point
          thresh[j] <- threshold
          proposal[j] <- tmp
        } 
      }
    }
    ## check to see if current iteration is large enough to adapt
    if (i >= burn + adapt.length){
      if (i %% 500 == 0 & i < adapt.stop){
        ## create space for block.cov 
        block.cov <- vector('list', length(block))
        for (b in 1:length(block)){
          ## iterate through list of parameters to be block sampled and compute empirical correlation structures
          #blocks = block[[b]]
          #print(b)
          block.cov[[b]] <- cov(matrix(pars[burn:i, (block[[b]]) ], ncol = length(block[[b]])))
          if (any(is.na(cor(block.cov[[b]])))){
            block.cov[[b]] <- cov(matrix(pars[(burn-1000):i, (block[[b]]) ], ncol = length(block[[b]])))
          }
          #print(block.cov[[b]])
          #print(blocks[[b]])
          #print(block.cov[[b]])
          block.cov[[b]] <- block.cov[[b]] / cov.adjust[b]
          
        }
      }
      for (j in (1:length(parnames))[-unlist(block)]){
        if (transform[j] == 'log'){
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            p_star <- new_par
            tmp <- rtruncnorm(1, a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])
            prop.correction <- log(dtruncnorm(log(new_par[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm(log(p_star[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(new_par[j]), sd = sd.proposal[j]))
          } else{
            p_star <- new_par
            tmp <- rnorm(1, mean = log(p_star[j]), sd = sd[j])
            prop.correction <- 0
          }
          p_star[j] <- exp(tmp)
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
        }
        if (transform[j] == 'none'){
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            p_star <- new_par
            tmp <- rtruncnorm(1, a = (plower[j]), b = (pupper[j]), mean = (p_star[j]), sd = sd.proposal[j])
            prop.correction <- log(dtruncnorm((new_par[j]), a = (plower[j]), b = (pupper[j]), mean = (p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm((p_star[j]), a = (plower[j]), b = (pupper[j]), mean = (new_par[j]), sd = sd.proposal[j]))
          } else{
            p_star <- new_par
            tmp <- rnorm(1, mean = (p_star[j]), sd = sd[j])
            prop.correction <- 0
          }
          p_star[j] <- (tmp)
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
        } 
      }
      for (b in 1:length(block)){
        blocks <- block[[b]]
        p_star <- new_par
        if (all(is.na(pupper[blocks]))){
          if (method == 'gibbs'){
            start.value = pars[i, blocks]
            burn.in.samples = 500
            tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], sigma = block.cov[[b]], algorithm = method, start.value = start.value, 
                            burn.in.samples = burn.in.samples)
          } else{
            tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], sigma = block.cov[[b]], algorithm = method)
          }
        } else{
          if (method == 'gibbs'){
            start.value = pars[i, blocks]
            burn.in.samples = 500
            #print(block.cov[[b]])
            tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]], algorithm = method, start.value = start.value, 
                            burn.in.samples = burn.in.samples)
          } else{
            tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]], algorithm = method)
          }
        }
        p_star[blocks] <- tmp
        lh_new <- lh_fun(pars = p_star, ...)
        if (any(is.na(pupper[blocks]))){
          prop.correction <- log(dtmvnorm(new_par[blocks], lower = plower[blocks], mean = p_star[blocks], sigma = block.cov[[b]])) - 
            log(dtmvnorm(x = p_star[blocks], mean = new_par[blocks], lower = plower[blocks], sigma = block.cov[[b]]))
        } else{
          if (is.nan(log(dtmvnorm(new_par[blocks], lower = plower[blocks], upper = pupper[blocks], mean = p_star[blocks], sigma = block.cov[[b]])) - 
                     log(dtmvnorm(x = p_star[blocks], mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]])))){
            prop.correction <- 0
          } else{
            prop.correction <- log(dtmvnorm(new_par[blocks], lower = plower[blocks], upper = pupper[blocks], mean = p_star[blocks], sigma = block.cov[[b]])) - 
              log(dtmvnorm(x = p_star[blocks], mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]])) 
          }
        }
        #print(prop.correction)
        threshold <- r*(lh_new - lh_old) + prop.correction
        if (log(runif(1,0,1)) < threshold){
          lh_old <- lh_new
          new_accept[k,blocks] <- 1
          new_par <- p_star
        }
      }
    }
  }
  pars <- rbind(pars, new_par)
  acceptance <- rbind(acceptance, colMeans(new_accept))
  return(list(pars = pars, acceptance = acceptance, lh_old = lh_old, cov = block.cov, sd = sd.proposal, thresh = thresh, p_star = proposal))
}

updateMCMC2_mw <- function(pars, parnames = colnames(pars), burn, lh_fun, block = NULL, i, upper, lower, transform = rep('none', length(parnames)),
                        adapt.length = 500, adapt.stop = 50000, sd.proposal, lh_old, acceptance, thin = 4, block.cov = NULL, check.accept = 3000, r = 1,
                        cov.adjust = rep(1, length(block)), method = "gibbs", fixed = rep('false', length(colnames(pars))),  ...){
  ## Ensure that par is of class matrix
  if (!(is.null(dim(pars)))){
    new_par <- pars[i,]
  } else{
    new_par <- pars
  }
  if (is.null(parnames)) parnames <- names(pars)
  if (is.null(parnames)) parnames <- rep(0, length(upper))
  new_accept <- matrix(0, ncol = length(parnames), nrow = thin)
  ## Initializes storage for log threshold values and proposed points
  thresh <- rep(0,length(parnames))
  proposal <- rep(0,length(parnames))
  for (k in 1:thin){  
    ## Repeats loop for thin interations
    if (i < burn + adapt.length){
      ## Checks to see if adaptive sampling is to performed
      for (j in 1:length(parnames)){
        #print(j)
        ## Loops over MCMC parameters
        if (i >= check.accept & i%%500 == 0){
          ## Tunes proposal standard deviations if wanted
          if (mean(acceptance[(i-499):i,j]) < .2){
            sd.proposal[j] <- sd.proposal[j] / 1.1
          } 
          if (mean(acceptance[(i-499):i,j]) > .4){
            sd.proposal[j] <- sd.proposal[j] * 1.1
          } 
        }
        ## Samples for log transformed variables
        if (transform[j] == 'log' & fixed[j] == 'false'){
          ## check if upper or lower bounds are NA
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            ## propose point
            #print(j)
            p_star <- new_par
            tmp <- rtruncnorm(1, a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])
            ## compute proposal correction term
            prop.correction <- log(dtruncnorm(log(new_par[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm(log(p_star[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(new_par[j]), sd = sd.proposal[j]))
          } else{
            ## if neither upper or lower are NA
            ## propose point
            p_star <- new_par
            tmp <- rnorm(1, mean = log(p_star[j]), sd = sd.proposal[j])
            prop.correction <- 0
          }
          ## Transform proposed point back
          p_star[j] <- exp(tmp)
          #print(p_star[j])
          lh_new <- lh_fun(pars = p_star, ...)
          #print(lh_new)
          threshold <- r*(lh_new - lh_old) + prop.correction
          #print(threshold)
          
          ## acceptance step
          #print(j)
          #print(lh_new)
          #print(prop.correction)
          #print(lh_old)
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
          ## store threshold value and proposed point
          thresh[j] <- threshold
          proposal[j] <- exp(tmp)
        }
        if (transform[j] == 'none' & fixed[j] == 'false'){
          ## check if upper or lower bounds are NA
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            ## propose point
            p_star <- new_par
            tmp <- rtruncnorm(1, a = plower[j], b = pupper[j], mean = p_star[j], sd = sd.proposal[j])
            ## compute proposal correction term
            prop.correction <- log(dtruncnorm(new_par[j], a = plower[j], b = pupper[j], mean = log(p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm(p_star[j], a = plower[j], b = pupper[j], mean = new_par[j], sd = sd.proposal[j]))
          } else{
            ## if neither upper or lower are NA
            ## propose point
            p_star <- new_par
            tmp <- rnorm(1, mean = p_star[j], sd = sd.proposal[j])
            prop.correction <- 0
          }
          ## Transform proposed point back
          p_star[j] <- tmp
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          #print(threshold)
          ## acceptance step
          #print(j)
          #print(lh_new)
          #print(prop.correction)
          #print(lh_old)
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
          ## store threshold value and proposed point
          thresh[j] <- threshold
          proposal[j] <- tmp
        } 
      }
    }
    ## check to see if current iteration is large enough to adapt
    if (i >= burn + adapt.length){
      if (i %% 500 == 0 & i < adapt.stop){
        ## create space for block.cov 
        block.cov <- vector('list', length(block))
        for (b in 1:length(block)){
          ## iterate through list of parameters to be block sampled and compute empirical correlation structures
          #blocks = block[[b]]
          #print(b)
          block.cov[[b]] <- cov(matrix(pars[burn:i, (block[[b]]) ], ncol = length(block[[b]])))
          if (any(is.na(cor(block.cov[[b]])))){
            block.cov[[b]] <- cov(matrix(pars[(i-2000):i, (block[[b]]) ], ncol = length(block[[b]])))
          }
          #print(block.cov[[b]])
          #print(blocks[[b]])
          #print(block.cov[[b]])
          block.cov[[b]] <- block.cov[[b]] / cov.adjust[b]
          
        }
      }
      for (j in (1:length(parnames))[-unlist(block)]){
        if (transform[j] == 'log' & fixed[j] == 'false'){
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            p_star <- new_par
            tmp <- rtruncnorm(1, a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])
            prop.correction <- log(dtruncnorm(log(new_par[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm(log(p_star[j]), a = log(plower[j]), b = log(pupper[j]), mean = log(new_par[j]), sd = sd.proposal[j]))
          } else{
            p_star <- new_par
            tmp <- rnorm(1, mean = log(p_star[j]), sd = sd[j])
            prop.correction <- 0
          }
          p_star[j] <- exp(tmp)
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
        }
        if (transform[j] == 'none' & fixed[j] == 'false'){
          if (!(is.na(upper[j])) & !(is.na(lower[j]))){
            p_star <- new_par
            tmp <- rtruncnorm(1, a = (plower[j]), b = (pupper[j]), mean = (p_star[j]), sd = sd.proposal[j])
            prop.correction <- log(dtruncnorm((new_par[j]), a = (plower[j]), b = (pupper[j]), mean = (p_star[j]), sd = sd.proposal[j])) - 
              log(dtruncnorm((p_star[j]), a = (plower[j]), b = (pupper[j]), mean = (new_par[j]), sd = sd.proposal[j]))
          } else{
            p_star <- new_par
            tmp <- rnorm(1, mean = (p_star[j]), sd = sd[j])
            prop.correction <- 0
          }
          p_star[j] <- (tmp)
          lh_new <- lh_fun(pars = p_star, ...)
          threshold <- r*(lh_new - lh_old) + prop.correction
          if (log(runif(1,0,1)) < threshold){
            lh_old <- lh_new
            new_accept[k,j] <- 1
            new_par <- p_star
          }
        } 
      }
      for (b in 1:length(block)){
        blocks <- block[[b]]
        p_star <- new_par
        if (all(is.na(pupper[blocks]))){
          if (method == 'gibbs'){
            start.value = pars[i, blocks]
            burn.in.samples = 500
            tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], sigma = block.cov[[b]], algorithm = method, start.value = start.value, 
                            burn.in.samples = burn.in.samples)
          } else{
            tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], sigma = block.cov[[b]], algorithm = method)
          }
        } else{
          if (method == 'gibbs'){
            start.value = pars[i, blocks]
            burn.in.samples = 500
            #print(block.cov[[b]])
            tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]], algorithm = method, start.value = start.value, 
                            burn.in.samples = burn.in.samples)
          } else{
            tmp <- rtmvnorm(1, mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]], algorithm = method)
          }
        }
        #print(tmp)
        p_star[blocks] <- tmp
        lh_new <- lh_fun(pars = p_star, ...)
        if (any(is.na(pupper[blocks]))){
          #print(log(dtmvnorm(new_par[blocks], lower = plower[blocks], mean = p_star[blocks], sigma = block.cov[[b]])))
          prop.correction <- log(dtmvnorm(new_par[blocks], lower = plower[blocks], mean = p_star[blocks], sigma = block.cov[[b]])) - 
            log(dtmvnorm(x = p_star[blocks], mean = new_par[blocks], lower = plower[blocks], sigma = block.cov[[b]]))
        } else{
          if (is.nan(log(dtmvnorm(new_par[blocks], lower = plower[blocks], upper = pupper[blocks], mean = p_star[blocks], sigma = block.cov[[b]])) - 
                     log(dtmvnorm(x = p_star[blocks], mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]])))){
            prop.correction <- 0
          } else{
            prop.correction <- log(dtmvnorm(new_par[blocks], lower = plower[blocks], upper = pupper[blocks], mean = p_star[blocks], sigma = block.cov[[b]])) - 
              log(dtmvnorm(x = p_star[blocks], mean = new_par[blocks], lower = plower[blocks], upper = pupper[blocks], sigma = block.cov[[b]])) 
          }
        }
        #print(prop.correction)
        threshold <- r*(lh_new - lh_old) + prop.correction
        if (log(runif(1,0,1)) < threshold){
          lh_old <- lh_new
          new_accept[k,blocks] <- 1
          new_par <- p_star
        }
      }
    }
  }
  pars <- rbind(pars, new_par)
  acceptance <- rbind(acceptance, colMeans(new_accept))
  return(list(pars = pars, acceptance = acceptance, lh_old = lh_old, cov = block.cov, sd = sd.proposal, thresh = thresh, p_star = proposal))
}
