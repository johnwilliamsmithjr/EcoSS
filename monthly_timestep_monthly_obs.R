## monthly_timestep_monthly_obs.R
## perform state space MCMC on DALECev model with monthly timestep model
## and monthly observed data

source('./setup/data_sim/monthly_timestep_monthly_obs_datasim.R')

source('./setup/MCMC/monthly_timestep_monthly_obs_MCMC.R')

##source('./monthly_timestep_monthly_obs_setup.R')
## time process
t1 <- Sys.time()

for (i in 2:chain_length){
  ## Update process parameters and precision parameters
  ## Function returns new likelihood, new parameters, threshold (if wanted), the proposed points (if wanted),
  ## and acceptance rate. Thin controls the rate at which samples are saved
  if (i == burn){
    sd[6] <- sd[6] / 100
    sd[9] <- sd[9] / 50
    sd[10] <- sd[10] / 15
    sd[11] <- sd[11] / 50
  }
  param_update <- updateMCMC2(pars = pars, burn = burn, lh_fun = ss_monthly_ll, i = i-1,
                              upper = pupper, lower = plower, transform = tf,
                              sd.proposal = sd, lh_old = old, block = block, acceptance = accept, 
                              C_obs = C_obs, C = C_samples[,,i-1], prec_obs = prec_obs, init_mean = C.0[i-1,],
                              init_var = init_var, ndays = ndays, G = G_obs, drivers = drivers_ev, init_ev = init_ev,
                              fluxes = fluxes_obs, block.cov = cov.block, thin = 3, r  = r, 
                              Ca = Ca, lat_e = lat_e, yearday = yearday, rad = rad, N = nday, Nit_e = Nit_e, method = 'gibbs', 
                              cov.adjust = c(1,1,1), preds = preds)
  pars <- param_update$pars
  old <- param_update$lh_old
  sd <- param_update$sd
  accept <- param_update$acceptance
  ## function can be set to block sample certain parameters after burn in
  if (!(is.null(param_update$cov))) cov.block <- param_update$cov
  ## Latent State Updating
  ## Function returns new sample for current iteration for both latent states and initial conditions
  p <- pars[i, 1:11]
  prec_add_sample <- pars[i, 12:16]
  preds <- getDALECpredictions2(p_e = p, C = C_samples[,,i-1], LAI = preds$Cfast[,1] / (LMA_e), maxt = maxt, mint = mint,
                                Ca = Ca, lat_e = lat_e, yearday = yearday, Nit_e = Nit_e, rad = rad, N = nday, G_obs = G_obs, 
                                ind = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6))-1, 
                                init_mean = C.0[i-1,], native_step = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6)))
  Cpredi <- preds$Cslow
  for (k in 1:5){
    prec_add_sample[k] <- rgamma(1, r*nmonths/2, 
                                 sum((C_samples[,k, i-1] - Cpredi[,k])^2)*.5*r)
  }
  pars[i, 12:16] <- prec_add_sample
  
  old <- ss_monthly_ll(C_obs, C = C_samples[,,i-1], C.0[i-1,], init_var, prec_obs = prec_obs, pars = pars[i,], fluxes = fluxes_obs,
                       G = G_obs, ndays = ndays, drivers = drivers_ev, init_ev = init_ev, 
                       Ca = Ca, lat_e = lat_e, yearday = yearday, rad = rad, N = nday, Nit_e = Nit_e, preds = preds, r = r)
  A <- preds$A
  b <- preds$b
  
  ## Update latent states
  LSU <- updateLatentStates3(A = A, b = b, pars = pars[i,,drop = FALSE], Samples = C_samples[,,i-1], IC = C.0[i-1,], init_mean = init_mean,
                             initial_prec = prec_init, data = C_obs, add_prec = prec_add_sample, obs_prec = prec_obs)
  C_samples[,,i] <- LSU$newsamp
  C.0[i,] <- LSU$newIC
  ## Update G
  buildCp <- getDALECpredictions2(p_e = p, C = C_samples[,,i], LAI = preds$Cfast[,1] / (LMA_e), maxt = maxt, mint = mint,
                                  Ca = Ca, lat_e = lat_e, yearday = yearday, Nit_e = Nit_e, rad = rad, N = nday, G_obs = G_obs, 
                                  ind = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6))-1, 
                                  init_mean = C.0[i,], native_step = cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 6)))
  Gpred = buildCp$G
  rm(buildCp)
  old <- ss_monthly_ll(C_obs, C = C_samples[,,i-1], init_mean, init_var, prec_obs = prec_obs, pars = pars[i,], fluxes = fluxes_obs,
                       G = G_obs, ndays = ndays, drivers = drivers_ev, init_ev = init_ev,
                       Ca = Ca, lat_e = lat_e, yearday = yearday, rad = rad, N = nday, Nit_e = Nit_e, preds = preds, r = r)
  ## Verbose updates, including trace plots, time elapsed, and acceptance rates over the last 500 iterations
  if (i %% 100 == 0 & update == TRUE){
    print(paste0('Acceptance over last 100 iterations ',colMeans2(accept[(i-99):i,])))
    print(pars[i,])
    for (var in 1:16){
      plot(pars[1:i, var], type = 'l', main = colnames(pars)[var])
    }
    print(Sys.time() - t1)
  }
}
