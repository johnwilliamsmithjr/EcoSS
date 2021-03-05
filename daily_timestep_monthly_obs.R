## daily_timestep_monthly_obs.R
## Perform MCMC for daily timestep latent state model 
## with monthly observations

## generate data and initialize MCMC
source('./daily_timestep_month_obs_setup.R')

t1 <- Sys.time()
for (i in 2:chain_length){
  ## thin samples before burn-in
  if (i <= burn){
    thin = 4
  } else{
    thin = 4
  }
  ## Update process parameters and precision parameters
  ## Function returns new likelihood, new parameters, threshold (if wanted), the proposed points (if wanted),
  ## and acceptance rate. Thin controls the rate at which samples are saved
  param_update <- updateMCMC(pars = pars, burn = 4000, lh_fun = cpp_ss_likelihood_ex, i = i-1,
                             upper = pupper, lower = plower, transform = rep('log', 16),
                             sd.proposal = sd, lh_old = old, block = block, acceptance = accept, 
                             Cobs = C_obs, C = C_samples[,,i-1], var_obs = var_obs, init_mean = init_mean,
                             init_sd = init_var, N = nday, G = G[i-1,], G_e = G_obs, drivers = drivers_ev, 
                             LAI = pmax(C_samples[,1,i-1]/LMA_e, .1), fluxes = fluxes_obs, block.cov = cov.block, thin = thin
                             , cov.adjust = c(1,1,1,1,1), method = 'gibbs')
  ## update parameters, likelihood, proposal standard deviations, acceptance
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
  ## get predictions
  buildCpredi <- Cpred(p, C_samples[, , i-1], LAI_ = pmax(.1, C_samples[,1,i-1]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, nday, G[i-1,])
  Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                  buildCpredi$Clitpred, buildCpredi$Csompred)
  
  ## update process precisions with Gibbs updates
  for (k in 1:5){
    prec_add_sample[k] <- rgamma(1, nday/2, 
                                 sum((C_samples[,k, i-1] - Cpredi[,k])^2)*.5)
  }
  
  ## update likelihood
  pars[i, 12:16] <- prec_add_sample
  old <- cpp_ss_likelihood_ex(Cobs = C_obs, C = C_samples[,,i-1], var_obs = var_obs, pars = pars[i,], init_mean = init_mean,
                              init_sd = init_var, N = nday, G = G[i-1,], G_e = G_obs, drivers = drivers_ev, LAI = pmax(C_samples[,1,i-1]/LMA_e),
                              fluxes = fluxes_obs)
  ## extract A_i, b_i from process model
  A <- Afun2(pars = pars[i,,drop = FALSE], Samples = C_samples[,,i-1], LAI = pmax(C_samples[,1,i-1]/LMA_e), 
             maxt = maxt, mint = mint, Ca = Ca, lat_e = lat_e, yearday = yearday, Nit_e = Nit_e, rad = rad, N = nday, G = G[i-1,])
  b <- bfun2(pars = pars[i,,drop = FALSE], Samples = C_samples[,,i-1], LAI = pmax(C_samples[,1,i-1]/LMA_e), 
             maxt = maxt, mint = mint, Ca = Ca, lat_e = lat_e, yearday = yearday, Nit_e = Nit_e, rad = rad, N = nday, G = G[i-1,])
  
  ## Update latent states
  LSU <- updateLatentStates3(A = A, b = b, pars = pars[i,,drop = FALSE], Samples = C_samples[,,i-1], IC = C.0[i-1,], init_mean = init_mean,
                             initial_prec = prec_init, data = C_obs, add_prec = prec_add_sample, obs_prec = prec_obs)
  C_samples[,,i] <- LSU$newsamp
  C.0[i,] <- LSU$newIC
  ## Update G
  buildCp <- Cpred(p, C_samples[,,i], LAI_ = pmax(C_samples[,1,i]/LMA_e), max_t_ = maxt, min_t_ = mint, Ca_ = Ca,
                   lat_e_ = lat_e, yearday_ = yearday, Nit_e_ = Nit_e, rad_ = rad, N = nday, G_ = G[i-1,])
  Gpred = buildCp$G
  rm(buildCp)
  Gmeans <- (phi_g*Gpred + tau_g*G_obs) / (phi_g + tau_g)
  G_sd <- 1 / sqrt(phi_g + tau_g)
  G[i,] <- rnorm(nday, mean = Gmeans, sd = G_sd)
  
  ## update likelihood
  old <- cpp_ss_likelihood_ex(Cobs = C_obs, C = C_samples[,,i], var_obs = var_obs, pars = pars[i,], init_mean = init_mean,
                              init_sd = init_var, N = nday, G = G[i,], G_e = G_obs, drivers = drivers_ev, LAI = pmax(C_samples[,1,i]/LMA_e),
                              fluxes = fluxes_obs)
  ## Verbose updates, including trace plots, time elapsed, and acceptance rates over the last 500 iterations
  if (i %% 500 == 0 & update == TRUE){
    print(paste0('Acceptance over last 500 iterations ',colMeans2(accept[(i-499):i,])))
    print(pars[i,])
    for (var in 1:16){
      plot(pars[1:i, var], type = 'l', main = colnames(pars)[var])
    }
    print(Sys.time() - t1)
    if (i >= burn+500) print(effectiveSize(pars[burn:i, ]))
    print(old)
  }
}
