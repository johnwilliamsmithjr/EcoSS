cpp_ss_likelihood = function(Cobs, C, Sigma, Cpred, sd, init_mean, init_sd, N, G, G_e, G_pred, ind){
  LL = 0
  for (i in 1:5){
    LL = LL + lgammad(sd[i], 10, 1)
  }
  LL = LL + SSLL(Cobs_ <- Cobs, C_ <- C, Sigma_ <- Sigma, Cpred_ <- Cpred, sd_ <- sd, init_mean_ <- init_mean, 
                 init_sd_ <- init_sd, N <- N, G_ <- G, G_e, rep(.025, length(G_e)), g_tau_ <- rep(.1^2, length(G_e)),
                 Gpred_ <- G_pred, ind = ind)
  #LL = LL + dmvnorm(G, G_e, sigma = diag(.025, length(G_e)), log = T)
  return(LL)
}

Lr_ll <- function(p, C, Lr_obs){
  Lr_ll <- sum(dnorm(Lr_obs, p[7]*C[,3], sd = .001, log = TRUE))
  return(Lr_ll)
}

cpp_ss_likelihood2 = function(Cobs, C, Sigma, Cpred, sd, init_mean, init_sd, N, G, G_e, ind){
  LL = 0
  #for (i in 1:5){
  #  LL = LL + lgammad(sd[i], 10, 1)
  #}
  LL = LL + SSLL(Cobs, C, Sigma, Cpred, sd, init_mean, init_sd, N, G, G_e, rep(.025, length(G_e)), ind = ind)
  #LL = LL + dmvnorm(G, G_e, sigma = diag(.025, length(G_e)), log = T)
  return(LL)
}
