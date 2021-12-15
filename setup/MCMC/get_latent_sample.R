## Program: get_latent_sample.R
## Author: John W Smith Jr
## Date: February 25th, 2019
## Description: Generates draw from latent distribution in a state space model with linear autoregressive process model

get_latent_sample = function(at = NULL, bt = NULL, atp = NULL, btp = NULL, type, initial_prec = NULL, add_prec, obs_prec = NULL, data,
                             xtm1 = NULL, xtp1 = NULL, init_mean = NULL){
  ## Check for valid type argument
  if (!(type %in% c('first', 'middle', 'last'))){
    stop('Invalid type argument')
  }
  ## Samples for interior timepoints
  if (type == 'middle'){
    if (!(is.na(data))){
      # Compute mean and variance for non-missing data
      vrnc = 1 / (add_prec*(1 + (atp^2)) + obs_prec )
      mu = (add_prec*(at*xtm1 + bt) + add_prec*(atp*xtp1 - atp*btp) + obs_prec*data)*vrnc
      #print(vrnc)
      #print(mu)
      # sample from distribution and return
      sample = rnorm(1, mean = mu, sd = sqrt(vrnc))
      return(sample)
    } else {
      ## compute mean and variance for missing data
      vrnc = 1 / (add_prec*(1 + (atp^2)))
      mu = add_prec*(at*xtm1 + bt + xtp1*atp - atp*btp)*vrnc
      # sample from distribution and return
      sample = rnorm(1, mean = mu, sd = sqrt(vrnc))
      return(sample)
    }
  }
  ## Samples for terminal timepoint
  if (type == 'last'){
    if (!(is.na(data))){
      # Compute mean and variance for non-missing data
      vrnc = 1 / (add_prec + obs_prec)
      mu = (add_prec*(at*xtm1 + bt) + obs_prec*data)*vrnc
      # sample from distribution and return
      sample = rnorm(1, mean = mu, sd = sqrt(vrnc))
      return(sample)
    } else{
      # Compute mean and variance for missing data
      mu = at*xtm1 + bt
      vrnc = 1 / (obs_prec)
      # sample from distribution and return
      sample = rnorm(1, mean = mu, sd = sqrt(vrnc))
      return(sample)
    }
  } 
  ## Samples for first timepoint
  if (type == 'first'){
    if (!(is.na(data))){
      # Compute mean and variance for non-missing data
      vrnc = 1 / ((atp^2)*add_prec + initial_prec)
      mu = (add_prec*atp*xtp1 - add_prec*atp*btp + initial_prec*init_mean)*vrnc
      #print(vrnc)
      #print(mu)
      # Sample from distribution and return
      sample = rnorm(1, mean = mu, sd = sqrt(vrnc))
      return(sample)
    }
  } else{
    # Compute mean and variance for missing data
    vrnc = 1 / initial_prec
    mu = init_mean
    # sample from distribution and return
    sample = rnorm(1, mean = mu, sd = sqrt(vrnc))
    return(sample)
  }
}

