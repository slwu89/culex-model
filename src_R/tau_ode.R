# --------------------------------------------------------------------------------
#   ODEs describing change in delay duration
# --------------------------------------------------------------------------------

tau_diffeqn <- function(t, y, params){
  
  # state variables
  DE <- y[1] # tau_E(t)
  DL <- y[2] # tau_L(t)
  DP <- y[3] # tau_P(t)
  
  # (lagged) temperature
  temp <- temperature(t, params)
  temp_E <- temperature(t - DE, params)
  temp_L <- temperature(t - DL, params)
  temp_P <- temperature(t - DP, params)
  
  # (lagged) development
  larvae_maturation <- larvae_maturation_rate(temp,params)
  larvae_maturation_L <- larvae_maturation_rate(temp_L,params)
  
  egg_maturation <- egg_maturation_rate(temp,params)
  egg_maturation_E <- egg_maturation_rate(temp_E,params)
  
  pupae_maturation <- pupae_maturation_rate(temp,params)
  pupae_maturation_P <- pupae_maturation_rate(temp_P,params)
  
  # DDEs describing change in state duration
  dDEdt = 1 - egg_maturation/egg_maturation_E
  dDLdt = 1 - larvae_maturation/larvae_maturation_L
  dDPdt = 1 - pupae_maturation/pupae_maturation_P
  
  # DDE system
  du <- rep(NaN, 3)
  
  du[1] = dDEdt # tau_E(t)
  du[2] = dDLdt # tau_L(t)
  du[3] = dDPdt # tau_P(t)
  
  return(list(du))
}


# --------------------------------------------------------------------------------
#   external forcing
# --------------------------------------------------------------------------------

# temperature as modified cosine function
temperature <- function(t, pars){
  
  phi = pars[["phi"]] # PHASE
  lambda = pars[["lambda"]] # A
  mu = pars[["mu"]] # M
  gamma = pars[["gamma"]] # POWER
  
  temp = 0.0
  
  if(t < 0.0){
    temp = (mu - lambda) + lambda * 2.0 * (0.5 * (1.0 + cos(2.0 * pi * (0.0 - phi) / 365.0)))^gamma
  } else {
    temp = (mu - lambda) + lambda * 2.0 * (0.5 * (1.0 + cos(2.0 * pi * (t - phi) / 365.0)))^gamma
  }
  
  return(temp)
}


# --------------------------------------------------------------------------------
#   lifecycle stage progression rates
# --------------------------------------------------------------------------------

# g_E
egg_maturation_rate <- function(temp, pars){
  
  alpha_E = pars[["alpha_E"]] # ALPHA
  beta_E = pars[["beta_E"]] # BETA
  maturation_min = pars[["maturation_min"]]
  
  # calculate egg development rate
  if(temp < 0.0){
    egg_maturation = 0.016667
  } else {
    egg_maturation = alpha_E * (temp^beta_E)
  }
  
  if(egg_maturation < maturation_min){
    egg_maturation = maturation_min
  }
  return(egg_maturation)
}

# g_L
larvae_maturation_rate <- function(temp, pars){
  
  alpha_L = pars[["alpha_L"]] # ALPHA
  beta_L = pars[["beta_L"]] # BETA
  maturation_min = pars[["maturation_min"]]
  
  # calculate larvae development rate
  if(temp < 0.0){
    larvae_maturation = 0.016667
  } else {
    larvae_maturation = alpha_L * (temp^beta_L)
  }
  
  if(larvae_maturation < maturation_min){
    larvae_maturation = maturation_min
  }
  return(larvae_maturation)
}

# g_P
pupae_maturation_rate <- function(temp, pars){
  
  alpha_P = pars[["alpha_P"]] # ALPHA
  beta_P = pars[["beta_P"]] # BETA
  maturation_min = pars[["maturation_min"]]
  
  # calculate larvae development rate
  if(temp < 0.0){
    pupae_maturation = 0.016667
  } else {
    pupae_maturation = alpha_P * (temp^beta_P)
  }
  
  if(pupae_maturation < maturation_min){
    pupae_maturation = maturation_min
  }
  return(pupae_maturation)
}
