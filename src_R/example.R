library(Rcpp)
library(RcppArmadillo)
library(deSolve)

library(parallel)
library(data.table)
library(ggplot2)

rm(list=ls());gc()

Rcpp::sourceCpp("src_R/culex.cpp")
source("src_R/tau_ode.R")

# parameters from https://github.com/davewi13/Temperate-Mosquito-DDE/blob/master/Chapter%202%20DDE%20code.f90
parameters = list(
  # temp
  "phi" = 1.4, # PHASE
  "lambda" = 6.3, # A
  "mu" = 10.3, # M
  "gamma" = 1.21, # POWER
  # photoperiod
  "L" = 51, # latitude
  # oviposition
  "max_egg" = 200, # max egg raft size, R
  # gonotrophic cycle
  "q1" = 0.2024, # KG
  "q2" = 74.48, # QG
  "q3" = 0.2456, # BG
  # egg death
  "nu_0E" = 0.0157,  # U3
  "nu_1E" = 20.5,  # U4
  "nu_2E" = 7,  # U5
  # larvae death
  "nu_0L" = 0.0157,  # U3
  "nu_1L" = 20.5,  # U4
  "nu_2L" = 7,  # U5
  # pupae death
  "nu_0P" = 0.0157,  # U3
  "nu_1P" = 20.5,  # U4
  "nu_2P" = 7,  # U5
  # adult death
  "alpha_A" = 2.166e-8, # ALPHA
  "beta_A" = 4.483, # BETA
  # egg maturation
  "alpha_E" = 0.0022, # ALPHA
  "beta_E" = 1.77, # BETA
  # larvae maturation
  "alpha_L" = 0.00315, # ALPHA
  "beta_L" = 1.12, # BETA
  # pupae maturation
  "alpha_P" = 0.0007109, # ALPHA
  "beta_P" = 1.8865648, # BETA
  # predation on pupae
  "a" = 1,
  "h" = 0.002,
  "r" = 0.001,
  "V" = 200,
  # max values
  "death_max" = 1.0,
  "death_min_a" = 0.01,
  "gon_min" = 0.0333,
  "maturation_min" = 0.016667
)

# calculate the combined parameters for predation on larvae
calc_pred <- function(pars){
  pars[["p0"]] = pars[["r"]] / pars[["h"]]
  pars[["p1"]] = pars[["V"]] / (pars[["a"]] * pars[["h"]])
  return(pars)
}

parameters <- calc_pred(parameters)

# time horizon
dt <- 0.1
tmax <- 365*5
maxstep <- tmax/dt

# initial # of adults
A0 <- 12000

# integrate ODEs describing maturation delays

# delays at t=0
tau0 <- c("E"=0,"L"=0,"P"=0)
temp0 <- temperature(0, parameters)
tau0[1] = 1 / egg_maturation_rate(temp0, parameters) # tau_E
tau0[2] = 1 / larvae_maturation_rate(temp0, parameters) # tau_L
tau0[3] = 1 / pupae_maturation_rate(temp0, parameters) # tau_P

# integrate past the simulation end time by a comfortable amount
times <- seq(from=0.0,to=tmax+200,by=dt)
tau_ode <- deSolve::ode(y = tau0, times = times, func = tau_diffeqn, parms = parameters, method = "ode23")

tau_traj <- as.data.table(tau_ode)
tau_traj[, "step" := seq_along(times)]
# delays are now in units of time steps
tau_traj[, c("E", "L", "P") := lapply(.SD, function(x){as.integer(round(x = x/dt))}), .SDcols = c("E", "L", "P")]
tau_traj[, "time" := NULL]

# go from backward-looking maturation delays to forward-looking delays
tau_traj[, "step_forward_E" := as.integer(step - E)]
tau_traj[, "step_forward_L" := as.integer(step - L)]
tau_traj[, "step_forward_P" := as.integer(step - P)]
# the step_forward_X columns is the time step at which new X's need to wait
# tauX timesteps until going to the next stage. If the last row in any of
# them is < maxstep, we didn't integrate the ODEs far enough ahead.
tauE <- as.integer(tau_traj[step_forward_E %in% 1:maxstep, "E"][[1]])
tauL <- as.integer(tau_traj[step_forward_E %in% 1:maxstep, "L"][[1]])
tauP <- as.integer(tau_traj[step_forward_E %in% 1:maxstep, "P"][[1]])

# draw samples
out <- parallel::mclapply(X = 1:10, FUN = function(runid) {
  mod <- create_culex_stochastic(p = 1, tau_E = tauE, tau_L = tauL, tau_P = tauP, dt = dt)
  set_A_stochastic(mod = mod, A = A0)
  
  out <- matrix(data = NaN, nrow = tmax, ncol = 4)
  out_i <- 1L
  
  for (i in 1:maxstep) {
    step_culex_stochastic(mod = mod, parameters = parameters)
    if ((i-1) %% (1/dt) == 0) {
      out[out_i, 1] <- get_E_stochastic(mod = mod)
      out[out_i, 2] <- get_L_stochastic(mod = mod)
      out[out_i, 3] <- get_P_stochastic(mod = mod)
      out[out_i, 4] <- get_A_stochastic(mod = mod)
      out_i <- out_i + 1L
    }
  }
  
  out <- as.data.table(out)
  setnames(out, c("E", "L", "P", "A"))
  out[ , "run" := as.integer(runid)]
  out[ , "day" := 1:tmax]
  return(out)
})

out_sum <- do.call(rbind, out)
out_sum <- melt(out_sum, id.vars = c("run",'day'))

# plot
ggplot(out_sum[variable == "A", ]) +
  geom_line(aes(x=day,y=value,group=run),alpha=0.25) +
  theme_bw()

ggplot(data = out_sum) +
  geom_line(aes(x=day, y =value, group = run, color = variable), alpha = 0.35) +
  facet_wrap(. ~ variable, scales = "free") +
  theme_bw()
