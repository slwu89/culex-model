library(Rcpp)
library(RcppArmadillo)
library(deSolve)

library(parallel)
library(data.table)
library(ggplot2)
library(here)

rm(list=ls());gc()

Rcpp::sourceCpp(here("src_R/culex.cpp"))
source(here("src_R/tau_ode.R"))

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
A0 <- c(600, 0)

# movement between 2 patches
psi <- matrix(
  data = c(0.95, 0.05,
           0.025, 0.975),
  nrow = 2, ncol = 2, byrow = TRUE
)

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

# solve a determiinistic trajectory
mod <- create_culex_deterministic(p = 2, tau_E = tauE, tau_L = tauL, tau_P = tauP, dt = dt, psi = psi)
set_A_deterministic(mod = mod, A = A0)

out_d <- data.table(day = rep(1:tmax, each = 8), stage = rep(c("E","L","P","A"), tmax*2), patch = rep(c(1,1,1,1,2,2,2,2), tmax), value = NaN)
setkey(out_d, "day")
out_d[, "patch" := as.factor(patch)]
out_i <- 1L

for (i in 1:maxstep) {
  step_culex_deterministic(mod = mod, parameters = parameters)
  if ((i-1) %% (1/dt) == 0) {
    out_d[day == out_i & stage == "E", "value" := as.vector(get_E_deterministic(mod))]
    out_d[day == out_i & stage == "L", "value" := as.vector(get_L_deterministic(mod))]
    out_d[day == out_i & stage == "P", "value" := as.vector(get_P_deterministic(mod))]
    out_d[day == out_i & stage == "A", "value" := as.vector(get_A_deterministic(mod))]
    out_i <- out_i + 1L
  }
}

# draw samples from stochastic model
out <- parallel::mclapply(X = 1:10, FUN = function(runid) {
  mod <- create_culex_stochastic(p = 2, tau_E = tauE, tau_L = tauL, tau_P = tauP, dt = dt, psi = psi)
  set_A_stochastic(mod = mod, A = A0)

  out <- data.table(day = rep(1:tmax, each = 8), stage = rep(c("E","L","P","A"), tmax*2), patch = rep(c(1,1,1,1,2,2,2,2), tmax), value = NaN)
  setkey(out, "day")
  out[, "patch" := as.factor(patch)]
  out_i <- 1L

  for (i in 1:maxstep) {
    step_culex_stochastic(mod = mod, parameters = parameters)
    if ((i-1) %% (1/dt) == 0) {
      out[day == out_i & stage == "E", "value" := as.vector(get_E_stochastic(mod))]
      out[day == out_i & stage == "L", "value" := as.vector(get_L_stochastic(mod))]
      out[day == out_i & stage == "P", "value" := as.vector(get_P_stochastic(mod))]
      out[day == out_i & stage == "A", "value" := as.vector(get_A_stochastic(mod))]
      out_i <- out_i + 1L
    }
  }

  out[ , "run" := as.integer(runid)]
  return(out)
})


out_sum <- do.call(rbind, out)

ggplot(data = out_sum) +
  geom_line(aes(x=day, y =value, group = interaction(run, patch), color = stage, linetype = patch), alpha = 0.35) +
  geom_line(data = out_d, aes(x=day, y =value, color = stage, linetype = patch)) +
  facet_wrap(. ~ stage, scales = "free") +
  theme_bw()

ggplot(data = out_sum[stage == "A", ]) +
  geom_line(aes(x=day, y =value, group = interaction(run, patch), color = stage, linetype = patch), alpha = 0.35) +
  geom_line(data = out_d[stage == "A", ], aes(x=day, y =value, color = stage, linetype = patch)) +
  theme_bw()
