setwd("~/2020_hOUwie/")

source("hOUwieSimmap.R")
source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)
require(partitions)

load("sim_fits/cd-fit-generating_model_1.Rsave")
load("sim_data/cd-data-generating_model_1.Rsave")

full_data <- all_data[[11]]$pars[1]
model <- out[[11]][1]
full_data$pars
model[[1]]$RegimeMap$maps

out[[1]][[1]]$p
simulating_transition_rate <- unlist(lapply(all_data, function(x) x$pars[1]))
estimated_transition_rate <- unlist(lapply(out, function(x) x[[1]]$p[1]))

plot(x = simulating_transition_rate, y = estimated_transition_rate, xlim=c(0.2,0.8), ylim=c(0.2,0.8))
abline(lm(estimated_transition_rate~simulating_transition_rate))
abline(a = 0, b = 1, col = "red")

simulating_transition_number <- unlist(lapply(all_data, function(x)
  sum(unlist(lapply(x$simmap[[1]]$maps, length)) - 1)))
estimated_transition_number <- unlist(lapply(out, function(x) 
  sum(unlist(lapply(x[[1]]$RegimeMap$maps, length)) - 1)))

plot(x = simulating_transition_number, y = estimated_transition_number)
abline(lm(estimated_transition_number~simulating_transition_number))
abline(a = 0, b = 1, col = "red")
