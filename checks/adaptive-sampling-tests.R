#### #### #### #### #### #### #### #### #### #### #### #### 
# imports
#### #### #### #### #### #### #### #### #### #### #### #### 

setwd("~/2020_hOUwie/")

source("hOUwieNode.R")
require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)
require(partitions)

#### #### #### #### #### #### #### #### #### #### #### #### 
# Functions
#### #### #### #### #### #### #### #### #### #### #### #### 


#### #### #### #### #### #### #### #### #### #### #### #### 
# run
#### #### #### #### #### #### #### #### #### #### #### #### 

phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = 10) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

# generate a dataset for CD
Q <- matrix(c(-1,1,1,-1), 2, 2)
root.p <- c(1,0)
alpha <- c(1e-10,1e-10)
sigma.sq <- c(0.1, 10)
theta <- c(10, 10)
dat_cd <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta[1], theta)


# cd_fit <- hOUwie(phy = phy, data = dat_cd$data, rate.cat = 1, nSim = 25, time_slice = 1.1, discrete_model = "ER", continuous_model = "OUM", recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = FALSE, optimizer = "nlopt_ln")
# p <- c(6.751825e+00, 1.232609e-04, 6.931472e-01, 1.025741e+01, 1.288142e+01)
cd_fit <- hOUwie(phy = phy, data = dat_cd$data, rate.cat = 1, nSim = 10, time_slice = 1.1, discrete_model = "ER", continuous_model = "BMV", recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln")
bm1_fit <- hOUwie(phy = phy, data = dat_cd$data, rate.cat = 1, nSim = 10, time_slice = 1.1, discrete_model = "ER", continuous_model = "BM1", recon = FALSE, sample_tips = FALSE, sample_nodes = FALSE, adaptive_sampling = FALSE, optimizer = "nlopt_ln")


best_fit_mapping <- correct_map_edges(cd_fit$simmaps[[which.max(cd_fit$all_disc_liks + cd_fit$all_cont_liks)]])
true_mapping <- dat_cd$simmap
par(mfrow=c(1,2))
plot(true_mapping)
plot(best_fit_mapping)

simmap <- best_fit_mapping
Rate.mat <- rbind(alpha, sigma.sq, theta)
# cd_fit$expected_vals
OUwie.basic(simmap, dat_cd$data, TRUE, alpha = alpha, sigma.sq = sigma.sq, theta = theta, return.expected.vals = TRUE)

test_expect <- getOUExpectations(simmap, Rate.mat)
tip_value <- dat_cd$data$x[8]
expected_value <- test_expect$expected_means[8]
expected_variance <- test_expect$expected_variances[8]
dnorm(tip_value, expected_value, sqrt(expected_variance), log = TRUE)

c(dOU(15, 10, 1, 1, 10, 1, log = TRUE), dOU(15, 20, 1, 1, 20, 1, log = TRUE))

logstuff <- c(dOU(14.8, 10, 1, 1, 10, 4, log = TRUE), dOU(14.8, 20, 1, 1, 20, 4, log = TRUE))
exp(logstuff)/sum(exp(logstuff))


