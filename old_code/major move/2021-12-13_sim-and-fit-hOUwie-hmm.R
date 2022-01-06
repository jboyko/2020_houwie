setwd("~/2020_hOUwie/")

source("hOUwieNode.R")
source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)
require(partitions)

#### #### #### #### #### #### #### #### #### #### #### #### 
# prerequisites
#### #### #### #### #### #### #### #### #### #### #### #### 

nTip <- 200
alpha <- c(5,5)
sigma2 <- alpha/5
theta <- c(12,24)
rate <- .1

#### #### #### #### #### #### #### #### #### #### #### #### 
# the 22 2-state models we want to fit - put into list form
#### #### #### #### #### #### #### #### #### #### #### #### 

continuous_models_cd_ou <- getAllContinuousModelStructures(2, "OU")
continuous_models_cd_bm <- getAllContinuousModelStructures(2, "BM")
continuous_models_cd_bmou <- getAllContinuousModelStructures(2, "BMOU")
continuous_models_cid_ou <- continuous_models_cd_ou[,c(1,1,2,2),]
continuous_models_cid_bm <- continuous_models_cd_bm[,c(1,1,2,2),]
continuous_models_cid_bmou <- continuous_models_cd_bmou[,c(1,1,2,2),]
continuous_models_cd_ou <- lapply(seq(dim(continuous_models_cd_ou)[3]), function(x) continuous_models_cd_ou[,,x])
continuous_models_cd_bm <- lapply(seq(dim(continuous_models_cd_bm)[3]), function(x) continuous_models_cd_bm[,,x])
continuous_models_cd_bmou <- lapply(seq(dim(continuous_models_cd_bmou)[3]), function(x) continuous_models_cd_bmou[,,x])
continuous_models_cid_ou <- lapply(seq(dim(continuous_models_cid_ou)[3]), function(x) continuous_models_cid_ou[,,x])
continuous_models_cid_bm <- lapply(seq(dim(continuous_models_cid_bm)[3]), function(x) continuous_models_cid_bm[,,x])
continuous_models_cid_bmou <- lapply(1:dim(continuous_models_cid_bmou)[3], function(x) continuous_models_cid_bmou[,,x])
all_model_structures <- c(continuous_models_cd_bm, continuous_models_cid_bm[-1], continuous_models_cd_ou, continuous_models_cid_ou[-1], continuous_models_cd_bmou, continuous_models_cid_bmou)

# the 2 discrete models being evaluated
discrete_model_cd <- equateStateMatPars(getRateCatMat(2), 1:2)
discrete_model_cid <- getFullMat(list(discrete_model_cd, discrete_model_cd), getRateCatMat(2))
discrete_model_cid <- equateStateMatPars(discrete_model_cid, c(1,2,3,4))

#### #### #### #### #### #### #### #### #### #### #### #### 
# Functions
#### #### #### #### #### #### #### #### #### #### #### #### 

# a function which generates a parameter set given a continuous model
generateParameters <- function(continuous_model, alpha, sigma2, theta, discrete_model, rate){
  k.alpha <- length(unique(na.omit(continuous_model[1,])))
  k.sigma <- length(unique(na.omit(continuous_model[2,])))
  k.theta <- length(unique(na.omit(continuous_model[3,])))
  k.discrete <- max(discrete_model, na.rm=TRUE)
  # par_alpha <- runif(k.alpha, minAlpha, maxAlpha)
  # par_sigma <- runif(k.sigma, minSigma2, maxSigma2)
  # par_theta <- runif(k.theta, minTheta, maxTheta)
  par_alpha <- alpha[1:k.alpha]
  par_sigma <- sigma2[1:k.sigma]
  par_theta <- theta[1:k.theta]
  par_rates <- rate[1:k.discrete]
  p <- c(par_rates, par_alpha, par_sigma, par_theta)
  return(p)
}

# a function to fit a single continuous model given the full data
  singleFit <- function(full_data, phy, continuous_model, discrete_model_cd, discrete_model_cid, houwie_parameters){
    # make hidden states observed
    time_slice <- houwie_parameters[1,1]
    nSim <- houwie_parameters[1,2]
    data <- full_data$data
    data[data[,2]==3,2] <- 1
    data[data[,2]==4,2] <- 2
    # determine if the model includes hidden states or not
    if(dim(continuous_model)[2] == 2){
      discrete_model <- discrete_model_cd
      rate.cat <- 1
    }else{
      discrete_model <- discrete_model_cid
      rate.cat <- 2
    }
    # fit the houwie model
    fit <- hOUwie(phy = phy, data = data, rate.cat = rate.cat, nSim = nSim, time_slice = time_slice, discrete_model = discrete_model, continuous_model = continuous_model, recon = FALSE, sample_tips = FALSE)
    return(fit)
  }

makeData <- function(nTip, continuous_model, discrete_model, alpha, sigma2, theta, rate){
  phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
  phy <- drop.extinct(phy)
  phy$edge.length <- phy$edge.length/max(branching.times(phy))
  pars <- generateParameters(continuous_model, alpha, sigma2, theta, discrete_model, rate)
  full_data <- try(generateData(phy, discrete_model, continuous_model, pars))
  while(class(full_data) == "try-error"){
    full_data <- try(generateData(phy, discrete_model, continuous_model, pars))
  }
  return(full_data)
}

# a single iteration
singleRun <- function(i, iter, nSim){
  model_name <- paste0("M", i)
  cat("Begining", model_name, "...\n")
  # deciding what tree size we would like to use
  phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
  phy <- drop.extinct(phy)
  phy$edge.length <- phy$edge.length/max(branching.times(phy))
  # generate data
  generating_model <- all_model_structures[[i]]
  discrete_model <- list(discrete_model_cd, discrete_model_cid)[[ifelse(dim(generating_model)[2] == 2, 1, 2)]]
  discrete_model[is.na(discrete_model)] <- 0
  pars <- generateParameters(generating_model, alpha, sigma2, theta, discrete_model, rate)
  full_data <- try(generateData(phy, discrete_model, generating_model, pars, root.p = "1"))
  while(class(full_data) == "try-error"){
    full_data <- try(generateData(phy, discrete_model, generating_model, pars, root.p = "1"))
  }
  # mclapply over all model structures
  houwie_parameters <- matrix(c(1.1, nSim), 1, 2)
  model_res <- mclapply(all_model_structures, function(x) singleFit(full_data, phy, x, discrete_model_cd, discrete_model_cid, houwie_parameters), mc.cores = 1)
  names(model_res) <- unlist(lapply(houwie_parameters, function(x) paste(x, collapse = "_")))
  file_name_res <- paste0("hmm_help_houwie/", model_name, "-", nSim, "_", iter, ".Rsave")
  out <- list(simulated_data = full_data, model_res = model_res)
  save(out, file = file_name_res)
}


#### #### #### #### #### #### #### #### #### #### #### #### 
# run
#### #### #### #### #### #### #### #### #### #### #### #### 
hidden_state_OUM <- all_model_structures[[15]]
OUM_classic <- all_model_structures[[8]]
OU1_classic <- all_model_structures[[4]]

run_hmm_help <- function(nSim=100){
  sim_dat <- makeData(nTip, hidden_state_OUM, discrete_model_cid, alpha, sigma2, theta, rate)
  R2_count <- length(which(sim_dat$data[,2] == 3 | sim_dat$data[,2] == 4))
  while(R2_count > nTip*0.3 & R2_count < nTip*.7){
    sim_dat <- makeData(nTip, hidden_state_OUM, discrete_model_cid, alpha, sigma2, theta, rate)
    R2_count <- length(which(sim_dat$data[,2] == 3 | sim_dat$data[,2] == 4))
  }
  houwie_dat <- sim_dat$data
  houwie_dat[,2][houwie_dat[,2] == 3] <- 1
  houwie_dat[,2][houwie_dat[,2] == 4] <- 2
  phy <- sim_dat$simmap[[1]]
  print("running hidden state")
  out_1 <- hOUwie(phy = phy, data = houwie_dat, rate.cat = 2, nSim = nSim, time_slice = 1.1, discrete_model = discrete_model_cid, continuous_model = hidden_state_OUM, recon = FALSE, sample_tips = TRUE, sample_nodes = TRUE, optimizer = "nlopt_ln", n_starts = 10, ncores = 10)
  print("running observed state")
  out_2 <- hOUwie(phy = phy, data = houwie_dat, rate.cat = 1, nSim = nSim, time_slice = 1.1, discrete_model = discrete_model_cd, continuous_model = OUM_classic, recon = FALSE, sample_tips = TRUE, sample_nodes = TRUE, optimizer = "nlopt_ln", n_starts = 10, ncores = 10)
  print("running OU1")
  out_3 <- hOUwie(phy = phy, data = houwie_dat, rate.cat = 1, nSim = nSim, time_slice = 1.1, discrete_model = discrete_model_cd, continuous_model = OU1_classic, recon = FALSE, sample_tips = TRUE, sample_nodes = TRUE, optimizer = "nlopt_ln", n_starts = 10, ncores = 10)
  out <- list(OU1=out_3, CD=out_2, CID=out_1)
  return(out)
}

out <- list()
map_number <- c(100, 250, 500, 1000)
for(i in 1:length(map_number)){
  nmaps <- map_number[i]
  out[[i]] <- mclapply(1:20, function(x) run_hmm_help(nmaps), mc.cores = 5)
}

# nmap_100 <- mclapply(1:20, function(x) run_hmm_help(100), mc.cores = 5)
# nmap_200 <- mclapply(1:20, function(x) run_hmm_help(200), mc.cores = 5)
# nmap_1000 <- mclapply(1:20, function(x) run_hmm_help(1000), mc.cores = 5)
# nmap_10000 <- mclapply(1:20, function(x) run_hmm_help(10000), mc.cores = 20)
