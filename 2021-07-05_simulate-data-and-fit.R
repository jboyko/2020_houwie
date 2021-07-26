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

nCores <- 24
nTip <- 100
nSim <- 25
minAlpha = 0.5
maxAlpha = 4
minSigma2 = 0.5
maxSigma2 = 2
minTheta = 10
maxTheta = 20
minRate = 0.1
maxRate = 1

#### #### #### #### #### #### #### #### #### #### #### #### 
# the 40 2-state models we want to fit - put into list form
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
all_model_structures <- c(continuous_models_cd_bm, continuous_models_cd_bmou, continuous_models_cd_ou, continuous_models_cid_bm, continuous_models_cid_bmou, continuous_models_cid_ou)

# the 2 discrete models being evaluated
discrete_model_cd <- equateStateMatPars(getRateCatMat(2), 1:2)
discrete_model_cid <- getFullMat(list(discrete_model_cd, discrete_model_cd), discrete_model_cd)
discrete_model_cid <- equateStateMatPars(discrete_model_cid, c(1,2))

# deciding what tree size we would like to use
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

#### #### #### #### #### #### #### #### #### #### #### #### 
# Functions
#### #### #### #### #### #### #### #### #### #### #### #### 

# a function which generates a parameter set given a continuous model
generateParameters <- function(continuous_model, minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta, discrete_model, minRate, maxRate){
  k.alpha <- length(unique(na.omit(continuous_model[1,])))
  k.sigma <- length(unique(na.omit(continuous_model[2,])))
  k.theta <- length(unique(na.omit(continuous_model[3,])))
  k.discrete <- max(na.omit(discrete_model))
  par_alpha <- runif(k.alpha, minAlpha, maxAlpha)
  par_sigma <- runif(k.sigma, minSigma2, maxSigma2)
  par_theta <- runif(k.theta, minTheta, maxTheta)
  par_rates <- runif(k.discrete, minRate, maxRate)
  p <- c(par_rates, par_alpha, par_sigma, par_theta)
  return(p)
}

# a function to fit a single continuous model given the full data
singleFit <- function(full_data, continuous_model, discrete_model_cd, discrete_model_cid, nSim){
  # make hidden states observed
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
  fit <- hOUwie(phy = phy, data = data, rate.cat = rate.cat, nSim = nSim, discrete_model = discrete_model, continuous_model = continuous_model, recon = TRUE, nodes = "all")
  return(fit)
}


  

#### #### #### #### #### #### #### #### #### #### #### #### 
# CD models generate
#### #### #### #### #### #### #### #### #### #### #### #### 
for(iter in 1:5){
for(i in 1:12){
  model_name <- paste0("M", i)
  cat("Begining", model_name, "...\n")
  # generate data
  generating_model <- all_model_structures[[i]]
  discrete_model <- list(discrete_model_cd, discrete_model_cid)[[ifelse(dim(generating_model)[2] == 2, 1, 2)]]
  pars <- generateParameters(generating_model, minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta, discrete_model, minRate, maxRate)
  full_data <- try(generateData(phy, discrete_model, generating_model, pars))
  while(class(full_data) == "try-error"){
    full_data <- try(generateData(phy, discrete_model, generating_model, pars))
  }
  file_name_data <- paste0("sim_data/", model_name, "-gen_", nTip, "-nTip_", nSim, "-nMap_", iter, "-iter.Rsave")
  save(full_data, file = file_name_data)
  # mclapply over all model structures
  out <- mclapply(all_model_structures, function(x) singleFit(full_data, x, discrete_model_cd, discrete_model_cid, nSim), mc.cores = nCores)
  # tmp <- singleFit(full_data, all_model_structures[[1]], discrete_model_cd, discrete_model_cid, 10)
  # p = c(0.483848937,  0.007448614,  0.693147181,  0.284018414, 15.829464035)
  file_name_res <- paste0("sim_fits/", model_name, "-gen_", nTip, "-nTip_", nSim, "-nMap_", iter, "-iter.Rsave")
  save(out, file = file_name_res)
  full_data <- NULL
}
}
#### #### #### #### #### #### #### #### #### #### #### #### 
# CID models generate
#### #### #### #### #### #### #### #### #### #### #### #### 
for(iter in 1:5){
for(i in 13:24){
  model_name <- paste0("M", i)
  cat("Begining", model_name, "...\n")
  # generate data
  generating_model <- all_model_structures[[i]]
  discrete_model <- list(discrete_model_cd, discrete_model_cid)[[ifelse(dim(generating_model)[2] == 2, 1, 2)]]
  pars <- generateParameters(generating_model, minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta, discrete_model, minRate, maxRate)
  full_data <- try(generateData(phy, discrete_model, generating_model, pars))
  while(class(full_data) == "try-error"){
    full_data <- try(generateData(phy, discrete_model, generating_model, pars))
  }
  file_name_data <- paste0("sim_data/", model_name, "-gen_", nTip, "-nTip_", nSim, "-nMap_", iter, "-iter.Rsave")
  save(full_data, file = file_name_data)
  # mclapply over all model structures
  out <- mclapply(all_model_structures, function(x) singleFit(full_data, x, discrete_model_cd, discrete_model_cid, nSim), mc.cores = nCores)
  file_name_res <- paste0("sim_fits/", model_name, "-gen_", nTip, "-nTip_", nSim, "-nMap_", iter, "-iter.Rsave")
  save(out, file = file_name_res)
  full_data <- NULL
}
}
#### #### #### #### #### #### #### #### #### #### #### #### 
# both models generate
#### #### #### #### #### #### #### #### #### #### #### #### 
for(iter in 1:5){
for(i in 1:24){
  model_name <- paste0("M", i)
  cat("Begining", model_name, "...\n")
  # generate data
  generating_model <- all_model_structures[[i]]
  discrete_model <- list(discrete_model_cd, discrete_model_cid)[[ifelse(dim(generating_model)[2] == 2, 1, 2)]]
  pars <- generateParameters(generating_model, minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta, discrete_model, minRate, maxRate)
  full_data <- try(generateData(phy, discrete_model, generating_model, pars))
  while(class(full_data) == "try-error"){
    full_data <- try(generateData(phy, discrete_model, generating_model, pars))
  }
  file_name_data <- paste0("sim_data/", model_name, "-gen_", nTip, "-nTip_", nSim, "-nMap_", iter, "-iter.Rsave")
  save(full_data, file = file_name_data)
  # mclapply over all model structures
  out <- mclapply(all_model_structures, function(x) singleFit(full_data, x, discrete_model_cd, discrete_model_cid, nSim), mc.cores = nCores)
  file_name_res <- paste0("sim_fits/", model_name, "-gen_", nTip, "-nTip_", nSim, "-nMap_", iter, "-iter.Rsave")
  save(out, file = file_name_res)
  full_data <- NULL
}
}
#### #### #### #### #### #### #### #### #### #### #### #### 
# ASR was still in log lik mode - this is now corrected, but have fix here
#### #### #### #### #### #### #### #### #### #### #### #### 
# 
# setwd("~/2020_hOUwie/")
# source("hOUwieNode.R")
# 
# files_fits <- dir("sim_fits/", full.names = TRUE)
# data_fits <-  dir("data_fits/", full.names = TRUE)
# 
# for(i in 1:length(files_fits)){
#   print(i)
#   load(files_fits[i])
#   for(j in 1:length(out)){
#     if(class(out[[j]]) != "houwie"){
#       next
#     }
#     if(mean(rowSums(out[[j]]$recon)) == 1){
#       next
#     }
#     out[[j]]$recon[out[[j]]$recon == 0] <- -Inf
#     tmp <- t(apply(out[[j]]$recon, 1, function(x) x - max(x)))
#     out[[j]]$recon <- round(exp(tmp)/rowSums(exp(tmp)), 10)
#   }
#   save(out, file = files_fits[i])
# }
# 
# undebug(hOUwieRecon)
# test <- hOUwieRecon(out[[j]])


