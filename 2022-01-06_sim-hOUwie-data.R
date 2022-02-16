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

# a function which generates a parameter set given a continuous model
generateParameters <- function(continuous_model, alpha, sigma.sq, theta, discrete_model, rate){
  par_alpha <- alpha[continuous_model[1,]]
  par_alpha[is.na(par_alpha)] <- 1e-10
  par_sigma <- sigma.sq[continuous_model[2,] - min(continuous_model[2,]) + 1]
  par_sigma[is.na(par_sigma)] <- 1e-10
  par_theta <- theta[continuous_model[3,] - min(continuous_model[3,]) + 1]
  par_theta[is.na(par_theta)] <- 1e-10
  k.discrete <- max(discrete_model, na.rm=TRUE)
  par_rates <- rate[seq_len(k.discrete)]
  discrete_model[discrete_model == 0] <- NA
  discrete_model[is.na(discrete_model)] <- max(discrete_model, na.rm = TRUE) + 1
  Q <- matrix(0, dim(discrete_model)[1], dim(discrete_model)[1])
  Q[] <- c(par_rates, 0)[discrete_model]
  diag(Q) <- -rowSums(Q)
  p <- list(alpha = par_alpha, sigma.sq = par_sigma, theta = par_theta, Q = Q)
  return(p)
}


generateDataset <- function(nTip, continuous_model_cd, continuous_model_cid, discrete_model_cd, discrete_model_cid, alpha, sigma.sq, theta, rate, root.p){
  phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
  phy <- drop.extinct(phy)
  phy$edge.length <- phy$edge.length/max(branching.times(phy))
  # generate a dataset for CD
  cd_pars <- generateParameters(continuous_model_cd, alpha, sigma.sq, theta, discrete_model_cd, rate)
  dat_cd <- hOUwie.sim(phy, cd_pars$Q, root.p[1:2], cd_pars$alpha, cd_pars$sigma.sq, cd_pars$theta[1], cd_pars$theta)
  # ensure sampling of both observed states
  tip_count <- length(which(dat_cd$data[,2] == 2))
  while(!(tip_count > nTip*0.25 & tip_count < nTip*.75)){
    dat_cd <- hOUwie.sim(phy, cd_pars$Q, root.p[1:2], cd_pars$alpha, cd_pars$sigma.sq, cd_pars$theta[1], cd_pars$theta)
    tip_count <- length(which(dat_cd$data[,2] == 2))
  }
  # generate a dataset for CID
  cid_pars <- generateParameters(continuous_model_cid, alpha, sigma.sq, theta, discrete_model_cid, rate)
  dat_cid <- hOUwie.sim(phy, cid_pars$Q, root.p, cid_pars$alpha, cid_pars$sigma.sq, cid_pars$theta[1], cid_pars$theta)
  # ensure sampling of both hidden states
  tip_count_1b <- length(which(dat_cid$data[,2] == 3))
  tip_count_2b <- length(which(dat_cid$data[,2] == 4))
  tip_count  <- length(which(dat_cid$data[,2] == 3 | dat_cid$data[,2] == 4))
  while(!(tip_count > nTip*0.25 & tip_count < nTip*.75 & tip_count_1b > nTip*0.05 & tip_count_2b > nTip*0.05)){
    dat_cid <- hOUwie.sim(phy, cid_pars$Q, root.p, cid_pars$alpha, cid_pars$sigma.sq, cid_pars$theta[1], cid_pars$theta)
    tip_count_1b <- length(which(dat_cid$data[,2] == 3))
    tip_count_2b <- length(which(dat_cid$data[,2] == 4))
    tip_count <- length(which(dat_cid$data[,2] == 3 | dat_cid$data[,2] == 4))
  }
  out <- list(dat_cd=dat_cd, dat_cid=dat_cid)
  return(out)
}

save_datasets <- function(focal_model_type, model_names, all_model_structures, iter, nTip){
  dir_name <- paste0("simulated_data/", focal_model_type, "/", nTip, "/")
  dat_name <- paste0("dataset_", iter, ".Rsave")
  if(dat_name %in% dir(dir_name)){
    print("already a dataset with this name.")
    return(NULL)
  }else{
    print(paste0(dat_name, " is on its way..."))
  }
  focal_models <- sort(model_names[grep(paste0("_", focal_model_type, "$"), model_names)])
  continuous_model_cd <- all_model_structures[names(all_model_structures) == focal_models[1]][[1]]
  continuous_model_cid <- all_model_structures[names(all_model_structures) == focal_models[2]][[1]]
  focal_dat <- generateDataset(nTip, continuous_model_cd, continuous_model_cid, discrete_model_cd, discrete_model_cid, alpha, sigma.sq, theta, rate, root.p)
  file_name <- paste0(dir_name, dat_name)
  save(focal_dat, file = file_name)
}

#### #### #### #### #### #### #### #### #### #### #### #### 
# prerequisites
#### #### #### #### #### #### #### #### #### #### #### #### 
# maybe adjust houwie data
nTip <- 100
alpha <- c(3, 1.5)
sigma.sq <- c(0.35, 1)
theta <- c(2, 0.75)
root.p <- c(1,0,0,0)
theta0 <- theta[1]
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
model_names <- c("CID_BM1", "CD_BMV", "CID+_BMV", "CID_OU1", "CD_OUA", "CD_OUV", "CD_OUVA", "CD_OUM",
                 "CD_OUMA", "CD_OUMV", "CD_OUMVA", "CID+_OUA", "CID+_OUV", "CID+_OUVA", "CID+_OUM", "CID+_OUMA",
                 "CID+_OUMV", "CID+_OUMVA", "CD_OUBM1", "CD_OUBMV", "CID+_OUBM1", "CID+_OUBMV")
names(all_model_structures) <- model_names

# the 2 discrete models being evaluated
discrete_model_cd <- equateStateMatPars(getRateCatMat(2), 1:2)
discrete_model_cid <- getFullMat(list(discrete_model_cd, discrete_model_cd), getRateCatMat(2))
discrete_model_cid <- equateStateMatPars(discrete_model_cid, c(1,2,3,4))

#### #### #### #### #### #### #### #### #### #### #### #### 
# run
#### #### #### #### #### #### #### #### #### #### #### #### 

# models 1 and 4 are standard CID and should always be included
# models 2 and 3 pair as BMV
# models 5 and 12 pair as OUA
# models 6 and 13 pair as OUV
# models 7 and 14 pair as OUVA
# models 8 and 15 pair as OUM
# models 9 and 16 pair as OUMA
# models 10 and 17 pair as OUMV
# models 11 and 18 pair as OUMVA
# models 19 and 21 pair as OUBM1
# models 20 and 22 pair as OUBMV

# for each model type i want to generate a dataset consistent with the CD and one consistent with CID+
model_types <- c("BMV", "OUA", "OUV", "OUVA", "OUM", "OUMA", "OUMV", "OUMVA", "OUBM1", "OUBMV")
# model_types <- c("BMV", "OUV", "OUM", "OUMV")
# for each number of tips
ntips <- c(25, 100, 250)
# for a certain number of datasets
ntips <- 250
for(i in 1:length(model_types)){
  sapply(1:50, function(x) save_datasets(model_types[i], model_names, all_model_structures, x, ntips))
}

focal_model_type <- "OUA"
focal_models <- sort(model_names[grep(paste0("_", focal_model_type, "$"), model_names)])
continuous_model_cd <- all_model_structures[names(all_model_structures) == focal_models[1]][[1]]
continuous_model_cid <- all_model_structures[names(all_model_structures) == focal_models[2]][[1]]
focal_dat <- generateDataset(nTip, continuous_model_cd, continuous_model_cid, discrete_model_cd, discrete_model_cid, alpha, sigma.sq, theta, rate, root.p)



# 
# for(i in 1:length(model_types)){
#   print(i)
#   focal_model_type <- model_types[i]
#   focal_models <- sort(model_names[grep(paste0("_", focal_model_type, "$"), model_names)])
#   continuous_model_cd <- all_model_structures[names(all_model_structures) == focal_models[1]][[1]]
#   continuous_model_cid <- all_model_structures[names(all_model_structures) == focal_models[2]][[1]]
#   for(j in 1:length(ntips)){
#     nTip <- ntips[j]
#     focal_dat_list <- mclapply(1:iter, function(x) generateDataset(nTip, continuous_model_cd, continuous_model_cid, discrete_model_cd, discrete_model_cid, alpha, sigma.sq, theta, rate, root.p), mc.cores = 11)
#     for(k in 1:iter){
#       focal_dat <- focal_dat_list[[k]]
#       file_name <- paste0("simulated_data/", focal_model_type, "/", nTip, "/", "dataset_", k, ".Rsave")
#       save(focal_dat, file = file_name)
#     }
#   }
# }
# 
# # visualizing the data for each model type
# # model_type_list <- list()
# for(i in 1:length(model_types)){
#   focal_model_type <- model_types[i]
#   for(j in 1:length(ntips)){
#     nTip <- ntips[j]
#     # iter_list <- list()
#     for(k in 1:iter){
#       file_name <- paste0("simulated_data/", focal_model_type, "/", nTip, "/", "dataset_", k, ".Rsave")
#       load(file_name)
#       dat_cd <- focal_dat$dat_cd$data
#       dat_cid <- focal_dat$dat_cid$data
#       # iter_list <- c(iter_list, list(dat_cd=dat_cd, dat_cid=dat_cid))
#     }
#   }
# }
# 
# dat_cd <- focal_dat$dat_cd$data
# dat_cid <- focal_dat$dat_cid$data
# 
# 
# 
# continuous_model_cd <- all_model_structures[[8]]
# continuous_model_cid <- all_model_structures[[15]]
# OU1_classic <- all_model_structures[[4]]
# 
# out <- list()
# map_number <- c(100, 250, 500, 1000)
# for(i in 1:length(map_number)){
#   nmaps <- map_number[i]
#   out[[i]] <- mclapply(1:20, function(x) run_hmm_help(nmaps), mc.cores = 5)
# }

# nmap_100 <- mclapply(1:20, function(x) run_hmm_help(100), mc.cores = 5)
# nmap_200 <- mclapply(1:20, function(x) run_hmm_help(200), mc.cores = 5)
# nmap_1000 <- mclapply(1:20, function(x) run_hmm_help(1000), mc.cores = 5)
# nmap_10000 <- mclapply(1:20, function(x) run_hmm_help(10000), mc.cores = 20)
