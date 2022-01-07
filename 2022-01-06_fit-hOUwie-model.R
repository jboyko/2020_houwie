#### #### #### #### #### #### #### #### #### #### #### #### 
# imports
#### #### #### #### #### #### #### #### #### #### #### #### 

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
# Functions
#### #### #### #### #### #### #### #### #### #### #### #### 

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

run_hmm_help <- function(nSim=100){
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
# for each model type i want to generate a dataset consistent with the CD and one consistent with CID+
model_types <- c("BMV", "OUA", "OUV", "OUVA", "OUM", "OUMA", "OUMV", "OUMVA", "OUBM1", "OUBMV")
# for each number of tips
ntips <- c(25, 100, 250)

# the dataests come from a particular structure and modeling results will also follow the same saving protocol
for(i in 1:length(model_types)){
  focal_model_type <- model_types[i]
  focal_models <- sort(model_names[grep(paste0("_", focal_model_type, "$"), model_names)])
  continuous_model_cd <- all_model_structures[names(all_model_structures) == focal_models[1]][[1]]
  continuous_model_cid <- all_model_structures[names(all_model_structures) == focal_models[2]][[1]]
  for(j in 1:length(ntips)){
    nTip <- ntips[j]
    focal_dir <- paste0("simulated_data/", focal_model_type, "/", nTip)
    dataset_files <- dir(focal_dir, full.names = TRUE)
  }
}

dataset_file <- dataset_files[1]
load(dataset_file)

phy <- focal_dat$dat_cd$simmap[[1]]
dat_cd <- focal_dat$dat_cd$data
dat_cid <- focal_dat$dat_cid$data
dat_cid[dat_cid[,2] == 3,2] <- 1
dat_cid[dat_cid[,2] == 4,2] <- 2

tmp <- hOUwie(phy = phy, data = dat_cid, rate.cat = 2, nSim = 10, time_slice = 1.1, discrete_model = discrete_model_cid, continuous_model = continuous_model_cid, recon = FALSE, sample_tips = TRUE, sample_nodes = TRUE, optimizer = "nlopt_ln")

focal_dat$dat_cd$index.cor


