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

singleRun <- function(dataset_file, nmap=25, continuous_model_cd, continuous_model_cid, n_starts=5, n_cores=5){
  load(dataset_file)
  phy <- focal_dat$dat_cd$simmap
  dat_cd <- focal_dat$dat_cd$data
  dat_cid <- focal_dat$dat_cid$data
  dat_cid[dat_cid[,2] == 3,2] <- 1
  dat_cid[dat_cid[,2] == 4,2] <- 2
  ## CD data
  # run a bm1
  bm1_fit <- hOUwie(phy = phy, data = dat_cd, rate.cat = 1, nSim = nmap, time_slice = 1.1, discrete_model = discrete_model_cd, continuous_model = bm1_model, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln", n_starts = n_starts, ncores = n_cores)
  # run an ou1
  ou1_fit <- hOUwie(phy = phy, data = dat_cd, rate.cat = 1, nSim = nmap, time_slice = 1.1, discrete_model = discrete_model_cd, continuous_model = ou1_model, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln", n_starts = n_starts, ncores = n_cores)
  # run a cd 
  cd_fit <- hOUwie(phy = phy, data = dat_cd, rate.cat = 1, nSim = nmap, time_slice = 1.1, discrete_model = discrete_model_cd, continuous_model = continuous_model_cd, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln", n_starts = n_starts, ncores = n_cores)
  # run a cid + 
  cid_fit <- hOUwie(phy = phy, data = dat_cd, rate.cat = 2, nSim = nmap, time_slice = 1.1, discrete_model = discrete_model_cid, continuous_model = continuous_model_cid, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln", n_starts = n_starts, ncores = n_cores)
  cd_out <- list(bm1_fit=bm1_fit, ou1_fit=ou1_fit, cd_fit=cd_fit, cid_fit=cid_fit)
  
  ## CID data
  # run a bm1
  bm1_fit <- hOUwie(phy = phy, data = dat_cid, rate.cat = 1, nSim = nmap, time_slice = 1.1, discrete_model = discrete_model_cd, continuous_model = bm1_model, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln", n_starts = n_starts, ncores = n_cores)
  # run an ou1
  ou1_fit <- hOUwie(phy = phy, data = dat_cid, rate.cat = 1, nSim = nmap, time_slice = 1.1, discrete_model = discrete_model_cd, continuous_model = ou1_model, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln", n_starts = n_starts, ncores = n_cores)
  # run a cd 
  cd_fit <- hOUwie(phy = phy, data = dat_cid, rate.cat = 1, nSim = nmap, time_slice = 1.1, discrete_model = discrete_model_cd, continuous_model = continuous_model_cd, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln", n_starts = n_starts, ncores = n_cores)
  # run a cid + 
  cid_fit <- hOUwie(phy = phy, data = dat_cid, rate.cat = 2, nSim = nmap, time_slice = 1.1, discrete_model = discrete_model_cid, continuous_model = continuous_model_cid, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln", n_starts = n_starts, ncores = n_cores)
  cid_out <- list(bm1_fit=bm1_fit, ou1_fit=ou1_fit, cd_fit=cd_fit, cid_fit=cid_fit)
  
  out <- list(cd_out = cd_out, cid_out = cid_out)
  modelres_file <- gsub("simulated_data/", "simulated_fit/", dataset_file)
  modelres_file <- gsub(".Rsave", paste0("_nmap=", nmap, ".Rsave"), modelres_file)
  save(out, file = modelres_file)
  # return(out)
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
model_types <- c("BMV","OUV", "OUM", "OUMV")
# model_types <- c("OUM")

# for each number of tips
ntips <- c(25, 100, 250)

bm1_model <- all_model_structures[[1]]
ou1_model <- all_model_structures[[4]]

# the dataests come from a particular structure and modeling results will also follow the same saving protocol
for(j in 1:length(ntips)){
  for(i in 1:length(model_types)){
    focal_model_type <- model_types[i]
    focal_models <- sort(model_names[grep(paste0("_", focal_model_type, "$"), model_names)])
    continuous_model_cd <- all_model_structures[names(all_model_structures) == focal_models[1]][[1]]
    continuous_model_cid <- all_model_structures[names(all_model_structures) == focal_models[2]][[1]]
    nTip <- ntips[j]
    focal_dir <- paste0("simulated_data/", focal_model_type, "/", nTip)
    dataset_files <- dir(focal_dir, full.names = TRUE)
    mclapply(dataset_files, function(x) singleRun(x, 250, continuous_model_cd, continuous_model_cid, n_starts = 5, n_cores = 5), mc.cores = 11)
  }
}



