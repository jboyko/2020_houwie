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
require(viridis)


#### #### #### #### #### #### #### #### #### #### #### #### 
# dataset
#### #### #### #### #### #### #### #### #### #### #### #### 

d <- getData("empirical/Ericaceae_niche.csv")
t <- read.tree("empirical/Ericaceae_Schwery_etal_2015.tre")

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
discrete_model_cid <- equateStateMatPars(discrete_model_cid, c(1,2))

data <- d$dat.arid.se
phy <- keep.tip(t, data$sp)

time_slice <- c(20, 60, 120)
nSim <- c(50, 100, 200)
houwie_parameters <- expand.grid(time_slice, nSim)
houwie_parameters[1,1]

eric_results <- vector("list", dim(houwie_parameters)[1])
for(i in 1:dim(houwie_parameters)[1]){
  eric_results[[i]] <- mclapply(all_model_structures, function(x) quickRun(data, phy, x, discrete_model_cid, discrete_model_cd, nSim = houwie_parameters[i,2], time_slice = houwie_parameters[i,1]), mc.cores = 22)
}

# evaluating hybrid models
eric_results_hmm <- vector("list", dim(houwie_parameters)[1])
hybrid_model_structures <- list(getOUParamStructure("BMS", "three.point", FALSE, FALSE, 4),
                                getOUParamStructure("OUM", "three.point", FALSE, FALSE, 4),
                                getOUParamStructure("OUMV", "three.point", FALSE, FALSE, 4),
                                getOUParamStructure("OUMA", "three.point", FALSE, FALSE, 4),
                                getOUParamStructure("OUMVA", "three.point", FALSE, FALSE, 4))

eric_results_hmm <- mclapply(hybrid_model_structures, function(x) quickRun(data, phy, x, discrete_model_cid, discrete_model_cd, nSim = 10, time_slice = 120, n_starts=10), mc.cores = 5)

# for(i in 1:dim(houwie_parameters)[1]){
#   eric_results_hmm[[i]] <- mclapply(hybrid_model_structures, function(x) quickRun(data, phy, x, discrete_model_cd_hmm, discrete_model_cd, nSim = houwie_parameters[i,2], time_slice = houwie_parameters[i,1]), mc.cores = 22)
# }

# fit <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = 10, time_slice = 120, discrete_model = discrete_model_cid, continuous_model = hybrid_model_structures[[2]], recon = FALSE, mserr = "known")





