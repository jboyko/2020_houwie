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
hybrid_model_structures <- list(getOUParamStructure("BMS", "three.point", FALSE, FALSE, 4),
                                getOUParamStructure("OUM", "three.point", FALSE, FALSE, 4),
                                getOUParamStructure("OUMV", "three.point", FALSE, FALSE, 4),
                                getOUParamStructure("OUMA", "three.point", FALSE, FALSE, 4),
                                getOUParamStructure("OUMVA", "three.point", FALSE, FALSE, 4))

all_model_structures <- c(all_model_structures, hybrid_model_structures)
# the 2 discrete models being evaluated
discrete_model_cd <- getRateCatMat(2)
discrete_model_cd_2 <- equateStateMatPars(getRateCatMat(2), c(1,2))
discrete_model_cid <- getFullMat(list(discrete_model_cd, discrete_model_cd), getRateCatMat(2))
discrete_model_cid <- equateStateMatPars(discrete_model_cid, c(1,2,3,4))

## function for running houwie paralellized
quickRun <- function(data, phy, continuous_model, discrete_model_cid, discrete_model_cd, nSim=50, time_slice=20){
  if(dim(continuous_model)[2] == 2){
    discrete_model <- discrete_model_cd
    rate.cat <- 1
  }else{
    discrete_model <- discrete_model_cid
    rate.cat <- 2
  }
  if(dim(data)[2] == 4){
    mserr <- "known"
  }else{
    mserr <- "none"
  }
  fit <- try(hOUwie(phy = phy, data = data, rate.cat = rate.cat, nSim = nSim, time_slice = time_slice, discrete_model = discrete_model, continuous_model = continuous_model, recon = FALSE, mserr = mserr, n_starts = 10, optimizer="nlopt_ln", sample_tips = FALSE))
  return(fit)
}


data <- d$dat.arid.se
phy <- keep.tip(t, data$sp)

eric_results <- mclapply(all_model_structures, function(x) quickRun(data, phy, x, discrete_model_cid, discrete_model_cd, nSim = 50, time_slice = 120), mc.cores = 27)

lapply(all_model_structures, function(x) dim(x)[2] < 4)

eric_results_er <- mclapply(all_model_structures[unlist(lapply(all_model_structures, function(x) dim(x)[2] < 4))], function(x) quickRun(data, phy, x, discrete_model_cid, discrete_model_cd_2, nSim = 50, time_slice = 120), mc.cores = 27)


# eric_results <- list()
# for(i in 1:27){
#   eric_results[[i]] <- quickRun(data, phy, all_model_structures[[i]], discrete_model_cid, discrete_model_cd, nSim = 10, time_slice = 120)
# }

# dats <- list(d$dat.arid.se, d$dat.temp.se, d$dat.prec.se, d$dat.pet.se)
# all_eric_results <- list()
# for(i in 1:length(dats)){
#   all_eric_results[[i]] <- mclapply(all_model_structures, function(x) quickRun(dats[[i]], phy, x, discrete_model_cid, discrete_model_cd, nSim = 50, time_slice = 120), mc.cores = 27)
# }

# for(i in 1:dim(houwie_parameters)[1]){
#   eric_results_hmm[[i]] <- mclapply(hybrid_model_structures, function(x) quickRun(data, phy, x, discrete_model_cd_hmm, discrete_model_cd, nSim = houwie_parameters[i,2], time_slice = houwie_parameters[i,1]), mc.cores = 22)
# }

# fit <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = 10, time_slice = 120, discrete_model = discrete_model_cid, continuous_model = hybrid_model_structures[[2]], recon = FALSE, mserr = "known")





