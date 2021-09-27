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
# functions
#### #### #### #### #### #### #### #### #### #### #### #### 

  getData <- function(csv){
    dat <- read.csv(csv)
    dat.temp.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, temp = dat$mean_temp, se_temp = dat$within_sp_var_temp)
    dat.temp.se <- dat.temp.se[which(apply(dat.temp.se, 1, function(x) !any(is.na(x)))),]
    dat.temp <- data.frame(sp = dat$species, reg = dat$Fruit_type, temp = dat$mean_temp)
    dat.temp <- dat.temp[which(apply(dat.temp, 1, function(x) !any(is.na(x)))),]
    
    dat.prec.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, prec = dat$mean_prec, se_prec = dat$within_sp_var_prec)
    dat.prec.se <- dat.prec.se[which(apply(dat.prec.se, 1, function(x) !any(is.na(x)))),]
    dat.prec <- data.frame(sp = dat$species, reg = dat$Fruit_type, prec = dat$mean_prec)
    dat.prec <- dat.prec[which(apply(dat.prec, 1, function(x) !any(is.na(x)))),]
    
    dat.pet.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, pet = dat$mean_pet, se_pet = dat$within_sp_var_pet)
    dat.pet.se <- dat.pet.se[which(apply(dat.pet.se, 1, function(x) !any(is.na(x)))),]
    dat.pet <- data.frame(sp = dat$species, reg = dat$Fruit_type, pet = dat$mean_pet)
    dat.pet <- dat.pet[which(apply(dat.pet, 1, function(x) !any(is.na(x)))),]
    
    dat.arid.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, arid = dat$mean_aridity, se_temp = dat$within_sp_var_aridity)
    dat.arid.se <- dat.arid.se[which(apply(dat.arid.se, 1, function(x) !any(is.na(x)))),]
    dat.arid <- data.frame(sp = dat$species, reg = dat$Fruit_type, arid = dat$mean_aridity)
    dat.arid <- dat.arid[which(apply(dat.arid, 1, function(x) !any(is.na(x)))),]
    return(list(dat.temp = dat.temp,
                dat.temp.se = dat.temp.se,
                dat.prec = dat.prec,
                dat.prec.se = dat.prec.se,
                dat.pet = dat.pet,
                dat.pet.se = dat.pet.se,
                dat.arid = dat.arid,
                dat.arid.se = dat.arid.se
    ))
  }

quickRun <- function(data, phy, continuous_model, discrete_model_cid, discrete_model_cd){
  if(dim(continuous_model)[2] == 2){
    discrete_model <- discrete_model_cd
    rate.cat <- 1
  }else{
    discrete_model <- discrete_model_cid
    rate.cat <- 2
  }
  fit <- hOUwie(phy = phy, data = data, rate.cat = rate.cat, nSim = 50, time_slice = 20, discrete_model = discrete_model, continuous_model = continuous_model, recon = FALSE)
  return(fit)
}

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

data <- d$dat.arid
phy <- keep.tip(t, data$sp)

test2 <- mclapply(all_model_structures, function(x) quickRun(data, phy, x, discrete_model_cid, discrete_model_cd), mc.cores = 22)



