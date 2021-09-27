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
# models to be run
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

#### #### #### #### #### #### #### #### #### #### #### #### 
# prerequisites
#### #### #### #### #### #### #### #### #### #### #### #### 

phy <- read.tree("empirical/Brassicales_5g.tre")
phy <- keep.tip(phy, phy$tip.label[grep("Cleomaceae", phy$tip.label)])
phy$tip.label <- gsub("Cleomaceae_", "", phy$tip.label)
dat <- read.csv("empirical/cleome_data.csv")
tmp <- dat$c4
tmp[is.na(tmp)] <- 2
plot(phy, cex = 0.6, show.tip.label = FALSE)
tiplabels(pch = 16, col = viridis(3)[tmp+1])
legend("topleft", legend = c("C3", "C4", "Unknown"), col = viridis(3), pch=16)

#### #### #### #### #### #### #### #### #### #### #### #### 
# corhmm models
#### #### #### #### #### #### #### #### #### #### #### #### 

dat_corhmm <- dat[,1:2]
dat_corhmm$c4[is.na(dat_corhmm$c4)] <- "?"

er_mat <- equateStateMatPars(getRateCatMat(2), c(1,2))
ard_mat <- getRateCatMat(2)
er_2_mat <- getFullMat(list(er_mat,er_mat), er_mat)
ard_2_mat <- getFullMat(list(ard_mat,ard_mat), er_mat)
er_ard_2_mat <- getFullMat(list(er_mat,ard_mat), er_mat)

ER_1 <- corHMM(phy, dat_corhmm, 1, rate.mat = er_mat)
ARD_1 <- corHMM(phy, dat_corhmm, 1, rate.mat = ard_mat)
ER_2 <- corHMM(phy, dat_corhmm, 2, rate.mat = er_2_mat)
ARD_2 <- corHMM(phy, dat_corhmm, 2, rate.mat = ard_2_mat)
ER_ARD_2 <- corHMM(phy, dat_corhmm, 2, rate.mat = er_ard_2_mat)

#### #### #### #### #### #### #### #### #### #### #### #### 
# houwie models
#### #### #### #### #### #### #### #### #### #### #### #### 







