setwd("2020_houwie/")
source("hOUwieNode.R")
source("Utils.R")

# source("/space_2/jamesboyko/2020_hOUwie/hOUwieSimmap.R")
# source("/space_2/jamesboyko/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)
require(data.table)

# organize the data
# dat <- read.csv("Ericaceae_niche.csv")
dat <- read.csv("empirical/Ericaceae_niche.csv")
dat <- data.frame(sp = dat$species, reg = dat$Fruit_type, arid = dat$mean_aridity)
dat <- dat[which(apply(dat, 1, function(x) !any(is.na(x)))),]

# organize the phylgoeny
# phy <- read.tree("Ericaceae_Schwery_etal_2015.tre")
phy <- read.tree("empirical/Ericaceae_Schwery_etal_2015.tre")

# drop species not matching in the phylogeny or dataset
toDropPhy <- phy$tip.label[!phy$tip.label %in% dat$sp]
toDropDat <- dat$Sp[!dat$sp %in% phy$tip.label]
dat <- dat[!dat$sp %in% toDropDat,]
phy <- drop.tip(phy, toDropPhy)
dim(dat)[1] == length(phy$tip.label)
dat <- dat[match(phy$tip.label, dat$sp), ]

# visualize and attack
Tmax <- max(branching.times(phy))
plot(phy, show.tip.label = FALSE, x.lim = c(0, Tmax + 0.2 * Tmax))
start = Tmax + (0.005 * Tmax)
dat.prop <- (dat[,3] - min(dat[,3]))/max(dat[,3] - min(dat[,3]))
xadd = Tmax * 0.2
jitter <- 0.1 * xadd
cols = c("brown", "purple")
for(i in 1:length(dat.prop)){
  lines(list(x = c(start, start + (dat.prop[i] * xadd)), y = c(i,i)), 
        col = cols[ifelse(dat$reg[i] == "Dry", 1, 2)],
        lwd = 1)
}
legend("bottomleft", legend = c("Dry/ Abiotic","Fleshy/ Biotic"), pch=15, col = cols)
# the 2 discrete models being evaluated
discrete_model_cd <- getRateCatMat(2)
discrete_model_cid <- getFullMat(list(discrete_model_cd, discrete_model_cd), equateStateMatPars(getRateCatMat(2), c(1,2)))

# set up our continuous character models 
# variable brownian motion models set up to be character dependent, independent, and a mix of both
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
HYBBMS <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, 4)
HYBOUM <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, 4)
HYBOUMVA <- getOUParamStructure("OUMVA", "three.point", FALSE, FALSE, 4)
all_model_structures <- c(continuous_models_cd_bm, continuous_models_cid_bm[-1], continuous_models_cd_ou, continuous_models_cid_ou[-1], continuous_models_cd_bmou, continuous_models_cid_bmou, list(HYBBMS), list(HYBOUM), list(HYBOUMVA))

model_names <- c("CID_BM1", "CD_BMV", "CID+_BMV", "CID_OU1", "CD_OUA", "CD_OUV", "CD_OUVA", "CD_OUM",
                 "CD_OUMA", "CD_OUMV", "CD_OUMVA", "CID+_OUA", "CID+_OUV", "CID+_OUVA", "CID+_OUM", "CID+_OUMA",
                 "CID+_OUMV", "CID+_OUMVA", "CD_OUBM1", "CD_OUBMV", "CID+_OUBM1", "CID+_OUBMV", 
                 "HYB_BMS", "HYB_OUM", "HYB_OUMVA")
names(all_model_structures) <- model_names



singleRun <- function(phy, dat, index, all_model_structures, discrete_model_cd, discrete_model_cid){
  cont_model <- all_model_structures[[index]]
  if(dim(cont_model)[2] == 2){
    disc_model <- discrete_model_cd
    rate.cat <- 1
  }else{
    disc_model <- discrete_model_cid
    rate.cat <- 2
  }
  fit <- hOUwie(phy = phy, data = dat, rate.cat = rate.cat, nSim = 50, time_slice = max(branching.times(phy)) + 1, discrete_model = disc_model, continuous_model = cont_model, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln", n_starts = 3, ncores = 3)
  model.name <- names(all_model_structures)[[index]]
  save(fit, file = paste0("empirical_fit/FitSD=", model.name, ".Rsave"))
  return(fit)
}

model_list <- list()
for(i in 1:25){
  model_list[[i]] <- singleRun(phy, dat, 25, all_model_structures, discrete_model_cd, discrete_model_cid)
}

out <- mclapply(1:25, function(x) singleRun(phy, dat, x, all_model_structures, discrete_model_cd, discrete_model_cid), mc.cores = 25)


# load("empirical_fit/FitSD=CD_OUBM1.Rsave")
# pars <- fit$p
# 
# 
# refit <- hOUwie(phy, dat, 1, discrete_model = fit$discrete_model, continuous_model = fit$continuous_model, p = pars, nSim = 10, adaptive_sampling = FALSE)

getModelTable(list("BM1" = fitBM1, "OU1" = fitOU1, "CDBMS" = fitBMS, "CDOUM" = fitOUM, "CIDBMS" = fitCIDBMS, "CIDOUM" = fitCIDOUM, "HYBOUM" = fitHYBOUM))


res
model.fit <- res[[1]]
# devtools::install_github("bomeara/dentist")
require(dentist)

dentFuc1 <- function(par, model.fit){
  out <- hOUwie(phy = model.fit$phy, data = model.fit$data, rate.cat = model.fit$rate.cat, nSim = model.fit$nSim, index.cor = model.fit$index.cor, index.ou = model.fit$index.ou, p = par)
  return(-out$loglik)
}

getDentResults <- function(model.fit){
  best_neglnL <- -model.fit$loglik
  p.ou <- p.cor <- c()
  for(i in 1:max(model.fit$index.cor, na.rm=TRUE)){
    p.cor[i] <- na.omit(model.fit$solution.cor[model.fit$index.cor == i])[1]
  }
  for(j in 1:max(model.fit$index.ou, na.rm=TRUE)){
    p.ou[j] <- na.omit(model.fit$solution.ou[model.fit$index.ou == j])[1]
  }
  p <- c(p.cor, p.ou)
  dented_results <- dent_walk(par=p, fn=dentFuc1, best_neglnL=best_neglnL, delta=2, nsteps=1000, print_freq=250, model.fit=model.fit)
}

dent_res <- mclapply(res, function(x) getDentResults(x), mc.cores = length(res))


