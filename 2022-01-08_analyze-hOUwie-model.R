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
require(ggplot2)
require(gridExtra)
require(reshape)

#### #### #### #### #### #### #### #### #### #### #### #### 
# Functions
#### #### #### #### #### #### #### #### #### #### #### #### 

quickload <- function(fit_file){
  load(fit_file)
  aicwt_cd <- getModelTable(out$cd_out)$AICwt
  aicwt_cid <- getModelTable(out$cid_out)$AICwt
  names(aicwt_cd) <- rownames(getModelTable(out$cd_out))
  names(aicwt_cid) <- rownames(getModelTable(out$cid_out))
  aic_wt_vector <- rbind(aicwt_cd, aicwt_cid)
  rownames(aic_wt_vector) <- c("CD", "CID")
  return(aic_wt_vector)
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
bm1_model <- all_model_structures[[1]]
ou1_model <- all_model_structures[[4]]

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
ntips <- c(250)
nmap <- 25

# the dataests come from a particular structure and modeling results will also follow the same saving protocol
model_type_list <- list()
for(i in 1:length(model_types)){
  focal_model_type <- model_types[i]
  focal_models <- sort(model_names[grep(paste0("_", focal_model_type, "$"), model_names)])
  continuous_model_cd <- all_model_structures[names(all_model_structures) == focal_models[1]][[1]]
  continuous_model_cid <- all_model_structures[names(all_model_structures) == focal_models[2]][[1]]
  for(j in 1:length(ntips)){
    nTip <- ntips[j]
    focal_dir <- paste0("simulated_fit/", focal_model_type, "/", nTip)
    fit_files <- dir(focal_dir, full.names = TRUE)
    fit_files <- fit_files[grep(paste0("nmap=", nmap, ".Rsave"), fit_files)]
    model_type_list[[i]] <- melt(do.call(rbind, lapply(fit_files, quickload)))
    model_type_list[[i]]$X2 <- factor(model_type_list[[i]]$X2 , levels=c("bm1_fit", "ou1_fit", "cid_fit", "cd_fit"))
  }
}
names(model_type_list) <- model_types

big_df <- do.call(rbind, model_type_list)
big_df$gen_model <- gsub("\\..*", "", rownames(big_df))

ggplot(big_df, aes(x = X2, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  ylab("AICwt") + 
  xlab("Alternative model structures") +
  theme(legend.position = "none") +
  facet_wrap(~gen_model + X1, ncol = 4)


# fit_file <- fit_files[11]
# wt_table <- do.call(rbind, lapply(fit_files, quickload))
# colMeans(wt_table)

# a closer examination of only OUM 25 tips
# for each model type i want to generate a dataset consistent with the CD and one consistent with CID+
focal_model_type <- "OUM"
# for each number of tips
nTip <- c(25)
nmap <- 25

# the dataests come from a particular structure and modeling results will also follow the same saving protocol
getPlotQuick <- function(nTip, nmap){
  focal_models <- sort(model_names[grep(paste0("_", focal_model_type, "$"), model_names)])
  continuous_model_cd <- all_model_structures[names(all_model_structures) == focal_models[1]][[1]]
  continuous_model_cid <- all_model_structures[names(all_model_structures) == focal_models[2]][[1]]
  focal_dir <- paste0("simulated_fit/", focal_model_type, "/", nTip)
  fit_files <- dir(focal_dir, full.names = TRUE)
  fit_files <- fit_files[grep(paste0("nmap=", nmap, ".Rsave"), fit_files)]
  model_type_table <- melt(do.call(rbind, lapply(fit_files, quickload)))
  model_type_table$X2 <- factor(model_type_table$X2 , levels=c("bm1_fit", "ou1_fit", "cid_fit", "cd_fit"))
  my_plot <- ggplot(model_type_table, aes(x = X2, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    ylab("AICwt") + 
    xlab("Alternative model structures") +
    theme(legend.position = "none") +
    facet_wrap(~X1)
  return(my_plot)
}

p_25 <- getPlotQuick(25, 25)
p_100 <- getPlotQuick(25, 100)
p_250 <- getPlotQuick(25, 250)

grid.arrange(p_25, p_100, p_250)

# big_df <- do.call(rbind, model_type_list)
# big_df$gen_model <- gsub("\\..*", "", rownames(big_df))


nTip <- c(25)
nmap <- 25
focal_models <- sort(model_names[grep(paste0("_", focal_model_type, "$"), model_names)])
continuous_model_cd <- all_model_structures[names(all_model_structures) == focal_models[1]][[1]]
continuous_model_cid <- all_model_structures[names(all_model_structures) == focal_models[2]][[1]]
focal_dir <- paste0("simulated_fit/", focal_model_type, "/", nTip)
fit_files <- dir(focal_dir, full.names = TRUE)
fit_files <- fit_files[grep(paste0("nmap=", nmap, ".Rsave"), fit_files)]
data_files <- dir(paste0("simulated_data/", focal_model_type, "/", nTip), full.names = TRUE)




load(fit_files[1]); getModelTable(out$cid_out)
load(data_files[1]); plot(focal_dat$dat_cid$simmap)
table_cd <- getModelTable(out$cd_out)
table_cid <- getModelTable(out$cid_out)


cid_dat <- focal_dat$dat_cid$data
cid_dat[cid_dat$reg == 2, 2] <- 1
cid_dat[(cid_dat$reg == 3 | cid_dat$reg == 4), 2] <- 2
phy <- reorder(focal_dat$dat_cid$simmap, "pruningwise")


p <- out$cid_out$cid_fit$p
p <- c(.1, 5, 2.5, 12, 24)
res <- hOUwie(phy, cid_dat, 2, discrete_model = discrete_model_cid, continuous_model = all_model_structures$`CID+_OUM`, time_slice = 1.1, nSim = 100, p = p); res$loglik
res$simmaps <- lapply(res$simmaps, correct_map_edges)
class(res$simmaps) <- c("multiSimmap", "multiPhylo")

prop_diffs <- unlist(lapply(res$simmaps, function(x) getPropDiff(x, phy))) 
min(prop_diffs); mean(prop_diffs)
# par(mfrow=c(1,2))
plot(x = prop_diffs, y = res$all_cont_liks + res$all_disc_liks, xlim=c(0,1), xlab="Difference from true map", ylab="lnLik")

getPropDiff <- function(map_1, map_2){
  map_1_to_fill <- matrix(0, dim(map_2$mapped.edge)[1], dim(map_2$mapped.edge)[2])
  map_1_to_fill[,as.numeric(colnames(map_1$mapped.edge))] <- map_1$mapped.edge
  branch_diff <- map_1_to_fill - map_2$mapped.edge
  out <- sum(branch_diff[branch_diff > 0])/sum(map_2$mapped.edge)
  return(out)
}



res <- hOUwie(phy, cid_dat, 2, discrete_model = discrete_model_cid, continuous_model = all_model_structures$`CID+_OUM`, time_slice = 1.1, nSim = 1000)


cid_dat[cid_dat$reg == 3, 2] <- 1
cid_dat[cid_dat$reg == 4, 2] <- 1
cid_dat[(cid_dat$reg == 3 | cid_dat$reg == 4), 2] <- 2



for(i in 1:11){
  print(i)
  load(fit_files[i]); getModelTable(out$cid_out)
  load(data_files[i])
  
  table_cid <- getModelTable(out$cid_out)
  cid_dat <- focal_dat$dat_cid$data
  cid_dat[cid_dat$reg == 2, 2] <- 1
  cid_dat[(cid_dat$reg == 3 | cid_dat$reg == 4), 2] <- 2
  phy <- reorder(focal_dat$dat_cid$simmap, "pruningwise")
  # p <- out$cid_out$cid_fit$p
  p_1 <- c(.1, 5, 2.5, 12, 24)
  p_2 <- c(.1, 1e-5, 2.5, 12, 24)
  p_3 <- out$cid_out$cid_fit$p
  res_1 <- hOUwie(phy, cid_dat, 2, discrete_model = discrete_model_cid, continuous_model = all_model_structures$`CID+_OUM`, time_slice = 1.1, nSim = 100, p = p_1); res_1$loglik
  res_1$simmaps <- lapply(res_1$simmaps, correct_map_edges)
  class(res_1$simmaps) <- c("multiSimmap", "multiPhylo")
  res_2 <- hOUwie(phy, cid_dat, 2, discrete_model = discrete_model_cid, continuous_model = all_model_structures$`CID+_OUM`, time_slice = 1.1, nSim = 100, p = p_2); res_2$loglik
  res_2$simmaps <- lapply(res_2$simmaps, correct_map_edges)
  class(res_2$simmaps) <- c("multiSimmap", "multiPhylo")
  res_3 <- hOUwie(phy, cid_dat, 2, discrete_model = discrete_model_cid, continuous_model = all_model_structures$`CID+_OUM`, time_slice = 1.1, nSim = 100, p = p_3); res_3$loglik
  res_3$simmaps <- lapply(res_3$simmaps, correct_map_edges)
  class(res_3$simmaps) <- c("multiSimmap", "multiPhylo")
  
  prop_diffs_1 <- unlist(lapply(res_1$simmaps, function(x) getPropDiff(x, phy))) 
  prop_diffs_2 <- unlist(lapply(res_2$simmaps, function(x) getPropDiff(x, phy))) 
  prop_diffs_3 <- unlist(lapply(res_3$simmaps, function(x) getPropDiff(x, phy))) 
  min(prop_diffs_1); mean(prop_diffs_1)
  min(prop_diffs_2); mean(prop_diffs_2)
  min(prop_diffs_3); mean(prop_diffs_3)
  
  pdf(file = paste0("~/Desktop/tmp/prop_image_", i, ".pdf"), width = 10, height = 5)
  par(mfrow=c(1,3))
  ylim <- c(range(c(res_1$all_cont_liks + res_1$all_disc_liks, res_2$all_cont_liks + res_2$all_disc_liks, res_3$all_cont_liks + res_3$all_disc_liks)))
  plot(x = prop_diffs_1, y = res_1$all_cont_liks + res_1$all_disc_liks, xlim=c(0,1), xlab="Difference from true map", ylab="lnLik", ylim=ylim, main = "pars = gen")
  plot(x = prop_diffs_2, y = res_2$all_cont_liks + res_2$all_disc_liks, xlim=c(0,1), xlab="Difference from true map", ylab="lnLik", ylim=ylim, main = "pars = gen but alpha=1e-5")
  plot(x = prop_diffs_3, y = res_3$all_cont_liks + res_3$all_disc_liks, xlim=c(0,1), xlab="Difference from true map", ylab="lnLik", ylim=ylim, main = paste0("pars = mle & dAIC = ", table_cid$dAIC[4]))
  dev.off()
}
