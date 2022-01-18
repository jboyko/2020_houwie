# updated houwie results

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
getPropDiff <- function(map_1, map_2){
  map_1_to_fill <- matrix(0, dim(map_2$mapped.edge)[1], dim(map_2$mapped.edge)[2])
  map_1_to_fill[,as.numeric(colnames(map_1$mapped.edge))] <- map_1$mapped.edge
  branch_diff <- map_1_to_fill - map_2$mapped.edge
  out <- sum(branch_diff[branch_diff > 0])/sum(map_2$mapped.edge)
  return(out)
}

quickloadModel <- function(fit_file){
  load(fit_file)
  aicwt_cd <- getModelTable(out$cd_out)$AICwt
  aicwt_cid <- getModelTable(out$cid_out)$AICwt
  names(aicwt_cd) <- rownames(getModelTable(out$cd_out))
  names(aicwt_cid) <- rownames(getModelTable(out$cid_out))
  aic_wt_vector <- rbind(aicwt_cd, aicwt_cid)
  rownames(aic_wt_vector) <- c("CD", "CID")
  return(aic_wt_vector)
}

quickloadData <- function(fit_file){
  load(fit_file)
  cd_dat <- out$cd_out[[1]]$data
  cid_dat <- out$cid_out[[1]]$data
  both_dat <- list(cd_dat=cd_dat, cid_dat=cid_dat, phy=out$cid_out[[1]]$phy)
  return(both_dat)
}

plothOUwieData <- function(phy, data){
  phy_height <- max(branching.times(phy))
  max_x <- 1.25 * phy_height
  min_x <- 1.025 * phy_height
  scaled_data <- ((data[,3]-min(data[,3]))/(max(data[,3])-min(data[,3])) * (max_x - min_x)) + min_x
  plot.phylo(phy, show.tip.label = FALSE, x.lim = c(0, max_x))
  tiplabels(col = c("#d7191c","#2c7bb6")[data$reg], pch = c(15, 16)[data$reg])
  plot_data <- data.frame(1:length(phy$tip.label), 1:length(phy$tip.label), min_x, scaled_data, c("#d7191c","#2c7bb6")[data$reg])
  apply(plot_data, 1, function(x) lines(x = x[3:4], y = x[1:2], col = x[5], lwd = 1))
}

#### #### #### #### #### #### #### #### #### #### #### #### 
# Model specification
#### #### #### #### #### #### #### #### #### #### #### #### 
focal_model_type <- c("OUMV")
nTip <- 25
nmap <- 25
focal_dir <- paste0("simulated_fit/", focal_model_type, "/", nTip)
fit_files <- dir(focal_dir, full.names = TRUE)
fit_files <- fit_files[grep(paste0("nmap=", nmap, ".Rsave"), fit_files)]
# fit_file <- fit_files[1]

#### #### #### #### #### #### #### #### #### #### #### #### 
# Typical dataset
#### #### #### #### #### #### #### #### #### #### #### #### 
data_list <- lapply(fit_files, quickloadData)
i = 1
load(fit_files[i])
getModelTable(out$cid_out)
focal_cd_dat <- data_list[[i]]$cd_dat
focal_cid_dat <- data_list[[i]]$cid_dat
focal_phy <- data_list[[i]]$phy
par(mfrow=c(1,2))
plothOUwieData(focal_phy, focal_cd_dat)
plothOUwieData(focal_phy, focal_cid_dat)


cid_fit <- out$cid_out$cid_fit
getModelTable(out$cid_out)
plothOUwieData(focal_phy, focal_cid_dat)
sort(cid_fit$all_cont_liks)
res_discrete <- hOUwie(focal_phy, focal_cid_dat,2, discrete_model = cid_fit$discrete_model, continuous_model = cid_fit$continuous_model, time_slice = 1.1, nSim = 500, p = c(0.1, 5, 2.5, 5, 12, 24), sample_nodes = FALSE, sample_tips = FALSE)
res_tips <- hOUwie(focal_phy, focal_cid_dat,2, discrete_model = cid_fit$discrete_model, continuous_model = cid_fit$continuous_model, time_slice = 1.1, nSim = 500, p = c(0.1, 5, 2.5, 5, 12, 24), sample_nodes = FALSE, sample_tips = TRUE)
res_cherry <- hOUwie(focal_phy, focal_cid_dat,2, discrete_model = cid_fit$discrete_model, continuous_model = cid_fit$continuous_model, time_slice = 1.1, nSim = 500, p = c(0.1, 5, 2.5, 5, 12, 24), sample_nodes = TRUE, sample_tips = FALSE)
res_tips_cherry <- hOUwie(focal_phy, focal_cid_dat,2, discrete_model = cid_fit$discrete_model, continuous_model = cid_fit$continuous_model, time_slice = 1.1, nSim = 500, p = c(0.1, 5, 2.5, 5, 12, 24), sample_nodes = TRUE, sample_tips = TRUE)

xlim <- range(c(res_discrete$all_disc_liks + res_discrete$all_cont_liks, res_tips$all_disc_liks + res_tips$all_cont_liks, res_cherry$all_disc_liks + res_cherry$all_cont_liks, res_tips_cherry$all_disc_liks + res_tips_cherry$all_cont_liks))

par(mfrow=c(2,2))
hist(res_discrete$all_disc_liks + res_discrete$all_cont_liks, main = "discrete only", xlab = "joint likelihood", xlim = xlim)
hist(res_tips$all_disc_liks + res_tips$all_cont_liks, main = "tip sampled only", xlab = "joint likelihood", xlim = xlim)
hist(res_cherry$all_disc_liks + res_cherry$all_cont_liks, main = "cherry", xlab = "joint likelihood", xlim = xlim)
hist(res_tips_cherry$all_disc_liks + res_tips_cherry$all_cont_liks, main = "cherry and tip sampling", xlab = "joint likelihood", xlim = xlim)


joint_liks <- res_cherry$all_disc_liks + res_cherry$all_cont_liks
best_mappings_index <- order(joint_liks, decreasing = TRUE)

distance_from_best <- lapply(res_cherry$simmaps, function(x) try(getPropDiff(res_cherry$simmaps[[best_mappings_index[1]]], x)))
to_remove <- unlist(lapply(distance_from_best, function(x) class(x) == "try-error"))

plot(x = joint_liks[!to_remove], y = unlist(distance_from_best[!to_remove]), xlab = "joint lik", ylab = "dist")

getPropDiff(res_cherry$simmaps[[1]], res_cherry$simmaps[[1]])

joint_liks[best_mappings_index[1000]]




map_1 <- correct_map_edges(res_cherry$simmaps[best_mappings_index[1]][[1]])
map_2 <- correct_map_edges(res_cherry$simmaps[best_mappings_index[2]][[1]])
map_3 <- correct_map_edges(res_cherry$simmaps[best_mappings_index[3]][[1]])
map_4 <- correct_map_edges(res_cherry$simmaps[best_mappings_index[4]][[1]])

phy_height <- max(branching.times(map_1))
max_x <- 1.25 * phy_height
min_x <- 1.025 * phy_height

par(mfrow=c(2,2))
plotSimmap(map_1, fsize = 0.01, xlim = c(0, max_x))
plotSimmap(map_2, fsize = 0.01, xlim = c(0, max_x))
plotSimmap(map_3, fsize = 0.01, xlim = c(0, max_x))
plotSimmap(map_4, fsize = 0.01, xlim = c(0, max_x))

plot(best_map)
phy <- best_map
phy_height <- max(branching.times(phy))
max_x <- 1.25 * phy_height
min_x <- 1.025 * phy_height
scaled_data <- ((data[,3]-min(data[,3]))/(max(data[,3])-min(data[,3])) * (max_x - min_x)) + min_x
plotSimmap(phy, fsize = 0.01, xlim = c(0, max_x))
tiplabels(col = c("#d7191c","#2c7bb6")[data$reg], pch = c(15, 16)[data$reg])
plot_data <- data.frame(1:length(phy$tip.label), 1:length(phy$tip.label), min_x, scaled_data, c("#d7191c","#2c7bb6")[data$reg])
apply(plot_data, 1, function(x) lines(x = x[3:4], y = x[1:2], col = x[5], lwd = 1))




#### #### #### #### #### #### #### #### #### #### #### #### 
# Modeling description and results
#### #### #### #### #### #### #### #### #### #### #### #### 

model_table <- melt(do.call(rbind, lapply(fit_files, quickloadModel)))
model_table$X2 <- factor(model_table$X2 , levels=c("bm1_fit", "ou1_fit", "cid_fit", "cd_fit"))

ggplot(model_table, aes(x = X2, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  ylab("AICwt") + 
  xlab("Alternative model structures") +
  theme(legend.position = "none") +
  facet_wrap(~X1)

#### #### #### #### #### #### #### #### #### #### #### #### 
# Parameter evaluation
#### #### #### #### #### #### #### #### #### #### #### #### 

#### #### #### #### #### #### #### #### #### #### #### #### 
# Increasing the number of stochastic maps
#### #### #### #### #### #### #### #### #### #### #### #### 

#### #### #### #### #### #### #### #### #### #### #### #### 
# The per-map performance
#### #### #### #### #### #### #### #### #### #### #### #### 

#### #### #### #### #### #### #### #### #### #### #### #### 
# Increasing the number of tips (only up to 25 stochastic mappings)
#### #### #### #### #### #### #### #### #### #### #### #### 




