# imports and seed
set.seed(1985)
setwd("~/2020_houwie/")
source("hOUwieNode.R")

require(parallel)
require(data.table)
require(OUwie)
require(geiger)
require(ape) 
require(magick) 
require(expm) 
require(POUMM)
require(ggplot2)
require(reshape2)
require(corHMM)
require(gridExtra)

# functions
generateParameters <- function(continuous_model, alpha, sigma.sq, theta, discrete_model, rate, vector=FALSE){
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
  if(vector){
    n_p_alpha <- length(unique(na.omit(continuous_model[1,])))
    n_p_sigma <- length(unique(na.omit(continuous_model[2,])))
    n_p_theta <- length(unique(na.omit(continuous_model[3,])))
    p <- c(par_rates, alpha[seq_len(n_p_alpha)], sigma.sq[seq_len(n_p_sigma)], theta[seq_len(n_p_theta)])
  }else{
    p <- list(alpha = par_alpha, sigma.sq = par_sigma, theta = par_theta, Q = Q)
  }
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

# # parameters
# alpha <- c(3, 1.5)
# sigma.sq <- c(0.35, 1)
# theta <- c(2, 0.75)
# root.p <- c(1,0,0,0)
# theta0 <- theta[1]
# rate <- .1
# ntip <- 100
# 
# oum_model_structure <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, k = 2)[,c(1,1,2,2)]
# discrete_model <- equateStateMatPars(getFullMat(list(equateStateMatPars(getRateCatMat(2), 1:2), equateStateMatPars(getRateCatMat(2), 1:2)), equateStateMatPars(getRateCatMat(2), 1:2)), 1:3)
# pars <- generateParameters(oum_model_structure, alpha, sigma.sq, theta, discrete_model, rate)
# 
# phy <- sim.bdtree(b = .4, d = .2, stop = "taxa", n = ntip)
# phy <- drop.extinct(phy)
# phy$edge.length <- phy$edge.length/max(branching.times(phy))
# 
# oum_data <- hOUwie.sim(phy, pars$Q, c(1,0,0,0), pars$alpha, pars$sigma.sq, pars$theta[1], pars$theta)
# oum_data$data[oum_data$data[,2] == 3, 2] <- 1
# oum_data$data[oum_data$data[,2] == 4, 2] <- 2
# 
# # test_map <- reorder.simmap(oum_data$simmap, "postorder")
# 
# true_mapping_res <- hOUwie.fixed(list(oum_data$simmap), oum_data$data, 2, discrete_model, oum_model_structure, max(branching.times(oum_data$simmap))+1, sample_nodes = TRUE, p = c(.1, 3, .35, 2, .75))
# node_sampling_res <- hOUwie(oum_data$simmap, oum_data$data, 2, discrete_model, oum_model_structure, max(branching.times(oum_data$simmap))+1, 100, sample_nodes = TRUE, p = c(.1, 3, .35, 2, .75))
# no_node_sampling_res <- hOUwie(oum_data$simmap, oum_data$data, 2, discrete_model, oum_model_structure, max(branching.times(oum_data$simmap))+1, 100, sample_nodes = FALSE, p = c(.1, 3, .35, 2, .75))
# 
# true_mapping_res$loglik
# plot_data <- data.frame(node_sampling = node_sampling_res$all_cont_liks + node_sampling_res$all_disc_liks, no_node_sampling = no_node_sampling_res$all_cont_liks + no_node_sampling_res$all_disc_liks)
# plot_data <- melt(plot_data)
# plot_data$variable <- as.factor(plot_data$variable)
# 
# ggplot(plot_data, aes(x=value, fill=variable)) +
#   geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
#   # geom_density(alpha=0.7) +
#   geom_vline(aes(xintercept=true_mapping_res$loglik), linetype="dashed") +
#   theme_classic()

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
cd_model_structures <- all_model_structures[grep("CID", model_names)][-c(1,3)]
# discrete models
discrete_model_cd <- equateStateMatPars(getRateCatMat(2), 1:2)
discrete_model <- equateStateMatPars(getFullMat(list(discrete_model_cd, discrete_model_cd), getRateCatMat(2)), 1:4)




# parameters
phy <- sim.bdtree(b = .5, d = .25, stop = "taxa", n = 100) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

alpha <- c(3, 1.5)
sigma.sq <- c(0.35, 1)
theta <- c(2, 0.75)
root.p <- c(1,0,0,0)
theta0 <- theta[1]
rate <- .1
ntip <- 100

i <- 4
plot_list <- plot_data_list <- list()
for(i in 1:length(cd_model_structures)){
  focal_cont_model_structure <- cd_model_structures[[i]]
  model_name <- names(cd_model_structures)[i]
  pars <- generateParameters(focal_cont_model_structure, alpha, sigma.sq, theta, discrete_model, rate)
  houwie_data <- hOUwie.sim(phy, pars$Q, c(1,0,0,0), pars$alpha, pars$sigma.sq, pars$theta[1], pars$theta)
  houwie_data$data[houwie_data$data[,2] == 3, 2] <- 1
  houwie_data$data[houwie_data$data[,2] == 4, 2] <- 1
  p <- generateParameters(focal_cont_model_structure, alpha, sigma.sq, theta, discrete_model, rate, TRUE)
  
  true_mapping_res <- hOUwie.fixed(list(houwie_data$simmap), houwie_data$data, 2, discrete_model, focal_cont_model_structure, max(branching.times(houwie_data$simmap))+1, sample_nodes = TRUE, p = p)
  true_mapping_llik <- true_mapping_res$loglik
  node_sampling_res <- hOUwie(houwie_data$simmap, houwie_data$data, 2, discrete_model, focal_cont_model_structure, max(branching.times(houwie_data$simmap))+1, 100, sample_nodes = TRUE, adaptive_sampling = FALSE, p = p)
  no_node_sampling_res <- hOUwie(houwie_data$simmap, houwie_data$data, 2, discrete_model, focal_cont_model_structure, max(branching.times(houwie_data$simmap))+1, 100, sample_nodes = FALSE, adaptive_sampling = FALSE, p = p)
  node_adpt_sampling_res <- hOUwie(houwie_data$simmap, houwie_data$data, 2, discrete_model, focal_cont_model_structure, max(branching.times(houwie_data$simmap))+1, 100, sample_nodes = TRUE, adaptive_sampling = TRUE, p = p)
  # no_node_adpt_sampling_res <- hOUwie(houwie_data$simmap, houwie_data$data, 2, discrete_model, focal_cont_model_structure, max(branching.times(houwie_data$simmap))+1, 100, sample_nodes = FALSE, adaptive_sampling = TRUE, p = p)
  plot_data <- data.frame(
    cherry_conditional = node_sampling_res$all_cont_liks + node_sampling_res$all_disc_liks, 
    discrete_only = no_node_sampling_res$all_cont_liks + no_node_sampling_res$all_disc_liks,
    adaptive_sampling = node_adpt_sampling_res$all_cont_liks + node_adpt_sampling_res$all_disc_liks
    )
  
  plot_data_list[[i]] <- rbind(
    data.frame(type = "cherry_conditional",
               disc = node_sampling_res$all_disc_liks,
               cont = node_sampling_res$all_cont_liks,
               joint = node_sampling_res$all_cont_liks + node_sampling_res$all_disc_liks),
    data.frame(type = "discrete_only",
               disc = no_node_sampling_res$all_disc_liks,
               cont = no_node_sampling_res$all_cont_liks,
               joint = no_node_sampling_res$all_cont_liks + no_node_sampling_res$all_disc_liks),
    data.frame(type = "adaptive_sampling",
               disc = node_adpt_sampling_res$all_disc_liks,
               cont = node_adpt_sampling_res$all_cont_liks,
               joint = node_adpt_sampling_res$all_cont_liks + node_adpt_sampling_res$all_disc_liks))
    
  plot_data <- melt(plot_data)
  plot_data$variable <- as.factor(plot_data$variable)
  
  plot_list[[i]] <- ggplot(plot_data, aes(x=value, fill=variable)) +
    geom_histogram(alpha=0.66, color="black", position="identity") +
    scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) + 
    geom_vline(xintercept = true_mapping_llik, linetype="dashed") +
    theme_classic() +
    ggtitle(paste0(letters[i], ") ", model_name))
}

final_plot <- grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], nrow = 5)

ggsave(final_plot, filename = "figures/raw/compare_simmap_generation.pdf", height = 10, width = 10, units = "in")


plot_list_b <- list()
for(i in 1:length(plot_data_list)){
  model_name <- names(cd_model_structures)[i]
  plot_list_b[[i]] <- ggplot(plot_data_list[[i]], aes(x=disc, y = cont, fill = type)) +
    geom_point(size = 3, shape=21, alpha = 0.75, color = "black") +
    scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) + 
    ggtitle(paste0(letters[i], ") ", model_name)) +
    theme_classic()
}

final_plot_b <- grid.arrange(plot_list_b[[1]], plot_list_b[[2]], plot_list_b[[3]], plot_list_b[[4]], plot_list_b[[5]], plot_list_b[[6]], plot_list_b[[7]], plot_list_b[[8]], plot_list_b[[9]], plot_list_b[[10]], nrow = 5)

ggsave(final_plot_b, filename = "figures/raw/compare_simmap_generation_points.pdf", height = 10, width = 10, units = "in")



