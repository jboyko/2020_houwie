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
require(reshape2)
require(RColorBrewer)
require(gridExtra)
require(gtable)

#### #### #### #### #### #### #### #### #### #### #### #### 
# run
#### #### #### #### #### #### #### #### #### #### #### #### 

# the model names
model_numbers <- 1:22
all_model_names <- paste0("M", model_numbers, "_")

# the model strucutres being tested:
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

model_param_files <- dir("model_power_test/", full.names = TRUE) # simulate under model i, fit under all models
param_model_res <- lapply(all_model_names, function(x) load_model(x, model_param_files))
names(param_model_res) <- paste0("M", model_numbers) # each element is a specific model

# param_model_res is a list of 22 elements (one for each model)
# within each element of 22, there is a list of 40 (one for each simulated dataset)
# within each element of 40, there is a list of 2 (one complete simulated dataset, and another list of model fits)
# within the second element of 2, there is a list of 22 model fits to the simulated data (from element 1)

#### #### #### #### #### #### #### #### #### #### #### #### 
# the checking of variable theta models 
#### #### #### #### #### #### #### #### #### #### #### #### 
# define the non hmm models
# non_hmm_models <- c(1,2,4,5,6,7,8,9,10,11,19,20)
# 
# # these are the variable theta models
# M8 <- param_model_res$M8
# M9 <- param_model_res$M9
# M10 <- param_model_res$M10
# M11 <- param_model_res$M11
# 
# getModelTable(M8[[40]][[2]])
# 
# mod_avg_m8 <- lapply(M8, function(x) getModelAvgParams(x[[2]][non_hmm_models])) 
# 
# model_i <- 8
# 
# model_plots <- list()
# for(model_i in 8:11){
#   model_fits <- lapply(param_model_res[[model_i]], function(x) x[[2]][non_hmm_models])
#   simulated_data <- lapply(param_model_res[[model_i]], function(x) x$simulated_data$data[,2:3])
#   root_states <- lapply(param_model_res[[model_i]], function(x) getRootState(x$simulated_data$simmap[[1]]))
#   # expected_data <- lapply(model_fits, function(x) cbind(reg=x$data[,2], x=x$expected_vals))
#   mod_avg_pars <- lapply(param_model_res[[model_i]], function(x) getModelAvgParams(x[[2]][non_hmm_models])) 
#   expected_data <- lapply(mod_avg_pars, function(x) x$mod_avg_expc)
#   lapply(mod_avg_pars, function(x) x$mod_avg_cont[3,])
#   
#   cols <- c("#377eb8", "#e41a1c")
#   df_list <- plots <- list()
#   # subset data to be state 1 or 2
#   for(root_state in 1:2){
#     sim_dat <- simulated_data[unlist(root_states) == root_state]
#     model_fit <- model_fits[unlist(root_states) == root_state]
#     exp_dat <- expected_data[unlist(root_states) == root_state]
#     mod_avg_par <- mod_avg_pars[unlist(root_states) == root_state]
#     cart_max <-  30
#     title_name <- paste0("Root state ", root_state)
#     ## state 1 is the root
#     # plot the true data and expected data histogram
#     df <- as.data.frame(do.call(rbind, sim_dat))
#     df$xx <- unlist(exp_dat)
#     df[,1] <- as.factor(df[,1])
#     hist <- ggplot(df, aes(x=x, fill = reg)) +
#       geom_histogram(color="black", bins = 20) +
#       scale_fill_manual(values = cols) + 
#       coord_cartesian(xlim=c(0, cart_max)) +
#       theme_classic() +
#       theme(legend.position = "none", 
#             axis.title.y=element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),
#             axis.title.x=element_blank(), axis.text.x = element_text(color="white")) +
#       ggtitle(ifelse(root_state == 1, "Observed data", ""))
#     
#     hist_exp <- ggplot(df, aes(x=xx, fill = reg)) +
#       geom_histogram(color="black", bins = 20) +
#       scale_fill_manual(values = cols) + 
#       coord_cartesian(xlim=c(0, cart_max)) +
#       theme_classic() +
#       theme(legend.position = "none", 
#             axis.title.y=element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
#       xlab("Continuous trait value") +
#       ggtitle(ifelse(root_state == 1, "Model expectation", ""))
#     
#     # plotting the distribution of theta estimates
#     theta_table <- do.call(rbind, lapply(mod_avg_par, function(x) x$mod_avg_cont[3,]))
#     theta_estimates <- melt(theta_table)
#     # plot boxplot
#     bp <- ggplot(theta_estimates, aes(y = value, fill = Var2)) +
#       scale_fill_manual(values = cols) + 
#       geom_hline(yintercept=12, linetype="dashed", color = cols[1]) +
#       geom_hline(yintercept=24, linetype="dashed", color = cols[2]) + 
#       geom_boxplot() + 
#       coord_flip(ylim=c(0, cart_max)) +
#       theme_classic() +
#       theme(legend.position = "none", plot.title = element_text(hjust = 0.5, colour = cols[root_state], face="bold"),
#             axis.title.y=element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(), 
#             axis.text.x = element_text(color="white")) +
#       ylab("Model averaged theta estimate") +
#       ggtitle(title_name)
#     
#     plots[[root_state]] <- grid.arrange(bp, hist, hist_exp, heights=c(2,3,3))
#   }
#   
#   model_plots[[model_i-7]] <- grid.arrange(plots[[1]], plots[[2]], nrow=1)
# }
# 
# # four separate models
# grid.arrange(model_plots[[1]], model_plots[[2]], model_plots[[3]], model_plots[[4]])

# summarized over all four models
# model_i <- 1
# res <- param_model_res[[model_i]]

getModDf <- function(res){
  simulated_data <- lapply(res, function(x) x$simulated_data$data[,2:3])
  root_states <- lapply(res, function(x) getRootState(x$simulated_data$simmap[[1]]))
  mod_avg_pars <- lapply(res, function(x) getModelAvgParams(x[[2]][non_hmm_models], force=FALSE)) 
  expected_data <- lapply(mod_avg_pars, function(x) x$mod_avg_expc)
  mod_avg_theta <- lapply(mod_avg_pars, function(x) x$mod_avg_cont[3,])
  Df <- list()
  for(i in 1:length(simulated_data)){
    Df[[i]] <- cbind(simulated_data[[i]], expected_data[[i]], root_states[[i]], mod_avg_theta[[i]][1], mod_avg_theta[[i]][2])
  }
  out <- do.call(rbind, Df)
  colnames(out) <- c("reg", "x", "xx", "root_state", "theta_1", "theta_2")
  rownames(out) <- NULL
  return(out)
}

non_hmm_models <- c(1,2,4,5,6,7,8,9,10,11,19,20)
df_1 <- do.call(rbind, lapply(param_model_res[8:11], getModDf))
df_1$reg <- as.factor(df_1$reg)
cart_max <-  30
cols <- c("#377eb8", "#e41a1c")

hist_obs <- ggplot(df_1, aes(x=x, fill = reg)) +
  scale_fill_manual(values = cols) + 
  geom_histogram(color="black", bins = 20) +
  coord_cartesian(xlim=c(0, cart_max)) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.y=element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.text.x = element_text(color="white")) +
  ggtitle("b) Observed data") +
  facet_wrap(~root_state) +
  theme(strip.background = element_blank(),strip.text.x = element_blank())

hist_exp <- ggplot(df_1, aes(x=xx, fill = reg)) +
  scale_fill_manual(values = cols) + 
  geom_histogram(color="black", bins = 30) +
  coord_cartesian(xlim=c(0, cart_max)) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.y=element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  xlab("Continuous trait value") +
  ggtitle("c) Model expectation") +
  facet_wrap(~root_state) +
  theme(strip.background = element_blank(),strip.text.x = element_blank())


# plotting the distribution of theta estimates
theta_estimates <- melt(df_1[,5:6])
theta_estimates$root_state <- as.factor(rep(df_1[,4], 2))
# theta_estimates <- melt(df_1[,5:6][df_1[,4]==1,])
# plot boxplot
bp <- ggplot(theta_estimates, aes(y = value, fill = variable)) +
  scale_fill_manual(values = cols) + 
  geom_hline(yintercept=12, linetype="dashed", color = cols[1]) +
  geom_hline(yintercept=24, linetype="dashed", color = cols[2]) + 
  geom_boxplot() + 
  coord_flip(ylim=c(0, cart_max)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y=element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(), 
        axis.text.x = element_text(color="white"), axis.title.x=element_blank()) +
  ggtitle("a) Model averaged theta estimate") +
  facet_wrap(~root_state) +
  theme(strip.background = element_blank(),strip.text.x = element_blank())

grid.arrange(bp, hist_obs, hist_exp, heights=c(2,3,3))







