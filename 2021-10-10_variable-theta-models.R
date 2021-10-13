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

model_param_files <- dir("nuiss_par_test/", full.names = TRUE) # simulate under model i, fit under all models
param_model_res <- lapply(all_model_names, function(x) load_model(x, model_param_files))
names(param_model_res) <- paste0("M", model_numbers) # each element is a specific model

#### #### #### #### #### #### #### #### #### #### #### #### 
# the checking of variable theta models 
#### #### #### #### #### #### #### #### #### #### #### #### 
model_fits <- lapply(param_model_res$M8, function(x) getBestModel(x$model_res))
simulated_data <- lapply(param_model_res$M8, function(x) x$simulated_data$data[,2:3])
root_states <- lapply(param_model_res$M8, function(x) getRootState(x$simulated_data$simmap[[1]]))
expected_data <- lapply(model_fits, function(x) cbind(reg=x$data[,2], x=x$expected_vals))


cols <- c("#e6550d", "#756bb1")

# subset data to be state 1 or 2
root_state <- 2
sim_dat <- simulated_data[unlist(root_states) == root_state]
model_fit <- model_fits[unlist(root_states) == root_state]
exp_dat <- expected_data[unlist(root_states) == root_state]
cart_max <-  30

## state 1 is the root
# plot the main histogram
df <- as.data.frame(do.call(rbind, sim_dat))
df[,1] <- as.factor(df[,1])
hist <- ggplot(df, aes(x=x, fill = reg)) +
  geom_histogram(color="black") +
  scale_fill_manual(values = cols) + 
  coord_cartesian(xlim=c(0, cart_max)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y=element_blank())


# model fits with root state 1
mean_est <- melt(do.call(rbind, lapply(model_fit, function(x) x$solution.cont[3,])))
# plot boxplot
bp <- ggplot(mean_est, aes(y = value, fill = Var2)) +
  scale_fill_manual(values = cols) + 
  geom_hline(yintercept=12, linetype="dashed", color = cols[1]) +
  geom_hline(yintercept=24, linetype="dashed", color = cols[2]) + 
  geom_boxplot() + 
  coord_flip(ylim=c(0, cart_max)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# expected hist
df <- as.data.frame(do.call(rbind, exp_dat))
df[,1] <- as.factor(df[,1])
hist_exp <- ggplot(df, aes(x=x, fill = reg)) +
  geom_histogram(color="black") +
  scale_fill_manual(values = cols) + 
  coord_cartesian(xlim=c(0, cart_max)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y=element_blank())
  

grid.arrange(bp, hist, hist_exp, heights=c(1,3,3))



