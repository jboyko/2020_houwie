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
# the checking of CID, CD, CID+ (1,2,3)
#### #### #### #### #### #### #### #### #### #### #### #### 
# define the non hmm models
# model_class <- c("CID", "CD", "CID+")[c(1,2,1,1, rep(2, 7), rep(3, 7), 2, 2, 3, 3)]

# clean failed optimizations out of the equation
getAdjustedModTable <- function(res){
  res_i <- res[[2]]
  mod_table <- getModelTable(res_i)
  mod_names <- rownames(mod_table)
  if(diff(range(mod_table$AIC)) > 1e10){
    max_aic <- max(mod_table$AIC)
    res_i <- res_i[abs(mod_table$AIC - max_aic)  < 1e10]
    mod_names <- mod_names[abs(mod_table$AIC - max_aic)  < 1e10]
    mod_table <- getModelTable(res_i)
    rownames(mod_table) <- mod_names
  }
  return(mod_table)
}
# res <- param_model_res$M1[[2]]
# getAdjustedModTable(res)

# get the dAIC for a model run named by the model class
getdAIC_table <- function(obj){
  model_index <- as.numeric(rownames(obj))
  dAIC_table <- data.frame(model_class=model_class[model_index], model_number = model_index, dAIC=obj$dAIC)
  return(dAIC_table)
}
# mod_table_list <- lapply(param_model_res$M1, getAdjustedModTable)
# obj <- mod_table_list[[1]]
# getdAIC_table(obj)

# create a dataframe for plotting
model_index <- c(1,2,1,1, rep(2, 7), rep(3, 7), 2, 2, 3, 3)
model_class <- c("CID", "CD", "CID+")[model_index]

# for each model
aic_list <- list()
for(i in 1:length(param_model_res)){
  mod_table_list <- lapply(param_model_res[[i]], getAdjustedModTable)
  AIC_table_list <- lapply(mod_table_list, getdAIC_table)
  aic_table <- do.call(rbind, AIC_table_list)
  aic_table$simulating_model <- model_class[i]
  aic_table$simulating_model_number <- i
  aic_list[[i]] <- aic_table
}

full_aic_table <- do.call(rbind, aic_list)

ggplot(full_aic_table, aes(x = model_class, y = dAIC)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(0, 25)) + 
  facet_wrap(~simulating_model)

# finding the proportion best for each scenario
# CD simulating 
cd_table <- full_aic_table[full_aic_table$simulating_model == "CD",]
summary(as.factor(cd_table$model_class[cd_table$dAIC == 0]))
# CID simulating
cid_table <- full_aic_table[full_aic_table$simulating_model == "CID",]
summary(as.factor(cid_table$model_class[cid_table$dAIC == 0]))
# CID+ simulating
cid2_table <- full_aic_table[full_aic_table$simulating_model == "CID+",]
summary(as.factor(cid2_table$model_class[cid2_table$dAIC == 0]))

# single example model
# # define the non hmm models
# model_index <- c(1,2,1,1, rep(2, 7), rep(3, 7), 2, 2, 3, 3)
# model_class <- c("CID", "CD", "CID+")[model_index]
# # clean failed optimizations out of the equation
# mod_table_list <- lapply(param_model_res$M14, getAdjustedModTable)
# 
# #  get the dAIC for a model run named by the model class
# AIC_table_list <- lapply(mod_table_list, getdAIC_table)
# aic_table <- do.call(rbind, AIC_table_list)
# 
# # create a plot
# ggplot(aic_table, aes(x = model_class, y = dAIC)) +
#   geom_boxplot() + 
#   coord_cartesian(ylim=c(0, 150))




