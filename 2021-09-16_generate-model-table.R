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
require(viridis)

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

#### #### #### #### #### #### #### #### #### #### #### #### 
# run
#### #### #### #### #### #### #### #### #### #### #### #### 

# the model table will have a row for every model (of which there are 22), and 12 columns (model, type, k, par_a, par_s, par_t, err_a, err_s, err_t, power_1, power_2, power_3)
model_table <- as.data.frame(matrix(data = NA, nrow=22, ncol=15, dimnames = list(paste0("M", 1:22), c("type_1","type_2", "k","alpha_free","sigma_free","theta_free","alpha_sign_error","sigma_sign_error","theta_sign_error","alpha_rmse","sigma_rmse","theta_rmse", "prop_best","avg_AICwt","prop_CD"))))

# fill out the first column of the table, the type_1 is either BM, OU, or mixed
model_table[,1] <- c("OU", "mixed", "BM")[(unlist(lapply(all_model_structures, function(x) sum(is.na(x[1,]))/length(x[1,])))*2)+1]

# fill out the second column of the table, the type_2 is either CD or CID (non-HMM or HMM) with 2 exceptions
model_table[,2] <- c("CD", "CID")[as.numeric(unlist(lapply(all_model_structures, function(x) dim(x)[2] == 4)))+1]
cd_vecor <- model_table[,2]
# fill out the third column of the table, the number of parameters
k_par_disc <- c(1,3)[as.numeric(unlist(lapply(all_model_structures, function(x) dim(x)[2] == 4)))+1] 
k_par_cont <- unlist(lapply(all_model_structures, function(x) max(x, na.rm = TRUE)))
model_table[,3] <- k_par_disc + k_par_cont

# fill out the fourth column of the table, whether or not alpha is free
model_table[,4] <- c(FALSE, TRUE)[unlist(lapply(all_model_structures, function(x) length(unique(x[1,]))))]

# fill out the fifth column of the table, whether or not sigma is free
model_table[,5] <- c(FALSE, TRUE)[unlist(lapply(all_model_structures, function(x) length(unique(x[2,]))))]

# fill out the sixth column of the table, whether or not theta is free
model_table[,6] <- c(FALSE, TRUE)[unlist(lapply(all_model_structures, function(x) length(unique(x[3,]))))]

# adjust the number of CID models based on whether all the parameters are not fixed to be equal
model_table[which(apply(model_table[,c(4:6)], 1, function(x) all(x == FALSE))),2] <- "CID"

# fill out the 7:9 columns. measures of sign error. how often do we guess the direction of the parameters when they are free
model_param_files <- dir("nuiss_par_test/", full.names = TRUE) # simulate under model i, fit under all models

# load the results of all model runs
param_model_res <- lapply(all_model_names, function(x) load_model(x, model_param_files))
names(param_model_res) <- paste0("M", model_numbers) # each element is a specific model
# example of the strucure of param_model_res
# model 1: param_model_res[[1]] 
# model 1, iteration 1: param_model_res[[1]][[1]]
# model 1, iteration 1, simulating information: param_model_res[[1]][[1]]$simulated_data
# model 1, iteration 1, 9 fitted models : param_model_res[[1]][[1]]$model_res
# for each iteration we will only take the maximum likelihood (same model with different nuiss parameters, so we are treating this as random restarts)
sign_error_table <- do.call(rbind, lapply(param_model_res, function(x) getSignErrorList(x)))
# sign error doesn't mean much for HMMs since their states are fluid
sign_error_table[unlist(lapply(all_model_structures, function(x) dim(x)[2] == 4)),] <- NA
model_table[,7:9] <- 1- round(sign_error_table, 3)

# fill out the 10:12 columns. measures of mean squared error 
rmse_table <- do.call(rbind, lapply(param_model_res, function(x) getRMSEList(x)))
model_table[,10:12] <- round(rmse_table, 3)

# fill out the 13:15 columns. measures of power 
model_power_files <- dir("model_power_test/", full.names = TRUE) # simulate and fit under model i
model_power_res <- lapply(all_model_names, function(x) load_model(x, model_power_files))
names(model_power_res) <- paste0("M", model_numbers) # each element is a specific model
model_power_list <- list()
model_type_vector <- model_table[,2]
for(i in model_numbers){
  focal_model <- paste0("M", i)
  model_power_table <- do.call(rbind, lapply(model_power_res[[i]], function(x) getPowerMeasures(x, focal_model, model_type_vector)))
  model_power_table <- model_power_table[!apply(model_power_table, 1, function(x) any(is.na(x))),]
  model_power_list[[i]] <- data.frame(prop_best=sum(model_power_table[,1])/length(model_power_table[,1]), dAIC=mean(model_power_table[,2]), prop_cd=sum(model_power_table[,3]=="CD")/length(model_power_table[,3]))
}
model_table[,13:15] <- round(do.call(rbind, model_power_list), 3)

# fill out the 16 column. measures of expected value error 
param_model_res$M1[[1]]$simulated_data$data[,3]





write.csv(model_table, file = "~/Desktop/houwie_model_table.csv")

