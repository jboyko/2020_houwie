#### #### #### #### #### #### #### #### #### #### #### #### 
# functions
#### #### #### #### #### #### #### #### #### #### #### #### 
get_parameter_names_discrete <- function(discrete_index, continuous_index){
  all_names <- c("rate", "alpha", "sigma", "theta")
  n_p_disc <- max(discrete_index, na.rm = TRUE)
  n_p_alpha <- length(unique(na.omit(continuous_index[1,])))
  n_p_sigma <- length(unique(na.omit(continuous_index[2,])))
  n_p_theta <- length(unique(na.omit(continuous_index[3,])))
  tmp_names <- c(rep("rate", n_p_disc), rep("alpha", n_p_alpha), rep("sigma", n_p_sigma), rep("theta", n_p_theta))
  out_names <- tmp_names
  for(i in 1:4){
    out_names[tmp_names %in% all_names[i]] <- paste0(all_names[i], "_", seq_len(length(which(tmp_names %in% all_names[i]))))
  }
  return(out_names)
}

get_log_error_from_file <- function(input_file){
  load(input_file)
  par_names <- get_parameter_names_discrete(out$model_res[[1]]$index.disc, out$model_res[[1]]$index.cont)
  log_diff <- log(do.call(cbind, lapply(out$model_res, function(x) x$p))) - log(out$simulated_data$pars)
  rownames(log_diff) <- par_names
  return(log_diff)
}

get_plotting_data <- function(results_by_iter){
  results <- do.call(rbind, results_by_iter)
  melted_results <- melt(results)
  nuiss_pars <- apply(do.call(rbind, strsplit(as.character(melted_results[,2]), "_")), 2, as.numeric)
  real_pars <- apply(do.call(rbind, strsplit(as.character(melted_results[,1]), "_")), 2, as.factor)
  plotting_data <- data.frame(par = real_pars[,1], par_no = real_pars[,2], time_slice = nuiss_pars[,1], n_map = nuiss_pars[,2], log_diff = melted_results[,3])
  return(plotting_data)
}
# input_file <- model_files[175]
get_sign_error_from_file <- function(input_file){
  load(input_file)
  par_names <- get_parameter_names_discrete(out$model_res[[1]]$index.disc, out$model_res[[1]]$index.cont)
  alpha_index <- grep("alpha", par_names)
  sigma_index <- grep("sigma", par_names)
  theta_index <- grep("theta", par_names)
  if(length(alpha_index) > 1){
    true_diff_alpha <- diff(out$simulated_data$pars[alpha_index])
    est_diff_alpha <- unlist(lapply(out$model_res, function(x) diff(x$p[alpha_index])))
    alpha <- est_diff_alpha > 0
  }else{
    alpha <- NA
  }
  if(length(sigma_index) > 1){
    true_diff_sigma <- diff(out$simulated_data$pars[sigma_index])
    est_diff_sigma <- unlist(lapply(out$model_res, function(x) diff(x$p[sigma_index])))
    sigma <- est_diff_sigma > 0
  }else{
    sigma <- NA
  }
  if(length(theta_index) > 1){
    true_diff_theta <- diff(out$simulated_data$pars[theta_index])
    est_diff_theta <- unlist(lapply(out$model_res, function(x) diff(x$p[theta_index])))
    theta <- est_diff_theta > 0
  }else{
    theta <- NA
  }
  error_table <- data.frame(alpha, sigma, theta)
  return(error_table)
}


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


#### #### #### #### #### #### #### #### #### #### #### #### 
# run
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

model_files <- dir("nuiss_par_test/", full.names = TRUE)
model_numbers <- 1:22
all_model_names <- paste0("M", model_numbers, "_")

#### #### #### #### #### #### #### #### #### #### #### #### 
# assessing model performance based on nuissance parameters
#### #### #### #### #### #### #### #### #### #### #### #### 



#### #### #### #### #### #### #### #### #### #### #### #### 
# figure S1 - all log error from many model iterations
#### #### #### #### #### #### #### #### #### #### #### #### 
all_things <- list()
for(i in 1:length(all_model_names)){
  results_by_iter <- lapply(model_files[grep(all_model_names[i], model_files)], get_log_error_from_file)
  plotting_data <- melt(results_by_iter)[,-4]
  plotting_data$par <- as.factor(gsub("*_.*", "", plotting_data$Var1))
  plotting_data$par_no <- as.factor(gsub(".*_", "", plotting_data$Var1))
  plotting_data$n_map <- as.factor(gsub(".*_", "", plotting_data$Var2))
  plotting_data$time_slice <- as.factor(gsub("_.*", "", plotting_data$Var2))
  all_things[[i]] <- cbind(model = gsub("_", "", all_model_names[i]), plotting_data)
}
all_plotting_data <- do.call(rbind, all_things)
all_plotting_data$par <- factor(all_plotting_data$par, levels = c("rate", "alpha", "sigma", "theta"))
all_plotting_data$model <- factor(all_plotting_data$model, levels = paste0("M", model_numbers))

my_plot <- ggplot(all_plotting_data, aes(y = value, x = par_no, fill = par_no)) +
  geom_hline(aes(yintercept=0), colour="#C0C0C0", lwd=1, linetype = "dotted") + 
  geom_violin() +
  scale_fill_manual(values = viridis(3)[c(2,3,1)]) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ylab("log diff est from true") +
  facet_grid(model + n_map ~ par + time_slice, margins = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# my_plot
ggsave("model_parameter_estimation_per_model.pdf", my_plot, "pdf", path = "~/2020_houwie/figures/raw/", width = 12, height = 88, units = "in", limitsize = FALSE)


#### #### #### #### #### #### #### #### #### #### #### #### 
# Type S (sign): the test statistic is in the opposite direction of the true effect size, given that the statistic is statistically significant
#### #### #### #### #### #### #### #### #### #### #### #### 
all_things <- list()
for(i in 1:length(all_model_names)){
  results_by_iter <- lapply(model_files[grep(all_model_names[i], model_files)], get_sign_error_from_file)
  plotting_data <- do.call(rbind, results_by_iter)
  rownames(plotting_data) <- NULL
  all_things[[i]] <- apply(plotting_data, 2, sum)/dim(plotting_data)[1]
}
all_plotting_data <- do.call(rbind, all_things)
model_type <- c("MK", "HMM")[as.numeric(unlist(lapply(all_model_structures, function(x) dim(x)[2] == 4))) + 1]
all_plotting_data <- data.frame(model = gsub("_", "", all_model_names), model_type = model_type, all_plotting_data)
plot_sign_data <- melt(all_plotting_data)

# all_plotting_data$par <- factor(all_plotting_data$par, levels = c("rate", "alpha", "sigma", "theta"))
# all_plotting_data$model <- factor(all_plotting_data$model, levels = paste0("M", model_numbers))
plot_sign_data_MK <- plot_sign_data[plot_sign_data$model_type == "MK",] # sign error doesn't mean anything for HMMs
plot_sign_data_MK$model <- factor(plot_sign_data_MK$model, levels = unique(plot_sign_data[plot_sign_data$model_type == "MK",1]))
ggplot(data = plot_sign_data_MK, aes(x = variable, y = value, fill = variable)) + 
  geom_hline(aes(yintercept=0.5), colour="#C0C0C0", lwd=1, linetype = "dotted") + 
  geom_bar(stat = "identity") +
  theme_bw(base_size = 12) +
  facet_wrap(~model) 


#### #### #### #### #### #### #### #### #### #### #### #### 
# Type M (magnitude or exaggeration ratio): the test statistic in magnitude exaggerates the true effect size, given that the statistic is statistically significant.
#### #### #### #### #### #### #### #### #### #### #### #### 

all_things <- list()
for(i in 1:length(all_model_names)){
  results_by_iter <- lapply(model_files[grep(all_model_names[i], model_files)], get_log_error_from_file)
  plotting_data <- melt(results_by_iter)[,-4]
  plotting_data$Var1 <- as.factor(gsub("*_.*", "", plotting_data$Var1))
  # plotting_data <- do.call(rbind, results_by_iter)
  all_things[[i]] <- cbind(model = gsub("_", "", all_model_names[i]), plotting_data)
}
all_plotting_data <- do.call(rbind, all_things)
all_plotting_data$Var1 <- factor(all_plotting_data$Var1, levels = c("rate", "alpha", "sigma", "theta"))
all_plotting_data$model <- factor(all_plotting_data$model, levels = paste0("M", model_numbers))

ggplot(all_plotting_data, aes(y = Var2, x = Var1)) +
  geom_violin() +
  ylab("log diff est from true") +
  facet_wrap(~model) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

