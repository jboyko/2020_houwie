setwd("~/2020_hOUwie/")

source("hOUwieNode.R")
source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(expm)
require(ggplot2)
require(reshape2)

#### #### #### #### #### #### #### #### #### #### #### #### 
# functions
#### #### #### #### #### #### #### #### #### #### #### #### 

summarizeFile <- function(file){
  load(file)
  true_pars <- obj$sim$pars
  model_pars <- obj$fit$p
  names(true_pars) <- names(model_pars) <- c("rate", "alpha", "sigma_sq", "theta_1", "theta_2")
  optim_pars <- strsplit(gsub(".Rsave", "", gsub(".*//", "", file)), "-")[[1]]
  names(optim_pars) <- c("algorithm", "n_starts", "init_pars")
  run_time <- obj$fit$run_time
  units(run_time) <- "mins"
  model_pars
  rmse <- log(true_pars) - log(model_pars)
  out <- data.frame(algorithm = optim_pars[1], n_starts_init = paste0(optim_pars[2], "_", optim_pars[3]),
                    rate = rmse[1], alpha = rmse[2], sigma_sq = rmse[3], theta_1 = rmse[4], theta_2 = rmse[5]
                    , row.names = NULL)
  return(out)
}

#### #### #### #### #### #### #### #### #### #### #### #### 
# running
#### #### #### #### #### #### #### #### #### #### #### #### 

files <- dir("optim_test/", full.names = TRUE)
results <- do.call(rbind, lapply(files, summarizeFile))
melted_results <- melt(results)
colnames(melted_results)
ggplot(melted_results, aes(x = algorithm, y = value, fill = n_starts_init)) + 
  geom_boxplot() + 
  ylab("log diff from true") + 
  facet_wrap(~variable)


boxplot(results$run_time ~ as.factor(paste0(results$algorithm, "_", results$n_starts_init)), xlab = "algorithm", ylab = "run_time (min)")

