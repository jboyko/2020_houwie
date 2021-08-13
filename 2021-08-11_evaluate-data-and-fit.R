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
require(dplyr)
require(reshape)

files_fits <- dir("sim_fits/", full.names = TRUE)
files_data <- dir("sim_data/", full.names = TRUE)

i = 1
# for a particular fit file
fit_path <- files_fits[685]

load(fit_path)
# load the fit
load(fit_path) # out
model_file <- gsub(".*/", "", fit_path)
model_name <- strsplit(model_file, "-")[[1]][1]
names(out) <- paste0("M", 1:24)

# find and load the simulating data
data_path <- files_data[grep(model_file, files_data)]
load(data_path) # full_data
tmp <- out[[which(names(out) == model_name)]]

hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 1.1, 50, p = full_data$pars, recon=FALSE)

getComparisonList <- function(fit_path, model_average=FALSE){
  # load the fit
  load(fit_path) # out
  model_file <- gsub(".*/", "", fit_path)
  model_name <- strsplit(model_file, "-")[[1]][1]
  names(out) <- paste0("M", 1:24)
  
  # find and load the simulating data
  data_path <- files_data[grep(model_file, files_data)]
  load(data_path) # full_data
  
  # organize the generating data and parameters
  true_generating_parameters <- getModelParams(organizeGeneratingData(full_data))
  
  if(model_average){
    # organize the fit 
    class_houwie <- unlist(lapply(out, function(x) class(x) == "houwie"))
    if(!all(class_houwie)){
      cat("try-error detected in ", fit_path, "\n")
    }
    out <- out[class_houwie]
    parameter_list <- lapply(out, getModelParams)
    
    # model average the fits
    model_average_table <- getModelTable(out)
    AIC_wts <- model_average_table$AICwt
    model_averaged_parameters <- matrix(rowSums(mapply("*", parameter_list, AIC_wts)),dim(true_generating_parameters)[1],dim(true_generating_parameters)[2],dimnames=list(row.names(true_generating_parameters),colnames(true_generating_parameters)))
    
    # output a table which compares the two
    comparison_list <- list(true_generating_parameters = true_generating_parameters,
                            model_averaged_parameters = model_averaged_parameters,
                            model_average_table = model_average_table)
    return(comparison_list)
  }else{
    true_model_est <- out[[which(names(out) == model_name)]]
    log_diff <- colMeans(log(true_generating_parameters) - log(getModelParams(true_model_est)))
    return(log_diff)
  }
}

pars_list <- list()
tri_index <- lower.tri(matrix(NA, 100, 100))
models <- paste("M", 1:24, sep = "")
# models <- CIDx_models
for(i in 1:length(models)){
  model_name <- models[i]
  print(model_name)
  model_files <- files_fits[grep(paste(model_name, "-", sep=""), files_fits)]
  pars_list[[i]] <- do.call(rbind, lapply(model_files, getComparisonList))
}
names(pars_list) <- models

boxplot(pars_list$M3)

# sign_table <- melt(do.call(rbind, out), id="model")
# sign_table[,1] <- factor(sign_table[,1], levels = models)
# ggplot(sign_table, aes(x = variable, y = value)) + 
#   geom_boxplot() +
#   labs(x = "model", y = "PropCorrect")


# getPropCorrect <- function(focal_run){
#   diff_table_est <-  diff_table_true <- matrix(NA, 4950, 4)
#   count <- 1
#   true_par <- log(focal_run$true_generating_parameters)[1:100,]
#   est_par <- log(focal_run$model_averaged_parameters)[1:100,]
#   for(i in 1:99){
#     for(j in (i+1):100){
#       diff_table_true[count, ] <- true_par[i, ] - true_par[j, ]
#       diff_table_est[count, ] <- est_par[i, ] - est_par[j, ]
#       count <- count + 1
#     }
#   }
#   # true_sd <- apply(log(focal_run$true_generating_parameters)[1:100,], 2, sd)
#   # est_sd <- apply(log(focal_run$model_averaged_parameters)[1:100,], 2, sd)
#   focal_prop_correct <- 1 - colMeans(abs((diff_table_true >= 0) - (diff_table_est >= 0)))
#   # focal_prop_correct <- colMeans(sqrt((true_par - est_par)^2))
#   #focal_prop_correct <- colMeans(diff_table_true - diff_table_est)
#   names(focal_prop_correct) <- colnames(focal_run$true_generating_parameters)
#   return(focal_prop_correct)
# }
# 
# out <- list()
# for(i in 1:length(models)){
#   print(i)
#   focal_model <- pars_list[[i]]
#   prop_correct_model <- do.call(rbind, lapply(focal_model, getPropCorrect))
#   out[[i]] <- data.frame(model = models[i], prop_correct_model)
# }
# 
# sign_table <- melt(do.call(rbind, out), id="model")
# sign_table[,1] <- factor(sign_table[,1], levels = models)
# ggplot(sign_table, aes(x = variable, y = value)) + 
#   geom_boxplot() +
#   labs(x = "model", y = "PropCorrect")







