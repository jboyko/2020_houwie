
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

files_fits <- dir("sim_fits/", full.names = TRUE)
files_data <- dir("sim_data/", full.names = TRUE)

CD_models <- paste("M", c(3:12), sep = "")
CID_models <- paste("M", c(1,2,13,14), sep = "")
CIDx_models <- paste("M", c(15:24), sep = "")

i = 1
# for a particular fit file
fit_path <- files_fits[i]

getComparisonList <- function(fit_path){
  # load the fit
  load(fit_path) # out
  model_name <- gsub(".*/", "", fit_path)
  
  # find and load the simulating data
  data_path <- files_data[grep(model_name, files_data)]
  load(data_path) # full_data
  
  # organize the generating data and parameters
  true_generating_parameters <- getModelParams(organizeGeneratingData(full_data))
  
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
}

require(viridis)

model_name <- CD_models[1]
model_files <- files_fits[grep(model_name, files_fits)]

pars_list <- lapply(model_files, getComparisonList)
RMSE <- do.call(rbind, lapply(pars_list, function(x) colMeans(sqrt((x[[1]][1:100,] - x[[2]][1:100,])^2))))

plot(x = pars_list$true_generating_parameters[,column_number], 
     y = pars_list$model_averaged_parameters[,column_number], 
     ylim = c(0, 10), xlim = c(0, 10), col = cols[1], ylab = "estimated", xlab = "generating", main = colnames(pars_list$true_generating_parameters)[column_number], pch = 16)
for(i in 2:length(files_fits)){
  pars_list <- getComparisonList(files_fits[i])
  points(x = pars_list$true_generating_parameters[,column_number], 
         y = pars_list$model_averaged_parameters[,column_number],
         col = cols[i], pch = 16)
}





