
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

good_model_structures <- c(1,2,5,7,9,11,13,14,15,16,17,18,19,20,21,22,25,27,29,31,33,34,35,36,37,38,39,40)

files_fits <- dir("sim_fits/", full.names = TRUE)
files_data <- dir("sim_data/", full.names = TRUE)

# for a particular fit file
fit_path <- files_fits[5]

# load the fit
load(fit_path) # out
model_name <- gsub(".*/", "", fit_path)

# find and load the simulating data
data_path <- files_data[grep(model_name, files_data)]
load(data_path) # full_data

# organize the generating data and parameters
true_generating_parameters <- getModelParams(organizeGeneratingData(full_data))

# organize the fit 
out <- out[good_model_structures]
parameter_list <- lapply(out, getModelParams)

# model average the fits
model_average_table <- getModelTable(out)
AIC_wts <- model_average_table$AICwt
model_averaged_parameters <- matrix(rowSums(mapply("*", parameter_list, AIC_wts)),dim(organized_generating_data)[1],dim(organized_generating_data)[2],dimnames=list(row.names(organized_generating_data),colnames(organized_generating_data)))

# output a table which compares the two
comparison_list <- list(true_generating_parameters = true_generating_parameters,
                        model_averaged_parameters = model_averaged_parameters)








