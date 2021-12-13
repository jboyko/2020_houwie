#### #### #### #### #### #### #### #### #### #### #### #### 
# imports
#### #### #### #### #### #### #### #### #### #### #### #### 

setwd("~/2020_hOUwie/")

require(OUwie)
require(corHMM)

source("hOUwieNode.R")
source("Utils.R")

load("empirical/eric_results_local.Rsave")

# remove failed optimizations
eric_results <- eric_results[-c(9, 14)]
model_names <- paste0("M", 1:27)[-c(9, 14)]
model_table <- round(getModelTable(eric_results), 2)
rownames(model_table) <- model_names
write.csv(model_table, file = "~/Desktop/model_table_eric_seed.csv")

non_hmm_eric_results <- eric_results[unlist(lapply(eric_results, function(x) x$rate.cat == 1))]

mod_params <- getModelAvgParams(non_hmm_eric_results)

best_fit <- eric_results[[which(rownames(model_table) == "M11")]]

exp(mod_params$mod_avg_cont[3,])*0.0001


hist(mod_params$mod_avg_expc - best_fit$data[,3])

range(best_fit$data[,3])


getModelTable(eric_results[[2]][-12])
getModelTable(eric_results[[3]][-12])
getModelTable(eric_results[[4]])

eric_results[[3]][[21]]

