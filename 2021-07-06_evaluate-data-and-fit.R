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

# for a particular fit file
fit_path <- files_fits[1]

# load the fit
load(fit_path) # out
model_name <- gsub(".*/", "", fit_path)

# find and load the simulating data
data_path <- files_data[grep(model_name, files_data)]
load(data_path) # full_data

# organize the fit 
parameter_list <- lapply(out, getModelParams)
out[[1]]$recon

# organize the data

# output a table which compares the two

parameter_list <- table_list <- list()
for(i in 1:length(files_fits)){
  load(files_fits[i])
  load(files_data[i])
  names(out) <- paste("M", 1:40, sep = "")
  names_out <- unlist(lapply(out, class))
  out <- out[names_out == "houwie"]
  for(j in 1:length(out)){
    out[[j]]$recon[out[[j]]$recon == 0] <- -Inf
    out[[j]]$recon <- out[[j]]$recon - max(out[[j]]$recon)
    out[[j]]$recon <- round(exp(out[[j]]$recon)/rowSums(exp(out[[j]]$recon)), 10)
  }
  parameter_list[[i]] <- lapply(out, getModelParams)
  table_list[[i]] <- getModelTable(out)
}

recon <- hOUwieRecon(out[[11]])
recon_table <- recon[101:199,]
recon_table <- t(apply(recon_table, 1, function(x) exp(x)/sum(exp(x))))
nodelabels(thermo = recon_table)

out[[9]]$index.cont



full_data
names(out) <- paste("M", 1:40, sep = "")
# out <- out[unlist(lapply(out, function(x) class(x) != "try-error"))]
getModelTable(out)

houwie_obj <- out[[1]]






out[[1]]$DiscLik
unlist(lapply(out[21:40], function(x) x$DiscLik))
Q <- full_data$index.cor
Q[] <- c(full_data$pars[1:3], 0)[full_data$index.cor]
diag(Q) <- -rowSums(Q)
getMapProbability(full_data$simmap[[1]], Q)

cor_data <- full_data$data[,c(1,2)]
cor_data[cor_data[,2] == 3, 2] <- 1
cor_data[cor_data[,2] == 4, 2] <- 2

tmp <- full_data$index.cor
tmp[tmp ==4] <- 0 
rate_1 <- corHMM(full_data$simmap[[1]], data = cor_data, rate.cat = 1, model = "ER")
rate_2 <- corHMM(full_data$simmap[[1]], data = cor_data, rate.cat = 2, rate.mat = tmp)

model <- out[[26]]

model[[1]]$RegimeMap$maps

out[[1]][[1]]$p
simulating_transition_rate <- unlist(lapply(all_data, function(x) x$pars[1]))
estimated_transition_rate <- unlist(lapply(out, function(x) x[[1]]$p[1]))

plot(x = simulating_transition_rate, y = estimated_transition_rate, xlim=c(0.2,0.8), ylim=c(0.2,0.8))
abline(lm(estimated_transition_rate~simulating_transition_rate))
abline(a = 0, b = 1, col = "red")

simulating_transition_number <- unlist(lapply(all_data, function(x)
  sum(unlist(lapply(x$simmap[[1]]$maps, length)) - 1)))
estimated_transition_number <- unlist(lapply(out, function(x) 
  sum(unlist(lapply(x[[1]]$RegimeMap$maps, length)) - 1)))

plot(x = simulating_transition_number, y = estimated_transition_number)
abline(lm(estimated_transition_number~simulating_transition_number))
abline(a = 0, b = 1, col = "red")
