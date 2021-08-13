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

i = 1
# for a particular fit file
fit_path <- files_fits[grep("M21-", files_fits)][10]

load(fit_path)
# load the fit
load(fit_path) # out
model_file <- gsub(".*/", "", fit_path)
model_name <- strsplit(model_file, "-")[[1]][1]
names(out) <- paste0("M", 1:24)
# find and load the simulating data
data_path <- files_data[grep(model_file, files_data)]
load(data_path) # full_data
full_data
tmp <- out[[which(names(out) == model_name)]]
p <- tmp$p
p[1] <- full_data$pars[1]
res_true_mk <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 50, p = p, recon=FALSE)
p <- tmp$p
p[2] <- full_data$pars[2]
res_true_alpha <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 50, p = p, recon=FALSE)
p <- tmp$p
p[3] <- full_data$pars[3]
res_true_sigma <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 50, p = p, recon=FALSE)
p <- tmp$p
p[c(4,5)] <- full_data$pars[c(4,5)]
res_true_thetas <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 50, p = p, recon=FALSE)
p <- tmp$p
p[c(2,3)] <- full_data$pars[c(2,3)]
res_true_sigma_alpha <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 50, p = p, recon=FALSE)
p <- tmp$p
p[c(2,4,5)] <- full_data$pars[c(2,4,5)]
res_true_thetas_alpha <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 50, p = p, recon=FALSE)
p <- tmp$p
p[c(3,4,5)] <- full_data$pars[c(3,4,5)]
res_true_thetas_sigma <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 50, p = p, recon=FALSE)
p <- tmp$p
p[c(2,3,4,5)] <- full_data$pars[c(2,3,4,5)]
res_true_thetas_alpha_sigma <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 50, p = p, recon=FALSE)
res_true_all <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 500, p = full_data$pars, recon=FALSE) 


nodes <- unlist(lapply(full_data$simmap[[1]]$maps[sapply(101:199, function(x) which(x == full_data$simmap[[1]]$edge[,1])[1])], function(y) as.numeric(names(y)[1])))
tmp_phy <- full_data$simmap[[1]]
tmp_phy$node.label <- nodes
test_true_nodes <- OUwie(phy = tmp_phy, data = full_data$data, model = "OUM", algorithm = "three.point", simmap.tree = FALSE)
test_true_map <- OUwie(phy = full_data$simmap[[1]], data = full_data$data, model = "OUM", algorithm = "three.point", simmap.tree = TRUE)
OUwie.fixed(phy = tmp_phy, data = full_data$data, model = "OUM", algorithm = "three.point", alpha = rep(full_data$pars[3],4), sigma.sq = rep(full_data$pars[4],4), theta = c(rep(full_data$pars[5],2), rep(full_data$pars[6],2)))
OUwie.fixed(phy = full_data$simmap[[1]], data = full_data$data, model = "OUM", algorithm = "three.point", simmap.tree = TRUE, alpha = rep(full_data$pars[3],4), sigma.sq = rep(full_data$pars[4],4), theta = c(rep(full_data$pars[5],2), rep(full_data$pars[6],2)))


tbl <- list(optimized=tmp, 
            true_alpha = res_true_alpha, 
            true_sigma = res_true_sigma,
            true_thetas = res_true_thetas,
            true_alpha_sigma = res_true_sigma_alpha,
            true_alpha_theta = res_true_thetas_alpha, 
            true_sigma_theta = res_true_thetas_sigma,
            true_alpha_sigma_theta = res_true_thetas_alpha_sigma,
            true_all = res_true_all)

getModelTable(tbl)

res_nlopt <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, discrete_model_cid, tmp$continuous_model, 0.1, 10, recon=FALSE, optimizer = "nlopt", opts = list("algorithm"="NLOPT_GN_DIRECT_L", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5))
res_sann_rgh <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 0.1, 50, recon=FALSE, optimizer = "sann", opts = list(max.call=10000, smooth=FALSE))
res_sann_smth <- hOUwie(tmp$phy, tmp$data, tmp$rate.cat, tmp$discrete_model, tmp$continuous_model, 1.1, 50, recon=FALSE, optimizer = "sann", opts = list(max.call=1e4, smooth=TRUE))
