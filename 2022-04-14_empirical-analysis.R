setwd("2020_houwie/")
source("hOUwieNode.R")
source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)
require(data.table)
require(reshape2)
require(ggplot2)
require(ggplotify)
require(gridExtra)
require(ggtree)
require(aplot)

phy <- read.tree("empirical/Ericaceae_Schwery_etal_2015.tre")

model_files <- dir("empirical_fit/", full.names = TRUE)
model_files <- model_files[-grep("bootstrap_res", model_files)]
out <- list()
for(i in 1:length(model_files)){
  print(i)
  load(model_files[i])
  out[[i]] <- fit
}
names(out) <- gsub(".Rsave", "", gsub(".*=", "", model_files))
desired_order <- c("CID_BM1", "CID_OU1", paste0("CD_", c("BMV", "OUV", "OUA", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA", "OUBM1", "OUBMV")), paste0("CID+_", c("BMV", "OUV", "OUA", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA", "OUBM1", "OUBMV")), "HYB_BMS", "HYB_OUM", "HYB_OUMVA")
out <- out[match(desired_order, names(out))]

# OUwie:::getModelAvgParams(out)

empirical_model_table <- getModelTable(out)
empirical_model_table <- round(empirical_model_table, 2)
empirical_model_table <- cbind(do.call(rbind, strsplit(rownames(empirical_model_table), "_")),empirical_model_table)
colnames(empirical_model_table)[1:2] <- c("Model class", "Model type")
write.csv(empirical_model_table, file = "tables/empirical-model-table.csv", row.names = FALSE)

# model averaged parameters
avg_pars <- getModelAvgParams(out[empirical_model_table$AICwt > 0])
avg_pars$mod_avg_disc
avg_pars$mod_avg_cont
exp(avg_pars$mod_avg_cont[3,]) * 0.0001 

# stationaary variance
avg_pars$mod_avg_cont[2,1]/ (2 * avg_pars$mod_avg_cont[1,1])
avg_pars$mod_avg_cont[2,2]/ (2 * avg_pars$mod_avg_cont[1,2])

# half=life
log(2)/avg_pars$mod_avg_cont[1,1]
log(2)/avg_pars$mod_avg_cont[1,2]

# half-life proportion
log(2)/avg_pars$mod_avg_cont[1,1]/max(branching.times(phy))
log(2)/avg_pars$mod_avg_cont[1,2]/max(branching.times(phy))
sum(phy$edge.length)

# expected value versus observed value plots
phy <- out[[1]]$phy
true_dat <- out[[1]]$data[,3]
names(true_dat) <- out[[1]]$data[,1]

true_dat <- exp(true_dat) * 0.0001
expc_dat <- exp(avg_pars$mod_avg_expc) * 0.0001

# true_dat <- out[[1]]$data[,3]
# expc_dat <- avg_pars$mod_avg_expc

(mean(true_dat - expc_dat))
(range(true_dat - expc_dat)) 
hist(true_dat - expc_dat)

# model avg model fit
avg_pars$mod_avg_disc
avg_pars$mod_avg_cont

cont_model <- out$CD_OUMVA$index.cont
disc_model <- out$CD_OUMVA$index.disc
dat <- out$CD_OUMVA$data
dat$reg <- as.factor(dat$reg)
p <- c(na.omit(c(avg_pars$mod_avg_disc)), c(t(avg_pars$mod_avg_cont)))
model_avg_houwie_fit <- hOUwie(phy, dat, 1, disc_model, cont_model, nSim = 100, sample_nodes = TRUE, adaptive_sampling = TRUE, p = p, mserr = "known")


## parameteric bootstrap
mappings <- model_avg_houwie_fit$simmaps
joint_liks_per_map <- model_avg_houwie_fit$all_cont_liks + model_avg_houwie_fit$all_disc_liks
# map <- mappings[[1]]
sim_based_on_map <- function(map, pars){
  root_state <- getRootState(map)
  map <- correct_map_edges(map)
  map <- reorder(map, "cladewise")
  sim_dat <- OUwie.sim(phy = map, simmap.tree = TRUE, alpha = pars[1,], sigma.sq = pars[2,], theta0 = pars[3,root_state], theta = pars[3,])
  return(sim_dat)
}

sim_based_on_mappings <- function(mappings, pars, joint_liks_per_map){
  weights <- exp(joint_liks_per_map - max(joint_liks_per_map))/sum(exp(joint_liks_per_map - max(joint_liks_per_map)))
  sim_set <- lapply(mappings[weights > 0], function(x) sim_based_on_map(x, pars))
  all_values <- do.call(rbind, lapply(sim_set, function(x) x[,2]))
  weighted_values <- colSums(all_values * weights[weights > 0])
  names(weighted_values) <- mappings[[1]]$tip.label
  return(weighted_values)
}

refit_houwie <- function(continuous_value, houwie_model){
  dat <- houwie_model$data
  dat[,3] <- continuous_value
  fit <- hOUwie(model_avg_houwie_fit$phy, dat, rate.cat = 1, discrete_model = model_avg_houwie_fit$discrete_model, continuous_model = model_avg_houwie_fit$index.cont, time_slice = model_avg_houwie_fit$time_slice, nSim = model_avg_houwie_fit$nSim, ip = model_avg_houwie_fit$p, sample_nodes = TRUE, adaptive_sampling = TRUE, mserr = "known")
  return(fit)
}


# continuous_values <- mclapply(1:100, function(x) sim_based_on_mappings(mappings, avg_pars$mod_avg_cont, joint_liks_per_map), mc.cores = 50)
# 
# bootstrap_res <- mclapply(continuous_values, function(x) refit_houwie(x, model_avg_houwie_fit), mc.cores = 50)
# save(bootstrap_res)

## sumarizning the boot strap res
load("empirical_fit/bootstrap_res.Rsave")
bootstrap_res[unlist(lapply(bootstrap_res, function(x) x$solution.cont[3,2] < 6))]

bootstrap_estimates <- do.call(rbind, lapply(bootstrap_res, function(x) x$p))
bootstrap_estimates[,c(7,8)] <- exp(bootstrap_estimates[,c(7,8)]) * 0.0001
# apply(bootstrap_estimates, 2, sort)
generating_p <- model_avg_houwie_fit$p
generating_p[c(7,8)] <- exp(generating_p[c(7,8)]) * 0.0001
std_errors <- sqrt(rowMeans((t(bootstrap_estimates) - generating_p)^2))

# expc_dat <- exp(out$`FitSD=HYB_OUM`$expected_vals) * 0.0001
# cbind(exp(out$`FitSD=HYB_OUM`$expected_vals) * 0.0001, exp(out$`FitSD=CD_OUMVA`$expected_vals) * 0.0001)
expc_dat <- exp(model_avg_houwie_fit$expected_vals) * 0.0001
true_dat <- exp(model_avg_houwie_fit$data[,3]) * 0.0001

layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
layout(mat = layout.matrix,widths = c(2,1.5))
cols_per_state <- c("#a6611a", "#018571", "#dfc27d", "#80cdc1")
disc_dat <- out[[1]]$data[,2]
names(disc_dat) <- out[[1]]$data[,1]
cols <- cols_per_state[as.numeric(as.factor(disc_dat))]
cols_b <- cols_per_state[3:4][as.numeric(as.factor(disc_dat))]


Tmax <-  max(branching.times(phy))
Ntaxa <- length(phy$tip.label)
offset <- 0.5
# plot.phylo(phy, show.tip.label = FALSE, edge.width = 0.5)
map <- correct_map_edges(model_avg_houwie_fit$simmaps[[1]])
cols_map <- cols_per_state[1:2]
names(cols_map) <- 1:2

a <- ggtree(map) +
  ggtitle("a) Ericaceae Phylogeny")

plot_data <- data.frame(id = phy$tip.label, value = true_dat, fruit_type = as.factor(disc_dat))
b <- ggplot(plot_data, aes(x = id, y = value, color = fruit_type)) + 
  geom_col(width = .000001) + 
  scale_color_manual(values = cols_per_state[1:2]) +
  coord_flip() + 
  theme_tree2() + 
  ylab("Aridity Index") +
  ggtitle("b) Observed AI") +
  theme(legend.position='none')

ab <- b %>% insert_left(a, width = 3)  
ab <- as.grob(ab)
# a <- as.grob(~plotSimmap(map, fsize = .0001, colors = cols_map, lwd = .75))
# tiplabels(pch = 16, col = cols, cex = 0.5, offset = offset)
plot_dat <- cbind(FruitType = c("Dry", "Fleshy"), melt(exp(model_avg_houwie_fit$solution.cont[3,]) * 0.0001))
c <- ggplot(plot_dat, aes(x = FruitType, y = value, color = FruitType)) +
  scale_color_manual(values = cols_per_state[1:2]) + 
  geom_point(size = 5, shape = 21) +
  geom_errorbar(aes(ymax = value + std_errors[c(7,8)], ymin = value - std_errors[c(7,8)]), width = 0.2) + 
  ylab("P/PET") + 
  xlab("Fruit Type") + 
  ggtitle("c) Estimated optimum") + 
  theme_classic() +
  theme(legend.position = "none")

plot_dat <- cbind(FruitType = c("Dry", "Fleshy"), melt(model_avg_houwie_fit$solution.cont[2,]))
d <- ggplot(plot_dat, aes(x = FruitType, y = value, color = FruitType)) +
  scale_color_manual(values = cols_per_state[1:2]) + 
  geom_point(size = 5, shape = 21) +
  geom_errorbar(aes(ymax = value + std_errors[c(5,6)], ymin = value - std_errors[c(5,6)]), width = 0.2) + 
  ylab("(P/PET)^2") + 
  xlab("Fruit Type") + 
  ggtitle("d) Estimated sigma squared") + 
  theme_classic() +
  theme(legend.position = "none")

plot_dat <- cbind(FruitType = c("Dry", "Fleshy"), melt(model_avg_houwie_fit$solution.cont[1,]))
e <- ggplot(plot_dat, aes(x = FruitType, y = value, color = FruitType)) +
  scale_color_manual(values = cols_per_state[1:2]) + 
  geom_point(size = 5, shape = 21) +
  geom_errorbar(aes(ymax = value + std_errors[c(3,4)], ymin = value - std_errors[c(3,4)]), width = 0.2) + 
  ylab("(P/PET)/Time") + 
  xlab("Fruit Type") + 
  ggtitle("e) Estimated alpha") + 
  theme_classic() +
  theme(legend.position = "none")

estimate_plots <- grid.arrange(c, d, e)

final <- grid.arrange(ab, estimate_plots, nrow = 1, widths = c(3,1.25))
ggsave(filename = "figures/raw/empirical_figure.pdf", plot = final, height = 7, width = 10, units = "in")


# evaluating elevation

library(maptools) 
library(raster)

layer <- raster("empirical/wc2-5/alt.bil")
points <- read.csv("empirical/Ericaceae_allpoints.csv")


source("~/rangers/R/miscellaneous.R") 

points_raw <- points[,c(2,8,7)]
points_ai <- points[,c(2,6,8,7)]
all_biomes <- localityToBiome(points_raw)
summary_biomes <- getBiomes(all_biomes)

climate_by_point <- DataFromPoints(points_raw, layer)
climate_dataset_alt <- GetClimateSummStats(climate_by_point)

climate_dataset_ai <- GetClimateSummStats(points_ai)

cols <- c("brown", "green")[as.numeric(model_avg_houwie_fit$data[,2][match(gsub(" ", "_", climate_dataset_ai$species), model_avg_houwie_fit$data[,1])])]

plot(climate_dataset_ai$mean, climate_dataset_alt$mean, col = cols)





# b <- ggplot(plot_data, aes(x = sp, y = value)) +

xlims <- range(c(true_dat, expc_dat))
# segments(expc_dat, 1:Ntaxa, true_dat, 1:Ntaxa, lwd=0.5, col = cols)
for(i in 1:dim(many_vectors)[1]){
  points(exp(many_vectors[i,])*0.0001, 1:Ntaxa, pch=16, cex = 0.2, col = rgb(0,0,1,0.05))
}
points(true_dat, 1:Ntaxa, pch=16, cex = 0.25, col = cols)
points(expc_dat, 1:Ntaxa, pch=16, cex = 0.25, col = "black")




# model_avg_reconstruction <- hOUwieRecon(model_avg_houwie_fit, nodes = "internal")






generate_plot <- function(dat, many_vectors){
  random_species <- sample(dat[,1], 2, replace = FALSE)
  simulated_values <- many_vectors[,match(random_species, colnames(many_vectors))]
  true_values <- dat[match(random_species, dat[,1]), 3]
  plot(x = simulated_values[,1], y = simulated_values[,2], xlim = range(c(c(simulated_values, true_values))), ylim = range(c(c(simulated_values, true_values))), xlab = random_species[1], ylab = random_species[2])
  points(x = true_values[1], y = true_values[2], col="red", pch = 16)
}

par(mfrow=c(3,3))
sapply(1:9, function(x) generate_plot(dat, many_vectors))

plot(many_vectors[,3], many_vectors[,100])
points(dat[4,3], dat[100,3], col = "red")




a <- ggtree(map) +
  geom_tippoint() + 
  geom_tiplab(size = 1)

diff_from_optima = true_dat #- exp(model_avg_houwie_fit$solution.cont[3,][as.numeric(as.factor(disc_dat))]) * 0.0001
plot_data <- data.frame(id = phy$tip.label, value = diff_from_optima, fruit_type = as.factor(disc_dat))
b <- ggplot(plot_data, aes(x = id, y = value, color = fruit_type)) + 
  geom_col(width = .000001) + 
  scale_color_manual(values = cols_per_state[1:2]) +
  coord_flip() + 
  theme_tree2() + 
  ylab("Diff from optimum AI") +
  ggtitle("b) Diff from optimum AI") +
  theme(legend.position='none')

ab <- b %>% insert_left(a, width = 2)  
ab <- as.grob(ab)
ggsave(file = "~/Desktop/big-eric-phy.pdf", ab, height = 30, width = 20, units = "in")

model_avg_houwie_fit$data[model_avg_houwie_fit$data[,1] == "Epacris_lanuginosa",]
model_avg_houwie_fit$data[model_avg_houwie_fit$data[,1] == "Pentachondra_involucrata",]
model_avg_houwie_fit$data[model_avg_houwie_fit$data[,1] == "Dracophyllum_trimorphum",]

