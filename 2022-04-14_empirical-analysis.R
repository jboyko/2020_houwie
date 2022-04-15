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

phy <- read.tree("empirical/Ericaceae_Schwery_etal_2015.tre")

model_files <- dir("empirical_fit/", full.names = TRUE)
out <- list()
for(i in 1:length(model_files)){
  print(i)
  load(model_files[i])
  out[[i]] <- fit
}
names(out) <- gsub(".Rsave", "", dir("empirical_fit/"))
empirical_model_table <- getModelTable(out)
empirical_model_table <- round(empirical_model_table, 2)
write.csv(empirical_model_table, file = "tables/empirical-model-table.csv")

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

true_dat <- true_dat
expc_dat <- avg_pars$mod_avg_expc

(mean(true_dat - expc_dat))
(range(true_dat - expc_dat)) 
hist(true_dat - expc_dat)

# line - data


disc_dat <- out[[1]]$data[,2]
names(disc_dat) <- out[[1]]$data[,1]
cols <- c("brown", "pink")[as.numeric(as.factor(disc_dat))]


# # Tree
# plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
# legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
# title(main=paste0(group))
# axisPhylo()
# 
# # Temperature
# x <- 1:length(temp)
# plot(temp, x, xlim=range(c(temp-se_temp, temp+se_temp)),
#      pch=19, yaxt = "n", xlab="Temperature (C*10)", ylab="", frame.plot=T, cex=0.3, col=mode_cols)
# suppressWarnings(arrows(temp-se_temp, x, temp+se_temp, x, length=0.05, angle=90, code=3, lwd=0.1, col=mode_cols))
layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
layout(mat = layout.matrix,widths = c(2,1))

Tmax <-  max(branching.times(phy))
Ntaxa <- length(phy$tip.label)
offset <- 0.5
plot.phylo(phy, show.tip.label = FALSE, edge.width = 0.5)
tiplabels(pch = 16, col = cols, cex = 0.5, offset = offset)
legend("bottomleft", legend=c("Dry", "Fleshy"), pt.bg = c("brown", "pink"), pch=21, cex=1, bty = "n")
axisPhylo()

xlims <- range(c(true_dat, expc_dat))
plot(expc_dat, 1:Ntaxa, pch=19, yaxt="n",cex = 0.25, xlim = xlims, frame.plot = FALSE, xlab = "Aridity Index", ylab = "", col = "red")
suppressWarnings(arrows(expc_dat, 1:Ntaxa, true_dat, 1:Ntaxa, length=0.05, code=1, lwd=0.5, col = cols))
points(true_dat, 1:Ntaxa, pch=19, cex = 0.25, col = "blue")
points(expc_dat, 1:Ntaxa, pch=19, cex = 0.25, col = "red")
legend("bottomright", legend=c("Observed", "Expected"), pt.bg = c("blue", "red"), bty = "n", pch = 21)


hist(expc_dat)


avg_pars$mod_avg_disc
avg_pars$mod_avg_cont

cont_model <- out$`FitSD=CD_OUMVA`$index.cont
disc_model <- out$`FitSD=CD_OUMVA`$index.disc
dat <- out$`FitSD=CD_OUMVA`$data
dat$reg <- as.factor(dat$reg)
p <- c(na.omit(c(avg_pars$mod_avg_disc)), c(t(avg_pars$mod_avg_cont)))
model_avg_houwie_fit <- hOUwie(phy, dat, 1, disc_model, cont_model, nSim = 100, sample_nodes = TRUE, adaptive_sampling = TRUE, p = p, mserr = "known")

model_avg_reconstruction <- hOUwieRecon(model_avg_houwie_fit, nodes = c(350))


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
  sim_set <- lapply(mappings, function(x) sim_based_on_map(x, pars))
  weights <- exp(joint_liks_per_map - max(joint_liks_per_map))/sum(exp(joint_liks_per_map - max(joint_liks_per_map)))
  all_values <- do.call(rbind, lapply(sim_set, function(x) x[,2]))
  weighted_values <- colSums(all_values * weights)
  names(weighted_values) <- mappings[[1]]$tip.label
  return(weighted_values)
}


many_tips <- lapply(1:50, function(x) sim_based_on_mappings(stoch_mappings, avg_pars$mod_avg_cont, joint_liks_per_map))
many_vectors <- do.call(rbind, many_tips)

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






