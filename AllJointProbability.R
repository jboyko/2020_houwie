# script for generating all the possible underlying mappings and looking at joint probablity. using this we can look at the bias produced by looking only at the discrete mappings.
setwd("~/2020_houwie/")
require(geiger)
require(corHMM)
require(partitions)
require(expm)
require(MASS)
require(phytools)
source("hOUwieNode.R")

nTip <- 3
rate <- 1
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
Tmax <- max(branching.times(phy))
Q <- matrix(c(-rate, rate,rate,-rate),2,2)
nState <- dim(Q)[1]
data <- sim.char(phy, Q, 1, "discrete")[,,1]
data <- data.frame(sp = names(data), d = data)
data <- cbind(data, rnorm(nTip, 10))

phy <- reorder.phylo(phy, "pruningwise")
edge_liks_list <- getEdgeLiks(phy, data[,c(1,2)], 2, 1, 1.1)
nodes_to_fix <- unique(phy$edge[,1])
possible_internal_combinations <- expand.grid(rep(list(1:2), length(nodes_to_fix)))
disc_liks <- internode_maps <- houwie_simmaps <- list()

for(h in 1:dim(possible_internal_combinations)[1]){
  print(h)
  combo_h <- as.numeric(possible_internal_combinations[h,])
  edge_liks_list_i <- edge_liks_list
  for(i in 1:length(nodes_to_fix)){
    node_i <- nodes_to_fix[i]
    anc_edges_to_fix <- which(phy$edge[,1] == node_i)
    dec_edges_to_fix <- which(phy$edge[,2] == node_i)
    state_j <- combo_h[i]
    fix_vector <- numeric(2)
    fix_vector[state_j] <- 1
    for(k in dec_edges_to_fix){
      edge_liks_list_i[[k]][1,] <- fix_vector
    }
    for(k in anc_edges_to_fix){
      last_row <- dim(edge_liks_list_i[[k]])[1]
      edge_liks_list_i[[k]][last_row,] <- fix_vector
    }
  }
  root_state <- numeric(2)
  root_state[combo_h[length(combo_h)]] <- 1
  tmp <- getInternodeMap(phy, Q, edge_liks_list_i, root_state, c(0.5, 0.5), 1)
  disc_liks[[h]] <- tmp[[1]]$llik
  internode_maps[[h]] <- tmp[[1]]$map
  houwie_simmaps[[h]] <- getMapFromSubstHistory(list(tmp[[1]]$map), phy)[[1]]
}

log(sum(exp(unlist(disc_liks))))
corHMM(phy, data[,c(1,2)], 1, model = "ER", p = rate)$loglik


plotting_simmaps <- houwie_simmaps
for(h in 1:length(plotting_simmaps)){
  simmaps <- plotting_simmaps[[h]]
  for(i in 1:length(simmaps$maps)){
    names_i <- names(simmaps$maps[[i]])
    if(length(names_i) > 1){
      if(names_i[1] == names_i[2]){
        simmaps$maps[[i]] <- sum(simmaps$maps[[i]])
        names(simmaps$maps[[i]]) <- names_i[1]
      }
    }
  }
  plotting_simmaps[[h]] <- simmaps
}

par(mfrow=c(2,2))
for(i in 1:length(plotting_simmaps)){
  plot(plotting_simmaps[[i]])
}

joint_probability_table <- matrix(NA, length(internode_maps), 2, dimnames = list(1:length(houwie_simmaps), c("llik_disc", "llik_cont")))
for(i in 1:length(houwie_simmaps)){
  llik_discrete <- unlist(disc_liks)[i]
  llik_continuous <- OUwie.basic(houwie_simmaps[[i]], data, simmap.tree=TRUE, alpha = c(0.1,1), sigma.sq = c(1,1), theta = c(10,10))
  joint_probability_table[i,] <- c(llik_discrete, llik_continuous)
}

log(sum(exp(joint_probability_table[,1]))) + log(0.5)
corHMM(phy, data[,c(1,2)], 1, model = "ER", p = rate)



