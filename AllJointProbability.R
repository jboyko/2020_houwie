# script for generating all the possible underlying mappings and looking at joint probablity. using this we can look at the bias produced by looking only at the discrete mappings.
getAllJointProbs <- function(phy, data){
  phy <- reorder.phylo(phy, "pruningwise")
  edge_liks_list <- getEdgeLiks(phy, data[,c(1,2)], 2, 1, 1.1)
  nodes_to_fix <- unique(phy$edge[,1])
  possible_internal_combinations <- expand.grid(rep(list(1:2), length(nodes_to_fix)))
  disc_liks <- internode_maps <- houwie_simmaps <- list()
  conditional_probs <- getConditionalInternodeLik(phy, Q, edge_liks_list)
  for(h in 1:dim(possible_internal_combinations)[1]){
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
    tmp <- getInternodeMap(phy, Q, edge_liks_list_i, root_state, c(0.5,0.5), 1)
    disc_liks[[h]] <- tmp[[1]]$llik
    internode_maps[[h]] <- tmp[[1]]$map
    houwie_simmaps[[h]] <- getMapFromSubstHistory(list(tmp[[1]]$map), phy)[[1]]
  }
  
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
  joint_probability_table <- matrix(NA, length(internode_maps), 3, dimnames = list(1:length(houwie_simmaps), c("disc", "cont", "total")))
  for(i in 1:length(houwie_simmaps)){
    llik_discrete <- unlist(disc_liks)[i]
    llik_continuous <- OUwie.basic(houwie_simmaps[[i]], data, simmap.tree=TRUE, alpha = alpha, sigma.sq = sig2, theta = theta)
    joint_probability_table[i,] <- c(llik_discrete, llik_continuous, llik_discrete + llik_continuous)
  }
  return(list(plotting_simmaps = plotting_simmaps, joint_probability_table = joint_probability_table))
}

setwd("~/2020_houwie/")
require(geiger)
require(corHMM)
require(OUwie)
require(partitions)
require(expm)
require(MASS)
require(phytools)
source("hOUwieNode.R")

nTip <- 5
rate <- 1
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
Tmax <- max(branching.times(phy))
Q <- matrix(c(-rate, rate,rate,-rate),2,2)
alpha = c(1,1)
sig2 = c(1,5)
theta0 = 5
theta = c(5,10)
houwie_data <- hOUwie.sim(phy, Q, c(1,0), alpha, sig2 = sig2, theta0, theta)
data <- houwie_data$data
true_map <- houwie_data$simmap[[1]]
plot(true_map)

Q <- Q/100
out <- getAllJointProbs(phy, data)

plotting_simmaps <- out$plotting_simmaps
data <- as.data.frame(out$joint_probability_table)

par(mfrow=c(4,4))
for(i in 1:length(plotting_simmaps)){
  plotSimmap(plotting_simmaps[[i]], fsize = 0.1, outline = TRUE, ylim = c(0, nTip + 2))
  llik_i <- round(data[i,], 3)
  name_i <- names(llik_i)
  llik_text <- paste0(name_i[1], ": ", llik_i[1], "\n",
                      name_i[2], ": ", llik_i[2], "\n",
                      name_i[3], ": ", llik_i[3])
  text(x = 0, y = 7, llik_text, adj = c(0, 1), font = 3)
  # text(x = 0.55, y = 5.8, name_text, adj = 1)
}

dev.off(); plot(data$disc, data$cont, cex=30+data[,3], pch=20)
text(label=rownames(data), data[,1], data[,2], col="red")


library(plyr)
library(ggplot2)
data$disc_prob <- exp(data$disc)
data$disc_weight <- data$disc_prob/sum(data$disc_prob)
results <- data.frame()
for (nmodels in sequence(nrow(data)-1)) {
  for (rep in sequence(20)) {
    for (type in sequence(2)) {
      if(type==1) {
        focal_rows <- sample(sequence(nrow(data)), size=nmodels, prob=data$disc_weight, replace=FALSE)
      } else {
        focal_rows <- sample(sequence(nrow(data)), size=nmodels, prob=NULL, replace=FALSE)
      }
      results <- plyr::rbind.fill(results, data.frame(nmodel=nmodels, type=ifelse(type==1,"discrete-weighted", "flat-weighted"), total_lnl=log(sum(exp(data$total[focal_rows])))))
    }
  }
}
ggplot(results, aes(x=nmodel, y=total_lnl)) + geom_point(aes(colour=factor(type), alpha=0.3)) + facet_wrap(~type)


