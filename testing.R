# establish the Mk likelihood

require(corHMM)
require(ape)
require(geiger)
require(expm)

avgLogs <- function(x){
  x <- x[!is.na(x)]
  out <- max(x) + log(sum(exp(x - max(x))))
  return(out)
}

phy <- rcoal(3)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
nTip <- length(phy$tip.label)
dat <- c("t1" = 1, "t2" = 2, "t3" = 1)

tmp <- fitDiscrete(phy, dat)
lik <- tmp$lik

L_z <- array(data=-Inf, dim = c(2, 2, dim(phy$edge)[1]))

mm = rbind(c(-1,1), c(1,-1))
for(tip_index in 1:nTip){
  focal_edge_index <- match(tip_index, phy$edge[,2])
  time_edge <- phy$edge.length[focal_edge_index]
  edge_Mk <- matrix(0, 2, 2)
  for(i in 1:2){
    Mk_vec <- rep(0, 2)
    Mk_vec[i] <- 1
    edge_Mk[i,] <- log(c(expm(mm * time_edge) %*% Mk_vec))
  }
  for(i in 1:2){
    Mk_llik <- edge_Mk[i,dat[tip_index]]
    L_z[i,dat[tip_index],focal_edge_index] <- Mk_llik
  }
}

pruningwise.index <- unique(phy$edge[reorder(phy, "pruningwise", index.only = TRUE), 1])
for(node_index in pruningwise.index[-length(pruningwise.index)]){
  focal_anc_index <- which(phy$edge[,2] %in% node_index)
  focal_dec_index <- which(phy$edge[,1] %in% node_index)
  time_edge <- phy$edge.length[focal_anc_index]
  edge_Mk <- matrix(0, 2, 2)
  for(i in 1:2){
    Mk_vec <- rep(0, 2)
    Mk_vec[i] <- 1
    edge_Mk[i,] <- log(c(expm(mm * time_edge) %*% Mk_vec))
  }
  for(i in 1:2){
    for(j in 1:2){
      Mk_llik <- edge_Mk[i,j]
      # L_z[i,j,focal_anc_index] <- Mk_llik + sum(apply(L_z[j,,focal_dec_index], 2, avgLogs))
      L_z[i,j,focal_anc_index] <- Mk_llik + sum(apply(L_z[j,,focal_dec_index], 2, max))
    }
  }
}


focal_dec_index <- which(phy$edge[,1] %in% (nTip+1))
L_root <- vector("numeric", 2)
for(k in 1:2){
  #P_k <- root.p[state_k]
  #L_root[k] <- log(P_k) + sum(apply(L_z[k,,focal_dec_index], 2, max))
  L_root[k] <- sum(apply(L_z[k,,focal_dec_index], 2, max))
  # L_root[k] <- sum(apply(L_z[k,,focal_dec_index], 2, avgLogs))
  
}
# print(L_root)
# llik <- max(L_root)
llik <- max(L_root) + log(sum(exp(L_root - max(L_root))))

lik(1, "flat")

# we're good on the Mk end...


# establish the 
