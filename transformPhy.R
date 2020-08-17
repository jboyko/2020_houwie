# This set of functions is responsible for transformiing the branchlengths of an input phylogeny so that when we use the three point algorith to calculate the quadratic products and determinent the interpsecific variance-covariance reflects that assumed by the particular OU model.

# Under "OUStationaryRoot", t is transformed to exp(-2 alpha (T-t)), where T is the mean root-to-tip distance. Under "OUFixedRroot", t is transformed to exp(-2 alpha (T-t)) * (1-exp(-2 alpha t)) per Ho and Ane (2014). I.e. if t (distance from the root to the begining of an epoch) is further away from the root than average (T) then the branch length is increased by a factor alpha.

# functions
# gets the path from a vertex to the root as an index of the edge matrix
getPathToRoot <- function(phy, tip){
  nTip <- length(phy$tip.label)
  root <- nTip + 1
  path <- 0
  count <- 1
  while(tip != root){
    tip.ind <- which(phy$edge[,2] == tip)
    path <- c(path, tip.ind)
    count <- count + 1
    tip <- phy$edge[tip.ind,1]
  }
  path <- path[-1]
  return(path)
}

# transforms the phylogeny
transformPhy <- function(phy, pars){
  # phy must be of class simmap
  nTip <- length(phy$tip.label)
  RootAge <- max(branching.times(phy))
  NodeAges <- branching.times(phy)[phy$edge[,1] - nTip]
  ModMap <- Map <- phy$maps
  D <- V_Tilde <- numeric(dim(phy$edge)[1])
  for(i in 1:dim(phy$edge)[1]){
    # evaluate the map for this particular edge and calculate the tipward variance
    NodeAge_i <- NodeAges[i]
    DistRoot_i <- RootAge - NodeAge_i
    Map_i <- Map[[i]]
    # the age of epoch j starts at the node age
    # EpochAge_j <- NodeAge_i
    Dist_rootward <- DistRoot_i
    w <- v <- 0
    for(j in 1:length(Map_i)){
      # distance the root of epoch j starts at the node distance and ends at node dist + epoch length
      Dist_tipward <- Dist_rootward + Map_i[j]
      # the length of the epoch is scaled by the alpha parameter of that epoch
      Sigma_j <- pars[,2][match(names(Map_i)[j], rownames(pars))]
      Alpha_j <- pars[,3][match(names(Map_i)[j], rownames(pars))]
      # calculate the descendent distance from the root based on a fixed root distribution
      tmp <- Sigma_j * (exp(2 * Alpha_j * Dist_tipward) - exp(2 * Alpha_j * Dist_rootward))/2/Alpha_j
      v <- v + tmp
      ModMap[[i]][j] <- tmp
      w <- w + (Alpha_j * (Dist_tipward - Dist_rootward)) 
      # The new distance from nodes
      Dist_rootward <- Dist_tipward
    }
    V_Tilde[i] <- v
    D[i] <- w
  }
  # calculates the diagonal matrix for each tip i
  DiagWt <- numeric(nTip)
  names(DiagWt) <- phy$tip.label
  for(i in 1:nTip){
    DiagWt[i] <- exp(-sum(D[getPathToRoot(phy, i)]))
  }
  
  phy$edge.length <- V_Tilde
  phy$maps <- ModMap
  obj <- list(tree = phy,
              diag = DiagWt)
  return(obj)
}


# run
require(corHMM)
require(phytools)
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]
MK <- corHMM(phy, data, 1)
phy <- MK$phy
data <- MK$data
model <- MK$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)
states <- MK$states
tip.states <- MK$tip.states
simmap <- makeSimmap(tree=phy, tip.states=tip.states, states=states, model=model, nSim=1, nCores=1)
plotSimmap(simmap[[1]])

# comparing to phylolm
# we need our alpha to be half the sigma to scale the branches in the same way
pars <- matrix(0, 3, 3, dimnames = list(c(1,2,3), c("opt", "sig", "alp")))
pars[,2] <- c(0.02,0.02,0.02)
pars[,3] <- c(0.01,0.01,0.01)

require(phylolm)
Q = rTrait(n=2,MK$phy)
y = rTrait(n=1,MK$phy)
P = cbind(1,y)

tre1 = transformPhy(simmap[[1]], pars)
tre2 = transf.branch.lengths(phy=MK$phy, model="OUfixedRoot", parameters = list(alpha=0.01))
three.point.compute(tre1$tree, P, Q, tre1$diag)
three.point.compute(tre2$tree, P, Q, tre2$diagWeight)


# comparing to OUwie and finding the likelihood
trait <- cbind(MK$data.legend, y)
OUwieRes <- OUwie(phy,trait,model=c("OUMVA"),root.station=FALSE, simmap.tree = TRUE)
debug(OUwie:::varcov.ou)

sum(V_Tilde[getPathToRoot(phy, 1)])

