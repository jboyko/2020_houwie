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

getOUExpectation <- function(phy, pars){
  # requirements
  RootAge <- max(branching.times(phy))
  NodeAges <- branching.times(phy)[phy$edge[,1] - nTip]
  Map <- phy$maps
  nTip <- length(phy$tip.label)
  RootOptim <- pars[1,][match(phy$node.label[1], rownames(pars))]
  X <- numeric(nTip)
  # calulate the expectation
  for(i in 1:nTip){
    # get the tip to root path for tip i
    Path_i <- getPathToRoot(phy, i)
    # get the root to tip mapping for path i
    Map_i <- rev(unlist(Map[Path_i]))
    # for each segment of the map track the expectation
    Dist_rootward <- 0
    TipOpt_i <- AncWt <- 0
    for(j in 1:length(Map_i)){
      Dist_tipward <- Dist_rootward + Map_i[j]
      Optim_j <- pars[,1][match(names(Map_i)[j], rownames(pars))]
      Alpha_j <- pars[,3][match(names(Map_i)[j], rownames(pars))]
      TipOpt_i <- TipOpt_i + Optim_j * (exp(Alpha_j * Dist_tipward) - exp(Alpha_j * Dist_rootward))
      AncWt <- AncWt + (Alpha_j * (Dist_tipward - Dist_rootward)) 
      Dist_rootward <- Dist_tipward
    }
    AncWt <- exp(-AncWt)
    X[i] <- AncWt * RootOptim + AncWt * TipOpt_i
  }
  names(X) <- phy$tip.label
  return(X)
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
simmap <- makeSimmap(tree=phy, data=data, model=model, rate.cat = 1, nSim=1, nCores=1)
#plotSimmap(simmap[[1]])

#### comparing to phylolm
# we need our alpha to be half the sigma to scale the branches in the same way
pars <- matrix(0, 3, 3, dimnames = list(c(1,2,3), c("opt", "sig", "alp")))
pars[,2] <- c(0.02,0.02,0.02)
pars[,3] <- c(0.01,0.01,0.01)

require(phylolm)
y = rTrait(n=1,MK$phy)
Q = y
Q[1:length(y)] = rep(0, length(y))

tre1 = transformPhy(simmap[[1]], pars)
tre2 = transf.branch.lengths(phy=simmap[[1]], model="OUfixedRoot", parameters = list(alpha=0.01))
three.point.compute(tre1$tree, y, Q, tre1$diag)
three.point.compute(tre2$tree, y, Q, tre2$diagWeight)


# comparing phylolm, OUwie, and finding the likelihood
require(OUwie)
trait <- cbind(MK$data.legend, y[as.character(MK$data.legend[,1])])
PhyloLMDat <- trait[,3]
names(PhyloLMDat) <- trait[,1]
sig <- c(0.02,0.02,0.02)
alp <- c(0.01,0.01,0.01)
opt <- c(0, 0, 0)
pars <- matrix(c(opt,sig,alp), 3, 3, dimnames = list(c(1,2,3), c("opt", "sig", "alp")))
nTip <- length(simmap[[1]]$tip.label)


tre1 <- transformPhy(simmap[[1]], pars)
X <- rep(0, length(PhyloLMDat))
names(X) <- names(PhyloLMDat)
comp <- three.point.compute(tre1$tree, PhyloLMDat, X, tre1$diag)


OUwie.fixed(simmap[[1]],trait,model=c("OUM"), simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alp,sigma.sq=sig,theta=opt)
#undebug(OU1d.loglik)
OU1d.loglik(trait=PhyloLMDat, phy=simmap[[1]], model="OUfixedRoot", parameters=list(ancestral.state=opt[1], alpha=alp[1],sigma2=sig[1], optimal.value=opt[1]))
# from phylolm:
-as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + Ntip(phy) * log(sig[1]/2/alp[1]) + 2 * alp[1]/sig[1] * (comp$PP - 2 * comp$QP + comp$QQ))/2
# my expectation would be...
-as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + comp$PP)/2


# adding in multiple alpha and multiple sigma. now only ouwie and 3 point
y = rTrait(n=1,MK$phy)
trait <- cbind(MK$data.legend, y[as.character(MK$data.legend[,1])])
PhyloLMDat <- trait[,3]
names(PhyloLMDat) <- trait[,1]
sig <- c(0.32,0.07,0.02)
alp <- c(0.05,0.01,0.02)
opt <- c(0, 0, 0)
pars <- matrix(c(opt,sig,alp), 3, 3, dimnames = list(c(1,2,3), c("opt", "sig", "alp")))
nTip <- length(simmap[[1]]$tip.label)

tre1 <- transformPhy(simmap[[1]], pars)
X <- as.numeric(opt[MK$data.legend[,2]])
names(X) <- names(PhyloLMDat)
comp <- three.point.compute(tre1$tree, PhyloLMDat, X, tre1$diag)

OUwie.fixed(simmap[[1]],trait,model=c("OUM"), simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alp,sigma.sq=sig,theta=opt)$loglik
# my expectation would be...
-as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + comp$PP)/2


# allowing for different optima means incorporating a different expectation (X) and incorporating some of the other quadratic products
y = rTrait(n=1,MK$phy)
trait <- cbind(MK$data.legend, y[as.character(MK$data.legend[,1])])
PhyloLMDat <- trait[,3]
names(PhyloLMDat) <- trait[,1]
sig <- c(0.32,0.07,0.02)
alp <- c(0.05,0.01,0.02)
opt <- c(1, 0, 1)
pars <- matrix(c(opt,sig,alp), 3, 3, dimnames = list(c(1,2,3), c("opt", "sig", "alp")))
nTip <- length(simmap[[1]]$tip.label)

tre1 <- transformPhy(simmap[[1]], pars)
X <- getOUExpectation(transformPhy(phy, pars)$tree, pars)
#X <- getOUExpectation(simmap[[1]], pars)
comp <- three.point.compute(tre1$tree, PhyloLMDat, X, tre1$diag)

OUwie.fixed(simmap[[1]],trait,model=c("OUM"), simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alp,sigma.sq=sig,theta=opt)$loglik
-as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + comp$PP)/2
-as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2

