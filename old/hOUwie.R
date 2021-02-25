# functions
## gets the path from a vertex to the root as an index of the edge matrix
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

## takes a node based reconstruction and returns a map (identical to a map from simmap)
getMapFromNode <- function(phy, tipstates, nodestates, shift.point){
  Map <- vector("list", dim(phy$edge)[1])
  Data <- c(tipstates, nodestates)
  NodeStates <- cbind(Data[phy$edge[,1]], Data[phy$edge[,2]])
  for(i in 1:dim(phy$edge)[1]){
    from <- as.character(NodeStates[i,1])
    to <- as.character(NodeStates[i,2])
    if(from == to){
      tmp <- phy$edge.length[i]
      names(tmp) <- from
      Map[[i]] <- tmp
    }else{
      shift.time <- shift.point * phy$edge.length[i]
      tmp <- c(shift.time, phy$edge.length[i] - shift.time)
      names(tmp) <- c(from, to)
      Map[[i]] <- tmp
    }
  }
  return(Map)
}

## transforms the phylogeny based on a set of paramaters and a simmap
transformPhy <- function(phy, map, pars){
  # phy must be of class simmap
  nTip <- length(phy$tip.label)
  RootAge <- max(branching.times(phy))
  NodeAges <- branching.times(phy)[phy$edge[,1] - nTip]
  ModMap <- Map <- map
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

## the likelihood function
getOULik <- function(phy, y, X, pars){
  # transform the phylogeny based on params
  tre <- transformPhy(phy, pars)
  # use the transformed phylogeny for the three point algorithm
  comp <- three.point.compute(tre$tree, y, X, tre$diag)
  # calculate the likelihood
  lik <- -as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
  return(lik)
}



# run/test node based
data(tworegime)
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
plot(tree)
nodelabels(pch=21, bg=select.reg)
trait[1:5,]
phy <- tree
tipstates <- trait[,2]
nodestates <- phy$node.label
shift.point <- 0.5

map <- getMapFromNode(phy, tipstates, nodestates, 0.5)
transformPhy(phy, map, pars)


# run/test simmap
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



## BASIC USAGE
require(phylolm)
require(OUwie)

# sim trait data
y = rTrait(n=1,MK$phy)
trait <- cbind(MK$data.legend, y[as.character(MK$data.legend[,1])])
PhyloLMDat <- trait[,3]
names(PhyloLMDat) <- trait[,1]

# params
sig <- c(0.32,0.07,0.02)
alp <- c(0.05,0.01,0.02)
opt <- c(1, -3, 10)
pars <- matrix(c(opt,sig,alp), 3, 3, dimnames = list(c(1,2,3), c("opt", "sig", "alp")))

# OUwie result
OUwieFixed <- OUwie.fixed(simmap[[1]],trait,model=c("OUMVA"), simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alp,sigma.sq=sig,theta=opt)

# use OUwie weights to calculate expectation
X <- colSums(t(OUwieFixed$regime.weights[-(1:3),1:3]) * opt)
# compare liks
getOULik(simmap[[1]], y, X, pars)
OUwieFixed$loglik





#### comparing to phylolm
# we need our alpha to be half the sigma to scale the branches in the same way
pars <- matrix(0, 3, 3, dimnames = list(c(1,2,3), c("opt", "sig", "alp")))
pars[,2] <- c(0.02,0.02,0.02)
pars[,3] <- c(0.01,0.01,0.01)

y = rTrait(n=1,MK$phy)
Q = y
Q[1:length(y)] = rep(0, length(y))

tre1 = transformPhy(simmap[[1]], pars)
tre2 = transf.branch.lengths(phy=simmap[[1]], model="OUfixedRoot", parameters = list(alpha=0.01))
three.point.compute(tre1$tree, y, Q, tre1$diag)
three.point.compute(tre2$tree, y, Q, tre2$diagWeight)


# comparing phylolm, OUwie, and finding the likelihood
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
opt <- c(1, -3, 10)
pars <- matrix(c(opt,sig,alp), 3, 3, dimnames = list(c(1,2,3), c("opt", "sig", "alp")))
nTip <- length(simmap[[1]]$tip.label)

OUwieFixed <- OUwie.fixed(simmap[[1]],trait,model=c("OUMVA"), simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alp,sigma.sq=sig,theta=opt)

tre1 <- transformPhy(simmap[[1]], pars)
#X <- getOUExpectation(transformPhy(phy, pars)$tree, pars)
#X <- getOUExpectation(simmap[[1]], pars)
X <- colSums(t(OUwieFixed$regime.weights[-(1:3),1:3]) * opt)
comp <- three.point.compute(tre1$tree, PhyloLMDat, X, tre1$diag)

OUwie.fixed(simmap[[1]],trait,model=c("OUMVA"), simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alp,sigma.sq=sig,theta=opt)
-as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2


# testing OUwie fixed
require(OUwie)
data(tworegime)
alpha=c(0.5632459,0.1726052)
sigma.sq=c(0.1064417,0.3461386)
theta=c(1.678196,0.4185894)

OUwieFixedInvert <- OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, shift.point=0.123, algorithm="invert")
#undebug(OUwie.fixed)
OUwieFixed3Point <- OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, shift.point=0.123, algorithm="three.point")
OUwieFixedInvert$loglik
OUwieFixed3Point$loglik

comparison <- identical(round(as.numeric(OUwieFixedInvert$loglik),5), round(as.numeric(OUwieFixed3Point$loglik),5))
comparison



# testing OUwie fixed OU1
require(OUwie)
data(tworegime)
alpha=c(0.5632459, 0.5632459)
sigma.sq=c(0.1064417, 0.1064417)
theta=c(1.678196, 1.678196)

OUwieFixedInvert <- OUwie.fixed(tree,trait,model=c("OU1"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, shift.point=0.123, algorithm="invert")
#undebug(OUwie.fixed)
OUwieFixed3Point <- OUwie.fixed(tree,trait,model=c("OU1"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, shift.point=0.123, algorithm="three.point")
OUwieFixedInvert$loglik
OUwieFixed3Point$loglik

comparison <- identical(round(as.numeric(OUwieFixedInvert$loglik),5), round(as.numeric(OUwieFixed3Point$loglik),5))
comparison


# testing the speed
require(OUwie)
require(corHMM)
require(geiger)

TestTime <- function(nTaxa){
  
  phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTaxa)
  phy <- drop.extinct(phy)
  q <- list(rbind(c(-5, 5), c(5, -5))/sum(phy$edge.length))
  d <- sim.char(phy, q, model = "discrete", n = 1)[,,1]
  c <- sim.char(phy, 5/sum(phy$edge.length), model = "BM", 1)[,,1]
  recon <- apply(round(ancRECON(phy, data.frame(sp = names(d), d = d), c(q[[1]][c(2,3)]), "marginal", 1)$lik.anc.states), 1, function(x) which(x == 1))
  phy$node.label <- recon
  
  tree <- phy
  trait <- data.frame(sp = names(d), reg = d, trt = c)
  alpha=c(sd(d),sd(d)/2)
  sigma.sq=c(var(d)/2, var(d))
  theta=c(mean(d) - sd(d), mean(d) + sd(d))
  
  InvertTime <- system.time(OUwieFixedInvert <- OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, shift.point=0.123, algorithm="invert"))
  #undebug(OUwie.fixed)
  ThrPntTime <- system.time(OUwieFixed3Point <- OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, shift.point=0.123, algorithm="three.point"))
  comparison <- identical(round(as.numeric(OUwieFixedInvert$loglik),5), round(as.numeric(OUwieFixed3Point$loglik),5))
  
  res <- c(InvertTime[3], ThrPntTime[3], comparison)
  names(res) <- c("InvertTime", "ThreeTime", "Match")
  return(res)
}


Times <- seq(from = 10, to = 2010, by = 50)
TimeTable <- matrix(0, length(Times), 3)
for(i in 1:length(Times)){
  TimeTable[i,] <- TestTime(Times[i])
}


plot(x = Times, y = TimeTable[,1], ylim = c(0, max(TimeTable)), pch = 16, col = "orange", xlab = "# of Tips", ylab = "Elapsed Time (s)")
points(x = Times, y = TimeTable[,2], pch = 16, col = "purple")
legend("topleft", legend = c("Invert", "3Pt"), pch = 16, col = c("orange", "purple"))

InvertTime
ThrPntTime



comparison <- identical(round(as.numeric(OUwieFixedInvert$loglik),5), round(as.numeric(OUwieFixed3Point$loglik),5))
comparison






require(OUwie)
data(tworegime)
set.seed(42)
alpha=c(0.358939, 0.3589399)
sigma.sq=c(0.5197486, 0.5197486)
theta=c( 1.3301447, 1.3301447)
OU1Invert <- OUwie.fixed(tree, trait, model=c("OU1"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="invert")

####FAILS IF RUNNING AS OU1:
OU13Point <- OUwie.fixed(tree, trait, model=c("OU1"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="three.point")

comparison <- identical(round(as.numeric(OU1Invert$loglik),5), round(as.numeric(OU13Point$loglik),5))
comparison




require(OUwie)
data(tworegime)
set.seed(42)
sigma.sq=c(0.5197486, 0.5197486*2)
#theta=c( 1.3264832, 1.3264832)

OU1Invert <- OUwie.fixed(tree, trait, model=c("BMS"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, sigma.sq=sigma.sq, shift.point=0.5, algorithm="invert")
theta=c(OU1Invert$theta[1,1], OU1Invert$theta[1,1])
####FAILS IF RUNNING AS OU1:
OU13Point <- OUwie.fixed(tree, trait, model=c("BMS"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, sigma.sq=sigma.sq, theta=theta, shift.point=0.5, algorithm="three.point")

OU1Invert$loglik
OU13Point$loglik

comparison <- identical(round(as.numeric(OU1Invert$loglik),5), round(as.numeric(OU13Point$loglik),5))
comparison

