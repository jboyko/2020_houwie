require(expm)
#### #### #### #### #### #### #### #### #### #### #### #### 
# functions
#### #### #### #### #### #### #### #### #### #### #### #### 

getBranchProb <- function(Q, bl, x0, xt, n){
  if(n == 0){
    path_i <- c(x0, xt)
    llik_vec <- getPathProb(path_i, bl, Q)
    names(llik_vec) <- paste0(path_i, collapse="_")
    return(llik_vec)
  }
  t <- bl/(n+1)
  p_mat <- expm(Q * t)
  combo.list <- list()
  for(i in 1:n){
    combo.list[[i]] <- c(1,2)
  }
  combos <- expand.grid(combo.list)
  llik_vec <- numeric(dim(combos)[1])
  for(i in 1:dim(combos)[1]){
    path_i <- c(x0, as.numeric(combos[i,]), xt)
    llik_vec[i] <- getPathProb(path_i, bl, Q)
    names(llik_vec)[i] <- paste0(path_i, collapse="_")
  }
  return(llik_vec)
}
# 
getPathProb <- function(path, bl, Q){
  section_length <- bl/(length(path)-1)
  p_mat <- expm(Q * section_length)
  path_cp <- path
  P <- vector("numeric", length(path)-1)
  for(i in 1:(length(path)-1)){
    P[i] <- p_mat[path_cp[1],path_cp[2]]
    path_cp <- path_cp[-1]
  }
  return(sum(log(P)))
}

#### #### #### #### #### #### #### #### #### #### #### #### 
# evaluation of branhces
#### #### #### #### #### #### #### #### #### #### #### #### 

Q <- matrix(c(-1,1,1,-1), 2, 2)
x0 = 1
xt = 1
n = 1
bl = 1

c(getPathProb(c(1,1), bl, Q))

c(getPathProb(c(1,1,1), bl, Q),
  getPathProb(c(1,2,1), bl, Q))

c(getPathProb(c(1,1,1,1), bl, Q),
  getPathProb(c(1,2,1,1), bl, Q),
  getPathProb(c(1,1,2,1), bl, Q),
  getPathProb(c(1,2,2,1), bl, Q))

getBranchProb(Q, bl, x0 = 1, xt = 2, n=2)

sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=0)))
sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=1)))
sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=2)))
sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=3)))

#### #### #### #### #### #### #### #### #### #### #### #### 
### what if we are talking about a node with 2 descendents rather than a single branch?
#### #### #### #### #### #### #### #### #### #### #### #### 

# Node = 0, Dec_1 = 0, Dec_2 = 1
c(expm(Q) %*% c(1,0))[1] * c(expm(Q) %*% c(0,1))[1]
# path prob
exp(getPathProb(c(1,1), bl, Q) + getPathProb(c(1,2), bl, Q))
# branch prob no change
exp(getBranchProb(Q, 1, x0 = 1, xt = 1, n=0) + getBranchProb(Q, 1, x0 = 1, xt = 2, n=0))
# branch prob 2 changes
# dec_1
D1 <- getBranchProb(Q, 1, x0 = 1, xt = 1, n=2)
# dec_2
D2 <- getBranchProb(Q, 1, x0 = 1, xt = 2, n=2)
combos <- expand.grid(1:4, 1:4)
P <- numeric(dim(combos)[1])
for(i in 1:dim(combos)[1]){
  combo_i <- as.numeric(combos[i,])
  P[i] <- D1[combo_i[1]] + D2[combo_i[2]]
}
sum(exp(P))

#### #### #### #### #### #### #### #### #### #### #### #### 
## let's expand out to an entire tree now
#### #### #### #### #### #### #### #### #### #### #### #### 
setwd("~/2020_houwie/")
require(geiger)
require(corHMM)
require(partitions)
require(expm)
require(MASS)
source("hOUwieNode.R")

nTip <- 5
rate <- 1
n <- 200
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
Tmax <- max(branching.times(phy))
Q <- matrix(c(-rate, rate,rate,-rate),2,2)
nState <- dim(Q)[1]
data <- sim.char(phy, Q, 1, "discrete")[,,1]
data <- data.frame(sp = names(data), d = data)

phy <- reorder.phylo(phy, "pruningwise")
edge_liks_list <- getEdgeLiks(phy, data, 2, 1, 0.25)

# edge_liks_list[[1]][2,] <- c(0,1)
# edge_liks_list[[7]][4,] <- c(0,1) 
# edge_liks_list[[8]][4,] <- c(0,1)

conditional_probs <- getConditionalInternodeLik(phy, Q, edge_liks_list, "yang")
internode_maps <- getInternodeMap(phy, Q, conditional_probs$edge_liks_list, conditional_probs$root_state, 1)
simmaps <- getMapFromSubstHistory(internode_maps, phy)
times_per_edge <- unlist(lapply(internode_maps[[1]], function(x) x[1]))
p_mats_per_edge <- lapply(times_per_edge, function(x) expm(Q * x))
llik_discrete_A <- unlist(lapply(simmaps, function(x) getMapProb(x, Q, c(0.5,0.5), p_mats_per_edge)))
# llik_discrete_B <- unlist(lapply(simmaps, function(x) getMapProb(x, Q, c(0.5,0.5), NULL)))
ouwie_res <- lapply(simmaps, function(x) OUwie.basic(x,cbind(data, rnorm(nTip, 10)), simmap.tree=TRUE, alpha = c(1,1), sigma.sq = c(1,1), theta = c(10,10)))

log(sum(exp(llik_discrete_A)))
corHMM(phy, data, 1, model = "ER", p = rate)



#### #### #### #### #### #### #### #### #### #### #### #### 
## let's test the hOUwie functions now
#### #### #### #### #### #### #### #### #### #### #### #### 
setwd("~/2020_houwie/")
require(geiger)
require(corHMM)
require(partitions)
require(expm)
require(MASS)
require(treeplyr)
require(OUwie)
require(parallel)

source("hOUwieNode.R")
source("Utils.R")

# import the data
scales <- read.csv("empirical/pnh-ms-master/R/data/scales2009.csv")
scales <- scales[!scales$species == "Carlia_fusca",] # for simplicity remove the only mixed forager
tetTree <- read.tree("empirical/pnh-ms-master/R/data/tetrapods.tre")
tetTree$tip.label[tetTree$tip.label == "Eumeces_schneideri"] <- "Eumeces_fasciatus"
tetTree <- multi2di(tetTree)
tetTree$edge.length[tetTree$edge.length==0] <- .Machine$double.eps
td <- make.treedata(tetTree, scales)
td <- reorder(td, "postorder")
phy <- td$phy
dat <- data.frame(sp = td$phy[['tip.label']], FM = td[['foraging.mode']], PE = td[['predator.escape']], FOG = td[['FG.frac']])

# generate the discrete model
getStateMat4Dat(dat[c(1:3)])$legend
StateDepA <- equateStateMatPars(getStateMat4Dat(dat[c(1:3)])$rate.mat, 1:12)
ParamDepA <- equateStateMatPars(getRateCatMat(2), 1:2)
DiscModelA <- getFullMat(list(StateDepA, StateDepA), ParamDepA)
DiscModelB <- equateStateMatPars(DiscModelA, 1:2)

# generate the continuous model
HYBBMS <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(DiscModelA)[1])
SWBMS <- CDBMS <- CIDBMS <- HYBBMS
CDBMS[2,] <- c(1:5,1:5)
CDBMS[3,] <- 6
CIDBMS[2,] <- c(1,1,1,1,1,2,2,2,2,2)
CIDBMS[3,] <- 3
# SWBMS[2,] <- c(1,2,1,3)
# SWBMS[3,] <- 4
HYBOUM <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(DiscModelA)[1])
PEOUM <- FMOUM <- SWOUM <- CDOUM <- CIDOUM <- HYBOUM
CDOUM[3,] <- c(3:7,3:7)
CIDOUM[3,] <- c(3,3,3,3,3,4,4,4,4,4)
FMOUM[3,] <- c(3,3,4,4,4,3,3,4,4,4)
PEOUM[3,] <- c(3,4,5,3,4,3,4,5,3,4)
SWOUM[3,] <- c(3,3,4,4,4,3,3,5,5,5)


fit.BM1_A <- hOUwie(phy, dat, 2, discrete_model = DiscModelB, continuous_model = "BM1", nSim = 10, root.p = c(rep(1,10)), p = c(1.218069e-03, 9.711943e-05, 4.450969e-04, 2.268529e+00))


