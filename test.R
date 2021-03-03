## functions
simulateData <- function(phy, no.of.trans, eta, gamma, quiet=FALSE){
  phy$edge.length <- phy$edge.length/max(branching.times(phy)) # following Cresler et al. (2015) we set T to be 1
  
  q <- no.of.trans/sum(phy$edge.length) # the expected number of trans gives our Mk rate
  Q = matrix(c(-q,q,q,-q), 2, 2)
  
  root.p = c(0,0) # we will sample the root with equal probability of either state
  root.p[sample(c(1,2), 1, prob =c(0.5, 0.5))] <- 1
  
  alpha = eta # alpha is eta/T and T is fixed at one
  sigma.sq = gamma # sigma sq is gamma * 1 / 1
  theta = c(4, 5) # following Cresler et al. (2015) we set delta theta to be 1
  theta0 = rnorm(1, theta[which(root.p == 1)], sqrt(sigma.sq[which(root.p == 1)]/2*alpha[which(root.p == 1)])) # following Cresler et al. (2015) we sample the root theta from the stationary distribution matchiing the root state
  full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
  obs.no.trans <- sum(unlist(lapply(full.data$simmap[[1]]$maps, function(x) length(x) - 1)))
  if(!quiet){
    cat("The observed number of transitions was found to be", obs.no.trans, "\n")
  }
  return(full.data)
}


# testing
#
require(corHMM)
require(OUwie)
require(parallel)
require(geiger)
source("~/2020_hOUwie/hOUwie.R")

nTip <- 100

prop.trans <- 0.1 # the expected number of markov transitions in terms of br proportion
eta <- c(1,1) # a dimensionless measure of selection opportunity (alpha*T)
gamma <- c(1,1) # a dimensionless measure noise intensity (sigma*T^(1/2)/deltaTheta)

phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) # this will be scaled to H=1
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
no.of.trans <- prop.trans * dim(phy$edge)[1]

data <- simulateData(phy, no.of.trans = no.of.trans, eta = eta, gamma = gamma)[[1]]

model.ou <- "OUM"
model.cor <- "ER"

# index.cor <- corHMM:::getFullMat(list(corHMM:::getStateMat4Dat(data[,c(1,2)], "ER")$rate.mat, corHMM:::getStateMat4Dat(data[,c(1,2)], "ER")$rate.mat), corHMM:::getStateMat4Dat(data[,c(1,2)], "ER")$rate.mat)
index.cor <- corHMM:::getStateMat4Dat(data[,c(1,2)], model.cor)$rate.mat
index.ou <- getOUParamStructure(model.ou, "three.point", FALSE, FALSE, dim(index.cor)[2])
p = c(1.5, 1, 1, 3, 4)

test <- hOUwie(phy = phy, data = data, rate.cat = 1, 
               index.cor = index.cor, root.p="yang", lb.cor=1e-5, ub.cor=1,
               index.ou = index.ou, root.station=FALSE, get.root.theta=FALSE, lb.ou=c(1e-5,1e-5,1e-1), ub.ou=c(2, 10, 10),
               nSim = 10, weighted = FALSE, opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="10", "ftol_rel"=.Machine$double.eps^0.5), p = p)

test2 <- fitNonCensored(phy = phy, data = data, rate.cat = 1, 
               model.cor = model.cor, root.p="yang", lb.cor=1e-5, ub.cor=1,
               model.ou = model.ou, root.station=FALSE, get.root.theta=FALSE, nSim = 10)


test
test2



debug(hOUwie)
debug(fitNonCensored)

source("~/2020_hOUwie/hOUwie.R")
test
# 