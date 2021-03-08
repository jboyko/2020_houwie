# simulation tests for optimizing the optimization of hOUwie 
simulateData <- function(phy, q, alpha, sigma.sq, quiet=FALSE){
  phy$edge.length <- phy$edge.length/max(branching.times(phy)) # following Cresler et al. (2015) we set T to be 1
  Q = matrix(c(-q,q,q,-q), 2, 2)
  
  root.p = c(0,0) # we will sample the root with equal probability of either state
  root.p[sample(c(1,2), 1, prob =c(0.5, 0.5))] <- 1
  
  theta = c(4, 5) # following Cresler et al. (2015) we set delta theta to be 1 (tempor at 2)
  theta0 = rnorm(1, theta[which(root.p == 1)], sqrt(sigma.sq[which(root.p == 1)]/2*alpha[which(root.p == 1)])) # following Beaulieu et al (2012) we sample the root theta from the stationary distribution matchiing the root state
  full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
  obs.no.trans <- sum(unlist(lapply(full.data$simmap[[1]]$maps, function(x) length(x) - 1)))
  if(!quiet){
    cat("The observed number of transitions was found to be", obs.no.trans, "\n")
  }
  return(full.data)
}

# get a phylogeny rescaled to a height of 1
source("~/2020_hOUwie/hOUwie.R")
#source("../hOUwie.R")
require(OUwie)
require(corHMM)
require(parallel)

# simulate data
nTip <- 250
q <- 2.5 # the expected number of markov transitions in terms of br proportion
alpha <- c(0.5,0.5) 
sigma.sq <- c(0.5,0.5) 
set.seed(1985)
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) # this will be scaled to H=1
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
data.houwie <- simulateData(phy, q, alpha, sigma.sq)[[1]]
model.cor="ER"
model.ou="OUM"
index.ou=NULL
rate.cat=1
root.station=FALSE
get.root.theta=FALSE
# evaluate 
n.evaluations <- 500 # the accuracy of our mean and sd estimates

hOUwie.dat <- organizeHOUwieDat(data.houwie, "none")
nObs <- length(hOUwie.dat$ObservedTraits)
#reorder phy
phy <- reorder(phy, "pruningwise")
# a way to speed up the three.point function
tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
Tmax <- max(branching.times(phy))

null.cor <- FALSE
if(is.null(model.cor)){
  model.cor <- "ER"
  null.cor <- TRUE
}
model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
index.ou <- getOUParamStructure(model.ou, "three.point", root.station, get.root.theta, dim(model.set.final$Q)[1])


est.pars <- c(2, 0.5, 0.5, 4, 5)
out <- hOUwie.dev.tmp(log(est.pars), 
                      phy=phy, rate.cat=rate.cat,
                      data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p="yang", 
                      data.ou=hOUwie.dat$data.ou, index.ou=index.ou, 
                      algorithm="three.point", mserr="none",
                      nSim=1000, nCores=1, tip.paths=tip.paths, weighted=FALSE)

hist(out$lik.vec)


par(mfrow=c(2,2))
d <- sample(out$lik.vec, size = 1000)
true.mean <- mean(d)
hist(d, main = "1000/1000", xlab = paste0("sd = ", round(sd(d), 3)))
abline(v = true.mean, col = "darkgreen", lwd = 2)
d <- sample(out$lik.vec, size = 900)
hist(d, main = "900/1000", xlab = paste0("sd = ", round(sd(d), 3)))
abline(v = true.mean, col = "darkgreen", lwd = 2)
d <- sample(out$lik.vec, size = 500)
hist(d, main = "500/1000", xlab = paste0("sd = ", round(sd(d), 3)))
abline(v = true.mean, col = "darkgreen", lwd = 2)
d <- sample(out$lik.vec, size = 100)
hist(d, main = "100/1000", xlab = paste0("sd = ", round(sd(d), 3)))
abline(v = true.mean, col = "darkgreen", lwd = 2)


