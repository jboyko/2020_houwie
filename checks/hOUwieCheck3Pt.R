# simulation tests for optimizing the optimization of hOUwie 
simulateData <- function(phy, q, alpha, sigma.sq, quiet=FALSE){
  # phy$edge.length <- phy$edge.length/max(branching.times(phy)) # following Cresler et al. (2015) we set T to be 1
  Q = matrix(c(-q,q,q,-q), 2, 2)
  
  root.p = c(0,0) # we will sample the root with equal probability of either state
  root.p[sample(c(1,2), 1, prob =c(0.5, 0.5))] <- 1
  
  theta = c(5, 10) # following Cresler et al. (2015) we set delta theta to be 1 (tempor at 2)
  theta0 = rnorm(1, theta[which(root.p == 1)], sqrt(sigma.sq[which(root.p == 1)]/2*alpha[which(root.p == 1)])) # following Cresler et al. (2015) we sample the root theta from the stationary distribution matchiing the root state
  full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
  obs.no.trans <- sum(unlist(lapply(full.data$simmap[[1]]$maps, function(x) length(x) - 1)))
  if(!quiet){
    cat("The observed number of transitions was found to be", obs.no.trans, "\n")
  }
  return(full.data)
}

nLikEval <- function(p.subset, nEval, nSimmaps, tree, trait){
  p <- c(mk = NA, alpha = NA, sig2 = NA, opt1 = 4, opt2 = 5)
  p[1] <- p.subset[1]
  p[2] <- p.subset[2]
  p[3] <- p.subset[3]
  lliks <- sapply(1:nEval, function(x) hOUwie(phy = tree, data = trait, rate.cat = 1, model.cor = "ER", model.ou = "OUM", p = p, nSim = nSimmaps, quiet = FALSE)$loglik)
  return(lliks)
}

# get a phylogeny rescaled to a height of 1
source("~/2020_hOUwie/hOUwieSimmap.R")
#source("hOUwie.R")
require(OUwie)
require(corHMM)
require(parallel)
require(lhs)
require(expm)

# simulate data
nTip <- 100
q <- 0.25 # the expected number of markov transitions in terms of br proportion
alpha <- c(1,1) 
sigma.sq <- c(0.5,0.5) 
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) # this will be scaled to H=1
phy <- drop.extinct(phy)
Tmax <- max(branching.times(phy))
print(Tmax)
# phy$edge.length <- phy$edge.length/max(branching.times(phy))
data.houwie <- simulateData(phy, q, alpha, sigma.sq)[[1]]

# set up a range of params to explore based on latin hypercube
n.par.combos <- 100 # the resolution of the heat map
RNGkind("L'Ecuyer-CMRG")
set.seed(1985)
X <- randomLHS(n = n.par.combos, k = 3)
X <- X * 1
X[,1] <- X[,1] * 2
colnames(X) <- c("p.mk", "alpha", "sig2")
# pairs(X)
par.set <- vector("list", dim(X)[1])
for(i in 1:dim(X)[1]){
  par.set[[i]] <- X[i,]
}

# evaluate 
n.evaluations <- 1 # the accuracy of our mean and sd estimates
n.simmaps <- 10

out <- list()
for(i in 1:n.par.combos){
  out[[i]] <- nLikEval(par.set[[i]], n.evaluations, n.simmaps, phy, data.houwie)
}

debug(OUwie.basic)

hOUwie(phy = phy, data = data.houwie, rate.cat = 1, model.cor = "ER", model.ou = "OUM", p = c(1, 8, 1, 5, 10), nSim = n.simmaps, quiet = FALSE)

# test <- hOUwie(phy = phy, data = data.houwie, rate.cat = 1, 
#                index.cor = index.cor, root.p="yang", lb.cor=1e-5, ub.cor=1,
#                index.ou = index.ou, root.station=FALSE, get.root.theta=FALSE, lb.ou=c(1e-5,1e-5,1e-1), ub.ou=c(2, 10, 10),
#                nSim = 10, weighted = FALSE)




out.sd <- unlist(lapply(out, function(x) sd( x[!x == -Inf])))
out.mu <- unlist(lapply(out, function(x) mean( x[!x == -Inf])))
true.par <- c("p.mk" = q, "alpha" = alpha[1], "sig2" = sigma.sq[1])

res <- do.call(rbind, par.set)
res <- cbind(res, out.mu, out.sd)
colnames(res)[c(4,5)] <- c("mean", "sd")
res <- res[apply(res, 1, function(x) !any(is.na(x))),]
res <- as.data.frame(res)

