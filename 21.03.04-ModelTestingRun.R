# a function that given a model structure and params, generates a hOUwie dataset.
generateData <- function(phy, index.cor, index.ou, pars, quiet=FALSE){
  phy$edge.length <- phy$edge.length/max(branching.times(phy)) # ensure tree height is 1
  k.cor <- max(index.cor, na.rm = TRUE) # number of corhmm params
  k.ou <- max(index.ou, na.rm = TRUE) # number of ouwie params
  
  if((k.cor + k.ou) != length(pars)){
    stop("Length of pars does not match index matrices.", call. = FALSE)
  }
  
  p.mk <- pars[1:k.cor]
  p.ou <- pars[(k.cor+1):length(pars)]
  
  # organize corhmm params
  index.cor[index.cor == 0] <- NA
  index.cor[is.na(index.cor)] <- max(index.cor, na.rm = TRUE) + 1
  Q <- matrix(0, dim(index.cor)[1], dim(index.cor)[1])
  Q[] <- c(p.mk, 0)[index.cor]
  diag(Q) <- -rowSums(Q)
  root.p = rep(0, dim(Q)[1]) # we will sample the root with equal probability of either state
  root.p[sample(1:dim(Q)[1], 1)] <- 1
  
  # organize ou params
  Rate.mat <- matrix(1, 3, dim(index.cor)[2])
  index.ou[is.na(index.ou)] <- max(index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.ou]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  theta0 = rnorm(1, theta[which(root.p == 1)], sqrt(sigma.sq[which(root.p == 1)]/2*alpha[which(root.p == 1)])) # following Beaulieu et al (2012) we sample the root theta from the stationary distribution matchiing the root state
  full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
  obs.no.trans <- sum(unlist(lapply(full.data$simmap[[1]]$maps, function(x) length(x) - 1)))
  if(!quiet){
    cat("The observed number of transitions was found to be", obs.no.trans, "\n")
  }
  simulators <- list(index.cor = index.cor, index.ou = index.ou, pars = pars)
  return(c(full.data, simulators))
}

SingleModelTestRun <- function(phy, index.cor, index.ou, model.cor, model.ou, pars, nmaps, iter){
  cat("Begining single run...\n")
  data.houwie <- generateData(phy, index.cor, index.ou, pars)
  cat("\nBegining two step...\n")
  TwoStepFit <- fitTwoStep(phy=phy, data=data.houwie$data, rate.cat=1, 
                           model.cor=model.cor, root.p="yang", lb.cor=1e-3, ub.cor=10,
                           model.ou=model.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none")
  cat("\nBegining non cens...\n")
  NonCensFit <- fitNonCensored(phy=phy, data=data.houwie$data, rate.cat=1, 
                               model.cor=model.cor, root.p="yang", lb.cor=1e-3, ub.cor=10,
                               model.ou=model.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none",
                               nSim=nmaps, nCores=1)
  cat("\nBegining houwie...\n")
  hOUwieFit <- hOUwie(phy=phy, data=data.houwie$data, rate.cat=1,
                      model.cor=model.cor, root.p="yang", lb.cor=1e-3, ub.cor=10,
                      model.ou=model.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none",
                      nSim=nmaps, nCores=1)
  
  cat("\nBegining truth\n")
  TruMapFit <- OUwie(phy=data.houwie$simmap[[1]], data=data.houwie$data,
                    model=model.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none",
                    simmap.tree=TRUE, scaleHeight=FALSE, algorithm="three.point", diagn=FALSE)
  
  obj <- list(data.houwie = data.houwie,
              TwoStepFit = TwoStepFit,
              NonCensFit  = NonCensFit,
              hOUwieFit = hOUwieFit,
              TruMapFit = TruMapFit)
  cat("Done.\n")
  
  file.name <- paste0("Fit=", model.cor, "_", model.ou, "-", "Pars=", paste(pars, collapse = "_"), "-Iter=", iter, ".Rsave")
  print(file.name)
  save(obj, file = file.name)
  return(obj)
}

# get a phylogeny rescaled to a height of 1
source("~/2020_hOUwie/hOUwie.R")
#source("../../../hOUwie.R")
require(OUwie)
require(corHMM)
require(parallel)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ER BMS
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) # this will be scaled to H=1
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

index.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
index.ou <- getOUParamStructure("BMS", "three.point", TRUE, TRUE, dim(index.cor)[1])
nmaps <- 100

# test <- SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", pars, nmaps, 1)
# test$TwoStepFit
# test$NonCensFit
# test$hOUwieFit

# run within the folder of the fitting model and a subfolder of simulating params
pars <- c(1, 1, 10, 5, 5)
mclapply(101:500, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", pars, nmaps, x), mc.cores = 15)
pars <- c(2, 1, 10, 5, 5)
mclapply(101:500, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", pars, nmaps, x), mc.cores = 15)
pars <- c(4, 1, 10, 5, 5)
mclapply(101:500, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", pars, nmaps, x), mc.cores = 15)
pars <- c(8, 1, 10, 5, 5)
mclapply(101:500, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", pars, nmaps, x), mc.cores = 15)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ER OUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) # this will be scaled to H=1
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

index.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
index.ou <- getOUParamStructure("OUM", "three.point", TRUE, TRUE, dim(index.cor)[1])
pars <- c(4, 2, 0.5, 5, 10)
nmaps <- 100

# test <- SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMV", pars, nmaps, 1)

# for server
pars <- c(1, 2, 0.5, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", pars, nmaps, x), mc.cores = 15)
pars <- c(2, 2, 0.5, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", pars, nmaps, x), mc.cores = 15)
pars <- c(4, 2, 0.5, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", pars, nmaps, x), mc.cores = 15)
pars <- c(8, 2, 0.5, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", pars, nmaps, x), mc.cores = 15)




## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ER OUMV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) # this will be scaled to H=1
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

index.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
index.ou <- getOUParamStructure("OUMV", "three.point", TRUE, TRUE, dim(index.cor)[1])
pars <- c(1, 2, 0.5, 0.5, 5, 10)
nmaps <- 100

# test <- SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMV", pars, nmaps, 1)

# for server
pars <- c(1, 2, 0.5, 0.5, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMV", pars, nmaps, x), mc.cores = 15)
pars <- c(1, 2, 0.5, 1, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMV", pars, nmaps, x), mc.cores = 15)
pars <- c(1, 2, 0.5, 2, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMV", pars, nmaps, x), mc.cores = 10)
pars <- c(1, 2, 0.5, 5, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMV", pars, nmaps, x), mc.cores = 10)

