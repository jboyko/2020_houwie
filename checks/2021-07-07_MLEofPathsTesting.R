setwd("~/2020_hOUwie/")

source("hOUwieSimmap.R")
source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)
require(partitions)


#### #### #### #### #### #### #### #### #### #### #### #### 
# functions
#### #### #### #### #### #### #### #### #### #### #### #### 
optMap <- function(p, phy, dat, nMap, type = "conditional"){
  model <- getRateCatMat(2)
  model[model > 0] <- p
  diag(model) <- -rowSums(model)
  simmap <- makeSimmap(tree=phy, data=dat, model=model, rate.cat=1, nSim=nMap, nCores=1)
  maps <- lapply(simmap, function(x) x$maps)
  llik <- max(unlist(lapply(maps, function(x) getMapProbability(x, model))))
  if(type == "joint"){
    MK <- corHMM(phy, dat, 1, model = "ER", p = p)
    llik <- llik + MK$loglik
  }
  cat("\r", llik)
  return(-llik)
}

getLiksByMapping <- function(phy, dat, nMap, type = "conditional"){
  opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
  out <- nloptr(x0 = 0.5, eval_f = optMap, lb = 1e-5, ub = 0.9, opts = opts, phy = phy, nMap = nMap, dat = dat, type = type)
  return(out)
}
#### #### #### #### #### #### #### #### #### #### #### #### 
# prerequisites
#### #### #### #### #### #### #### #### #### #### #### #### 
nCores <- 40
nTip <- 100
minRate = 0.25
maxRate = 0.75

#### #### #### #### #### #### #### #### #### #### #### #### 
# the phylogeny
#### #### #### #### #### #### #### #### #### #### #### #### 
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

#### #### #### #### #### #### #### #### #### #### #### #### 
# do the things
#### #### #### #### #### #### #### #### #### #### #### #### 
maxIter <- 80

singleRun <- function(phy, nMaps, minRate, maxRate){
  rate <- runif(1, min = minRate, max = maxRate)
  Q <- matrix(c(-rate, rate,rate,-rate),2,2)
  data <- simCharacterHistory(phy, Q, c(0.5,0.5))$TipStates
  dat <- data.frame(sp = names(data), d = data)
  corHMM_res <- corHMM(phy, dat, 1, model = "ER")
  Mapping_res <- list()
  for(i in 1:length(nMaps)){
    cat("\n\nWorking on nMap", i, "of length", nMaps[i], "\n")
    Mapping_res[[i]] <- getLiksByMapping(phy, dat, nMaps[i])
  }
  fit_summary <- data.frame(true_rate = rate, corhmm_rate = corHMM_res$solution[1,2], mapping_rate = unlist(lapply(Mapping_res, function(x) x$solution)), nMaps = nMaps)
  cat("\n")
  return(list(fit_summary = fit_summary,
              data = dat,
              corHMM_res = corHMM_res,
              Mapping_res = Mapping_res))
}

result <- mclapply(1:40, function(x) singleRun(phy, c(10,50,100,500), minRate, maxRate), mc.cores = 40)




require(corHMM)
require(nloptr)

##run corhmm
MK <- corHMM(phy, dat, 1, model = "ER")
model <- MK$solution
# optimize the map, only a few maps for speed
out <- nloptr(x0 = 0.1, eval_f = optMap, lb = 1e-5, ub = 0.5, opts = opts, phy = phy, nMap = 10, dat = dat)

# the results
print(c(Mapping_rate = out$solution, corHMM_rate = MK$solution[2,1]))

QA <- MK$solution
diag(QA) <- 0
diag(QA) <- -rowSums(QA)

QB <- QA
QB[QB > 0] <- out$solution

# what if we use the corHMM optimization for the likelihood
simmapA <- makeSimmap(tree=phy, data=dat, model=QA, rate.cat=1, nSim=1000, nCores=1)
MapProbsA <- unlist(lapply(simmapA, function(x) getMapProbability(x$maps, QA)))
# what if we use the mapping optimization for the likelihood
simmapB <- makeSimmap(tree=phy, data=dat, model=QB, rate.cat=1, nSim=1000, nCores=1)
MapProbsB <- unlist(lapply(simmapB, function(x) getMapProbability(x$maps, Q)))

print(c(Mapping_rate_lik = max(MapProbsB), corHMM_rate_lik = max(MapProbsA)))









