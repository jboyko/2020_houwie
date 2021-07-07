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
optMap <- function(p, phy, dat, nMap, type = "max"){
  model <- getRateCatMat(2)
  model[model > 0] <- p
  diag(model) <- -rowSums(model)
  simmap <- makeSimmap(tree=phy, data=dat, model=model, rate.cat=1, nSim=nMap, nCores=1)
  maps <- lapply(simmap, function(x) x$maps)
  llik <- unlist(lapply(maps, function(x) getMapProbability(x, model)))
  if(type == "max"){
    llik <- max(llik)
  }
  if(type == "sum"){
    llik <- max(llik) + log(sum(exp(llik - max(llik)))) - log(nMap)
  }
  cat("\r", llik)
  return(-llik)
}

getLiksByMapping <- function(phy, dat, nMap, type = "max"){
  opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
  out <- nloptr(x0 = 0.5, eval_f = optMap, lb = 1e-5, ub = 0.9, opts = opts, phy = phy, nMap = nMap, dat = dat, type = type)
  return(out)
}
#### #### #### #### #### #### #### #### #### #### #### #### 
# prerequisites
#### #### #### #### #### #### #### #### #### #### #### #### 
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
# a small aside looking at variance
#### #### #### #### #### #### #### #### #### #### #### #### 

# rate <- runif(1, min = minRate, max = maxRate)
# Q <- matrix(c(-rate, rate,rate,-rate),2,2)
# data <- simCharacterHistory(phy, Q, c(0.5,0.5))$TipStates
# dat <- data.frame(sp = names(data), d = data)
# simmap <- makeSimmap(tree=phy, data=dat, model=Q, rate.cat=1, nSim=10, nCores=1)
# maps <- lapply(simmap, function(x) x$maps)
# llik_10 <- unlist(lapply(maps, function(x) getMapProbability(x, Q)))
# simmap <- makeSimmap(tree=phy, data=dat, model=Q, rate.cat=1, nSim=100, nCores=1)
# maps <- lapply(simmap, function(x) x$maps)
# llik_100 <- unlist(lapply(maps, function(x) getMapProbability(x, Q)))
# simmap <- makeSimmap(tree=phy, data=dat, model=Q, rate.cat=1, nSim=500, nCores=1)
# maps <- lapply(simmap, function(x) x$maps)
# llik_500 <- unlist(lapply(maps, function(x) getMapProbability(x, Q)))
# 
# c(Map10_lik = max(llik_10), Map10_var = var(llik_10),
#   Map100_lik = max(llik_100), Map100_var = var(llik_100),
#   Map500_lik = max(llik_500), Map500_var = var(llik_500))

#### #### #### #### #### #### #### #### #### #### #### #### 
# do the things
#### #### #### #### #### #### #### #### #### #### #### #### 
singleRun <- function(phy, nMaps, minRate, maxRate){
  rate <- runif(1, min = minRate, max = maxRate)
  Q <- matrix(c(-rate, rate,rate,-rate),2,2)
  data <- simCharacterHistory(phy, Q, c(0.5,0.5))$TipStates
  dat <- data.frame(sp = names(data), d = data)
  corHMM_res <- corHMM(phy, dat, 1, model = "ER")
  Mapping_res_B <- Mapping_res_A <- list()
  for(i in 1:length(nMaps)){
    cat("\n\nWorking on nMapA", i, "of length", nMaps[i], "\n")
    Mapping_res_A[[i]] <- getLiksByMapping(phy, dat, nMaps[i], type = "max")
  }
  for(i in 1:length(nMaps)){
    cat("\n\nWorking on nMapB", i, "of length", nMaps[i], "\n")
    Mapping_res_B[[i]] <- getLiksByMapping(phy, dat, nMaps[i], type = "sum")
  }
  Mapping_res <- c(Mapping_res_A, Mapping_res_B)
  corhmm_rate <- corHMM_res$solution[1,2]
  mapping_rate <- unlist(lapply(Mapping_res, function(x) x$solution))
  estimated_rate <- c(corhmm_rate, mapping_rate)
  nStochMaps <- c(NA, nMaps, nMaps)
  type <- c("corhmm", rep("max", length(nMaps)), rep("sum", length(nMaps)))
  fit_summary <- data.frame(true_rate = rate, estimated_rate = estimated_rate, nMaps = nStochMaps, type = type)
  cat("\n")
  return(list(fit_summary = fit_summary,
              data = dat,
              corHMM_res = corHMM_res,
              Mapping_res = Mapping_res))
}

result <- mclapply(1:40, function(x) singleRun(phy, nMaps = c(10,50,100,500), minRate, maxRate), mc.cores = 40)







require(viridis)
load("discrete_liks_check.Rsave")
par(mfrow=c(1,2))
plot(x = discrete_liks_check$true_rate, y = discrete_liks_check$corhmm_rate, xlim = c(0,1), ylim = c(0,1), xlab = "true rate", ylab = "corhmm rate")
abline(a = 0, b = 1, col = "red")
abline(lm(discrete_liks_check$corhmm_rate ~ discrete_liks_check$true_rate))

plot(x = discrete_liks_check$true_rate, y = discrete_liks_check$mapping_rate, xlim = c(0,1), ylim = c(0,1), xlab = "true rate", ylab = "mapping rate")
abline(a = 0, b = 1, col = "red")
abline(lm(discrete_liks_check$mapping_rate ~ discrete_liks_check$true_rate))


abline(a = 0, b = 1)



