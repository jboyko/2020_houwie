source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")
# source("/space_2/jamesboyko/2020_hOUwie/hOUwie.R")
# source("/space_2/jamesboyko/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)

####### ####### ####### ####### ####### ####### ####### ####### ####### 
# models to fit 
####### ####### ####### ####### ####### ####### ####### ####### ####### 
# CD BM1
fit.CDBM1.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CDBM1.ou <- getOUParamStructure("BM1", "three.point", FALSE, FALSE, dim(fit.CDBM1.cor)[1])
# CD BMS
fit.CDBMS.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.CDBMS.cor)[1])
# CD OU1
fit.CDOU1.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CDOU1.ou <- getOUParamStructure("OU1", "three.point", FALSE, FALSE, dim(fit.CDOU1.cor)[1])
# CD OUM
fit.CDOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.CDOUM.cor)[1])
# fit CID BMS
fit.CIDBMS.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CIDBMS.cor <- getFullMat(list(fit.CIDBMS.cor, fit.CIDBMS.cor), fit.CIDBMS.cor)
fit.CIDBMS.cor <- equateStateMatPars(fit.CIDBMS.cor, c(1,2,3))
fit.CIDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.CIDBMS.cor)[1])
fit.CIDBMS.ou[2,] <- c(1,1,2,2)
fit.CIDBMS.ou[3,] <- 3
# fit CID OUM
fit.CIDOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CIDOUM.cor <- getFullMat(list(fit.CIDOUM.cor, fit.CIDOUM.cor), fit.CIDOUM.cor)
fit.CIDOUM.cor <- equateStateMatPars(fit.CIDOUM.cor, c(1,2,3))
fit.CIDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.CIDOUM.cor)[1])
fit.CIDOUM.ou[3,] <- c(3,3,4,4)
# fit Hyb BMS
fit.HYBBMS.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.HYBBMS.cor <- getFullMat(list(fit.HYBBMS.cor, fit.HYBBMS.cor), fit.HYBBMS.cor)
fit.HYBBMS.cor <- equateStateMatPars(fit.HYBBMS.cor, c(1,2,3))
fit.HYBBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.HYBBMS.cor)[1])
# fit Hyb OUM
fit.HYBOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.HYBOUM.cor <- getFullMat(list(fit.HYBOUM.cor, fit.HYBOUM.cor), fit.HYBOUM.cor)
fit.HYBOUM.cor <- equateStateMatPars(fit.HYBOUM.cor, c(1,2,3))
fit.HYBOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.HYBOUM.cor)[1])

####### ####### ####### ####### ####### ####### ####### ####### ####### 
# simulate under one and fit under all others
####### ####### ####### ####### ####### ####### ####### ####### ####### 
index.cor <- list(fit.CDBM1.cor, fit.CDBMS.cor, fit.CDOU1.cor, fit.CDOUM.cor, 
                  fit.CIDBMS.cor, fit.CIDOUM.cor, fit.HYBBMS.cor, fit.HYBOUM.cor)
index.ou <- list(fit.CDBM1.ou, fit.CDBMS.ou, fit.CDOU1.ou, fit.CDOUM.ou, 
                 fit.CIDBMS.ou, fit.CIDOUM.ou, fit.HYBBMS.ou, fit.HYBOUM.ou)
rate.cats <- c(CDBM1=1,CDBMS=1,CDOU1=1,CDOUM=1,
               CIDBMS=2,CIDOUM=2,HYBBMS=2,HYBOUM=2)
pars.list <- list(CDBM1 = c(1, 5, 10),
                  CDBMS = c(1, 1, 10, 10),  # mk, sigma1, sigma2, theta
                  CDOU1 = c(1, 10, 5, 10),
                  CDOUM = c(1, 10, 5, 5, 10),  # mk, alpha, sigma, theta1, theta2
                  CIDBMS = c(1, 1, 10, 10), # mk, sigma1, sigma2, theta
                  CIDOUM = c(1, 10, 5, 5, 10),  # mk, alpha, sigma, theta1, theta2
                  HYBBMS = c(1, 0.1, 1, 5, 10, 10), # mk, sigma1, sigma2, sigma3, sigma4, theta
                  HYBOUM = c(1, 10, 5, 15, 10, 5, 20)) # mk, alpha, sigma, theta1, theta2, theta3, theta4

# prelims
n.iter <- 84
nmaps <- 100
ncores <- 84
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

for(i in 1:length(index.cor)){
  model.name <- names(pars.list)[i]
  pars <- pars.list[[i]]
  sim.index.cor <- index.cor[[i]]
  sim.index.ou <- index.ou[[i]]
  sim.rate.cat <- rate.cats[i]
  mclapply(1:n.iter, function(x) ModelSetRun(phy, pars, model.name, sim.index.cor, sim.index.ou, sim.rate.cat, index.cor, index.ou, rate.cats, nmaps, x), mc.cores = ncores)
}

singleSummary <- function(Rsave){
  load(Rsave)
  cat("\nBegining", Rsave, "\n")
  SumSet <- summarize.hOUwie.set(out[[2]], prune = FALSE, nSim = 100)
  file.name <- gsub("-Pars", "-Summary-Pars", Rsave)
  save(SumSet, file = file.name)
}

for(i in dir()){
  Rsaves <- dir(i, full.names = TRUE)
  nCores <- length(Rsaves)
  mclapply(Rsaves, function(x) singleSummary(x), mc.cores = nCores)
}


####### ####### ####### ####### ####### ####### ####### ####### ####### 
## For fitting to a super model
####### ####### ####### ####### ####### ####### ####### ####### ####### 

sim.HYBOUMVA.cor <- getRateCatMat(2)
rate.cat.mat <- equateStateMatPars(getRateCatMat(3), 1:6)
sim.index.cor <- getFullMat(list(sim.HYBOUMVA.cor, sim.HYBOUMVA.cor, sim.HYBOUMVA.cor), rate.cat.mat)
sim.index.ou <- getOUParamStructure("OUMVA", "three.point", FALSE, FALSE, dim(sim.index.cor)[1])

n.iter <- 84
nmaps <- 250
ncores <- 84
nTip <- 500
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
p.mk <- unique(sim.index.cor[sim.index.cor > 0])/sum(p.mk <- unique(sim.index.cor[sim.index.cor > 0])) * 5
pars <- c(p.mk,1,2,1e-5,1e-5,1e-5, 4, # R1 is OUM
          2,1,10,4,4,6, # R2 is BMS
          5,10,5,10,5,15) # R3 is BM + OUMA
model.name <- "HYBOUMVA"
sim.rate.cat <- 3
index.cor <- list(fit.CDBMS.cor, fit.CDOUM.cor, fit.CIDBMS.cor, fit.CIDOUM.cor, fit.HYBBMS.cor, fit.HYBOUM.cor)
index.ou <- list(fit.CDBMS.ou, fit.CDOUM.ou, fit.CIDBMS.ou, fit.CIDOUM.ou, fit.HYBBMS.ou, fit.HYBOUM.ou)
rate.cats <- c(1,1,2,2,2,2)
mclapply(1:n.iter, function(x) ModelSetRun(phy, pars, model.name, sim.index.cor, sim.index.ou, sim.rate.cat, index.cor, index.ou, rate.cats, nmaps, x), mc.cores = ncores)



####### ####### ####### ####### ####### ####### ####### ####### ####### 
# Simulate just under CID OUM and HYBOUM since that had distinguishability issues
####### ####### ####### ####### ####### ####### ####### ####### ####### 

# start by varrying the signal in the OU process, next try to adjust the discrete process adn get more HRM transitions

index.cor <- list(fit.CDBMS.cor, fit.CDOUM.cor, fit.CIDBMS.cor, fit.CIDOUM.cor, fit.HYBBMS.cor, fit.HYBOUM.cor)
index.ou <- list(fit.CDBMS.ou, fit.CDOUM.ou, fit.CIDBMS.ou, fit.CIDOUM.ou, fit.HYBBMS.ou, fit.HYBOUM.ou)
rate.cats <- c(1,1,2,2,2,2)
pars.list <- list(CIDOUM = c(2, 1, 0.1, 5, 20))  # mk, alpha, sigma, theta1, theta2
sim.index.cor.all <- index.cor[c(4,6)]
sim.index.ou.all <- index.ou[c(4,6)]
sim.rate.cats <- c(2,2)
# prelims
n.iter <- 84
nmaps <- 250
ncores <- 84
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

for(i in 1:length(sim.index.cor.all)){
  model.name <- names(pars.list)[i]
  pars <- pars.list[[i]]
  sim.index.cor <- sim.index.cor.all[[i]]
  sim.index.ou <- sim.index.ou.all[[i]]
  sim.rate.cat <- sim.rate.cats[i]
  mclapply(1:n.iter, function(x) ModelSetRun(phy, pars, model.name, sim.index.cor, sim.index.ou, sim.rate.cat, index.cor, index.ou, rate.cats, nmaps, x), mc.cores = ncores)
}


singleSummary <- function(Rsave){
  load(Rsave)
  cat("\nBegining", Rsave, "\n")
  SumSet <- summarize.hOUwie.set(out[[2]], prune = FALSE, nSim = 100)
  file.name <- gsub("-Pars", "-Summary-Pars", Rsave)
  save(SumSet, file = file.name)
}

for(i in dir()){
  Rsaves <- dir(i, full.names = TRUE)
  nCores <- length(Rsaves)
  mclapply(Rsaves, function(x) singleSummary(x), mc.cores = nCores)
}



# old code below

####### ####### ####### ####### ####### ####### ####### ####### ####### 
# can we detect hidden states when they're present?
####### ####### ####### ####### ####### ####### ####### ####### ####### 
# from previous work we know that we can detect hidden states when there are differences in the state dependent process with discrete, but what about OU differences?
sim.index.cor <- fit.HYBBMS.cor
sim.index.ou <- fit.HYBBMS.ou
sim.rate.cat <- 2

fit.index.cor <- list(fit.CDBMS.cor, fit.CDOUM.cor, fit.CIDBMS.cor, fit.CIDOUM.cor, fit.HYBBMS.cor, fit.HYBOUM.cor)
fit.index.ou <- list(fit.CDBMS.ou, fit.CDOUM.ou, fit.CIDBMS.ou, fit.CIDOUM.ou, fit.HYBBMS.ou, fit.HYBOUM.ou)
fit.rate.cat <- c(1,1,2,2,2,2)

mclapply(1:n.iter, function(x) ModelSetRun(phy, pars, "HYB_BMS", sim.index.cor, sim.index.ou, sim.rate.cat, fit.index.cor, fit.index.ou, fit.rate.cat, nmaps, x), mc.cores = ncores)


####### ####### ####### ####### ####### ####### ####### ####### ####### 
# can we avoid detecting hidden states when they're absent?
####### ####### ####### ####### ####### ####### ####### ####### ####### 
# there are often problems in PCMs where overly complex models are selected. one of the advantages of hOUwie is that it explictly models the process of regimes and is time heterogenious. this may assist it avoid the large-p-small-n problem, but we should make sure it doesn't always select overly complex models. 
sim.index.cor <- fit.CDBMS.cor
sim.index.ou <- fit.CDBMS.ou
sim.rate.cat <- 1

fit.index.cor <- list(fit.CDBMS.cor, fit.CDOUM.cor, fit.CIDBMS.cor, fit.CIDOUM.cor, fit.HYBBMS.cor, fit.HYBOUM.cor)
fit.index.ou <- list(fit.CDBMS.ou, fit.CDOUM.ou, fit.CIDBMS.ou, fit.CIDOUM.ou, fit.HYBBMS.ou, fit.HYBOUM.ou)
fit.rate.cat <- c(1,1,2,2,2,2)

pars <- c(1,
          1, 10,
          5)

mclapply(1:n.iter, function(x) ModelSetRun(phy, pars, "CD_BMS", sim.index.cor, sim.index.ou, sim.rate.cat, fit.index.cor, fit.index.ou, fit.rate.cat, nmaps, x), mc.cores = ncores)

####### ####### ####### ####### ####### ####### ####### ####### ####### 
# can we detect character dependence when it's present?
####### ####### ####### ####### ####### ####### ####### ####### ####### 
# a common question when using hOUwie could be does my focal trait associate with my ecological variable? in this case we should be able to recover character dependence 
sim.index.cor <- fit.CDOUM.cor
sim.index.ou <- fit.CDOUM.ou
sim.rate.cat <- 1

fit.index.cor <- list(fit.CDBMS.cor, fit.CDOUM.cor, fit.CIDBMS.cor, fit.CIDOUM.cor, fit.HYBBMS.cor, fit.HYBOUM.cor)
fit.index.ou <- list(fit.CDBMS.ou, fit.CDOUM.ou, fit.CIDBMS.ou, fit.CIDOUM.ou, fit.HYBBMS.ou, fit.HYBOUM.ou)
fit.rate.cat <- c(1,1,2,2,2,2)

pars <- c(1,
          1, 10,
          5)

mclapply(1:n.iter, function(x) ModelSetRun(phy, pars, "HYB_OUM", sim.index.cor, sim.index.ou, sim.rate.cat, fit.index.cor, fit.index.ou, fit.rate.cat, nmaps, x), mc.cores = ncores)


####### ####### ####### ####### ####### ####### ####### ####### ####### 
# can we avoid detecting character dependence when it's absent?
####### ####### ####### ####### ####### ####### ####### ####### ####### 
sim.index.cor <- fit.CIDOUM.ou
sim.index.ou <- fit.CIDOUM.ou
sim.rate.cat <- 1

fit.index.cor <- list(fit.CDBMS.cor, fit.CDOUM.cor, fit.CIDBMS.cor, fit.CIDOUM.cor, fit.HYBBMS.cor, fit.HYBOUM.cor)
fit.index.ou <- list(fit.CDBMS.ou, fit.CDOUM.ou, fit.CIDBMS.ou, fit.CIDOUM.ou, fit.HYBBMS.ou, fit.HYBOUM.ou)
fit.rate.cat <- c(1,1,2,2,2,2)

pars <- c(1,
          1, 10,
          5)

mclapply(1:n.iter, function(x) ModelSetRun(phy, pars, "HYB_OUM", sim.index.cor, sim.index.ou, sim.rate.cat, fit.index.cor, fit.index.ou, fit.rate.cat, nmaps, x), mc.cores = ncores)


####### ####### ####### ####### ####### ####### ####### ####### ####### 
# how does hOUwie perform when the model is outside the original set?
####### ####### ####### ####### ####### ####### ####### ####### ####### 
# simulate both CD and CID 


ModelSetRun <- function(phy, pars, model.name, sim.index.cor, sim.index.ou, sim.rate.cat, fit.index.cor, fit.index.ou, fit.rate.cat, nmaps, iter){
  nModels <- length(fit.index.ou)
  nState <- dim(sim.index.cor)[1]
  cat("Begining single run for", nModels, "models...\n")
  data.houwie <- generateData(phy, sim.index.cor, sim.index.ou, pars)
  while(length(unique(data.houwie$data[,2])) != nState){
    data.houwie <- generateData(phy, index.cor, index.ou, pars)
  }
  if(sim.rate.cat == 2){
    data <- data.houwie$data
    data[data[,2] == 3, 2] <- 1
    data[data[,2] == 4, 2] <- 2
    data[,2] <- as.factor(as.numeric(data[,2]))
  }else{
    data <- data.houwie$data
  }
  obj <- vector("list", nModels)
  for(j in 1:nModels){
    obj[[j]] <- hOUwie(phy=phy, data=data, rate.cat=fit.rate.cat[j],
                       index.cor=fit.index.cor[[j]], root.p="yang",
                       index.ou=fit.index.ou[[j]], root.station=FALSE, get.root.theta=FALSE, mserr="none",
                       nSim=nmaps, nCores=1)
  }
  out <- list(simulating.data = data.houwie,
              model.fits = obj)
  file.name <- paste0("/space_2/jamesboyko/2020_hOUwie/ModelTesting/", model.name, "/Sim=", model.name, "-", "Pars=", paste(pars, collapse = "_"), "-Iter=", iter, ".Rsave")
  save(out, file = file.name)
  cat("Done.\n")
  return(out)
}

