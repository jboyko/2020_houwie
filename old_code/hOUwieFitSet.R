source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")

source("hOUwie.R")
source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(gridExtra)
require(geiger)

fitModelset <- function(phy, data){
  CID.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  CID.cor <- getFullMat(list(CID.cor, CID.cor), CID.cor)
  CID.cor <- equateStateMatPars(CID.cor, c(1,2,3))
  CIDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(CID.cor)[1])
  CIDOUM.ou[3,] <- c(3,3,4,4)
  CIDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(CID.cor)[1])
  CIDBMS.ou[2,] <- c(1,1,2,2)
  CIDBMS.ou[3,] <- c(3,3,3,3)
  
  BM1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1", all.roots = TRUE)
  OU1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1", all.roots = TRUE)
  OUM <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", all.roots = TRUE)
  BMS <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", all.roots = TRUE)
  CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CIDOUM.ou, all.roots = TRUE)
  CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CIDBMS.ou, all.roots = TRUE)
  
  ModelFits <- list(BM1 = BM1, OU1 = OU1, BMS = BMS, OUM = OUM, CIDBMS = CIDBMS, CIDOUM = CIDOUM)
  ModelTable <- getModelTable(ModelFits)
  return(list(ModelFits = ModelFits,
              ModelTable = ModelTable))
}

getDatasets <- function(model, nDatum){
  if(model == "BM1"){
    fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
    fit.ou <- getOUParamStructure("BM1", "three.point", FALSE, FALSE, dim(fit.cor)[1])
    pars = c(1, 5, 10)  # mk, sigma, theta1
  }
  if(model == "OU1"){
    fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
    fit.ou <- getOUParamStructure("OU1", "three.point", FALSE, FALSE, dim(fit.cor)[1])
    pars = c(1, 3, 3, 8)  # mk, alpha, sigma, theta1
  }
  if(model == "BMS"){
    fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
    fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
    pars = c(1, 1, 10, 8)  # mk, sigma1, sigma2, theta
  }
  if(model == "OUM"){
    fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
    fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
    pars = c(1, 3, 3, 3, 8)  # mk, alpha, sigma, theta1, theta2
  }
  if(model == "CIDBMS"){
    fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
    fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
    fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
    fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
    fit.ou[2,] <- c(1,1,2,2)
    fit.ou[3,] <- c(3,3,3,3)
    pars = c(1, 1, 10, 10)  # mk, sigma1, sigma2, theta
  }
  if(model == "CIDOUM"){
    fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
    
    fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
    fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
    fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
    fit.ou[3,] <- c(3,3,4,4)
    pars = c(1, 3, 3, 2, 10)  # mk, alpha, sigma, theta1, theta2
  }
  data.houwie <- lapply(1:nDatum, function(x) generateData(phy, fit.cor, fit.ou, pars))
  for(i in 1:length(data.houwie)){
    data.houwie[[i]]$data[data.houwie[[i]]$data[,2]==3,2] <- 1
    data.houwie[[i]]$data[data.houwie[[i]]$data[,2]==4,2] <- 2
  }
  return(data.houwie)
}



nTip <- 50
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
nSim <- 50
nCores <- 50

models <- c("BM1", "OU1", "BMS", "OUM", "CIDBMS", "CIDOUM")
# models <- c("CIDBMS", "CIDOUM")


for(i in 1:length(models)){
  FullData <- getDatasets(models[i], nSim)
  Datum <- lapply(FullData, function(x) x$data)
  ModelSetFits <- mclapply(Datum, function(x) fitModelset(phy, x), mc.cores = nCores)
  out <- list(Datasets = FullData, ModelSetFits = ModelSetFits)
  save(out, file = paste0(models[i], "-", nTip, "-hOUwieFit.Rsave"))
}

ModelFits <- out$ModelSetFits[!unlist(lapply(out$ModelSetFits, class)) == "try-error"]
colSums(do.call(rbind, lapply(ModelFits, function(x) x$ModelTable[,4])) == 0)

getMkParVec <- function(ModelFits){
  ModelFits <- ModelFits$ModelFits
  models <- c("BM1", "OU1", "BMS", "OUM", "CIDBMS", "CIDOUM")
  out.vec <- unlist(lapply(ModelFits, function(x) x$solution.cor[1,2]))
  names(out.vec) <- models
  return(out.vec)
}

countZeros <- function(x){
  return(length(which(x == 0)))
}



source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")

# source("hOUwie.R")
# source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(gridExtra)
require(geiger)

load("~/2020_hOUwie/ModelFits/CIDOUM-hOUwieFit.Rsave")
ModelFits <- out$ModelSetFits[!unlist(lapply(out$ModelSetFits, class)) == "try-error"]
DataSets <- out$Datasets[!unlist(lapply(out$ModelSetFits, class)) == "try-error"]

cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
i = 34
j = 6
ModelFits[[i]]$ModelTable
mod <- ModelFits[[i]]$ModelFits[[j]]
mod$collapse <- TRUE
mod$p <- getPVecFromModel(mod)
# mod$p <- DataSets[[i]]$pars
mod$nBins <- 10
recon.mod <- hOUwieRecon(hOUwie.model = mod)
# for(i in 1:length(DataSets)){
#   file.name <- paste0("~/2020_hOUwie/dataImages/DataSet-", i, ".png")
#   png(file.name, width = 960, 480)
  par(mfrow=c(1,2))
  cols<-setNames(c("black","grey", "black", "grey"),
                 c("1","2","3","4"))
  plotDataSet(DataSets[[i]]); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
  nodelabels(pch = 16, col = cols[recon.mod$NodeStates[,1]])
  cols<-setNames(c("black","black", "grey", "grey"),
                 c("1","2","3","4"))
  plotDataSet(DataSets[[i]]); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
  nodelabels(pch = 16, col = cols[recon.mod$NodeStates[,1]])
  
#   grid.table(round(ModelFits[[i]]$ModelTable, 2))
#   dev.off()
# }

# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(DataSets[[i]]); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)


# for the presentaTion
# NonCIDModelTables <- lapply(ModelFits, function(x) x$ModelTable[c(1,2,4),1:3])
# for(i in 1:length(NonCIDModelTables)){
#   dAIC <- NonCIDModelTables[[i]]$AIC - min(NonCIDModelTables[[i]]$AIC)
#   AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
#   NonCIDModelTables[[i]] <- cbind(NonCIDModelTables[[i]], dAIC = dAIC, AICwt = AICwt)
# }
# 
# boxplot(do.call(rbind, lapply(NonCIDModelTables, function(x) x$AICwt)))


