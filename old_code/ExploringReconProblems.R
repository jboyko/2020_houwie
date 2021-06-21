source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(gridExtra)
require(geiger)
require(dentist)

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

hOUwie(mod$phy, mod$data, mod$rate.cat, mod$nSim, mod$nBins, mod$collapse, "joint", mod$model.cor, mod$index.cor, mod$root.p, mod$lb.cor, mod$ub.cor, mod$model.ou, mod$index.ou, p = mod$p)

# lik surface
dentFuc1 <- function(par, mod){
  out <- hOUwie(mod$phy, mod$data, mod$rate.cat, mod$nSim, mod$nBins, mod$collapse, "joint", mod$model.cor, mod$index.cor, mod$root.p, mod$lb.cor, mod$ub.cor, mod$model.ou, mod$index.ou, p = par, split.LnLiks = FALSE)
  return(-out$loglik)
}

names(mod$p) <- c("mk", "alpha", "sigma", "theta1", "theta2")
dented_results <- dent_walk(par=mod$p, fn=dentFuc1, best_neglnL=59.01998, delta=4, nsteps=500, print_freq=250, mod=mod)
print(dented_results)
plot(dented_results)


# examining simulating parameters to check for optimization errors
TrueCIDOUMFits <- list()
for(i in 1:length(ModelFits)){
  mod <- ModelFits[[i]]$ModelFits[[6]]
  mod$collapse <- TRUE
  mod$nBins <- 10
  mod$p <- DataSets[[i]]$pars
  TrueCIDOUMFits[[i]] <- hOUwie(mod$phy, mod$data, mod$rate.cat, mod$nSim, mod$nBins, mod$collapse, "joint", mod$model.cor, mod$index.cor, mod$root.p, mod$lb.cor, mod$ub.cor, mod$model.ou, mod$index.ou, p = mod$p)
  
}
i = 26
j = 6
ModelFits[[i]]$ModelTable
mod <- ModelFits[[i]]$ModelFits[[j]]
mod$collapse <- TRUE
# mod$p <- getPVecFromModel(mod)
mod$p <- DataSets[[i]]$pars

# creating a model table with split likelihoods 
BM1Fits <- CIDBMSFits <- CIDOUMFits <- TrueCIDOUMFits <- list()
for(i in 1:length(ModelFits)){
  mod <- ModelFits[[i]]$ModelFits[[6]]
  mod$collapse <- TRUE
  mod$nBins <- 10
  mod$p <- DataSets[[i]]$pars
  TrueCIDOUMFits[[i]] <- hOUwie(mod$phy, mod$data, mod$rate.cat, mod$nSim, mod$nBins, mod$collapse, "joint", mod$model.cor, mod$index.cor, mod$root.p, mod$lb.cor, mod$ub.cor, mod$model.ou, mod$index.ou, p = mod$p, split.LnLiks = TRUE)
  
  mod <- ModelFits[[i]]$ModelFits[[6]]
  mod$collapse <- TRUE
  mod$nBins <- 10
  mod$p <- getPVecFromModel(mod)
  CIDOUMFits[[i]] <- hOUwie(mod$phy, mod$data, mod$rate.cat, mod$nSim, mod$nBins, mod$collapse, "joint", mod$model.cor, mod$index.cor, mod$root.p, mod$lb.cor, mod$ub.cor, mod$model.ou, mod$index.ou, p = mod$p, split.LnLiks = TRUE)
  
  mod <- ModelFits[[i]]$ModelFits[[5]]
  mod$collapse <- TRUE
  mod$nBins <- 10
  mod$p <- getPVecFromModel(mod)
  CIDBMSFits[[i]] <- hOUwie(mod$phy, mod$data, mod$rate.cat, mod$nSim, mod$nBins, mod$collapse, "joint", mod$model.cor, mod$index.cor, mod$root.p, mod$lb.cor, mod$ub.cor, mod$model.ou, mod$index.ou, p = mod$p, split.LnLiks = TRUE)
  
  mod <- ModelFits[[i]]$ModelFits[[1]]
  mod$collapse <- TRUE
  mod$nBins <- 10
  mod$p <- getPVecFromModel(mod)
  BM1Fits[[i]] <- hOUwie(mod$phy, mod$data, mod$rate.cat, mod$nSim, mod$nBins, mod$collapse, "joint", mod$model.cor, mod$index.cor, mod$root.p, mod$lb.cor, mod$ub.cor, mod$model.ou, mod$index.ou, p = mod$p, split.LnLiks = TRUE)
  
}

ToSaveCSV <- list()
for(i in 1:length(BM1Fits)){
  tmp <- list(BM1 = BM1Fits[[i]],
       CIDBMS = CIDBMSFits[[i]],
       CIDOUM = CIDOUMFits[[i]],
       TruCIDOUM = TrueCIDOUMFits[[i]])
  tmpTable <- round(getModelTable(tmp), 3)
  tmpOULik <- unlist(lapply(tmp, function(x) x$OULnLik))
  tmpMkLik <- unlist(lapply(tmp, function(x) x$MkLnLik))
  tmpTable <- cbind(Model = rownames(tmpTable), tmpTable[,1:2], OUlnLik = round(tmpOULik, 3), MklnLik = round(tmpMkLik, 3), tmpTable[,3:5])
  ToSaveCSV[[i]] <- rbind(tmpTable, "", "")
  rownames(ToSaveCSV[[i]]) <- NULL
}



out.obj <- do.call(rbind, ToSaveCSV)
write.csv(out.obj, file = "~/Desktop/ModelTables.csv")
