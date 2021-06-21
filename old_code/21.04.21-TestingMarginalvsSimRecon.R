source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")
#source("../../hOUwie.R")
require(phytools)
require(OUwie)
require(corHMM)
require(parallel)
require(reshape2)
require(ggplot2)
require(gridExtra)

getTipStatesFromMap <- function(Map, nState){
  BranchMaps <- Map$maps[Map$edge[,2] <= length(Map$tip.label)]
  Descendents <- Map$edge[,2][Map$edge[,2] <= length(Map$tip.label)]
  SpeciesByBranch <- Map$tip.label[Descendents]
  TipStates <- as.numeric(unlist(lapply(BranchMaps, function(x) names(x)[length(x)])))
  names(TipStates) <- SpeciesByBranch
  TipStates <- TipStates[Map$tip.label]
  TipStatesInMatrix <- matrix(0, length(TipStates), nState, dimnames = list(Map$tip.label))
  for(i in 1:length(TipStates)){
    TipStatesInMatrix[i,TipStates[i]] <- 1
  }
  return(TipStatesInMatrix)
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under HYBOUM 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we detect character dependence when it's present?
# can we avoid detecting hidden states when they're absent?

Rsaves <- dir("~/2020_hOUwie/ModelTesting/HYBOUM//", full.names = TRUE)
Rsaves_Summary <- Rsaves[grep("-Summary-", Rsaves)]
Rsaves_Fits <- Rsaves[-grep("-Summary-", Rsaves)]

SearchString <- "-Iter=21.Rsave"
Rsaves_Summary_i <- Rsaves_Summary[grep(SearchString, Rsaves_Summary)]
Rsaves_Fits_i <- Rsaves_Fits[grep(SearchString, Rsaves_Fits)]
load(Rsaves_Summary_i) #SumSet
load(Rsaves_Fits_i) #out

ModelFit <- out$model.fits[[3]]
MarginalRecon <- SumSet$model.recon[[3]]

SquaredDiff <- c()

nMap = 100
Maps <- makeSimmap(tree = ModelFit$phy, data = ModelFit$hOUwie.dat$data.cor, rate.cat = 2, model = ModelFit$solution.cor, nSim = nMap)
TipStatesPerMap <- lapply(Maps, function(x) getTipStatesFromMap(x, 4))
AvgTipStatesOverMaps <- Reduce("+", TipStatesPerMap)/length(TipStatesPerMap)
SquaredDiff <- rbind(SquaredDiff, cbind(nMap, rowMeans(abs(MarginalRecon - AvgTipStatesOverMaps))))


nMap = 250
Maps <- makeSimmap(tree = ModelFit$phy, data = ModelFit$hOUwie.dat$data.cor, rate.cat = 2, model = ModelFit$solution.cor, nSim = nMap)
TipStatesPerMap <- lapply(Maps, function(x) getTipStatesFromMap(x, 4))
AvgTipStatesOverMaps <- Reduce("+", TipStatesPerMap)/length(TipStatesPerMap)
SquaredDiff <- rbind(SquaredDiff, cbind(nMap, rowMeans(abs(MarginalRecon - AvgTipStatesOverMaps))))


nMap = 500
Maps <- makeSimmap(tree = ModelFit$phy, data = ModelFit$hOUwie.dat$data.cor, rate.cat = 2, model = ModelFit$solution.cor, nSim = nMap)
TipStatesPerMap <- lapply(Maps, function(x) getTipStatesFromMap(x, 4))
AvgTipStatesOverMaps <- Reduce("+", TipStatesPerMap)/length(TipStatesPerMap)
SquaredDiff <- rbind(SquaredDiff, cbind(nMap, rowMeans(abs(MarginalRecon - AvgTipStatesOverMaps))))


nMap = 1000
Maps <- makeSimmap(tree = ModelFit$phy, data = ModelFit$hOUwie.dat$data.cor, rate.cat = 2, model = ModelFit$solution.cor, nSim = nMap)
TipStatesPerMap <- lapply(Maps, function(x) getTipStatesFromMap(x, 4))
AvgTipStatesOverMaps <- Reduce("+", TipStatesPerMap)/length(TipStatesPerMap)
SquaredDiff <- rbind(SquaredDiff, cbind(nMap, rowMeans(abs(MarginalRecon - AvgTipStatesOverMaps))))
colnames(SquaredDiff) <- c("nMap", "Diff")
dat <- data.frame(nMap = as.factor(SquaredDiff[,1]), Diff = as.numeric(SquaredDiff[,2]))

ggplot(data = dat, aes(x = nMap, y = Diff)) +
  labs(x = "Number of stochastic maps", y = "Mean difference between\n true marginal and simmap average") + 
  theme_linedraw() +
  geom_boxplot()


