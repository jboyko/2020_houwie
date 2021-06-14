source("~/2020_hOUwie/hOUwieSimmap.R")
source("~/2020_hOUwie/Utils.R")

# source("/space_2/jamesboyko/2020_hOUwie/hOUwieSimmap.R")
# source("/space_2/jamesboyko/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)

# organize the data
# dat <- read.csv("Ericaceae_niche.csv")
dat <- read.csv("~/2021_SeedDispersal/trait_data/Ericaceae_niche.csv")
dat <- data.frame(sp = dat$species, reg = dat$Fruit_type, arid = dat$mean_aridity)
dat <- dat[which(apply(dat, 1, function(x) !any(is.na(x)))),]

# organize the phylgoeny
# phy <- read.tree("Ericaceae_Schwery_etal_2015.tre")
phy <- read.tree("~/2021_SeedDispersal/trees/Ericaceae_Schwery_etal_2015.tre")

# drop species not matching in the phylogeny or dataset
toDropPhy <- phy$tip.label[!phy$tip.label %in% dat$sp]
toDropDat <- dat$Sp[!dat$sp %in% phy$tip.label]
dat <- dat[!dat$sp %in% toDropDat,]
phy <- drop.tip(phy, toDropPhy)
dim(dat)[1] == length(phy$tip.label)
dat <- dat[match(phy$tip.label, dat$sp), ]

# set up our discrete character models
# discrete model A allows for transitions between carnivores and herbivores
StateDepA <- equateStateMatPars(getStateMat4Dat(dat[c(1,2)])$rate.mat, 1:2)
ParamDepA <- equateStateMatPars(getRateCatMat(2), 1:2)
DiscModelA <- getFullMat(list(StateDepA, StateDepA), ParamDepA)

# visualize and attack
Tmax <- max(branching.times(phy))
plot(phy, show.tip.label = FALSE, x.lim = c(0, Tmax + 0.2 * Tmax))
start = Tmax + (0.005 * Tmax)
dat.prop <- (dat[,3] - min(dat[,3]))/max(dat[,3] - min(dat[,3]))
xadd = Tmax * 0.2
jitter <- 0.1 * xadd
cols = c("brown", "purple")
for(i in 1:length(dat.prop)){
  lines(list(x = c(start, start + (dat.prop[i] * xadd)), y = c(i,i)), 
        col = cols[ifelse(dat$reg[i] == "Dry", 1, 2)],
        lwd = 1)
}
legend("bottomleft", legend = c("Dry/ Abiotic","Fleshy/ Biotic"), pch=15, col = cols)

# set up our continuous character models 
# variable brownian motion models set up to be character dependent, independent, and a mix of both
HYBBMS <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(DiscModelA)[1])
CDBMS <- CIDBMS <- HYBBMS
CDBMS[2,] <- c(1,2,1,2)
CDBMS[3,] <- 3
CIDBMS[2,] <- c(1,1,2,2)
CIDBMS[3,] <- 3
# OUM models set up to be character dependent, independent, and a mix of both
HYBOUM <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(DiscModelA)[1])
CDOUM <- CIDOUM <- HYBOUM
CDOUM[3,] <- c(3,4,3,4)
CIDOUM[3,] <- c(3,3,4,4)

HYBOUM <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(DiscModelA)[1])
CDOUMV <- CIDOUMV <- HYBOUM
CDOUMV[2,] <- c(2,3,2,3)
CDOUMV[3,] <- c(4,5,4,5)
CIDOUMV[3,] <- c(2,2,3,3)
CIDOUMV[3,] <- c(4,4,5,5)
CDOUM_ABIOTIC <- CDOUM_BIOTIC <- HYBOUM
CDOUM_BIOTIC[3,] <- c(3,4,3,5)
CDOUM_ABIOTIC[3,] <- c(3,4,5,4)


# fit the models
nSim <- 500
mserr <- "none"
dual = FALSE
collapse <- TRUE
root.station <- FALSE
get.root.theta <- FALSE
lb.disc = 0.00001
ub.disc = NULL
lb.cont = NULL
ub.cont = c(5, 10, 20)
opts = NULL


singleRun <- function(index){
  if(index == 1){
    model.name = "BM1"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, model.ou = "BM1", weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 2){
    model.name = "OU1"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, model.ou = "OU1", weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 3){
    model.name = "CDBMS"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CDBMS, weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 4){
    model.name = "CDOUM"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CDOUM, weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 5){
    model.name = "CIDBMS"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CIDBMS, weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 6){
    model.name = "CIDOUM"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CIDOUM, weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 7){
    model.name = "HYBBMS"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = HYBBMS, weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 8){
    model.name = "HYBOUM"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = HYBOUM, weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 9){
    model.name = "CDOUMV"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CDOUMV, weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 10){
    model.name = "CIDOUMV"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CIDOUMV, weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 11){
    model.name = "CDOUM_BIOTIC"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CDOUM_BIOTIC, weighted = TRUE, lb.cor = lb.disc)
  }
  if(index == 12){
    model.name = "CDOUM_ABIOTIC"
    fit <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CDOUM_ABIOTIC, weighted = TRUE, lb.cor = lb.disc)
  }
  save(fit, file = paste0("FitSD=", model.name, ".Rsave"))
  return(fit)
}

res <- mclapply(1:12, function(x) singleRun(x), mc.cores = 12)


getModelTable(list("BM1" = fitBM1, "OU1" = fitOU1, "CDBMS" = fitBMS, "CDOUM" = fitOUM, "CIDBMS" = fitCIDBMS, "CIDOUM" = fitCIDOUM, "HYBOUM" = fitHYBOUM))
