source("~/2020_hOUwie/hOUwieSimmap.R")
source("~/2020_hOUwie/Utils.R")

# source("/space_2/jamesboyko/2020_hOUwie/hOUwieSimmap.R")
# source("/space_2/jamesboyko/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)

# organize the data
# species with 80% or more of their diet being plant matter are herbivores, 20% or less are carnivore, else omnivore
data <- read.csv("~/2020_hOUwie/empirical/mammals/MammalDataSmall.csv")
data <- read.csv("MammalDataSmall.csv")
dat <- data.frame(Sp = gsub(" ", "_", data$Species), Diet = data$Main_food, Mass = data$Body_mass_g)
dat$Diet <- as.character(dat$Diet)
VegiScore <- data$Seeds + data$Fungi + data$Flowers.Gum + data$Roots.Tubers + data$Green.Plants + data$Fruit
dat[VegiScore >= 80, 2] <- "Herbivore"
dat[VegiScore <= 20, 2] <- "Carnivore"
dat[VegiScore < 80 & VegiScore > 20 , 2] <- "Omnivore"
dat$Mass <- log(dat$Mass)

# organize the phylgoeny
MCC_100 <- read.nexus("~/2020_hOUwie/empirical/mammals/tree-pruner-bdbf3821-fb61-418b-a45a-d363261c3c8f/output.nex")
MCC_100 <- read.nexus("output.nex")
phy <- MCC_100[[1]]
toDropPhy <- phy$tip.label[!phy$tip.label %in% dat$Sp]
toDropDat <- dat$Sp[!dat$Sp %in% phy$tip.label]

# drop species not matching in the phylogeny or dataset
dat <- dat[!dat$Sp %in% toDropDat,]
phy <- drop.tip(phy, toDropPhy)
dim(dat)[1] == length(phy$tip.label)

# randomly resolve the phylogeny to be completely bifurcating
phy <- multi2di(phy)

# set up our discrete character models
# discrete model A allows for transitions between carnivores and herbivores
StateDepA <- equateStateMatPars(getStateMat4Dat(dat[c(1,2)])$rate.mat, 1:2)
ParamDepA <- equateStateMatPars(getRateCatMat(2), 1:2)
DiscModelA <- getFullMat(list(StateDepA, StateDepA), ParamDepA)

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

# fit the models
nSim <- 50
mserr <- "none"
dual = FALSE
collapse <- TRUE
root.station <- FALSE
get.root.theta <- FALSE
lb.disc = 0.001
ub.disc = NULL
lb.cont = NULL
ub.cont = c(5, 10, 20)
opts = NULL

fitBM1 <- hOUwie(phy, dat, 2, nSim = nSim, index.cor = DiscModelA, model.ou = "BM1", weighted = TRUE, lb.cor = lb.disc)
fitOU1 <- hOUwie(phy, dat, 2, nSim = nSim, index.cor = DiscModelA, model.ou = "OU1", weighted = TRUE, lb.cor = lb.disc)
fitBMS <- hOUwie(phy, dat, 2, nSim = nSim, index.cor = DiscModelA, index.ou = CDBMS, weighted = TRUE, lb.cor = lb.disc)
fitOUM <- hOUwie(phy, dat, 2, nSim = nSim, index.cor = DiscModelA, index.ou = CDOUM, weighted = TRUE, lb.cor = lb.disc)
fitCIDOUM <- hOUwie(phy, dat, 2, nSim = nSim, index.cor = DiscModelA, index.ou = CIDOUM, weighted = TRUE, lb.cor = lb.disc)
fitHYBOUM <- hOUwie(phy, dat, 2, nSim = nSim, index.cor = DiscModelA, index.ou = HYBOUM, weighted = TRUE, lb.cor = lb.disc)
fitCIDBMS <- hOUwie(phy, dat, 2, nSim = nSim, index.cor = DiscModelA, index.ou = CIDBMS, weighted = TRUE, lb.cor = lb.disc)

getModelTable(list("BM1" = fitBM1, "OU1" = fitOU1, "CDBMS" = fitBMS, "CDOUM" = fitOUM, "CIDBMS" = fitCIDBMS, "CIDOUM" = fitCIDOUM, "HYBOUM" = fitHYBOUM))




