setwd("~/2020_hOUwie/empirical/pnh-ms-master/R/R")
source("~/2020_hOUwie/hOUwieSimmap.R")
source("~/2020_hOUwie/Utils.R")

library(treeplyr)
require(OUwie)
require(corHMM)
require(parallel)

# import the data
scales <- read.csv("../data/scales2009.csv")
scales <- scales[!scales$species == "Carlia_fusca",] # for simplicity remove the only mixed forager
tetTree <- read.tree("../data/tetrapods.tre")
tetTree$tip.label[tetTree$tip.label == "Eumeces_schneideri"] <- "Eumeces_fasciatus"
tetTree <- multi2di(tetTree)
tetTree$edge.length[tetTree$edge.length==0] <- .Machine$double.eps
td <- make.treedata(tetTree, scales)
td <- reorder(td, "postorder")
phy <- td$phy
dat <- data.frame(sp = td$phy[['tip.label']], reg = td[['foraging.mode']], FG = td[['FG.frac']])

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
        col = cols[ifelse(dat$reg[i] == "AF", 1, 2)],
        lwd = 1)
}
legend("topleft", legend = c("AF","SW"), pch=15, col = cols)

# generate the discrete model
getStateMat4Dat(dat[c(1,2)])$legend
StateDepA <- equateStateMatPars(getStateMat4Dat(dat[c(1,2)])$rate.mat, 1:2)
ParamDepA <- equateStateMatPars(getRateCatMat(2), 1:2)
DiscModelA <- getFullMat(list(StateDepA, StateDepA), ParamDepA)
DiscModelB <- equateStateMatPars(DiscModelA, 1:2)

# generate the continuous model
HYBBMS <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(DiscModelA)[1])
SWBMS <- CDBMS <- CIDBMS <- HYBBMS
CDBMS[2,] <- c(1,2,1,2)
CDBMS[3,] <- 3
CIDBMS[2,] <- c(1,1,2,2)
CIDBMS[3,] <- 3
SWBMS[2,] <- c(1,2,1,3)
SWBMS[3,] <- 4
HYBOUM <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(DiscModelA)[1])
SWOUM <- CDOUM <- CIDOUM <- HYBOUM
CDOUM[3,] <- c(3,4,3,4)
CIDOUM[3,] <- c(3,3,4,4)
SWOUM[3,] <- c(3,4,3,5)


nSim <- 100

fit.BM1 <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, model.ou = "BM1", weighted = TRUE)
fit.OU1 <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, model.ou = "OU1", weighted = TRUE)
fit.CD <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CDOUM, weighted = TRUE)
fit.CID <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CIDOUM, weighted = TRUE)
fit.CID <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CIDOUM, weighted = TRUE)
fit.SWA <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = SWOUM, weighted = TRUE)
fit.SWB <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, index.ou = SWOUM, weighted = TRUE)
getModelTable(list(BM1 = fit.BM1, OU1 = fit.OU1, CD = fit.CD, CID = fit.CID, SWA = fit.SWA, SWB = fit.SWB))

# Interestingly, the differ- ent strategies are often quite divergent, with active foragers (AFs) relying on frequent but often slow locomotion, whereas ambush predators move rarely, but rapidly (Miles et al. 2007 and refer- ences therein), resulting in phenotypes shaped by selection for energetic efficiency and higher stamina, versus high-speed burst locomotion (Miles et al. 2007 and references therein, Pruitt 2010).



