setwd("~/2020_hOUwie/")

source("hOUwieSimmap.R")
source("Utils.R")

require(treeplyr)
require(OUwie)
require(corHMM)
require(parallel)

# import the data
scales <- read.csv("empirical/pnh-ms-master/R/data/scales2009.csv")
scales <- scales[!scales$species == "Carlia_fusca",] # for simplicity remove the only mixed forager
tetTree <- read.tree("empirical/pnh-ms-master/R/data/tetrapods.tre")
tetTree$tip.label[tetTree$tip.label == "Eumeces_schneideri"] <- "Eumeces_fasciatus"
tetTree <- multi2di(tetTree)
tetTree$edge.length[tetTree$edge.length==0] <- .Machine$double.eps
td <- make.treedata(tetTree, scales)
td <- reorder(td, "postorder")
phy <- td$phy
dat <- data.frame(sp = td$phy[['tip.label']], FM = td[['foraging.mode']], PE = td[['predator.escape']], FOG = td[['FG.frac']])

# visualize and attack
# Tmax <- max(branching.times(phy))
# plot(phy, show.tip.label = FALSE, x.lim = c(0, Tmax + 0.2 * Tmax))
# start = Tmax + (0.005 * Tmax)
# dat.prop <- (dat[,3] - min(dat[,3]))/max(dat[,3] - min(dat[,3]))
# xadd = Tmax * 0.2
# jitter <- 0.1 * xadd
# cols = c("brown", "purple")
# for(i in 1:length(dat.prop)){
#   lines(list(x = c(start, start + (dat.prop[i] * xadd)), y = c(i,i)), 
#         col = cols[ifelse(dat$reg[i] == "AF", 1, 2)],
#         lwd = 1)
# }
# legend("topleft", legend = c("AF","SW"), pch=15, col = cols)

# generate the discrete model
getStateMat4Dat(dat[c(1:3)])$legend
StateDepA <- equateStateMatPars(getStateMat4Dat(dat[c(1:3)])$rate.mat, 1:12)
ParamDepA <- equateStateMatPars(getRateCatMat(2), 1:2)
DiscModelA <- getFullMat(list(StateDepA, StateDepA), ParamDepA)
DiscModelB <- equateStateMatPars(DiscModelA, 1:2)

# generate the continuous model
HYBBMS <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(DiscModelA)[1])
SWBMS <- CDBMS <- CIDBMS <- HYBBMS
CDBMS[2,] <- c(1:5,1:5)
CDBMS[3,] <- 6
CIDBMS[2,] <- c(1,1,1,1,1,2,2,2,2,2)
CIDBMS[3,] <- 3
# SWBMS[2,] <- c(1,2,1,3)
# SWBMS[3,] <- 4
HYBOUM <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(DiscModelA)[1])
PEOUM <- FMOUM <- SWOUM <- CDOUM <- CIDOUM <- HYBOUM
CDOUM[3,] <- c(3:7,3:7)
CIDOUM[3,] <- c(3,3,3,3,3,4,4,4,4,4)
FMOUM[3,] <- c(3,3,4,4,4,3,3,4,4,4)
PEOUM[3,] <- c(3,4,5,3,4,3,4,5,3,4)
SWOUM[3,] <- c(3,3,4,4,4,3,3,5,5,5)


nSim <- 2

quickRun <- function(index){
  if(index == 1){
    fit.BM1_A <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, model.ou = "BM1")
    return(fit.BM1_A)
  }
  if(index == 2){
    fit.BM1_B <- hOUwie(phy, dat, 1, nSim, model.cor = "ER", model.ou = "BM1")
    return(fit.BM1_B)
  }
  if(index == 3){
    fit.OU1_A <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, model.ou = "OU1")
    return(fit.OU1_A)
  }
  if(index == 4){
    fit.OU1_B <- hOUwie(phy, dat, 1, nSim, model.cor = "ER", model.ou = "OU1")
    return(fit.OU1_B)
  }
  if(index == 5){
    fit.CD_A <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CDOUM)
    return(fit.CD_A)
  }
  if(index == 6){
    fit.CD_B <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, index.ou = CDOUM)
    return(fit.CD_B)
  }
  if(index == 7){
    fit.CID_A <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CIDOUM)
    return(fit.CID_A)
  }
  if(index == 8){
    fit.CID_B <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, index.ou = CIDOUM)
    return(fit.CID_B)
  }
  if(index == 9){
    fit.FM_A <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = FMOUM)
    return(fit.FM_A)
  }
  if(index == 10){
    fit.FM_B <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, index.ou = FMOUM)
    return(fit.FM_B)
  }
  if(index == 11){
    fit.PE_A <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = PEOUM)
    return(fit.PE_A)
  }
  if(index == 12){
    fit.PE_B <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, index.ou = PEOUM)
    return(fit.PE_B)
  }
  if(index == 13){
    fit.SW_A <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = SWOUM)
    return(fit.SW_A)
  }
  if(index == 14){
    fit.SW_B <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, index.ou = SWOUM)
    return(fit.SW_B)
  }
  if(index == 15){
    fit.CDBMS_A <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CDBMS)
    return(fit.CDBMS_A)
  }
  if(index == 16){
    fit.CDBMS_B <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, index.ou = CDBMS)
    return(fit.CDBMS_B)
  }
  if(index == 17){
    fit.CIDBMS_A <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelA, index.ou = CIDBMS)
    return(fit.CIDBMS_A)
  }
  if(index == 18){
    fit.CIDBMS_B <- hOUwie(phy, dat, 2, nSim, index.cor = DiscModelB, index.ou = CIDBMS)
    return(fit.CIDBMS_B)
  }
}

results <- mclapply(1:18, function(x) quickRun(x), mc.cores = 18)


names(results) <- c("BM1_A", "BM1_B", "OU1_A", "OU1_B", "CDOUM_A", "CDOUM_B", "CIDOUM_A", "CIDOUM_B", "FMOUM_A", "FMOUM_B", "PEOUM_A", "PEOUM_B", "SWOUM_A", "SWOUM_B", "CDBMS_A", "CDBMS_B", "CIDBMS_A", "CIDBMS_B")
getModelTable(results)
unlist(lapply(results, function(x) x$DiscLik))
unlist(lapply(results, function(x) x$ContLik))

# Interestingly, the differ- ent strategies are often quite divergent, with active foragers (AFs) relying on frequent but often slow locomotion, whereas ambush predators move rarely, but rapidly (Miles et al. 2007 and refer- ences therein), resulting in phenotypes shaped by selection for energetic efficiency and higher stamina, versus high-speed burst locomotion (Miles et al. 2007 and references therein, Pruitt 2010).
load("empirical/scales2009/scales-results.Rsave")
require(viridis)
getModelTable(results)
cols <- viridis(5)
names(cols) <- 1:5
par(mfrow=c(1,2))
plot(results[[6]]$RegimeMap, col = cols, lwd = 2)
legend("topleft", legend = getStateMat4Dat(dat[c(1:3)])$legend, col = cols, pch = 16)
plot(results[[12]]$RegimeMap, col = cols)

# do our results make sense with what has been said in the past?
# can we detect where the regime shift occurs?
# can we detect a CID model if it were true?
# do we suffer from the single shift along a clade issue?
results[[6]]
results[[12]]




