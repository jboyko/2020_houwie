# comparison of BMS and BM1
source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")
require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)

#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
# check the joint recon of pupko in corhmm and houwie
#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 

nTip <- 100
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(1, 2, 1, 2, 10)  # mk, alpha, sigma, theta1, theta2

data.houwie <- generateData(phy, fit.cor, fit.ou, pars)

p = 2
corHMMMod <- corHMM(data.houwie$simmap[[1]], data.houwie$data[,c(1,2)], 1, model = "ER", p = p)
corHMMASR <- ancRECON(data.houwie$simmap[[1]], data.houwie$data[,c(1,2)], p, method=c("joint"),
                  1, ntraits=NULL, rate.mat=NULL, 
                  model="ER", root.p=NULL, get.likelihood=FALSE, get.tip.states = FALSE, collapse = TRUE)

hOUwieMod <- hOUwie(data.houwie$simmap[[1]], data.houwie$data, 1, 10, TRUE, "joint", model.cor = "ER", model.ou = "OUM", p = pars)
hOUwieASR <- hOUwieRecon(hOUwieMod)

# check if the recons match
CombinedASR <- cbind(hOUwieASR$NodeStates[,1], corHMMASR$lik.anc.states)
all(apply(CombinedASR, 1, function(x) x[1] == x[2]))

# check if the Mk logliks match
corHMMMod$loglik == hOUwieMod$MkLnLik

# visualize the difference
par(mfrow=c(1,2))
cols<-setNames(c("black","grey"),
               c("1","2"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1", "2"), pch=16, col = cols)
nodelabels(pch=16, col = cols[corHMMASR$lik.anc.states])
cols<-setNames(c("black","grey"),
               c("1","2"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1", "2"), pch=16, col = cols)
nodelabels(pch=16, col = cols[hOUwieASR$NodeStates[,1]])

#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
# check the branchwise average of an OU process
#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
z_1 <- 10
z0_1 <- 5
t_1 <- 1
alpha_1 = 2
sigma_1 = 1
theta_1 = 10

z_2 <- 5
z0_2 <- 7.5
t_2 <- 1
alpha_2 = 2
sigma_2 = 1
theta_2 = 10


OU_llik <- dOU(z, z0, t, alpha, theta, sigma, TRUE)





#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
# next test
#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 


data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
data[,2] <- as.numeric(data[,2])
data$reg <- c("Y", "N")[data$reg]
data <- data.frame(sp = data$sp, X =data$reg, Y = data$reg, x = data$x)

nBins <- 10
rate.cat <- 2

test4 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1")
test4$objective
exp(test4$solution)
p <- exp(test4$solution)
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1", p = p)

test3 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
test3$objective
exp(test3$solution)
p <- c(p[1], p[2], p[2], p[3])
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)

# hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = exp(test3$solution))
BMS.recon <- hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = exp(test3$solution))

