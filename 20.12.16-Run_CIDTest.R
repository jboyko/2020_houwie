# This script will test for an improvement of bias/variance for parameters when doing corHMM then OUwie vs. hOUwie formulated following Revell (2011)

# imports
require(corHMM)
require(OUwie)
require(parallel)
require(geiger)

# generate data 
nSim <- 1
p.mk <- 1
nTip <- 100
root.p = c(1, 0, 0, 0)
theta = c(3, 3, 10, 10)
theta0 = 3
alpha = c(1, 1, 1, 1)
sig2= c(3, 3, 1, 1)

phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip)
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
R1 <- R2 <- equateStateMatPars(getRateCatMat(2), list(c(1,2)))
FullMat <- getFullMat(list(R1, R2), R1)
rates <- c(1, 2, 0.5)
Q = FullMat
Q[!Q==0] = rates[FullMat]
diag(Q) <- -rowSums(Q)

singleRun <- function(){
  data <- OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta, 1)
  TrueMap <- data[[2]][[1]]
  data.true <- data[[1]]
  data.obs <- data[[1]]
  data.obs[,2][data.obs[,2] == 3] <- 1
  data.obs[,2][data.obs[,2] == 4] <- 2
  # fit the data to TrueMap, then corHMM + OUwie, then hOUwie
  fitOU_tru <- OUwie(TrueMap, data.true, model = "OUMV", simmap.tree = TRUE, scaleHeight = FALSE, algorithm = "three.point")
  fitMK <- corHMM(phy, data.obs[,c(1,2)], 2, rate.mat = FullMat, get.tip.states = TRUE)
  tmp <- apply(fitMK$tip.states, 1, function(x) sample(1:4, 1, prob = x))
  while(!all(1:4 %in% tmp)){
    tmp <- apply(fitMK$tip.states, 1, function(x) sample(1:4, 1, prob = x))
  }
  data.est <- data.true 
  data.est[,2] <- tmp
  simmaps <- makeSimmap(phy, data.est[,c(1,2)], fitMK$solution, 1, nSim = nSim)
  fitOU_est <- lapply(simmaps, function(x) OUwie(x, data.est, model = "OUMV", simmap.tree = TRUE, scaleHeight = FALSE, algorithm = "three.point"))
  res_hOUwie <- OUwie:::hOUwie(phy, data.obs, rate.cat=2, rate.mat = FullMat, model.ou = "OUMV", nSim = nSim)
  res_hOUwie.wt <- OUwie:::hOUwie(phy, data.obs, rate.cat = 2, rate.mat = FullMat, model.ou = "OUMV", nSim = nSim, weighted = TRUE)
  return(list(fitMK = fitMK,
              fitOU_tru = fitOU_tru,
              fitOU_est = fitOU_est,
              res_hOUwie = res_hOUwie,
              res_hOUwie.wt = res_hOUwie.wt))
}

# for server
# LiamTest <- mclapply(1:100, function(x) singleRun(), mc.cores = 50)

# for local testing
# test <- lapply(1:3, function(x) singleRun())

# analyze for bias

# transition rate of 1
  LiamTest <- LiamTest[unlist(lapply(LiamTest, function(x) class(x))) != "try-error"]
  sig2_TrueMapEst <- do.call(rbind, lapply(LiamTest, function(x) x$fitOU_tru$solution[2,]))
  sig2_StocMapEst <- lapply(LiamTest, function(x) do.call(rbind, (lapply(x$fitOU_est, function(y) y$solution[2,]))))
  sig2_StocMapEst <- do.call(rbind, sig2_StocMapEst)
  sig2_hOUwieTrue <- do.call(rbind, lapply(LiamTest, function(x) exp(x$res_hOUwie$solution)[2:3]))
  sig2_hOUwieWt <- do.call(rbind, lapply(LiamTest, function(x) exp(x$res_hOUwie.wt$solution)[2:3]))
  
  q <- paste("q", p.mk, sep="_")
df.Q1_Sig2_1 <- rbind(
  data.frame(par = q, fit = "TrueMap", val = sig2_TrueMapEst[,1]),
  data.frame(par = q, fit = "StocMap", val = sig2_StocMapEst[,1]),
  data.frame(par = q, fit = "hOUwTru", val = sig2_hOUwieTrue[,1]),
  data.frame(par = q, fit = "hOuwWt", val = sig2_hOUwieWt[,1])
)
df.Q1_Sig2_2 <- rbind(
  data.frame(par = q, fit = "TrueMap", val = sig2_TrueMapEst[,2]),
  data.frame(par = q, fit = "StocMap", val = sig2_StocMapEst[,2]),
  data.frame(par = q, fit = "hOUwTru", val = sig2_hOUwieTrue[,2]),
  data.frame(par = q, fit = "hOuwWt", val = sig2_hOUwieWt[,2])
)

file.name <- paste(format(Sys.time(), "%y_%m_%d"), "-Test_q", p.mk, "_hOUwieVsLiam.R", sep = "")
save(df.Q1_Sig2_1, df.Q1_Sig2_2, file = file.name)

load("~/2020_hOUwie/BoykoTest.Rsave")
require(viridis)
require(ggplot2)

cols <- viridis(4)

ggplot(df.Q1_Sig2_1, aes(x=par, y=val, fill=fit)) + 
  labs(x = "Transition Rate", y = expression('estimated '*sigma[1]^2)) +
  scale_fill_manual(values=cols) + 
  theme(text = element_text(size = 20)) + 
  geom_boxplot()

ggplot(df.Q1_Sig2_2, aes(x=par, y=val, fill=fit)) + 
  labs(x = "Transition Rate", y = expression('estimated '*sigma[2]^2)) +
  scale_fill_manual(values=cols) + 
  theme(text = element_text(size = 20)) + 
  geom_boxplot()





