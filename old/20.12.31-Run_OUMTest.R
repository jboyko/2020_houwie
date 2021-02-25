# This script will test for an improvement of bias/variance for OUM

# imports
require(corHMM)
require(OUwie)
require(parallel)
require(geiger)

# generate data 
nSim <- 100
p.mk <- 1
nTip <- 100
root.p = c(1, 0)
theta = c(3, 8)
theta0 = 3
alpha = c(2, 2)
sig2= c(2, 2)

phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip)
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
Q = matrix(c(-p.mk,p.mk,p.mk,-p.mk), 2, 2)

singleRun <- function(){
  data <- OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta, 1)
  TrueMap <- data[[2]][[1]]
  data <- data[[1]]
  # fit the data to corHMM then OUwie and then hOUwie
  fitMK <- corHMM(phy, data[,c(1,2)], 1, model = "ER")
  fitOU_tru <- OUwie(TrueMap, data, model = "OUM", simmap.tree = TRUE, scaleHeight = FALSE, algorithm = "three.point")
  simmaps <- makeSimmap(phy, data[,c(1,2)], fitMK$solution, 1, nSim = nSim)
  fitOU_est <- lapply(simmaps, function(x) OUwie(x, data, model = "OUM", simmap.tree = TRUE, scaleHeight = FALSE, algorithm = "three.point"))
  res_hOUwie <- OUwie:::hOUwie(phy, data, 1, model.cor = "ER", model.ou = "OUM", nSim = nSim)
  res_hOUwie.wt <- OUwie:::hOUwie(phy, data, 1, model.cor = "ER", model.ou = "OUM", nSim = nSim, weighted = TRUE)
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
file.name <- paste(format(Sys.time(), "%y.%m.%d"), "LiamTest_OUMTest.Rsave", sep = "")
save(LiamTest, file = file.name)

LiamTest <- LiamTest[unlist(lapply(LiamTest, function(x) class(x))) != "try-error"]
# table order: alpha1, alpha2, sigma1, sigma2, theta1, theta2
# OU True 
OUTrueTab <- do.call(rbind, lapply(LiamTest, function(x) c(t(x$fitOU_tru$solution))))
colnames(OUTrueTab) <- c("alpha1", "alpha2", "sigma1", "sigma2", "theta1", "theta2")
# OU Est
OUEstTab <- do.call(rbind, lapply(LiamTest, function(x) colMeans(do.call(rbind, lapply(x$fitOU_est, function(y) c(t(y$solution)))))))
colnames(OUEstTab) <- c("alpha1", "alpha2", "sigma1", "sigma2", "theta1", "theta2")
# hOUwie Wt
hOUWtTab <- do.call(rbind, lapply(LiamTest, function(x) exp(x$res_hOUwie.wt$solution)))
colnames(hOUWtTab) <- c("mk", "alpha", "sigma", "theta1", "theta2")
# hOUwie True
hOUTruTab <- do.call(rbind, lapply(LiamTest, function(x) exp(x$res_hOUwie$solution)))
colnames(hOUWtTab) <- c("mk", "alpha", "sigma", "theta1", "theta2")
  

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

file.name <- paste(format(Sys.time(), "%y.%m.%d"), "LiamTest_OUMTest.R", sep = "")
save(df.Q1_Sig2_1, df.Q1_Sig2_2, file = file.name)

load("~/2020_hOUwie/20_12_17-Test_q1_OUMTest.R")
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





