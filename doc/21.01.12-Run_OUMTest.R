# This script will test for an improvement of bias/variance for OUM

# imports
require(corHMM)
require(OUwie)
require(parallel)
require(geiger)

# generate data 
nSim <- 100
p.mk <- 8
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
file.name <- paste(format(Sys.time(), "%y.%m.%d"), "_OUM_MK", p.mk, "_Test.Rsave", sep = "")
save(LiamTest, file = file.name)

LiamTest <- LiamTest[unlist(lapply(LiamTest, function(x) class(x))) != "try-error"]

# table order: alpha1, alpha2, sigma1, sigma2, theta1, theta2
# OU True 
OUTrueTab <- do.call(rbind, lapply(LiamTest, function(x) c(t(x$fitOU_tru$solution))))
OUTrueTab <- OUTrueTab[,c(1,3,5,6)]
OUTrueTab <- as.data.frame(OUTrueTab)
OUTrueTab$model <- "OUTrue"
OUTrueTab <- cbind(p.mk, OUTrueTab)
colnames(OUTrueTab) <- c("mk", "alpha", "sigma","theta1", "theta2", "model")

# OU Est
OUEstTab <- do.call(rbind, lapply(LiamTest, function(x) colMeans(do.call(rbind, lapply(x$fitOU_est, function(y) c(t(y$solution)))))))
OUEstTab <- OUEstTab[,c(1,3,5,6)]
OUEstTab <- as.data.frame(OUEstTab)
OUEstTab$model <- "OUEst"
OUEstTab <- cbind(unlist(lapply(LiamTest, function(x) x$fitMK$solution[1,2])), OUEstTab)
colnames(OUEstTab) <- c("mk", "alpha", "sigma","theta1", "theta2", "model")

# hOUwie Wt
hOUWtTab <- do.call(rbind, lapply(LiamTest, function(x) exp(x$res_hOUwie.wt$solution)))
colnames(hOUWtTab) <- c("mk", "alpha", "sigma", "theta1", "theta2")
hOUWtTab <- as.data.frame(hOUWtTab)
hOUWtTab$model <- "hOUWt"

# hOUwie True
hOUTruTab <- do.call(rbind, lapply(LiamTest, function(x) exp(x$res_hOUwie$solution)))
colnames(hOUTruTab) <- c("mk", "alpha", "sigma", "theta1", "theta2")
hOUTruTab <- as.data.frame(hOUTruTab)
hOUTruTab$model <- "hOUTru"

OUMTestTable <- rbind(OUTrueTab, OUEstTab, hOUWtTab, hOUTruTab)

file.name <- paste(format(Sys.time(), "%y.%m.%d"), "_OUM_MK", p.mk, "_TestTable.Rsave", sep = "")
save(OUMTestTable, file = file.name)

# load("~/2020_hOUwie/")
# require(viridis)
# require(ggplot2)
# 
# cols <- viridis(4)
# 
# ggplot(df.Q1_Sig2_1, aes(x=par, y=val, fill=fit)) + 
#   labs(x = "Transition Rate", y = expression('estimated '*sigma[1]^2)) +
#   scale_fill_manual(values=cols) + 
#   theme(text = element_text(size = 20)) + 
#   geom_boxplot()
# 
# ggplot(df.Q1_Sig2_2, aes(x=par, y=val, fill=fit)) + 
#   labs(x = "Transition Rate", y = expression('estimated '*sigma[2]^2)) +
#   scale_fill_manual(values=cols) + 
#   theme(text = element_text(size = 20)) + 
#   geom_boxplot()





