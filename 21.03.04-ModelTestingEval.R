organizeSimulators <- function(simulators){
  k.cor <- max(simulators$index.cor, na.rm = TRUE) - 1 # number of corhmm params
  k.ou <- max(simulators$index.ou, na.rm = TRUE) - 1 # number of ouwie params
  p.mk <- simulators$pars[1:k.cor]
  p.ou <- simulators$pars[(k.cor+1):length(simulators$pars)]
  simulators$index.cor[simulators$index.cor == 0] <- NA
  simulators$index.cor[is.na(simulators$index.cor)] <- max(simulators$index.cor, na.rm = TRUE) + 1
  Q <- matrix(0, dim(simulators$index.cor)[1], dim(simulators$index.cor)[1])
  Q[] <- c(p.mk, 0)[simulators$index.cor]
  diag(Q) <- -rowSums(Q)
  Rate.mat <- matrix(1, 3, dim(simulators$index.cor)[2])
  simulators$index.ou[is.na(simulators$index.ou)] <- max(simulators$index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[simulators$index.ou]
  rownames(Rate.mat) <- c("alpha", "sigma.sq", "theta")
  return(list(pars.cor = Q,
              pars.ou = Rate.mat))
}

getParDifferences <- function(simiulators, fit){
  
  fit$solution.ou
  simulators$index.ou
  
  k.cor <- max(simulators$index.cor, na.rm = TRUE) - 1 # number of corhmm params
  k.ou <- max(simulators$index.ou, na.rm = TRUE) - 1 # number of ouwie params
  p.mk <- simulators$pars[1:k.cor]
  p.ou <- simulators$pars[(k.cor+1):length(simulators$pars)]
  simulators$index.cor[simulators$index.cor == 0] <- NA
  simulators$index.cor[is.na(simulators$index.cor)] <- max(simulators$index.cor, na.rm = TRUE) + 1
  Q <- matrix(0, dim(simulators$index.cor)[1], dim(simulators$index.cor)[1])
  Q[] <- c(p.mk, 0)[simulators$index.cor]
  diag(Q) <- -rowSums(Q)
  Rate.mat <- matrix(1, 3, dim(simulators$index.cor)[2])
  simulators$index.ou[is.na(simulators$index.ou)] <- max(simulators$index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[simulators$index.ou]
  rownames(Rate.mat) <- c("alpha", "sigma.sq", "theta")
  
}

# general measure of error
# ou and cor specific error
# param specific error

source("~/2020_hOUwie/hOUwie.R")
#source("../../hOUwie.R")
require(OUwie)
require(corHMM)
require(parallel)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ER BMS
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


sub.folders <- dir("~/2020_hOUwie/ModelTesting/ER_BMS", full.names = TRUE)
MSE.tables <- Bias.tables <- Var.tables <- vector("list", length(sub.folders))

par(mfrow=c(2,2))
for(i in 1:length(sub.folders)){
  Rsaves <- dir(sub.folders[i], full.names = TRUE)
  big.obj <- list()
  for(j in 1:length(Rsaves)){
    load(Rsaves[j])
    big.obj[[j]] <- obj
  }
  truepars <- organizeSimulators(obj$data.houwie)
  
  ## sigma squared ratio
  sig2.TwoStep <- do.call(rbind, lapply(big.obj, function(x) x$TwoStepFit$solution.ou[2,]))
  sig2.NonCens <- do.call(rbind, lapply(big.obj, function(x) x$NonCensFit$solution.ou[2,]))
  sig2.hOUwie <- do.call(rbind, lapply(big.obj, function(x) x$hOUwieFit$solution.ou[2,]))
  sig2.TruMap <- do.call(rbind, lapply(big.obj, function(x) x$TruMapFit$solution[2,]))
  
  sig2.ratio.table <- data.frame(TwoStep=sig2.TwoStep[,2]/sig2.TwoStep[,1],
                                 NonCens=sig2.NonCens[,2]/sig2.NonCens[,1],
                                 hOUwie=sig2.hOUwie[,2]/sig2.hOUwie[,1],
                                 TrueMap=sig2.TruMap[,2]/sig2.TruMap[,1])
  sig2.ratio.table <- sig2.ratio.table[!apply(sig2.ratio.table, 1, function(x) any(x > 100)),]
  
  sig2.resid.table <- data.frame(TwoStepSig1=sig2.TwoStep[,1] - truepars$pars.ou[2,1],
                                 TwoStepSig2=sig2.TwoStep[,2] - truepars$pars.ou[2,2],
                                 NonCensSig1=sig2.NonCens[,1] - truepars$pars.ou[2,1],
                                 NonCensSig2=sig2.NonCens[,2] - truepars$pars.ou[2,2],
                                 hOUwieSig1=sig2.hOUwie[,1] - truepars$pars.ou[2,1],
                                 hOUwieSig2=sig2.hOUwie[,2] - truepars$pars.ou[2,2],
                                 TrueMapSig1=sig2.TruMap[,1] - truepars$pars.ou[2,1],
                                 TrueMapSig2=sig2.TruMap[,2] - truepars$pars.ou[2,2])
  
  
  ## mean squared error (MSE)
  ## E(thet_hat - theta)^2
  MSE.tables[[i]] <- data.frame(
    MSE.Sig1 = c(mean((sig2.TwoStep[,1] - truepars$pars.ou[2,1])^2),
                mean((sig2.NonCens[,1] - truepars$pars.ou[2,1])^2),
                mean((sig2.hOUwie[,1] - truepars$pars.ou[2,1])^2),
                mean((sig2.TruMap[,1] - truepars$pars.ou[2,1])^2)),
    MSE.Sig2 = c(mean((sig2.TwoStep[,2] - truepars$pars.ou[2,2])^2),
                mean((sig2.NonCens[,2] - truepars$pars.ou[2,2])^2),
                mean((sig2.hOUwie[,2] - truepars$pars.ou[2,2])^2),
                mean((sig2.TruMap[,2] - truepars$pars.ou[2,2])^2)),
    row.names = c("TwoStep", "NoneCens", "hOUwie", "TrueMap")
  )
  
  ## calculate the bias
  ## E(theta_hat) - theta
  ## E(theta_hat - theta) (i.e. the expectation of the error)
  Bias.tables[[i]] <- data.frame(
    B.Sig1 = c(mean(sig2.TwoStep[,1]) - truepars$pars.ou[2,1],
                mean(sig2.NonCens[,1]) - truepars$pars.ou[2,1],
                mean(sig2.hOUwie[,1]) - truepars$pars.ou[2,1],
                mean(sig2.TruMap[,1]) - truepars$pars.ou[2,1]),
    B.Sig2 = c(mean(sig2.TwoStep[,2]) - truepars$pars.ou[2,2],
                mean(sig2.NonCens[,2]) - truepars$pars.ou[2,2],
                mean(sig2.hOUwie[,2]) - truepars$pars.ou[2,2],
                mean(sig2.TruMap[,2]) - truepars$pars.ou[2,2]),
    row.names = c("TwoStep", "NoneCens", "hOUwie", "TrueMap")
  )
  
  ## calculate the variance 
  Var.tables[[i]] <- data.frame(
    Var.Sig1 = c(var(sig2.TwoStep[,1]),
                var(sig2.NonCens[,1]),
                var(sig2.hOUwie[,1]),
                var(sig2.TruMap[,1])),
    Var.Sig2 = c(var(sig2.TwoStep[,2]),
                var(sig2.NonCens[,2]),
                var(sig2.hOUwie[,2]),
                var(sig2.TruMap[,2])),
    row.names = c("TwoStep", "NoneCens", "hOUwie", "TrueMap")
  )
  
  boxplot(sig2.ratio.table,
          ylab = expression(paste(sigma [2] ^2 / sigma [1] ^2 )),
          xlab = "Model used to estimate 500 datasets",
          ylim = c(0, 40),
          main = paste("q =", truepars$pars.cor[1,2]))
  abline(h = truepars$pars.ou[2,2]/truepars$pars.ou[2,1], col = "red")
}


ErrTable <- rbind(
  cbind(rate=1, cbind(MSE.tables[[1]], Bias.tables[[1]], Var.tables[[1]])),NA,
  cbind(rate=2, cbind(MSE.tables[[2]], Bias.tables[[2]], Var.tables[[2]])),NA,
  cbind(rate=4, cbind(MSE.tables[[3]], Bias.tables[[3]], Var.tables[[3]])),NA,
  cbind(rate=8, cbind(MSE.tables[[4]], Bias.tables[[4]], Var.tables[[4]]))
)

round(ErrTable, 3)


## theta estimate
theta.TwoStep <- do.call(rbind, lapply(big.obj, function(x) x$TwoStepFit$solution.ou[3,]))
theta.NonCens <- do.call(rbind, lapply(big.obj, function(x) x$NonCensFit$solution.ou[3,]))
theta.hOUwie <- do.call(rbind, lapply(big.obj, function(x) x$hOUwieFit$solution.ou[3,]))
theta.TruMap <- do.call(rbind, lapply(big.obj, function(x) x$TruMapFit$solution[3,]))

theta.table <- data.frame(TwoStep=theta.TwoStep[,1],
                          NonCens=theta.NonCens[,1],
                          hOUwie=theta.hOUwie[,1],
                          TrueMap=theta.TruMap[,1])

boxplot(theta.table, 
        ylab = expression(paste(theta [estimate])),
        xlab = "Model used to estimate 500 datasets",
        main = paste("q =", truepars$pars.cor[1,2]))
abline(h = truepars$pars.ou[3,1], col = "red")


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ER OUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

sub.folders <- dir("~/2020_hOUwie/ModelTesting/ER_OUM", full.names = TRUE)
MSE.tables <- Bias.tables <- Var.tables <- vector("list", length(sub.folders))

Rsaves <- dir(sub.folders[i], full.names = TRUE)
big.obj <- list()
for(j in 1:length(Rsaves)){
  load(Rsaves[j])
  big.obj[[j]] <- obj
}
truepars <- organizeSimulators(obj$data.houwie)

## sigma squared ratio
theta.TwoStep <- do.call(rbind, lapply(big.obj, function(x) x$TwoStepFit$solution.ou[3,]))
theta.NonCens <- do.call(rbind, lapply(big.obj, function(x) x$NonCensFit$solution.ou[3,]))
theta.hOUwie <- do.call(rbind, lapply(big.obj, function(x) x$hOUwieFit$solution.ou[3,]))
theta.TruMap <- do.call(rbind, lapply(big.obj, function(x) x$TruMapFit$solution[3,]))

theta.ratio.table <- data.frame(TwoStep=theta.TwoStep[,2]/theta.TwoStep[,1],
                                NonCens=theta.NonCens[,2]/theta.NonCens[,1],
                                hOUwie=theta.hOUwie[,2]/theta.hOUwie[,1],
                                TrueMap=theta.TruMap[,2]/theta.TruMap[,1])

theta.ratio.table <- theta.ratio.table[!apply(theta.ratio.table, 1, function(x) any(x > 100)),]

theta.resid.table <- data.frame(TwoStepTheta1=theta.TwoStep[,1] - truepars$pars.ou[2,1],
                                TwoStepTheta2=theta.TwoStep[,2] - truepars$pars.ou[2,2],
                                NonCensTheta1=theta.NonCens[,1] - truepars$pars.ou[2,1],
                                NonCensTheta2=theta.NonCens[,2] - truepars$pars.ou[2,2],
                                hOUwieTheta1=theta.hOUwie[,1] - truepars$pars.ou[2,1],
                                hOUwieTheta2=theta.hOUwie[,2] - truepars$pars.ou[2,2],
                                TrueMapTheta1=theta.TruMap[,1] - truepars$pars.ou[2,1],
                                TrueMapTheta2=theta.TruMap[,2] - truepars$pars.ou[2,2])

## mean squared error (MSE)
## E(thet_hat - theta)^2
MSE.tables[[i]] <- data.frame(
  MSE.Theta1 = c(mean((theta.TwoStep[,1] - truepars$pars.ou[3,1])^2),
                 mean((theta.NonCens[,1] - truepars$pars.ou[3,1])^2),
                 mean((theta.hOUwie[,1] - truepars$pars.ou[3,1])^2),
                 mean((theta.TruMap[,1] - truepars$pars.ou[3,1])^2)),
  MSE.Theta2 = c(mean((theta.TwoStep[,2] - truepars$pars.ou[3,2])^2),
                 mean((theta.NonCens[,2] - truepars$pars.ou[3,2])^2),
                 mean((theta.hOUwie[,2] - truepars$pars.ou[3,2])^2),
                 mean((theta.TruMap[,2] - truepars$pars.ou[3,2])^2)),
  row.names = c("TwoStep", "NoneCens", "hOUwie", "TrueMap")
)

## calculate the bias
## E(theta_hat) - theta
## E(theta_hat - theta) (i.e. the expectation of the error)
Bias.tables[[i]] <- data.frame(
  B.Theta1 = c(mean(theta.TwoStep[,1]) - truepars$pars.ou[2,1],
               mean(theta.NonCens[,1]) - truepars$pars.ou[2,1],
               mean(theta.hOUwie[,1]) - truepars$pars.ou[2,1],
               mean(theta.TruMap[,1]) - truepars$pars.ou[2,1]),
  B.Theta2 = c(mean(theta.TwoStep[,2]) - truepars$pars.ou[2,2],
               mean(theta.NonCens[,2]) - truepars$pars.ou[2,2],
               mean(theta.hOUwie[,2]) - truepars$pars.ou[2,2],
               mean(theta.TruMap[,2]) - truepars$pars.ou[2,2]),
  row.names = c("TwoStep", "NoneCens", "hOUwie", "TrueMap")
)

## calculate the variance 
Var.tables[[i]] <- data.frame(
  Var.Theta1 = c(var(theta.TwoStep[,1]),
                 var(theta.NonCens[,1]),
                 var(theta.hOUwie[,1]),
                 var(theta.TruMap[,1])),
  Var.Theta2 = c(var(theta.TwoStep[,2]),
                 var(theta.NonCens[,2]),
                 var(theta.hOUwie[,2]),
                 var(theta.TruMap[,2])),
  row.names = c("TwoStep", "NoneCens", "hOUwie", "TrueMap")
)

boxplot(theta.ratio.table,
        ylab = expression(paste(theta [2]/ theta [1])),
        xlab = "Model used to estimate 500 datasets",
        ylim = c(0, 20),
        main = paste("q =", truepars$pars.cor[1,2]))
abline(h = truepars$pars.ou[3,2]/truepars$pars.ou[3,1], col = "red")

