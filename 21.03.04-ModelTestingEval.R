g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

getBigObj <- function(sub.folder){
  Rsaves <- dir(sub.folder, full.names = TRUE)
  big.obj <- list()
  for(j in 1:length(Rsaves)){
    load(Rsaves[j])
    big.obj[[j]] <- obj
  }
  return(big.obj)
}

getTrueParsVector <- function(truepars, index.ou){
  if(any(index.ou[1,] > index.ou[2,])){
    index.ou[1,][index.ou[1,] > index.ou[2,]] <- 0
  }
  pars.count <- c(alpha=0, sigma=0, theta=0)
  out.vec <- c()
  for(i in 1:max(index.ou, na.rm = TRUE)){
    index_i <- index.ou == i
    name.par <- rownames(index_i)[apply(index_i, 1, any)]
    pars.count[apply(index_i, 1, any)] <- pars.count[apply(index_i, 1, any)] + 1
    name.par <- paste(name.par, pars.count[apply(index_i, 1, any)], sep = "_")
    vec_i <- truepars$pars.ou[index_i][1]
    names(vec_i) <- name.par
    out.vec <- c(out.vec, vec_i)
  }
  return(out.vec)
}

getParTable <- function(index.ou, big.obj){
  # determine how many of each parameter we have
  if(any(index.ou[1,] > index.ou[2,])){
    index.ou[1,][index.ou[1,] > index.ou[2,]] <- 0
  }
  pars.count <- c(alpha=0, sigma=0, theta=0)
  table <- c()
  for(i in 1:max(index.ou, na.rm = TRUE)){
    index_i <- index.ou == i
    name.par <- rownames(index_i)[apply(index_i, 1, any)]
    pars.count[apply(index_i, 1, any)] <- pars.count[apply(index_i, 1, any)] + 1
    name.par <- paste(name.par, pars.count[apply(index_i, 1, any)], sep = "_")
    table_i <- rbind(
      data.frame(par = name.par, model = "TwoStepFit", 
                 value = unlist(lapply(big.obj, function(x) 
                   x$TwoStepFit$solution.ou[index_i][1]))),
      data.frame(par = name.par, model = "NonCensFit", 
                 value = unlist(lapply(big.obj, function(x) 
                   x$NonCensFit$solution.ou[index_i][1]))),
      data.frame(par = name.par, model = "hOUwieFit", 
                 value = unlist(lapply(big.obj, function(x) 
                   x$hOUwieFit$solution.ou[index_i][1]))),
      data.frame(par = name.par, model = "TruMapFit", 
                 value = unlist(lapply(big.obj, function(x) 
                   x$TruMapFit$solution[index_i][1])))
    )
    table <- rbind(table, table_i)
  }
  return(table)
}

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
require(ggplot2)
require(gridExtra)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## All models
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# sub.folders <- dir("~/2020_hOUwie/ModelTesting/ER_BMS", full.names = TRUE)

sub.folders <- dir("~/2020_hOUwie/ModelTesting/ER_OUM", full.names = TRUE)

# sub.folders <- dir("~/2020_hOUwie/ModelTesting/ER_OUMV", full.names = TRUE)

sub.folder <- sub.folders[1]
count <- 1
p.resids <- p.RMSE <- list()
for(sub.folder in sub.folders){
  big.obj <- getBigObj(sub.folder)
  truepars <- organizeSimulators(big.obj[[1]]$data.houwie)
  index.ou <- big.obj[[1]]$data.houwie$index.ou
  table.resids <- table.pars <- getParTable(index.ou, big.obj)
  sim.pars <- getTrueParsVector(truepars, index.ou)
  
  for(j in 1:length(sim.pars)){
    table.resids[names(sim.pars)[j] == table.resids[,1], 3] <- table.resids[names(sim.pars)[j] == table.resids[,1], 3] - sim.pars[j]
  }
  table.resids$model <- factor(table.resids$model, levels = c("TwoStepFit", "NonCensFit", "hOUwieFit", "TruMapFit"))
  RMSE <- aggregate(table.resids$value, by=list(table.resids$model, table.resids$par), function(x) sqrt(mean(x^2)))
  SE <- aggregate(table.resids$value, by=list(table.resids$model, table.resids$par), function(x) sd(x)/length(x))
  RMSE.table <- data.frame(Model = RMSE[,1], Param = RMSE[,2], Value = RMSE[,3], SE = SE[,3])
  
  Model <- paste(big.obj[[1]]$hOUwieFit$model.cor, big.obj[[1]]$hOUwieFit$model.ou, sep="+")
  Model <- paste("Model:", Model)
  Pars <- paste(sim.pars, collapse = ", ")
  Pars <- paste("Generating Pars:", Pars)
  Rate <- truepars$pars.cor[2,1]
  Rate <- paste("Mk Rate:", Rate)
  
  pd <- position_dodge(0.1) # move them .05 to the left and right

  # p.RMSE[[count]] <- 
    ggplot(RMSE.table, aes(x=Param, y=Value, colour=Model)) +
    geom_errorbar(aes(ymin=Value-SE, ymax=Value+SE), width=.1) + 
      ylim(-10,10) +
    geom_point()
  
  
  # resid graph
  p.resids[[count]] <- ggplot(table.resids, aes(x=par, y=value, fill=model)) + 
    theme_linedraw() +
    labs(x = "Parameter", y = "RMSE", title = Model, subtitle = Pars, caption = Rate) + 
    theme(legend.position="bottom") +
    ylim(-10,10) +
    scale_fill_brewer() +
    geom_boxplot()
  
  count <- count + 1
}

legend <- g_legend(p.resids[[1]])

grid.arrange(arrangeGrob(p.resids[[1]] + theme(legend.position="none"),
                         p.resids[[2]] + theme(legend.position="none") + labs(title = ""),
                         p.resids[[3]] + theme(legend.position="none") + labs(title = ""),
                         p.resids[[4]] + theme(legend.position="none") + labs(title = ""),
                         nrow=2),
             legend, nrow=2, heights=c(10, 1))


grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]])



MSE.table # mean squared resids
Bias.table # mean resids
Var.table # var of params









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

par(mfrow=c(2,2))
for(i in 1:4){
  
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

theta.resid.table <- data.frame(TwoStepTheta1=theta.TwoStep[,1] - truepars$pars.ou[3,1],
                                TwoStepTheta2=theta.TwoStep[,2] - truepars$pars.ou[3,2],
                                NonCensTheta1=theta.NonCens[,1] - truepars$pars.ou[3,1],
                                NonCensTheta2=theta.NonCens[,2] - truepars$pars.ou[3,2],
                                hOUwieTheta1=theta.hOUwie[,1] - truepars$pars.ou[3,1],
                                hOUwieTheta2=theta.hOUwie[,2] - truepars$pars.ou[3,2],
                                TrueMapTheta1=theta.TruMap[,1] - truepars$pars.ou[3,1],
                                TrueMapTheta2=theta.TruMap[,2] - truepars$pars.ou[3,2])

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
  B.Theta1 = c(mean(theta.TwoStep[,1]) - truepars$pars.ou[3,1],
               mean(theta.NonCens[,1]) - truepars$pars.ou[3,1],
               mean(theta.hOUwie[,1]) - truepars$pars.ou[3,1],
               mean(theta.TruMap[,1]) - truepars$pars.ou[3,1]),
  B.Theta2 = c(mean(theta.TwoStep[,2]) - truepars$pars.ou[3,2],
               mean(theta.NonCens[,2]) - truepars$pars.ou[3,2],
               mean(theta.hOUwie[,2]) - truepars$pars.ou[3,2],
               mean(theta.TruMap[,2]) - truepars$pars.ou[3,2]),
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

}

MSE.tables





