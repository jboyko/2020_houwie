# general measure of error
# ou and cor specific error
# param specific error

source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")
#source("../../hOUwie.R")
require(phytools)
require(OUwie)
require(corHMM)
require(parallel)
require(ggplot2)
require(gridExtra)

# models to fit 
# fit Mk BMS, fit Mk OUM
# fit CD BMS, fit CD OUM
# fit CID BMS, fit CID OUM

# can we detect hidden states when they're present?
# can we avoid detecting hidden states when they're absent?
# can we detect character dependence when it's present?
# can we avoid detecting character dependence when it's absent?
# how does hOUwie perform when the model is highly complex and outside of a reasonable range?

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## BMS 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## BMS Mk rate variation

sub.folders <- dir("~/2020_hOUwie/ModelSetTesting/ER_BMS/100", full.names = TRUE)
sub.folders <- dir("~/2020_hOUwie/ModelSetTesting/ER_BMS/250", full.names = TRUE)

big.obj <- getBigObj(sub.folders)
p.RMSE <- plotRSME(big.obj, "resid")

legend <- g_legend(p.RMSE[[1]] + theme(legend.position = "bottom"))
grid.arrange(
  arrangeGrob(
    p.RMSE[[1]] + theme(legend.position="none"),
    p.RMSE[[2]] + theme(legend.position="none") + labs(title = ""),
    p.RMSE[[3]] + theme(legend.position="none") + labs(title = ""),
    #p.RMSE[[4]] + theme(legend.position="none") + labs(title = ""),
    nrow=1),
  legend, nrow=2, heights=c(10, 1)
)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## OUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## OUM alpha:sigma ratio
sub.folders <- dir("~/2020_hOUwie/ModelSetTesting/ER_OUM/100/", full.names = TRUE)[c(2,3,5)]
#sub.folders <- dir("~/2020_hOUwie/ModelSetTesting/ER_OUM/250/", full.names = TRUE)

big.obj <- getBigObj(sub.folders)
p.RMSE <- plotRSME(big.obj, "RMSE")
legend <- g_legend(p.RMSE[[1]] + theme(legend.position = "bottom"))
grid.arrange(
  arrangeGrob(
    p.RMSE[[1]] + theme(legend.position="none") + ylim(c(0,30)),
    p.RMSE[[2]] + theme(legend.position="none") + labs(title = "") + ylim(c(0,30)),
    p.RMSE[[3]] + theme(legend.position="none") + labs(title = "") + ylim(c(0,30)),
    #p.RMSE[[4]] + theme(legend.position="none") + labs(title = ""),
    nrow=1),
  legend, nrow=2, heights=c(10, 1)
)

# closer examinations
big.obj <- getBigObj(sub.folders)[[1]]
TrueTheta1 <- big.obj[[1]]$data.houwie$pars[4]
TrueTheta2 <- big.obj[[1]]$data.houwie$pars[5]

E.Theta1 <- unlist(lapply(big.obj, function(x) x$TruMapFit$solution[3,1]))
E.Theta2 <- unlist(lapply(big.obj, function(x) x$TruMapFit$solution[3,2]))

ASR <- unlist(lapply(big.obj, function(x) 
  length(x$data.houwie$simmap[[1]]$maps[[1]])[1]))

big.obj[[37]]$data.houwie$simmap[[1]]$maps[[1]]
plotSimmap(big.obj[[28]]$data.houwie$simmap[[1]])

PErrTheta1 = abs(E.Theta1 - TrueTheta1)/TrueTheta1 * 100
PErrTheta2 = abs(E.Theta2 - TrueTheta2)/TrueTheta2 * 100

Model <- c(rep("%ErrorrTheta_1", length(PErrTheta1)), rep("%ErrorrTheta_2", length(PErrTheta2)))
dat <- data.frame(RootState = as.factor(ASR), Model = Model, PercErr = c(PErrTheta1, PErrTheta2))
MeanError <- aggregate(dat$PercErr, by=list(dat$RootState, dat$Model), function(x) mean(x))

ggplot(MeanError, aes(x=Group.1, y=x, fill=Group.2)) +
  labs(x = "Observed Root State", y = "Percent Error (%)") +
  theme_linedraw() +
  geom_point(shape=21, size = 4) 


ggplot(dat, aes(x=RootState, y=PercErr, fill=Model)) +
  labs(x = "Observed Root State", y = "Percent Error (%)") +
  theme_linedraw() +
  geom_boxplot()


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CID HMM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
sub.folders <- dir("~/2020_hOUwie/ModelSetTesting/HMM_OUM/500/", full.names = TRUE)

big.obj <- getBigObj2(sub.folders)
nTip <- dim(big.obj[[1]][[1]]$simulating.data$data)[1]

big.obj[[1]][[1]]$model.fits[[1]]$index.cor
big.obj[[1]][[1]]$model.fits[[2]]$index.cor
big.obj[[1]][[1]]$model.fits[[3]]$index.cor
big.obj[[1]][[1]]$model.fits[[4]]$index.cor

big.obj[[1]][[1]]$model.fits[[1]]$index.ou
big.obj[[1]][[1]]$model.fits[[2]]$index.ou
big.obj[[1]][[1]]$model.fits[[3]]$index.ou
big.obj[[1]][[1]]$model.fits[[4]]$index.ou

model.fits <- lapply(big.obj[[1]], function(x) x$model.fits)
loglik <- do.call(rbind, (lapply(model.fits, function(x) simplify2array(lapply(x, "[[", "loglik")))))
par.count <- do.call(rbind, (lapply(model.fits, function(x) simplify2array(lapply(x, "[[", "param.count")))))
AIC <- do.call(rbind, (lapply(model.fits, function(x) simplify2array(lapply(x, "[[", "AIC")))))
AICc <- do.call(rbind, (lapply(model.fits, function(x) simplify2array(lapply(x, "[[", "AICc")))))

boxplot(AICc)

AICcTable <- do.call(rbind, lapply(model.fits, function(x) unlist(lapply(x, function(y) y$AICc))))


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## parsimony
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
Rsaves <- dir("~/2020_hOUwie/ModelSetTesting/parsimony/", full.names = TRUE)
big.obj <- list()
for(j in 1:length(Rsaves)){
  load(Rsaves[j])
  big.obj[[j]] <- obj
}

nTip <- dim(big.obj[[1]][[1]]$data)[1]
p.RMSE <- plotRSME(big.obj, 1, "RMSE")


big.obj[[1]]$hOUwieFitPars
big.obj[[1]]$hOUwieFit
big.obj[[1]]$TruMapFit

model.fits <- lapply(big.obj, function(x) x$hOUwieFitPars)
loglik <- do.call(rbind, (lapply(model.fits, function(x) simplify2array(lapply(x, "[[", "loglik")))))
par.count <- do.call(rbind, (lapply(model.fits, function(x) simplify2array(lapply(x, "[[", "param.count")))))
AIC <- do.call(rbind, (lapply(model.fits, function(x) simplify2array(lapply(x, "[[", "AIC")))))
AICc <- do.call(rbind, (lapply(model.fits, function(x) simplify2array(lapply(x, "[[", "AICc")))))

boxplot(AICc)

AICcTable <- do.call(rbind, lapply(model.fits, function(x) unlist(lapply(x, function(y) y$AICc))))


