source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")
#source("../../hOUwie.R")
require(phytools)
require(OUwie)
require(corHMM)
require(parallel)
require(reshape2)
require(ggplot2)
require(gridExtra)

# can we detect hidden states when they're present?
# can we avoid detecting hidden states when they're absent?
# can we detect character dependence when it's present?
# can we avoid detecting character dependence when it's absent?
# how does hOUwie perform when the model is highly complex and outside of a reasonable range?

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under CDBMS 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we detect character dependence when it's present?
# can we avoid detecting hidden states when they're absent?
Models <- c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
nTip <- 100
nMap <- 100

plots.aic <- plots.resid <- list()
for(i in 1:length(Models)){
  Model <- Models[i]
  Rsaves <- getRsaves(Model, nTip, nMap)
  ResidTables <- AICTables <- list()
  for(j in 1:length(Rsaves)){
    load(Rsaves[j])
    AICTables[[j]] <- out$model.fits$model.table$dAICc
    TruePars <- getMapWeightedPars(out$simulating.data, "Root")
    ResidTables[[j]] <- TruePars[,c(2,3,4)] - out$model.fits$model.pars[,c(1,2,3)]
  }
  AICTables <- do.call(rbind, AICTables)
  ResidTables <- do.call(rbind, ResidTables)
  colnames(AICTables) <-  c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
  MeltedAICwtTable <- melt(AICTables)
  MeltedResidTables <- melt(ResidTables)
  plots.aic[[i]] <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
    theme_linedraw() +
    labs(x = "Fitting Model", y = "dAIC for a given dataset", title = paste0("Simulated under ", Model)) +
    theme(legend.position="bottom") +
    scale_fill_brewer() +
    geom_boxplot()
  plots.resid[[i]] <- ggplot(MeltedResidTables, aes(x=variable, y=value)) +
    theme_linedraw() +
    labs(x = "Parameter", y = "Residual Error", title = paste0("Simulated under ", Model)) +
    theme(legend.position="bottom") +
    scale_fill_brewer() +
    geom_boxplot()
}

grid.arrange(plots.aic[[1]],plots.aic[[2]],plots.aic[[3]],plots.aic[[4]],plots.aic[[5]],plots.aic[[6]])
grid.arrange(plots.resid[[1]],plots.resid[[2]],plots.resid[[3]],plots.resid[[4]],plots.resid[[5]],plots.resid[[6]])





Rsaves <- getRsaves("CIDOUM", nTip, nMap)[1:82]
Rsaves <- getRsaves("CIDOUM", nTip, 250)
ResidTables <- AICTables <- list()
for(j in 1:length(Rsaves)){
  load(Rsaves[j])
  #AICTables[[j]] <- getModelTable(out$model.fits$model.fits, "BIC")$dBIC
  AICTables[[j]] <- getModelTable(out$model.fits, "BIC")$dBIC
  # AICTables[[j]] <- out$model.fits$model.table$dAICc
  TruePars <- getMapWeightedPars(out$simulating.data, "Root")
  ResidTables[[j]] <- out$model.fits$model.pars[,c(1,2,3)] - TruePars[,c(2,3,4)]
}
AICTables <- do.call(rbind, AICTables)
ResidTables <- do.call(rbind, ResidTables)
colnames(AICTables) <-  c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
MeltedAICwtTable <- melt(AICTables)
MeltedResidTables <- melt(ResidTables)
ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
  theme_linedraw() +
  labs(x = "Fitting Model", y = "dAIC for a given dataset", title = paste0("Simulated under ", "HYBOUM")) +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()
ggplot(MeltedResidTables, aes(x=variable, y=value)) +
  theme_linedraw() +
  labs(x = "Parameter", y = "Residual Error", title = paste0("Simulated under ", Model)) +
  theme(legend.position="bottom") +
  ylim(-10, 10) + 
  scale_fill_brewer() +
  geom_violin()


# out$model.fits$model.pars
# out$simulating.data$pars
# 
# CDBMS <- getTotalErrorRate(Rsaves, "Raw")
# p.CDBMS <- ggplot(CDBMS, aes(x=Var2, y=value, fill = Var1)) +
#   theme_linedraw() +
#   labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under CDBMS") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()
# 
# tmp <- CDBMS[CDBMS$P_sigma.sq < 20, ]
# ggplot(tmp) + 
#   geom_point(aes(x = E_sigma.sq, y = P_sigma.sq, color = Type), position = position_dodge(width = 0.1))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under CDOUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we detect character dependence when it's present?
# can we avoid detecting hidden states when they're absent?

model <- "CDOUM"
Rsaves <- dir(paste0("~/2020_hOUwie/ModelTesting/", model), full.names = TRUE)
Rsaves <- Rsaves[grep("nMap", Rsaves)]

AICTables <- list()
for(i in 1:length(Rsaves)){
  load(Rsaves[i])
  AICTables[[i]] <- out$model.fits$model.table$AICcwt
}
AICTables <- do.call(rbind, AICTables)
colnames(AICTables) <-  c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
MeltedAICwtTable <- melt(AICTables)

p.CDOUM <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
  theme_linedraw() +
  labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = paste0("Simulated under ", model)) +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

# CDOUM <- getTotalErrorRate(Rsaves, "Raw")
# 
# p.CDOUM <- ggplot(CDOUM, aes(x=Var2, y=value, fill = Var1)) +
#   theme_linedraw() +
#   labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under CDOUM") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under CIDBMS
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we avoid detecting character dependence when it's absent?
# can we detect hidden states when they're present?

model <- "CIDBMS"
Rsaves <- dir(paste0("~/2020_hOUwie/ModelTesting/", model), full.names = TRUE)
Rsaves <- Rsaves[grep("nMap", Rsaves)]

AICTables <- list()
for(i in 1:length(Rsaves)){
  load(Rsaves[i])
  AICTables[[i]] <- out$model.fits$model.table$AICcwt
}
AICTables <- do.call(rbind, AICTables)
colnames(AICTables) <-  c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
MeltedAICwtTable <- melt(AICTables)

p.CIDBMS <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
  theme_linedraw() +
  labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = paste0("Simulated under ", model)) +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

# CIDBMS <- getTotalErrorRate(Rsaves)
# p.CIDBMS <- ggplot(CIDBMS, aes(x=Var2, y=value, fill = Var1)) +
#   theme_linedraw() +
#   labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under CIDBMS") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under CIDOUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we avoid detecting character dependence when it's absent?
# can we detect hidden states when they're present?

model <- "CIDOUM"
Rsaves <- dir(paste0("~/2020_hOUwie/ModelTesting/", model), full.names = TRUE)
Rsaves <- Rsaves[grep("nMap", Rsaves)]

AICTables <- list()
for(i in 1:length(Rsaves)){
  load(Rsaves[i])
  AICTables[[i]] <- out$model.fits$model.table$AICcwt
}
AICTables <- do.call(rbind, AICTables)
colnames(AICTables) <-  c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
MeltedAICwtTable <- melt(AICTables)

p.CIDOUM <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
  theme_linedraw() +
  labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = paste0("Simulated under ", model)) +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

# CIDOUM <- getTotalErrorRate(Rsaves)
# p.CIDOUM <- ggplot(CIDOUM, aes(x=Var2, y=value, fill = Var1)) +
#   theme_linedraw() +
#   labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under CIDOUM") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under HYBBMS
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we detect both the character dependent and independent effects?

model <- "HYBBMS"
Rsaves <- dir(paste0("~/2020_hOUwie/ModelTesting/", model), full.names = TRUE)
Rsaves <- Rsaves[grep("nMap", Rsaves)]

AICTables <- list()
for(i in 1:length(Rsaves)){
  load(Rsaves[i])
  AICTables[[i]] <- out$model.fits$model.table$AICcwt
}
AICTables <- do.call(rbind, AICTables)
colnames(AICTables) <-  c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
MeltedAICwtTable <- melt(AICTables)

p.HYBBMS <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
  theme_linedraw() +
  labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = paste0("Simulated under ", model)) +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

# HYBBMS <- getTotalErrorRate(Rsaves)
# p.HYBBMS <- ggplot(HYBBMS, aes(x=Var2, y=value, fill = Var1)) +
#   theme_linedraw() +
#   labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under HYBBMS") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under HYBOUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we detect both the character dependent and independent effects?

model <- "HYBOUM"
Rsaves <- dir(paste0("~/2020_hOUwie/ModelTesting/", model), full.names = TRUE)
Rsaves <- Rsaves[grep("nMap", Rsaves)]

AICTables <- list()
for(i in 1:length(Rsaves)){
  load(Rsaves[i])
  AICTables[[i]] <- out$model.fits$model.table$AICcwt
}
AICTables <- do.call(rbind, AICTables)
colnames(AICTables) <-  c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
MeltedAICwtTable <- melt(AICTables)

p.HYBOUM <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
  theme_linedraw() +
  labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = paste0("Simulated under ", model)) +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

# HYBOUM <- getTotalErrorRate(Rsaves)
# p.HYBOUM <- ggplot(HYBOUM, aes(x=Var2, y=value, fill = Var1)) +
#   theme_linedraw() +
#   labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under HYBOUM") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()
# 
grid.arrange(p.CDBMS, p.CDOUM, p.CIDBMS, p.CIDOUM, p.HYBBMS, p.HYBOUM)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under HYBOUMVA
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we detect both the character dependent and independent effects?

Rsaves <- dir("~/2020_hOUwie/ModelTesting/HYBOUMVA/", full.names = TRUE)
models <- c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")

# AICTables <- lapply(Rsaves, getAICParamTablefromRsave)
# AICcWtTable <- do.call(rbind, lapply(AICTables, function(x) x$AICcWt))
# colnames(AICcWtTable) <- models
# MeltedAICwtTable <- melt(AICcWtTable)
# 
# p.HYBOUMVA <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
#   theme_linedraw() +
#   labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = "Simulated under HYBOUMVA") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()
# 
# test <- summarize.hOUwie.set(out[[2]])

HYBOUMVA <- getTotalErrorRate(Rsaves)
p.HYBOUMVA <- ggplot(HYBOUMVA, aes(x=Var2, y=value, fill = Var1)) +
  theme_linedraw() +
  labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under HYBOUMVA") +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

