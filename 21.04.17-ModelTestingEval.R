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

Rsaves <- dir("~/2020_hOUwie/ModelTesting/CDBMS/", full.names = TRUE)

# AICTables <- AICTables[unlist(lapply(AICTables, length))!=0]
# AICcWtTable <- do.call(rbind, lapply(AICTables, function(x) x$AICcwt))
# colnames(AICcWtTable) <- models
# MeltedAICwtTable <- melt(AICcWtTable)


# p.CDBMS <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
#   theme_linedraw() +
#   labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = "Simulated under CDBMS") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()

CDBMS <- getTotalErrorRate(Rsaves, "Raw")
p.CDBMS <- ggplot(CDBMS, aes(x=Var2, y=value, fill = Var1)) +
  theme_linedraw() +
  labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under CDBMS") +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

tmp <- CDBMS[CDBMS$P_sigma.sq < 20, ]
ggplot(tmp) + 
  geom_point(aes(x = E_sigma.sq, y = P_sigma.sq, color = Type), position = position_dodge(width = 0.1))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under CDOUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we detect character dependence when it's present?
# can we avoid detecting hidden states when they're absent?

Rsaves <- dir("~/2020_hOUwie/ModelTesting/CDOUM/", full.names = TRUE)
models <- c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")

AICTables <- lapply(Rsaves, getAICParamTablefromRsave)
AICcWtTable <- do.call(rbind, lapply(AICTables, function(x) x$AICcWt))
colnames(AICcWtTable) <- models
MeltedAICwtTable <- melt(AICcWtTable)

# p.CDOUM <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
#   theme_linedraw() +
#   labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = "Simulated under CDOUM") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()

CDOUM <- getTotalErrorRate(Rsaves, "Raw")

p.CDOUM <- ggplot(CDOUM, aes(x=Var2, y=value, fill = Var1)) +
  theme_linedraw() +
  labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under CDOUM") +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under CIDBMS
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we avoid detecting character dependence when it's absent?
# can we detect hidden states when they're present?

Rsaves <- dir("~/2020_hOUwie/ModelTesting/CIDBMS/", full.names = TRUE)
# Rsaves <- Rsaves[-grep("-Summary-", Rsaves)]
models <- c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")

AICTables <- lapply(Rsaves, getAICParamTablefromRsave)
unlist(lapply(AICTables, function(x) sum(x$AICcWt)))
AICcWtTable <- do.call(rbind, lapply(AICTables, function(x) x$AICcWt))
colnames(AICcWtTable) <- models
MeltedAICwtTable <- melt(AICcWtTable)

# p.CIDBMS <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
#   theme_linedraw() +
#   labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = "Simulated under CIDBMS") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()

CIDBMS <- getTotalErrorRate(Rsaves)
p.CIDBMS <- ggplot(CIDBMS, aes(x=Var2, y=value, fill = Var1)) +
  theme_linedraw() +
  labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under CIDBMS") +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under CIDOUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we avoid detecting character dependence when it's absent?
# can we detect hidden states when they're present?

Rsaves <- dir("~/2020_hOUwie/ModelTesting/CIDOUM/", full.names = TRUE)
models <- c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")

AICTables <- lapply(Rsaves, getAICParamTablefromRsave)
AICcWtTable <- do.call(rbind, lapply(AICTables, function(x) x$AICcWt))
colnames(AICcWtTable) <- models
MeltedAICwtTable <- melt(AICcWtTable)

# p.CIDOUM <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
#   theme_linedraw() +
#   labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = "Simulated under CIDOUM") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()

CIDOUM <- getTotalErrorRate(Rsaves)
p.CIDOUM <- ggplot(CIDOUM, aes(x=Var2, y=value, fill = Var1)) +
  theme_linedraw() +
  labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under CIDOUM") +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under HYBBMS
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we detect both the character dependent and independent effects?

Rsaves <- dir("~/2020_hOUwie/ModelTesting/HYBBMS/", full.names = TRUE)
models <- c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")

AICTables <- lapply(Rsaves, getAICParamTablefromRsave)
AICcWtTable <- do.call(rbind, lapply(AICTables, function(x) x$AICcWt))
colnames(AICcWtTable) <- models
MeltedAICwtTable <- melt(AICcWtTable)

# p.HYBBMS <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
#   theme_linedraw() +
#   labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = "Simulated under HYBBMS") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()

HYBBMS <- getTotalErrorRate(Rsaves)
p.HYBBMS <- ggplot(HYBBMS, aes(x=Var2, y=value, fill = Var1)) +
  theme_linedraw() +
  labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under HYBBMS") +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Simulate under HYBOUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# can we detect both the character dependent and independent effects?

Rsaves <- dir("~/2020_hOUwie/ModelTesting/HYBOUM/", full.names = TRUE)
models <- c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")

AICTables <- lapply(Rsaves, getAICParamTablefromRsave)
AICcWtTable <- do.call(rbind, lapply(AICTables, function(x) x$AICcWt))
colnames(AICcWtTable) <- models
MeltedAICwtTable <- melt(AICcWtTable)

# p.HYBOUM <- ggplot(MeltedAICwtTable, aes(x=Var2, y=value)) +
#   theme_linedraw() +
#   labs(x = "Fitting Model", y = "AICcWt for a given dataset", title = "Simulated under HYBOUM") +
#   theme(legend.position="bottom") +
#   scale_fill_brewer() +
#   geom_boxplot()

HYBOUM <- getTotalErrorRate(Rsaves)
p.HYBOUM <- ggplot(HYBOUM, aes(x=Var2, y=value, fill = Var1)) +
  theme_linedraw() +
  labs(x = "Parameter", y = "Mean Residual Error", title = "Simulated under HYBOUM") +
  theme(legend.position="bottom") +
  scale_fill_brewer() +
  geom_boxplot()

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

