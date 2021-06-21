## Comparison of different causal model structures in phylogenetic comparative methods

## Load packages
library(bayou)
library(phylolm)

## We have three traits. One focal trait Y is an effect of both X and Z. However, the 
## graph of causal relationships between them differ, as well as whether there are 
## phylogenetic effects. 

## Define true causal parameters: 
Beta_XY <- runif(1, -1, 1)
Beta_ZY <- runif(1, -1, 1)
Beta_ZX <-Beta_XZ <- runif(1, -1, 1)
 
rse <- 0.5
ntips <- 500
fitLMs <- function(x, y, z){
  res <- list()
  #OLS xz
  res$ols.xz <- summary(lm(y~x+z))
  #OLS just x
  res$ols.x <- summary(lm(y~x))
  #PGLS xz
  res$pgls.xz <- summary(phylolm(y~x+z, phy=tree, model="BM"))
  #PGLS x
  res$pgls.x <- summary(phylolm(y~x, phy=tree, model="BM"))
  return(res)
}

## Case Set I. X -> Y <- Z
### Case IA. Phylogenetic effects of X. 
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- rnorm(ntips, 0, 1)
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z + rnorm(ntips, 0, rse)
.res <- fitLMs(.X, .Y, .Z)
.res$ols.xz
.res$pgls.xz

### Case IB. Phylogenetic effects in X and Z
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z + rnorm(ntips, 0, rse)
pairs(data.frame(Y=.Y, X=.X, Z=.Z))
.res <- fitLMs(.X, .Y, .Z)
.res

### Case IC. Phylogenetic effects in X, Z and Y
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z + sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
pairs(data.frame(Y=.Y, X=.X, Z=.Z))
.res <- fitLMs(.X, .Y, .Z)
.res

### Case ID. Phylogenetic effects in Z
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- rnorm(ntips, 0, 1) #sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z# + rnorm(ntips, 0, rse) #sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
.df <-data.frame(Y=.Y, X=.X, Z=.Z) 
#apply(.df, 2, phylosig, tree=tree, method="lambda")
pairs(.df)
.res <- fitLMs(.X, .Y, .Z)
.res
Beta_XY

### Case ID. Phylogenetic effects in Y
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- rnorm(ntips, 0, 1) #sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- rnorm(ntips, 0, 1) #sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z + sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
.df <-data.frame(Y=.Y, X=.X, Z=.Z) 
apply(.df, 2, phylosig, tree=tree, method="lambda")
pairs(.df)
.res <- fitLMs(.X, .Y, .Z)
.res


## Case II. X -> Y <- Z <- X
### Case IIA. Phylogenetic effects of X. 
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- .X * Beta_XZ + rnorm(ntips, 0, rse) #sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z + rnorm(ntips, 0, rse) #sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
pairs(data.frame(Y=.Y, X=.X, Z=.Z))
.res <- fitLMs(.X, .Y, .Z)
.res

### Case IIB Phylogenetic signal in X and Y
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- .X * Beta_XZ + rnorm(ntips, 0, rse) #sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z + sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
pairs(data.frame(Y=.Y, X=.X, Z=.Z))
.res <- fitLMs(.X, .Y, .Z)
.res

### Case IIC Phylogenetic signal in X, Z
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- .X * Beta_XZ + sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1] #sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z + rnorm(ntips, 0, rse) #sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
pairs(data.frame(Y=.Y, X=.X, Z=.Z))
.res <- fitLMs(.X, .Y, .Z)
.res

### Case IID Phylogenetic signal in X, Y
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- .X * Beta_XZ + rnorm(ntips, 0, rse) #sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z + sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
pairs(data.frame(Y=.Y, X=.X, Z=.Z))
.res <- fitLMs(.X, .Y, .Z)
.res

### Case IIE Phylogenetic Signal in Z
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.X <- rnorm(ntips, 0, rse) #sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.X <- scale(.X)
.Z <- .X * Beta_XZ + sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.Y <- Beta_XY*.X + Beta_ZY*.Z + rnorm(ntips, 0, rse) #sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
pairs(data.frame(Y=.Y, X=.X, Z=.Z))
.res <- fitLMs(.X, .Y, .Z)
.res



### Case III Z -> X -> Y <- Z
### Case IIIA. Phylogenetic signal in Z
tree <- sim.bdtree(stop="taxa", n=ntips)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
.Z <- sim.char(tree, par=1, model="BM", nsim=1, root=0)[,,1]
.Z <- scale(.Z)
.X <- .Z * Beta_ZX + rnorm(ntips, 0, rse)
.X <- scale(.X)
.Y <- Beta_XY*.X + Beta_ZY*.Z + sim.char(tree, par=rse^2, model="BM", nsim=1, root=0)[,,1]
pairs(data.frame(Y=.Y, X=.X, Z=.Z))
.res <- fitLMs(.X, .Y, .Z)
.res


### On the Felsenstein Tree:
ntips=100
rse=0.5
felsTree <- function(n, stemL=1, lambda=1){
  phy1 <- reorder(sim.bdtree(b=1, d=0, stop="taxa", n=n), "postorder")
  phy2 <- reorder(sim.bdtree(b=1, d=0, stop="taxa", n=n), "postorder")
  phy1 <- rescale(phy1, model="lambda", lambda)
  phy2 <- rescale(phy1, model="lambda", lambda)
  phy1$edge.length <- phy1$edge.length/max(branching.times(phy1))
  phy2$edge.length <- phy2$edge.length/max(branching.times(phy2))
  o <- list(phy1=phy1, phy2=phy2)
  phy2$tip.label <- paste("s", (n+1):(2*n), sep="")
  phy <- list(edge=NULL, Nnode=(n-1)*2 + 1, tip.label=c(phy1$tip.label, phy2$tip.label), edge.length=NULL)
  phy1$edge[phy1$edge > n] <- phy1$edge[phy1$edge > n]+(n+1)
  phy2$edge[phy2$edge > n] <- phy2$edge[phy2$edge > n]+2*n
  phy2$edge[phy2$edge <= n] <- phy2$edge[phy2$edge <=n]+n
  
  phy$edge <- rbind(phy1$edge, phy2$edge, c((2*n+1), (2*n+n+1)), c((2*n+1), (2*n+2)))
  phy$edge.length <- c(phy1$edge.length, phy2$edge.length, rep(stemL, 2))
  class(phy) <- class(o$phy1)
  phy <- reorder(phy, "postorder")
  plot(phy)
  return(phy)
}
OLS1 <- OLS2 <- PGLS1 <- PGLS2 <- NULL
for(i in 1:2000){
  phy <- felsTree(ntips/2, stemL=9, lambda=0)
  phy$edge.length <- phy$edge.length/(max(branching.times(phy)))

  ## BS(X) -> MR(Y) <- BT(Z) <- BS(X)
  X <- sim.char(phy, par=1, model="BM", nsim=1, root=0)[,,1]
  Z <- rnorm(ntips, X, rse)
  X <- X + rnorm(ntips, 0, rse)
  Z <- ifelse(Z>median(Z), 1, 0)
  Y <- rnorm(ntips, 1*X+-3*Z, rse)
  df1 <- data.frame(Y=Y, X=X, Z=Z)
  #pairs(df1)

  ## BT(Z) -> BS(X) -> MR(Y) <- BT(Z)
  Z <- sim.char(phy, par=1, model="BM", nsim=1, root=0)[,,1]
  X <- rnorm(ntips, Z, rse)
  Z <- ifelse(Z>median(Z), 1, 0)
  plot(phy, show.tip.label=FALSE)
  tiplabels(pch=21, bg=Z)
  Y <- rnorm(ntips, 1*X+-3*Z, rse)
  df2 <- data.frame(Y=Y, X=X, Z=Z)
  #pairs(df2)

  pgls2 <- phylolm(df2$Y~df2$X, phy=phy, model="BM")
  pgls1 <- phylolm(df1$Y~df1$X, phy=phy, model="BM")
  ols2 <- lm(df2$Y~df2$X)
  ols1 <- lm(df1$Y~df1$X)
  OLS1 <- rbind(OLS1, ols1$coef)
  OLS2 <- rbind(OLS2, ols2$coef)
  PGLS1 <- rbind(PGLS1, pgls1$coef)
  PGLS2 <- rbind(PGLS2, pgls2$coef)

  par(mfrow=c(1,2))
  plot(df1$X, df1$Y)
  abline(ols1, col="blue", lty=2)
  abline(pgls1, col="red", lty=2)
  ols1$coef; pgls1$coef

  plot(df2$X, df2$Y)
  abline(ols2, col="blue", lty=2)
  abline(pgls2, col="red", lty=2)
  ols2$coef; pgls2$coef
  apply(df1, 2, sd)
  apply(df2, 2, sd)
}

hist(OLS1[,2], breaks=seq(-2,2,0.2))
hist(OLS2[,2], add=TRUE, col="red", breaks=seq(-2,2,0.2))

hist(PGLS1[,2], breaks=seq(-2,2,0.2))
hist(PGLS2[,2], add=TRUE, col="red", breaks=seq(-2,2,0.2))
par(mfrow=c(1,2))
plot(OLS1[,2], PGLS1[,2])
plot(OLS2[,2], PGLS2[,2])

quantile(OLS1[,2], c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(OLS2[,2], c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(PGLS1[,2], c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(PGLS2[,2], c(0.025, 0.25, 0.5, 0.75, 0.975))


ntips=40
phy <- felsTree(ntips/2, stemL=5, lambda=0)
par(mfrow=c(2,2))
phy$edge.length <- phy$edge.length/(max(branching.times(phy)))
set.seed(6)
Z <- sim.char(phy, par=1, model="BM", nsim=1, root=0)[,,1]
X <- rnorm(ntips, Z, rse)
Z <- ifelse(Z>median(Z), 1, 0)
Y <- rnorm(ntips, -1*X+3*Z, rse)
df1 <- data.frame(Y=Y, X=X, Z=Z)

png("PGLSOLS2.png", width=900, height=575)
par(mfrow=c(1,2))
par(mar=c(2.5,4,3,2))
plot(phy, show.tip.label=FALSE, type="fan")
text(-1,1.25, label="C.", cex=2)
tiplabels(pch=21, bg=df1$Z, cex=2)
legend(-1.05, -1.025, legend=c("Migratory", "Non-migratory"), pch=21, pt.bg=c(1,0), cex=1.5)
par(mar=c(5,6,4,2))
plot(df1$X, df1$Y, pch=21, bg=df1$Z, xlab="Body Size (B)", ylab="Species Abundance (N)", cex.lab=2, cex=1.2)
text(-1, 3.5, label="D.", cex=2)
pgls1 <- phylolm(df1$Y~df1$X, phy=phy, model="BM")
ols1 <- lm(df1$Y~df1$X)
legend(-1.2, -.85, legend=c("OLS", "PGLS"), lty=c(2,3), col=c("black", "gray75"), cex=1.5, lwd=3)

abline(ols1, col="black", lty=2, lwd=3)
abline(pgls1, col="gray75", lty=3, lwd=3)
dev.off()



