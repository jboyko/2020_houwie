# calculates tip edge stats based step 1 of Ho and Ane's (2014) algorithm
getTipStats <- function(t, y, x){
  logV = log(t)       #aka logd
  p = 1/t             # aka vec11
  uY = y
  uX = x
  Q = (uX %*% uY)/t
  return(list(p = p, uY = uY, uX = uX, Q = Q, logV = logV))
}

getNodeStats <- function(time, DescList){
  # descendent subtree values
  ps <- unlist(lapply(DescList, function(x) x$p))
  logVs <- unlist(lapply(DescList, function(x) x$logV))
  uYs <- unlist(lapply(DescList, function(x) x$uY))
  uXs <- lapply(DescList, function(x) x$uX)
  Qs <- lapply(DescList, function(x) x$Q)
  # combination for node values
  pA <- sum(ps)
  ws <- ps/pA
  logV = sum(logVs) + log(1 + time*pA)
  p = pA/(1 + time*pA)
  uY = sum(ws * uYs)
  uX = rowSums(sapply(1:length(ws), function(x) ws[x] * uXs[[x]]))
  Q = Reduce("+", Qs) - QsSub(time, pA, uX, uY)
  return(list(p = p, uY = uY, uX = uX, Q = Q, logV = logV))
}

QsSub <- function(t, pa, ux, uy){
  obj <- t(((t * pa^2)/(1 + t*pa)) %*% ux) %*% uy
  return(obj)
}

getThreePoint <- function(phy, X, Y){
  generations <- hisse:::FindGenerations(phy)
  nTip <- Ntip(phy)
  
  dec <- phy$edge[,2]
  tip.index <- match(1:nTip, dec)
  tip_times <- phy$edge.length[tip.index]
  TipList <- lapply(1:nTip, function(x) getTipStats(tip_times[x], y = Y[x], x = as.matrix(X[x,])))
  
  DescList <- vector("list", dim(phy$edge)[1])
  for( i in 1:length(tip.index)){
    DescList[[tip.index[i]]] <- TipList[[i]]
  }
  
  for(i in 1:(length(generations)-1)){
    for(j in 1:length(generations[[i]])){
      anc.index <- which(phy$edge[,1] %in% generations[[i]][j])
      dec.index <- which(phy$edge[,2] %in% generations[[i]][j])
      anc.time <- phy$edge.length[dec.index]
      DescList[[dec.index]] <- getNodeStats(anc.time, DescList[anc.index])
    }
  }
  
  root <- generations[[length(generations)]]
  root.index <- which(phy$edge[,1] %in% root)
  res <- getNodeStats(0, DescList[root.index])
  return(res)
}

require(hisse)
require(phylolm)
phy = rtree(10)
Y = rTrait(n=1,phy)
p = rTrait(n=1,phy)
X = cbind(1,p)

getThreePoint(phy, X, Y)
three.point.compute(phy, X, Y)


# T1 <- getTipStats(t[2], Y[1], as.matrix(X[1,]))
# T2 <- getTipStats(t[3], Y[2], as.matrix(X[2,]))
# 
# DescList <- list(T1, T2)
# 
# N1 <- getNodeStats(t[1], DescList)
# T3 <- getTipStats(t[4], Y[3], as.matrix(X[3,]))
# 
# DecList <- list(N1, T3)
# getNodeStats(0, DecList)
# 
# 
