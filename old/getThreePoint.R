# calculates tip edge stats based step 1 of Ho and Ane's (2014) algorithm
getTipStats <- function(t, y, x){
  logV = log(t)
  p = 1/t
  uY = y
  uX = x
  Qxy = (uX %*% y)/t
  Qyy = (uY %*% y)/t
  Qxx = (uX %*% t(x))/t
  return(list(p = p, uY = uY, uX = uX,
              Qxy = Qxy, Qyy = Qyy, Qxx = Qxx, 
              logV = logV))
}

# combines descendent values
getNodeStats <- function(time, DescList){
  # descendent subtree values
  ps <- unlist(lapply(DescList, function(x) x$p))
  logVs <- unlist(lapply(DescList, function(x) x$logV))
  uYs <- unlist(lapply(DescList, function(x) x$uY)) # currently univariate, could be changed here
  uXs <- lapply(DescList, function(x) as.matrix(x$uX))
  Qxy <- lapply(DescList, function(x) as.matrix(x$Qxy))
  Qyy <- lapply(DescList, function(x) as.matrix(x$Qyy))
  Qxx <- lapply(DescList, function(x) as.matrix(x$Qxx))
  # combination for node values
  pA <- sum(ps)
  ws <- ps/pA
  logV = sum(logVs) + log(1 + time*pA)
  p = pA/(1 + time*pA)
  uY = sum(ws * uYs) 
  if(dim(uXs[[1]])[1] > 1){
    uX = rowSums(sapply(1:length(ws), function(x) ws[x] * uXs[[x]]))
  }else{
    uX = sum(ws * unlist(uXs))
  }
  Qxy = Reduce("+", Qxy) - QsSub(time, pA, t(uX), uY)
  Qxx = Reduce("+", Qxx) - QsSub(time, pA, t(uX), t(uX))
  Qyy = Reduce("+", Qyy) - QsSub(time, pA, uY, uY)
  return(list(p = p, uY = uY, uX = uX, 
              Qxy = Qxy, Qyy = Qyy, Qxx = Qxx,
              logV = logV))
}

QsSub <- function(t, pa, ux, uy){
  obj <- t(((t * pa^2)/(1 + t*pa)) %*% ux) %*% uy
  return(obj)
}

getThreePoint <- function(phy, X, Y){
  generations <- corHMM:::FindGenerations(phy)
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

# testing
# 
# require(hisse)
# require(phylolm)
# 
# phy = sim.bdtree(b = 1, d = 0, stop = "taxa", n = 10)
# Y = rTrait(n=1,phy)
# p = rTrait(n=1,phy)
# X = cbind(p)
# 
# getThreePoint(phy, X, Y)
# three.point.compute(phy, X, Y)
# 
