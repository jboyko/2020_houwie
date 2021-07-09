require(expm)

getBranchProb <- function(Q, bl, x0, xt, n){
  if(n == 0){
    path_i <- c(x0, xt)
    llik_vec <- getPathProb(path_i, bl, Q)
    names(llik_vec) <- paste0(path_i, collapse="_")
    return(llik_vec)
  }
  t <- bl/(n+1)
  p_mat <- expm(Q * t)
  combo.list <- list()
  for(i in 1:n){
    combo.list[[i]] <- c(1,2)
  }
  combos <- expand.grid(combo.list)
  llik_vec <- numeric(dim(combos)[1])
  for(i in 1:dim(combos)[1]){
    path_i <- c(x0, as.numeric(combos[i,]), xt)
    llik_vec[i] <- getPathProb(path_i, bl, Q)
    names(llik_vec)[i] <- paste0(path_i, collapse="_")
  }
  return(llik_vec)
}

getPathProb <- function(path, bl, Q){
  section_length <- bl/(length(path)-1)
  p_mat <- expm(Q * section_length)
  path_cp <- path
  P <- vector("numeric", length(path)-1)
  for(i in 1:(length(path)-1)){
    P[i] <- p_mat[path_cp[1],path_cp[2]]
    path_cp <- path_cp[-1]
  }
  return(sum(log(P)))
}

Q <- matrix(c(-1,1,1,-1), 2, 2)
x0 = 1
xt = 1
n = 1
bl = 1

c(getPathProb(c(1,1), bl, Q))

c(getPathProb(c(1,1,1), bl, Q),
  getPathProb(c(1,2,1), bl, Q))

c(getPathProb(c(1,1,1,1), bl, Q),
  getPathProb(c(1,2,1,1), bl, Q),
  getPathProb(c(1,1,2,1), bl, Q),
  getPathProb(c(1,2,2,1), bl, Q))

getBranchProb(Q, bl, x0 = 1, xt = 2, n=2)

sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=0)))
sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=1)))
sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=2)))
sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=3)))

Q = matrix(c(-1,1,.5,-.5), 2, 2)

sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=0)))
sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=1)))
sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=2)))
sum(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=3)))



getBranchProb(Q*100, bl, x0 = 1, xt = 1, n=0)
getBranchProb(Q*100, bl, x0 = 1, xt = 1, n=1)
getBranchProb(Q*100, bl, x0 = 1, xt = 1, n=2)
round(exp(getBranchProb(Q, bl, x0 = 1, xt = 1, n=3))/sum(exp(getBranchProb(Q*100, bl, x0 = 1, xt = 1, n=3))), 3)


