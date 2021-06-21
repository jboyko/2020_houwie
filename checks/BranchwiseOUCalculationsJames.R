# # calculate a contrast 
# getC_i <- function(lambda, x_l, x_r, v_l, v_r){
#   eta_i <- getEta_i(lambda, v_l, v_r)
#   c_i <- (eta_i * exp(-lambda * v_r) * x_l) - (eta_i * exp(-lambda * v_l) * x_r)
#   return(c_i)
# }
# 
# # calculate the intermediate
# getW_i <- function(lambda, x_l, x_r, v_l, v_r){
#   T_i <- getT_i(lambda, v_l, v_r)
#   A_i <- getA_i(lambda, T_i, v_r)
#   B_i <- getA_i(lambda, T_i, v_l)
#   w_i <- (A_i * x_l) + (B_i * x_r)
#   return(w_i)
# }
# 
# # calculate eta_i
# getEta_i <- function(lambda, v_l, v_r){
#   r_value <- exp(-2 * lambda * v_r)
#   l_value <- exp(-2 * lambda * v_l)
#   b_value <- -2 * exp(-2 * lambda * (v_r + v_l))
#   denom <- r_value + l_value + b_value
#   out <- sqrt(2/denom)
#   return(out)
# }
# 
# # get t, a pre-req of a_i
# getT_i <- function(lambda, v_l, v_r){
#   A <- exp(2 * lambda * v_l) + exp(2 * lambda * v_r) - 2
#   B <- 1 - exp(-2 * lambda * (v_r + v_l))
#   out <- 1/sqrt(A * B)
#   return(out)
# }
# 
# # calculate a_i
# getA_i <- function(lambda, t, v){
#   out <- t * (exp(lambda * v) - exp(-lambda * v))
#   return(out)
# }
# 
# # calculate the distance 
# getDeltaD <- function(lambda, v_l, v_r){
#   numer <- exp(2 * lambda * (v_r + v_l)) - 1
#   denom <- exp(2 * lambda * v_l) + exp(2 * lambda * v_r) - 2
#   out <- 1/lambda * log(sqrt(numer/denom))
#   return(out)
# }

getContrasts <- function(phy, dat, alpha, sigma2){
  phy <- reorder(phy, "pruningwise")
  x.anc <- vector("numeric", max(phy$edge))
  names(x.anc) <- 1:max(phy$edge)
  x.anc[1:length(phy$tip.label)] <- dat[phy$tip.label]
  weights <- contrasts <- numeric(length(phy$tip.label) - 1)
  distances <- phy$edge.length
  # for each internal node Wi
  count <- 1
  for(W_i in unique(phy$edge[,1])){
    edge_index <- which(phy$edge[,1] == W_i)
    dec_index <- phy$edge[edge_index,2]
    # define two descendent trait values as xl and xr
    dec_values <- x.anc[dec_index]
    # define two descendent distances from Wi as vl and vr
    dec_distance <- distances[edge_index]
    # calculate the variance of ci 
    variances_i <- varOU(dec_distance, alpha, sqrt(sigma2))
    # calculate contrast ci as a difference between xl and xr
    contrasts[count] <- (dec_values[1] * exp(-alpha * variances_i[1])) - (dec_values[2] * exp(-alpha * variances_i[2]))
    # calculate an intermediate value wi as a weighted average between xl and xr
    weights[count] <- sum(1/variances_i * dec_values)/sum(1/variances_i)
    # assign the new ancestral state value
    x.anc[W_i] <- weights[count]
    # assign a new distance with variance accumulated due to estimating the ancestral state
    D_i <- prod(variances_i)/(sum(variances_i))
    distances[phy$edge[,2] == W_i] <- distances[phy$edge[,2] == W_i] + D_i
    count <- count +1
  }
  return(contrasts)
}


# test the method 
require(geiger)
require(phylolm)
require(POUMM)

phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = 10) 
phy <- drop.extinct(phy)

phy = rtree(10)
alpha = 1
sigma2 = 1
sigma2_error = 0
ancestral.state = 0
optimal.value = 0
parameters = list(ancestral.state=ancestral.state, alpha=alpha,
                  sigma2=sigma2,sigma2_error=sigma2_error,
                  optimal.value=optimal.value)
dat = rTrait(n = 1, phy, model = "OU", parameters = parameters)
PossibleContrasts <- getContrasts(phy, dat, alpha, sigma2)
c("PhylmFixed" = OU1d.loglik(dat, phy, "OUfixedRoot", parameters = parameters),
  "PhylmRandm" = OU1d.loglik(dat, phy, "OUrandomRoot", parameters = parameters),
  "Contrast" = sum(dnorm(PossibleContrasts, 0, sqrt(sigma2/2*alpha), log = TRUE)))




