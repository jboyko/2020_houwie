require(partitions)

getAllContinuousModelStructures <- function(k, type = "OU"){
  # index.mat <- matrix(0, 3, k, dimnames = list(c("alpha", "sigma.sq", "theta"), c(1:k)))
  # we want all unique combinations of a parameter. then we can add a single all same
  # how many combinations are there of 1:k numbers? 
  potential_combos <- apply(partitions:::setparts(k), 2, function(x) paste(x, collapse="_"))
  # this technically isn't all the possible alpha combinations, but for sim purposes we're fine.
  if(type == "OU"){
    alpha.combos <- potential_combos
  }
  if(type == "BM"){
    alpha.combos <- paste(rep(0, k), collapse="_")
  }
  if(type == "BMOU"){
    needed_numerals <- 1:((2^k)-2)
    alpha.combos <- apply(sapply(needed_numerals, function(x) as.numeric(intToBits(x)[1:k])), 2, function(x) paste(x, collapse="_")) # currently doesn't allow for BM mixed with OUA
  }
  sigma.sq.combos <- potential_combos
  theta.combos <- potential_combos
  all_combos <- expand.grid(list(alpha.combos, sigma.sq.combos, theta.combos))
  index_mats <- array(NA, c(3, k, dim(all_combos)[1]), dimnames = list(c("alpha", "sigma.sq", "theta"), c(1:k)))
  for(i in 1:dim(all_combos)[1]){
    alpha_i <- as.numeric(unlist(strsplit(as.character(all_combos[i,1]), "_")))
    alpha_i[alpha_i == 0] <- NA
    sigma_i <- max(c(0, alpha_i), na.rm = TRUE) + as.numeric(unlist(strsplit(as.character(all_combos[i,2]), "_")))
    theta_i <- max(sigma_i) + as.numeric(unlist(strsplit(as.character(all_combos[i,3]), "_")))
    index_mats[,,i] <- rbind(alpha_i, sigma_i, theta_i)
  }
  return(index_mats)
}

# the 2-state models we want to fit
continuous_models_cd_ou <- getAllContinuousModelStructures(2, "OU")
continuous_models_cd_bm <- getAllContinuousModelStructures(2, "BM")
continuous_models_cd_bmou <- getAllContinuousModelStructures(2, "BMOU")


independence_values <- apply(continuous_models_cd_ou, 3, function(x) c(x[1,] - min(x[1,]), 
                                                                       x[2,] - min(x[2,]), 
                                                                       x[3,] - min(x[3,])))

dist(t(independence_values))


# expected_differences <- cbind(
#   alpha = abs(independence_values[1,] - independence_values[2,]),
#   sigma = abs(independence_values[3,] - independence_values[4,]),
#   theta = abs(independence_values[5,] - independence_values[6,]))
# 
# p <- sapply(1:1000, function(x) runif(k, 0, 1))
# p <- t(lhs:::improvedLHS(1000, 6))
# observed_differences <- cbind(alpha = abs(p[1,] - p[2,]), 
#                               sigma = abs(p[3,] - p[4,]), 
#                               theta = abs(p[5,] - p[6,]))
# 
# 
# par(mfrow=c(1,3))
# hist(observed_differences[,1])
# hist(observed_differences[,2])
# hist(observed_differences[,3])
# 
# pairs(observed_differences)
