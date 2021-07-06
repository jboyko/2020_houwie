setwd("~/2020_hOUwie/")

source("hOUwieSimmap.R")
source("Utils.R")

require(partitions)

# the 2-state models we want to fit
continuous_models_cd_ou <- getAllContinuousModelStructures(2, "OU")
continuous_models_cd_bm <- getAllContinuousModelStructures(2, "BM")
continuous_models_cd_bmou <- getAllContinuousModelStructures(2, "BMOU")

k <- max(continuous_models_cd_ou[,,8], na.rm = TRUE)
p <- runif(k, 0, 1)

independence_values <- apply(continuous_models_cd_ou, 3, function(x) c(x[1,] - min(x[1,]) + 1, 
                                                                       x[2,] - min(x[2,]) + 1, 
                                                                       x[3,] - min(x[3,]) + 1))

expected_differences <- cbind(
  alpha = abs(independence_values[1,] - independence_values[2,]),
  sigma = abs(independence_values[3,] - independence_values[4,]),
  theta = abs(independence_values[5,] - independence_values[6,]))

p <- sapply(1:1000, function(x) runif(k, 0, 1))
p <- t(lhs:::improvedLHS(1000, 6))
observed_differences <- cbind(alpha = abs(p[1,] - p[2,]), 
                              sigma = abs(p[3,] - p[4,]), 
                              theta = abs(p[5,] - p[6,]))


par(mfrow=c(1,3))
hist(observed_differences[,1])
hist(observed_differences[,2])
hist(observed_differences[,3])

pairs(observed_differences)
