require(POUMM)
require(expm)

getOUProbBranch <- function(tip_value, states, pars, bl){
  # states are read from tip to root states[1] = tip value, states[2] = node value
  times <- c(0, bl/2, bl)
  alphas <- c(pars[1,states])
  sigma2s <- c(pars[2,states])
  thetas <- c(pars[3,states])
  values <- c(thetas[1], thetas)
  # finding the expected value
  tip_weight <- exp(-((alphas[1] * (times[2] - times[1])) + (alphas[2] * (times[3] - times[2]))))
  theta_1_weight <- tip_weight * (exp(alphas[1] * times[2]) - exp(alphas[1] * times[1]))
  theta_2_weight <- tip_weight * (exp(alphas[2] * times[3]) - exp(alphas[2] * times[2]))
  weights <- c(tip_weight, theta_1_weight, theta_2_weight)/sum(c(tip_weight, theta_1_weight, theta_2_weight))
  expected_value <- sum(values * weights)
  # finding the variance 
  var_weight <- exp(-((2 * alphas[1] * (times[2] - times[1])) + (2 * alphas[2] * (times[2] - times[1]))))
  var_1 <- sigma2s[1]/(2*alphas[1]) * (exp(2 * alphas[1] * times[2]) - exp(2 * alphas[1] * times[1]))
  var_2 <- sigma2s[2]/(2*alphas[2]) * (exp(2 * alphas[2] * times[3]) - exp(2 * alphas[2] * times[2]))
  variance <- sum(c(var_1, var_2)) * var_weight
  
  loglik <- dnorm(tip_value, expected_value, sqrt(variance), log = TRUE)
  return(loglik)
}

getJointProbBranch <- function(states, tip_value, pars, bl, P_mat){
  coninuous_prob <- getOUProbBranch(tip_value, states, pars, bl)
  discrete_prob <- log(P_mat[states[1], states[2]])
  joint_prob <- sum(coninuous_prob, discrete_prob)
  return(joint_prob)
}

getJointBranchMatrix <- function(possible_internal, possible_external, tip_value, pars, bl, P_mat){
  cond_matrix <- matrix(NA, dim(P_mat)[1], dim(P_mat)[2])
  for(i in possible_internal){
    for(j in possible_external){
      cond_matrix[i,j] <- getJointProbBranch(c(i,j), tip_value, pars, bl, P_mat)
    }
  }
  colnames(cond_matrix) <- rownames(cond_matrix) <- c(1:max(possible_internal))
  return(cond_matrix)
}

alpha_1 <- 0.1
alpha_2 <- 0.1
sigma2_1 <- 0.5
sigma2_2 <- 0.5
theta_1 <- 10
theta_2 <- 10
tip_value <- 7.5
rate <- 1

n_states <- 2
states <- c(2,1)
bl = 1
pars <- matrix(c(alpha_1, alpha_2, sigma2_1, sigma2_2, theta_1, theta_2), 3, 2, byrow = TRUE, dimnames = c(list(c("alpha", "sigma2", "theta"), c("1", "2"))))
Q <- matrix(c(-rate,rate,rate,-rate), 2, 2)
P_mat <- expm(Q * bl)


coninuous_prob <- getOUProbBranch(tip_value, states, pars, bl)
discrete_prob <- log(P_mat[states[1], states[2]])
joint_prob <- sum(coninuous_prob, discrete_prob)
#  dOU(tip_value, 10, 1, 0.1, 10, sqrt(0.5), log = TRUE)
getJointProbBranch(states, tip_value, pars, bl, P_mat)


possible_internal <- c(1,2)
possible_external <- c(1,2)


# internal nodes are i, external nodes are j
getJointBranchMatrix(possible_internal, possible_external, tip_value, pars, bl, P_mat)




