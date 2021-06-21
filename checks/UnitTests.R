# comparison of BMS and BM1
source("~/2020_hOUwie/hOUwieEM.R")
source("~/2020_hOUwie/Utils.R")
require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)

#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
# check the probability of a particular path vs all possible paths
#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 

Q <- matrix(c(-1,1,1,-1), 2, 2)
time_ij <- 1
state_i <- 1
state_j <- 2
shift.point <- 0.5

path <- c(0.75, 1)
names(path) <- c(1, 2)

probPathOld(path, Q)
probPath(path, Q)


manualP <- dexp(0.5, 1) * (1 - pexp(0.5 * 1, 2))

#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
# check the probability of a particular path vs all possible paths
#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 

# Browse[6]> path
# 1         2         4 
# 0.4574998 0.4574998 0.4574998 
# Browse[6]> print(Q)
# (1,R1)      (2,R1)      (1,R2)      (2,R2)
# (1,R1) -0.03341306  0.01670653  0.01670653  0.00000000
# (2,R1)  0.01670653 -0.03341306  0.00000000  0.01670653
# (1,R2)  0.01670653  0.00000000 -0.03341306  0.01670653
# (2,R2)  0.00000000  0.01670653  0.01670653 -0.03341306



