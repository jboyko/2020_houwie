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

Q <- matrix(c(-2,1,2,-1), 2, 2)
time_ij <- 1
state_i <- 1
state_j <- 2
shift.point <- 0.33

pathP <- probPath(Q, state_i, state_i, time_ij, shift.point)
manualP <- dexp(shift.point * time_ij, 2) * (1 - pexp(time_ij - shift.point * time_ij, 2))

#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
# check the probability of a particular path vs all possible paths
#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 

