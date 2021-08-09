source("~/2020_hOUwie/hOUwieNode.R")
source("~/2020_hOUwie/Utils.R")
require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)

#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
# check that the probability of all mappings is the same as corhmm 
#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = 5) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
Q <- matrix(c(-1,1,1,-1),2,2)
data <- data.frame(sp = phy$tip.label, reg = c(1,1,1,2,2), x = c(5,5,5,10,10))
houwie_disc <- getAllJointProbs(phy, data, 1, 1.1, Q, c(4,4), c(1,1), c(5,10))
corhmm_disc <- corHMM(phy, data[,c(1,2)], 1, model = "ER", p = 1, root.p = c(0.5,0.5))
round(log(sum(exp(houwie_disc[,1]))),5) == round(corhmm_disc$loglik,5)

#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 
# next test goes here
#### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### #### #### ### ### 

