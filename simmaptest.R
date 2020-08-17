require(corHMM)
require(OUwie)

data(tworegime)

#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
plot(tree)
nodelabels(pch=21, bg=select.reg)

dat <- data.frame(sp = tree$tip.label, d = trait[,2])
MK <- corHMM(phy = tree, rate.cat = 1, dat)
phy <- MK$phy
data <- MK$data
model <- MK$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)
states <- MK$states
tip.states <- MK$tip.states
phy$node.label <- NULL
## run get simmap (can be plotted using phytools)
simmap <- makeSimmap(tree=phy, tip.states=tip.states, states=states, model=model, 
                     nSim=1, nCores=1)

OURes <- OUwie(simmap[[1]], trait, model=c("OUM"), simmap.tree = TRUE)
undebug(OUwie:::weight.mat)

debug(OUwie)


