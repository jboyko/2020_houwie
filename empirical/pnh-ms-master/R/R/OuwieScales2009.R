library(OUwie)

td <- reorder(td, "postorder")
simpars <- bayou::priorSim(prior.phryno, td$phy, nsim=1)
simpars$pars[[1]]$alpha <- 0.15
simpars$pars[[1]]$sig2 <- 0.001
simpars$pars[[1]]$theta <- c(0.6, 0.3)
simdat <- bayou::dataSim(simpars$pars[[1]], model="OU", td$phy)

ouwie.PE <- bayou2OUwie(H.PE, td$phy, simdat$dat)
ouwie.FM <- bayou2OUwie(H.FM, td$phy, simdat$dat)
ouwie.FMPE <- bayou2OUwie(H.FMPE, td$phy, simdat$dat)
ouwie.phryno <- bayou2OUwie(simpars$pars[[1]], td$phy, simdat$dat)

ouwFit.PE <- OUwie(ouwie.PE$tree, ouwie.PE$dat, model="OUM")
ouwFit.FM <- OUwie(ouwie.FM$tree, ouwie.FM$dat, model="OUM")
ouwFit.FMPE <- OUwie(ouwie.FMPE$tree, ouwie.FMPE$dat, model="OUM")
ouwFit.phryno <- OUwie(ouwie.phryno$tree, ouwie.phryno$dat, model="OUM")


ouwFit.PE$AICc
ouwFit.FM$AICc
ouwFit.FMPE$AICc
ouwFit.phryno$AICc