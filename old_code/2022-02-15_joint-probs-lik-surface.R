setwd("~/2020_houwie/")
require(geiger)
require(parallel)
require(corHMM)
require(OUwie)
require(partitions)
require(expm)
require(MASS)
require(phytools)
require(ggtree)
require(ggplot2)
require(reshape2)
require(gridExtra)
require(ggplotify)
source("hOUwieNode.R")

# balenced_tree <- read.tree(text = "((S1:0.5,S2:0.5):0.5,(S3:0.5,S4:0.5):0.5);")
phy <- read.tree(text = "(((S1:0.33, S2:0.33):0.33,S3:0.66):0.34,S4:1);")
plot(phy)

phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = 20) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

#### #### #### #### #### #### #### #### #### 
# OUM CD
#### #### #### #### #### #### #### #### #### 
# generating parameters
p = c(.1,1,.1,10,20)
Q = matrix(c(-p[1],p[1],p[1],-p[1]),2,2)
# consistent dataset with equal distribution of discrete characters
set.seed(51444)
dat_cd <- hOUwie.sim(phy, Q, c(1,0), p[c(2,2)], p[c(3,3)], p[4], p[4:5])
plot(dat_cd$simmap)
# fit under true parameters
res <- hOUwie(phy, dat_cd$data, 1, "ER", "OUM", 1.1, 20, p = p, sample_nodes = FALSE)
# examining the likelihood surface
# free_pars <- matrix(corHMM:::GetLHSPoints(0.1, 20, 500), 250, 2)
free_pars <- as.data.frame(expand.grid(list(1:11, 6:25)))
joint_res <- apply(free_pars, 1, function(x) hOUwie(phy, dat_cd$data, 1, "ER", "OUM", 1.1, 20, p = c(p[1:3], x), sample_nodes = FALSE))
joint_lliks <- lapply(joint_res, function(x) x$loglik)
llik_mapping <- data.frame(theta_1 = free_pars[,1], theta_2 = free_pars[,2], joint_lliks = (unlist(joint_lliks)))
plot_theta <- ggplot(llik_mapping, aes(theta_1, theta_2, z = joint_lliks)) + 
  geom_contour_filled() +
  geom_vline(xintercept=10, linetype="dashed", alpha = 0.4, colour = "red") +
  geom_hline(yintercept=20, linetype="dashed", alpha = 0.4, colour = "red")



#### #### #### #### #### #### #### #### #### 
# OUV CD
#### #### #### #### #### #### #### #### #### 
p = c(.1,1,.1,10,10)
Q = matrix(c(-p[1],p[1],p[1],-p[1]),2,2)
# consistent dataset with equal distribution of discrete characters
set.seed(1992)
dat_cd <- hOUwie.sim(phy, Q, c(1,0), p[c(2,2)], p[c(3,4)], p[4], p[c(5,5)])
plot(dat_cd$simmap)
# fit under true parameters
cont_model <- getAllContinuousModelStructures(2)[,,3]
res <- hOUwie(phy, dat_cd$data, 1, "ER", cont_model, 1.1, 20, p = p, sample_nodes = FALSE)
# examining the likelihood surface
# free_pars <- matrix(corHMM:::GetLHSPoints(0.1, 20, 500), 250, 2)
free_pars <- as.data.frame(expand.grid(list(seq(from = 0.01, to = 1, length.out = 20), seq(from = 1, to = 20, length.out = 20))))
joint_res <- apply(free_pars, 1, function(x) hOUwie(phy, dat_cd$data, 1, "ER", cont_model, 1.1, 20, p = c(p[1:2], x, p[5]), sample_nodes = FALSE))
joint_lliks <- lapply(joint_res, function(x) x$loglik)
llik_mapping <- data.frame(sigma_1 = free_pars[,1], sigma_2 = free_pars[,2], joint_lliks = log(-unlist(joint_lliks)))
plot_sigma <- ggplot(llik_mapping, aes(sigma_1, sigma_2, z = joint_lliks)) + 
  geom_contour_filled() +
  geom_vline(xintercept=.1, linetype="dashed", alpha = 0.4, colour = "red") +
  geom_hline(yintercept=10, linetype="dashed", alpha = 0.4, colour = "red")


#### #### #### #### #### #### #### #### #### 
# OUA CD
#### #### #### #### #### #### #### #### #### 
p = c(.1,.1,10,.1,10)
Q = matrix(c(-p[1],p[1],p[1],-p[1]),2,2)
# consistent dataset with equal distribution of discrete characters
# set.seed(1992)
dat_cd <- hOUwie.sim(phy, Q, c(1,0), p[c(2,3)], p[c(4,4)], p[5], p[c(5,5)])
plot(dat_cd$simmap)
# fit under true parameters
cont_model <- getAllContinuousModelStructures(2)[,,2]
res <- hOUwie(phy, dat_cd$data, 1, "ER", cont_model, 1.1, 200, p = p, sample_nodes = FALSE)
# examining the likelihood surface
# free_pars <- matrix(corHMM:::GetLHSPoints(0.1, 20, 500), 250, 2)
free_pars <- as.data.frame(expand.grid(list(seq(from = 0.01, to = 1, length.out = 20), seq(from = 1, to = 20, length.out = 20))))
joint_res <- apply(free_pars, 1, function(x) hOUwie(phy, dat_cd$data, 1, "ER", cont_model, 1.1, 20, p = c(p[1], x, p[4:5]), sample_nodes = FALSE))
joint_lliks <- lapply(joint_res, function(x) x$loglik)
llik_mapping <- data.frame(alpha_1 = free_pars[,1], alpha_2 = free_pars[,2], joint_lliks = unlist(joint_lliks))
plot_alpha <- ggplot(llik_mapping, aes(alpha_1, alpha_2, z = joint_lliks)) + 
  geom_contour_filled() +
  geom_vline(xintercept=.1, linetype="dashed", alpha = 0.4, colour = "red") +
  geom_hline(yintercept=10, linetype="dashed", alpha = 0.4, colour = "red")

head(llik_mapping[sort(llik_mapping$joint_lliks, decreasing = TRUE, index=TRUE)$ix,])
grid.arrange(plot_theta, plot_sigma, plot_alpha)

