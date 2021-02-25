require(OUwie)
data(tworegime)
set.seed(42)
Tmax <- max(branching.times(tree))

# BM1 
sig2_A=c(0.4669113, 0.4669113)
sig2_B=sig2_A*Tmax
tree_B <- tree_A <- tree
tree_B$edge.length <- tree_B$edge.length/Tmax

Unscaled <- OUwie.fixed(tree_A, trait, model=c("BM1"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, sigma.sq=sig2_A, shift.point=0.5, algorithm="invert")
Scaled <- OUwie.fixed(tree_B, trait, model=c("BM1"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, sigma.sq=sig2_B, shift.point=0.5, algorithm="invert")

comparison <- identical(round(as.numeric(Unscaled$loglik),3), round(as.numeric(Scaled$loglik),3))
comparison


# BMS 
sig2_A=c(0.2424788, 0.7007112)
sig2_B=sig2_A*Tmax
tree_B <- tree_A <- tree
tree_B$edge.length <- tree_B$edge.length/Tmax

Unscaled <- OUwie.fixed(tree_A, trait, model=c("BMS"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, sigma.sq=sig2_A, shift.point=0.5, algorithm="invert")
Scaled <- OUwie.fixed(tree_B, trait, model=c("BMS"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, sigma.sq=sig2_B, shift.point=0.5, algorithm="invert")

comparison <- identical(round(as.numeric(Unscaled$loglik),3), round(as.numeric(Scaled$loglik),3))
comparison


# OU1
alpha_A=c(0.358939, 0.3589399)
sig2_A=c(0.5197486, 0.5197486)
theta_B=theta_A=c(1.3301447, 1.3301447)
alpha_B=alpha_A*Tmax
sig2_B=sig2_A*Tmax
tree_B <- tree_A <- tree
tree_B$edge.length <- tree_B$edge.length/Tmax

Unscaled <- OUwie.fixed(tree_A, trait, model=c("OU1"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_A, sigma.sq=sig2_A,theta=theta_A, shift.point=0.5, algorithm="invert")
Scaled <- OUwie.fixed(tree_B, trait, model=c("OU1"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_B, sigma.sq=sig2_B,theta=theta_B, shift.point=0.5, algorithm="invert")
comparison <- identical(round(as.numeric(Unscaled$loglik),5), round(as.numeric(Scaled$loglik),5))
comparison


# OUM
alpha_A=c(1.3916589, 1.3916589)
sig2_A=c(0.6545502, 0.6545502)
theta_B=theta_A=c(1.6751330, 0.4424138)
alpha_B=alpha_A*Tmax
sig2_B=sig2_A*Tmax
tree_B <- tree_A <- tree
tree_B$edge.length <- tree_B$edge.length/Tmax

Unscaled <- OUwie.fixed(tree_A, trait, model=c("OUM"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_A, sigma.sq=sig2_A,theta=theta_A, shift.point=0.5, algorithm="invert")
Scaled <- OUwie.fixed(tree_B, trait, model=c("OUM"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_B, sigma.sq=sig2_B,theta=theta_B, shift.point=0.5, algorithm="invert")
comparison <- identical(round(as.numeric(Unscaled$loglik),5), round(as.numeric(Scaled$loglik),5))
comparison


# OUMV
alpha_A=c(1.7110818, 1.711082)
sig2_A=c(0.3517019, 1.076479)
theta_B=theta_A=c(1.676894, 0.5563541)
alpha_B=alpha_A*Tmax
sig2_B=sig2_A*Tmax
tree_B <- tree_A <- tree
tree_B$edge.length <- tree_B$edge.length/Tmax

Unscaled <- OUwie.fixed(tree_A, trait, model=c("OUMV"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_A, sigma.sq=sig2_A,theta=theta_A, shift.point=0.5, algorithm="invert")
Scaled <- OUwie.fixed(tree_B, trait, model=c("OUMV"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_B, sigma.sq=sig2_B,theta=theta_B, shift.point=0.5, algorithm="invert")
comparison <- identical(round(as.numeric(Unscaled$loglik),5), round(as.numeric(Scaled$loglik),5))
comparison


# OUMA
alpha_A=c(1.6501816, 1.0294487)
sig2_A=c(0.7082462, 0.7082462)
theta_B=theta_A=c(1.6765718, 0.1516105)
alpha_B=alpha_A*Tmax
sig2_B=sig2_A*Tmax
tree_B <- tree_A <- tree
tree_B$edge.length <- tree_B$edge.length/Tmax

Unscaled <- OUwie.fixed(tree_A, trait, model=c("OUMA"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_A, sigma.sq=sig2_A,theta=theta_A, shift.point=0.5, algorithm="invert")
Scaled <- OUwie.fixed(tree_B, trait, model=c("OUMA"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_B, sigma.sq=sig2_B,theta=theta_B, shift.point=0.5, algorithm="invert")
comparison <- identical(round(as.numeric(Unscaled$loglik),5), round(as.numeric(Scaled$loglik),5))
comparison


# OUMVA
alpha_A=c(3.0793193, 0.6060786)
sig2_A=c(0.4735485, 1.7049102)
theta_B=theta_A=c(1.68189033, -1.032546)
alpha_B=alpha_A*Tmax
sig2_B=sig2_A*Tmax
tree_B <- tree_A <- tree
tree_B$edge.length <- tree_B$edge.length/Tmax

Unscaled <- OUwie.fixed(tree_A, trait, model=c("OUMA"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_A, sigma.sq=sig2_A,theta=theta_A, shift.point=0.5, algorithm="invert")
Scaled <- OUwie.fixed(tree_B, trait, model=c("OUMA"), simmap.tree=FALSE, scaleHeight=FALSE, clade=NULL, alpha=alpha_B, sigma.sq=sig2_B,theta=theta_B, shift.point=0.5, algorithm="invert")
comparison <- identical(round(as.numeric(Unscaled$loglik),5), round(as.numeric(Scaled$loglik),5))
comparison
