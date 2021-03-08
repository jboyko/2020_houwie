# get a phylogeny rescaled to a height of 1
source("~/2020_hOUwie/hOUwie.R")
#source("hOUwie.R")
require(OUwie)
require(corHMM)
require(parallel)


setwd("~/2020_hOUwie/likConverge/")
Rsaves <- dir("~/2020_hOUwie/likConverge/")
rates <- unique(unlist(lapply(strsplit(Rsaves, "-"), function(x) x[3])))
maps <- sort(unique(unlist(lapply(strsplit(Rsaves, "-"), function(x) x[2]))))[c(3,2,5,1,4)]


i = j = 1
fin <- c()
par(mfrow=c(1,3))
for(i in 1:length(rates)){
  ToLoad <- Rsaves[grep(rates[i], Rsaves)]
  dat <- NA
  for(j in 1:length(maps)){
    load(ToLoad[grep(maps[j], ToLoad)])
    # dat <- cbind(dat, unlist(lapply(out, function(x) x$loglik)))
    dat <- cbind(dat, unlist(out))
  }
  dat <- dat[,-1]
  colnames(dat) <- gsub("simmaps", "", maps)
  fin <- rbind(fin, dat)
  boxplot(dat, main = rates[i], xlab = "no. of simmaps", ylab = "llik")
}


boxplot(fin)
dat
