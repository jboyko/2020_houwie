# this script will analyze whether hOUwie can preform better than independent corHMM > stoch Map > OUwie

require(viridis)
require(ggplot2)

setwd("~/2020_hOUwie/OUMTestResults//")

load("21.01.15-OUMTestTable.R")

ggplot(OUMTestTable, aes(x=model, y=alpha)) + 
  labs(x = "Transition Rate", y = expression('estimated '*alpha)) +
  theme(text = element_text(size = 20)) + 
  geom_boxplot() +
  geom_hline(yintercept = 2)

ggplot(OUMTestTable, aes(x=model, y=sigma)) + 
  labs(x = "Transition Rate", y = expression('estimated '*sigma^2)) +
  theme(text = element_text(size = 20)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 2)

ggplot(OUMTestTable, aes(x=model, y=theta1)) + 
  labs(x = "Transition Rate", y = expression('estimated '*theta[1])) +
  theme(text = element_text(size = 20)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 3)

ggplot(OUMTestTable, aes(x=model, y=theta2)) + 
  labs(x = "Transition Rate", y = expression('estimated '*theta[2])) +
  theme(text = element_text(size = 20)) + 
  geom_boxplot()



jpeg(filename = "~/2020_hOUwie/doc/hOUwieVsStoch.jpg")
par(mfrow=c(2,2))
ggplot(df.Sig2_ratio, aes(x=par, y=val, fill=fit)) + 
  labs(x = "Transition Rate", y = expression('estimated '*sigma[2]^2/sigma[1]^2)) +
  scale_fill_manual(values=cols) + 
  theme(text = element_text(size = 20)) + 
  geom_boxplot()
dev.off()



