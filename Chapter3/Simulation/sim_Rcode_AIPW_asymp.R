
# Doubly robust inference for non-probability surveys
# author: Ali Rafei
# E-mail: arafei@umich.edu
# Date  : 12-06-2021
# Paper : https://arxiv.org/abs/2101.07456

# This set of codes replicate the Simulation study by Chen et al (2019)

#rm(list=ls())
setwd("C:\\Users\\arafei\\Desktop\\Thesis\\Simulation\\Chapter III\\Chen_sim\\Frequentist")

library(survey)
library(MASS)

source("sim_func_chen_asymp_IPSW.R") # Read the main DR adjustment function

set.seed(01012021)

sigma=1; # The variance of the ouctome
nSim <- 5000 # Number of simulation iteration
N <- 1000000 # population size
nR <- 100 # 
nA <- 100

# Initial auxiliary varialles
Z1 <- rbinom(N, 1, 0.5)
Z2 <- runif(N, 0, 2)
Z3 <- rexp(N, 1)
Z4 <- rchisq(N, 4)
pop <- data.frame(Z1=Z1, Z2=Z2, Z3=Z3, Z4=Z4) # Create the population
pop$Z5 <- log(pop$Z4)

# Derive the final auxiliary variales
pop$X1 <- pop$Z1
pop$X2 <- pop$Z2+0.3*pop$Z1
pop$X3 <- pop$Z3+0.2*(pop$X1+pop$X2)
pop$X4 <- pop$Z4+0.1*(pop$X1+pop$X2+pop$X3)
pop$X5 <- log(pop$Z4)*(pop$X1+pop$X2+pop$X3)
pop$e_y <- rnorm(N, 0, 1)

# Function finding sigma for a given rho
rtcor <- function(x, e, rho){
uniroot(function(t){cor(x+t*e, x)-rho}, interval=c(0, 20))$root
}

# Create ouctomes with rho = c(0.2, 0.5, 0.8)
sigma <- rtcor(2+pop$X1+pop$X2+pop$X3+pop$X4, pop$e_y, rho=0.2)
pop$Y1 <- 2+pop$X1+pop$X2+pop$X3+pop$X4+sigma*pop$e_y
sigma <- rtcor(2+pop$X1+pop$X2+pop$X3+pop$X4, pop$e_y, rho=0.5)
pop$Y2 <- 2+pop$X1+pop$X2+pop$X3+pop$X4+sigma*pop$e_y
sigma <- rtcor(2+pop$X1+pop$X2+pop$X3+pop$X4, pop$e_y, rho=0.8)
pop$Y3 <- 2+pop$X1+pop$X2+pop$X3+pop$X4+sigma*pop$e_y

#pop$Y2 <- rbinom(N, 1, exp(-8+pop$X1+pop$X2+pop$X3+pop$X4)/(1+exp(-8+pop$X1+pop$X2+pop$X3+pop$X4)))

cor(pop$Y3, 2+pop$X1+pop$X2+pop$X3+pop$X4)
summary(pop$Y1)
summary(pop$Y2)
summary(pop$Y3)

# Function finding intercept for the propensity model in S_R
rtrat <- function(x){
uniroot(function(t){max(t+x)/min(t+x)-50}, interval=c(0, 20))$root
}

c <- rtrat(pop$X3)
pop$Z0 <- c+pop$X3
max(pop$Z0)/min(pop$Z0)

# Generate piR, the selection probabilities associated with S_R
pop$piR <- nR*pop$Z0/sum(pop$Z0) 
sum(pop$piR)

# Function finding intercept for the propensity model in S_A
rtexp1 <- function(x){
uniroot(function(t){mean(exp(t+x)/(1+exp(t+x)))-nA/N}, interval=c(-20, 20))$root
}

gamma1 <- rtexp1(0.1*pop$X1+0.2*pop$X2+0.1*pop$X3+0.2*pop$X4)

# Generate piA, the selection probabilities associated with S_A
pop$piA <- exp(gamma1+0.1*pop$X1+0.2*pop$X2+0.1*pop$X3+0.2*pop$X4)/(1+exp(gamma1+0.1*pop$X1+0.2*pop$X2+0.1*pop$X3+0.2*pop$X4))
sum(pop$piR)
sum(pop$piA)

# Generate and store the selection indicators in S_R and S_A across simulation iterations
ZZ_R <- list()
ZZ_A <- list()
for(i in 1:nSim){
	ZZ_R[[i]] <- (1:N)[as.logical(rbinom(N, 1, pop$piR))]
	ZZ_A[[i]] <- (1:N)[as.logical(rbinom(N, 1, pop$piA))]
}

#write.csv(pop, paste0("sim_pop", ".csv"), col.names=T, sep=",")
#pop <- read.csv(paste0("sim_pop", ".csv"), header=T, sep=",")
#pop$X <- NULL

# Store the true FPMs
mu1 <- mean(pop$Y1)
mu2 <- mean(pop$Y2)
mu3 <- mean(pop$Y2)

#saveRDS(ZZ_R, file=paste0("smp_idx0_", nR, "_", ".rds"))
#saveRDS(ZZ_A, file=paste0("smp_idx1_", nA, "_", ".rds"))
#ZZ_R <- readRDS(file=paste0("smp_idx0_", nR, "_", ".rds"))
#ZZ_A <- readRDS(file=paste0("smp_idx1_", nA, "_", ".rds"))


###############################################################

# Set number of parameters
b <- 6
p <- 4+4*b

# res saves the outputs for each outcomes across the iterations
res1 <- as.data.frame(matrix(NA, nSim, 7*p))
res2 <- as.data.frame(matrix(NA, nSim, 7*p))
res3 <- as.data.frame(matrix(NA, nSim, 7*p))

# restructure the res contents
names(res1) <- c(paste0("mu", 1:p), paste0("se", 1:p), paste0("bias", 1:p), paste0("rmse", 1:p), paste0("ci_ll", 1:p), paste0("ci_ul", 1:p), paste0("cov", 1:p))
names(res2) <- c(paste0("mu", 1:p), paste0("se", 1:p), paste0("bias", 1:p), paste0("rmse", 1:p), paste0("ci_ll", 1:p), paste0("ci_ul", 1:p), paste0("cov", 1:p))
names(res3) <- c(paste0("mu", 1:p), paste0("se", 1:p), paste0("bias", 1:p), paste0("rmse", 1:p), paste0("ci_ll", 1:p), paste0("ci_ul", 1:p), paste0("cov", 1:p))

pcor <- wcor <- matrix(NA, nSim, b-1)

for(i in 1:nSim){
  
  #
	s_R <- pop[ZZ_R[[i]], ]
	s_A <- pop[ZZ_A[[i]], ]

	s_R$wght <- 1/s_R$piR
	s_A$wght <- 1/s_A$piR
	s_A$wght0 <- 1/s_A$piA
	s_R$z <- 0
	s_A$z <- 1
	smp <- rbind(s_A[, c("X1", "X2", "X3", "X4", "Z1", "Z2", "Z3", "Z4", "X5", "z", "wght")], s_R[, c("X1", "X2", "X3", "X4", "Z1", "Z2", "Z3", "Z4", "X5", "z", "wght")])

	# Run the AIPW_KH() function for different outcomes across different combination of model specification PS: TRUE, FALSE; PM: TRUE FALSE
	out10 <- AIPW_KH(y=c(s_A$Y1, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(T, T), family="gaussian")
	out11 <- AIPW_KH(y=c(s_A$Y1, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(T, F), family="gaussian")
	out12 <- AIPW_KH(y=c(s_A$Y1, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(F, T), family="gaussian")
	out13 <- AIPW_KH(y=c(s_A$Y1, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(F, F), family="gaussian")	
	out20 <- AIPW_KH(y=c(s_A$Y2, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(T, T), family="gaussian")
	out21 <- AIPW_KH(y=c(s_A$Y2, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(T, F), family="gaussian")
	out22 <- AIPW_KH(y=c(s_A$Y2, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(F, T), family="gaussian")
	out23 <- AIPW_KH(y=c(s_A$Y2, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(F, F), family="gaussian")	
	out30 <- AIPW_KH(y=c(s_A$Y3, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(T, T), family="gaussian")
	out31 <- AIPW_KH(y=c(s_A$Y3, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(T, F), family="gaussian")
	out32 <- AIPW_KH(y=c(s_A$Y3, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(F, T), family="gaussian")
	out33 <- AIPW_KH(y=c(s_A$Y3, rep(NA, nrow(s_R))), x=smp, deltaB=smp$z, sw=smp$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), spec=c(F, F), family="gaussian")	

	# Estimating unweighted mean and 95%CI	
	res1$mu1[i] <- mean(s_R$Y1, na.rm=T)
	res1$se1[i] <- sd(s_R$Y1, na.rm=T)/sqrt(nR)
	res1$ci_ll1[i] <- res1$mu1[i] - qnorm(0.975)*res1$se1[i]
	res1$ci_ul1[i] <- res1$mu1[i] + qnorm(0.975)*res1$se1[i]

	# Estimating fully weighted mean and 95%CI
	res1$mu2[i] <- weighted.mean(s_R$Y1, w=s_R$wght, na.rm=T)
	res1$se2[i] <- weighted.se(s_R$Y1, s_R$wght, na.rm=T)
	res1$ci_ll2[i] <- res1$mu2[i] - qnorm(0.975)*res1$se2[i]
	res1$ci_ul2[i] <- res1$mu2[i] + qnorm(0.975)*res1$se2[i]

	# Estimating unweighted mean and 95%CI	
	res1$mu3[i] <- mean(s_A$Y1, na.rm=T)
	res1$se3[i] <- sd(s_A$Y1, na.rm=T)/sqrt(nA)
	res1$ci_ll3[i] <- res1$mu3[i] - qnorm(0.975)*res1$se3[i]
	res1$ci_ul3[i] <- res1$mu3[i] + qnorm(0.975)*res1$se3[i]

	# Estimating fully weighted mean and 95%CI
	res1$mu4[i] <- weighted.mean(s_A$Y1, w=s_A$wght0, na.rm=T)
	res1$se4[i] <- weighted.se(s_A$Y1, s_A$wght0, na.rm=T)
	res1$ci_ll4[i] <- res1$mu4[i] - qnorm(0.975)*res1$se4[i]
	res1$ci_ul4[i] <- res1$mu4[i] + qnorm(0.975)*res1$se4[i]

	res1[i, paste0("mu", 5:p)] <- c(out10$Mu, out11$Mu, out12$Mu, out13$Mu)
	res1[i, paste0("se", 5:p)] <- c(out10$SE, out11$SE, out12$SE, out13$SE)

	# Estimating unweighted mean and 95%CI	
	res2$mu1[i] <- mean(s_R$Y2, na.rm=T)
	res2$se1[i] <- sd(s_R$Y2, na.rm=T)/sqrt(nR)
	res2$ci_ll1[i] <- res2$mu1[i] - qnorm(0.975)*res2$se1[i]
	res2$ci_ul1[i] <- res2$mu1[i] + qnorm(0.975)*res2$se1[i]

	# Estimating fully weighted mean and 95%CI
	res2$mu2[i] <- weighted.mean(s_R$Y2, w=s_R$wght, na.rm=T)
	res2$se2[i] <- weighted.se(s_R$Y2, s_R$wght, na.rm=T)
	res2$ci_ll2[i] <- res2$mu2[i] - qnorm(0.975)*res2$se2[i]
	res2$ci_ul2[i] <- res2$mu2[i] + qnorm(0.975)*res2$se2[i]

	# Estimating unweighted mean and 95%CI	
	res2$mu3[i] <- mean(s_A$Y2, na.rm=T)
	res2$se3[i] <- sd(s_A$Y2, na.rm=T)/sqrt(nA)
	res2$ci_ll3[i] <- res2$mu3[i] - qnorm(0.975)*res2$se3[i]
	res2$ci_ul3[i] <- res2$mu3[i] + qnorm(0.975)*res2$se3[i]

	# Estimating fully weighted mean and 95%CI
	res2$mu4[i] <- weighted.mean(s_A$Y2, w=s_A$wght0, na.rm=T)
	res2$se4[i] <- weighted.se(s_A$Y2, s_A$wght0, na.rm=T)
	res2$ci_ll4[i] <- res2$mu4[i] - qnorm(0.975)*res2$se4[i]
	res2$ci_ul4[i] <- res2$mu4[i] + qnorm(0.975)*res2$se4[i]

	res2[i, paste0("mu", 5:p)] <- c(out20$Mu, out21$Mu, out22$Mu, out23$Mu)
	res2[i, paste0("se", 5:p)] <- c(out20$SE, out21$SE, out22$SE, out23$SE)

	# Estimating unweighted mean and 95%CI	
	res3$mu1[i] <- mean(s_R$Y3, na.rm=T)
	res3$se1[i] <- sd(s_R$Y3, na.rm=T)/sqrt(nR)
	res3$ci_ll1[i] <- res3$mu1[i] - qnorm(0.975)*res3$se1[i]
	res3$ci_ul1[i] <- res3$mu1[i] + qnorm(0.975)*res3$se1[i]

	# Estimating fully weighted mean and 95%CI
	res3$mu2[i] <- weighted.mean(s_R$Y3, w=s_R$wght, na.rm=T)
	res3$se2[i] <- weighted.se(s_R$Y3, s_R$wght, na.rm=T)
	res3$ci_ll2[i] <- res3$mu2[i] - qnorm(0.975)*res3$se2[i]
	res3$ci_ul2[i] <- res3$mu2[i] + qnorm(0.975)*res3$se2[i]

	# Estimating unweighted mean and 95%CI	
	res3$mu3[i] <- mean(s_A$Y3, na.rm=T)
	res3$se3[i] <- sd(s_A$Y3, na.rm=T)/sqrt(nA)
	res3$ci_ll3[i] <- res3$mu3[i] - qnorm(0.975)*res3$se3[i]
	res3$ci_ul3[i] <- res3$mu3[i] + qnorm(0.975)*res3$se3[i]

	# Estimating fully weighted mean and 95%CI
	res3$mu4[i] <- weighted.mean(s_A$Y3, w=s_A$wght0, na.rm=T)
	res3$se4[i] <- weighted.se(s_A$Y3, s_A$wght0, na.rm=T)
	res3$ci_ll4[i] <- res3$mu4[i] - qnorm(0.975)*res3$se4[i]
	res3$ci_ul4[i] <- res3$mu4[i] + qnorm(0.975)*res3$se4[i]

	res3[i, paste0("mu", 5:p)] <- c(out30$Mu, out31$Mu, out32$Mu, out33$Mu)
	res3[i, paste0("se", 5:p)] <- c(out30$SE, out31$SE, out32$SE, out33$SE)

	wcor[i, ] <- t(cor(s_A$wght0, out10$RW))
	pcor[i, ] <- t(cor(s_A$piA, 1/out10$RW))
}

res1[, paste0("bias", 1:p)] <- res1[, paste0("mu", 1:p)]-mu1
res1[, paste0("rmse", 1:p)] <- res1[, paste0("bias", 1:p)]^2
res1[, paste0("ci_ll", 5:p)] <- res1[, paste0("mu", 5:p)]-qt(0.975, nA-10)*res1[, paste0("se", 5:p)]
res1[, paste0("ci_ul", 5:p)] <- res1[, paste0("mu", 5:p)]+qt(0.975, nA-10)*res1[, paste0("se", 5:p)]
res1[, paste0("cov", 1:p)] <- res1[, paste0("ci_ll", 1:p)]<mu1 & mu1<res1[, paste0("ci_ul", 1:p)]

res2[, paste0("bias", 1:p)] <- res2[, paste0("mu", 1:p)]-mu2
res2[, paste0("rmse", 1:p)] <- res2[, paste0("bias", 1:p)]^2
res2[, paste0("ci_ll", 5:p)] <- res2[, paste0("mu", 5:p)]-qt(0.975, nA-10)*res2[, paste0("se", 5:p)]
res2[, paste0("ci_ul", 5:p)] <- res2[, paste0("mu", 5:p)]+qt(0.975, nA-10)*res2[, paste0("se", 5:p)]
res2[, paste0("cov", 1:p)] <- res2[, paste0("ci_ll", 1:p)]<mu2 & mu2<res2[, paste0("ci_ul", 1:p)]

res3[, paste0("bias", 1:p)] <- res3[, paste0("mu", 1:p)]-mu3
res3[, paste0("rmse", 1:p)] <- res3[, paste0("bias", 1:p)]^2
res3[, paste0("ci_ll", 5:p)] <- res3[, paste0("mu", 5:p)]-qt(0.975, nA-10)*res3[, paste0("se", 5:p)]
res3[, paste0("ci_ul", 5:p)] <- res3[, paste0("mu", 5:p)]+qt(0.975, nA-10)*res3[, paste0("se", 5:p)]
res3[, paste0("cov", 1:p)] <- res3[, paste0("ci_ll", 1:p)]<mu3 & mu3<res3[, paste0("ci_ul", 1:p)]

tb <- as.data.frame(matrix(NA, p, 28))
names(tb) <- c("method", "bias_y1", "rMSE_y1", "cov_rate_y1", "se_ratio_y1", "MeanSE_y1", "TrueSE_y1", "Q5_y1", "Q95_y1", "biasSE_y1", "bias_y2", "rMSE_y2", "cov_rate_y2", "se_ratio_y2", "MeanSE_y2", "TrueSE_y2", "Q5_y2", "Q95_y2", "biasSE_y2", "bias_y3", "rMSE_y3", "cov_rate_y3", "se_ratio_y3", "MeanSE_y3", "TrueSE_y3", "Q5_y3", "Q95_y3", "biasSE_y3")
tb$method <- c("UW_PS", "TW_PS", "UW_NS", "TW_NS", paste(out10$Method, "TT", sep="_"), paste(out10$Method, "TF", sep="_"), paste(out10$Method, "FT", sep="_"), paste(out10$Method, "FF", sep="_"))

tb[, 2] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)mean(x, na.rm=T))/mu1
tb[, 3] <- sqrt(apply(res1[, paste0("rmse", 1:p)], 2, function(x)mean(x, na.rm=T)))
tb[, 4] <- 100*apply(res1[, paste0("cov", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 5] <- apply(res1[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))/apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 6] <- apply(res1[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 7] <- apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 8] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.05, na.rm=T))/mu1
tb[, 9] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.95, na.rm=T))/mu1
tb[, 10] <- 100*apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))/mu1

tb[, 11] <- 100*apply(res2[, paste0("bias", 1:p)], 2, function(x)mean(x, na.rm=T))/mu2
tb[, 12] <- sqrt(apply(res2[, paste0("rmse", 1:p)], 2, function(x)mean(x, na.rm=T)))
tb[, 13] <- 100*apply(res2[, paste0("cov", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 14] <- apply(res2[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))/apply(res2[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 15] <- apply(res2[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 16] <- apply(res2[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 17] <- 100*apply(res2[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.05, na.rm=T))/mu2
tb[, 18] <- 100*apply(res2[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.95, na.rm=T))/mu2
tb[, 19] <- 100*apply(res2[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))/mu2

tb[, 20] <- 100*apply(res3[, paste0("bias", 1:p)], 2, function(x)mean(x, na.rm=T))/mu3
tb[, 21] <- sqrt(apply(res3[, paste0("rmse", 1:p)], 2, function(x)mean(x, na.rm=T)))
tb[, 22] <- 100*apply(res3[, paste0("cov", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 23] <- apply(res3[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))/apply(res3[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 24] <- apply(res3[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 25] <- apply(res3[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 26] <- 100*apply(res3[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.05, na.rm=T))/mu3
tb[, 27] <- 100*apply(res3[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.95, na.rm=T))/mu3
tb[, 28] <- 100*apply(res3[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))/mu3

tb


write.csv(res1, paste0("res1_TT", "_", nR, "_", nA, "_IPSW.csv"), col.names=T, sep=",")
write.csv(res2, paste0("res2_TT", "_", nR, "_", nA, "_IPSW.csv"), col.names=T, sep=",")
write.csv(res3, paste0("res3_TT", "_", nR, "_", nA, "_IPSW.csv"), col.names=T, sep=",")
write.csv(tb, paste0("tb_TT", "_", nR, "_", nA, "_IPSW.csv"), col.names=T, sep=",")

colnames(wcor) <- out10$Model[1:(b-1)]
colnames(pcor) <- out10$Model[1:(b-1)]
wc <- apply(wcor, 2, function(x)mean(x, na.rm=T))
pc <- apply(pcor, 2, function(x)mean(x, na.rm=T))

write.csv(wc, paste0("wcor_TT", "_", nR, "_", nA, "_IPSW.csv"), col.names=T, sep=",")
write.csv(pc, paste0("pcor_TT", "_", nR, "_", nA, "_IPSW.csv"), col.names=T, sep=",")

