# Doubly robust inference for non-probability surveys
# author: Ali Rafei
# E-mail: arafei@umich.edu
# Date  : 12-06-2021
# Paper : https://arxiv.org/abs/2101.07456

# The following function provides doubly robust (DR) inference 
# for a non-probability survey (S_A) where a parallel probability
# survey (S_R) is available as the reference survey. In both
# samples, it is assumed a common set of auxiliary variables, X, 
# are observed, but the outcome, Y, is supposed to be observed 
# only in S_A. On the other hand, S_R comes with a set of
# sampling weights, while pseudo-weights in S_A are unoberved.

# The function AIPW_KH estimates the finite population mean for a 
# continuous or binary outcome [family=c('gaussian', 'binomial')] 
# based on the Augmented Inverse Propensity Weighting (AIPW) idea 
# proposed by Kim & Haziza (2014) that yields DR point and variance 
# estimates by solving the estimating equations associated with the 
# proponsity model and outcome models.

# Major limitations:
# 1. The same set of X should be used in propensity and outcome models.
# 2. Unique solutions in the joint EE are not guaranteed.
# 3. The sampling weights of S_R are calculable for units of S_A.

# The function uses three different propensity approaches:
# 1. PMLE proposed by Chen et al (2019)
# 2. PAPP proposed by Rafei et al (2020)
# 3. IPSW proposed by Wang et al (2020)

# Inputs of AIPW_KH:
# y is a vector of outcome; NA for the probability sample (n x 1)
# x is a matrix of covariates without intercept  (n x p)
# deltaA is a vector of the binary indicator of belonging to the nonprobability sample;
# i.e., 1 if the unit belongs to the nonprobability sample, and 0 otherwise  (n x 1)
# sw is a vector of weights: 1 for the unit in the nonprobability sample and the 
# design weight for the unit in the probability sample  (n x 1) family specifies 
# the outcome model
# family = "gaussian": a linear regression model for the continuous outcome
# family = "binomial": a logistic regression model for the binary outcome
# The sampling score is a logistic regression model for the probability of selection into the nonprobability sample given X
# The outcome model is a linear regression model for continuous outcome or a logistic regression model for binary outcome
# spec: Model specification status (TRUE, FALSE) used for the simulation
# Method = c('PAPP', 'PMLE', 'IPSW', 'PM', 'AIPW')
#
# return:
# Method: The name of the method used
# Mu: adjusted mean
# SE: adjusted SE
# RW: the associated estimated pseudo weights for units of S_A

# install.packages("rootSolve")  # @import MASS rootSolve ncvreg stats
library(rootSolve)  #needed for multiroot() to solve joint estimating equations
library(sas7bdat)
library(ggplot2)
library(gridExtra)
library(survey)
library(gtools)
library(betareg)
library(MASS)

AIPW_KH <- function(y, x, deltaA, sw, method=c("PAPP", "PMLE", "IPSW", "PM", "IPSW"), family=c('gaussian', 'binomial')){

# Model specification in the simulation  
	PS_f <- PM_f <- colnames(x)
	n<-nrow(x) #combined sample size
	pp <- length(PS_f) #number of covariates in PS model
	pm <- length(PM_f) #number of covariates in PM model
	x.ARp <- as.matrix(cbind(1, x[, PS_f])) # Pooled design matrix for PS
	x.ARm <- as.matrix(cbind(1, x[, PM_f])) # Pooled design martix for PM
	y.AR <- y # Pooled ouctome: NA for units of S_R
	y.AR[which(is.na(y.AR))] <- 0 # Set NAs to 0
	nA <- sum(deltaA) # Sample size in S_A
	nR <- sum(1-deltaA) # Sample size in S_R
	nAR <- nA + nR # Total sample size
	loc.A <- which(deltaA==1) # Indicators of being in S_A
	loc.R <- which(deltaA==0) # Indicators of being in S_R
	x.Rp <- x.ARp[loc.R,] # X of PS model in S_R
	x.Rm <- x.ARm[loc.R,] # X of PM model in S_R
	x.Ap <- x.ARp[loc.A,] # X of PS model in S_A
	x.Am <- x.ARm[loc.A,] # X of PM model in S_A
	y.A <- y.AR[loc.A] # Outcome in S_A
	N <- sum(sw[loc.R]) # Estimated population size
	sw.A <- sw[loc.A] # Sampling weights of S_R in S_A
	sw.R <- sw[loc.R] # Sampling weights of S_R in S_R

	nms <- out_m <- out_se <- pwght <- c() #empty oututs
	
	# Estimate propensity weighted estimate and associated asmptotic SE based on PAPP
	if("PAPP"%in%method){
		fit0 <- glm(z~., data=data.frame(X=x.ARp[,-1], z=deltaA), family=binomial(link="logit"))  # Fit P(Z_i=1|X)
		ei <- exp(predict(fit0, type="link"))
		weight <- sw.A/ei[loc.A]  # Estimate PAPP weights
		est.pdr <- weighted.mean(y.A, w=weight, na.rm=T) # Comput weighted estimates
		# Compute asymptotic SE using TSL method
		wd <- (y.A-est.pdr)*weight
		b <- (t(wd)%*%x.ARp[loc.A, ]/N)%*%ginv(t(x.ARp)%*%(x.ARp*ei/(N*((1+ei)^2))))
		A <- (1/N^2)*sum((1-1/weight)*(wd^2))
		B <- (t(wd/(1+ei[loc.A]))%*%x.ARp[loc.A, ]/(N^2))%*%t(b)
		C <- b%*%(t(x.ARp)%*%(x.ARp*ei/((N^2)*((1+ei)^2))))%*%t(b)
		se.pdr <- sqrt(A - 2*B + C)
		out_m <- c(out_m, est.pdr)
		out_se <- c(out_se, se.pdr)
		pwght <- cbind(pwght, weight)
		nms <- c(nms, "PAPP") # Store method name
	}
	if("PMLE"%in%method){
		W <- sw
		W[loc.A] <- 1
		PZi <- chen_glm(X=x.ARp[,-1], Y=deltaA, di=W)$predict  #Fit PS model based on Chen et al's PMLE method
		weight <- 1/PZi[loc.A] # Estimate PMLE weights
		# Compute asymptotic SE using TSL method 
		est.pdr <- weighted.mean(y.A, w=weight, na.rm=T)
		wd <- (y.A-est.pdr)*weight
		b <- (t((weight-1)*(y.A-est.pdr))%*%x.ARp[loc.A, ])%*%ginv(t(x.ARp[loc.R, ])%*%(x.ARp[loc.R, ]*sw.R*PZi[loc.R]*(1-PZi[loc.R])))
		D <- nR*cov(PZi[loc.R]*sw.R*x.ARp[loc.R, ])/(N^2)
		se.pdr <- sqrt((1/N^2)*sum((1-1/weight)*((wd-x.ARp[loc.A, ]%*%t(b))^2)) + b%*%D%*%t(b))
		out_m <- c(out_m, est.pdr)
		out_se <- c(out_se, se.pdr)
		pwght <- cbind(pwght, weight)
		nms <- c(nms, "PMLE")
	}
	if("IPSW"%in%method){
		W <- sw
		W[loc.A] <- 1
		PZi <- wglm(X=x.ARp[,-1], Y=deltaA, di=W)$predict
		weight <- (1-PZi[loc.A])/PZi[loc.A]
		est.pdr <- weighted.mean(y.A, w=weight, na.rm=T)
		wd <- (y.A-est.pdr)*weight
		b <- (t((est.pdr-y.A)/weight)%*%x.ARp[loc.A, ])%*%ginv(t(x.ARp[loc.R, ])%*%(x.ARp[loc.R, ]*PZi[loc.R]))
		D <- nR*cov(PZi[loc.R]*sw.R*x.ARp[loc.R, ])/(N^2)
		se.pdr <- sqrt((1/N^2)*sum((1-1/weight)*((wd-x.ARp[loc.A, ]%*%t(b))^2)) + b%*%D%*%t(b))
		out_m <- c(out_m, est.pdr)
		out_se <- c(out_se, se.pdr)
		pwght <- cbind(pwght, weight)
		nms <- c(nms, "IPSW")
	}
	if("PM"%in%method){
		if(family=="gaussian"){
			fit1 <- lm(y~., data=data.frame(x.Am[,-1], y=y.A))
			pred <- predict(fit1, newdata=data.frame(x.Rm[,-1]),type="response",se.fit = TRUE)
			m.A <- pred$fit
			est.pdr <- weighted.mean(m.A, w=sw.R, na.rm=T)
			c <- ginv(t(x.ARm[loc.A, ])%*%x.ARm[loc.A, ]/nA)%*%t(t(sw.R)%*%x.ARm[loc.R, ]/N)
			ve.pdr <- weighted.var(y0=m.A, w0=sw.R) + sum((fit1$residuals^2)*(x.ARm[loc.A, ]%*%c)^2)/(nA^2)
			se.pdr<-sqrt(ve.pdr)
		}else if(family=="binomial"){
			fit1 <- glm(y~., data=data.frame(x.Am[,-1], y=y.A), family=binomial(link="logit"))
			m.A <- predict(fit1, newdata=data.frame(x.Rm[,-1]), type="response")
			m.B <- predict(fit1, type="response")
			est.pdr <- weighted.mean(m.A, w=sw[loc.R], na.rm=T)
			c <- ginv(t(x.ARm[loc.A, ])%*%x.ARm[loc.A, ]/nA)%*%t(t(sw.R)%*%x.ARm[loc.R, ]/N)
			ve.pdr <- weighted.var(y0=m.A, w0=sw.R) + sum((fit1$residuals^2)*((x.ARm[loc.A, ]*m.B*(1-m.B))%*%c)^2)/(nA^2)
			se.pdr<-sqrt(ve.pdr)			
		}
		out_m <- c(out_m, est.pdr)
		out_se <- c(out_se, se.pdr)
		nms <- c(nms, "PM")
	}
	if("PAPP"%in%method & "AIPW"%in%method){
		if(family=="gaussian"){
			##################################################
			## joint estimating equations
			##################################################

			par0 <- rep(0, (pp + pm + 2))
			par0[1:(pp+1)] <- glm(z~., data=data.frame(X=x.ARp[,-1], z=deltaA))$coefficients
			par0[(pp+2):(pp + pm + 2)] <- lm(y~., data=data.frame(X=x.Am[,-1], y=y.A))$coefficients
			jeepar<-rootSolve::multiroot(Uee_papp,par0,deltaA=deltaA,x.ARp=x.ARp,x.ARm=x.ARm,y.AR=y.AR,sw=sw)$root
			jeepar1<-jeepar[ 1:(pp+1)]
			jeepar2<-jeepar[(pp+2):(pp + pm + 2)]
			lpi.est <- exp(as.vector(x.ARp %*% as.matrix(jeepar1)))
			piZ.est  <- lpi.est/(1+lpi.est)
			m.est <- as.vector(x.ARm%*% as.matrix(jeepar2))
			piB.est <- lpi.est/sw

			##################################################
			## doubly robust estimation for FPM
			##################################################
			est.pdr <- sum((y.AR-m.est)*deltaA/piB.est)/N + sum(m.est*(1-deltaA)*sw)/N
			m.A <- m.est[loc.R]
			sigmasqhat <- sum((y.AR[loc.A]-m.est[loc.A])^2)/(nA-pm-1)
			ve.pdr <- weighted.var(y0=m.A, w0=sw.R) + (sum(deltaA*(1-2*piB.est )/piB.est ^2) + (N )) * sigmasqhat/(N^2)
			se.pdr<-sqrt(ve.pdr)

		}else if(family=="binomial"){
			##################################################
			## joint estimating equations
			##################################################
			par0 <- rep(0, (pp + pm + 2))
			par0[1:(pp+1)] <- glm(z~., data=data.frame(X=x.ARp[,-1], z=deltaA), family=binomial(link="logit"))$coefficients
			par0[(pp+2):(pp+pm+2)] <- glm(y~., data=data.frame(X=x.Am[,-1], y=y.A), family=binomial(link="logit"))$coefficients
			jeepar<-rootSolve::multiroot(Uee_binary_papp,par0,deltaA=deltaA,x.ARp=x.ARp,x.ARm=x.ARm,y.AR=y.AR,sw=sw)$root
			jeepar1<-jeepar[ 1:(pp + 1)]
			jeepar2<-jeepar[(pp + 2):(pp + pm + 2)]
			lpi.est <- exp(as.vector(x.ARp %*% as.matrix(jeepar1)))
			piZ.est  <- lpi.est/(1+lpi.est)
			lm.est <- as.vector(x.ARm%*% as.matrix(jeepar2))
			m.est  <- expoit(lm.est)
			piB.est <- lpi.est/sw

			##################################################
			## doubly robust estimation for FPM
			##################################################
			est.pdr <- sum((y.AR-m.est)*deltaA/piB.est)/N+(sum(m.est*(1-deltaA)*sw))/N
			m.A <- m.est[loc.R]
			sigmasqhat <- m.est*(1-m.est)
			ve.pdr <- weighted.var(y0=m.A, w0=sw.R) + (sum(sigmasqhat*deltaA*(1-2*piB.est )/piB.est ^2) + sum(sigmasqhat[loc.R]/sw.R))/(N^2)
			se.pdr<-sqrt(ve.pdr)
		}
		out_m <- c(out_m, est.pdr)
		out_se <- c(out_se, se.pdr)
		pwght <- cbind(pwght, 1/piB.est[loc.A])
		nms <- c(nms, "AIPW_PAPP")
	}
	if("PMLE"%in%method & "AIPW"%in%method){
		if(family=="gaussian"){

			##################################################
			## joint eestimating equations
			##################################################
			par0 <- rep(0, (pp + pm + 2))
			W <- sw
			W[loc.A] <- 1
			par0[1:(pp+1)] <- chen_glm(X=x.ARp[,-1], Y=deltaA, di=W)$coef
			par0[(pp+2):(pp + pm + 2)] <- lm(y~., data=data.frame(X=x.Am[,-1], y=y.A))$coefficients
			jeepar<-rootSolve::multiroot(Uee_pmle,par0,deltaA=deltaA,x.ARp=x.ARp,x.ARm=x.ARm,y.AR=y.AR,sw=sw)$root
			jeepar1<-jeepar[ 1:(pp+1)]
			jeepar2<-jeepar[(pp+2):(pp + pm + 2)]
			lpi.est <- exp(as.vector(x.ARp %*% as.matrix(jeepar1)))
			piB.est  <- lpi.est/(1+lpi.est)
			m.est <- as.vector(x.ARm%*% as.matrix(jeepar2))

			##################################################
			## doubly robust estimation for FPM
			##################################################
			est.pdr <- sum((y.AR-m.est)*deltaA/piB.est)/N + sum(m.est*(1-deltaA)*sw)/N
			m.A <- m.est[loc.R]
			sigmasqhat <- sum((y.AR[loc.A]-m.est[loc.A])^2)/(nA-pm-1)
			ve.pdr <- weighted.var(y0=m.A, w0=sw.R) + (sum(deltaA*(1-2*piB.est )/piB.est ^2) + (N )) * sigmasqhat/(N^2)
			se.pdr<-sqrt(ve.pdr)

		}else if(family=="binomial"){

			##################################################
			## joint eestimating equations
			##################################################
			par0 <- rep(0, (pp + pm + 2))
			W <- sw
			W[loc.A] <- 1
			par0[1:(pp+1)] <- chen_glm(X=x.ARp[,-1], Y=deltaA, di=W)$coef
			par0[(pp+2):(pp+pm+2)] <- glm(y~., data=data.frame(X=x.Am[,-1], y=y.A), family=binomial(link="logit"))$coefficients
			jeepar<-rootSolve::multiroot(Uee_binary_pmle,par0,deltaA=deltaA,x.ARp=x.ARp,x.ARm=x.ARm,y.AR=y.AR,sw=sw)$root
			jeepar1<-jeepar[ 1:(pp + 1)]
			jeepar2<-jeepar[(pp + 2):(pp + pm + 2)]
			lpi.est <- exp(as.vector(x.ARp %*% as.matrix(jeepar1)))
			piB.est  <- lpi.est/(1+lpi.est)
			lm.est <- as.vector(x.ARm%*% as.matrix(jeepar2))
			m.est  <- expoit(lm.est)

			##################################################
			## doubly robust estimation for FPM
			##################################################
			est.pdr <- sum((y.AR-m.est)*deltaA/piB.est)/N+(sum(m.est*(1-deltaA)*sw))/N
			m.A <- m.est[loc.R]
			sigmasqhat <- m.est*(1-m.est)
			ve.pdr <- weighted.var(y0=m.A, w0=sw.R) + (sum(sigmasqhat*deltaA*(1-2*piB.est )/piB.est ^2) + sum(sigmasqhat[loc.R]/sw.R))/(N^2)
		}
		out_m <- c(out_m, est.pdr)
		out_se <- c(out_se, se.pdr)
		pwght <- cbind(pwght, 1/piB.est[loc.A])
		nms <- c(nms, "AIPW_PMLE")
	}
	colnames(pwght) <- nms[nms!="PM"]
	names(out_m) <- names(out_se) <- nms
	return(list(Method=nms, Mu=out_m, SE=out_se, LCL=out_m-1.96*out_se, UCL=out_m+1.96*out_se, RW=pwght))
}


####################################################################################
## function list
####################################################################################

# joint score equation for alpha and beta of continuous outcome under PAPP
Uee_papp <- function(par, deltaA, x.ARp, x.ARm, y.AR, sw){
	pp<-ncol(x.ARp)-1
	pm<-ncol(x.ARm)-1
	nAR <- length(deltaA)
	alpha <- par[1:(pp+1)]
	beta <- par[(pp+2):(pp + pm + 2)]
	N<- sum(sw[deltaA==0])
	lpiZ <- exp(x.ARp%*%alpha)
	piZ <- lpiZ/(1+lpiZ)
	deltaA <- 1-deltaA
	y.AR[which(is.na(y.AR))] <- 0
	resAB <- (y.AR-x.ARm%*%beta)
	piZ <- as.vector(piZ)
	sw <- as.vector(sw)
	resAB <- as.vector(resAB)
	c(apply(x.ARm*(deltaA/piZ-1)*sw,2,sum), apply(x.ARp*deltaA*(1/piZ-1)*resAB*sw,2,sum))/N
}

# joint score equation for alpha and beta of continuous outcome under PMLE
Uee_pmle <- function(par, deltaA, x.ARp, x.ARm, y.AR, sw){
	pp<-ncol(x.ARp)-1
	pm<-ncol(x.ARm)-1
	nAR <- length(deltaA)
	alpha <- par[1:(pp+1)]
	beta <- par[(pp+2):(pp + pm + 2)]
	N<- sum(sw[deltaA==0])
	lpiB <- exp(x.ARp%*%alpha)
	piB <- lpiB/(1+lpiB)
	deltaA <- 1-deltaA
	y.AR[which(is.na(y.AR))] <- 0
	resAB <- (y.AR-x.ARm%*%beta)
	piB <- as.vector(piB)
	sw <- as.vector(sw)
	resAB <- as.vector(resAB)
	c(apply(x.ARm*deltaA/piB-x.ARm*deltaA*sw,2,sum), apply(x.ARp*deltaA*(1/piB-1)*resAB,2,sum))/nAR
}

# joint score equation for alpha and beta of binary outcome under PAPP
Uee_binary_papp <- function(par,deltaA,x.ARp,x.ARm,y.AR,sw){
	pp<-ncol(x.ARp)-1
	pm<-ncol(x.ARm)-1
	if(is.null(pp))pp=1
	if(is.null(pm))pm=1
	nAR<-length(deltaA)
	alpha<-par[1:(pp+1)]
	beta <-par[(pp+2):(pp+pm+2)]
	N<- sum(sw[deltaA==0])
	lpiZ<-x.ARp%*%alpha
	piZ <-exp(lpiZ)/(1+exp(lpiZ))
	deltaA<-1-deltaA
	mAB<-expoit(x.ARm%*%beta)
	resAB<-(y.AR- mAB)
	piZ<-as.vector(piZ)
	sw<-as.vector(sw)
	resAB<-as.vector(resAB)
	mAB<-as.vector(mAB)
	c(apply(x.ARm*mAB*(1-mAB)*(deltaA/piZ-1)*sw,2,sum), apply(x.ARp*deltaA*(1/piZ-1)*resAB*sw, 2, sum))/N
}

# joint score equation for alpha and beta of binary outcome under PMLE
Uee_binary_pmle <- function(par,deltaA,x.ARp,x.ARm,y.AR,sw){
	pp<-ncol(x.ARp)-1
	pm<-ncol(x.ARm)-1
	if(is.null(pp))pp=1
	if(is.null(pm))pm=1
	nAR<-length(deltaA)
	alpha<-par[1:(pp+1)]
	beta <-par[(pp+2):(pp+pm+2)]
	N<- sum(sw[deltaA==0])
	lpiB<-x.ARp%*%alpha
	piB <-exp(lpiB)/(1+exp(lpiB))
	deltaA<-1-deltaA
	mAB<-expoit(x.ARm%*%beta)
	resAB<-(y.AR- mAB)
	piB<-as.vector(piB)
	sw<-as.vector(sw)
	resAB<-as.vector(resAB)
	mAB<-as.vector(mAB)
	c(apply(x.ARm*deltaA/piB*mAB*(1-mAB) -x.ARm*deltaA*sw*mAB*(1-mAB),2,sum), apply(x.ARp*deltaA*(1/piB-1)*resAB,2,sum))/nAR
}

# Logistic function
expoit<-function(x){
	exp(x)/(1+exp(x))
}

# Function fitting weighted logistic regression based on Chen et al's method
chen_glm <- function(X, Y, di, threshold = 1e-10, max_iter = 100){
	X <- as.matrix(X)
	n <- nrow(X)
	np <- ncol(X) + 1
	X <- cbind(rep(1, n), X)
	calc_p <- function(X, beta){
		beta <- as.vector(beta)
		return(1 / (1+ exp(-X%*%beta)))
	}  
	beta <- rep(0, np)
	diff <- 10000 
	iter_count <- 0
	while(diff > threshold){
		pi <- as.vector(calc_p(X, beta))
		U <- t(X)%*%Y-t(di[Y==0]*X[Y==0, ])%*%pi[Y==0]
		WW <- pi[Y==0]*(1-pi[Y==0])*di[Y==0]
		H <- t(X[Y==0, ])%*%(WW*X[Y==0,])
		beta_change <- ginv(H) %*% U
		beta <- beta + beta_change
		diff <- sum(beta_change^2)
		iter_count <- iter_count + 1
		if(iter_count > max_iter) {
			stop("Not converging...")
		}
	} 
	list(coef=as.numeric(beta), predict=as.numeric(calc_p(X, beta)))
}

# Regular weighted logistic regression model based on PMLE 
wglm <- function(X, Y, di, threshold = 1e-10, max_iter = 100){
	X <- as.matrix(X)
	n <- nrow(X)
	np <- ncol(X) + 1
	X <- cbind(rep(1, n), X)
	calc_p <- function(X, beta){
		beta <- as.vector(beta)
		return(1 / (1+ exp(-X%*%beta)))
	}
	beta <- rep(0, np)
	diff <- 10000 
	iter_count <- 0
	while(diff > threshold){
		pi <- as.vector(calc_p(X, beta))
		U <- t(di*X)%*%(Y - pi)
		H <- t(X)%*%(di*pi*(1-pi)*X)
		beta_change = ginv(H) %*% U
		beta <- beta + beta_change
		diff <- sum(beta_change^2)
		iter_count <- iter_count + 1
		if(iter_count > max_iter) {
			stop("Not converging...")
		}
	} 
	list(coef=as.numeric(beta), predict=as.numeric(calc_p(X, beta)))
}


# Weighted SE based on TSL
weighted.se <- function(y0, w0, na.rm=T){
	if(na.rm){
		w0 <- w0[!is.na(y0)]
		y0 <- y0[!is.na(y0)]	
	}
	n <- length(y0)
	x <- sum(y0 * w0)
	z <- sum(w0)
	sqrt((var(y0*w0) + ((x/z)^2)*var(w0) - 2*(x/z)*cov(y0*w0, w0))*n/(z^2))
}

# Weighted variance based on TSL
weighted.var <- function(y0, w0, na.rm=T){
	if(na.rm){
		w0 <- w0[!is.na(y0)]
		y0 <- y0[!is.na(y0)]	
	}
	n <- length(y0)
	x <- sum(y0 * w0)
	z <- sum(w0)
	(var(y0*w0) + ((x/z)^2)*var(w0) - 2*(x/z)*cov(y0*w0, w0))*n/(z^2)
}
