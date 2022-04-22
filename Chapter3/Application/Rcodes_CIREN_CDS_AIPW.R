

#########################################################
#   Date: 04/12/2022
#   Author: Ali Rafei
#   Email: arafei@umich.edu
#########################################################


# Goal: To adjust for selection bias in
# Crash Injury Research Engineering System (CIREN) data
# Using Crashworthiness Data System (CDS) as the benchmark

# The two datasets are already harmonized and available
# in "CIREN_cleaed.cav" and "CDS_cleaned.csv"
# Outcomes: injury rate in: Head, Abdomen, Thorax, Lower_extermity, Spine
# Covariates: Sex, Age_g, Race, Ethn, modelYear_g, vehMake, 
#		  vehType, Role, seatRow, resUse, bagDep, lightCon, 
#		  defLoc, Rollover, hospDay_g, damDist, crashType, 
#		  direct, Fatal, deltaV, injScale

# Note that for simplicity, we ignore the clusters/strata of CDS

# Set the working directory
setwd("D:\\My Thesis\\Chapter III\\CIREN")

# Read the user-defined functions
source("func_AIPW.R")

# Read cleaned harmonized data of CIREN and CDS
# CIREN dataset involves covariates + outcomes
# CDS dataset involves covariates + outcome + sampling weights 
CIREN <- read.csv("CIREN_cleaned.csv", header=T, sep=",")
CDS   <- read.csv("CDS_cleaned.csv", header=T, sep=",")

# Set the weights in CIREN to 1 temp according to Valliant & Dever 2011
CIREN$wght <- 1

# Create the indicator (Z_i) of being in CIREN
CIREN$Z <- 1
CDS$Z <- 0

# Append the two datasets
CIREN_CDS <- rbind(CIREN, CDS)

# list of covariates for different type of modeling 
# picked based on prior stepwise(AIC) variable selection
# For propensity modeling p(Z_i=1 | X_i)
nXz <- c("Age_g", "modelYear_g", "vehMake", "hospDay_g", "resUse", "bagDep", "deltaV", "injScale")
# For outcome modeling 
nXy_head <- c("Age_g", "vehType", "defLoc", "injScale", "Fatal")
nXy_abd <- c("Age_g", "vehMake", "vehType", "resUse", "defLoc", "damDist", "deltaV", "injScale", "Fatal")
nXy_thrx <- c("Sex", "Age_g", "modelYear_g", "bagDep", "Rollover", "hospDay_g", "damDist", "deltaV", "injScale")
nXy_lext <- c("Sex", "Age_g", "modelYear_g", "resUse", "Rollover", "hospDay_g", "damDist", "deltaV", "injScale", "Fatal")
nXy_spin <- c("Age_g", "modelYear_g", "seatRow", "Rollover", "hospDay_g", "damDist", "injScale")

# drop item-level missingness in covariates
CIREN_CDS1 <- CIREN_CDS[complete.cases(CIREN_CDS[, nXz]), ]

# Creating pseudo-weights based on
# 1. IPSW: Valliant & Dever (2011)
fit1 <- wglm(X=data.frame(model.matrix(~., data=CIREN_CDS1[, nXz])[, -1]), Y=CIREN_CDS1$Z, di=CIREN_CDS1$wght)
CIREN_CDS1$pw_IPSW <- (1-fit1$predict)/fit1$predict

# 2. PMLE: Chen et al (2011)
fit2 <- chen_glm(X=data.frame(model.matrix(~., data=CIREN_CDS1[, nXz])[, -1]), Y=CIREN_CDS1$Z, di=CIREN_CDS1$wght)
CIREN_CDS1$pw_PMLE <- 1/fit2$predict

# 3. PAPP: Elliott & Valliant (2017)
# Step 1: modeling sampling weights
CIREN_CDS1$Pi <- 1/CIREN_CDS1$wght
fit3 <- betareg(as.formula(paste("Pi", paste(nXz, collapse = " + "), sep = " ~ ")), data=CIREN_CDS1[CIREN_CDS1$Z==0, ], link="logit")
CIREN_CDS1$Pi[CIREN_CDS1$Z==1] <- predict(fit3, newdata=CIREN_CDS1[CIREN_CDS1$Z==1, nXz])
CIREN_CDS1$wght[CIREN_CDS1$Z==1] <- 1/CIREN_CDS1$Pi[CIREN_CDS1$Z==1]
# Step 1: modeling Z_i
fit4 <- glm(as.formula(paste("Z", paste(nXz, collapse = " + "), sep = " ~ ")), data=CIREN_CDS1, family=binomial(link=logit))
CIREN_CDS1$Pz <- fit4$fitted.values
CIREN_CDS1$pw_PAPP <- CIREN_CDS1$wght*(1-CIREN_CDS1$Pz)/CIREN_CDS1$Pz

# Estimate population total
N <- round(sum(CIREN_CDS1$wght[CIREN_CDS1$Z==0]))

# Normalizing estimated pseudo-weights
CIREN_CDS1$pw_PAPP[CIREN_CDS1$Z==1] <- N*CIREN_CDS1$pw_PAPP[CIREN_CDS1$Z==1]/sum(CIREN_CDS1$pw_PAPP[CIREN_CDS1$Z==1], na.rm=T)
CIREN_CDS1$pw_IPSW[CIREN_CDS1$Z==1] <- N*CIREN_CDS1$pw_IPSW[CIREN_CDS1$Z==1]/sum(CIREN_CDS1$pw_IPSW[CIREN_CDS1$Z==1], na.rm=T)
CIREN_CDS1$pw_PMLE[CIREN_CDS1$Z==1] <- N*CIREN_CDS1$pw_PMLE[CIREN_CDS1$Z==1]/sum(CIREN_CDS1$pw_PMLE[CIREN_CDS1$Z==1], na.rm=T)

#draw the box-plot of estimated pseudo-weights in CIREN
mm <- factor(c(rep(0, nrow(CIREN_CDS1[CIREN_CDS1$Z==0, ])), rep(1:2, rep(nrow(CIREN_CDS1[CIREN_CDS1$Z==1, ]), 2))), levels=0:2, labels=c("CDS_WGHT", c("CIREN_PAPP", "CIREN_PMLE")))
dfw <- data.frame(Method=mm, wght=c(CIREN_CDS1$wght[CIREN_CDS1$Z==0], CIREN_CDS1$pw_PAPP[CIREN_CDS1$Z==1], CIREN_CDS1$pw_IPSW[CIREN_CDS1$Z==1]))
ggplot(dfw, aes(x=Method, y=log(wght), fill=Method)) + 
	geom_boxplot() +
	theme_bw(base_size = 15)+
	labs(y="Weights (log scale)", x="Method")+
	ggtitle("(b)")+
	theme(legend.position="none", plot.title = element_text(hjust = 0.5))


#######################################################
#								 	#
#         distribution of covariates in CIREN/CDS	#
#									#
#######################################################

#set the sampling design
svydn_CDS <- svydesign(ids=~1, strata=NULL, weights=~wght, data=CIREN_CDS1[CIREN_CDS1$Z==0, ], nest=T)
svydn_CIREN <- svydesign(ids=~1, weights=~pw_PAPP, data=CIREN_CDS1[CIREN_CDS1$Z==1, ])


#draw the grouped bar plots
df1 <- 100*c(svymean(~factor(Sex), design=svydn_CIREN, na.rm=T), svymean(~factor(Sex), design=svydn_CDS, na.rm=T))
dff1 <- data.frame(Gender=factor(c(0, 1, 0, 1), levels=0:1, labels=c("Male", "Female")), 
		     Study=factor(c(0, 0, 1, 1), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov1 <- ggplot(data=dff1, aes(x=Gender, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Gender")+
		ylim(0, 100)+
		labs(x="", y="percent (%)")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(Age_g), design=svydn_CIREN, na.rm=T), svymean(~factor(Age_g), design=svydn_CDS, na.rm=T))
dff2 <- data.frame(AgeGroup=factor(c(0:3, 0:3), levels=0:3, labels=c("16-19", "20-39", "40-64", "65+")), 
		     Study=factor(c(rep(0, 4), rep(1, 4)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov2 <- ggplot(data=dff2, aes(x=AgeGroup, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Age group (yrs)")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(Race), design=svydn_CIREN, na.rm=T), svymean(~factor(Race), design=svydn_CDS, na.rm=T))
dff3 <- data.frame(Race=factor(c(0:2, 0:2), levels=0:2, labels=c("White", "Black", "Other")), 
		     Study=factor(c(rep(0, 3), rep(1, 3)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov3 <- ggplot(data=dff3, aes(x=Race, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Race")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(Ethn), design=svydn_CIREN, na.rm=T), svymean(~factor(Ethn), design=svydn_CDS, na.rm=T))
dff4 <- data.frame(Ethn=factor(c(0:1, 0:1), levels=0:1, labels=c("Non-Hispanic", "Hispanic")), 
		     Study=factor(c(rep(0, 2), rep(1, 2)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov4 <- ggplot(data=dff4, aes(x=Ethn, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Ethnicity")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(modelYear_g), design=svydn_CIREN, na.rm=T), svymean(~factor(modelYear_g), design=svydn_CDS, na.rm=T))
dff5 <- data.frame(modelYear=factor(c(0:3, 0:3), levels=0:3, labels=c("<2004", "2004-2007", "2008-2011", "2012+")), 
		     Study=factor(c(rep(0, 4), rep(1, 4)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov5 <- ggplot(data=dff5, aes(x=modelYear, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Model year")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(vehMake), design=svydn_CIREN, na.rm=T), svymean(~factor(vehMake), design=svydn_CDS, na.rm=T))
dff6 <- data.frame(carMake=factor(c(0, 1, 2, 3, 0, 1, 2, 3), levels=0:3, labels=c("American", "Japanese", "Korean", "Other")), 
		     Study=factor(c(1, 1, 1, 1, 2, 2, 2, 2), levels=1:2, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov6 <- ggplot(data=dff6, aes(x=carMake, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Vehicle make")+
		ylim(0, 100)+
		labs(x="", y="percent (%)")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(vehType), design=svydn_CIREN, na.rm=T), svymean(~factor(vehType), design=svydn_CDS, na.rm=T))
dff7 <- data.frame(vehType=factor(c(1:4, 1:4), levels=1:4, labels=c("Car", "SUV", "Van", "TRUCK")), 
		     Study=factor(c(1, 1, 1, 1, 2, 2, 2, 2), levels=1:2, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov7 <- ggplot(data=dff7, aes(x=vehType, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Vehicle type")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(Role), design=svydn_CIREN, na.rm=T), svymean(~factor(Role), design=svydn_CDS, na.rm=T))
dff8 <- data.frame(Role=factor(c(0:1, 0:1), levels=0:1, labels=c("Passenger", "Driver")), 
		     Study=factor(c(rep(0, 2), rep(1, 2)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov8 <- ggplot(data=dff8, aes(x=Role, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Occupant role")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))


df1 <- 100*c(svymean(~factor(seatRow), design=svydn_CIREN, na.rm=T), svymean(~factor(seatRow), design=svydn_CDS, na.rm=T))
dff9 <- data.frame(seatRow=factor(c(0:1, 0:1), levels=0:1, labels=c("Front", "Rear")), 
		     Study=factor(c(rep(0, 2), rep(1, 2)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov9 <- ggplot(data=dff9, aes(x=seatRow, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Seating row")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(lightCon), design=svydn_CIREN, na.rm=T), svymean(~factor(lightCon), design=svydn_CDS, na.rm=T))
dff10 <- data.frame(lightCon=factor(c(0:3, 0:3), levels=0:3, labels=c("Daylight", "Dark but lighted", "dark", "dawn/dusk")), 
		     Study=factor(c(rep(0, 4), rep(1, 4)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov10 <- ggplot(data=dff10, aes(x=lightCon, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Light condition")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(defLoc), design=svydn_CIREN, na.rm=T), svymean(~factor(defLoc), design=svydn_CDS, na.rm=T))
dff11 <- data.frame(defLoc=factor(c(0:3, 0:3), levels=0:3, labels=c("Front", "Left", "Right", "Top")), 
		     Study=factor(c(rep(0, 4), rep(1, 4)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov11 <- ggplot(data=dff11, aes(x=defLoc, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Deformation location")+
		ylim(0, 100)+
		labs(x="", y="percent (%)")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(hospDay_g), design=svydn_CIREN, na.rm=T), svymean(~factor(hospDay_g), design=svydn_CDS, na.rm=T))
dff12 <- data.frame(hospDay_g=factor(c(0:4, 0:4), levels=0:4, labels=c("0", "1-3", "4-7", "8-14", "15+")), 
		     Study=factor(c(rep(0, 5), rep(1, 5)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov12 <- ggplot(data=dff12, aes(x=hospDay_g, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Days hospitalized")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(damDist), design=svydn_CIREN, na.rm=T), svymean(~factor(damDist), design=svydn_CDS, na.rm=T))
dff13 <- data.frame(damDist=factor(c(0:3, 0:3), levels=0:3, labels=c("Wide", "Narrow", "Corner", "Rollover")), 
		     Study=factor(c(rep(0, 4), rep(1, 4)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov13 <- ggplot(data=dff13, aes(x=damDist, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Damage distribution")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(Rollover), design=svydn_CIREN, na.rm=T), svymean(~factor(Rollover), design=svydn_CDS, na.rm=T))
dff14 <- data.frame(Rollover=factor(c(0, 1, 0, 1), levels=0:1, labels=c("Not", "Yes")), 
		     Study=factor(c(0, 0, 1, 1), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov14 <- ggplot(data=dff14, aes(x=Rollover, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=3.5,   position=position_dodge(width=0.9))+
		ggtitle("Rollover status")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(injScale), design=svydn_CIREN, na.rm=T), svymean(~factor(injScale), design=svydn_CDS, na.rm=T))
dff15 <- data.frame(injScale=factor(c(0:3, 0:3), levels=0:3, labels=c("Moderate", "Serious", "Severe", "Critical")), 
		     Study=factor(c(rep(0, 4), rep(1, 4)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov15 <- ggplot(data=dff15, aes(x=injScale, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=5,   position=position_dodge(width=0.9))+
		ggtitle("Maximum Abbr Injury Scale")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(crashType), design=svydn_CIREN, na.rm=T), svymean(~factor(crashType), design=svydn_CDS, na.rm=T))
dff16 <- data.frame(crashType=factor(c(0:2, 0:2), levels=0:2, labels=c("Frontal", "Side", "Rollover")), 
		     Study=factor(c(rep(0, 3), rep(1, 3)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov16 <- ggplot(data=dff16, aes(x=crashType, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=5,   position=position_dodge(width=0.9))+
		ggtitle("Crash Type")+
		ylim(0, 100)+
		labs(x="", y="percent (%)")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(resUse), design=svydn_CIREN, na.rm=T), svymean(~factor(resUse), design=svydn_CDS, na.rm=T))
dff17 <- data.frame(resUse=factor(c(0:1, 0:1), levels=0:1, labels=c("No", "Yes")), 
		     Study=factor(c(rep(0, 2), rep(1, 2)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov17 <- ggplot(data=dff17, aes(x=resUse, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=5,   position=position_dodge(width=0.9))+
		ggtitle("Seatbelt Use")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(bagDep), design=svydn_CIREN, na.rm=T), svymean(~factor(resUse), design=svydn_CDS, na.rm=T))
dff18 <- data.frame(bagDep=factor(c(0:1, 0:1), levels=0:1, labels=c("No", "Yes")), 
		     Study=factor(c(rep(0, 2), rep(1, 2)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov18 <- ggplot(data=dff18, aes(x=bagDep, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=5,   position=position_dodge(width=0.9))+
		ggtitle("Airbag Deployed")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(Fatal), design=svydn_CIREN, na.rm=T), svymean(~factor(Fatal), design=svydn_CDS, na.rm=T))
dff19<- data.frame(Fatal=factor(c(0:1, 0:1), levels=0:1, labels=c("No", "Yes")), 
		     Study=factor(c(rep(0, 2), rep(1, 2)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov19 <- ggplot(data=dff19, aes(x=Fatal, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=5,   position=position_dodge(width=0.9))+
		ggtitle("Fatal")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

df1 <- 100*c(svymean(~factor(deltaV), design=svydn_CIREN, na.rm=T), svymean(~factor(deltaV), design=svydn_CDS, na.rm=T))
dff20<- data.frame(deltaV=factor(c(0:2, 0:2), levels=0:2, labels=c("Minor", "Moderate", "Severe")), 
		     Study=factor(c(rep(0, 3), rep(1, 3)), levels=0:1, labels=c("CIREN", "CDS")),
		     Frequency=as.numeric(df1))
pcov20 <- ggplot(data=dff20, aes(x=deltaV, y=Frequency, fill=Study)) +
		geom_bar(stat="identity", position="dodge")+
		geom_text(aes(label=round(Frequency, 1)), vjust=-0.3, color="black", size=5,   position=position_dodge(width=0.9))+
		ggtitle("deltaV")+
		ylim(0, 100)+
		labs(x="", y="")+
		theme_light(base_size = 15)+
		theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

pdf("cat_covs_PAPP_bar_minX.pdf", width=45, height=20)
grid.arrange(pcov1, pcov2, pcov3, pcov4, pcov5, pcov6, pcov7, pcov8, pcov9, pcov10, pcov11, pcov12, pcov13, pcov14, pcov15, pcov16, pcov17, pcov18, pcov19, pcov20, ncol=5, nrow=4, heights=rep(4, 4), widths=c(rep(5, 4), 7))
dev.off()


#######################################################
#									                                    #
#                    AIPW estimation			            #
#									                                    #
#######################################################

CIREN_CDS2 <- CIREN_CDS1[complete.cases(CIREN_CDS1[, union(nXz, nXy_head)]), ]
out10_head <- AIPW_KH(y=CIREN_CDS2$Head, x=model.matrix(~., data=CIREN_CDS2[, nXy_head])[, -1], deltaA=CIREN_CDS2$Z, sw=CIREN_CDS2$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), family='binomial')

CIREN_CDS2 <- CIREN_CDS1[complete.cases(CIREN_CDS1[, union(nXz, nXy_abd)]), ]
out10_abd  <- AIPW_KH(y=CIREN_CDS2$Abd , x=model.matrix(~., data=CIREN_CDS2[, nXy_abd ])[, -1], deltaA=CIREN_CDS2$Z, sw=CIREN_CDS2$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), family='binomial')

CIREN_CDS2 <- CIREN_CDS1[complete.cases(CIREN_CDS1[, union(nXz, nXy_thrx)]), ]
out10_thrx <- AIPW_KH(y=CIREN_CDS2$Thrx, x=model.matrix(~., data=CIREN_CDS2[, nXy_thrx])[, -1], deltaA=CIREN_CDS2$Z, sw=CIREN_CDS2$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), family='binomial')

CIREN_CDS2 <- CIREN_CDS1[complete.cases(CIREN_CDS1[, union(nXz, nXy_lext)]), ]
out10_lext <- AIPW_KH(y=CIREN_CDS2$Lext, x=model.matrix(~., data=CIREN_CDS2[, nXy_lext])[, -1], deltaA=CIREN_CDS2$Z, sw=CIREN_CDS2$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), family='binomial')

CIREN_CDS2 <- CIREN_CDS1[complete.cases(CIREN_CDS1[, union(nXz, nXy_spin)]), ]
out10_spin <- AIPW_KH(y=CIREN_CDS2$Spin, x=model.matrix(~., data=CIREN_CDS2[, nXy_spin])[, -1], deltaA=CIREN_CDS2$Z, sw=CIREN_CDS2$wght, method=c("PAPP", "PMLE", "IPSW", "PM", "AIPW"), family='binomial')

#################################################

nxn <- c("Head", "Abd", "Thrx", "Lext", "Spin")

svydn_CDSW <- svydesign(ids=~1, weights=~wght, data=CIREN_CDS1[CIREN_CDS1$Z==0, ])
svydn_CDSU <- svydesign(ids=~1, weights=~1, data=CIREN_CDS1[CIREN_CDS1$Z==0, ])
svydn_CIRENU <- svydesign(ids=~1, weights=~1, data=CIREN_CDS1[CIREN_CDS1$Z==1, ])

CDSU_m <- svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSU, na.rm=T)
CDSU_ci <- confint(svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSU, na.rm=T))
CDSW_m <- svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSW, na.rm=T)
CDSW_ci <- confint(svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSW, na.rm=T))

CIRENU_m <- svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CIRENU, na.rm=T)
CIRENU_ci <- confint(svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CIRENU, na.rm=T))

result <- data.frame(matrix(NA, 9, 16))
names(result) <- c("Method", nxn, paste(nxn, "ll", sep="_"), paste(nxn, "ul", sep="_"))
result$Method <- factor(c(0, 1, 0, 2, 3, 4, 5, 6, 7), levels=0:7, labels=c("Naive", "FW", "PAPP", "PMLE", "IPSW", "PM", "AIPW_PAPP", "AIPW_PMLE")) 
result$Robust <- factor(c(0, 0, 0, 1, 1, 1, 1, 2, 2), levels=0:2, labels=c("unweighted", "non-robust", "robust"))
result$Study <-  factor(c(0, 0, 1, 1, 1, 0, 2, 2, 2), levels=0:2, labels=c("CDS", "CIREN", "CIREN/CDS"))
result[1, 2:6] <- svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSU, na.rm=T)
result[1, 7:11] <- confint(svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSU, na.rm=T))[, 1]
result[1, 12:16] <- confint(svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSU, na.rm=T))[, 2]
result[2, 2:6] <- svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSW, na.rm=T)
result[2, 7:11] <- confint(svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSW, na.rm=T))[, 1]
result[2, 12:16] <- confint(svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CDSW, na.rm=T))[, 2]
result[3, 2:6] <- svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CIRENU, na.rm=T)
result[3, 7:11] <- confint(svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CIRENU, na.rm=T))[, 1]
result[3, 12:16] <- confint(svymean(~Head+Abd+Thrx+Lext+Spin, design=svydn_CIRENU, na.rm=T))[, 2]

result[4:9, 2] <- out10_head$Mu
result[4:9, 3] <- out10_abd$Mu
result[4:9, 4] <- out10_thrx$Mu
result[4:9, 5] <- out10_lext$Mu
result[4:9, 6] <- out10_spin$Mu

result[4:9, 7] <- out10_head$LCL
result[4:9, 8] <- out10_abd$LCL
result[4:9, 9] <- out10_thrx$LCL
result[4:9, 10] <- out10_lext$LCL
result[4:9, 11] <- out10_spin$LCL

result[4:9, 12] <- out10_head$UCL
result[4:9, 13] <- out10_abd$UCL
result[4:9, 14] <- out10_thrx$UCL
result[4:9, 15] <- out10_lext$UCL
result[4:9, 16] <- out10_spin$UCL

result1 <- result
result1[1:9, 2:16] <- 100*result1[1:9, 2:16]

# Plot the adjusted estimate

library(ggplot2)
library(gridExtra)
out1 <- ggplot(data = result1[-(2), ], aes(Method, Head, shape=Study, colour=Robust)) +
    	  geom_point(size = 6, position=position_dodge(width=0.5)) +
	  geom_rect(data = NULL, aes(xmin = 0, xmax = 8, ymin = result1[2, 7], ymax = result1[2, 12]), alpha=0.05, color = NA)+
	  geom_hline(yintercept = result1[2, 2], colour="red", linetype="dashed", size=1) +
    	  geom_errorbar(
        aes(ymin = Head_ll, ymax = Head_ul),
        width = 0.2, size=1.5,
        linetype = "solid",
        position=position_dodge(width=0.5)) +
	  labs(x="", y="Percent (%)") +
	  ylim(-1, 50) +
	  ggtitle("Head") +
    	  theme_bw(base_size = 30) + 
	  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), plot.title = element_text(hjust = 0.5), legend.position="none")

out2 <- ggplot(data = result1[-(2), ], aes(Method, Abd, shape=Study, colour=Robust)) +
    	  geom_point(size = 6, position=position_dodge(width=0.5)) +
	  geom_rect(data = NULL, aes(xmin = 0, xmax = 8, ymin = result1[2, 8], ymax = result1[2, 13]), alpha=0.05, color = NA)+
	  geom_hline(yintercept = result1[2, 3], colour="red", linetype="dashed", size=1) +
    	  geom_errorbar(
        aes(ymin = Abd_ll, ymax = Abd_ul),
        width = 0.2, size=1.5,
        linetype = "solid",
        position=position_dodge(width=0.5)) +
	  labs(x="", y="", color="Robustness") +
	  ylim(-1, 50) +
	  #guides(shape=FALSE)+
	  ggtitle("Abdomen") +
    	  theme_bw(base_size = 30) + 
	  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust = 0.5))

out3 <- ggplot(data = result1[-(2), ], aes(Method, Thrx, shape=Study, colour=Robust)) +
    	  geom_point(size = 6, position=position_dodge(width=0.5)) +
	  geom_rect(data = NULL, aes(xmin = 0, xmax = 8, ymin = result1[2, 9], ymax = result1[2, 14]), alpha=0.05, color = NA)+
	  geom_hline(yintercept = result1[2, 4], colour="red", linetype="dashed", size=1) +
    	  geom_errorbar(
        aes(ymin = Thrx_ll, ymax = Thrx_ul),
        width = 0.2, size=1.5,
        linetype = "solid",
        position=position_dodge(width=0.5)) +
	  labs(x="Method", y="Percent (%)") +
	  ylim(-1, 50) +
	  ggtitle("Thorax") +
    	  theme_bw(base_size = 30) + 
	  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x = element_text(angle = 20, hjust = 1))

out4 <- ggplot(data = result1[-(2), ], aes(Method, Lext, shape=Study, colour=Robust)) +
    	  geom_point(size = 6, position=position_dodge(width=0.5)) +
	  geom_rect(data = NULL, aes(xmin = 0, xmax = 8, ymin = result1[2, 10], ymax = result1[2, 15]), alpha=0.05, color = NA)+
	  geom_hline(yintercept = result1[2, 5], colour="red", linetype="dashed", size=1) +
    	  geom_errorbar(
        aes(ymin = Lext_ll, ymax = Lext_ul),
        width = 0.2, size=1.5,
        linetype = "solid",
        position=position_dodge(width=0.5)) +
	  labs(x="Method", y="") +
	  ylim(-1, 50) +
	  ggtitle("Lower extreminty") +
    	  theme_bw(base_size = 30) + 
	  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x = element_text(angle = 20, hjust = 1))

out5 <- ggplot(data = result1[-(2), ], aes(Method, Spin, shape=Study, colour=Robust)) +
    	  geom_point(size = 6, position=position_dodge(width=0.5)) +
	  geom_rect(data = NULL, aes(xmin = 0, xmax = 8, ymin = result1[2, 11], ymax = result1[2, 16]), alpha=0.05, color = NA)+
	  geom_hline(yintercept = result1[2, 6], colour="red", linetype="dashed", size=1) +
    	  geom_errorbar(
        aes(ymin = Spin_ll, ymax = Spin_ul),
        width = 0.2, size=1.5,
        linetype = "solid",
        position=position_dodge(width=0.5)) +
	  labs(x="Method", y="", color="Robustness") +
	  #guides(shape=FALSE)+
	  ylim(-1, 50) +
	  ggtitle("Spine") +
    	  theme_bw(base_size = 30) + 
	  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1))

pdf("out_11_p3.pdf", width=32, height=15)
grid.arrange(out1, out2, out3, out4, ncol=2, nrow=2, heights=c(7, 7.3), widths=c(8, 9.2))
dev.off()
