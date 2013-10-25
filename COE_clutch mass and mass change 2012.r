##########################################################################
#
#	CLUTCH MASS and MASS CHANGE DURING INCUBATION analysis - Nome 2012
#
# Name - Samantha Franks
# 23 Oct 2013
#
##########################################################################


rm(list=ls()) #clear all previous
library(NCStats) # for Subset function
library(nlme)

### Data import and setup

d1 <- read.csv("R/sandpiper/COE Nome 2012.csv",header=T)
nest<-read.csv("R/sandpiper/nest data all sites 2008-2010 and 2012_clean.csv",header=T)
band<-read.csv("R/sandpiper/ALL banding breeding 2008-2010 and 2012.csv", header=T)
nest2012 <- Subset(nest, year=="2012")
band2012 <- Subset(band, Year=="2012" & fullnest != "12nome_no_nest")

d1$Band <- as.factor(d1$Band)

morpho <- with(band2012,data.frame(fullnest,Band,Flag,Sex,Age,Wing = Wing.mm,Culmen = Culmen.mm, Tarsus = Diagonal.tarsus.mm, Weight = Weight.g))

uniquebands <- unique(morpho)

# merge morphometrics with COE dataset
merge1 <- merge(uniquebands, d1, by.x="Band", by.y="Band", all.x = T, all.y = T)

# merge with nest dataset
merge2 <- merge(merge1,nest2012,by.x="fullnest",by.y="fullnest",all.x=T,all.y=T)

COE <- with(merge2, data.frame(fullnest, INIT.JULIAN, fate, failed, Flag = Flag.x, Band, Sex = Sex.x, dC, dN, Wing, Culmen, Tarsus, Weight, NumberEggs, LengthEgg1, LengthEgg2, LengthEgg3, LengthEgg4, WidthEgg1, WidthEgg2, WidthEgg3, WidthEgg4, AverageEggLength, AverageEggWidth, AverageChickTarsus, AverageChickCulmen, Response = AverageEggVolume.mL, EarlyWeight, LateWeight, WeightChange))

COE$NumberEggs <- as.factor(COE$NumberEggs)

coe <- Subset(COE, dC != "NA")

# create principal component of body size using wing and tarsus (PC1)
pc <- coe[,c("Wing","Tarsus")]
model <- prcomp(pc,scale=TRUE)
summary(model)
# biplot(model)
pc1<-predict(model)[,1] # predicted values of PC1

coe2 <- data.frame(coe,BodySize = pc1)

female <- Subset(coe2, Sex=="F" & dC != "NA")
male <- Subset(coe2, Sex=="M" & dC != "NA")

###-------------------------------------------------------------###
#	Analysis 1 - clutch mass vs female initiation date, dC, dN, and body size
###-------------------------------------------------------------###

dat <- female

modellist <- list()  # set up blank list to hold model formulas
modellist[[1]] <- formula(Response ~ INIT.JULIAN)
modellist[[2]] <- formula(Response ~ INIT.JULIAN + dC)
modellist[[3]] <- formula(Response ~ INIT.JULIAN + dN)
modellist[[4]] <- formula(Response ~ INIT.JULIAN + dC + dN)
modellist[[5]] <- formula(Response ~ dC)
modellist[[6]] <- formula(Response ~ dN)
modellist[[7]] <- formula(Response ~ dC + dN)
modellist[[8]] <- formula(Response ~ 1)
modellist[[9]] <- formula(Response ~ BodySize)
modellist[[10]] <- formula(Response ~ INIT.JULIAN + BodySize)
modellist[[11]] <- formula(Response ~ BodySize + dC)
modellist[[12]] <- formula(Response ~ BodySize + dN)
modellist[[13]] <- formula(Response ~ BodySize + dC + dN)
modellist[[14]] <- formula(Response ~ BodySize + INIT.JULIAN + dC + dN)
modellist[[15]] <- formula(Response ~ BodySize + INIT.JULIAN + dC)
modellist[[16]] <- formula(Response ~ BodySize + INIT.JULIAN + dN)
modellist[[17]] <- formula(Response ~ BodySize + NumberEggs)
modellist[[18]] <- formula(Response ~ BodySize + INIT.JULIAN + NumberEggs)
modellist[[19]] <- formula(Response ~ BodySize + NumberEggs + dC)
modellist[[20]] <- formula(Response ~ BodySize + NumberEggs + dN)
modellist[[21]] <- formula(Response ~ BodySize + NumberEggs + dC + dN)
modellist[[22]] <- formula(Response ~ BodySize + INIT.JULIAN + NumberEggs + dC)
modellist[[23]] <- formula(Response ~ BodySize + INIT.JULIAN + NumberEggs + dN)
modellist[[24]] <- formula(Response ~ BodySize + INIT.JULIAN + NumberEggs + dC + dN)
modellist[[25]] <- formula(Response ~ NumberEggs + INIT.JULIAN)
modellist[[26]] <- formula(Response ~ NumberEggs + INIT.JULIAN + dC)
modellist[[27]] <- formula(Response ~ NumberEggs + INIT.JULIAN + dN)
modellist[[28]] <- formula(Response ~ NumberEggs + INIT.JULIAN + dC + dN)
modellist[[29]] <- formula(Response ~ NumberEggs + dC)
modellist[[30]] <- formula(Response ~ NumberEggs + dN)
modellist[[31]] <- formula(Response ~ NumberEggs + dC + dN)
modellist[[32]] <- formula(Response ~ NumberEggs)


globalmodel <- lm(Response ~ INIT.JULIAN + dC + dN + BodySize + NumberEggs, data=dat)

ModelList<-modellist
coefnames <- names(globalmodel$coef)
modeloutput<-list()	# blank model output list

# AIC function
# a function to produce relevant calculations for an AIC table that is used in conjunction with a list containing the model names
calculate.AIC<-function(aictable,modellist) {
modelnames<-modellist
delta.aic<-aictable$AIC.c-min(aictable$AIC.c)
lik.aic<-exp(-delta.aic/2)
aic.w<-lik.aic/(sum(lik.aic))
aic.table<-data.frame(modelnames,AIC.table,delta.aic,lik.aic,aic.w)
}

# blank Coefficient Matrix
CoefMatrix <- matrix(NA, nrow=length(ModelList), ncol=length(coefnames), dimnames = list(c(1:length(ModelList)), coefnames))

# blank SE Matrix
SEMatrix <- matrix(NA, nrow=length(ModelList), ncol=length(coefnames), dimnames = list(c(1:length(ModelList)), coefnames))

# blank AIC matrix
AIC.table<-matrix(NA, nrow=length(ModelList), ncol=6, dimnames = list(c(1:length(ModelList)), c("r.squared","n.obs","df","-2loglik","AIC","AIC.c")))

# loop for calculating model output, filling CoefMatrix and SEMatrix

for(i in 1:length(ModelList)){

modeloutput[[i]]<-lm(ModelList[[i]], data=dat)
m<-modeloutput[[i]]

# fill in row "i" of the AIC table
aic<-AIC(m)
df<-length(m$coef) + 1 # +1 for the variance
n.obs<-nrow(dat)
AIC.table[i,"r.squared"]<-summary(m)$adj.r.squared
AIC.table[i,"n.obs"]<-n.obs
AIC.table[i,"df"]<-df
AIC.table[i,"-2loglik"]<- -2*logLik(m)
AIC.table[i,"AIC"]<-aic
aic.c<-aic + (2*df*(df+1))/(n.obs-df-1)
AIC.table[i,"AIC.c"]<-aic.c

# fill in row "i", column "j" of the Coefficient Matrix
for (j in 1:length(coefnames)) {
	
if (coefnames[j] == "(Intercept)" & length(names(m$coef)) == 1) {CoefMatrix[i,j] <- summary(m)$coef[coefnames[j],1]} else {
CoefMatrix[i,j] <- summary(m)$coef[,1][coefnames[j]]
# CoefMatrix[i,j] <- summary(m)$coef[coefnames[j],1] # gives subscript out of bounds error if coefnames[j] is not one of the parameters in the model
# CoefMatrix[i,j] <- summary(m)$coef[,1][coefnames[j]] gives NA if model is the NULL model (only 1 parameter, the intercept) because in that case, summary(m)$coef is a vector (1-dimension only) rather than a dataframe (2-dimensions)
# code a special case for the NULL (intercept only) model, all other models with other parameters get the "else" statement
} # close if else statement

# fill in row "i" of the Standard Error Matrix
if (coefnames[j] == "(Intercept)" & length(names(m$coef)) == 1) {SEMatrix[i,j] <- summary(m)$coef[coefnames[j],2]} else {
SEMatrix[i,j] <- summary(m)$coef[,2][coefnames[j]]
} # close if else statement

} # close inner "j" loop

} # close outer "i" loop

# calculate deltaAICs, likelihoods, and AIC weights
AIC.output<-data.frame(AIC.table)
aic<-calculate.AIC(AIC.output,as.character(ModelList))
#for (i in 4:9) {
# aic[,i]<-round(aic[,i], digits=3)
# }
aic.ordered<-aic[rev(order(aic$aic.w)),]

# see AIC table
aic.ordered[,c(2,5,6,7,8,9,10)] <- round(aic.ordered[,c(2,5,6,7,8,9,10)], digits = 3)
print(aic.ordered)

########## RESULTS ##########

# Female egg volume is strongly positively correlated with winter dC values. Model with dC has AICw 0.692, r^2 = 0.402. Next closest model includes dC + initiation date, AICw = 0.128.  Higher dC values = larger eggs
# in contrast, running the same analysis with the male of the pair, dC comes out as top model at AICw = 0.33, but Null model is next at AICw = 0.16 and deltaAIC = 1.5.  Also, r squared for male dC model is low, r^2 = 0.169

#############################

out <- c("female clutch mass vs dC dN body size and initiation date")
# write AIC table output to .csv files
setwd("R/output/Chapter 4/COE MS") # working directory for output
write.csv(aic.ordered,paste("AIC output_",out,".csv",sep=""),row.names=FALSE)

# write model output to .csv files

sink(paste("model summaries_",out,".doc",sep=""))
for (i in 1:length(modellist)) {
print(paste("MODEL ", i, sep = ""))
print(ModelList[[i]])
print(summary(modeloutput[[i]]))
print(paste("###################################################################"))
}
sink()

#--------------------------------------------------
# calculate parameter likelihoods, weighted estimates, and unconditional SEs and CIs
#--------------------------------------------------

CoefMatrix[is.na(CoefMatrix)] <- 0
SEMatrix[is.na(SEMatrix)] <- 0

### parameter likelihoods

parlik<-rep(0,length(coefnames))

# if parameter j is included in the model i (ie. there is a non-zero value in the CoefMatrix for model i), add model i's AIC weight from the dataframe "aic" to the value of parlik[j], otherwise, add "0"

for (j in 1:length(coefnames)) {

for (i in 1:length(modeloutput)) {

if (CoefMatrix[i,j] != 0) {parlik[j] <- parlik[j] + aic[i,"aic.w"]} else {parlik[j] <- parlik[j]}

}
}

### weighted parameter estimates

weightedpar<-rep(NA,length(coefnames))

for (j in 1:length(coefnames)) {

weightedpar[j] <- sum(CoefMatrix[,j] * aic[,"aic.w"])

}

### unconditional SEs and CIs
# AICw * sqrt((SE^2) + (coef - weightedpar)^2 ))

uncondSE<-rep(NA,length=ncol(SEMatrix))


for (j in 1:length(coefnames)) {

uncondSE[j] <- sum(aic[,"aic.w"] * sqrt((SEMatrix[,j]^2) + (CoefMatrix[,j] - weightedpar[j])^2))
}

CI<-1.96 * uncondSE

lowerCI<-weightedpar - CI
upperCI<-weightedpar + CI

### parameter table

par.table <- data.frame(coefnames, parlik, weightedpar, uncondSE, lowerCI, upperCI)
par.table.new <- data.frame(coefnames,round(par.table[,2:6], digits=3)) # values rounded to 3 decimal places
print(par.table.new)

write.csv(par.table.new,paste("parameter estimates_",out,".csv",sep=""),row.names=FALSE)

setwd("/Users/samantha/Documents/Sam's Stuff")	# return working directory to your usual default

###-------------------------------------------------------------###
#	Analysis 2 - incubation mass  change vs dC, dN, Sex with RE(nest)
###-------------------------------------------------------------###

library(AICcmodavg)

# MODEL SELECTION ANALYSIS - with random intercept of Nest

dat <- with(coe2,data.frame(Response = WeightChange, Sex, dC, dN, fullnest), na.rm=TRUE)
dat <- na.omit(dat)

# set up blank matrices and lists

m1<-lme(Response~1, random = ~1 | fullnest, data=dat, method = "ML")
m2<-lme(Response~Sex + dC, random = ~1 | fullnest, data=dat, method = "ML")
m3<-lme(Response~Sex + dN, random = ~1 | fullnest, data=dat, method = "ML")
m4<-lme(Response~Sex + dC + dN, random = ~1 | fullnest, data=dat, method = "ML")
m5<-lme(Response~Sex, random = ~1 | fullnest, data=dat, method = "ML")

modellist<-list(m1,m2,m3,m4,m5)
modnames<-c("Intercept","Sex + dC","Sex + dN", "Sex + dC + dN","Sex")

dNtablesort<-aictab.lme(modellist,modnames,sort=TRUE)
dNtablesort

dNtable <- aictab.lme(modellist,modnames,sort=FALSE)
dNtable

write.csv(dNtablesort,"R/output/Chapter 4/COE MS/incubation mass change_aic.csv",row.names=FALSE)

########## RESULTS ##########

# No relationship between incubation mass change and dC/dN

#############################

###-------------------------------------------------------------###
#	Analysis 3 - early incubation body condition vs dC, dN, Sex
###-------------------------------------------------------------###

# change Response to EarlyWeight or LateWeight for early vs late incubation mass
dat <- with(coe2, data.frame(Response = LateWeight, dC, dN, Sex, fullnest, BodySize))

modellist <- list()  # set up blank list to hold model formulas
modellist[[1]] <- formula(Response ~ 1)
modellist[[2]] <- formula(Response ~ BodySize)
modellist[[3]] <- formula(Response ~ BodySize + dC)
modellist[[4]] <- formula(Response ~ BodySize + dN)
modellist[[5]] <- formula(Response ~ BodySize + dC + dN)

globalmodel <- lm(Response ~ dC + dN + BodySize, data=dat)

ModelList<-modellist
coefnames <- names(globalmodel$coef)
modeloutput<-list()	# blank model output list

# AIC function
# a function to produce relevant calculations for an AIC table that is used in conjunction with a list containing the model names
calculate.AIC<-function(aictable,modellist) {
modelnames<-modellist
delta.aic<-aictable$AIC.c-min(aictable$AIC.c)
lik.aic<-exp(-delta.aic/2)
aic.w<-lik.aic/(sum(lik.aic))
aic.table<-data.frame(modelnames,AIC.table,delta.aic,lik.aic,aic.w)
}

# blank Coefficient Matrix
CoefMatrix <- matrix(NA, nrow=length(ModelList), ncol=length(coefnames), dimnames = list(c(1:length(ModelList)), coefnames))

# blank SE Matrix
SEMatrix <- matrix(NA, nrow=length(ModelList), ncol=length(coefnames), dimnames = list(c(1:length(ModelList)), coefnames))

# blank AIC matrix
AIC.table<-matrix(NA, nrow=length(ModelList), ncol=6, dimnames = list(c(1:length(ModelList)), c("r.squared","n.obs","df","-2loglik","AIC","AIC.c")))

# loop for calculating model output, filling CoefMatrix and SEMatrix

for(i in 1:length(ModelList)){

modeloutput[[i]]<-lm(ModelList[[i]], data=dat)
m<-modeloutput[[i]]

# fill in row "i" of the AIC table
aic<-AIC(m)
df<-length(m$coef) + 1 # +1 for the variance
n.obs<-nrow(dat)
AIC.table[i,"r.squared"]<-summary(m)$adj.r.squared
AIC.table[i,"n.obs"]<-n.obs
AIC.table[i,"df"]<-df
AIC.table[i,"-2loglik"]<- -2*logLik(m)
AIC.table[i,"AIC"]<-aic
aic.c<-aic + (2*df*(df+1))/(n.obs-df-1)
AIC.table[i,"AIC.c"]<-aic.c

# fill in row "i", column "j" of the Coefficient Matrix
for (j in 1:length(coefnames)) {
	
if (coefnames[j] == "(Intercept)" & length(names(m$coef)) == 1) {CoefMatrix[i,j] <- summary(m)$coef[coefnames[j],1]} else {
CoefMatrix[i,j] <- summary(m)$coef[,1][coefnames[j]]
# CoefMatrix[i,j] <- summary(m)$coef[coefnames[j],1] # gives subscript out of bounds error if coefnames[j] is not one of the parameters in the model
# CoefMatrix[i,j] <- summary(m)$coef[,1][coefnames[j]] gives NA if model is the NULL model (only 1 parameter, the intercept) because in that case, summary(m)$coef is a vector (1-dimension only) rather than a dataframe (2-dimensions)
# code a special case for the NULL (intercept only) model, all other models with other parameters get the "else" statement
} # close if else statement

# fill in row "i" of the Standard Error Matrix
if (coefnames[j] == "(Intercept)" & length(names(m$coef)) == 1) {SEMatrix[i,j] <- summary(m)$coef[coefnames[j],2]} else {
SEMatrix[i,j] <- summary(m)$coef[,2][coefnames[j]]
} # close if else statement

} # close inner "j" loop

} # close outer "i" loop

# calculate deltaAICs, likelihoods, and AIC weights
AIC.output<-data.frame(AIC.table)
aic<-calculate.AIC(AIC.output,as.character(ModelList))
#for (i in 4:9) {
# aic[,i]<-round(aic[,i], digits=3)
# }
aic.ordered<-aic[rev(order(aic$aic.w)),]

# see AIC table
aic.ordered[,c(2,5,6,7,8,9,10)] <- round(aic.ordered[,c(2,5,6,7,8,9,10)], digits = 3)
print(aic.ordered)

########## RESULTS ##########

# No relationship between early incuation mass or late incubation mass and dC/dN
# Mass is mainly a function of principal component of body size (wing and tarsus)

#############################

out <- c("female clutch mass vs dC dN body size and initiation date")
# write AIC table output to .csv files
setwd("R/output/Chapter 4/COE MS") # working directory for output
write.csv(aic.ordered,paste("AIC output_",out,".csv",sep=""),row.names=FALSE)

# write model output to .csv files

sink(paste("model summaries_",out,".doc",sep=""))
for (i in 1:length(modellist)) {
print(paste("MODEL ", i, sep = ""))
print(ModelList[[i]])
print(summary(modeloutput[[i]]))
print(paste("###################################################################"))
}
sink()
