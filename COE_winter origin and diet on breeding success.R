##########################################################################
#
# Carry-over effects - relationship of winter origin and diet/habitat on nest success
# Samantha Franks
# 16 July 2013
# 23 Oct 2013
#
#	with revised assignments using sex-specific priors and corrected isotope data Oct 2013
#
##########################################################################

rm(list=ls()) #clear all previous

library(lattice) # plotting functions
library(NCStats)
#library(faraway)
#library(nlme)
#library(lme4)
#library(MuMIn)
library(AICcmodavg)

#------------------------------------------------------------
# AIC function
#------------------------------------------------------------

# a function to produce relevant calculations for an AIC table that is used in conjunction with a list containing the model names
calculate.AIC.2<-function(aictable,modellist) {
modelnames<-modellist
delta.aic<-aictable$AIC.c-min(aictable$AIC.c)
lik.aic<-exp(-delta.aic/2)
aic.w<-lik.aic/(sum(lik.aic))
aic.table<-data.frame(modelnames,AIC.table,delta.aic,lik.aic,aic.w)
print(aic.table)
}


#------------------------------------------
# Load data
#------------------------------------------

dd<-read.csv("R/sandpiper/merged breeding nest banding isotope and assignment data_20120108.csv", header=TRUE)
str(dd)
wesa<-Subset(dd,Species=="WESA") # WESA only
str(wesa)
head(wesa)

# remove incomplete records (NAs)
t1<-na.omit(wesa)
t2<-data.frame(t1)

# remove duplicate records
t3<-unique(t2)

# adults only
complete<-Subset(t3,Revised.Age=="ASY")
str(complete)
colnames(complete)<-c("Site","Year","Latitude","Longitude","Band","fullnest","INIT.JULIAN","Species","Recap.btwn.years","Revised.Age","Age","Sex","Wing.mm","Culmen.mm","Tarsus.mm","Weight.g","dN.feathers","dC.feathers","dD.feathers","Origin")

# write datafile of complete records

write.csv(complete,"R/sandpiper/merged breeding nest banding and isotope data_complete records.csv",row.names=F)

# identification of unknowns - changed in all master data sheets to females
# 2008_228156951 is female based on pair dynamics (other bird is definitely male)
# 2010_140165960 is female based on culmen (mistakenly coded as U)
# 2010_140165965 is female based on other morphometrics besides culmen (very large wing and tarsus!)
# known<-Subset(complete,Sex!="U") # subset out males and females

# set Year as a factor
complete$Year<-factor(complete$Year)	

# sample size by SITE and YEAR
n<-aggregate(list(N=complete$fullnest),list(Site=complete$Site,Year=complete$Year),length)
n<-reshape(n,idvar=c("Site"),timevar="Year",direction="wide")
n[is.na(n)] <- 0
n

#-------------------------------------------------
# remove renests for each site and year (late outlier Julian dates
#-------------------------------------------------

### NOME
dotplot(INIT.JULIAN~Year,data=complete[complete$Site=="Nome",]) # 2 late outliers in 2008 (>145), 1 late outlier in 2009 (>157), 4 late outliers in 2010 (>155)


# 2008
# using observations from Nome in 2008, reverse order them based on INIT.JULIAN
# structure for ordering is dd[order(dd$columnname),] - in this case, dd is complete, where Site = Nome and Year = 2008
head(complete[complete$Site=="Nome"  & complete$Year=="2008",][rev(order(complete[complete$Site=="Nome" & complete$Year=="2008",]$INIT.JULIAN)),])	# dates are 154 and 147
firstnests<-Subset(complete,Band!="2008_228156938" & Band!="2008_228156948")

# 2009
head(complete[complete$Site=="Nome"  & complete$Year=="2009",][rev(order(complete[complete$Site=="Nome" & complete$Year=="2009",]$INIT.JULIAN)),])	# dates are 160
firstnests<-Subset(firstnests,Band!="2009_197100520")

# 2010
head(complete[complete$Site=="Nome"  & complete$Year=="2010",][rev(order(complete[complete$Site=="Nome" & complete$Year=="2010",]$INIT.JULIAN)),])	# dates are 158-163
firstnests<-Subset(firstnests,fullnest!="10nome_10wesabs13" & fullnest!="10nome_10wesanh17" & fullnest!="10nome_10wesaij10" & fullnest!="10nome_10wesaij11")

dotplot(INIT.JULIAN~Year,data=firstnests[firstnests$Site=="Nome",])

### Espenberg
dotplot(INIT.JULIAN~Year,data=complete[complete$Site=="Espenberg",]) # 0 outliers

### Barrow
dotplot(INIT.JULIAN~Year,data=complete[complete$Site=="Barrow",]) # 0 outliers, sample size from 2009 is too small (4)

### Kotzebue
dotplot(INIT.JULIAN~Year,data=complete[complete$Site=="Kotzebue",]) # 0 outliers, small sample size

### Russia
dotplot(INIT.JULIAN~Year,data=complete[complete$Site=="Russia",]) # 0 outliers

### YKD
dotplot(INIT.JULIAN~Year,data=complete[complete$Site=="YKD",]) # potentially outliers at day=160 or later, but will leave in for now

### Cape Krusenstern
dotplot(INIT.JULIAN~Year,data=complete[complete$Site=="Cape Krusenstern",]) # 0 outliers

# sample size by SITE and YEAR
n<-aggregate(list(N=firstnests$fullnest),list(Site=firstnests$Site,Year=firstnests$Year),length)
n<-reshape(n,idvar=c("Site"),timevar="Year",direction="wide")
n[is.na(n)] <- 0
n

#--------------------------------------
# Add variable "Distance from breeding site" to dataframe
#--------------------------------------

distance <- read.csv("R/sandpiper/distance breeding to wintering regions.csv",header=TRUE)

dist <- rep(NA,nrow(firstnests))

worigin <- levels(firstnests$Origin)
borigin <- levels(firstnests$Site)


for (i in 1:nrow(firstnests)) {
	
	for (j in 1:length(worigin)){
		for (k in 1:length(borigin)){
			if (firstnests[i,"Origin"]==worigin[j] & firstnests[i,"Site"]==borigin[k]) {dist[i] <- distance$Distance[distance$Origin==worigin[j] & distance$Breeding==borigin[k]]} else {k <- k+1}
			}
			j <- j+1
			}
}

append <- data.frame(firstnests,Distance=dist)

# confidence in assignment of birds with nest initiation dates and isotope values
assign<-read.csv("R/sandpiper/summary all breeding assignments_with sex specific priors_20111228.csv",header=TRUE)
assign.id<-paste(assign$Year,"_",assign$Band,sep="")
assign<-cbind(assign.id,assign)

merge1 <- merge(append,assign,by.x="Band",by.y="assign.id")
correctrate <- data.frame(merge1[,c(1:21,38:42)])
summary(correctrate$F.C50)

# > 50% = 168/172 = 98%
# > 70% = 128/172 = 74%

#--------------------------------------
# Remove outlier dC values
#--------------------------------------

temp <- Subset(append, dC.feathers > -20)

#--------------------------------------
# Add a YearSite term which combines year and site into one categorical variable
#
# Add variable standardized nest initiation date to dataset
# Standardized in relation to the mean for each site/year
#
#--------------------------------------


YearSite <- as.factor(paste(temp$Year,temp$Site,sep=""))
newnests <- data.frame(temp,YearSite)

meannests <- aggregate(list(MeanDate=newnests$INIT.JULIAN),list(YearSite=newnests$YearSite),mean)

position <- match(newnests$YearSite,meannests$YearSite)

meannest <- meannests$MeanDate[position]

standard <- newnests$INIT.JULIAN-meannest

stan.nests <- data.frame(newnests,meannest,standard)

# test for collinearity between dN and Distance

library(car)
m <- lm(standard~dN.feathers + Distance, data=stan.nests)
vif(m)

#--------------------------------------
#
# Add variable standardized dN to dataset
# Standardized in relation to the mean for each BREEDING site/year
#
#--------------------------------------


#meannestsdN <- aggregate(list(MeandN=stan.nests$dN.feathers),list(YearSite=stan.nests$YearSite),mean)

#position <- match(stan.nests$YearSite,meannestsdN$YearSite)

#meannestdN <- meannestsdN$MeandN[position]

#standard.dN <- stan.nests$dN.feathers-meannestdN

#stan.nests.dN <- data.frame(stan.nests,meannestdN,standard)


#--------------------------------------
# Standardized dN by WINTER site
#--------------------------------------

# dN for each individual standardized relative to the mean for its region of WINTER origin

meannestsdN <- aggregate(list(MeandN=stan.nests$dN.feathers),list(Winter=stan.nests$Origin),mean)

position <- match(stan.nests$Origin,meannestsdN$Winter)

meannestdN <- meannestsdN$MeandN[position]

standard.dN <- stan.nests$dN.feathers-meannestdN

stan.nests.dN <- data.frame(stan.nests,meannestdN,standard.dN)

#--------------------------------------
# Standardized dC by WINTER site
#--------------------------------------

# dC for each individual standardized relative to the mean for its region of WINTER origin

meannestsdC <- aggregate(list(MeandC=stan.nests$dC.feathers),list(Winter=stan.nests$Origin),mean)

position <- match(stan.nests$Origin,meannestsdC$Winter)

meannestdC <- meannestsdC$MeandC[position]

standard.dC <- stan.nests$dC.feathers-meannestdC

stan.nests.dNdC <- data.frame(stan.nests.dN,meannestdC,standard.dC)


###############################################################
#--------------------------------------
# Nest initiation by Year and Site - Figure 2a
#--------------------------------------

m <- aov(INIT.JULIAN ~ Year + Site, data=stan.nests)
anova(m)
TukeyHSD(m,"Year")
TukeyHSD(m,"Site")
par(mar=c(5,10,4,2))
plot(TukeyHSD(m,"Year"),las=1,cex.axis=0.6)
par(mar=c(5,10,4,2))
plot(TukeyHSD(m,"Site"),las=1,cex.axis=0.6)

m <- aov(INIT.JULIAN ~ YearSite, data=stan.nests)
anova(m)
TukeyHSD(m)
par(mar=c(5,10,4,2))
plot(TukeyHSD(m),las=1,cex.axis=0.6)

m <- lm(INIT.JULIAN ~ Year + Site, data=stan.nests)
summary(m)

#-------------------------------------
# dataset of unique nests only to calculate means and sd of initiation dates and to draw boxplot
#-------------------------------------

nestsonly <- with(stan.nests,data.frame(fullnest,INIT.JULIAN,YearSite))

unique.nests <- unique(nestsonly)

###Stat = column which has the value (eg. C or N value)
###Category = Site (or Sex, Spp, etc)
catmeans<-function(Stat, Category){
    aggregate(list(Mean=Stat),list(Region=Category), mean,na.rm = TRUE)}
catsd<-function(Stat, Category){
    aggregate(list(SD=Stat),list(Region=Category), sd,na.rm = TRUE)}

init.mean<-catmeans(unique.nests$INIT.JULIAN, unique.nests$YearSite)
init.sd<-catsd(unique.nests$INIT.JULIAN, unique.nests$YearSite)

init<-data.frame(Region=init.mean$Region,Mean=init.mean$Mean,SD=init.sd$SD)


tiff("R/output/Chapter 4/Figure 2.tiff",res=150,height=21,width=21,units="cm")
par(mfrow=c(2,2), mar=c(5,4,2,1))

#tiff("R/output/Chapter 4/Figure 2a.tiff",res=150,height=15,width=15,units="cm")
boxplot(INIT.JULIAN~YearSite,data=unique.nests,xaxt="n",ylim=c(133,180),ylab="",las=1)
axis(1,at=c(1:10),labels=F)
title(ylab="Nest initiation date", line=3)
samplesize <- c(9,8,3,7,19,11,27,8,7,21)
text(1:10,133,samplesize,cex=0.8)
text(1:10,par("usr")[3]-1.7, srt = 35, adj = 1,labels = levels(stan.nests$YearSite),xpd = TRUE,cex=0.8)
mtext("a)", side=3, line = 0.25, adj = 0, cex=1.2)
#dev.off()

#--------------------------------------
# dC and dN by Year and Site - Figure 2b and 2c
#--------------------------------------

#tiff("R/output/Chapter 4/Figure 2b.tiff",res=150,height=15,width=15,units="cm")
boxplot(dC.feathers~YearSite,data=stan.nests,xaxt="n",ylab="",las=1,ylim=c(-21,-5))
axis(1,at=c(1:10),labels=F)
title(ylab=(expression(bold(paste(delta^13,"C"," (‰)")))), line=2.5)
samplesize <- summary(stan.nests$YearSite)
text(1:10,-21,samplesize,cex=0.8)
text(1:10,par("usr")[3]-0.6, srt = 35, adj = 1,labels = levels(stan.nests$YearSite),xpd = TRUE,cex=0.8)
mtext("b)", side=3, line = 0.25, adj = 0, cex=1.2)
#dev.off()

#tiff("R/output/Chapter 4/Figure 2c.tiff",res=150,height=15,width=15,units="cm")
boxplot(dN.feathers~YearSite,data=stan.nests,xaxt="n",ylab="",las=1,ylim=c(6,20))
axis(1,at=c(1:10),labels=F)
title(ylab=(expression(bold(paste(delta^15,"N"," (‰)")))), line=2.5)
samplesize <- summary(stan.nests$YearSite)
text(1:10,6,samplesize,cex=0.8)
text(1:10,par("usr")[3]-0.6, srt = 35, adj = 1,labels = levels(stan.nests$YearSite),xpd = TRUE,cex=0.8)
mtext("c)", side=3, line = 0.25, adj = 0, cex=1.2)
dev.off()

#--------------------------------------
# Nest initiation by Year and Site
#--------------------------------------

minnests <- aggregate(list(MinDate=newnests$INIT.JULIAN),list(YearSite=newnests$YearSite),min)
maxnests <- aggregate(list(MaxDate=newnests$INIT.JULIAN),list(YearSite=newnests$YearSite),max)
range1 <- merge(meannests,minnests,by.x="YearSite")
range2 <- merge(range1,maxnests,by.x="YearSite")
diffnests <- with(range2,MaxDate-MinDate)
diffnests2 <- with(range2,MaxDate-MeanDate)
diffnests3 <- with(range2, MeanDate-MinDate)
range <- cbind(range2,diffnests, diffnests2,diffnests3)
range	

WeekEarly <- stan.nests$meannest-5
WeekLate <- stan.nests$meannest+5

hist(stan.nests$standard)
withinweek <- stan.nests$INIT.JULIAN > WeekEarly & stan.nests$INIT.JULIAN < WeekLate
summary(withinweek)
# 103 out of 172 (60%) individuals nested within 1 week of the site mean initiation date
# 130/172 (75%) individuals nested within a 10 day period around the mean initiation date
# 143/172 (83%) individuals nested within a 2 week period around the mean initiation date

#--------------------------------------
# Male vs Female migration distance
#--------------------------------------

pairs <- with(stan.nests,data.frame(fullnest,Sex,Origin,Distance))

newpairs <- reshape(pairs,v.names=c("Origin","Distance"),idvar="fullnest", timevar = c("Sex"), direction="wide")

completepairs <- na.omit(newpairs)

par(mar=c(5,4.5,4,2))
plot(Distance.M~Distance.F, data=completepairs, pch=16, las=1, xlab="", ylab="", xlim=c(3800,10100),ylim=c(3800,10100))
title(xlab="Female migration distance (km)", line=2.5)
title(ylab="Male migration distance (km)", line=3.5)

#--------------------------------------
# Standardized dN by winter site
#--------------------------------------

# dN for each individual standardized relative to the mean for its region of winter origin

##############################################
##############################################		MODELS
##############################################

#--------------------------------------
# Models and Analysis - effect of dN, dC and MIGRATION DISTANCE on nest initiation - with standardized initiation dates and standardized dN and dC by WINTER site, YearSite
#--------------------------------------


iso.name <- c("standardized by winter site") # change to desired isotope

dat<-with(stan.nests.dNdC,data.frame(Response=standard,YearSite,Sex,Origin,Distance, dN = standard.dN, dC=standard.dC))


female <- Subset(dat,Sex=="F")
male <- Subset(dat,Sex=="M")


### MALES

dat <- male

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Distance)
modellist[[2]]<-formula(Response~dN)
modellist[[3]]<-formula(Response~Distance + dN)
modellist[[4]]<-formula(Response~Distance*dN)
modellist[[5]]<-formula(Response~NULL)
modellist[[6]]<-formula(Response~dC)
modellist[[7]]<-formula(Response~Distance + dC)
modellist[[8]]<-formula(Response~Distance*dC)
modellist[[9]]<-formula(Response~dC + dN)
modellist[[10]]<-formula(Response~Distance + dC + dN)
modellist[[11]]<-formula(Response~Distance*dN*dC)
modellist[[12]]<-formula(Response~dC*dN)


# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_males",sep="")
wd <- c("R/output/Chapter 4/breeding/males")
source("R/scripts/AIC script revised.r")




### FEMALES

dat <- female

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Distance)
modellist[[2]]<-formula(Response~dN)
modellist[[3]]<-formula(Response~Distance + dN)
modellist[[4]]<-formula(Response~Distance*dN)
modellist[[5]]<-formula(Response~NULL)
modellist[[6]]<-formula(Response~dC)
modellist[[7]]<-formula(Response~Distance + dC)
modellist[[8]]<-formula(Response~Distance*dC)
modellist[[9]]<-formula(Response~dC + dN)
modellist[[10]]<-formula(Response~Distance + dC + dN)
modellist[[11]]<-formula(Response~Distance*dN*dC)
modellist[[12]]<-formula(Response~dC*dN)


# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_females",sep="")
wd <- c("R/output/Chapter 4/breeding/females")
source("R/scripts/AIC script revised.r")

##############################################
##############################################
##############################################


#--------------------------------------
# Models and Analysis - effect of dN, dC and MIGRATION DISTANCE on nest initiation - with standardized initiation dates and standardized dN by breeding site, no YearSite
#--------------------------------------


### Candidate Model set; Response = INIT.JULIAN
# Distance
# dN feathers
# Distance + dN feathers
# Distance*dN feathers

# create dataset with only the variables of interest from the main dataset "firstnests"

iso.name <- c("standard dN") # change to desired isotope

dat<-with(stan.nests.dN,data.frame(Response=standard,YearSite,Sex,Origin,Distance, dN = standard.dN, dC=dC.feathers))


female <- Subset(dat,Sex=="F")
male <- Subset(dat,Sex=="M")


### MALES

dat <- male

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Distance)
modellist[[2]]<-formula(Response~dN)
modellist[[3]]<-formula(Response~Distance + dN)
modellist[[4]]<-formula(Response~Distance*dN)
modellist[[5]]<-formula(Response~NULL)
modellist[[6]]<-formula(Response~dC)
modellist[[7]]<-formula(Response~Distance + dC)
modellist[[8]]<-formula(Response~Distance*dC)
modellist[[9]]<-formula(Response~dC + dN)
modellist[[10]]<-formula(Response~Distance + dC + dN)
modellist[[11]]<-formula(Response~Distance*dN*dC)
modellist[[12]]<-formula(Response~dC*dN)


# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_males",sep="")
wd <- c("R/output/Chapter 4/breeding/males")
source("R/scripts/AIC script revised.r")




### FEMALES

dat <- female

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Distance)
modellist[[2]]<-formula(Response~dN)
modellist[[3]]<-formula(Response~Distance + dN)
modellist[[4]]<-formula(Response~Distance*dN)
modellist[[5]]<-formula(Response~NULL)
modellist[[6]]<-formula(Response~dC)
modellist[[7]]<-formula(Response~Distance + dC)
modellist[[8]]<-formula(Response~Distance*dC)
modellist[[9]]<-formula(Response~dC + dN)
modellist[[10]]<-formula(Response~Distance + dC + dN)
modellist[[11]]<-formula(Response~Distance*dN*dC)
modellist[[12]]<-formula(Response~dC*dN)


# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_females",sep="")
wd <- c("R/output/Chapter 4/breeding/females")
source("R/scripts/AIC script revised.r")

##############################################
##############################################
##############################################

#--------------------------------------
# Models and Analysis - effect of dN, dC and MIGRATION DISTANCE on nest initiation - test for subset of sites (largest sample sizes), with site x isotope interaction
#--------------------------------------

newsites <- Subset(stan.nests, YearSite=="2009Nome" | YearSite=="2010Nome" | YearSite=="2009YKD")

iso.name <- c("N C distance") # change to desired isotope

dat<-with(newsites,data.frame(Response=INIT.JULIAN,YearSite,Sex,Origin,Distance, dN = dN.feathers, dC=dC.feathers))

female <- Subset(dat,Sex=="F")
male <- Subset(dat,Sex=="M")



### MALES

dat <- male

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Distance + YearSite)
modellist[[2]]<-formula(Response~dN + YearSite)
modellist[[3]]<-formula(Response~Distance + dN + YearSite)
modellist[[4]]<-formula(Response~Distance*dN*YearSite)
modellist[[5]]<-formula(Response~NULL)
modellist[[6]]<-formula(Response~dC + YearSite)
modellist[[7]]<-formula(Response~Distance + dC + YearSite)
modellist[[8]]<-formula(Response~Distance*dC*YearSite)
modellist[[9]]<-formula(Response~dC + dN + YearSite)
modellist[[10]]<-formula(Response~Distance + dC + dN + YearSite)
modellist[[11]]<-formula(Response~Distance*dN*dC*YearSite)
modellist[[12]]<-formula(Response~Distance*YearSite)
modellist[[13]]<-formula(Response~dN*YearSite)
modellist[[14]]<-formula(Response~dC*YearSite)


# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC*YearSite,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_males_sitesubset",sep="")
wd <- c("R/output/Chapter 4/breeding/males")
source("R/scripts/AIC script revised.r")


### FEMALES

dat <- female

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Distance + YearSite)
modellist[[2]]<-formula(Response~dN + YearSite)
modellist[[3]]<-formula(Response~Distance + dN + YearSite)
modellist[[4]]<-formula(Response~Distance*dN*YearSite)
modellist[[5]]<-formula(Response~NULL)
modellist[[6]]<-formula(Response~dC + YearSite)
modellist[[7]]<-formula(Response~Distance + dC + YearSite)
modellist[[8]]<-formula(Response~Distance*dC*YearSite)
modellist[[9]]<-formula(Response~dC + dN + YearSite)
modellist[[10]]<-formula(Response~Distance + dC + dN + YearSite)
modellist[[11]]<-formula(Response~Distance*dN*dC*YearSite)
modellist[[12]]<-formula(Response~Distance*YearSite)
modellist[[13]]<-formula(Response~dN*YearSite)
modellist[[14]]<-formula(Response~dC*YearSite)


# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC*YearSite,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_females_sitesubset",sep="")
wd <- c("R/output/Chapter 4/breeding/females")
source("R/scripts/AIC script revised.r")


##############################################
##############################################
##############################################


#--------------------------------------
# Models and Analysis - effect of dN, dC and MIGRATION DISTANCE on nest initiation - with standardized initiation dates, no YearSite
#
#
#	**** this is the analysis used in thesis chapter #5
#--------------------------------------


### Candidate Model set; Response = INIT.JULIAN
# Distance
# dN feathers
# Distance + dN feathers
# Distance*dN feathers

# create dataset with only the variables of interest from the main dataset "firstnests"

iso.name <- c("N C distance_20120314") # change to desired isotope

dat<-with(stan.nests,data.frame(Response=standard,YearSite,Sex,Origin,Distance, dN = dN.feathers, dC=dC.feathers))

dat <- data.frame(dat,YearSite)


female <- Subset(dat,Sex=="F")
male <- Subset(dat,Sex=="M")


### MALES

dat <- male

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Distance)
modellist[[2]]<-formula(Response~dN)
modellist[[3]]<-formula(Response~Distance + dN)
modellist[[4]]<-formula(Response~Distance*dN)
modellist[[5]]<-formula(Response~NULL)
modellist[[6]]<-formula(Response~dC)
modellist[[7]]<-formula(Response~Distance + dC)
modellist[[8]]<-formula(Response~Distance*dC)
modellist[[9]]<-formula(Response~dC + dN)
modellist[[10]]<-formula(Response~Distance + dC + dN)
modellist[[11]]<-formula(Response~Distance*dN*dC)
modellist[[12]]<-formula(Response~dC*dN)

# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_males",sep="")
wd <- c("R/output/Chapter 4/breeding/males")
source("R/scripts/AIC script revised.r")

#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------


p <- par.table$weightedpar

median.dist <- median(dat$Distance)
median.dC <- median(dat$dC)
median.dN <- median(dat$dN)
int <- p[1]
N.max <- max(dat$dN)
N.min <- min(dat$dN)
D.max <- max(dat$Distance)
D.min <- min(dat$Distance)
C.max <- max(dat$dC)
C.min <- min(dat$dC)

#parameters p[5:13] are YearSite parameters
#int #dist #dN #dC #dist:dN #dist:dC #dN:dC #dist:dN:dC
line.dN <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*median.dist*x) + (p[6]*median.dist*median.dC) + (p[7]*x*median.dC) + (p[8]*median.dist*x*median.dC))

line.D <- function(x)(int + (p[2]*x) + (p[3]*median.dN) + (p[4]*median.dC) + (p[5]*x*median.dN) + (p[6]*x*median.dC) + (p[7]*median.dN*median.dC) + (p[8]*x*median.dN*median.dC))

line.dC <- function(x)(int + (p[2]*median.dist) + (p[3]*median.dN) + (p[4]*x) + (p[5]*median.dist*median.dN) + (p[6]*median.dist*x) + (p[7]*median.dN*x) + (p[8]*median.dist*median.dN*x))

tiff("R/output/Chapter 4/Figure 3.tiff",res=150,height=15,width=25,units="cm")
par(mfrow=c(2,3),mar=c(5,5,2,1))

plot(Response~Distance, data=dat, type="p", pch=16, ylab = "Standardized nest initiation date", xlab="", las=1, cex.axis=1.2, cex.lab=1.2)
lines(seq(D.min,D.max,100),line.D(seq(D.min,D.max,100)),col="black",lty=3,lwd=2)
mtext(c("a) Males"),side=3,line=0.5,font=2,cex=0.8,adj=0,las=1)

plot(Response~dC, data=dat, type="p", pch=16, las=1, ylab="", xlab="", cex.axis=1.2, cex.lab=1.2)
lines(seq(C.min,C.max,0.1),line.dC(seq(C.min,C.max,0.1)),col="black",lty=3,lwd=2)
mtext(c("b)"),side=3,line=0.5,font=2,cex=0.8,adj=0,las=1)

plot(Response~dN, data=dat, type="p", pch=16, las=1, ylab="", xlab="", cex.axis=1.2, cex.lab=1.2)
lines(seq(N.min,N.max,0.1),line.dN(seq(N.min,N.max,0.1)),col="black",lty=3,lwd=2)
mtext(c("c)"),side=3,line=0.5,font=2,cex=0.8,adj=0,las=1)



### FEMALES

dat <- female

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Distance)
modellist[[2]]<-formula(Response~dN)
modellist[[3]]<-formula(Response~Distance + dN)
modellist[[4]]<-formula(Response~Distance*dN)
modellist[[5]]<-formula(Response~NULL)
modellist[[6]]<-formula(Response~dC)
modellist[[7]]<-formula(Response~Distance + dC)
modellist[[8]]<-formula(Response~Distance*dC)
modellist[[9]]<-formula(Response~dC + dN)
modellist[[10]]<-formula(Response~Distance + dC + dN)
modellist[[11]]<-formula(Response~Distance*dN*dC)
modellist[[12]]<-formula(Response~dC*dN)

# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_females",sep="")
wd <- c("R/output/Chapter 4/breeding/females")
source("R/scripts/AIC script revised.r")

#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------


p <- par.table$weightedpar

median.dist <- median(dat$Distance)
median.dC <- median(dat$dC)
median.dN <- median(dat$dN)
int <- p[1]
N.max <- max(dat$dN)
N.min <- min(dat$dN)
D.max <- max(dat$Distance)
D.min <- min(dat$Distance)
C.max <- max(dat$dC)
C.min <- min(dat$dC)

#parameters p[5:13] are YearSite parameters
#int #dist #dN #dC #site(change to appropriate parameter depending on site) #dist:dN #dist:dC #dN:dC #dist:dN:dC
line.dN <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*median.dist*x) + (p[6]*median.dist*median.dC) + (p[7]*x*median.dC) + (p[8]*median.dist*x*median.dC))

line.D <- function(x)(int + (p[2]*x) + (p[3]*median.dN) + (p[4]*median.dC) + (p[5]*x*median.dN) + (p[6]*x*median.dC) + (p[7]*median.dN*median.dC) + (p[8]*x*median.dN*median.dC))

line.dC <- function(x)(int + (p[2]*median.dist) + (p[3]*median.dN) + (p[4]*x) + (p[5]*median.dist*median.dN) + (p[6]*median.dist*x) + (p[7]*median.dN*x) + (p[8]*median.dist*median.dN*x))

plot(Response~Distance, data=dat, type="p", pch=16, ylab = "Standardized nest initiation date", xlab = "Migration distance (km)", las=1, cex.axis=1.2, cex.lab=1.2)
lines(seq(D.min,D.max,100),line.D(seq(D.min,D.max,100)),col="black",lty=3,lwd=2)
mtext(c("d) Females"),side=3,line=0.5,font=2,cex=0.8,adj=0,las=1)

plot(Response~dC, data=dat, type="p", pch=16, ylab="", xlab = (expression(bold(paste(delta^13,"C"," (‰)")))), las=1, cex.axis=1.2, cex.lab=1.2)
lines(seq(C.min,C.max,0.1),line.dC(seq(C.min,C.max,0.1)),col="black",lty=3,lwd=2)
mtext(c("e)"),side=3,line=0.5,font=2,cex=0.8,adj=0,las=1)

plot(Response~dN, data=dat, type="p", pch=16, ylab="", xlab = (expression(bold(paste(delta^15,"N"," (‰)")))), las=1, cex.axis=1.2, cex.lab=1.2)
lines(seq(N.min,N.max,0.1),line.dN(seq(N.min,N.max,0.1)),col="black",lty=3,lwd=2)
mtext(c("f)"),side=3,line=0.5,font=2,cex=0.8,adj=0,las=1)

dev.off()



#############################################################################################################
##############################################
##############################################

#--------------------------------------
# Models and Analysis - effect of dN, dC and MIGRATION DISTANCE on nest initiation - with YearSite as covariate
#--------------------------------------

# include YearSite as fixed effect covariate, no random effects

### Candidate Model set; Response = INIT.JULIAN
# Year + Site + Sex
# Distance + Year + Site + Sex
# dN feathers + Year + Site + Sex
# Distance + dN feathers + Year + Site + Sex
# Distance*dN feathers + Year + Site + Sex

# create dataset with only the variables of interest from the main dataset "firstnests"

iso.name <- c("N C distance") # change to desired isotope

dat<-with(stan.nests,data.frame(Response=INIT.JULIAN,YearSite,Sex,Origin,Distance, dN = dN.feathers, dC=dC.feathers))

dat <- data.frame(dat,YearSite)


female <- Subset(dat,Sex=="F")
male <- Subset(dat,Sex=="M")




### MALES

dat <- male

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~YearSite)
modellist[[2]]<-formula(Response~Distance + YearSite)
modellist[[3]]<-formula(Response~dN + YearSite)
modellist[[4]]<-formula(Response~Distance + dN + YearSite)
modellist[[5]]<-formula(Response~Distance*dN + YearSite)
modellist[[6]]<-formula(Response~NULL)
modellist[[7]]<-formula(Response~dC + YearSite)
modellist[[8]]<-formula(Response~Distance + dC + YearSite)
modellist[[9]]<-formula(Response~Distance*dC + YearSite)
modellist[[10]]<-formula(Response~dC + dN + YearSite)
modellist[[11]]<-formula(Response~Distance + dC + dN + YearSite)
modellist[[12]]<-formula(Response~Distance*dN*dC + YearSite)


# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC + YearSite,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_males",sep="")
wd <- c("R/output/Chapter 4/breeding/males")
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------


p <- par.table$weightedpar

median.dist <- median(dat$Distance)
median.dC <- median(dat$dC)
int <- p[1]
x.max <- max(dat$dN)
x.min <- min(dat$dN)

#parameters p[5:13] are YearSite parameters
#int #dist #dN #dC #site(change to appropriate parameter depending on site) #dist:dN #dist:dC #dN:dC #dist:dN:dC
line.ES08 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.NO08 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*1) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.BA09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*1) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.KO09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*1) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.NO09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*1) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.RU09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*1) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.YK09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*1) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.BA10 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*1) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.CK10 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*1) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.NO10 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*1) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

par(c(1,2))
plot(Response~dN, data=dat, type="n", ylab = "Nest initiation date", xlab = "dN")

lines(seq(x.min,x.max,1),line.ES08(seq(x.min,x.max,1)),col="black",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.NO08(seq(x.min,x.max,1)),col="grey",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.BA09(seq(x.min,x.max,1)),col="blue",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.KO09(seq(x.min,x.max,1)),col="red",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.NO09(seq(x.min,x.max,1)),col="orange",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.RU09(seq(x.min,x.max,1)),col="darkgreen",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.YK09(seq(x.min,x.max,1)),col="lightgreen",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.BA10(seq(x.min,x.max,1)),col="purple",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.CK10(seq(x.min,x.max,1)),col="pink",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.NO10(seq(x.min,x.max,1)),col="brown",lty=1,lwd=2)

#par(mar=c(5,4,4,10))
#plot(Response~dN,data=dat,col=c(1:7)[dat$Site],pch=c(15:17)[dat$Year])
#legend(20,150,levels(dat$Site),col=c(1:7),pch=16,xpd=TRUE)

#m1 <- modeloutput[[3]]

#xyplot(predict(m1, type="response")~dN, groups = YearSite, data=dat, type="l", col = rainbow(10))


############################################################
############################################################


### FEMALES

# set up blank matrices and lists

### FEMALES

dat <- female

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~YearSite)
modellist[[2]]<-formula(Response~Distance + YearSite)
modellist[[3]]<-formula(Response~dN + YearSite)
modellist[[4]]<-formula(Response~Distance + dN + YearSite)
modellist[[5]]<-formula(Response~Distance*dN + YearSite)
modellist[[6]]<-formula(Response~NULL)
modellist[[7]]<-formula(Response~dC + YearSite)
modellist[[8]]<-formula(Response~Distance + dC + YearSite)
modellist[[9]]<-formula(Response~Distance*dC + YearSite)
modellist[[10]]<-formula(Response~dC + dN + YearSite)
modellist[[11]]<-formula(Response~Distance + dC + dN + YearSite)
modellist[[12]]<-formula(Response~Distance*dN*dC + YearSite)


# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC + YearSite,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_females",sep="")
wd <- c("R/output/Chapter 4/breeding/females") # set desired directory for AIC output
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------

p <- par.table$weightedpar

median.dist <- median(dat$Distance)
median.dC <- median(dat$dC)
int <- p[1]
x.max <- max(dat$dN)
x.min <- min(dat$dN)

#parameters p[5:13] are YearSite parameters
#int #dist #dN #dC #site(change to appropriate parameter depending on site) #dist:dN #dist:dC #dN:dC #dist:dN:dC
line.ES08 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.NO08 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*1) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.BA09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*1) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.KO09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*1) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.NO09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*1) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.RU09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*1) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.YK09 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*1) + (p[11]*0) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.BA10 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*1) + (p[12]*0) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.CK10 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*1) + (p[13]*0) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))

line.NO10 <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*0) + (p[6]*0) + (p[7]*0) + (p[8]*0) + (p[9]*0) + (p[10]*0) + (p[11]*0) + (p[12]*0) + (p[13]*1) + (p[14]*median.dist*x) + (p[15]*median.dist*median.dC) + (p[16]*x*median.dC) + (p[17]*median.dist*x*median.dC))


plot(Response~dN, data=dat, type="n", ylab = "Nest initiation date", xlab = "dN")

lines(seq(x.min,x.max,1),line.ES08(seq(x.min,x.max,1)),col="black",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.NO08(seq(x.min,x.max,1)),col="grey",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.BA09(seq(x.min,x.max,1)),col="blue",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.KO09(seq(x.min,x.max,1)),col="red",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.NO09(seq(x.min,x.max,1)),col="orange",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.RU09(seq(x.min,x.max,1)),col="darkgreen",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.YK09(seq(x.min,x.max,1)),col="lightgreen",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.BA10(seq(x.min,x.max,1)),col="purple",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.CK10(seq(x.min,x.max,1)),col="pink",lty=1,lwd=2)
lines(seq(x.min,x.max,1),line.NO10(seq(x.min,x.max,1)),col="brown",lty=1,lwd=2)


#--------------------------------------
# Models and Analysis - effect of dN and ORIGIN on nest initiation
#--------------------------------------

# include Site, Sex, and Year as fixed effects, no random effects

### Candidate Model set; Response = INIT.JULIAN
# Year + Site + Sex
# Origin + Year + Site + Sex
# dN feathers + Year + Site + Sex
# Origin + dN feathers + Year + Site + Sex
# Origin*dN feathers + Year + Site + Sex

# create dataset with only the variables of interest from the main dataset "firstnests"

iso.name <- c("N") # change to desired isotope

dat<-with(append,data.frame(Response=INIT.JULIAN,Site,Year,Sex,Origin,Distance,isotope = dN.feathers))	# change "isotope" column to desired isotope

female <- Subset(dat,Sex=="F")
male <- Subset(dat,Sex=="M")

### POOLED SEXES

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Year + Site + Sex)
modellist[[2]]<-formula(Response~Origin + Year + Site + Sex)
modellist[[3]]<-formula(Response~isotope + Year + Site + Sex)
modellist[[4]]<-formula(Response~Origin + isotope + Year + Site + Sex)
modellist[[5]]<-formula(Response~Origin*isotope + Year + Site + Sex)
modellist[[6]]<-formula(Response~isotope*Sex + Year + Site)
modellist[[7]]<-formula(Response~Origin*Sex + Year + Site)
modellist[[8]]<-formula(Response~Year + Site)
modellist[[9]]<-formula(Response~isotope + Year + Site)
modellist[[10]]<-formula(Response~Origin + Year + Site)
modellist[[11]]<-formula(Response~Origin + isotope + Year + Site)
modellist[[12]]<-formula(Response~Origin*isotope + Year + Site)
modellist[[13]]<-formula(Response~Origin*isotope*Sex + Year + Site)
modellist[[14]]<-formula(Response~NULL)

# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Origin*isotope*Sex + Year + Site,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------

out <- paste(iso.name,"_pooled",sep="")
wd <- c("R/output/Chapter 4/breeding/")
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------

par(mar=c(5,4,4,10))
plot(Response~isotope,data=dat,col=c(1:7)[dat$Site],pch=c(15:17)[dat$Year])
legend(20,150,levels(dat$Site),col=c(1:7),pch=16,xpd=TRUE)

############################################################
############################################################


### MALES

dat <- male

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Year + Site)
modellist[[2]]<-formula(Response~Origin + Year + Site)
modellist[[3]]<-formula(Response~isotope + Year + Site)
modellist[[4]]<-formula(Response~Origin + isotope + Year + Site)
modellist[[5]]<-formula(Response~Origin*isotope + Year + Site)
modellist[[6]]<-formula(Response~NULL)

# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Origin*isotope + Year + Site,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_males",sep="")
wd <- c("R/output/Chapter 4/breeding/males")
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------

par(mar=c(5,4,4,10))
plot(Response~isotope,data=dat,col=c(1:7)[dat$Site],pch=c(15:17)[dat$Year])
legend(20,150,levels(dat$Site),col=c(1:7),pch=16,xpd=TRUE)



############################################################
############################################################


### FEMALES

# set up blank matrices and lists

### FEMALES

dat <- female

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Year + Site)
modellist[[2]]<-formula(Response~Origin + Year + Site)
modellist[[3]]<-formula(Response~isotope + Year + Site)
modellist[[4]]<-formula(Response~Origin + isotope + Year + Site)
modellist[[5]]<-formula(Response~Origin*isotope + Year + Site)
modellist[[6]]<-formula(Response~NULL)

# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Origin*isotope + Year + Site,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------

out <- paste(iso.name,"_females",sep="")
wd <- c("R/output/Chapter 4/breeding/females") # set desired directory for AIC output
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------

par(mar=c(5,4,4,10))
plot(Response~isotope,data=dat,col=c(1:7)[dat$Site],pch=c(15:17)[dat$Year])
legend(20,150,levels(dat$Site),col=c(1:7),pch=16,xpd=TRUE)

###############################################################

#--------------------------------------
# Models and Analysis - effect of dN and MIGRATION DISTANCE on nest initiation
#--------------------------------------

# include Site, Sex, and Year as fixed effects, no random effects

### Candidate Model set; Response = INIT.JULIAN
# Year + Site + Sex
# Distance + Year + Site + Sex
# dN feathers + Year + Site + Sex
# Distance + dN feathers + Year + Site + Sex
# Distance*dN feathers + Year + Site + Sex

# create dataset with only the variables of interest from the main dataset "firstnests"

iso.name <- c("N distance") # change to desired isotope

dat<-with(append,data.frame(Response=INIT.JULIAN,Site,Year,Sex,Origin,Distance,isotope = dN.feathers))	# change "isotope" column to desired isotope

female <- Subset(dat,Sex=="F")
male <- Subset(dat,Sex=="M")

### POOLED SEXES

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Year + Site + Sex)
modellist[[2]]<-formula(Response~Distance + Year + Site + Sex)
modellist[[3]]<-formula(Response~isotope + Year + Site + Sex)
modellist[[4]]<-formula(Response~Distance + isotope + Year + Site + Sex)
modellist[[5]]<-formula(Response~Distance*isotope + Year + Site + Sex)
modellist[[6]]<-formula(Response~isotope*Sex + Year + Site)
modellist[[7]]<-formula(Response~Distance*Sex + Year + Site)
modellist[[8]]<-formula(Response~Year + Site)
modellist[[9]]<-formula(Response~isotope + Year + Site)
modellist[[10]]<-formula(Response~Distance + Year + Site)
modellist[[11]]<-formula(Response~Distance + isotope + Year + Site)
modellist[[12]]<-formula(Response~Distance*isotope + Year + Site)
modellist[[13]]<-formula(Response~Distance*isotope*Sex + Year + Site)
modellist[[14]]<-formula(Response~NULL)

# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*isotope*Sex + Year + Site,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------

out <- paste(iso.name,"_pooled",sep="")
wd <- c("R/output/Chapter 4/breeding/")
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------

par(mar=c(5,4,4,10))
plot(Response~isotope,data=dat,col=c(1:7)[dat$Site],pch=c(15:17)[dat$Year])
legend(20,150,levels(dat$Site),col=c(1:7),pch=16,xpd=TRUE)

############################################################
############################################################


### MALES

dat <- male

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Year + Site)
modellist[[2]]<-formula(Response~Distance + Year + Site)
modellist[[3]]<-formula(Response~isotope + Year + Site)
modellist[[4]]<-formula(Response~Distance + isotope + Year + Site)
modellist[[5]]<-formula(Response~Distance*isotope + Year + Site)
modellist[[6]]<-formula(Response~NULL)

# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*isotope + Year + Site,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_males",sep="")
wd <- c("R/output/Chapter 4/breeding/males")
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------

par(mar=c(5,4,4,10))
plot(Response~isotope,data=dat,col=c(1:7)[dat$Site],pch=c(15:17)[dat$Year])
legend(20,150,levels(dat$Site),col=c(1:7),pch=16,xpd=TRUE)



############################################################
############################################################


### FEMALES

# set up blank matrices and lists

### FEMALES

dat <- female

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Year + Site)
modellist[[2]]<-formula(Response~Distance + Year + Site)
modellist[[3]]<-formula(Response~isotope + Year + Site)
modellist[[4]]<-formula(Response~Distance + isotope + Year + Site)
modellist[[5]]<-formula(Response~Distance*isotope + Year + Site)
modellist[[6]]<-formula(Response~NULL)

# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*isotope + Year + Site,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_females",sep="")
wd <- c("R/output/Chapter 4/breeding/females") # set desired directory for AIC output
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------

par(mar=c(5,4,4,10))
plot(Response~isotope,data=dat,col=c(1:7)[dat$Site],pch=c(15:17)[dat$Year])
legend(20,150,levels(dat$Site),col=c(1:7),pch=16,xpd=TRUE)

###############################################################

#--------------------------------------
# Models and Analysis - effect of dN, dC and MIGRATION DISTANCE on nest initiation
#--------------------------------------

# include Site, Sex, and Year as fixed effects, no random effects

### Candidate Model set; Response = INIT.JULIAN
# Year + Site + Sex
# Distance + Year + Site + Sex
# dN feathers + Year + Site + Sex
# Distance + dN feathers + Year + Site + Sex
# Distance*dN feathers + Year + Site + Sex

# create dataset with only the variables of interest from the main dataset "firstnests"

iso.name <- c("N C distance") # change to desired isotope

dat<-with(append,data.frame(Response=INIT.JULIAN,Site,Year,Sex,Origin,Distance, dN = dN.feathers, dC=dC.feathers))


female <- Subset(dat,Sex=="F")
male <- Subset(dat,Sex=="M")


############################################################
############################################################


### MALES

dat <- male

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Year + Site)
modellist[[2]]<-formula(Response~Distance + Year + Site)
modellist[[3]]<-formula(Response~dN + Year + Site)
modellist[[4]]<-formula(Response~Distance + dN + Year + Site)
modellist[[5]]<-formula(Response~Distance*dN + Year + Site)
modellist[[6]]<-formula(Response~NULL)
modellist[[7]]<-formula(Response~dC + Year + Site)
modellist[[8]]<-formula(Response~Distance + dC + Year + Site)
modellist[[9]]<-formula(Response~Distance*dC + Year + Site)
modellist[[10]]<-formula(Response~dC + dN + Year + Site)
modellist[[11]]<-formula(Response~Distance + dC + dN + Year + Site)
modellist[[12]]<-formula(Response~Distance*dN*dC + Year + Site)


# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC + Year + Site,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_males",sep="")
wd <- c("R/output/Chapter 4/breeding/males")
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------

par(mar=c(5,4,4,10))
plot(Response~dN,data=dat,col=c(1:7)[dat$Site],pch=c(15:17)[dat$Year])
legend(20,150,levels(dat$Site),col=c(1:7),pch=16,xpd=TRUE)

par(mar=c(5,4,4,10))
plot(Response~dN,data=dat,col=c(1:3)[dat$Year])
legend(20,150,levels(dat$Year),col=c(1:3), pch=16,xpd=TRUE)

m1 <- modeloutput[[3]]

# group years and sites as unique levels of explanatory variable year.site
# create new variable that concatenates year and site

y <- m1$coef[1] + 

year.site <- as.factor(paste(dat$Year,dat$Site,sep=""))

new.dat <- cbind(dat,year.site)

m2 <- lm(Response ~ dN + year.site, data=new.dat)
xyplot(predict(m2, type="response")~dN, groups = year.site, data=new.dat, type="l", col = rainbow(10))



############################################################
############################################################


### FEMALES

# set up blank matrices and lists

### FEMALES

dat <- female

# set up blank matrices and lists

modellist<-list()
modellist[[1]]<-formula(Response~Year + Site)
modellist[[2]]<-formula(Response~Distance + Year + Site)
modellist[[3]]<-formula(Response~dN + Year + Site)
modellist[[4]]<-formula(Response~Distance + dN + Year + Site)
modellist[[5]]<-formula(Response~Distance*dN + Year + Site)
modellist[[6]]<-formula(Response~NULL)
modellist[[7]]<-formula(Response~dC + Year + Site)
modellist[[8]]<-formula(Response~Distance + dC + Year + Site)
modellist[[9]]<-formula(Response~Distance*dC + Year + Site)
modellist[[10]]<-formula(Response~dC + dN + Year + Site)
modellist[[11]]<-formula(Response~Distance + dC + dN + Year + Site)
modellist[[12]]<-formula(Response~Distance*dN*dC + Year + Site)

# creating vector with all possible parameters (from global model)
globalmodel<-lm(Response~Distance*dN*dC + Year + Site,data=dat)

#------------------------------------------
# The AIC analysis
#------------------------------------------
out <- paste(iso.name,"_females",sep="")
wd <- c("R/output/Chapter 4/breeding/females") # set desired directory for AIC output
source("R/scripts/AIC script revised.r")


#-------------------------------------------------
# Graph Initiation vs dN.feathers by Site and Year
#-------------------------------------------------

par(mar=c(5,4,4,10))
plot(Response~isotope,data=dat,col=c(1:7)[dat$Site],pch=c(15:17)[dat$Year])
legend(20,150,levels(dat$Site),col=c(1:7),pch=16,xpd=TRUE)

m <- modeloutput[[3]]

xyplot(predict(m, type="response")~dN, groups=Site, data=dat, type="l")



