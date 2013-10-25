##########################################################################
# Merge breeding banding, nest, and isotope data
# Name - Samantha Franks
# 16 July 2013
# 14 Aug 2013
# 23 Oct 2013
#
#	Modified to reflect new assignments with sex-specific priors and corrected isootpe values Oct 2013
##########################################################################

### Analyses
# - initation date vs distance, dN, dC (2008-10 only)
# - clutch mass vs dN. dC (females 2012)
# - incubation mass change vs dN, dC (2012, control for nest)

### Datasets required
# - banding, iso, nest data 2008-10
# - banding, iso, nest, COE 2009-2012 Nome only
# - banding, nest, COE 2012

rm(list=ls()) #clear all previous
library(NCStats) # for Subset function

nest<-read.csv("R/sandpiper/nest data all sites 2008-2010 and 2012_clean.csv",header=T)
iso<-read.csv("R/sandpiper/C N and H isotope data breeding_revised.csv", header=T)
band<-read.csv("R/sandpiper/ALL banding breeding 2008-2010 and 2012.csv", header=T)
assign<-read.csv("R/sandpiper/summary all breeding assignments_with sex specific priors_revised_20131021.csv",header=TRUE)
COE <- read.csv("R/sandpiper/COE Nome 2012.csv",header=T)

nest$year <- as.factor(nest$year)
band$Year <- as.factor(band$Year)
iso$Year <- as.factor(iso$Year)
assign$Year <- as.factor(assign$Year)

#---------------------------
# Give band numbers unique ids for each year
#---------------------------

band.id<-paste(band$Year,"_",band$Band,sep="")
band<-cbind(band.id,band)

iso.id<-paste(iso$Year,"_",iso$Band,sep="")
iso<-cbind(iso.id,iso)

assign.id<-paste(assign$Year,"_",assign$Band,sep="")
assign<-cbind(assign.id,assign)

COE.id <- paste("2012_", COE$Band, sep="")
COE <- data.frame(Band = COE.id, COE)

#---------------------------
# Create datasets for relevant analyses
#---------------------------

### Dataset 1 - banding, iso, nest data 2008-10

# remove 2012 individuals from this dataset
band1 <- Subset(band, Year != "2012")
nest1 <- Subset(nest, year != "2012")

# merge banding 2008-10 with isotope data
merge1<-merge(band1,iso,by.x="band.id",by.y="iso.id",all.x=T,all.y=T)
band.iso <- with(merge1, data.frame(Site=Site.x, Year=Year.x, Latitude, Longitude, Band=band.id, fullnest, Species=Species.x, Recap.btwn.years, Revised.Age, Age=Age.x, Sex=Sex.x, Wing.mm, Culmen.mm, Tarsus.mm=Corrected.long.tarsus.mm, Weight.g, dN=dN.feathers, dC=dC.feathers, dD=dD.feathers))

# remove all individuals without any isotope data from dataset (will not be part of analysis)
BI0810 <- Subset(band.iso, dC != "NA")

# merge banding/iso data with assignment data	
merge2<-merge(BI0810,assign,by.x="Band",by.y="assign.id",all.x=T,all.y=T)
assign.band.iso <- with(merge2, data.frame(Site=Site.x, Year=Year.x, Latitude=Latitude.x, Longitude=Longitude.x, Band, fullnest, Species, Recap.btwn.years, Revised.Age,Age=Age.x, Sex=Sex.x, Wing.mm, Culmen.mm, Tarsus.mm, Weight.g, dN, dC, dD, Origin=Region.assigned))

ABI0810 <- Subset(assign.band.iso, Revised.Age != "SY")	# remove SYs from dataset

# merge banding/iso/assignment data with nesting data
merge3<-merge(nest1,ABI0810,by.x="fullnest",by.y="fullnest",all.x=T,all.y=T)

nest.assign.band.iso <- with(merge3, data.frame(Site1=site, Year1=year, Site2=Site, Year2=Year, Latitude, Longitude, Band, fullnest, INIT.JULIAN, failed, fail.cause, Species, Recap.btwn.years, Revised.Age, Age, Sex, Wing.mm, Culmen.mm, Tarsus.mm, Weight.g, dN, dC, dD, Origin))

fulldata <- Subset(nest.assign.band.iso, INIT.JULIAN != "NA")
fulldata2 <- Subset(fulldata, dN != "NA")

write.csv(nest.assign.band.iso,"R/sandpiper/merged breeding nest banding isotope assignment data 2008-10_20131023.csv",row.names=FALSE)


#####################################
#####################################
#####################################
#####################################
#####################################

### Dataset 2 - banding, iso, nest, COE data - Nome 2009-2012

### STEP 1 - create 2009-10 dataset
nest2 <- Subset(nest, site == "Nome")
nestNome0910 <- Subset(nest2, year == "2009" | year == "2010")
band2 <- Subset(band, Site == "Nome" & Species == "WESA")
bandNome0910 <- Subset(band2, Year == "2009" | Year == "2010")
iso2 <- Subset(iso, Site == "Nome")
isoNome0910 <- Subset(iso2, Year == "2009" | Year == "2010")

# merge banding with isotope data
merge1<-merge(bandNome0910, isoNome0910, by.x="band.id", by.y="iso.id", all.x=T, all.y=T)

band.iso.Nome <- with(merge1, data.frame(Site=Site.x, Year=Year.x, Latitude, Longitude, Band=band.id, fullnest=fullnest.x, Species=Species.x, Recap.btwn.years, Revised.Age, Age=Age.x, Sex=Sex.x, Wing.mm, Culmen.mm, Tarsus.mm=Corrected.long.tarsus.mm, Weight.g, dN.feathers, dC.feathers, dD.feathers))

# remove any observation that doesn't have an isotope measurement or is an SY by isotope age
BINome <- Subset(band.iso.Nome, dN.feathers != "NA" & Revised.Age != "SY")

# merge banding/iso data with nesting data
merge2<-merge(nestNome0910, BINome, by.x="fullnest", by.y="fullnest", all.x=T, all.y=T)

nest.band.iso.Nome <- with(merge2, data.frame(Year1=year, Year2=Year, Band, fullnest, INIT.JULIAN, failed, fail.cause, Recap.btwn.years, Revised.Age, Age, Sex, Wing.mm, Culmen.mm, Tarsus.mm, Weight.g, dN.feathers, dC.feathers))

# remove any observation that doesn't have both an initiation date AND an isotope measurement
NBINome <- Subset(nest.band.iso.Nome, dN.feathers != "NA" & INIT.JULIAN !="NA")

### Step 2 - create 2012 dataset

COENome <- with(COE, data.frame(COE.id, Nest.ID, Fullnest.ID, Band, dN, dC, Sex))
nestNome12 <- Subset(nest2, year == "2012")
bandNome12 <- Subset(band2, Year == "2012")

# merge banding with isotope data
merge1 <- merge(bandNome12, COENome, by.x="band.id", by.y="Band", all.x=T, all.y=T)

bandCOENome <- with(merge1, data.frame(Site, Year, fullnest, Fullnest.ID, Band=band.id, Age, Sex=Sex.y, Wing.mm, Culmen.mm, Tarsus.mm=Diagonal.tarsus.mm, Weight.g, dN, dC))

BCNome <- bandCOENome[-40,]

# merge banding/iso data with nest data
merge2 <- merge(nestNome12, BCNome, by.x = "fullnest", by.y = "fullnest", all.x = T, all.y = T)

NBCNome <- with(merge2, data.frame(Year, Site, fullnest, Band, INIT.JULIAN, failed, fail.cause, Age, Sex, Wing.mm, Culmen.mm, Tarsus.mm, Weight.g, dN, dC))

# merge2[,c("fullnest","Band","INIT.JULIAN","dN")]

### Step 3 - merge 2009-10 AND 2012 datasets

# data variables are: year, site, fullnest, band, init.julian, failed, fail.cause, age/revised age, sex, wing, culmen, tarsus, weight, dN, dC


#######
NBINome2 <- Subset(NBINome, Revised.Age != "SY")
NBINome3 <- Subset(NBINome2, INIT.JULIAN != "NA")
NBINome4 <- with(NBINome3, data.frame(Year = Year1, Band, fullnest, INIT.JULIAN, failed, fail.cause, Recap = Recap.btwn.years, Revised.Age, Age, Sex, Wing.mm, Culmen.mm, Tarsus.mm, Weight.g, dN = dN.feathers, dC = dC.feathers))

# final merge - COE 2012 data relevant columns

merge3 <- merge(COENome, NBINome, by.x = "Fullnest.ID", by.y = "fullnest", all.x = T, all.y = T)

# final modifications and extraction of relevant rows with full data in Excel
write.csv(merge3, "R/output/Chapter 4/COE MS/Nome 2009-2012 data.csv", row.names = F)


### 2012 dataset

COENome <- with(COE, data.frame(COE.id, Nest.ID, Fullnest.ID, Band, dN, dC, Sex))

### merge 2009-10 and 2012 datasets