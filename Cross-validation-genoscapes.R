library(geosphere)
library(fields)
library(MASS)
library(raster)
library(move)
library(vioplot)
library(igraph)
library(sp)
library(rgeos)
library(dggridR)
library(rworldmap)
library(mapplots)
library(ebirdst)
library(car)
library(ade4)

setwd("/Volumes/Marius_SSD/American-Flyway/Connectivity_NAbirds/Redistribution-model")



##  Common Yellowthroat  ##

# Load probability surfaces obtained using OriGen
breeding.surfaces.origen <- readRDS("Data/Genoscapes/CommonYellowthroat/COYE_BreedingSurfaces.rds")

# Load empirical location of individuals sampled on the breeding ground
empirical.breeding.locations <- read.table("Data/Genoscapes/CommonYellowthroat/COYE_breeding.loc", header=T)

# Compute distance between empirical and predicted breeding locations
distance.error.COYE <- vector()
for(k in 3:length(colnames(breeding.surfaces.origen))){
	if(length(which(as.character(empirical.breeding.locations$Sample) == colnames(breeding.surfaces.origen)[k])) > 0){
		ee <- empirical.breeding.locations[which(as.character(empirical.breeding.locations$Sample) == colnames(breeding.surfaces.origen)[k]),][,c(5,6)] 
		ss <- breeding.surfaces.origen[which(breeding.surfaces.origen[,k] == max(breeding.surfaces.origen[,k])),1:2]
		distance.error.COYE[k] <- rdist.earth(ee,ss, miles=F)
	}	
}
distance.error.COYE <- distance.error.COYE[-c(1,2)]




##  Willow Flycatcher  ##

# Load probability surfaces obtained using OriGen
breeding.surfaces.origen <- readRDS("Data/Genoscapes/WillowFlycatcher/WIFL_BreedingSurfaces.rds")

# Load empirical location of individuals sampled on the breeding ground
empirical.breeding.locations <- read.table("Data/Genoscapes/WillowFlycatcher/WIFL_breeding.loc", header=T)

# Compute distance between empirical and predicted breeding locations
distance.error.WIFL <- vector()
for(k in 3:length(colnames(breeding.surfaces.origen))){
	if(length(which(as.character(empirical.breeding.locations$ID) == colnames(breeding.surfaces.origen)[k])) > 0){
		ee <- empirical.breeding.locations[which(as.character(empirical.breeding.locations$ID) == colnames(breeding.surfaces.origen)[k]),][,c(5,6)] 
		ss <- breeding.surfaces.origen[which(breeding.surfaces.origen[,k] == max(breeding.surfaces.origen[,k])),1:2]
		distance.error.WIFL[k] <- rdist.earth(ee,ss, miles=F)
	}	
}
distance.error.WIFL <- distance.error.WIFL[-c(1,2)]




##  Yellow Warbler  ##

# Load probability surfaces obtained using OriGen
breeding.surfaces.origen <- readRDS("Data/Genoscapes/YellowWarbler/BreedingSurfaces.rds")

# Load empirical location of individuals sampled on the breeding ground
wintering.assignments <- readRDS("Data/Genoscapes/YellowWarbler/YWAR_cleaned.rds")
empirical.breeding.locations <- wintering.assignments[which(wintering.assignments$Stage=="Breeding"),]

# Compute distance between empirical and predicted breeding locations
distance.error.YEWA <- vector()
for(k in 3:length(colnames(breeding.surfaces.origen))){
	if(length(which(empirical.breeding.locations$Sample == colnames(breeding.surfaces.origen)[k])) > 0){
		ee <- empirical.breeding.locations[which(empirical.breeding.locations$Sample == colnames(breeding.surfaces.origen)[k]),][,c(3,4)] 
		ss <- breeding.surfaces.origen[which(breeding.surfaces.origen[,k] == max(breeding.surfaces.origen[,k])),1:2]
		distance.error.YEWA[k] <- rdist.earth(ee,ss, miles=F)
	}	
}
distance.error.YEWA <- distance.error.YEWA[-c(1,2)]




##  Wilson's Warbler  ##

# Load probability surfaces obtained using OriGen
breeding.surfaces.origen <- readRDS("Data/Genoscapes/WilsonsWarbler/WIWA_BreedingSurfaces.rds")

# Load empirical location of individuals sampled on the breeding ground
empirical.breeding.locations <- read.table("Data/Genoscapes/WilsonsWarbler/WIWA_breeding.loc", header=T)

# Compute distance between empirical and predicted breeding locations
distance.error.WIWA <- vector()
for(k in 3:length(colnames(breeding.surfaces.origen))){
	if(length(which(as.character(empirical.breeding.locations$ID) == colnames(breeding.surfaces.origen)[k])) > 0){
		ee <- empirical.breeding.locations[which(as.character(empirical.breeding.locations$ID) == colnames(breeding.surfaces.origen)[k]),][,c(5,6)] 
		ss <- breeding.surfaces.origen[which(breeding.surfaces.origen[,k] == max(breeding.surfaces.origen[,k])),1:2]
		distance.error.WIWA[k] <- rdist.earth(ee,ss, miles=F)
	}	
}
distance.error.WIWA <- distance.error.WIWA[-c(1,2)]




##  Common Loon  ##

# Load probability surfaces obtained using OriGen
breeding.surfaces.origen <- readRDS("Data/Genoscapes/CommonLoon/COLO_BreedingSurfaces.rds")

# Load empirical location of individuals sampled on the breeding ground
empirical.breeding.locations <- read.table("Data/Genoscapes/CommonLoon/COLO_breeding.loc", header=T)

# Compute distance between empirical and predicted breeding locations
distance.error.COLO <- vector()
for(k in 3:length(colnames(breeding.surfaces.origen))){
	if(length(which(as.character(empirical.breeding.locations$ID) == colnames(breeding.surfaces.origen)[k])) > 0){
		ee <- empirical.breeding.locations[which(as.character(empirical.breeding.locations$ID) == colnames(breeding.surfaces.origen)[k]),][,c(5,6)] 
		ss <- breeding.surfaces.origen[which(breeding.surfaces.origen[,k] == max(breeding.surfaces.origen[,k])),1:2]
		distance.error.COLO[k] <- rdist.earth(ee,ss, miles=F)
	}	
}
distance.error.COLO <- distance.error.COLO[-c(1,2)]



##  FIGURE S12: results of cross-validation  ##

par(mfrow=c(2,3), mar=c(3,3,3,1), mgp=c(2,1,0), bg="white")

hist(distance.error.COYE, xlab="Error distance (km)", col="light grey", breaks=20, border="grey", main="Common Yellowthroat")
abline(v=mean(distance.error.COYE, na.rm=T))
hist(distance.error.WIFL, xlab="Error distance (km)", col="light grey", breaks=20, border="grey", main="Willow Flycatcher")
abline(v=mean(distance.error.WIFL, na.rm=T))
hist(distance.error.YEWA, xlab="Error distance (km)", col="light grey", breaks=20, border="grey", main="Yellow Warbler")
abline(v=mean(distance.error.YEWA, na.rm=T))
hist(distance.error.WIWA, xlab="Error distance (km)", col="light grey", breaks=20, border="grey", main="Wilson's Warbler")
abline(v=mean(distance.error.WIWA, na.rm=T))
hist(distance.error.COLO, xlab="Error distance (km)", col="light grey", breaks=20, border="grey", main="Common Loon")
abline(v=mean(distance.error.COLO, na.rm=T))



