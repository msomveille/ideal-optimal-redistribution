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


##  Declare objects ##

nullModel.r.all <- list()
nullModel.dist.r.all <- list()
nullModel.abund.r.all <- list()
ORSIM.r <- vector()


##  Run null models for the six species shown in Fig S1:  ##


spp_sel1 <- c(1,3,4,6,9,13)
spp_sel <- spp_sel_all[spp_sel1]

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r


	#map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	#plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	#plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	#for(i in 1:length(empiricalData_breedingHexagons3)){
	#			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][winter.destinations.nulldist[i],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
	#			plot(spl, add=T, col="black", lwd=1.5, cex=2)
	#}


for(j in 2:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)	
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
	
}









##  Figure S2  ##

spp_sel1 <- c(20,26,28,22,16,15)
spp_sel <- spp_sel_all[spp_sel1]

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	breeding.regions <- spTransform(breeding.regions, sr)
	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/YellowWarbler/YWAR_cleaned.rds")
	wintering.assignments <- wintering.assignments[which(wintering.assignments$Stage=="Wintering"),]
	breeding.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/YellowWarbler/BreedingSurfaces.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/YellowWarbler/WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$Sample == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$Sample == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	wintering_origen <- wintering_origen[-toRemove]
	breeding_origen <- breeding_origen[-toRemove]

	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

	ptsNB <- as.matrix(wintering.assignments[,3:4][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
	
j=2
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/CommonYellowthroat/COYE_WinterProbSurfaces.rds")
	winter.data <- read.csv("Redistribution-model/Genoscapes/CommonYellowthroat/WinterData.csv")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- winter.data[c(5,6)]

	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

	ptsNB <- winter.data[c(5,6)]
	colnames(ptsNB) <- c("x","y")
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	

j=3
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/WilsonsWarbler/WIWA.WinterAssignment.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/WilsonsWarbler/WIWA_WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]

	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}	
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	if(length(toRemove)>0){
		wintering_origen <- wintering_origen[-toRemove]
		breeding_origen <- breeding_origen[-toRemove]
	}
	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))
	ptsNB <- as.matrix(wintering.assignments[,6:5][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	

j=4
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	trackingData <- read.csv("Redistribution-model/Tracking-data/Vermivora/vermNBlatWKbk.csv", header=T)
	trackingData$depLON <- -trackingData$depLON
	trackingData$nbLON <- -trackingData$nbLON
	ptsBR <- as.matrix(trackingData[,3:2][which(trackingData$sp == 2),]) # 1 = gowwar ; 2 = buwwar
	ptsNB <- as.matrix(trackingData[,4:5][which(trackingData$sp == 2),])
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
	
for(j in 5:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
}




	
	
##  Figure S3  ##

spp_sel1 <- c(27,11,19,10,12,14)
spp_sel <- spp_sel_all[spp_sel1]

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/WillowFlycatcher/winter-migrant-assignments.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/WillowFlycatcher/WIFL_WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	if(length(toRemove)>0){
		wintering_origen <- wintering_origen[-toRemove]
		breeding_origen <- breeding_origen[-toRemove]
	}
	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

	ptsNB <- as.matrix(wintering.assignments[,13:12][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r

for(j in 2:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
}






##  Figure S4  ##

spp_sel1 <- c(7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]
	
j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]

	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/CommonLoon/COLO.WinteringAssignment.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/CommonLoon/COLO_WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	if(length(toRemove)>0){
		wintering_origen <- wintering_origen[-toRemove]
		breeding_origen <- breeding_origen[-toRemove]
	}
	ptsBR1 <- as.matrix(wintering.surfaces.origen[,1:2][breeding_origen,])
	sptsBR <- SpatialPoints(ptsBR1)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))
	ptsNB1 <- as.matrix(wintering.assignments[,5:4][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB1)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons_gen <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons_gen <- unlist(winteringHexagons[toKeep])
	ptsBR2 <- vector()
	ptsNB2 <- vector()
	bandingData <- read.csv("Redistribution-model/Banding-data/COLO_combine.csv", header=T)
	for(i in 1:length(bandingData$GISBLONG)){
		if(bandingData$B_SEASON[i] == "B"){
			ptsBR2 = rbind(ptsBR2, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
		}else{
			ptsBR2 = rbind(ptsBR2, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
		}	
	}	
	for(i in 1:length(bandingData$GISBLONG)){
		if(bandingData$R_SEASON[i] == "B"){
			ptsNB2 = rbind(ptsNB2, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
		}else{
		ptsNB2 = rbind(ptsNB2, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
		}	
	}
	dfBR = data.frame(a = 1:length(ptsBR2[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR2, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons_band <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	dfNB = data.frame(a = 1:length(ptsNB2[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB2, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons_band <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])

	empiricalData_breedingHexagons <- c(empiricalData_breedingHexagons_gen, empiricalData_breedingHexagons_band)
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	empiricalData_nonbreedingHexagons <- c(empiricalData_nonbreedingHexagons_gen, empiricalData_nonbreedingHexagons_band)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])
	ptsBR <- rbind(ptsBR1, ptsBR2)
	ptsNB <- rbind(ptsNB1, ptsNB2)

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}

	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r

for(j in 2:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[i,] > 0),]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
}



spp_sel1 <- c(1,3,4,6,9,13,20,26,28,22,16,15,27,11,19,10,12,14,7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]


EMD.r2 <- EMD.r[spp_sel]
nullModel.r.all2 <- nullModel.r.all[spp_sel]
nullModel.abund.r.all2 <- nullModel.abund.r.all[spp_sel]
nullModel.dist.r.all2 <- nullModel.dist.r.all[spp_sel]

pval.nullModel <- vector()
pval.nullModel.abund <- vector()
pval.nullModel.dist <- vector()
for(i in 1:length(EMD.r2)){
	pval.nullModel[i] <- length(which(nullModel.r.all2[[i]] > EMD.r2[i])) / length(nullModel.r.all2[[i]])
	pval.nullModel.abund[i] <- length(which(nullModel.abund.r.all2[[i]] > EMD.r2[i]))/length(nullModel.abund.r.all2[[i]])
	pval.nullModel.dist[i] <- length(which(nullModel.dist.r.all2[[i]] > EMD.r2[i]))/length(nullModel.dist.r.all2[[i]])
}


lapply(nullModel.r.all2, function(x) length(which(x)))




















#####    LONGITUDINAL PATTERNS  ##### 


nullModel.r.all <- list()
nullModel.dist.r.all <- list()
nullModel.abund.r.all <- list()
EMD.r <- vector()


##  Figure S1  ##

spp_sel1 <- c(1,3,4,6,9,13)
spp_sel <- spp_sel_all[spp_sel1]

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r


	#map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	#plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	#plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	#for(i in 1:length(empiricalData_breedingHexagons3)){
	#			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][winter.destinations.nulldist[i],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
	#			plot(spl, add=T, col="black", lwd=1.5, cex=2)
	#}


for(j in 2:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
	
}









##  Figure S2  ##

spp_sel1 <- c(20,26,28,22,16,15)
spp_sel <- spp_sel_all[spp_sel1]

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	breeding.regions <- spTransform(breeding.regions, sr)
	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/YellowWarbler/YWAR_cleaned.rds")
	wintering.assignments <- wintering.assignments[which(wintering.assignments$Stage=="Wintering"),]
	breeding.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/YellowWarbler/BreedingSurfaces.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/YellowWarbler/WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$Sample == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$Sample == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	wintering_origen <- wintering_origen[-toRemove]
	breeding_origen <- breeding_origen[-toRemove]

	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

	ptsNB <- as.matrix(wintering.assignments[,3:4][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
	
j=2
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/CommonYellowthroat/COYE_WinterProbSurfaces.rds")
	winter.data <- read.csv("Redistribution-model/Genoscapes/CommonYellowthroat/WinterData.csv")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- winter.data[c(5,6)]

	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

	ptsNB <- winter.data[c(5,6)]
	colnames(ptsNB) <- c("x","y")
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	

j=3
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/WilsonsWarbler/WIWA.WinterAssignment.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/WilsonsWarbler/WIWA_WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]

	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}	
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	if(length(toRemove)>0){
		wintering_origen <- wintering_origen[-toRemove]
		breeding_origen <- breeding_origen[-toRemove]
	}
	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))
	ptsNB <- as.matrix(wintering.assignments[,6:5][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	

j=4
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	trackingData <- read.csv("Redistribution-model/Tracking-data/Vermivora/vermNBlatWKbk.csv", header=T)
	trackingData$depLON <- -trackingData$depLON
	trackingData$nbLON <- -trackingData$nbLON
	ptsBR <- as.matrix(trackingData[,3:2][which(trackingData$sp == 2),]) # 1 = gowwar ; 2 = buwwar
	ptsNB <- as.matrix(trackingData[,4:5][which(trackingData$sp == 2),])
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
	
for(j in 5:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
}




	
	
##  Figure S3  ##

spp_sel1 <- c(27,11,19,10,12,14)
spp_sel <- spp_sel_all[spp_sel1]

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/WillowFlycatcher/winter-migrant-assignments.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/WillowFlycatcher/WIFL_WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	if(length(toRemove)>0){
		wintering_origen <- wintering_origen[-toRemove]
		breeding_origen <- breeding_origen[-toRemove]
	}
	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

	ptsNB <- as.matrix(wintering.assignments[,13:12][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r

for(j in 2:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
}






##  Figure S4  ##

spp_sel1 <- c(7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]
	
j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]

	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/CommonLoon/COLO.WinteringAssignment.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/CommonLoon/COLO_WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	if(length(toRemove)>0){
		wintering_origen <- wintering_origen[-toRemove]
		breeding_origen <- breeding_origen[-toRemove]
	}
	ptsBR1 <- as.matrix(wintering.surfaces.origen[,1:2][breeding_origen,])
	sptsBR <- SpatialPoints(ptsBR1)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))
	ptsNB1 <- as.matrix(wintering.assignments[,5:4][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB1)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons_gen <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons_gen <- unlist(winteringHexagons[toKeep])
	ptsBR2 <- vector()
	ptsNB2 <- vector()
	bandingData <- read.csv("Redistribution-model/Banding-data/COLO_combine.csv", header=T)
	for(i in 1:length(bandingData$GISBLONG)){
		if(bandingData$B_SEASON[i] == "B"){
			ptsBR2 = rbind(ptsBR2, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
		}else{
			ptsBR2 = rbind(ptsBR2, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
		}	
	}	
	for(i in 1:length(bandingData$GISBLONG)){
		if(bandingData$R_SEASON[i] == "B"){
			ptsNB2 = rbind(ptsNB2, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
		}else{
		ptsNB2 = rbind(ptsNB2, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
		}	
	}
	dfBR = data.frame(a = 1:length(ptsBR2[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR2, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons_band <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	dfNB = data.frame(a = 1:length(ptsNB2[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB2, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons_band <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])

	empiricalData_breedingHexagons <- c(empiricalData_breedingHexagons_gen, empiricalData_breedingHexagons_band)
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	empiricalData_nonbreedingHexagons <- c(empiricalData_nonbreedingHexagons_gen, empiricalData_nonbreedingHexagons_band)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])
	ptsBR <- rbind(ptsBR1, ptsBR2)
	ptsNB <- rbind(ptsNB1, ptsNB2)

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r

for(j in 2:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),1][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}
	
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
}


nullModel.r.all.longitude <- nullModel.r.all
nullModel.dist.r.all.longitude <- nullModel.dist.r.all
nullModel.abund.r.all.longitude <- nullModel.abund.r.all
EMD.r.longitude <- EMD.r


spp_sel1 <- c(1,3,4,6,9,13,20,26,28,22,16,15,27,11,19,10,12,14,7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]


EMD.r.longitude2 <- EMD.r.longitude[spp_sel]
nullModel.r.all.longitude2 <- nullModel.r.all.longitude[spp_sel]
nullModel.abund.r.all.longitude2 <- nullModel.abund.r.all.longitude[spp_sel]
nullModel.dist.r.all.longitude2 <- nullModel.dist.r.all.longitude[spp_sel]

pval.nullModel <- vector()
pval.nullModel.abund <- vector()
pval.nullModel.dist <- vector()
for(i in 1:length(EMD.r.longitude2)){
	pval.nullModel[i] <- length(which(nullModel.r.all.longitude2[[i]] > EMD.r.longitude2[i])) / length(nullModel.r.all.longitude2[[i]])
	pval.nullModel.abund[i] <- length(which(nullModel.abund.r.all.longitude2[[i]] > EMD.r.longitude2[i]))/length(nullModel.abund.r.all.longitude2[[i]])
	pval.nullModel.dist[i] <- length(which(nullModel.dist.r.all.longitude2[[i]] > EMD.r.longitude2[i]))/length(nullModel.dist.r.all.longitude2[[i]])
}


lapply(nullModel.r.all2, function(x) length(which(x)))

spp_name[spp_sel]












#####    LATITUDINAL PATTERNS  ##### 


nullModel.r.all <- list()
nullModel.dist.r.all <- list()
nullModel.abund.r.all <- list()
EMD.r <- vector()


##  Figure S1  ##

spp_sel1 <- c(1,3,4,6,9,13)
spp_sel <- spp_sel_all[spp_sel1]

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r


	#map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	#plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	#plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	#for(i in 1:length(empiricalData_breedingHexagons3)){
	#			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][winter.destinations.nulldist[i],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
	#			plot(spl, add=T, col="black", lwd=1.5, cex=2)
	#}


for(j in 2:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
	
}









##  Figure S2  ##

spp_sel1 <- c(20,26,28,22,16,15)
spp_sel <- spp_sel_all[spp_sel1]

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	breeding.regions <- spTransform(breeding.regions, sr)
	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/YellowWarbler/YWAR_cleaned.rds")
	wintering.assignments <- wintering.assignments[which(wintering.assignments$Stage=="Wintering"),]
	breeding.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/YellowWarbler/BreedingSurfaces.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/YellowWarbler/WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$Sample == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$Sample == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	wintering_origen <- wintering_origen[-toRemove]
	breeding_origen <- breeding_origen[-toRemove]

	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

	ptsNB <- as.matrix(wintering.assignments[,3:4][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
	
j=2
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/CommonYellowthroat/COYE_WinterProbSurfaces.rds")
	winter.data <- read.csv("Redistribution-model/Genoscapes/CommonYellowthroat/WinterData.csv")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- winter.data[c(5,6)]

	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

	ptsNB <- winter.data[c(5,6)]
	colnames(ptsNB) <- c("x","y")
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	

j=3
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/WilsonsWarbler/WIWA.WinterAssignment.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/WilsonsWarbler/WIWA_WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]

	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}	
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	if(length(toRemove)>0){
		wintering_origen <- wintering_origen[-toRemove]
		breeding_origen <- breeding_origen[-toRemove]
	}
	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))
	ptsNB <- as.matrix(wintering.assignments[,6:5][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r

j=4
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	trackingData <- read.csv("Redistribution-model/Tracking-data/Vermivora/vermNBlatWKbk.csv", header=T)
	trackingData$depLON <- -trackingData$depLON
	trackingData$nbLON <- -trackingData$nbLON
	ptsBR <- as.matrix(trackingData[,3:2][which(trackingData$sp == 2),]) # 1 = gowwar ; 2 = buwwar
	ptsNB <- as.matrix(trackingData[,4:5][which(trackingData$sp == 2),])
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
	
for(j in 5:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
}




	
	
##  Figure S3  ##

spp_sel1 <- c(27,11,19,10,12,14)
spp_sel <- spp_sel_all[spp_sel1]

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/WillowFlycatcher/winter-migrant-assignments.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/WillowFlycatcher/WIFL_WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	if(length(toRemove)>0){
		wintering_origen <- wintering_origen[-toRemove]
		breeding_origen <- breeding_origen[-toRemove]
	}
	ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
	sptsBR <- SpatialPoints(ptsBR)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

	ptsNB <- as.matrix(wintering.assignments[,13:12][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons <- unlist(winteringHexagons[toKeep])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons)
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r

for(j in 2:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
}






##  Figure S4  ##

spp_sel1 <- c(7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]
	
j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]

	wintering.assignments <- readRDS("Redistribution-model/Genoscapes/CommonLoon/COLO.WinteringAssignment.rds")
	wintering.surfaces.origen <- readRDS("Redistribution-model/Genoscapes/CommonLoon/COLO_WinterProbSurfaces.rds")

	# Get hexagons for genoscapes breeding and wintering destinations 
	breeding_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
	}
	breeding_origen <- breeding_origen[-c(1,2)]
	wintering_origen <- vector()
	for(k in 3:length(colnames(wintering.surfaces.origen))){
		if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
			wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
		}else{
			wintering_origen[k] <- NA
		}
	}
	wintering_origen <- wintering_origen[-c(1,2)]
	toRemove <- which(is.na(wintering_origen) == T)
	if(length(toRemove)>0){
		wintering_origen <- wintering_origen[-toRemove]
		breeding_origen <- breeding_origen[-toRemove]
	}
	ptsBR1 <- as.matrix(wintering.surfaces.origen[,1:2][breeding_origen,])
	sptsBR <- SpatialPoints(ptsBR1)
	proj4string(sptsBR) <- proj4string(hexgrid3_stem)
	breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
	breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))
	ptsNB1 <- as.matrix(wintering.assignments[,5:4][wintering_origen,])
	sptsNB <- SpatialPoints(ptsNB1)
	proj4string(sptsNB) <- proj4string(hexgrid3_stem)
	winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
	winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

	toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
	empiricalData_breedingHexagons_gen <- unlist(breedingHexagons[toKeep])
	empiricalData_nonbreedingHexagons_gen <- unlist(winteringHexagons[toKeep])
	ptsBR2 <- vector()
	ptsNB2 <- vector()
	bandingData <- read.csv("Redistribution-model/Banding-data/COLO_combine.csv", header=T)
	for(i in 1:length(bandingData$GISBLONG)){
		if(bandingData$B_SEASON[i] == "B"){
			ptsBR2 = rbind(ptsBR2, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
		}else{
			ptsBR2 = rbind(ptsBR2, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
		}	
	}	
	for(i in 1:length(bandingData$GISBLONG)){
		if(bandingData$R_SEASON[i] == "B"){
			ptsNB2 = rbind(ptsNB2, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
		}else{
		ptsNB2 = rbind(ptsNB2, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
		}	
	}
	dfBR = data.frame(a = 1:length(ptsBR2[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR2, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons_band <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	dfNB = data.frame(a = 1:length(ptsNB2[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB2, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons_band <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])

	empiricalData_breedingHexagons <- c(empiricalData_breedingHexagons_gen, empiricalData_breedingHexagons_band)
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	empiricalData_nonbreedingHexagons <- c(empiricalData_nonbreedingHexagons_gen, empiricalData_nonbreedingHexagons_band)
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])
	ptsBR <- rbind(ptsBR1, ptsBR2)
	ptsNB <- rbind(ptsNB1, ptsNB2)

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r

for(j in 2:length(spp_sel)){

	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	ptsBR <- vector()
	ptsNB <- vector()
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Redistribution-model/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$B_SEASON[i] == "B"){
				ptsBR = rbind(ptsBR, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsBR = rbind(ptsBR, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
		for(i in 1:length(bandingData$GISBLONG)){
			if(bandingData$R_SEASON[i] == "B"){
				ptsNB = rbind(ptsNB, c(bandingData$GISBLONG[i], bandingData$GISBLAT[i]))
			}else{
				ptsNB = rbind(ptsNB, c(bandingData$GISRLONG[i], bandingData$GISRLAT[i]))
			}	
		}
	}
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Redistribution-model/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}
	dfBR = data.frame(a = 1:length(ptsBR[,1]))
	sp.ptsBR <- SpatialPointsDataFrame(ptsBR, dfBR)
	proj4string(sp.ptsBR) <- sr
	empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
	empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
	dfNB = data.frame(a = 1:length(ptsNB[,1]))
	sp.ptsNB <- SpatialPointsDataFrame(ptsNB, dfNB)
	proj4string(sp.ptsNB) <- sr
	empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])
	empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])

	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]
	distanceMat.empirical <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}
	winter.destinations <- list()
	EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),2][which(EMD_flows_sel[i,] > 0)]
		if(length(which(EMD_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- mean(winter.destinations[[i]])
		}
	}
	distanceMat.simulatedBR <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulatedBR <- rbind(distanceMat.simulatedBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulatedBR <- distanceMat.simulatedBR[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}
	
	EMD.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
		
	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	for(k in 1:100){
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
		distanceMat.nullBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullBR <- rbind(distanceMat.nullBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.nullBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
		distanceMat.nullabundBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabundBR <- rbind(distanceMat.nullabundBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabundBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		
		
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x))))[1])
		distanceMat.nulldistBR <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldistBR <- rbind(distanceMat.nulldistBR, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldistBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r
}


nullModel.r.all.latitude <- nullModel.r.all
nullModel.dist.r.all.latitude <- nullModel.dist.r.all
nullModel.abund.r.all.latitude <- nullModel.abund.r.all
EMD.r.latitude <- EMD.r


spp_sel1 <- c(1,3,4,6,9,13,20,26,28,22,16,15,27,11,19,10,12,14,7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]


EMD.r.latitude2 <- EMD.r.latitude[spp_sel]
nullModel.r.all.latitude2 <- nullModel.r.all.latitude[spp_sel]
nullModel.abund.r.all.latitude2 <- nullModel.abund.r.all.latitude[spp_sel]
nullModel.dist.r.all.latitude2 <- nullModel.dist.r.all.latitude[spp_sel]

pval.nullModel <- vector()
pval.nullModel.abund <- vector()
pval.nullModel.dist <- vector()
for(i in 1:length(EMD.r.latitude2)){
	pval.nullModel[i] <- length(which(nullModel.r.all.latitude2[[i]] > EMD.r.latitude2[i])) / length(nullModel.r.all.latitude2[[i]])
	pval.nullModel.abund[i] <- length(which(nullModel.abund.r.all.latitude2[[i]] > EMD.r.latitude2[i]))/length(nullModel.abund.r.all.latitude2[[i]])
	pval.nullModel.dist[i] <- length(which(nullModel.dist.r.all.latitude2[[i]] > EMD.r.latitude2[i]))/length(nullModel.dist.r.all.latitude2[[i]])
}


lapply(nullModel.r.all2, function(x) length(which(x)))

cbind(spp_name[spp_sel], EMD.r.latitude2, pval.nullModel, pval.nullModel.abund, pval.nullModel.dist)

