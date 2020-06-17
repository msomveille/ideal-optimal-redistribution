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


##  Construct an hexagon grid covering the Americas  ##

hexgrid <- dgconstruct(projection="ISEA", topology="HEXAGON", res=7, metric=T)
hexgrid_center <- dgSEQNUM_to_GEO(hexgrid, 1:21872) # 65612 / 590492
hexgrid_centroids <- cbind(hexgrid_center$lon_deg, hexgrid_center$lat_deg)
hex_sel <- which(hexgrid_centroids[,1] < -30 & hexgrid_centroids[,1] > -170 & hexgrid_centroids[,2] > -60 & hexgrid_centroids[,2] < 80)
hexgrid2_stem <- dgcellstogrid(hexgrid, hex_sel, frame=F,wrapcells=TRUE)
hexgrid2_stem_centroids <- matrix(unlist(lapply(hexgrid2_stem@polygons, function(x) x@labpt)), byrow=T, ncol=2)
newmap <- getMap(resolution = "low")
newmap <- spTransform(newmap, proj4string(hexgrid2_stem))
newmap@data$world <- rep(1,length(newmap@data$SOVEREIGNT))
newmap <- gBuffer(newmap, byid=TRUE, width=0)
newmap <- gUnaryUnion(newmap, id=newmap@data$world)
hexgrid3_stem <- gIntersection(hexgrid2_stem, newmap, byid=T)
hexgrid3_stem_centroids <- matrix(unlist(lapply(hexgrid3_stem@polygons, function(x) x@labpt)), byrow=T, ncol=2)
sr <- proj4string(hexgrid3_stem)


##  Lists of species names for the various datasets  ##

spp_name <- c("Wood Thrush", "Yellow Warbler", "Tree Swallow", "Swainson's Thrush", "Wilson's Warbler", "Common Yellowthroat", "Canada Warbler", "American Redstart", "Hermit Thrush", "Kentucky Warbler", "American Kestrel", "American Robin", "Anna's Hummingbird", "Common Loon", "Burrowing Owl", "Painted Bunting", "Western Tanager", "Willow Flycatcher", "American Goldfinch", "Brown headed Cowbird", "Brown Thrasher", "Common Grackle", "Purple Finch", "Red winged Blackbird", "White throated Sparrow", "Blue winged Warbler", "Golden winged Warbler", "Ovenbird", "Prothonotary Warbler", "Osprey", "Barn Swallow", "Golden crowned Sparrow", "Gray Catbird", "Blackpoll Warbler")
spp <- c("woothr", "yelwar", "treswa", "swathr", "wlswar", "comyel", "canwar", "amered", "herthr", "kenwar", "amekes", "amerob", "annhum", "comloo", "burowl", "paibun", "westan", "wilfly", "amegfi", "bnhcow", "brnthr", "comgra", "purfin", "rewbla", "whtspa", "buwwar", "gowwar", "ovenbi1", "prowar", "osprey", "barswa", "gocspa", "grycat", "bkpwar")
spp_banding <- c("WOTH", "YEWA", "TRES", "SWTH", NA, "COYE", NA, "AMRE", "HETH", NA, "AMKE", "AMRO", NA, "COLO", "BUOW", "PABU", NA, "WIFL", "AMGO", "BHCO", "BRTH", "COGR", "PUFI", "RWBL", "WTSP", NA, NA, "OVEN", NA, "OSPR", "BARS", "GCSP", "GRCA", NA)
spp_tracking <- c("wood thrush", NA, NA, "swainson's thrush", NA, NA, NA, NA, "hermit thrush", NA, NA, NA, NA, NA, NA, NA, NA ,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "ovenbird", "prothonotary warbler", "osprey", "barn swallow", NA, "gray catbird", "blackpoll warbler")


# Selected species for the analysis
spp_sel_all <- c(1,3,4,9,11,12,14,15,19,20,21,22,23,24,25,28,30,31,33,2,16,26,27,32,34,6,18,5)



##  Declare objects ##

nullModel.r.all <- list()
nullModel.dist.r.all <- list()
nullModel.abund.r.all <- list()
ORSIM.r <- vector()
nullModel.r.all.longitude <- list()
nullModel.dist.r.all.longitude <- list()
nullModel.abund.r.all.longitude <- list()
ORSIM.r.longitude <- vector()
nullModel.r.all.latitude <- list()
nullModel.dist.r.all.latitude <- list()
nullModel.abund.r.all.latitude <- list()
ORSIM.r.latitude <- vector()

##  Run null models for the six species shown in Fig S1: Wood Thrush, Swainson's Thrush, Hermit Thrush, American Robin, American Goldfinch, and Purple Finch  ##

# Selected species
spp_sel1 <- c(1,3,4,6,9,13)
spp_sel <- spp_sel_all[spp_sel1]

for(j in 1:length(spp_sel)){

	##  Load species seasonal abundance distributions  ##
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

	##  Load ORSIM output (simulated migratory connectivity)  ##
	ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

	##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species  ##
	ptsBR <- vector()
	ptsNB <- vector()
	# Banding data
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Data/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
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
	# Tracking data
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Data/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}

	##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
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

	##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

	##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical <- vector() 
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	
	##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	
	##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}

	##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
	winter.destinations <- list()
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
		if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}

	##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	
	##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	
	##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulated <- distanceMat.simulated[-aaa,]
		distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
		distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}

	##  Compute the Mantel correlation coefficients	##
	ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}

	##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
	

	##  Compute the null models  ##	

	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	nullModel.r.longitude <- vector()
	nullModel.dist.r.longitude <- vector()
	nullModel.abund.r.longitude <- vector()
	nullModel.r.latitude <- vector()
	nullModel.dist.r.latitude <- vector()
	nullModel.abund.r.latitude <- vector()
	
	for(k in 1:1000){
	
		##  Null model 1: randomizing migratory destinations  ##
	
		# Sampling winter destinations at random across the wintering ground
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
		##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
		# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
		##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
		# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
	nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
	nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
	nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
	nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
	nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
	nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude
}


##  Plot simulated conectivity from the null model
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][winter.destinations.nulldist[i],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
	plot(spl, add=T, col="black", lwd=1.5, cex=2)
}





##  Run null models for the six species shown in Fig S1: Yellow Warbler, Common Yellowthroat, Wilson's Warbler, Blue winged Warbler, Ovenbird, and White-throated Sparrow  ##

# Selected species
spp_sel1 <- c(20,26,28,22,16,15)
spp_sel <- spp_sel_all[spp_sel1]

# Yellow Warbler

j = 1

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]
	
##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species (genetic data for yellow warbler)  ##
wintering.assignments <- readRDS("Data/Genoscapes/YellowWarbler/YWAR_cleaned.rds")
wintering.assignments <- wintering.assignments[which(wintering.assignments$Stage=="Wintering"),]
breeding.surfaces.origen <- readRDS("Data/Genoscapes/YellowWarbler/BreedingSurfaces.rds")
wintering.surfaces.origen <- readRDS("Data/Genoscapes/YellowWarbler/WinterProbSurfaces.rds") 
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

##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
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
	
##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
ptsBR <- ptsBR[toKeep,]
ptsNB <- ptsNB[toKeep,]

##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical <- vector() 
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
}
	
##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
}
	
##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
}

##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
winter.destinations <- list()
ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
for(i in 1:length(empiricalData_breedingHexagons3)){
	winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
	if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
		winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
	}
}

##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
}
	
##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
aaa <- which(lapply(winter.destinations, length) == 0)
if(length(aaa)>0){
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
	distanceMat.simulated <- distanceMat.simulated[-aaa,]
	distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
	distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
	distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
}

##  Compute the Mantel correlation coefficients	##
ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
}
	
##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
}

##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
}
	

##  Compute the null models  ##	

nullModel.r <- vector()
nullModel.dist.r <- vector()
nullModel.abund.r <- vector()
nullModel.r.longitude <- vector()
nullModel.dist.r.longitude <- vector()
nullModel.abund.r.longitude <- vector()
nullModel.r.latitude <- vector()
nullModel.dist.r.latitude <- vector()
nullModel.abund.r.latitude <- vector()

for(k in 1:1000){
	
	##  Null model 1: randomizing migratory destinations  ##
	
	# Sampling winter destinations at random across the wintering ground
	winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)

	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
	##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
	# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
	winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
	##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
	# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
	winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
}
nullModel.r.all[[spp_sel[j]]] <- nullModel.r
nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude
	
	
	
# Commong Yellowthroat

j = 2

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]
	
##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species (genetic data for common yellowthroat)  ##
wintering.surfaces.origen <- readRDS("Data/Genoscapes/CommonYellowthroat/COYE_WinterProbSurfaces.rds")
winter.data <- read.csv("Data/Genoscapes/CommonYellowthroat/WinterData.csv")
breeding_origen <- vector()
for(k in 3:length(colnames(wintering.surfaces.origen))){
	breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
}
breeding_origen <- breeding_origen[-c(1,2)]
wintering_origen <- winter.data[c(5,6)]

##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
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

##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
ptsBR <- ptsBR[toKeep,]
ptsNB <- ptsNB[toKeep,]

##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical <- vector() 
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
}
	
##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
}
	
##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
}

##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
winter.destinations <- list()
ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
for(i in 1:length(empiricalData_breedingHexagons3)){
	winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
	if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
		winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
	}
}

##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
}
	
##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
aaa <- which(lapply(winter.destinations, length) == 0)
if(length(aaa)>0){
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
	distanceMat.simulated <- distanceMat.simulated[-aaa,]
	distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
	distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
	distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
}

##  Compute the Mantel correlation coefficients	##
ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
}
	
##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
}

##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
}
	

##  Compute the null models  ##	

nullModel.r <- vector()
nullModel.dist.r <- vector()
nullModel.abund.r <- vector()
nullModel.r.longitude <- vector()
nullModel.dist.r.longitude <- vector()
nullModel.abund.r.longitude <- vector()
nullModel.r.latitude <- vector()
nullModel.dist.r.latitude <- vector()
nullModel.abund.r.latitude <- vector()

for(k in 1:1000){
	
	##  Null model 1: randomizing migratory destinations  ##
	
	# Sampling winter destinations at random across the wintering ground
	winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)

	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
	##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
	# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
	winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
	##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
	# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
	winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
}
nullModel.r.all[[spp_sel[j]]] <- nullModel.r
nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude


# Wilson's Warbler

j = 3

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species (genetic data for common yellowthroat)  ##
wintering.assignments <- readRDS("Data/Genoscapes/WilsonsWarbler/WIWA.WinterAssignment.rds")
wintering.surfaces.origen <- readRDS("Data/Genoscapes/WilsonsWarbler/WIWA_WinterProbSurfaces.rds")
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

##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
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

##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##	
toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
ptsBR <- ptsBR[toKeep,]
ptsNB <- ptsNB[toKeep,]


##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical <- vector() 
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
}
	
##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
}
	
##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
}

##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
winter.destinations <- list()
ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
for(i in 1:length(empiricalData_breedingHexagons3)){
	winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
	if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
		winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
	}
}

##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
}
	
##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
aaa <- which(lapply(winter.destinations, length) == 0)
if(length(aaa)>0){
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
	distanceMat.simulated <- distanceMat.simulated[-aaa,]
	distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
	distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
	distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
}

##  Compute the Mantel correlation coefficients	##
ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
}
	
##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
}

##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
}
	

##  Compute the null models  ##	

nullModel.r <- vector()
nullModel.dist.r <- vector()
nullModel.abund.r <- vector()
nullModel.r.longitude <- vector()
nullModel.dist.r.longitude <- vector()
nullModel.abund.r.longitude <- vector()
nullModel.r.latitude <- vector()
nullModel.dist.r.latitude <- vector()
nullModel.abund.r.latitude <- vector()

for(k in 1:1000){
	
	##  Null model 1: randomizing migratory destinations  ##
	
	# Sampling winter destinations at random across the wintering ground
	winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)

	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
	##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
	# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
	winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
	##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
	# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
	winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
}
nullModel.r.all[[spp_sel[j]]] <- nullModel.r
nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude



# Blue winged Warbler

j = 4

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species  ##
trackingData <- read.csv("Data/Tracking-data/Vermivora/vermNBlatWKbk.csv", header=T)
trackingData$depLON <- -trackingData$depLON
trackingData$nbLON <- -trackingData$nbLON
ptsBR <- as.matrix(trackingData[,3:2][which(trackingData$sp == 2),]) # 1 = gowwar ; 2 = buwwar
ptsNB <- as.matrix(trackingData[,4:5][which(trackingData$sp == 2),])


##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
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

##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
ptsBR <- ptsBR[toKeep,]
ptsNB <- ptsNB[toKeep,]

##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical <- vector() 
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
}
	
##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
}
	
##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
}

##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
winter.destinations <- list()
ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
for(i in 1:length(empiricalData_breedingHexagons3)){
	winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
	if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
		winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
	}
}

##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
}
	
##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
aaa <- which(lapply(winter.destinations, length) == 0)
if(length(aaa)>0){
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
	distanceMat.simulated <- distanceMat.simulated[-aaa,]
	distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
	distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
	distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
}

##  Compute the Mantel correlation coefficients	##
ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
}
	
##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
}

##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
}
	

##  Compute the null models  ##	

nullModel.r <- vector()
nullModel.dist.r <- vector()
nullModel.abund.r <- vector()
nullModel.r.longitude <- vector()
nullModel.dist.r.longitude <- vector()
nullModel.abund.r.longitude <- vector()
nullModel.r.latitude <- vector()
nullModel.dist.r.latitude <- vector()
nullModel.abund.r.latitude <- vector()

for(k in 1:1000){
	
	##  Null model 1: randomizing migratory destinations  ##
	
	# Sampling winter destinations at random across the wintering ground
	winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)

	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
	##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
	# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
	winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
	##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
	# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
	winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
}
nullModel.r.all[[spp_sel[j]]] <- nullModel.r
nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude



# Ovenbird and White-throated Sparrow
	
for(j in 5:length(spp_sel)){

	##  Load species seasonal abundance distributions  ##
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

	##  Load ORSIM output (simulated migratory connectivity)  ##
	ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

	##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species  ##
	ptsBR <- vector()
	ptsNB <- vector()
	# Banding data
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Data/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
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
	# Tracking data
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Data/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}

	##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
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

	##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

		##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical <- vector() 
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	
	##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	
	##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}

	##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
	winter.destinations <- list()
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
		if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}

	##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	
	##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	
	##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulated <- distanceMat.simulated[-aaa,]
		distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
		distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}

	##  Compute the Mantel correlation coefficients	##
	ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}

	##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
	

	##  Compute the null models  ##	

	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	nullModel.r.longitude <- vector()
	nullModel.dist.r.longitude <- vector()
	nullModel.abund.r.longitude <- vector()
	nullModel.r.latitude <- vector()
	nullModel.dist.r.latitude <- vector()
	nullModel.abund.r.latitude <- vector()
	
	for(k in 1:1000){
	
		##  Null model 1: randomizing migratory destinations  ##
	
		# Sampling winter destinations at random across the wintering ground
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
		##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
		# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
		##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
		# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
	nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
	nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
	nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
	nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
	nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
	nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude
}




	
##  Run null models for the six species shown in Fig S1: Willow Flycatcher, Brown Trasher, Gray Catbird, Brown-headed Cowbird, Common Grackle, adn Red-winged Blackbird  ##

# Selected species
spp_sel1 <- c(27,11,19,10,12,14)
spp_sel <- spp_sel_all[spp_sel1]

# Willow Flycatcher

j = 1

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species  ##
wintering.assignments <- readRDS("Data/Genoscapes/WillowFlycatcher/winter-migrant-assignments.rds")
wintering.surfaces.origen <- readRDS("Data/Genoscapes/WillowFlycatcher/WIFL_WinterProbSurfaces.rds")
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

##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
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

##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
ptsBR <- ptsBR[toKeep,]
ptsNB <- ptsNB[toKeep,]

##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical <- vector() 
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
}
	
##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
}
	
##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
}

##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
winter.destinations <- list()
ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
for(i in 1:length(empiricalData_breedingHexagons3)){
	winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
	if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
		winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
	}
}

##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
}
	
##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
aaa <- which(lapply(winter.destinations, length) == 0)
if(length(aaa)>0){
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
	distanceMat.simulated <- distanceMat.simulated[-aaa,]
	distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
	distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
	distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
}

##  Compute the Mantel correlation coefficients	##
ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
}
	
##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
}

##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
}
	

##  Compute the null models  ##	

nullModel.r <- vector()
nullModel.dist.r <- vector()
nullModel.abund.r <- vector()
nullModel.r.longitude <- vector()
nullModel.dist.r.longitude <- vector()
nullModel.abund.r.longitude <- vector()
nullModel.r.latitude <- vector()
nullModel.dist.r.latitude <- vector()
nullModel.abund.r.latitude <- vector()

for(k in 1:1000){
	
	##  Null model 1: randomizing migratory destinations  ##
	
	# Sampling winter destinations at random across the wintering ground
	winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)

	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
	##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
	# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
	winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
	##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
	# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
	winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
}
nullModel.r.all[[spp_sel[j]]] <- nullModel.r
nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude


# Brown Trasher, Gray Catbird, Brown-headed Cowbird, Common Grackle, adn Red-winged Blackbird

for(j in 2:length(spp_sel)){

	##  Load species seasonal abundance distributions  ##
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

	##  Load ORSIM output (simulated migratory connectivity)  ##
	ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

	##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species  ##
	ptsBR <- vector()
	ptsNB <- vector()
	# Banding data
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Data/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
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
	# Tracking data
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Data/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}

	##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
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

	##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

	##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical <- vector() 
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	
	##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	
	##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}

	##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
	winter.destinations <- list()
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
		if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}

	##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	
	##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	
	##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulated <- distanceMat.simulated[-aaa,]
		distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
		distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}

	##  Compute the Mantel correlation coefficients	##
	ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}

	##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
	

	##  Compute the null models  ##	

	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	nullModel.r.longitude <- vector()
	nullModel.dist.r.longitude <- vector()
	nullModel.abund.r.longitude <- vector()
	nullModel.r.latitude <- vector()
	nullModel.dist.r.latitude <- vector()
	nullModel.abund.r.latitude <- vector()
	
	for(k in 1:1000){
	
		##  Null model 1: randomizing migratory destinations  ##
	
		# Sampling winter destinations at random across the wintering ground
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
		##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
		# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
		##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
		# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
	nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
	nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
	nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
	nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
	nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
	nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude
}





##  Run null models for the six species shown in Fig S1: Common Loon, Tree Swallow, Barn Swallow, American Kestrel, Burrowing Owl, and Osprey  ##

# Selected species
spp_sel1 <- c(7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]
	
# Common Loon

j = 1

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species from genetic data  ##
wintering.assignments <- readRDS("Data/Genoscapes/CommonLoon/COLO.WinteringAssignment.rds")
wintering.surfaces.origen <- readRDS("Data/Genoscapes/CommonLoon/COLO_WinterProbSurfaces.rds")
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

##  Get hexagons containing the empirical breeding and wintering locations of individuals from genetic data  ## 
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

##  Get empirical migration data for the species from banding data  ##
ptsBR2 <- vector()
ptsNB2 <- vector()
bandingData <- read.csv("Data/Banding-data/COLO_combine.csv", header=T)
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

##  Get hexagons containing the empirical breeding and wintering locations of individuals from banding data  ## 
dfBR = data.frame(a = 1:length(ptsBR2[,1]))
sp.ptsBR <- SpatialPointsDataFrame(ptsBR2, dfBR)
proj4string(sp.ptsBR) <- sr
empiricalData_breedingHexagons_band <- over(sp.ptsBR, hexgrid3_stem[which(summer.abundance>0),])
dfNB = data.frame(a = 1:length(ptsNB2[,1]))
sp.ptsNB <- SpatialPointsDataFrame(ptsNB2, dfNB)
proj4string(sp.ptsNB) <- sr
empiricalData_nonbreedingHexagons_band <- over(sp.ptsNB, hexgrid3_stem[which(winter.abundance>0),])

##  Combine hexagons containing the empirical breeding and wintering locations of individuals from both genetic and banding data  ##
empiricalData_breedingHexagons <- c(empiricalData_breedingHexagons_gen, empiricalData_breedingHexagons_band)
empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
empiricalData_nonbreedingHexagons <- c(empiricalData_nonbreedingHexagons_gen, empiricalData_nonbreedingHexagons_band)
empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])
ptsBR <- rbind(ptsBR1, ptsBR2)
ptsNB <- rbind(ptsNB1, ptsNB2)

##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
ptsBR <- ptsBR[toKeep,]
ptsNB <- ptsNB[toKeep,]

##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical <- vector() 
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
}
	
##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
}
	
##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
}

##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
winter.destinations <- list()
ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
for(i in 1:length(empiricalData_breedingHexagons3)){
	winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
	if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
		winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
	}
}

##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
}
	
##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
distanceMat.simulated.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
}
	
aaa <- which(lapply(winter.destinations, length) == 0)
if(length(aaa)>0){
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
	distanceMat.simulated <- distanceMat.simulated[-aaa,]
	distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
	distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
	distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
}

##  Compute the Mantel correlation coefficients	##
ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
}
	
##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.longitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
}

##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
distanceMat.simulated.all.latitude <- vector()
for(i in 1:length(empiricalData_breedingHexagons3)){
	distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
}
	

##  Compute the null models  ##	

nullModel.r <- vector()
nullModel.dist.r <- vector()
nullModel.abund.r <- vector()
nullModel.r.longitude <- vector()
nullModel.dist.r.longitude <- vector()
nullModel.abund.r.longitude <- vector()
nullModel.r.latitude <- vector()
nullModel.dist.r.latitude <- vector()
nullModel.abund.r.latitude <- vector()

for(k in 1:1000){
	
	##  Null model 1: randomizing migratory destinations  ##
	
	# Sampling winter destinations at random across the wintering ground
	winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)

	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.null.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
	##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
	# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
	winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nulldist.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
	##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
	# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
	winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
	# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
	}
		
	# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
	}
		
	# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
	distanceMat.nullabund.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
	}
		
	# Compute the Mantel correlation coefficients
	nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
}
nullModel.r.all[[spp_sel[j]]] <- nullModel.r
nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude


# Tree Swallow, Barn Swallow, American Kestrel, Burrowing Owl, and Osprey	

for(j in 2:length(spp_sel)){

	##  Load species seasonal abundance distributions  ##
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

	##  Load ORSIM output (simulated migratory connectivity)  ##
	ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

	##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species  ##
	ptsBR <- vector()
	ptsNB <- vector()
	# Banding data
	if(is.na(spp_banding[spp_sel[j]])==F){
		bandingData <- read.csv(paste("Data/Banding-data/", spp_banding[spp_sel[j]], "_combine.csv", sep=""), header=T)
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
	# Tracking data
	if(is.na(spp_tracking[spp_sel[j]])==F){
		trackingData <- read.csv("Data/Tracking-data/data.csv", header=T)
		trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[j]]),]
		ptsBR <- rbind(ptsBR, cbind(trackingData_sel$blon, trackingData_sel$blat))
		ptsNB <- rbind(ptsNB, cbind(trackingData_sel$wlon1, trackingData_sel$wlat1))
	}

	##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
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

	##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
	toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
	empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
	empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
	ptsBR <- ptsBR[toKeep,]
	ptsNB <- ptsNB[toKeep,]

		##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical <- vector() 
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , 	hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
	}
	
	##  Compute pairwise longitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1][empiricalData_nonbreedingHexagons3] )
	}
	
	##  Compute pairwise latitudinal difference between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
	distanceMat.empirical.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){	
		distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2][empiricalData_nonbreedingHexagons3] )
	}

	##  Get wintering sites (hexagon) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate  ##
	winter.destinations <- list()
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3,]
	for(i in 1:length(empiricalData_breedingHexagons3)){
		winter.destinations[[i]] <- hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel[i,] > 0),]
		if(length(which(ORSIM_flows_sel[i,] > 0)) > 1){
			winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
		}
	}

	##  Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
	}
	
	##  Compute pairwise longitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.longitude <- rbind(distanceMat.simulated.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	
	##  Compute pairwise latitudinal difference between empirical breeding sites and corresponding simulated wintering sites  ##
	distanceMat.simulated.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.latitude <- rbind(distanceMat.simulated.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - unlist(winter.destinations) )
	}
	
	aaa <- which(lapply(winter.destinations, length) == 0)
	if(length(aaa)>0){
		empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons3[-aaa]
		distanceMat.simulated <- distanceMat.simulated[-aaa,]
		distanceMat.simulated.longitude <- distanceMat.simulated.longitude[-aaa,]
		distanceMat.simulated.latitude <- distanceMat.simulated.latitude[-aaa,]
		distanceMat.empirical <- distanceMat.empirical[-aaa,-aaa]
	}

	##  Compute the Mantel correlation coefficients	##
	ORSIM.r[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	ORSIM.r.longitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
	ORSIM.r.latitude[spp_sel[j]] <- mantel.rtest(as.dist(distanceMat.simulated.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),] ) / 1000)
	}
	
	##  Compute pariwise longitudinal difference between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all.longitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),1] )
	}

	##  Compute pariwise latitudinal difference between empirical breeding sites and all wintering sites 	##	
	distanceMat.simulated.all.latitude <- vector()
	for(i in 1:length(empiricalData_breedingHexagons3)){
		distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[which(winter.abundance>0),2] )
	}
	

	##  Compute the null models  ##	

	nullModel.r <- vector()
	nullModel.dist.r <- vector()
	nullModel.abund.r <- vector()
	nullModel.r.longitude <- vector()
	nullModel.dist.r.longitude <- vector()
	nullModel.abund.r.longitude <- vector()
	nullModel.r.latitude <- vector()
	nullModel.dist.r.latitude <- vector()
	nullModel.abund.r.latitude <- vector()
	
	for(k in 1:1000){
	
		##  Null model 1: randomizing migratory destinations  ##
	
		# Sampling winter destinations at random across the wintering ground
		winter.destinations.null <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T)
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null <- rbind(distanceMat.null, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.null,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null.longitude <- rbind(distanceMat.null.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.null.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.null.latitude <- rbind(distanceMat.null.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.null,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.r[k] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
		
		##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
	
		# Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
		winter.destinations.nulldist <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist <- rbind(distanceMat.nulldist, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nulldist,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist.longitude <- rbind(distanceMat.nulldist.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nulldist.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nulldist.latitude <- rbind(distanceMat.nulldist.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nulldist,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.dist.r[k] <- mantel.rtest(as.dist(distanceMat.nulldist, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.dist.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.dist.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nulldist.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	
	
		##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
	
		# Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
		winter.destinations.nullabund <- sample(which(winter.abundance>0), length(empiricalData_breedingHexagons3), replace=T, prob=winter.abundance[which(winter.abundance>0)])
	
		# Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund <- rbind(distanceMat.nullabund, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[winter.destinations.nullabund,] ) / 1000)
		}
		
		# Compute pairwise longitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund.longitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund.longitude <- rbind(distanceMat.nullabund.longitude, hexgrid3_stem_centroids[which(summer.abundance>0),1][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,1] )
		}
		
		# Compute pairwise latitudinal difference between empirical breeding sites and corresponding sampled wintering sites
		distanceMat.nullabund.latitude <- vector()
		for(i in 1:length(empiricalData_breedingHexagons3)){
			distanceMat.nullabund.latitude <- rbind(distanceMat.nullabund.latitude, hexgrid3_stem_centroids[which(summer.abundance>0),2][empiricalData_breedingHexagons3[i]] - hexgrid3_stem_centroids[winter.destinations.nullabund,2] )
		}
		
		# Compute the Mantel correlation coefficients
		nullModel.abund.r[k] <- mantel.rtest(as.dist(distanceMat.nullabund, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
		nullModel.abund.r.longitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
		nullModel.abund.r.latitude[k] <- mantel.rtest(as.dist(distanceMat.nullabund.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
	}
	nullModel.r.all[[spp_sel[j]]] <- nullModel.r
	nullModel.abund.r.all[[spp_sel[j]]] <- nullModel.abund.r
	nullModel.dist.r.all[[spp_sel[j]]] <- nullModel.dist.r	
	nullModel.r.all.longitude[[spp_sel[j]]] <- nullModel.r.longitude
	nullModel.abund.r.all.longitude[[spp_sel[j]]] <- nullModel.abund.r.longitude
	nullModel.dist.r.all.longitude[[spp_sel[j]]] <- nullModel.dist.r.longitude	
	nullModel.r.all.latitude[[spp_sel[j]]] <- nullModel.r.latitude
	nullModel.abund.r.all.latitude[[spp_sel[j]]] <- nullModel.abund.r.latitude
	nullModel.dist.r.all.latitude[[spp_sel[j]]] <- nullModel.dist.r.latitude
}




##  Comparing the performance of ORSIM to the null models  ##

# Select all species
spp_sel1 <- c(1,3,4,6,9,13,20,26,28,22,16,15,27,11,19,10,12,14,7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]

ORSIM.r2 <- ORSIM.r[spp_sel]
nullModel.r.all2 <- nullModel.r.all[spp_sel]
nullModel.abund.r.all2 <- nullModel.abund.r.all[spp_sel]
nullModel.dist.r.all2 <- nullModel.dist.r.all[spp_sel]

# Compute the p-values (i.e. fraction of null simulations with a Mantel correlation higher than ORSIM)
pval.nullModel <- vector()
pval.nullModel.abund <- vector()
pval.nullModel.dist <- vector()
for(i in 1:length(ORSIM.r2)){
	pval.nullModel[i] <- length(which(nullModel.r.all2[[i]] > ORSIM.r2[i])) / length(nullModel.r.all2[[i]])
	pval.nullModel.abund[i] <- length(which(nullModel.abund.r.all2[[i]] > ORSIM.r2[i]))/length(nullModel.abund.r.all2[[i]])
	pval.nullModel.dist[i] <- length(which(nullModel.dist.r.all2[[i]] > ORSIM.r2[i]))/length(nullModel.dist.r.all2[[i]])
}


# For longitudinal patterns of migration

ORSIM.r.longitude2 <- ORSIM.r.longitude[spp_sel]
nullModel.r.all.longitude2 <- nullModel.r.all.longitude[spp_sel]
nullModel.abund.r.all.longitude2 <- nullModel.abund.r.all.longitude[spp_sel]
nullModel.dist.r.all.longitude2 <- nullModel.dist.r.all.longitude[spp_sel]

# Compute the p-values (i.e. fraction of null simulations with a Mantel correlation higher than ORSIM)
pval.nullModel.longitude <- vector()
pval.nullModel.abund.longitude <- vector()
pval.nullModel.dist.longitude <- vector()
for(i in 1:length(ORSIM.r.longitude2)){
	pval.nullModel.longitude[i] <- length(which(nullModel.r.all.longitude2[[i]] > ORSIM.r.longitude2[i])) / length(nullModel.r.all.longitude2[[i]])
	pval.nullModel.abund.longitude[i] <- length(which(nullModel.abund.r.all.longitude2[[i]] > ORSIM.r.longitude2[i]))/length(nullModel.abund.r.all.longitude2[[i]])
	pval.nullModel.dist.longitude[i] <- length(which(nullModel.dist.r.all.longitude2[[i]] > ORSIM.r.longitude2[i]))/length(nullModel.dist.r.all.longitude2[[i]])
}


# For latitudinal patterns of migration

ORSIM.r.latitude2 <- ORSIM.r.latitude[spp_sel]
nullModel.r.all.latitude2 <- nullModel.r.all.latitude[spp_sel]
nullModel.abund.r.all.latitude2 <- nullModel.abund.r.all.latitude[spp_sel]
nullModel.dist.r.all.latitude2 <- nullModel.dist.r.all.latitude[spp_sel]

# Compute the p-values (i.e. fraction of null simulations with a Mantel correlation higher than ORSIM)
pval.nullModel.latitude <- vector()
pval.nullModel.abund.latitude <- vector()
pval.nullModel.dist.latitude <- vector()
for(i in 1:length(ORSIM.r.latitude2)){
	pval.nullModel.latitude[i] <- length(which(nullModel.r.all.latitude2[[i]] > ORSIM.r.latitude2[i])) / length(nullModel.r.all.latitude2[[i]])
	pval.nullModel.abund.latitude[i] <- length(which(nullModel.abund.r.all.latitude2[[i]] > ORSIM.r.latitude2[i]))/length(nullModel.abund.r.all.latitude2[[i]])
	pval.nullModel.dist.latitude[i] <- length(which(nullModel.dist.r.all.latitude2[[i]] > ORSIM.r.latitude2[i]))/length(nullModel.dist.r.all.latitude2[[i]])
}

