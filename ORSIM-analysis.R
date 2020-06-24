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




##  Load and save species seasonal abundance distribution (estimated from eBird data using spatio-tempororal exploratory models)  ##

sp_path <- ebirdst_download(spp[1])
abunds <- load_raster("abundance_seasonal", path = sp_path)
writeRaster(abunds[[1]], paste("Data/STEMs/stems_", spp[i], "_BR.tif", sep=""), overwrite = T)
writeRaster(abunds[[2]], paste("Data/STEMs/stems_", spp[i], "_NB.tif", sep=""), overwrite = T)
srr <- proj4string(abunds[[1]])
hexgrid3_stem_2 <- spTransform(hexgrid3_stem, srr)  # Project the hexagonal grid using the same projection than the eBird STEM rasters
for(i in 2:length(spp)){
	sp_path <- ebirdst_download(spp[i])
	abunds <- load_raster("abundance_seasonal", path = sp_path)
	writeRaster(abunds[[1]], paste("Data/STEMs/stems_", spp[i], "_BR.tif", sep=""), overwrite = T)
	writeRaster(abunds[[2]], paste("Data/STEMs/stems_", spp[i], "_NB.tif", sep=""), overwrite = T)
	print(i)
}


##  Convert species seasonal abundance rasters to values on the grid of hexagons  ##

for(i in 1:length(spp)){

	# Breeding season
	abunds_BR <- raster(paste("Data/STEMs/stems_", spp[i], "_BR.tif", sep=""))
	abundance.hex <- extract(abunds_BR, hexgrid3_stem_2, fun=mean, na.rm=T)
	abundance.hex[which(is.na(abundance.hex)==T)] <- 0
	write.csv(abundance.hex, paste("Data/STEMs/seasonalAbundance_", spp[i], "_BR.csv", sep=""), row.names=F)
	
	# Non-breeding season
	abunds_NB <- raster(paste("Data/STEMs/stems_", spp[i], "_NB.tif", sep=""))
	abundance.hex <- extract(abunds_NB, hexgrid3_stem_2, fun=mean, na.rm=T)
	abundance.hex[which(is.na(abundance.hex)==T)] <- 0
	write.csv(abundance.hex, paste("Data/STEMs/seasonalAbundance_", spp[i], "_NB.csv", sep=""), row.names=F)
}


##  Compute matrix of pairwise distance between every pair of hexagons on the grid  ##
distance.matrix <- vector()
for(i in 1:length(hexgrid3_stem_centroids[,1])){	
	distance.matrix <- rbind(distance.matrix, distHaversine(hexgrid3_stem_centroids[i,] , hexgrid3_stem_centroids)/1000)
}
write.table(distance.matrix, "ideal-optimal-redistribution/distanceMatrix.csv", row.names=F, col.names=F, sep=";")



##  Execute python script ORSIM.py on the command line to run the ORSIM model and generate simulated migratory connectivity for a given species  ##


##  FIGURE 1: Ideal Optimal Redistribution Simulator (using ovenbird as an example)  ##

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv("Data/STEMs/seasonalAbundance_ovenbi1_BR.csv", header=F)[,1]
winter.abundance <- read.csv("Data/STEMs/seasonalAbundance_ovenbi1_NB.csv", header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv("ORSIM-outputs/ORSIMresults_ovenbi1.csv", header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

##  Plot the top two panels  ##

pdf("Manuscript/Figures/Fig1a.pdf", width=8, height=3)
par(mfrow=c(1,2), mar=c(0.1,0.1,0.1,0.1), mgp=c(2,1,0), bg="white")

# (a) Seasonal relative abundance
map("world", fill=T, col="light grey", border="light grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -30), ylim=c(0,70))
rbPal.BR <- colorRampPalette(c("yellow2","red3"))
rbPal.NB <- colorRampPalette(c("skyblue","dark blue"))
datcol.BR <- rbPal.BR(5)[as.numeric(cut(summer.abundance, breaks=seq(range(summer.abundance, na.rm=T)[1], range(summer.abundance, na.rm=T)[2], range(summer.abundance, na.rm=T)[2]/5)))]
datcol.BR[which(summer.abundance == range(summer.abundance, na.rm=T)[2])] <- "red3"
datcol.BR[which(summer.abundance == "NaN" | is.na(summer.abundance) == T | summer.abundance == 0)] <- "light grey"
datcol.NB <- rbPal.NB(5)[as.numeric(cut(winter.abundance, breaks=seq(range(winter.abundance, na.rm=T)[1], range(winter.abundance, na.rm=T)[2], range(winter.abundance, na.rm=T)[2]/5)))]
datcol.NB[which(winter.abundance == range(winter.abundance, na.rm=T)[2])] <- "dark blue"
datcol.NB[which(winter.abundance == "NaN" | is.na(winter.abundance) == T | winter.abundance == 0)] <- "light grey"
datcol <- datcol.BR
datcol[which(datcol.BR == "light grey" & datcol.NB != "light grey")] <- datcol.NB[which(datcol.BR == "light grey" & datcol.NB != "light grey")]
plot(hexgrid3_stem, col=datcol, border=datcol, add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][251,], col=datcol[which(summer.abundance>0)][251], border="black", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][125,], col=datcol[which(winter.abundance>0)][125], border="black", add=T)	
spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][125,], hexgrid3_stem_centroids[which(summer.abundance>0),][251,] )),ID="a")))
plot(spl, add=T, col="black", lwd=2, lty=2)
legend("bottomleft", inset=0.06, bg="white", box.col="white", c("> 4.14","3.11–4.14","2.07–3.11","1.04–2.07", "< 1.04"), fill=rev(rbPal.BR(5)), cex=0.8, border=rev(rbPal.BR(5)), title="Breeding\nabundance")
title(main="(a) Seasonal relative abundance", line=0.7, cex.main=1)

# (b) Bipartite graph of potential movement
plot(1, 2, pch=21, cex=2, col=NULL, bg="yellow2", xlab="", ylab="", xlim=c(0.5,7.5), ylim=c(0.5,2.5), axes=F)
for(i in 1:7){
	for(j in 1:7){
		segments(i, 1, j, 2, col="grey")
	}
}
segments(3, 1, 3, 2, col="black", lwd=2.5, lty=2)
points(1, 2, pch=21, cex=3, col= NULL, bg="yellow2")
points(2, 2, pch=21, cex=3, col= NULL, bg="orange")
points(3, 2, pch=21, cex=3, col= "black", bg="red")
points(4, 2, pch=21, cex=3, col= NULL, bg="red3")
points(5, 2, pch=21, cex=3, col= NULL, bg="orange")
points(6, 2, pch=21, cex=3, col= NULL, bg="yellow2")
points(7, 2, pch=21, cex=3, col= NULL, bg="yellow2")

points(1, 1, pch=21, cex=3, col= NULL, bg="skyblue")
points(2, 1, pch=21, cex=3, col= NULL, bg="skyblue")
points(3, 1, pch=21, cex=3, col= "black", bg="skyblue3")
points(4, 1, pch=21, cex=3, col= NULL, bg="skyblue3")
points(5, 1, pch=21, cex=3, col= NULL, bg="skyblue4")
points(6, 1, pch=21, cex=3, col= NULL, bg="dark blue")
points(7, 1, pch=21, cex=3, col= NULL, bg="skyblue3")

mtext(expression(paste("(b"["i"],", k"["b"["i"]],")",sep="")), side=3, at=3, line=-2.5)
mtext(expression(paste("(w"["i"],", k"["w"["i"]],")",sep="")), side=1, at=3, line=-2.2)
mtext(expression("d"["ij"]), side=1, at=2.8, line=-6)

title(main="(b) Bipartite graph of potential movement", line=0.7, cex.main=1)
dev.off() 

##  Plot the bottom two panels  ##

pdf("Manuscript/Figures/Fig1b.pdf", width=8, height=3)
par(mfrow=c(1,2), mar=c(0.1,0.1,0.1,0.1), mgp=c(2,1,0), bg="white")

# (d) Simulated migratory connectivity
map("world", fill=T, col="light grey", border="light grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -30), ylim=c(0,70))
rbPal.BR <- colorRampPalette(c("yellow2","red3"))
rbPal.NB <- colorRampPalette(c("skyblue","dark blue"))
datcol.BR <- rbPal.BR(5)[as.numeric(cut(summer.abundance, breaks=seq(range(summer.abundance, na.rm=T)[1], range(summer.abundance, na.rm=T)[2], range(summer.abundance, na.rm=T)[2]/5)))]
datcol.BR[which(summer.abundance == range(summer.abundance, na.rm=T)[2])] <- "red3"
datcol.BR[which(summer.abundance == "NaN" | is.na(summer.abundance) == T | summer.abundance == 0)] <- "light grey"
datcol.NB <- rbPal.NB(5)[as.numeric(cut(winter.abundance, breaks=seq(range(winter.abundance, na.rm=T)[1], range(winter.abundance, na.rm=T)[2], range(winter.abundance, na.rm=T)[2]/5)))]
datcol.NB[which(winter.abundance == range(winter.abundance, na.rm=T)[2])] <- "dark blue"
datcol.NB[which(winter.abundance == "NaN" | is.na(winter.abundance) == T | winter.abundance == 0)] <- "light grey"
datcol <- datcol.BR
datcol[which(datcol.BR == "light grey" & datcol.NB != "light grey")] <- datcol.NB[which(datcol.BR == "light grey" & datcol.NB != "light grey")]
plot(hexgrid3_stem, col=datcol, border=datcol, add=T)
sampled.links <- which(EMD_flows>0, arr.ind=T)[sample(1:length(which(EMD_flows>0)), 100, replace=F),]
for(i in 1:length(sampled.links[,1])){
	spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][sampled.links[i,2],], hexgrid3_stem_centroids[which(summer.abundance>0),][sampled.links[i,1],] )),ID="a")))
	plot(spl, add=T, col="black", lwd=EMD_flows[sampled.links[i,1], sampled.links[i,2]]*250)
}
legend("bottomleft", inset=0.06, bg="white", box.col="white", c("> 0.52","0.39–0.52","0.26–0.39","0.13–0.26", "< 0.13"), fill=rev(rbPal.NB(5)), cex=0.8, border=rev(rbPal.NB(5)), title="Winter\nabundance")
title(main="(d) Simulated migratory connectivity", line=0.7, cex.main=1)

# (c) Bipartite graph of realised movement

plot(1, 2, pch=21, cex=2, col=NULL, bg="yellow2", xlab="", ylab="", xlim=c(0.5,7.5), ylim=c(0.5,2.5), axes=F)

segments(1, 1, 2, 2, col="black", lwd=1)
segments(2, 1, 2, 2, col="black", lwd=1)
segments(3, 1, 1, 2, col="black", lwd=1)
segments(3, 1, 3, 2, col="black", lwd=2)
segments(4, 1, 3, 2, col="black", lwd=2)
segments(4, 1, 5, 2, col="black", lwd=1)
segments(5, 1, 4, 2, col="black", lwd=3)
segments(6, 1, 4, 2, col="black", lwd=5)
segments(6, 1, 5, 2, col="black", lwd=2)
segments(7, 1, 6, 2, col="black", lwd=1)
segments(7, 1, 7, 2, col="black", lwd=1)

points(1, 2, pch=21, cex=3, col= NULL, bg="yellow2")
points(2, 2, pch=21, cex=3, col= NULL, bg="orange")
points(3, 2, pch=21, cex=3, col= "black", bg="red")
points(4, 2, pch=21, cex=3, col= NULL, bg="red3")
points(5, 2, pch=21, cex=3, col= NULL, bg="orange")
points(6, 2, pch=21, cex=3, col= NULL, bg="yellow2")
points(7, 2, pch=21, cex=3, col= NULL, bg="yellow2")

points(1, 1, pch=21, cex=3, col= NULL, bg="skyblue")
points(2, 1, pch=21, cex=3, col= NULL, bg="skyblue")
points(3, 1, pch=21, cex=3, col= "black", bg="skyblue3")
points(4, 1, pch=21, cex=3, col= NULL, bg="skyblue3")
points(5, 1, pch=21, cex=3, col= NULL, bg="skyblue4")
points(6, 1, pch=21, cex=3, col= NULL, bg="dark blue")
points(7, 1, pch=21, cex=3, col= NULL, bg="skyblue3")

mtext(expression(paste("(b"["i"],", k"["b"["i"]],")",sep="")), side=3, at=3, line=-2.5)
mtext(expression(paste("(w"["i"],", k"["w"["i"]],")",sep="")), side=1, at=3, line=-2.2)
mtext(expression("f"["ij"]), side=1, at=2.8, line=-6)

title(main="(c) Bipartite graph of realised movement", line=0.7, cex.main=1)

dev.off()







##  FIGURE 2: Empirical versus simulated migratory connectivity (for three speices: ovenbird, purple finch, yellow warbler)  ##
	

# Select species to map
spp_sel1 <- c(16,13,20)
spp_sel <- spp_sel_all[spp_sel1]

# Plot the figure
pdf("Manuscript/Figures/Fig2_mantel.pdf", width=12, height=10)
par(mfrow=c(3,3), mar=c(0.1,0.1,0.1,0.1), mgp=c(2,1,0), bg="white")

# Ovenbird

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[1]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[1]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[1]], ".csv", sep=""), header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species (banding and tracking data for ovenbird)  ##
ptsBR <- vector()
ptsNB <- vector()
# Banding data
if(is.na(spp_banding[spp_sel[1]])==F){
	bandingData <- read.csv(paste("Data/Banding-data/", spp_banding[spp_sel[1]], "_combine.csv", sep=""), header=T)
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
if(is.na(spp_tracking[spp_sel[1]])==F){
	trackingData <- read.csv("Data/Tracking-data/data.csv", header=T)
	trackingData_sel <- trackingData[which(trackingData$species == spp_tracking[spp_sel[1]]),]
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
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[1]], side=1, line=-1.3, at=-135, cex=1.3)
mtext("(a)", side=3, line=0.5, at=-175, cex=1.3)
title("Empirical connectivity", line=1, cex.main=2)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2)
		}
	}	
}
mtext("(b)", side=3, line=0.5, at=-175, cex=1.3)
title("Simulated connectivity", line=1, cex.main=2)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical")
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext("(c)", side=3, line=0.5, at=700, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=4400, cex=1.4)
title("Simulated vs. empirical", line=1, cex.main=2)	
		
	
# Purple Finch

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[2]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[2]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[2]], ".csv", sep=""), header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species (banding data for purple finch)  ##
ptsBR <- vector()
ptsNB <- vector()
# Banding data
if(is.na(spp_banding[spp_sel[2]])==F){
	bandingData <- read.csv(paste("Data/Banding-data/", spp_banding[spp_sel[2]], "_combine.csv", sep=""), header=T)
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
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[2]], side=1, line=-1.3, at=-135, cex=1.3)
mtext("(d)", side=3, line=0.5, at=-175, cex=1.3)
title("Empirical connectivity", line=1, cex.main=2)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2)
		}
	}	
}
mtext("(e)", side=3, line=0.5, at=-175, cex=1.3)
title("Simulated connectivity", line=1, cex.main=2)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical")
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext("(f)", side=3, line=0.5, at=700, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=4400, cex=1.4)
title("Simulated vs. empirical", line=1, cex.main=2)		
	
	
# Yellow Warbler

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[3]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[3]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[3]], ".csv", sep=""), header=FALSE)
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
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[3]], side=1, line=-1.3, at=-135, cex=1.3)
mtext("(g)", side=3, line=0.5, at=-175, cex=1.3)
title("Empirical connectivity", line=1, cex.main=2)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2)
		}
	}	
}
mtext("(h)", side=3, line=0.5, at=-175, cex=1.3)
title("Simulated connectivity", line=1, cex.main=2)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical")
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext("(i)", side=3, line=0.5, at=700, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=4400, cex=1.4)
title("Simulated vs. empirical", line=1, cex.main=2)	
	
dev.off()	







##  Statistical tests for the effect of various variables on the Mantel correlation coefficient  ##


# Select species
spp_sel1 <- c(1,3,4,6,9,13,20,26,28,22,16,15,27,11,19,10,12,14,7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]

# Load Mantel correlation coefficients resulting from ORSIM for all species
load("ideal-optimal-redistribution/nullSimulations.RData")
ORSIM.r.sel <- ORSIM.r[spp_sel1]

# Test the effect of data source on the results 
dataSource <- c("T", "T", "BT", "B", "B", "B", "G", "G", "G", "T", "T", "B", "G", "B", "B", "B", "B", "B", "BG", "B", "T", "B", "B", "BT") # T: tracking data; B: banding data; G: genetic data; BT: mix of banding and tracking data; BG: mix of banding and genetic data
summary(aov(R2 ~ dataSource))

# Test the effect of sample size on the results 
sampleSize <- c(102, 22, 10, 872, 211, 358, 178, 91, 232, 17, 53, 18, 261, 94, 29, 530, 2737, 482, 161, 19, 26, 278, 12, 413) # Number of individuals used to map migratory connectivity
summary(lm(R2 ~ sampleSize)) 

# Compute seasonal range size, seasonal range overlap and seasonal range size difference
rangeSizeBR <- vector()
rangeSizeNB <- vector()
rangeOverlap <- vector()
for(j in 1:length(spp_sel)){
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	rangeSizeBR[j] <- length(which(summer.abundance > 0))
	rangeSizeNB[j] <- length(which(winter.abundance > 0))
	rangeOverlap[j] <- length(which(summer.abundance > 0 & winter.abundance > 0)) / length(which(summer.abundance > 0 | winter.abundance > 0))
}
rangeSizeDifference <- rangeSizeNB / rangeSizeBR

summary(lm(R2 ~ rangeSizeBR))  # Test the effect of breeding range size
summary(lm(R2 ~ rangeSizeNB))  # Test the effect of wintering range size
summary(lm(R2 ~ rangeSizeDifference))  # Test the effect of range size difference
summary(lm(R2 ~ rangeOverlap))  # Test the effect of range overlap







##  FIGURE S2: Empirical versus simulated migratory connectivity for 6 species: Wood Thrush, Swainson's Thrush, Hermit Thrush, American Robin, American Goldfinch, and Purple Finch  ##


# Selected species
spp_sel1 <- c(1,3,4,6,9,13)
spp_sel <- spp_sel_all[spp_sel1]

# Plot the figure
pdf("Manuscript/Figures/FigS2.pdf", width=12, height=20)
par(mfrow=c(6,3), mar=c(2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", "(q)", "(r)")


# Wood Thrush

j = 1

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
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)
title("Empirical connectivity", line=1, cex.main=2.2)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2)
		}
	}	
}
mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)
title("Simulated connectivity", line=1, cex.main=2.2)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
title("Simulated vs. empirical", line=1, cex.main=2.2)	

	
# Swainson's Thrush, Hermit Thrush, American Robin, American Goldfinch, and Purple Finch

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

	##  Compute the Mantel correlation coefficient	##
	ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	##  Plot the figure for this species  ##
	# Empirical migratory connectivity
	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
	for(i in 1:length(ptsBR[,1])){
		spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
		plot(spl, add=T, col="black", lwd=1.5)
	}
	mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
	mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	# Simulated migratory connectivity
	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
		if(sum(ORSIM_flows_sel)>0){
			for(k in 1:length(which(ORSIM_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	# Simulated versus empirical
	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)	
}
dev.off()	
	
	
	


##  FIGURE S3: Empirical versus simulated migratory connectivity for 6 species: Yellow Warbler, Common Yellowthroat, Wilson's Warbler, Blue winged Warbler, Ovenbird, and White-throated Sparrow  ##


# Selected species
spp_sel1 <- c(20,26,28,22,16,15)
spp_sel <- spp_sel_all[spp_sel1]
	
# Plot the figure
pdf("Manuscript/Figures/FigS3.pdf", width=12, height=20)
par(mfrow=c(6,3), mar=c(2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", "(q)", "(r)")


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
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[j]], side=1, line=-1.3, at=-135, cex=1.3)
mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)
title("Empirical connectivity", line=1, cex.main=2.2)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2.2)
		}
	}	
}
mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)
title("Simulated connectivity", line=1, cex.main=2)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,10000), ylim=c(0,10000), cex.lab=1.5)
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=8000, cex=1.4)
title("Simulated vs. empirical", line=1, cex.main=2.2)	
	
	
	
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
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2)
		}
	}	
}
mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,8000), ylim=c(0,8000), cex.lab=1.5)
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=6500, cex=1.4)



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
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2)
		}
	}	
}
mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,10000), ylim=c(0,10000), cex.lab=1.5)
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=8000, cex=1.4)


	
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
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)
title("Empirical connectivity", line=1, cex.main=2.2)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2)
		}
	}	
}
mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)
title("Simulated connectivity", line=1, cex.main=2.2)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
title("Simulated vs. empirical", line=1, cex.main=2.2)
	
	
	
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

	##  Compute the Mantel correlation coefficient	##
	ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	##  Plot the figure for this species  ##
	# Empirical migratory connectivity
	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
	for(i in 1:length(ptsBR[,1])){
		spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
		plot(spl, add=T, col="black", lwd=1.5)
	}
	mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
	mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	# Simulated migratory connectivity
	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
		if(sum(ORSIM_flows_sel)>0){
			for(k in 1:length(which(ORSIM_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	# Simulated versus empirical
	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)	
}
dev.off()		
	
	
	
	
	
	
##  FIGURE S4: Empirical versus simulated migratory connectivity for 6 species: Willow Flycatcher, Brown Trasher, Gray Catbird, Brown-headed Cowbird, Common Grackle, adn Red-winged Blackbird  ##

# Selected species
spp_sel1 <- c(27,11,19,10,12,14)
spp_sel <- spp_sel_all[spp_sel1]
	
# Plot the figure
pdf("Manuscript/Figures/FigS4.pdf", width=12, height=20)
par(mfrow=c(6,3), mar=c(2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", "(q)", "(r)")


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
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)
title("Empirical connectivity", line=1, cex.main=2.2)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2)
		}
	}	
}
mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)
title("Simulated connectivity", line=1, cex.main=2.2)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
title("Simulated vs. empirical", line=1, cex.main=2.2)	
	
	
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

	##  Compute the Mantel correlation coefficient	##
	ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	##  Plot the figure for this species  ##
	# Empirical migratory connectivity
	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
	for(i in 1:length(ptsBR[,1])){
		spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
		plot(spl, add=T, col="black", lwd=1.5)
	}
	mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
	mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	# Simulated migratory connectivity
	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
		if(sum(ORSIM_flows_sel)>0){
			for(k in 1:length(which(ORSIM_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	# Simulated versus empirical
	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)	
}
dev.off()





##  FIGURE S5: Empirical versus simulated migratory connectivity for 6 species: Common Loon, Tree Swallow, Barn Swallow, American Kestrel, Burrowing Owl, and Osprey  ##

# Selected species
spp_sel1 <- c(7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]
	
# Plot the figure
pdf("Manuscript/Figures/FigS5.pdf", width=12, height=20)
par(mfrow=c(6,3), mar=c(2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", "(q)", "(r)")


# Common Loon

j = 1

##  Load species seasonal abundance distributions  ##
summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]

##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species from genetic data  ##
#wintering.assignments <- readRDS("Data/Genoscapes/CommonLoon/COLO.WinteringAssignment.rds")
#wintering.surfaces.origen <- readRDS("Data/Genoscapes/CommonLoon/COLO_WinterProbSurfaces.rds")
#breeding_origen <- vector()
#for(k in 3:length(colnames(wintering.surfaces.origen))){
#	breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
#}
#breeding_origen <- breeding_origen[-c(1,2)]
#wintering_origen <- vector()
#for(k in 3:length(colnames(wintering.surfaces.origen))){
#	if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
#		wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
#	}else{
#		wintering_origen[k] <- NA
#	}
#}
#wintering_origen <- wintering_origen[-c(1,2)]
#toRemove <- which(is.na(wintering_origen) == T)
#if(length(toRemove)>0){
#	wintering_origen <- wintering_origen[-toRemove]
#	breeding_origen <- breeding_origen[-toRemove]
#}

##  Get hexagons containing the empirical breeding and wintering locations of individuals from genetic data  ## 
#ptsBR1 <- as.matrix(wintering.surfaces.origen[,1:2][breeding_origen,])
#sptsBR <- SpatialPoints(ptsBR1)
#proj4string(sptsBR) <- proj4string(hexgrid3_stem)
#breedingPoints_inHexagons <- gContains(hexgrid3_stem[which(summer.abundance>0),], sptsBR, byid=T)
#breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))

#ptsNB1 <- as.matrix(wintering.assignments[,5:4][wintering_origen,])
#sptsNB <- SpatialPoints(ptsNB1)
#proj4string(sptsNB) <- proj4string(hexgrid3_stem)
#winteringPoints_inHexagons <- gContains(hexgrid3_stem[which(winter.abundance>0),], sptsNB, byid=T)
#winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))

#toKeep <- which(lapply(breedingHexagons, length) > 0 & lapply(winteringHexagons, length) > 0)
#empiricalData_breedingHexagons_gen <- unlist(breedingHexagons[toKeep])
#empiricalData_nonbreedingHexagons_gen <- unlist(winteringHexagons[toKeep])

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
#empiricalData_breedingHexagons <- c(empiricalData_breedingHexagons_gen, empiricalData_breedingHexagons_band)
#empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
#empiricalData_nonbreedingHexagons <- c(empiricalData_nonbreedingHexagons_gen, empiricalData_nonbreedingHexagons_band)
#empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])
#ptsBR <- rbind(ptsBR1, ptsBR2)
#ptsNB <- rbind(ptsNB1, ptsNB2)
empiricalData_breedingHexagons <- empiricalData_breedingHexagons_band
empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
empiricalData_nonbreedingHexagons <- empiricalData_nonbreedingHexagons_band
empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])
ptsBR <- ptsBR2
ptsNB <- ptsNB2

##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
ptsBR <- ptsBR[toKeep,]
ptsNB <- ptsNB[toKeep,]

##  Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical dataset  ##
distanceMat.empirical <- vector() 
for(i in 1:length(empiricalData_breedingHexagons3)){	
	distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine( hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] , hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3,] ) / 1000)
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

##  Compute the Mantel correlation coefficient	##
ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
##  Plot the figure for this species  ##
# Empirical migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(ptsBR[,1])){
	spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)
title("Empirical connectivity", line=1, cex.main=2.2)

# Simulated migratory connectivity
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
	if(sum(ORSIM_flows_sel)>0){
		for(k in 1:length(which(ORSIM_flows_sel > 0))){
			spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
			plot(spl, add=T, col="black", lwd=1.5, cex=2)
		}
	}	
}
mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)
title("Simulated connectivity", line=1, cex.main=2.2)

# Simulated versus empirical
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
axis(side=1)
axis(side=2)
abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
title("Simulated vs. empirical", line=1, cex.main=2.2)	

	

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

	##  Compute the Mantel correlation coefficient	##
	ORSIM.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	##  Plot the figure for this species  ##
	# Empirical migratory connectivity
	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-40,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
	for(i in 1:length(ptsBR[,1])){
		spl = SpatialLines(list(Lines(Line(rbind(ptsBR[i,], ptsNB[i,])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
		plot(spl, add=T, col="black", lwd=1.5)
	}
	mtext(spp_name[spp_sel[j]], side=1, line=-4, at=-135, cex=1.3)
	mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	# Simulated migratory connectivity
	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-40,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(ORSIM_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		ORSIM_flows_sel <- ORSIM_flows[empiricalData_breedingHexagons3[i],]
		if(sum(ORSIM_flows_sel)>0){
			for(k in 1:length(which(ORSIM_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(ORSIM_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	# Simulated versus empirical
	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,12000), ylim=c(0,12000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(ORSIM.rM,3))), side=1, line=-1.3, at=9700, cex=1.4)	
}
dev.off()









##  FIGURE S6: distribution of links from breeding sites in simulated migratory connectivity ##
	
# Select species	
spp_sel1 <- c(1,3,4,6,9,13,20,26,28,22,16,15,27,11,19,10,12,14,7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]

pdf("Manuscript/Figures/FigS6.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")
for(j in 1:length(spp_sel1)){
	if(is.na(spp_sel1[j])==T){
		plot.new()
	}else{
		##  Load species seasonal abundance distributions  ##
		summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
		winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

		##  Load ORSIM output (simulated migratory connectivity)  ##
		ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
		ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]
		
		links <- apply(ORSIM_flows, 1, function(x) length(which(x > 0)))
		links <- links[links>0]
		maxLinks <- max(links)
		
		hist(links, main=spp_name[spp_sel[j]], xlab="Number of links", ylab="Frequency", xlim=c(0.5,10.5), breaks=seq(0.5,maxLinks+1,1), col="grey")
		mtext(paste("max = ", maxLinks, sep=""), side=3, line=-3, at=7, cex=0.8)
	}
}
dev.off()






##  FIGURE S7: QQ-Plot for exponential distribution for the distribution of links from breeding sites  ##

pdf("Manuscript/Figures/FigS7.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")
for(j in 1:length(spp_sel1)){
	##  Load species seasonal abundance distributions  ##
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

	##  Load ORSIM output (simulated migratory connectivity)  ##
	ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	links <- apply(ORSIM_flows, 1, function(x) length(which(x > 0)))
	links <- links[links>0]
	
	qqPlot(links, distribution="exp", id=F, pch=20, cex=1.1, lwd=1.5, col.lines="red", main=spp_name[spp_sel[j]])
}
dev.off()






##  FIGURE S8: distribution of links from wintering sites in simulated migratory connectivity ##

pdf("Manuscript/Figures/FigS8.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")
for(j in 1:length(spp_sel1)){
	if(is.na(spp_sel1[j])==T){
		plot.new()
	}else{
		##  Load species seasonal abundance distributions  ##
		summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
		winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

		##  Load ORSIM output (simulated migratory connectivity)  ##
		ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
		ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]
	
		links <- apply(ORSIM_flows, 2, function(x) length(which(x > 0)))
		links <- links[links>0]
		maxLinks <- max(links)
		
		hist(links, main=spp_name[spp_sel[j]], xlab="Number of links", ylab="Frequency", xlim=c(0.5,15.5), breaks=seq(0.5,maxLinks+1,1), col="grey")
		mtext(paste("max = ", maxLinks, sep=""), side=3, line=-3, at=7, cex=0.8)
	}
}
dev.off()






##  FIGURE S9: QQ-Plot for exponential distribution for the distribution of links from wintering sites  ##

pdf("Manuscript/Figures/FigS9.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")
for(j in 1:length(spp_sel1)){
	##  Load species seasonal abundance distributions  ##
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

	##  Load ORSIM output (simulated migratory connectivity)  ##
	ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]
		
	links <- apply(ORSIM_flows, 2, function(x) length(which(x > 0)))
	links <- links[links>0]
	
	qqPlot(links, distribution="exp", id=F, pch=20, cex=1.1, lwd=1.5, col.lines="red", main=spp_name[spp_sel[j]])
}
dev.off()


	
	
	
	
	
##  Correlations between site-level migratory connectivity and local relative abundance  ##

glm.results.BR <- list()
glm.results.NB <- list()
for(j in 1:length(spp_sel1)){
	##  Load species seasonal abundance distributions  ##
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

	##  Load ORSIM output (simulated migratory connectivity)  ##
	ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	linksBR <- apply(ORSIM_flows, 1, function(x) length(which(x > 0)))
	linksNB <- apply(ORSIM_flows, 2, function(x) length(which(x > 0)))

	glm.results.BR[[j]] <- summary(glm(linksBR ~ summer.abundance[which(summer.abundance>0)]))$coefficients[2,c(1,4)]
	glm.results.NB[[j]] <- summary(glm(linksNB ~ winter.abundance[which(winter.abundance>0)]))$coefficients[2,c(1,4)]
}





##  FIGURE S10: species' relative abundance distribution  ##

pdf("Manuscript/Figures/FigS10.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")	
for(j in 1:length(spp_sel1)){
	
	##  Load species seasonal abundance distributions  ##
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

	map("world", fill=T, col="light grey", border="light grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -30), ylim=c(-40,70))
	rbPal.BR <- colorRampPalette(c("yellow2","red3"))
	rbPal.NB <- colorRampPalette(c("skyblue","dark blue"))
	datcol.BR <- rbPal.BR(6)[as.numeric(cut(summer.abundance, breaks=seq(range(summer.abundance, na.rm=T)[1], range(summer.abundance, na.rm=T)[2], range(summer.abundance, na.rm=T)[2]/6)))]
	datcol.BR[which(summer.abundance == range(summer.abundance, na.rm=T)[2])] <- "red3"
	datcol.BR[which(summer.abundance == "NaN" | is.na(summer.abundance) == T | summer.abundance == 0)] <- "light grey"
	datcol.NB <- rbPal.NB(6)[as.numeric(cut(winter.abundance, breaks=seq(range(winter.abundance, na.rm=T)[1], range(winter.abundance, na.rm=T)[2], range(winter.abundance, na.rm=T)[2]/6)))]
	datcol.NB[which(winter.abundance == range(winter.abundance, na.rm=T)[2])] <- "dark blue"
	datcol.NB[which(winter.abundance == "NaN" | is.na(winter.abundance) == T | winter.abundance == 0)] <- "light grey"
	datcol <- datcol.NB
	datcol[which(datcol.NB == "light grey" & datcol.BR != "light grey")] <- datcol.BR[which(datcol.NB == "light grey" & datcol.BR != "light grey")]
	plot(hexgrid3_stem, col=datcol, border=datcol, add=T)
	legend("bottomleft", inset=0.06, bg="white", box.col="white", legend=round(rev(seq(range(summer.abundance, na.rm=T)[1], range(summer.abundance, na.rm=T)[2], range(summer.abundance, na.rm=T)[2]/6))[1:6],2), fill=rev(rbPal.BR(6)), cex=0.85, border=rev(rbPal.BR(6)), title="Relative\nabundance")
	legend("bottomleft", inset=c(0.3,0.06), bg="white", box.col="white", legend=round(rev(seq(range(winter.abundance, na.rm=T)[1], range(winter.abundance, na.rm=T)[2], range(winter.abundance, na.rm=T)[2]/6))[1:6],2), fill=rev(rbPal.NB(6)), cex=0.85, border=rev(rbPal.NB(6)))
	title(main=spp_name[spp_sel[j]], line=0.3)
}
dev.off()




##  FIGURE S11: simulated site-level migratory connectivity  ##

pdf("Manuscript/Figures/FigS11.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")	
for(j in 1:length(spp_sel1)){
	
	##  Load species seasonal abundance distributions  ##
	summer.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Data/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]

	##  Load ORSIM output (simulated migratory connectivity)  ##
	ORSIM_results <- read.csv(paste("ORSIM-outputs/ORSIMresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	ORSIM_flows <- ORSIM_results[which(summer.abundance>0), which(winter.abundance>0)]
	
	linksBR <- apply(ORSIM_flows, 1, function(x) length(which(x > 0)))
	linksNB <- apply(ORSIM_flows, 2, function(x) length(which(x > 0)))	
	
	map("world", fill=T, col="light grey", border="light grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -30), ylim=c(-40,70))
	rbPal.BR <- colorRampPalette(c("yellow2","red3"))
	rbPal.NB <- colorRampPalette(c("skyblue","dark blue"))
	datcol.BR <- rep("light grey", length(summer.abundance))
	datcol.BR[which(summer.abundance>0)] <- rbPal.BR(6)[as.numeric(cut(linksBR, breaks=c(0.5,1.5,2.5,5.5,10.5,20.5,100)))]
	datcol.BR[which(summer.abundance>0)][which(linksBR == 0)] <- "dark grey"
	datcol.NB <- rep("light grey", length(winter.abundance))
	datcol.NB[which(winter.abundance>0)] <- rbPal.NB(6)[as.numeric(cut(linksNB, breaks=c(0.5,1.5,2.5,5.5,10.5,20.5,100)))]
	#datcol.NB[which(winter.abundance>0)] <- rbPal.NB(5)[as.numeric(cut(linksNB, breaks=seq(range(linksNB, na.rm=T)[1], range(linksNB, na.rm=T)[2], range(linksNB, na.rm=T)[2]/5)))]
	datcol.NB[which(winter.abundance>0)][which(linksNB == range(linksNB, na.rm=T)[2])] <- "dark blue"
	datcol.NB[which(winter.abundance>0)][which(linksNB == 0)] <- "dark grey"
	datcol <- datcol.NB
	datcol[which(datcol.NB == "light grey" & datcol.BR != "light grey")] <- datcol.BR[which(datcol.NB == "light grey" & datcol.BR != "light grey")]
	plot(hexgrid3_stem, col=datcol, border=datcol, add=T)
	if(j==1){
		legend("bottomleft", inset=0.06, bg="white", box.col="white", c("> 20", "11–20","6–10","3–5","2", "1"), fill=rev(rbPal.BR(6)), cex=1, border=rev(rbPal.BR(6)), title="Number\nof links")
	}
	if(j==5){
		legend("bottomleft", inset=0.06, bg="white", box.col="white", c("> 20", "11–20","6–10","3–5","2", "1"), fill=rev(rbPal.NB(6)), cex=1, border=rev(rbPal.NB(6)), title="Number\nof links")
	}
	title(main=spp_name[spp_sel[j]], line=0.3)
}
dev.off()

