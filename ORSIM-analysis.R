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

##  Get empirical migration data (location of breeding and wintering sites of individuals) for the species (banding and tracking data for ovenbird)  ##
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







##  FIGURE S1:  ##

spp_sel1 <- c(1,3,4,6,9,13)
spp_sel <- spp_sel_all[spp_sel1]
	
pdf("Redistribution-model/EMD-outputs/FigS1_mantel.pdf", width=12, height=20)
par(mfrow=c(6,3), mar=c(2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

let <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", "(q)", "(r)")

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
	
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared	
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]
	
	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)
	title("Empirical connectivity", line=1, cex.main=2.2)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)
	title("Simulated connectivity", line=1, cex.main=2.2)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
	title("Simulated vs. empirical", line=1, cex.main=2.2)	



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
	
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]

	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)

}
dev.off()	
	
	
	


##  Figure S2  ##

spp_sel1 <- c(20,26,28,22,16,15)
spp_sel <- spp_sel_all[spp_sel1]
	
pdf("Redistribution-model/EMD-outputs/FigS2_mantel.pdf", width=12, height=20)
par(mfrow=c(6,3), mar=c(2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

j=1
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	
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
	
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared	
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]
	
	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)
	title("Empirical connectivity", line=1, cex.main=2.2)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)
	title("Simulated connectivity", line=1, cex.main=2.2)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,10000), ylim=c(0,10000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
	title("Simulated vs. empirical", line=1, cex.main=2.2)	
	
	
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
		
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared	
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]
	
	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,8000), ylim=c(0,8000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
		

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
		
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared	
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]
	
	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,10000), ylim=c(0,10000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
		

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
	
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared	
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]
	
	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
	
	
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
		
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]

	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)	
	
}
dev.off()	
	
	
	
	
	
	
##  Figure S3  ##

spp_sel1 <- c(27,11,19,10,12,14)
spp_sel <- spp_sel_all[spp_sel1]
	
pdf("Redistribution-model/EMD-outputs/FigS3_mantel.pdf", width=12, height=20)
par(mfrow=c(6,3), mar=c(2,0.1,0.1,0.1), mgp=c(1.9,0.8,0), bg="white")

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
	
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared	
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]
	
	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)
	title("Empirical connectivity", line=1, cex.main=2.2)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)
	title("Simulated connectivity", line=1, cex.main=2.2)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
	title("Simulated vs. empirical", line=1, cex.main=2.2)
	

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
	
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]

	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
}
dev.off()





##  Figure S4  ##

spp_sel1 <- c(7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]
	
pdf("Redistribution-model/EMD-outputs/FigS4_mantel.pdf", width=12, height=20)
par(mfrow=c(6,3), mar=c(2,0.1,0.1,0.1), mgp=c(1.9,0.8,0), bg="white")

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
	
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared	
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]
	
	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)
	title("Empirical connectivity", line=1, cex.main=2.2)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)
	title("Simulated connectivity", line=1, cex.main=2.2)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=5700, cex=1.4)
	title("Simulated vs. empirical", line=1, cex.main=2.2)
	

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
		
	#EMD.r2 <- summary(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)))$r.squared
	#nullModel.r2 <- null.r2[[spp_sel1[j]]]

	EMD.rM <- mantel.rtest(as.dist(distanceMat.simulatedBR, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
	
	# Plots
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
	mtext(let[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.3)

	map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-40,70))
	plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
	plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
	plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[empiricalData_breedingHexagons3,], 2, sum) > 0),], col="blue", border="blue", add=T)
	for(i in 1:length(empiricalData_breedingHexagons3)){
		EMD_flows_sel <- EMD_flows[empiricalData_breedingHexagons3[i],]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.5, cex=2)
			}
		}	
	}
	mtext(let[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.3)

	plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulatedBR), pch=20, cex=0.8, axes=F, main=NULL, ylab="Simulated", xlab="Empirical", xlim=c(0,12000), ylim=c(0,12000), cex.lab=1.5)
	axis(side=1)
	axis(side=2)
	abline(lm(as.vector(distanceMat.simulatedBR) ~ as.vector(distanceMat.empirical)), col="red", lwd=2)
	mtext(let[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.3)
	mtext(bquote('r'["M"]*' = '* .(round(EMD.rM,3))), side=1, line=-1.3, at=9700, cex=1.4)

}
dev.off()










##  Fig S5: distribution of links from breeding sites in simulated migratory connectivity ##
	
spp_sel1 <- c(1,3,4,6,9,13,20,26,28,22,16,15,27,11,19,10,12,14,7,2,18,5,8,17)
spp_sel <- spp_sel_all[spp_sel1]

pdf("Redistribution-model/EMD-outputs/FigS5.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")
for(j in 1:length(spp_sel1)){
	if(is.na(spp_sel1[j])==T){
		plot.new()
	}else{
		summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
		winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
		EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
		EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
		links <- apply(EMD_flows, 1, function(x) length(which(x > 0)))
		links <- links[links>0]
		maxLinks <- max(links)
		hist(links, main=spp_name[spp_sel[j]], xlab="Number of links", ylab="Frequency", xlim=c(0.5,10.5), breaks=seq(0.5,maxLinks+1,1), col="grey")
		mtext(paste("max = ", maxLinks, sep=""), side=3, line=-3, at=7, cex=0.8)
	}
}
dev.off()


##  Fig S6: QQ-Plot for exponential distribution for the distribution of links from breeding sites  ##

pdf("Redistribution-model/EMD-outputs/FigS6.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")
for(j in 1:length(spp_sel1)){
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	links <- apply(EMD_flows, 1, function(x) length(which(x > 0)))
	links <- links[links>0]
	qqPlot(links, distribution="exp", id=F, pch=20, cex=1.1, lwd=1.5, col.lines="red", main=spp_name[spp_sel[j]])
}
dev.off()


##  Fig S7: distribution of links from wintering sites in simulated migratory connectivity ##

pdf("Redistribution-model/EMD-outputs/FigS7.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")
for(j in 1:length(spp_sel1)){
	if(is.na(spp_sel1[j])==T){
		plot.new()
	}else{
		summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
		winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
		EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
		EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
		links <- apply(EMD_flows, 2, function(x) length(which(x > 0)))
		links <- links[links>0]
		maxLinks <- max(links)
		hist(links, main=spp_name[spp_sel[j]], xlab="Number of links", ylab="Frequency", xlim=c(0.5,15.5), breaks=seq(0.5,maxLinks+1,1), col="grey")
		mtext(paste("max = ", maxLinks, sep=""), side=3, line=-3, at=7, cex=0.8)
	}
}
dev.off()


##  Fig S8: QQ-Plot for exponential distribution for the distribution of links from wintering sites  ##

pdf("Redistribution-model/EMD-outputs/FigS8.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")
for(j in 1:length(spp_sel1)){
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	links <- apply(EMD_flows, 2, function(x) length(which(x > 0)))
	links <- links[links>0]
	qqPlot(links, distribution="exp", id=F, pch=20, cex=1.1, lwd=1.5, col.lines="red", main=spp_name[spp_sel[j]])
}
dev.off()


	
##  Relative abundance versus degree   ##

#pdf("Redistribution-model/EMD-outputs/FigS10.pdf", width=9, height=12)
#par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")
glm.results.BR <- list()
glm.results.NB <- list()
for(j in 1:length(spp_sel1)){
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	linksBR <- apply(EMD_flows, 1, function(x) length(which(x > 0)))
	linksNB <- apply(EMD_flows, 2, function(x) length(which(x > 0)))
	#plot(winter.abundance[which(winter.abundance>0)], links, pch=20)		
	#abline(glm(links ~ winter.abundance[which(winter.abundance>0)]), col="red")
	glm.results.BR[[j]] <- summary(glm(linksBR ~ summer.abundance[which(summer.abundance>0)]))$coefficients[2,c(1,4)]
	glm.results.NB[[j]] <- summary(glm(linksNB ~ winter.abundance[which(winter.abundance>0)]))$coefficients[2,c(1,4)]
}
#dev.off()


##  Figs: maps of relative abundance and degree spatial distribution  ##

pdf("Redistribution-model/EMD-outputs/FigS11.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")	
for(j in 1:length(spp_sel1)){
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
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

pdf("Redistribution-model/EMD-outputs/FigS12.pdf", width=9, height=12)
par(mfrow=c(6,4), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")	
for(j in 1:length(spp_sel1)){
	summer.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_BR.csv", sep=""), header=F)[,1]
	winter.abundance <- read.csv(paste("Redistribution-model/STEMs/seasonalAbundance_", spp[spp_sel[j]], "_NB.csv", sep=""), header=F)[,1]
	EMD_results <- read.csv(paste("Redistribution-model/EMD-outputs/EMDresults_", spp[spp_sel[j]], ".csv", sep=""), header=FALSE)
	EMD_flows <- EMD_results[which(summer.abundance>0), which(winter.abundance>0)]
	linksBR <- apply(EMD_flows, 1, function(x) length(which(x > 0)))
	linksNB <- apply(EMD_flows, 2, function(x) length(which(x > 0)))	
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


legend("bottomleft", inset=0.06, bg="white", box.col="white", c("> 4.14","3.11–4.14","2.07–3.11","1.04–2.07", "< 1.04"), fill=rev(rbPal.BR(5)), cex=0.8, border=rev(rbPal.BR(5)), title="Breeding\nabundance")
mtext("(a) Seasonal relative abundance", side=3, at=-100, line=0.7, cex=1.25)

	
	
	
	
	
	
	
	
	
	




# Simulated from wintering sites
summer.destinations <- list()
EMD_flows_sel <- EMD_flows[,empiricalData_nonbreedingHexagons3]
for(i in 1:length(empiricalData_nonbreedingHexagons3)){
	summer.destinations[[i]] <- hexgrid3_stem_centroids[which(summer.abundance>0),][which(EMD_flows_sel[,i] > 0),]
	if(length(which(EMD_flows_sel[,i] > 0)) > 1){
		summer.destinations[[i]] <- apply(summer.destinations[[i]], 2, mean)
	}
	if(length(which(EMD_flows_sel[,i] > 0)) == 0){
		summer.destinations[[i]] <- NA
	}
}
empiricalData_nonbreedingHexagons4 <- empiricalData_nonbreedingHexagons3[-which(is.na(summer.destinations) == T)]
distanceMat.simulatedNB <- vector()
for(i in 1:length(empiricalData_nonbreedingHexagons4)){
	distanceMat.simulatedNB <- rbind(distanceMat.simulatedNB, distHaversine( hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons4[i],] , matrix(unlist(summer.destinations[-which(is.na(summer.destinations) == T)]), ncol=2, byrow=T) ) / 1000)
}


## PLOTS
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	spl = SpatialLines(list(Lines(Line(rbind( winter.destinations[[i]], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
	plot(spl, add=T, col="black", lwd=1.5)
}
for(i in 1:length(empiricalData_breedingHexagons3)){
	spl = SpatialLines(list(Lines(Line(rbind( winter.destinations[[i]], hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],] )),ID="a")))
	plot(spl, add=T, col="black", lwd=1.5)
}

hexgrid3_stem_centroids[which(summer.abundance>0),][empiricalData_breedingHexagons3[i],]


map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(empiricalData_breedingHexagons3),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(empiricalData_nonbreedingHexagons3),], col="blue", border="blue", add=T)
for(i in 1:length(empiricalData_breedingHexagons3)){
	spl = SpatialLines(list(Lines(Line(rbind( summer.destinations[[i]], hexgrid3_stem_centroids[which(winter.abundance>0),][empiricalData_nonbreedingHexagons3[i],] )),ID="a")))
	plot(spl, add=T, col="black", lwd=1.5)
}









##  Plot empirical and simulated data by regions  

# Empirical data
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
bandingData_breedingHexagons2_neighs <- list()
for(i in 1:length(bandingData_breedingHexagons2)){
	bandingData_breedingHexagons2_neighs[[i]] <- c(bandingData_breedingHexagons2[i], which(gTouches(hexgrid3_stem[which(summer.abundance>0),][bandingData_breedingHexagons2[i]], hexgrid3_stem[which(summer.abundance>0),], byid=T) == TRUE))
}
bandingData_nonbreedingHexagons2_neighs <- list()
for(i in 1:length(bandingData_nonbreedingHexagons2)){
	bandingData_nonbreedingHexagons2_neighs[[i]] <- c(bandingData_nonbreedingHexagons2[i], which(gTouches(hexgrid3_stem[which(winter.abundance>0),][bandingData_nonbreedingHexagons2[i]], hexgrid3_stem[which(winter.abundance>0),], byid=T) == TRUE))
}
plot(hexgrid3_stem[which(summer.abundance>0),][unique(unlist(bandingData_breedingHexagons2_neighs)),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(unlist(bandingData_nonbreedingHexagons2_neighs)),], col="blue", border="blue", add=T)
for(i in 1:length(bandingData$GISBLONG)){
	spl = SpatialLines(list(Lines(Line(rbind( c(bandingData$GISBLONG[i],bandingData$GISBLAT[i]), c(bandingData$GISRLONG[i],bandingData$GISRLAT[i]) )),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="black", lwd=1.5)
}
title("Empirical connectivity (regions)", line=1.5)

# Simulated data from breeding grounds
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][unique(unlist(bandingData_breedingHexagons2_neighs)),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[unique(unlist(bandingData_breedingHexagons2_neighs)),], 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(bandingData_breedingHexagons2)){
	for(j in 1:length(bandingData_breedingHexagons2_neighs[[i]])){
		EMD_flows_sel <- EMD_flows[bandingData_breedingHexagons2_neighs[[i]],][j,]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(summer.abundance>0),][bandingData_breedingHexagons2_neighs[[i]][j],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1.2)
			}
		}
	}	
}
title("Simulated connectivity from breeding regions", line=1.5)


# Simulated data from wintering grounds
map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][unique(unlist(bandingData_nonbreedingHexagons2_neighs)),], col="blue", border="blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][which(apply(EMD_flows[,unique(unlist(bandingData_nonbreedingHexagons2_neighs))], 1, sum) > 0),], col="red", border="red", add=T)
for(i in 1:length(bandingData_nonbreedingHexagons2)){
	for(j in 1:length(bandingData_nonbreedingHexagons2_neighs[[i]])){
		EMD_flows_sel <- EMD_flows[,bandingData_nonbreedingHexagons2_neighs[[i]]][,j]
		if(sum(EMD_flows_sel)>0){
			for(k in 1:length(which(EMD_flows_sel > 0))){
				spl = SpatialLines(list(Lines(Line(rbind( hexgrid3_stem_centroids[which(summer.abundance>0),][which(EMD_flows_sel > 0)[k],], hexgrid3_stem_centroids[which(winter.abundance>0),][bandingData_nonbreedingHexagons2_neighs[[i]][j],] )),ID="a")))
				plot(spl, add=T, col="black", lwd=1)
			}
		}
	}	
}
title("Simulated connectivity from wintering regions", line=1.5)









##  Compute distance between banding data and nearest simulated destinations  ##



minDist <- vector()
for(i in 1:length(bandingData_breedingHexagons3)){	
	EMD_flows_sel <- EMD_flows[bandingData_breedingHexagons2[i],]
	minDists1 <- vector()
	for(j in 1:length(which(bandingData_breedingHexagons3 == bandingData_breedingHexagons3[i]))){
		minDists1[j] <- min(distHaversine(hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[j,] > 0),], hexgrid3_stem_centroids[which(winter.abundance>0),][bandingData_nonbreedingHexagons3[which(bandingData_breedingHexagons3 == bandingData_breedingHexagons3[i])],]) / 1000)
	}
	EMD_flows_sel <- EMD_flows[,bandingData_nonbreedingHexagons2[i]]
	minDists2 <- vector()
	for(j in 1:length(which(bandingData_nonbreedingHexagons3 == bandingData_nonbreedingHexagons3[i]))){
		if(sum(EMD_flows_sel[,j])>0){
			minDists2[j] <- min(distHaversine(hexgrid3_stem_centroids[which(summer.abundance>0),][which(EMD_flows_sel[,j] > 0),], hexgrid3_stem_centroids[which(summer.abundance>0),][bandingData_breedingHexagons3[which(bandingData_nonbreedingHexagons3 == bandingData_nonbreedingHexagons3[i])],]) / 1000)
		}
	}
	minDist[i] <- min(c(minDists1, minDists2), na.rm=T)
}
	
	
	minDist[i] <- min(minDists)



# Between sites
minDist <- vector()
for(i in 1:length(bandingData_breedingHexagons2)){	
	EMD_flows_sel <- EMD_flows[which(bandingData_breedingHexagons == bandingData_breedingHexagons2[i]),]
	minDists <- vector()
	for(j in 1:length(which(bandingData_breedingHexagons == bandingData_breedingHexagons2[i]))){
		minDists[j] <- min(distHaversine(hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows_sel[j,] > 0),], hexgrid3_stem_centroids[which(winter.abundance>0),][bandingData_nonbreedingHexagons[which(bandingData_breedingHexagons == bandingData_breedingHexagons2[i])],]) / 1000)
	}
	minDist[i] <- min(minDists)
}
median(minDist, na.rm=T)


# Between regions









##  Plot connectivity  ##

map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="blue", border="blue", add=T)
for(i in 1:length(EMD_flows[,1])){
	if(length(which(EMD_flows[i,] > 0)) > 1){
		for(j in 1:length(which(EMD_flows[i,] > 0))){
			spl = SpatialLines(list(Lines(Line(rbind(hexgrid3_stem_centroids[which(summer.abundance>0),][i,] + sample(seq(-0.5,0.5,0.001),1), hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows[i,] > 0),][j,] + sample(seq(-0.5,0.5,0.001),1))),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
			plot(spl, add=T, col="black", lwd=0.25)
		}
	}
	if(length(which(EMD_flows[i,] > 0)) == 1){
		spl = SpatialLines(list(Lines(Line(rbind(hexgrid3_stem_centroids[which(summer.abundance>0),][i,] + sample(seq(-0.5,0.5,0.001),1), hexgrid3_stem_centroids[which(winter.abundance>0),][which(EMD_flows[i,] > 0),] + sample(seq(-0.5,0.5,0.001),1))),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
		plot(spl, add=T, col="black", lwd=0.25)
	}
}




#pdf(paste("Redistribution-model/YellowWarbler/outputs/RegionalConnectivity_new_", alpha, "_", beta, ".pdf", sep=""), width=8, height=8)
pdf("Redistribution-model/YellowWarbler/outputs/RegionalConnectivity_greatcircledistance_new.pdf", width=8, height=8)

par(mfrow=c(2,2), mar=c(0.5,0.5,0.5,0.5), mgp=c(2,1,0), bg="white")

# South west
hexsBR <- which(hexgrid3_stem_centroids[,1][which(summer.abundance>0)] <= -110 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] <= 40 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] > 25)
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="dark blue", border="dark blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][hexsBR], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[hexsBR,], 2, sum) > 0),], col="orange", border="orange", add=T)

# North west
hexsBR <- which(hexgrid3_stem_centroids[,1][which(summer.abundance>0)] <= -120 & hexgrid3_stem_centroids[,1][which(summer.abundance>0)] > -140 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] > 45 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] <= 60)
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="dark blue", border="dark blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][hexsBR], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[hexsBR,], 2, sum) > 0),], col="orange", border="orange", add=T)

# Central
hexsBR <- which(hexgrid3_stem_centroids[,1][which(summer.abundance>0)] <= -80 & hexgrid3_stem_centroids[,1][which(summer.abundance>0)] > -100 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] <= 45 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] > 25)
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="dark blue", border="dark blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][hexsBR], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[hexsBR,], 2, sum) > 0),], col="orange", border="orange", add=T)

# East
hexsBR <- which(hexgrid3_stem_centroids[,1][which(summer.abundance>0)] > -75 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] > 25)
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="dark blue", border="dark blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][hexsBR], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[hexsBR,], 2, sum) > 0),], col="orange", border="orange", add=T)

dev.off()

## Plot migration flow

# Compute least cost paths between the pairs of breeding and non-breeding sites predicted by the EMD model 
EMD_flows2 <- list()
leastCostPaths <- list()
for(i in 1:length(EMD_flows[,1])){
	leastCostPaths[[i]] <- shortest_paths(g, land[which(summer.abundance > 0)][i], to=land[which(winter.abundance > 0)])$vpath
	leastCostPaths[[i]] <- leastCostPaths[[i]][which(EMD_flows[i,] > 0)]
	EMD_flows2[[i]] <- EMD_flows[i,][which(EMD_flows[i,] > 0)]
}

# Create a table of hexagons crossed during migration and the abundance of birds crossing it
df1 <- data.frame()
for(k in 1:length(EMD_flows2)){
	EMD_flows3 <- list()
	for(i in 1:length(EMD_flows2[[k]])){
		EMD_flows3[[i]] <- unlist(rep(EMD_flows2[[k]][i], length(leastCostPaths[[k]][[i]])))
	}
	df2 <- data.frame()
	for(j in 1:length(unique(unlist(leastCostPaths[[k]])))){
		df2[j,1] <- unique(unlist(leastCostPaths[[k]]))[j]
		df2[j,2] <- sum(unlist(EMD_flows3)[which(unlist(leastCostPaths[[k]]) == unique(unlist(leastCostPaths[[k]]))[j])])
	}
	df1 <- rbind(df1, df2)
}
df <- data.frame()
for(j in 1:length(unique(df1[,1]))){
	df[j,1] <- unique(df1[,1])[j]
	df[j,2] <- sum(df1[,2][which(df1[,1] == df[j,1])])
}


# Plot migration trajectories for single breeding site
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
k = 350
for(i in 1:length(leastCostPaths[[k]])){
	points(hexgrid2_stem_centroids[leastCostPaths[[k]][[i]],1], hexgrid2_stem_centroids[leastCostPaths[[k]][[i]],2], type='l', col="black")	
	#plot(hexgrid2_stem[leastCostPaths[[k]][[i]],], col="yellow", border="yellow", add=T)
	plot(hexgrid2_stem[leastCostPaths[[k]][[i]][1],], col="red", border="red", add=T)
	plot(hexgrid2_stem[leastCostPaths[[k]][[i]][length(leastCostPaths[[k]][[i]])],], col="blue", border="blue", add=T)
}

# Plot migration flow across entire species
pdf(paste("Redistribution-model/YellowWarbler/outputs/MigrationFlows_", alpha, "_", beta, ".pdf", sep=""), width=4, height=4)
par(mfrow=c(1,1), mar=c(0.5,0.5,0.5,0.5), mgp=c(2,1,0), bg="white")
rbPal <- colorRampPalette(c("yellow2","red3"))
datcol <- rbPal(6)[cut(df[,2], breaks=c(0, 0.01, 0.2, 0.5, 0.10, 0.25, 1))]
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid2_stem[df[,1]], col=datcol, border=datcol, add=T)
dev.off()









# Plot connectivity by regions (for EMD)  –– Tree Swallow

# South west
hexsBR <- which(hexgrid3_stem_centroids[,1][which(summer.abundance>0)] <= -120 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] <= 55 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] > 25)
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="dark blue", border="dark blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][hexsBR], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[hexsBR,], 2, sum) > 0),], col="orange", border="orange", add=T)

# North west
hexsBR <- which(hexgrid3_stem_centroids[,1][which(summer.abundance>0)] <= -110 & hexgrid3_stem_centroids[,1][which(summer.abundance>0)] > -180 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] > 60 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] <= 80)
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="dark blue", border="dark blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][hexsBR], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[hexsBR,], 2, sum) > 0),], col="orange", border="orange", add=T)

# Central
hexsBR <- which(hexgrid3_stem_centroids[,1][which(summer.abundance>0)] <= -90 & hexgrid3_stem_centroids[,1][which(summer.abundance>0)] > -105 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] <= 65 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] > 45)
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="dark blue", border="dark blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][hexsBR], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[hexsBR,], 2, sum) > 0),], col="orange", border="orange", add=T)

# East
hexsBR <- which(hexgrid3_stem_centroids[,1][which(summer.abundance>0)] > -80 & hexgrid3_stem_centroids[,2][which(summer.abundance>0)] > 25)
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="dark blue", border="dark blue", add=T)
plot(hexgrid3_stem[which(summer.abundance>0),][hexsBR], col="red", border="red", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),][which(apply(EMD_flows[hexsBR,], 2, sum) > 0),], col="orange", border="orange", add=T)


##  GENOSCAPES  ##

genoscapes <- readOGR("Redistribution-model/Genoscapes/","BirdGenoscape_SppData", verbose=FALSE)
map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(genoscapes[38,], add=T)


which(as.character(genoscapes@data$SCINAME) == "Setophaga petechia")



## EMD in R ##

summer.abundance[which(summer.abundance == "NaN")] <- 0
winter.abundance[which(winter.abundance == "NaN")] <- 0
BR = cbind(summer.abundance, hexgrid3_stem_centroids)
NB = cbind(winter.abundance, hexgrid3_stem_centroids)
distanceMatrix <- rdist.earth(hexgrid3_stem_centroids, hexgrid3_stem_centroids, miles=F)
emd.results <- emdr(BR, NB, flows=T)

map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(0,70))
plot(hexgrid3_stem[which(summer.abundance>0),], col="dark blue", border="dark blue", add=T)
plot(hexgrid3_stem[which(winter.abundance>0),], col="light blue", border="light blue", add=T)

plot(hexgrid3_stem[unique(attr(emd.results, "flows")[[1]]),], col="red", border="red", add=T)
plot(hexgrid3_stem[unique(attr(emd.results, "flows")[[2]]),], col="orange", border="orange", add=T)

for(i in 1:length(attr(emd.results, "flows")[[1]])){
	spl = SpatialLines(list(Lines(Line(rbind(hexgrid3_stem_centroids[attr(emd.results, "flows")[[1]][i],], hexgrid3_stem_centroids[attr(emd.results, "flows")[[2]][i],])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	plot(spl, add=T, col="red", lwd=0.1)
}































plot(hexgrid3_stem_centroids[,1], hexgrid3_stem_centroids[,2], pch=20, cex=1.3, ylim=c(0,90), xlim=c(-150,-40), col="grey")
points(hexgrid3_stem_centroids[,1][which(summer.abundance>0)], hexgrid3_stem_centroids[,2][which(summer.abundance>0)], pch=20, cex=1.3, ylim=c(0,90), xlim=c(-150,-40), col="dark blue")
points(hexgrid3_stem_centroids[,1][which(winter.abundance>0)], hexgrid3_stem_centroids[,2][which(winter.abundance>0)], pch=20, cex=1.3, ylim=c(0,90), xlim=c(-150,-40), col="light blue")
points(hexgrid3_stem_centroids[,1][which(summer.abundance>0)][k], hexgrid3_stem_centroids[,2][which(summer.abundance>0)][k], pch=20, cex=2, ylim=c(0,90), xlim=c(-150,-40), col="red")
points(hexgrid3_stem_centroids[,1][which(winter.abundance>0)][hexs.NB[which(hexs.BR == k)]], hexgrid3_stem_centroids[,2][which(winter.abundance>0)][hexs.NB[which(hexs.BR == k)]], pch=20, cex=2, ylim=c(0,90), xlim=c(-150,-40), col="orange")
k=k+10




# Compute average intensity and bearing in each hexagon

bearingMatrix <- read.csv("Redistribution-model/bearingMatrix_stem.csv", sep=" ", header=F)
bearingMatrix <- bearingMatrix[which(summer.abundance > 0),][, which(winter.abundance > 0)]

hex.migrationDirection <- list()
for(i in 1:length(hexs.BR)){
	spl = SpatialLines(list(Lines(Line(rbind(hexgrid3_stem_centroids[which(summer.abundance > 0),][hexs.BR[i],], hexgrid3_stem_centroids[which(winter.abundance > 0),][hexs.NB[i],])),ID="a")), proj4string = CRS(proj4string(hexgrid3_stem)))
	#plot(spl, add=T, col="red")
	intersec <- gIntersects(spl, hexgrid3_stem, byid=T)
	hex.migrationDirection[which(intersec==T)] = lapply(hex.migrationDirection[which(intersec==T)], function(x) append(x,bearingMatrix[hexs.BR[i], hexs.NB[i]]))
}

# Migration intensity
hex.migrationDirection.intensity <- unlist(lapply(hex.migrationDirection, length))
hex.migrationDirection.intensity[(length(hex.migrationDirection.intensity)+1):length(hexgrid3_stem_centroids[,1])] <- 0
rbPal <- colorRampPalette(c("yellow2","red3"))
datcolRes <- rbPal(6)[as.numeric(cut(hex.migrationDirection.intensity, breaks=seq(0,600,100)))] #seq(0,500,100) seq(0,375,75) c(0,150,300,450,600,760)
datcolRes[which(hex.migrationDirection.intensity == 0)] <- "grey"
plot(hexgrid3_stem_centroids[,1], hexgrid3_stem_centroids[,2], pch=20, cex=1.3, ylim=c(0,90), xlim=c(-150,-40), col=datcolRes)

# Migration direction
get.u <- function(l,a){
	return(l*sin((pi/180)*a))
}
get.v <- function(l,a){
	return(l*cos((pi/180)*a))
}
hex.migrationDirection.u <- list()
hex.migrationDirection.v <- list()
for(i in 1:length(hex.migrationDirection)){
	hex.migrationDirection.u[[i]] <- get.u(hex.migrationDirection.intensity[i], hex.migrationDirection[[i]])
	hex.migrationDirection.v[[i]] <- get.v(hex.migrationDirection.intensity[i], hex.migrationDirection[[i]])
}
hex.migrationDirection.u.mean <- unlist(lapply(hex.migrationDirection.u, function(x){mean(x,na.rm=T)}))
hex.migrationDirection.v.mean <- unlist(lapply(hex.migrationDirection.v, function(x){mean(x,na.rm=T)}))
hex.migrationDirection.u.mean[(length(hex.migrationDirection.u.mean)+1):length(hexgrid3_stem_centroids[,1])] <- NA
hex.migrationDirection.v.mean[(length(hex.migrationDirection.v.mean)+1):length(hexgrid3_stem_centroids[,1])] <- NA


# Spring migration

map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-170, -30), ylim=c(10,70))
plot(isea3h7_hexgrid_WH[Z.albicollis_BR,], col="orange", border="orange", add=T)
plot(isea3h7_hexgrid_WH[Z.albicollis_NB,], col="light blue", border="light blue", add=T)
for(i in 1:length(isea3h7_hexgrid_WH)){
	if(hex.migrationDirection.intensity > 0){arrow.plot(envData$LONGITUDE[i], envData$LATITUDE[i], u=hex.migrationDirection.u.mean[i], v=hex.migrationDirection.v.mean[i], arrow.ex= hex.migrationDirection.intensity[i]/10000, lwd=2, length=0.02)}
}


# Fall migration

map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-170, -30), ylim=c(10,70))
plot(hexgrid3_stem[which(summer.abundance > 0),], col="orange", border="orange", add=T)
plot(hexgrid3_stem[which(winter.abundance > 0),], col="light blue", border="light blue", add=T)
for(i in 1:length(hexgrid3_stem_centroids)){
	if(hex.migrationDirection.intensity[i] > 0){arrow.plot(hexgrid3_stem_centroids[,1][i], hexgrid3_stem_centroids[,2][i], u=-hex.migrationDirection.u.mean[i], v=-hex.migrationDirection.v.mean[i], arrow.ex= hex.migrationDirection.intensity[i]/10000, lwd=2, length=0.02)}
}












##  Compute pairwise distance and bearing along rhumb lines  ##
#distanceMatrix <- matrix(nrow=length(hexgrid3_centroids[,1]), ncol=length(hexgrid3_centroids[,1]))
#bearingMatrix <- matrix(nrow=length(hexgrid3_centroids[,1]), ncol=length(hexgrid3_centroids[,1]))
#for(i in 1:length(hexgrid3_centroids[,1])){	
#	distanceMatrix[i,] <- apply(hexgrid3_centroids, 1, function(x) (distRhumb(x, hexgrid3_centroids[i,])/1000))	
#	bearingMatrix[i,] <- apply(hexgrid3_centroids, 1, function(x) bearingRhumb(x, hexgrid3_centroids[i,]))
#}
#write.table(distanceMatrix, "Redistribution-model/distanceMatrix_stem.csv", row.names=F, col.names=F)
#write.table(bearingMatrix, "Redistribution-model/bearingMatrix_stem.csv", row.names=F, col.names=F)






# Compute average migration intensity and bearing in each hexagon

# Parameters
beta = 0.05 #0.05 	# Parameter associated with the level of stochasticity (smaller value = more stochasticity)
alpha1 = 0.1 	# Parameter associated with cost of elevation
alpha2 = 160		# Parameter associated with cost of crossing water

migration.steps <- matrix(0,nrow=length(hexas.winter), ncol=length(hexas.summer))
neigList <- apply(transitionMat, 1, function(x) which(x>0))

elevationAmericas <- elevationAmericas.mat[1,]
elevationAmericasList <- lapply(neigList, function(x) elevationAmericas[x])
water <- water.mat[1,]
waterList <- lapply(neigList, function(x) water[x])

hexgrid3_centroids <- hexgrid2_centroids[which(water == 0),]

migration.tracks <- list()
for(k in 1:length(hexs.BR)){
		
	destination_hexagon <- hexas.winter[hexs.NB[k]] #7593 #11906 # land[5300]
	bearingsToDestinationHex = apply(hexgrid2_centroids, 1, function(x) bearingRhumb(x, hexgrid3_centroids[destination_hexagon,]))
	
	# Transition matrix with transition probabilities based on attraction to destination hexagon, elevation and water costs
	transitionList <- list()
	for(i in 1:length(hexgrid2_centroids[,1])){
		bearingsNeighbours <- bearingRhumb(hexgrid2_centroids[i,], hexgrid2_centroids[neigList[[i]],])
		bearingsNeighbours <- abs(bearingsNeighbours - bearingsToDestinationHex[i])
		bearingsNeighbours <- sapply(bearingsNeighbours, function(x) ifelse(x>180, 360-x, x))
		transitionList[[i]] <- exp(-beta * (bearingsNeighbours + (alpha1*elevationAmericasList[[i]]) + (alpha2*waterList[[i]])))
	}
		
	starting.hex <- hexas.summer[hexs.BR[k]] #4620 #1898 #2302 #2726 # land[1500]
	hex.population <- rep(0,length(hexgrid2_centroids[,1]))
	hex.population[which(water == 0)][starting.hex] = 10
	occupied.hexs.list <- list()
	t=2
	occupied.hexs.list[[1]] <- matrix(10)
	row.names(occupied.hexs.list[[1]]) <- which(water == 0)[starting.hex]
	while(hex.population[which(water == 0)][destination_hexagon] < 5 & t < 1000){
		occupied.hexs <- which(hex.population>0)
		occupied.hexs.prelist <- list()
		hex.population2 <- rep(0,length(hexgrid2_centroids[,1]))
		for(i in 1:length(occupied.hexs)){
			
			if(occupied.hexs[i] != which(water == 0)[destination_hexagon]){
				#hex.neighbours <- which(transitionMat2[occupied.hexs[i],] > 0)
				hex.neighbours <- neigList[[occupied.hexs[i]]]
				#hex.neighbours.selected <- sample(hex.neighbours, hex.population[occupied.hexs[i]], replace=T, prob= transitionMat2[occupied.hexs[i],][hex.neighbours]/sum(transitionMat2[occupied.hexs[i],][hex.neighbours]))
				hex.neighbours.selected <- sample(hex.neighbours, hex.population[occupied.hexs[i]], replace=T, prob= transitionList[[occupied.hexs[i]]]/sum(transitionList[[occupied.hexs[i]]]))
				for(j in 1:length(as.matrix(table(hex.neighbours.selected))[,1])){
					hex.population2[as.numeric(row.names(as.matrix(table(hex.neighbours.selected)))[j])] = hex.population2[as.numeric(row.names(as.matrix(table(hex.neighbours.selected)))[j])] + as.matrix(table(hex.neighbours.selected))[j,1]
				}
				occupied.hexs.prelist[[i]] <- as.matrix(table(hex.neighbours.selected))
			}else{
				hex.population2[occupied.hexs[i]] <- hex.population2[occupied.hexs[i]] + hex.population[occupied.hexs[i]]
				occupied.hexs.prelist[[i]] <- matrix(hex.population[occupied.hexs[i]])
				row.names(occupied.hexs.prelist[[i]]) <- occupied.hexs[i]
			}		
		}
		occupied.hexs.list[[t]] <- occupied.hexs.prelist[[1]]
		if(i>1){
			for(h in 2:i){
				occupied.hexs.list[[t]] <- rbind(occupied.hexs.list[[t]], occupied.hexs.prelist[[h]])		
			}
		}
		hex.population <- hex.population2
		#hex.population[as.numeric(row.names(occupied.hexs.list[[t-1]]))] <- 0
		t=t+1
	}
	migration.tracks[[k]] <- occupied.hexs.list
	#migration.steps[hw,hs] <- t
	#print(hs)
}


















## OLD STUFF


# Summer -- 07 June to 03 Aug   ##02 May to 21 August  # YelWar: 14 June – 13 July

weekly.abundance.summer <- list()
weekly.abundance.summer.hex <- list()
for(i in 23:30){ #18:34
	#weekly.abundance.summer[[i]] <- raster(paste("stem/grycat/abundance_umean/", list.files("stem/grycat/abundance_umean/")[i], sep=""))	
	weekly.abundance.summer[[i]] <- abunds[[i]]
	weekly.abundance.summer[[i]] <- projectRaster(weekly.abundance.summer[[i]], crs=sr)
	weekly.abundance.summer.hex[[i]] <- extract(weekly.abundance.summer[[i]], hexgrid3_stem, fun=mean, na.rm=T)
}
#summer.abundance <- (weekly.abundance.summer.hex[[18]] + weekly.abundance.summer.hex[[19]] + weekly.abundance.summer.hex[[20]] + weekly.abundance.summer.hex[[21]] + weekly.abundance.summer.hex[[22]] + weekly.abundance.summer.hex[[23]] + weekly.abundance.summer.hex[[24]] + weekly.abundance.summer.hex[[25]] + weekly.abundance.summer.hex[[26]] + weekly.abundance.summer.hex[[27]] + weekly.abundance.summer.hex[[28]] + weekly.abundance.summer.hex[[29]] + weekly.abundance.summer.hex[[30]] + weekly.abundance.summer.hex[[31]] + weekly.abundance.summer.hex[[32]] + weekly.abundance.summer.hex[[33]] + weekly.abundance.summer.hex[[34]]) / 17
summer.abundance <- (weekly.abundance.summer.hex[[23]] + weekly.abundance.summer.hex[[24]] + weekly.abundance.summer.hex[[25]] + weekly.abundance.summer.hex[[26]] + weekly.abundance.summer.hex[[27]] + weekly.abundance.summer.hex[[28]] + weekly.abundance.summer.hex[[29]] + weekly.abundance.summer.hex[[30]]) / 8


# Winter -- 14 Novembre to 28 February # YelWar: 9 Nov – 1 Mar

weekly.abundance.winter <- list()
weekly.abundance.winter.hex <- list()
for(i in c(1:8,46:52)){
	#weekly.abundance.winter[[i]] <- raster(paste("stem/grycat/abundance_umean/", list.files("stem/grycat/abundance_umean/")[i], sep=""))
	weekly.abundance.winter[[i]] <- abunds[[i]]
	weekly.abundance.winter[[i]] <- projectRaster(weekly.abundance.winter[[i]], crs=sr)
	weekly.abundance.winter.hex[[i]] <- extract(weekly.abundance.winter[[i]], hexgrid3_stem, fun=mean, na.rm=T)
}
winter.abundance <- (weekly.abundance.winter.hex[[1]] + weekly.abundance.winter.hex[[2]] + weekly.abundance.winter.hex[[3]] + weekly.abundance.winter.hex[[4]] + weekly.abundance.winter.hex[[5]] + weekly.abundance.winter.hex[[6]] + weekly.abundance.winter.hex[[7]] + weekly.abundance.winter.hex[[8]] + weekly.abundance.winter.hex[[46]] + weekly.abundance.winter.hex[[47]] + weekly.abundance.winter.hex[[48]] + weekly.abundance.winter.hex[[49]] + weekly.abundance.winter.hex[[50]] + weekly.abundance.winter.hex[[51]] + weekly.abundance.winter.hex[[52]]) / 15

seasonal.abundance <- cbind(summer.abundance, winter.abundance)


plot(NULL,xlim=c(-170,-30), ylim=c(-60,80))
plot(weekly.abundance.winter[[1]], add=T)
plot(newmap, add=T)


rbPal <- colorRampPalette(c("yellow2","red3"))
datcol <- rbPal(5)[as.numeric(cut(summer.abundance, breaks=seq(range(summer.abundance, na.rm=T)[1], range(summer.abundance, na.rm=T)[2], range(summer.abundance, na.rm=T)[2]/5)))]
datcol[which(summer.abundance == "NaN" | is.na(summer.abundance) == T | summer.abundance == 0)] <- "grey"
plot(hexgrid3_stem, col=datcol, border=datcol)

points(hexgrid3_centroids[which(summer.abundance>0),][232,1], hexgrid3_centroids[which(summer.abundance>0),][232,2], pch=20, col="blue")
points(hexgrid3_centroids[which(winter.abundance>0),][168,1], hexgrid3_centroids[which(winter.abundance>0),][168,2], pch=20, col="blue")


#write.table(seasonal.abundance, "seasonalAbundance_whtspa.csv", row.names=F, col.names=F)
seasonal.abundance <- read.csv("Redistribution-model/YellowWarbler/seasonalAbundance_yelwar.csv", sep=" ", header=F)
summer.abundance = seasonal.abundance[,1] 
winter.abundance = seasonal.abundance[,2] 


plot(hexgrid3_stem, col="white", border="white", bg='grey', xlim=c(-150,-40), ylim=c(15,70))
plot(hexgrid3_stem[which(winter.abundance>0.01)], add=T, col="light blue", border="light blue")
plot(hexgrid3_stem[which(summer.abundance>0.01)], add=T, col="orange2", border="orange2")


##  Compute pairwise distances between seasonally occupied hexagons

load("Redistribution-model/data_transitionMatrix.RData")


hexas.summer <- vector()
for(i in 1:length(which(summer.abundance > 0))){
	dd = as.matrix(dist(rbind(hexgrid3_stem_centroids[which(summer.abundance > 0.01),][i,], hexgrid2_stem_centroids[land,])))[1,(2:(length(land)+1))]
	hexas.summer[i] <- which.min(dd)	
	print(i)
}
hexas.winter <- vector()
for(i in 1:length(which(winter.abundance > 0))){
	dd = as.matrix(dist(rbind(hexgrid3_stem_centroids[which(winter.abundance > 0.01),][i,], hexgrid2_centroids[land,])))[1,(2:(length(land)+1))]
	hexas.winter[i] <- which.min(dd)	
	print(i)
}
write.csv(hexas.summer, "Redistribution-model/hexasSummer.csv")
write.csv(hexas.winter, "Redistribution-model/hexasWinter.csv")


which(summer.abundance > 0)
which(winter.abundance > 0)


hexas.summer <- read.csv("Redistribution-model/hexasSummer.csv", header=F)[,1]
hexas.winter <- read.csv("Redistribution-model/hexasWinter.csv", header=F)[,1]

















k = 225
a <- as.data.frame(table(hexs.NB[which(hexs.BR == k)]))
for(i in 1:length(a[,1])){
	if(a[i,2]<=5){
		points(hexgrid3_centroids[,1][which(winter.abundance != 0)[as.numeric(as.vector(a[,1]))[i]]], hexgrid3_centroids[,2][which(winter.abundance != 0)[as.numeric(as.vector(a[,1]))[i]]], pch=20, cex=1.3, ylim=c(0,90), xlim=c(-150,-40), col="orange")
	}else if(a[i,2]<=10){
		points(hexgrid3_centroids[,1][which(winter.abundance != 0)[as.numeric(as.vector(a[,1]))[i]]], hexgrid3_centroids[,2][which(winter.abundance != 0)[as.numeric(as.vector(a[,1]))[i]]], pch=20, cex=1.3, ylim=c(0,90), xlim=c(-150,-40), col="red")
	}else{
		points(hexgrid3_centroids[,1][which(winter.abundance != 0)[as.numeric(as.vector(a[,1]))[i]]], hexgrid3_centroids[,2][which(winter.abundance != 0)[as.numeric(as.vector(a[,1]))[i]]], pch=20, cex=1.3, ylim=c(0,90), xlim=c(-150,-40), col="brown4")
	}
}
points(hexgrid3_centroids[,1][which(summer.abundance != 0)[k]], hexgrid3_centroids[,2][which(summer.abundance != 0)[k]], pch=20, cex=1.3, ylim=c(0,90), xlim=c(-150,-40), col="green")

for(i in 1:length(hexs.BR)){
	d = rdist.earth(t(as.matrix(hexgrid3_centroids[which(summer.abundance != 0),][hexs.BR[i],])), t(as.matrix(hexgrid3_centroids[which(winter.abundance != 0),][hexs.NB[i],])))
	inter <- gcIntermediate(hexgrid3_centroids[which(summer.abundance != 0),][hexs.BR[i],], hexgrid3_centroids[which(winter.abundance != 0),][hexs.NB[i],], n=50, addStartEnd=T)
	if(d<1000){
		lines(inter, lwd=1, col="yellow")
	}else if(d<2000){
		lines(inter, lwd=1, col="orange")
	}else if(d<3000){
		lines(inter, lwd=1, col="red")
	}else{
		lines(inter, lwd=1, col="brown4")
	}
}

map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-170, -30), ylim=c(10,70))
i = 680
br.hexas = c(i, which(gTouches(hexgrid3[which(summer.abundance != 0)][i], hexgrid3[which(summer.abundance != 0)], byid=T) == TRUE))
plot(hexgrid3[which(summer.abundance != 0),][br.hexas], col="dark green", border="dark green", add=T)
a = vector()
for(k in 1:length(br.hexas)){
	a <- c(a,which(hexs.BR == br.hexas[k]))
}
aa <- as.data.frame(table(hexs.NB[a]))
rbPal <- colorRampPalette(c("yellow2","red3"))
datcol <- rbPal(5)[as.numeric(cut(aa[,2], breaks=seq(0, range(aa[,2], na.rm=T)[2]+1, (range(aa[,2], na.rm=T)[2]+1)/5)))]
plot(hexgrid3[which(winter.abundance != 0)][as.numeric(as.vector(aa[,1]))], col=datcol, border=datcol, add=T)





for(i in 1:length(aa[,1])){
	if(aa[i,2]<=5){
		plot(hexgrid3[which(winter.abundance != 0),][as.numeric(as.vector(aa[,1]))[i]]], col="yellow", border="yellow", add=T)
	}else if(a[i,2]<=10){
		plot(hexgrid3[which(winter.abundance != 0),][as.numeric(as.vector(aa[,1]))[i]]], col="orange", border="orange", add=T)
	}else{
		plot(hexgrid3[which(winter.abundance != 0),][as.numeric(as.vector(aa[,1]))[i]]], col="red", border="red", add=T)
	}
}





# Compute average intensity and bearing in each hexagon

bearingMatrix <- read.csv("Redistribution-model/bearingMatrix_stem.csv", sep=" ", header=F)
bearingMatrix <- bearingMatrix[which(summer.abundance > 0),][, which(winter.abundance > 0)]

hex.migrationDirection <- list()
for(i in 1:length(hexs.BR)){
	spl = SpatialLines(list(Lines(Line(rbind(hexgrid3_centroids[which(summer.abundance != 0),][hexs.BR[i],], hexgrid3_centroids[which(winter.abundance != 0),][hexs.NB[i],])),ID="a")), proj4string = CRS(proj4string(hexgrid3)))
	#plot(spl, add=T, col="red")
	intersec <- gIntersects(spl, hexgrid3, byid=T)
	hex.migrationDirection[which(intersec==T)] = lapply(hex.migrationDirection[which(intersec==T)], function(x) append(x,bearingMatrix[hexs.BR[i], hexs.NB[i]]))
}

# Migration intensity
hex.migrationDirection.intensity <- unlist(lapply(hex.migrationDirection, length))
hex.migrationDirection.intensity[(length(hex.migrationDirection.intensity)+1):length(hexgrid3_centroids[,1])] <- 0
rbPal <- colorRampPalette(c("yellow2","red3"))
datcolRes <- rbPal(6)[as.numeric(cut(hex.migrationDirection.intensity, breaks=seq(0,600,100)))] #seq(0,500,100) seq(0,375,75) c(0,150,300,450,600,760)
datcolRes[which(hex.migrationDirection.intensity == 0)] <- "grey"
plot(hexgrid3_centroids[,1], hexgrid3_centroids[,2], pch=20, cex=1.3, ylim=c(0,90), xlim=c(-150,-40), col=datcolRes)

# Migration direction
get.u <- function(l,a){
	return(l*sin((pi/180)*a))
}
get.v <- function(l,a){
	return(l*cos((pi/180)*a))
}
hex.migrationDirection.u <- list()
hex.migrationDirection.v <- list()
for(i in 1:length(hex.migrationDirection)){
	hex.migrationDirection.u[[i]] <- get.u(hex.migrationDirection.intensity[i], hex.migrationDirection[[i]])
	hex.migrationDirection.v[[i]] <- get.v(hex.migrationDirection.intensity[i], hex.migrationDirection[[i]])
}
hex.migrationDirection.u.mean <- unlist(lapply(hex.migrationDirection.u, function(x){mean(x,na.rm=T)}))
hex.migrationDirection.v.mean <- unlist(lapply(hex.migrationDirection.v, function(x){mean(x,na.rm=T)}))
hex.migrationDirection.u.mean[(length(hex.migrationDirection.u.mean)+1):length(hexgrid3_centroids[,1])] <- NA
hex.migrationDirection.v.mean[(length(hex.migrationDirection.v.mean)+1):length(hexgrid3_centroids[,1])] <- NA


# Spring migration

map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-170, -30), ylim=c(10,70))
plot(isea3h7_hexgrid_WH[Z.albicollis_BR,], col="orange", border="orange", add=T)
plot(isea3h7_hexgrid_WH[Z.albicollis_NB,], col="light blue", border="light blue", add=T)
for(i in 1:length(isea3h7_hexgrid_WH)){
	if(hex.migrationDirection.intensity > 0){arrow.plot(envData$LONGITUDE[i], envData$LATITUDE[i], u=hex.migrationDirection.u.mean[i], v=hex.migrationDirection.v.mean[i], arrow.ex= hex.migrationDirection.intensity[i]/10000, lwd=2, length=0.02)}
}


# Fall migration

map("world", fill=T, col="white", border="white", bg="dark grey", mar=c(0,0,0,0), xlim=c(-170, -30), ylim=c(10,70))
plot(hexgrid3[which(summer.abundance != 0),], col="orange", border="orange", add=T)
plot(hexgrid3[which(winter.abundance != 0),], col="light blue", border="light blue", add=T)
for(i in 1:length(hexgrid3_centroids)){
	if(hex.migrationDirection.intensity[i] > 0){arrow.plot(hexgrid3_centroids[,1][i], hexgrid3_centroids[,2][i], u=-hex.migrationDirection.u.mean[i], v=-hex.migrationDirection.v.mean[i], arrow.ex= hex.migrationDirection.intensity[i]/10000, lwd=2, length=0.02)}
}











## BANDING DATA

encounter.data <- read.csv("/Volumes/Marius_SSD/Banding_data/Encounter_files/Somveille_all_Rs_of_3_spp_201810110843.csv")
encounter.data.grycat <- encounter.data[which(encounter.data$E_SPECIES_NAME == "Gray Catbird"),]
#encounter.data.yelwar <- encounter.data[which(encounter.data$B_SPECIES_NAME == "Yellow Warbler"),]

encounter.data.grycat <- encounter.data.grycat[which(((encounter.data.grycat$BANDING_MONTH >=4 & encounter.data.grycat$BANDING_MONTH <= 9) & (encounter.data.grycat$ENCOUNTER_MONTH >= 10 | encounter.data.grycat$ENCOUNTER_MONTH <= 3)) | ((encounter.data.grycat$ENCOUNTER_MONTH >=4 & encounter.data.grycat$ENCOUNTER_MONTH <= 9) & (encounter.data.grycat$BANDING_MONTH >= 10 | encounter.data.grycat$BANDING_MONTH <= 3))),]

encounter.data.herthr <- encounter.data[which(encounter.data$E_SPECIES_NAME == "Hermit Thrush"),]

Hermit Thrush


encounter.data.yelwar <- encounter.data.yelwar[which(((encounter.data.yelwar$BANDING_MONTH >=4 & encounter.data.yelwar$BANDING_MONTH <= 9) & (encounter.data.yelwar$ENCOUNTER_MONTH >= 10 | encounter.data.yelwar$ENCOUNTER_MONTH <= 3) & encounter.data.yelwar$ENCOUNTER_MONTH <= 12) | ((encounter.data.yelwar$ENCOUNTER_MONTH >=4 & encounter.data.yelwar$ENCOUNTER_MONTH <= 9) & (encounter.data.yelwar$BANDING_MONTH >= 10 | encounter.data.yelwar$BANDING_MONTH <= 3) & encounter.data.yelwar$BANDING_MONTH <= 12)),]

cbind(encounter.data.yelwar$B_LAT_DECIMAL_DEGREES, encounter.data.yelwar$BANDING_MONTH, encounter.data.yelwar$E_LAT_DECIMAL_DEGREES, encounter.data.yelwar$ENCOUNTER_MONTH)

for(i in 1:length(encounter.data.yelwar$B_LON_DECIMAL_DEGREES)){
	inter <- gcIntermediate(cbind(encounter.data.yelwar$B_LON_DECIMAL_DEGREES[i], encounter.data.yelwar$B_LAT_DECIMAL_DEGREES[i]), cbind(encounter.data.yelwar$E_LON_DECIMAL_DEGREES[i], encounter.data.yelwar$E_LAT_DECIMAL_DEGREES[i]), n=50, addStartEnd=T)
	lines(inter, lwd=2, col="yellow")
	points(encounter.data.yelwar$B_LON_DECIMAL_DEGREES[i], encounter.data.yelwar$B_LAT_DECIMAL_DEGREES[i], pch=20, col="orange")
	points(encounter.data.yelwar$E_LON_DECIMAL_DEGREES[i], encounter.data.yelwar$E_LAT_DECIMAL_DEGREES[i], pch=20, col="orange")
}




