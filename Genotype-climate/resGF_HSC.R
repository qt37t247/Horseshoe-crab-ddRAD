## Load library
library(vcfR)
library(resGF)
library(raster)

## Obtain some environmental data
files <- list.files(pattern = ".asc$", recursive = F, full.names = T)
files_5026 <- list.files(path = "./2050_RCP26", pattern = ".asc$", recursive = F, full.names = T)
files_0026 <- list.files(path = "./2100_RCP26", pattern = ".asc$", recursive = F, full.names = T)
files_5085 <- list.files(path = "./2050_RCP85", pattern = ".asc$", recursive = F, full.names = T)
files_0085 <- list.files(path = "./2100_RCP85", pattern = ".asc$", recursive = F, full.names = T)

nows <- do.call("stack", lapply(files, function(fn){raster(fn)}))
e5026s <- do.call("stack", lapply(files_5026, function(fn){raster(fn)}))
e0026s <- do.call("stack", lapply(files_0026, function(fn){raster(fn)}))
e5085s <- do.call("stack", lapply(files_5085, function(fn){raster(fn)}))
e0085s <- do.call("stack", lapply(files_0085, function(fn){raster(fn)}))

library(RColorBrewer)
pal <- colorRampPalette(c("blue","red"))

###########################CR###############################
## create a genind object
CRv <- read.vcfR("CR0L.vcf")
CR = vcfR2genind(CRv)

## obtain spatial data
CR.coord  <- read.table(file = 'CR.coordinates', sep = '\t', header = T, row.names = 1, fill = T)
CR.coord.sp <- SpatialPointsDataFrame(CR.coord[,c('x','y')], proj4string=crs(nows), data=CR.coord)

library(rSDM)
CR.coord.sp <- points2nearestcell(CR.coord.sp, nows, layer=19)

# Gradient Forest
library(gradientForest)
# extract points from
CR.clim.points <- raster::extract(nows, CR.coord.sp)
CR.clim.points[is.na(CR.clim.points)] <- 0 #change three NA values of bathymetry to 0
#generates the PCNMs
library(vegan)
CR.pcnm <- pcnm(dist(CR.coord.sp@data))
CR.keep <- round(length(which(CR.pcnm$value > 0))/2)
CR.pcnm.keep <- scores(CR.pcnm)[,1:CR.keep]  #keep half of positive ones as suggested by some authors

CR.nbvar <- length(colnames(CR.clim.points))
CR.env.gf <- cbind(CR.clim.points[ , c(seq(1, CR.nbvar, by=1)) ], CR.pcnm.keep)
CR.maxLevel <- log2(0.368*nrow(CR.env.gf)/2)
CR.env.gf <- as.data.frame(CR.env.gf)
CR.snp <- as.data.frame(CR$tab)

# run gradient forest
colnames(CR.snp) <- paste("A",1:length(colnames(CR.snp)), sep = "")
rownames(CR.snp) <- rownames(CR.env.gf)
CR.gf <- gradientForest(cbind(CR.env.gf, CR.snp), predictor.vars=colnames(CR.env.gf), response.vars=colnames(CR.snp), ntree=500, maxLevel=CR.maxLevel, trace=T, corr.threshold=0.50, compact = T, nbin = 201)
# generate resistance surface
CR.single_r <- resGF(CR.gf, nows, save.image = T)
dev.off()
pdf("CR_final.pdf")
plot(CR.single_r, col = pal(15))
dev.off()

CR.single_r <- resGF(CR.gf, e5026s, save.image = F)
dev.off()
pdf("CR_5026e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()

CR.single_r <- resGF(CR.gf, e0026s, save.image = F)
dev.off()
pdf("CR_0026e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()

CR.single_r <- resGF(CR.gf, e5085s, save.image = F)
dev.off()
pdf("CR_5085e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()

CR.single_r <- resGF(CR.gf, e0085s, save.image = F)
dev.off()
pdf("CR_0085e.pdf")
plot(CR.single_r, col = pal(15))
dev.off()

###########################CR (excl. Sabah)###############################
## create a genind object
CRv <- read.vcfR("CR0Lx.vcf")
CR = vcfR2genind(CRv)

## obtain spatial data
CR.coord  <- read.table(file = 'CRx.coordinates', sep = '\t', header = T, row.names = 1, fill = T)
CR.coord.sp <- SpatialPointsDataFrame(CR.coord[,c('x','y')], proj4string=crs(nows), data=CR.coord)

library(rSDM)
CR.coord.sp <- points2nearestcell(CR.coord.sp, nows, layer=19)

# Gradient Forest
library(gradientForest)
# extract points from
CR.clim.points <- raster::extract(nows, CR.coord.sp)
CR.clim.points[is.na(CR.clim.points)] <- 0 #change three NA values of bathymetry to 0
#generates the PCNMs
library(vegan)
CR.pcnm <- pcnm(dist(CR.coord.sp@data))
CR.keep <- round(length(which(CR.pcnm$value > 0))/2)
CR.pcnm.keep <- scores(CR.pcnm)[,1:CR.keep]  #keep half of positive ones as suggested by some authors

CR.nbvar <- length(colnames(CR.clim.points))
CR.env.gf <- cbind(CR.clim.points[ , c(seq(1, CR.nbvar, by=1)) ], CR.pcnm.keep)
CR.maxLevel <- log2(0.368*nrow(CR.env.gf)/2)
CR.env.gf <- as.data.frame(CR.env.gf)
CR.snp <- as.data.frame(CR$tab)

# run gradient forest
colnames(CR.snp) <- paste("A",1:length(colnames(CR.snp)), sep = "")
rownames(CR.snp) <- rownames(CR.env.gf)
CR.gf <- gradientForest(cbind(CR.env.gf, CR.snp), predictor.vars=colnames(CR.env.gf), response.vars=colnames(CR.snp), ntree=500, maxLevel=CR.maxLevel, trace=T, corr.threshold=0.50, compact = T, nbin = 201)
# generate resistance surface
CR.single_r <- resGF(CR.gf, nows, save.image = T)
dev.off()
plot(CR.single_r, col = pal(15))

CR.single_r <- resGF(CR.gf, e5026s, save.image = F)
dev.off()
plot(CR.single_r, col = pal(15))

CR.single_r <- resGF(CR.gf, e0026s, save.image = F)
dev.off()
pdf("CRx_0026.pdf")
plot(CR.single_r, col = pal(15))
dev.off()

CR.single_r <- resGF(CR.gf, e5085s, save.image = F)
dev.off()
pdf("CRx_5085.pdf")
plot(CR.single_r, col = pal(15))
dev.off()


CR.single_r <- resGF(CR.gf, e0085s, save.image = F)
dev.off()
pdf("CRx_0085.pdf")
plot(CR.single_r, col = pal(15))
dev.off()

#####################################TG##############################
## create a genind object
TGv <- read.vcfR("TG0L.vcf")
TG = vcfR2genind(TGv)

## obtain spatial data
TG.coord  <- read.table(file = 'TG.coordinates', sep = '\t', header = T, row.names = 1, fill = T)
TG.coord.sp <- SpatialPointsDataFrame(TG.coord[,c('x','y')], proj4string=crs(nows), data=TG.coord)

library(rSDM)
TG.coord.sp <- points2nearestcell(TG.coord.sp, nows, layer=19)

# Gradient Forest
library(gradientForest)
# extract points from
TG.clim.points <- raster::extract(nows, TG.coord.sp)
TG.clim.points[is.na(TG.clim.points)] <- 0 #change four NA values of bathymetry to 0
#generates the PCNMs
library(vegan)
TG.pcnm <- pcnm(dist(TG.coord.sp@data))
TG.keep <- round(length(which(TG.pcnm$value > 0))/2)
TG.pcnm.keep <- scores(TG.pcnm)[,1:TG.keep]  #keep half of positive ones as suggested by some authors

TG.nbvar <- length(colnames(TG.clim.points))
TG.env.gf <- cbind(TG.clim.points[ , c(seq(1, TG.nbvar, by=1)) ], TG.pcnm.keep)
TG.maxLevel <- log2(0.368*nrow(TG.env.gf)/2)
TG.env.gf <- as.data.frame(TG.env.gf)
TG.snp <- as.data.frame(TG$tab)

# run gradient forest
colnames(TG.snp) <- paste("A",1:length(colnames(TG.snp)), sep = "")
rownames(TG.snp) <- rownames(TG.env.gf)
TG.gf <- gradientForest(cbind(TG.env.gf, TG.snp), predictor.vars=colnames(TG.env.gf), response.vars=colnames(TG.snp), ntree=500, maxLevel=TG.maxLevel, trace=T, corr.threshold=0.50, compact = T, nbin = 201)
# generate resistance surface
TG.single_r <- resGF(TG.gf, nows, save.image = T)
dev.off()
plot(TG.single_r, col = pal(15))

TG.single_r <- resGF(TG.gf, e5026s, save.image = F)
dev.off()
plot(TG.single_r, col = pal(15))

TG.single_r <- resGF(TG.gf, e0026s, save.image = F)
dev.off()
plot(TG.single_r, col = pal(15))

TG.single_r <- resGF(TG.gf, e5085s, save.image = F)
dev.off()
plot(TG.single_r, col = pal(15))

TG.single_r <- resGF(TG.gf, e0085s, save.image = F)
dev.off()
plot(TG.single_r, col = pal(15))


#####################################TT##############################
## create a genind object
TTv <- read.vcfR("TT0L.vcf")
TT = vcfR2genind(TTv)

## obtain spatial data
TT.coord  <- read.table(file = 'TT.coordinates', sep = '\t', header = T, row.names = 1, fill = T)
TT.coord.sp <- SpatialPointsDataFrame(TT.coord[,c('x','y')], proj4string=crs(nows), data=TT.coord)

library(rSDM)
TT.coord.sp <- points2nearestcell(TT.coord.sp, nows, layer=19)

# Gradient Forest
library(gradientForest)
# extract points from
TT.clim.points <- raster::extract(nows, TT.coord.sp)
#generates the PCNMs
library(vegan)
TT.pcnm <- pcnm(dist(TT.coord.sp@data))
TT.keep <- round(length(which(TT.pcnm$value > 0))/2)
TT.pcnm.keep <- scores(TT.pcnm)[,1:TT.keep]  #keep half of positive ones as suggested by some authors

TT.nbvar <- length(colnames(TT.clim.points))
TT.env.gf <- cbind(TT.clim.points[ , c(seq(1, TT.nbvar, by=1)) ], TT.pcnm.keep)
TT.maxLevel <- log2(0.368*nrow(TT.env.gf)/2)
TT.env.gf <- as.data.frame(TT.env.gf)
TT.snp <- as.data.frame(TT$tab)

# run gradient forest
colnames(TT.snp) <- paste("A",1:length(colnames(TT.snp)), sep = "")
rownames(TT.snp) <- rownames(TT.env.gf)
TT.gf <- gradientForest(cbind(TT.env.gf, TT.snp), predictor.vars=colnames(TT.env.gf), response.vars=colnames(TT.snp), ntree=500, maxLevel=TT.maxLevel, trace=T, corr.threshold=0.50, compact = T, nbin = 201)
# generate resistance surface
TT.single_r <- resGF(TT.gf, nows, save.image = T)
dev.off()
plot(TT.single_r, col = pal(15))

TT.single_r <- resGF(TT.gf, e5026s, save.image = F)
dev.off()
plot(TT.single_r, col = pal(15))

TT.single_r <- resGF(TT.gf, e0026s, save.image = F)
dev.off()
plot(TT.single_r, col = pal(15))

TT.single_r <- resGF(TT.gf, e5085s, save.image = F)
dev.off()
plot(TT.single_r, col = pal(15))

TT.single_r <- resGF(TT.gf, e0085s, save.image = F)
dev.off()
plot(TT.single_r, col = pal(15))
