library(gradientForest)
library(raster)
library(RColorBrewer)

#Load pop allele frequency
tt_snp <- read.csv("TT.alfreq", row.names = 1)
#Load bio-factor
tt_env <- read.csv("TT_env.csv", header = TRUE)
tt_env <- tt_env[, -1]
rownames(tt_env) <- rownames(tt_snp)#Match row names
#GF analysis
nSites <- dim(tt_snp)[1]
nSpecs <- dim(tt_snp)[2]
lev <- floor(log2(nSites*0.368/2))
lev
gf <- gradientForest(cbind(tt_env,tt_snp), predictor.vars=colnames(tt_env), 
                     response.vars=colnames(tt_snp), ntree = 500, maxLevel = lev, 
                     corr.threshold = 0.5, compact = T, nbin = 201)
gf
##Important variables:
##[1] BO21_chlorange_ss    BO21_salinitymean_ss BO21_temprange_bdmin BO21_temprange_ss   
##[5] BO21_tempmean_bdmin 

#Plot cumulative important
most_important <- names(importance(gf))[1:17]
plot(gf, plot.type = "O", legend = F)
plot(gf, plot.type = "C", imp.vars = most_important, show.species = F, common.scale = T)

#Simulate current genetic composition
tt_current_bio <- read.csv("current_bio.csv")
tt_current_bio <- tt_current_bio[,-1][complete.cases(tt_current_bio[,-1]), ]
imp.vars <- c("BO21_chlorange_ss","BO21_salinitymean_ss","BO21_temprange_bdmin","BO21_temprange_ss","BO21_tempmean_bdmin")
Trns_grid <- cbind(tt_current_bio[,c("x","y")], predict(gf,tt_current_bio[,imp.vars]))
PCs <- prcomp(Trns_grid[,imp.vars])
sgn <- sign(PCs$rotation["BO21_chlorange_ss",])
PCs$rotation <- sweep(PCs$rotation,2,sgn,"*")
PCs$x <- sweep(PCs$x,2,sgn,"*")
a1 <- PCs$x[,1]
a2 <- PCs$x[,2]
a3 <- PCs$x[,3]
r <- a1+a2
g <- -a2
b <- a3+a2-a1
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255

## plotting PCA
pdf("tt_env_PCA.pdf")
plot(PCs$x[,c("PC1","PC2")],pch=".", ylim = c(-0.1, 0.25), cex=3, asp=1, col=rgb(r,g,b, max = 255))
l.x <- PCs$rotation[,1]/10
l.y <- PCs$rotation[,2]/10
l.pos <- l.y # Create a vector of y axis coordinates
lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
hi <- which(l.y > 0) # Get variables on the top half
# Replace values in the vector
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")
text(l.x, l.y, labels=row.names(PCs$rotation), cex=0.6, pos=l.pos)
arrows(x0=0, x1=l.x, y0=0, y1=l.y, length=0.05)
dev.off()

## plotting transformed genotype-climate association across the distribution ranges
pdf("tt_env_PCA_sites.pdf")
plot(Trns_grid[,c("x","y")],pch=".", cex=3, asp=1, col=rgb(r,g,b, max = 255))
dev.off()

#Simulate future genetic composition
#Just use 2100-RCP8.5 as an example
fut_0085 <- read.csv("fut_0085.csv")
fut_0085 <- fut_0085[,-1][complete.cases(fut_0085[,-1]), ]
colnames(fut_0085) <- colnames(tt_current_bio)
Trns_grid_0085 <- cbind(fut_0085[,c("x","y")], predict(gf,fut_0085[,imp.vars]))
##Compute the GV value
library(dplyr)
library(tidyr)
euclidean <- function(a, b) sqrt(sum((a - b)^2))#Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_0085)){
  enc_dist <- euclidean(Trns_grid[i,imp.vars], Trns_grid_0085[i,imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GV_0085"
}
## this for loop takes long time
tt_gv_0085 <- cbind(fut_0085[,c("x","y")], df)
write.csv(tt_gv_0085, "tt_genomic_offset_0085.csv")

#Convert points of average value to raster
tt_gv_0085_raster <- rasterFromXYZ(tt_gv_0085)

pdf("tt_offset.pdf")
pal <- colorRampPalette(c("white","red"))
plot(tt_gv_0085_raster, col = pal(15))
dev.off()