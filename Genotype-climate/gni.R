#This script shows how to calculate the genome-niche index (gni) 
#load package
library(ABCoptim)
library(raster)
library(dplyr)
library(RColorBrewer)

#load data
files <- list.files(pattern = ".asc$", recursive = F)
filez <- list.files(pattern = ".csv$", recursive = F)

gni <- list()

for (i in 1:length(files)){
  
  gni[[i]] <- as.data.frame(rasterToPoints(stack(list(nsc=raster(files[i]), go=rasterFromXYZ(read.csv(filez[i])[,-1])))))

}


for (i in 1:length(gni)){
  
gni_read <- na.omit(gni[[i]])

#build function to normalise value
trans_num <- function(x, a, b){
  xmin <- min(x)
  xmax <- max(x)
  z <- ((b-a)*(x-xmin))/(xmax - xmin) + a
  return(z)
}

#normalise value beforehand
gni_read$lon_nsc <- trans_num(gni_read$nsc, 0.1, 0.9)%>%log()
gni_read$lon_go <- trans_num(gni_read$go, 0.1, 0.9)%>%log()

#build function to find the smallest value of min Y
fn_bee <- 
  function(k) {
    sum(((k*gni_read$lon_nsc+(1-k)*gni_read$lon_go)/gni_read$lon_nsc)+
          (gni_read$lon_nsc/(k*gni_read$lon_nsc+(1-k)*gni_read$lon_go))+
          ((k*gni_read$lon_nsc+(1-k)*gni_read$lon_go)/gni_read$lon_go)+
          (gni_read$lon_go/(k*gni_read$lon_nsc+(1-k)*gni_read$lon_go))-4)
  }

#ABC algorithm
abc_read <- abc_optim(0.1, fn_bee, FoodNumber = 20,
                      lb=0,ub=1,limit =50,
                      maxCycle=1000)

#calculate the gni index
gni_read$gni <- ((gni_read$nsc%>%trans_num(0.1,0.9))^abc_read[["par"]])*
  ((gni_read$go%>%trans_num(0.1,0.9))^(1-abc_read[["par"]]))

write.csv(gni_read, paste("gni_", filez[i], sep=""))
writeRaster(rasterFromXYZ(gni_read[,-c(3:6)]), paste("gni_", files[i]), "ascii")

pdf(paste(files[i], "_gni.pdf", sep=""))
pal <- colorRampPalette(c("lightgreen","darkred"))
plot(rasterFromXYZ(gni_read[,-c(3:6)]), col = pal(15), breaks = seq(0, 1, length.out=16))
dev.off()

}
