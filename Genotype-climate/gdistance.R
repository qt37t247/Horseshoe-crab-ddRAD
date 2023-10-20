library(poppr)
gen_data <- adegenet::read.genepop('CR0L.gen', ncode=3)
gen_dist <- poppr::diss.dist(gen_data)
gen_dis<-as.data.frame(as.matrix(dist(gen_dist)))
write.table(gen_dis, 'CR0L_gdist.txt')

gen_data <- adegenet::read.genepop('TG0L.gen', ncode=3)
gen_dist <- poppr::diss.dist(gen_data)
gen_dis<-as.data.frame(as.matrix(dist(gen_dist)))
write.table(gen_dis, 'TG0L_gdist.txt')

gen_data <- adegenet::read.genepop('TT0L.gen', ncode=3)
gen_dist <- poppr::diss.dist(gen_data)
gen_dis<-as.data.frame(as.matrix(dist(gen_dist)))
write.table(gen_dis, 'TT0L_gdist.txt')

library(raster)
sp_data <- read.csv("CR_loc.csv")
xy <- sp_data[, c(2,3)]
sp_data <- SpatialPointsDataFrame(coords = xy, data = sp_data, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
sp_dis <- spDists(sp_data, longlat = TRUE)
write.table(sp_dis, 'CR_dist.txt')

sp_data <- read.csv("TG_loc.csv")
xy <- sp_data[, c(2,3)]
sp_data <- SpatialPointsDataFrame(coords = xy, data = sp_data, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
sp_dis <- spDists(sp_data, longlat = TRUE)
write.table(sp_dis, 'TG_dist.txt')

sp_data <- read.csv("TT_loc.csv")
xy <- sp_data[, c(2,3)]
sp_data <- SpatialPointsDataFrame(coords = xy, data = sp_data, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
sp_dis <- spDists(sp_data, longlat = TRUE)
write.table(sp_dis, 'TT_dist.txt')
