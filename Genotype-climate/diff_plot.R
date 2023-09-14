library(raster)
library(RColorBrewer)

### CR

CR <- raster("Carcinoscorpius_rotundicauda.asc")
CR_2100_RCP85 <- raster("Carcinoscorpius_rotundicauda_2100_RCP85.asc")

CR_diff <- CR - CR_2100_RCP85

writeRaster(CR_diff, "CR_diff_2100_RCP85.asc", "ascii")

pdf("CR_diff_2100_RCP85.pdf")
pal <- colorRampPalette(c("blue","white","red"))
plot(CR_diff, col = pal(15), breaks = seq(-1, 1, length.out=16))
dev.off()


### TG

TG <- raster("Tachypleus_gigas.asc")
TG_2100_RCP85 <- raster("Tachypleus_gigas_2100_RCP85.asc")

TG_diff <- TG - TG_2100_RCP85

writeRaster(TG_diff, "TG_diff_2100_RCP85.asc", "ascii", overwrite=TRUE)

pdf("TG_diff_2100_RCP85.pdf")
pal <- colorRampPalette(c("blue","white","red"))
plot(TG_diff, col = pal(15), breaks = seq(-1, 1, length.out=16))
dev.off()


### TT

TT <- raster("Tachypleus_tridentatus.asc")
TT_2100_RCP85 <- raster("Tachypleus_tridentatus_2100_RCP85.asc")

TT_diff <- TT - TT_2100_RCP85

writeRaster(TT_diff, "TT_diff_2100_RCP85.asc", "ascii")

pdf("TT_diff_2100_RCP85.pdf")
pal <- colorRampPalette(c("blue","white","red"))
plot(TT_diff, col = pal(15), breaks = seq(-1, 1, length.out=16))
dev.off()

