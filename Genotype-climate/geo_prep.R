library(sdmpredictors)
library(raster)

## download raster layers for present day
hsc <- load_layers(c(
  "BO21_chlomean_ss",
  "BO21_chlorange_ss",
  "BO21_curvelmean_ss",
  "BO21_curvelrange_ss",
  "BO21_curvelmean_bdmin",
  "BO21_curvelrange_bdmin",
  "BO21_salinitymean_ss",
  "BO21_salinityrange_ss",
  "BO21_salinitymean_bdmin",
  "BO21_salinityrange_bdmin",
  "BO21_tempmean_ss",
  "BO21_temprange_ss",
  "BO21_tempmean_bdmin",
  "BO21_temprange_bdmin"
  ))

ahsc.ext <- extent(80, 140, -20, 40)
hsc.crop <- crop(hsc, ahsc.ext)
hsc_current <- as.data.frame(hsc.crop, xy=TRUE)
write.csv(hsc_current, "current_bio.csv")

## extact env for CR
cr_loc <- read.csv("CR_loc.csv")
xy <- cr_loc[, c(2,3)]
spdf <- SpatialPointsDataFrame(coords = xy, data = cr_loc, proj4string = crs(hsc))
temp <- extract(hsc.crop, spdf)

write.csv(temp, "CR_env.csv")

## extact env for TG
tg_loc <- read.csv("TG_loc.csv")
xy <- tg_loc[, c(2,3)]
spdf <- SpatialPointsDataFrame(coords = xy, data = tg_loc, proj4string = crs(hsc))
temp <- extract(hsc.crop, spdf)

write.csv(temp, "TG_env.csv")

## extact env for TT
tt_loc <- read.csv("TT_loc.csv")
xy <- tt_loc[, c(2,3)]
spdf <- SpatialPointsDataFrame(coords = xy, data = tt_loc, proj4string = crs(hsc))
temp <- extract(hsc.crop, spdf)

write.csv(temp, "TT_env.csv")


## download raster layers for 2050 RCP26
hsc_2050_26 <- load_layers(c(
#  "BO21_RCP26_2050_chlomean_ss",
#  "BO21_RCP26_2050_chlorange_ss",
#  "BO21_RCP26_2050_curvelmean_ss",
#  "BO21_RCP26_2050_curvelrange_ss",
#  "BO21_RCP26_2050_curvelmean_bdmin",
#  "BO21_RCP26_2050_curvelrange_bdmin",
#  "BO21_RCP26_2050_salinitymean_ss",
#  "BO21_RCP26_2050_salinityrange_ss",
#  "BO21_RCP26_2050_salinitymean_bdmin",
#  "BO21_RCP26_2050_salinityrange_bdmin",
  "BO21_RCP26_2050_tempmean_ss",
  "BO21_RCP26_2050_temprange_ss"
#  "BO21_RCP26_2050_tempmean_bdmin",
#  "BO21_RCP26_2050_temprange_bdmin"
))

hsc_2050_26.crop <- crop(hsc_2050_26, ahsc.ext)
hsc_5026 <- as.data.frame(hsc_2050_26.crop, xy=TRUE)
write.csv(hsc_5026, "fut_5026.csv")

## download raster layers for 2050 RCP45
hsc_2050_45 <- load_layers(c(
  "BO21_RCP45_2050_chlomean_ss",
  "BO21_RCP45_2050_chlorange_ss",
  "BO21_RCP45_2050_curvelmean_ss",
  "BO21_RCP45_2050_curvelrange_ss",
  "BO21_RCP45_2050_curvelmean_bdmin",
  "BO21_RCP45_2050_curvelrange_bdmin",
  "BO21_RCP45_2050_salinitymean_ss",
  "BO21_RCP45_2050_salinityrange_ss",
  "BO21_RCP45_2050_salinitymean_bdmin",
  "BO21_RCP45_2050_salinityrange_bdmin",
  "BO21_RCP45_2050_tempmean_ss",
  "BO21_RCP45_2050_temprange_ss",
  "BO21_RCP45_2050_tempmean_bdmin",
  "BO21_RCP45_2050_temprange_bdmin"
))

hsc_2050_45.crop <- crop(hsc_2050_45, ahsc.ext)
hsc_5045 <- as.data.frame(hsc_2050_45.crop, xy=TRUE)
write.csv(hsc_5045, "fut_5045.csv")

## download raster layers for 2050 RCP85
hsc_2050_85 <- load_layers(c(
  "BO21_RCP85_2050_chlomean_ss",
  "BO21_RCP85_2050_chlorange_ss",
  "BO21_RCP85_2050_curvelmean_ss",
  "BO21_RCP85_2050_curvelrange_ss",
  "BO21_RCP85_2050_curvelmean_bdmin",
  "BO21_RCP85_2050_curvelrange_bdmin",
  "BO21_RCP85_2050_salinitymean_ss",
  "BO21_RCP85_2050_salinityrange_ss",
  "BO21_RCP85_2050_salinitymean_bdmin",
  "BO21_RCP85_2050_salinityrange_bdmin",
  "BO21_RCP85_2050_tempmean_ss",
  "BO21_RCP85_2050_temprange_ss",
  "BO21_RCP85_2050_tempmean_bdmin",
  "BO21_RCP85_2050_temprange_bdmin"
))

hsc_2050_85.crop <- crop(hsc_2050_85, ahsc.ext)
hsc_5085 <- as.data.frame(hsc_2050_85.crop, xy=TRUE)
write.csv(hsc_5085, "fut_5085.csv")


## download raster layers for 2100 RCP26
hsc_2100_26 <- load_layers(c(
  "BO21_RCP26_2100_chlomean_ss",
  "BO21_RCP26_2100_chlorange_ss",
  "BO21_RCP26_2100_curvelmean_ss",
  "BO21_RCP26_2100_curvelrange_ss",
  "BO21_RCP26_2100_curvelmean_bdmin",
  "BO21_RCP26_2100_curvelrange_bdmin",
  "BO21_RCP26_2100_salinitymean_ss",
  "BO21_RCP26_2100_salinityrange_ss",
  "BO21_RCP26_2100_salinitymean_bdmin",
  "BO21_RCP26_2100_salinityrange_bdmin",
  "BO21_RCP26_2100_tempmean_ss",
  "BO21_RCP26_2100_temprange_ss",
  "BO21_RCP26_2100_tempmean_bdmin",
  "BO21_RCP26_2100_temprange_bdmin"
))

hsc_2100_26.crop <- crop(hsc_2100_26, ahsc.ext)
hsc_0026 <- as.data.frame(hsc_2100_26.crop, xy=TRUE)
write.csv(hsc_0026, "fut_0026.csv")

## download raster layers for 2100 RCP45
hsc_2100_45 <- load_layers(c(
  "BO21_RCP45_2100_chlomean_ss",
  "BO21_RCP45_2100_chlorange_ss",
  "BO21_RCP45_2100_curvelmean_ss",
  "BO21_RCP45_2100_curvelrange_ss",
  "BO21_RCP45_2100_curvelmean_bdmin",
  "BO21_RCP45_2100_curvelrange_bdmin",
  "BO21_RCP45_2100_salinitymean_ss",
  "BO21_RCP45_2100_salinityrange_ss",
  "BO21_RCP45_2100_salinitymean_bdmin",
  "BO21_RCP45_2100_salinityrange_bdmin",
  "BO21_RCP45_2100_tempmean_ss",
  "BO21_RCP45_2100_temprange_ss",
  "BO21_RCP45_2100_tempmean_bdmin",
  "BO21_RCP45_2100_temprange_bdmin"
))

hsc_2100_45.crop <- crop(hsc_2100_45, ahsc.ext)
hsc_0045 <- as.data.frame(hsc_2100_45.crop, xy=TRUE)
write.csv(hsc_0045, "fut_0045.csv")

## download raster layers for 2100 RCP85
hsc_2100_85 <- load_layers(c(
  "BO21_RCP85_2100_chlomean_ss",
  "BO21_RCP85_2100_chlorange_ss",
  "BO21_RCP85_2100_curvelmean_ss",
  "BO21_RCP85_2100_curvelrange_ss",
  "BO21_RCP85_2100_curvelmean_bdmin",
  "BO21_RCP85_2100_curvelrange_bdmin",
  "BO21_RCP85_2100_salinitymean_ss",
  "BO21_RCP85_2100_salinityrange_ss",
  "BO21_RCP85_2100_salinitymean_bdmin",
  "BO21_RCP85_2100_salinityrange_bdmin",
  "BO21_RCP85_2100_tempmean_ss",
  "BO21_RCP85_2100_temprange_ss",
  "BO21_RCP85_2100_tempmean_bdmin",
  "BO21_RCP85_2100_temprange_bdmin"
))

hsc_2100_85.crop <- crop(hsc_2100_85, ahsc.ext)
hsc_0085 <- as.data.frame(hsc_2100_85.crop, xy=TRUE)
write.csv(hsc_0085, "fut_0085.csv")

