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

tt_loc <- read.csv("TT_loc.csv")
xy <- tt_loc[, c(2,3)]
spdf <- SpatialPointsDataFrame(coords = xy, data = tt_loc, proj4string = crs(hsc))
temp <- extract(hsc.crop, spdf)

write.csv(temp, "TT_env.csv")

