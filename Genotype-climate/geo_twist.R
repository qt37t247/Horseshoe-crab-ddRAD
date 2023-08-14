xx <- load_layers(c(
  "BO21_RCP85_2100_chlomean_ss",
  "BO21_RCP85_2100_chlorange_ss",
  "BO21_RCP85_2100_curvelmean_ss",
  "BO21_RCP85_2100_curvelrange_ss"
#  "BO21_RCP85_2100_curvelmean_bdmin",
#  "BO21_RCP85_2100_curvelrange_bdmin",
#  "BO21_RCP85_2100_salinitymean_ss",
#  "BO21_RCP85_2100_salinityrange_ss",
#  "BO21_RCP85_2100_salinitymean_bdmin",
#  "BO21_RCP85_2100_salinityrange_bdmin",
#  "BO21_RCP85_2100_tempmean_ss",
#  "BO21_RCP85_2100_temprange_ss",
#  "BO21_RCP85_2100_tempmean_bdmin",
#  "BO21_RCP85_2100_temprange_bdmin"
))

xx_1 <- load_layers(c(
#  "BO21_RCP85_2100_chlomean_ss",
#  "BO21_RCP85_2100_chlorange_ss",
#  "BO21_RCP85_2100_curvelmean_ss",
#  "BO21_RCP85_2100_curvelrange_ss",
  "BO21_RCP85_2100_curvelmean_bdmin",
  "BO21_RCP85_2100_curvelrange_bdmin",
  "BO21_RCP85_2100_salinitymean_ss",
  "BO21_RCP85_2100_salinityrange_ss"
#  "BO21_RCP85_2100_salinitymean_bdmin",
#  "BO21_RCP85_2100_salinityrange_bdmin",
#  "BO21_RCP85_2100_tempmean_ss",
#  "BO21_RCP85_2100_temprange_ss",
#  "BO21_RCP85_2100_tempmean_bdmin",
#  "BO21_RCP85_2100_temprange_bdmin"
))

xx_2 <- load_layers(c(
#  "BO21_RCP85_2100_chlomean_ss",
#  "BO21_RCP85_2100_chlorange_ss",
#  "BO21_RCP85_2100_curvelmean_ss",
#  "BO21_RCP85_2100_curvelrange_ss",
#  "BO21_RCP85_2100_curvelmean_bdmin",
#  "BO21_RCP85_2100_curvelrange_bdmin",
#  "BO21_RCP85_2100_salinitymean_ss",
#  "BO21_RCP85_2100_salinityrange_ss",
  "BO21_RCP85_2100_salinitymean_bdmin",
  "BO21_RCP85_2100_salinityrange_bdmin",
  "BO21_RCP85_2100_tempmean_ss",
  "BO21_RCP85_2100_temprange_ss"
#  "BO21_RCP85_2100_tempmean_bdmin",
#  "BO21_RCP85_2100_temprange_bdmin"
))

xx_3 <- load_layers(c(
#  "BO21_RCP85_2100_chlomean_ss",
#  "BO21_RCP85_2100_chlorange_ss",
#  "BO21_RCP85_2100_curvelmean_ss",
#  "BO21_RCP85_2100_curvelrange_ss",
#  "BO21_RCP85_2100_curvelmean_bdmin",
#  "BO21_RCP85_2100_curvelrange_bdmin",
#  "BO21_RCP85_2100_salinitymean_ss",
#  "BO21_RCP85_2100_salinityrange_ss",
#  "BO21_RCP85_2100_salinitymean_bdmin",
#  "BO21_RCP85_2100_salinityrange_bdmin",
#  "BO21_RCP85_2100_tempmean_ss",
#  "BO21_RCP85_2100_temprange_ss",
  "BO21_RCP85_2100_tempmean_bdmin",
  "BO21_RCP85_2100_temprange_bdmin"
))

ahsc.ext <- extent(80, 140, -20, 40)

xx.crop <- crop(xx, ahsc.ext)
xxx <- as.data.frame(xx.crop, xy=TRUE)
xx_1.crop <- crop(xx_1, ahsc.ext)
xxx_1 <- as.data.frame(xx_1.crop, xy=TRUE)
xx_2.crop <- crop(xx_2, ahsc.ext)
xxx_2 <- as.data.frame(xx_2.crop, xy=TRUE)
xx_3.crop <- crop(xx_3, ahsc.ext)
xxx_3 <- as.data.frame(xx_3.crop, xy=TRUE)

xxx <- cbind(xxx, xxx_1[,-(1:2)])
xxx <- cbind(xxx, xxx_2[,-(1:2)])
xxx <- cbind(xxx, xxx_3[,-(1:2)])

write.csv(xxx, "fut_0085.csv")

rm(list=ls())