##################PCA######################################
library(SNPRelate)

vcf.fn <- "DATA.vcf"

snpgdsVCF2GDS(vcf.fn, "DATA.gds",  method="biallelic.only")

genofile <- openfn.gds("DATA.gds")

ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

pdf("DATA.pdf")  
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2], pch=19, col="yellow")
text(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],labels=ccm_pca$sample.id, cex=0.2)
dev.off()

tq1<-ccm_pca$eigenvect[,1]
tq2<-ccm_pca$eigenvect[,2]
tq3<-ccm_pca$eigenvect[,3]
tq4<-ccm_pca$eigenvect[,4]
pop<-ccm_pca$sample.id
tq<-cbind(tq1,tq2,tq3,tq4)
tq<-cbind(pop,tq)
tq<-as.data.frame(tq)

write.csv(tq,"pca_DATA.csv")
##I exported the PCA data to plot with excel for easier color manipulation

##################DAPC (some operations needed in R console)######################################
##Load R package
library(poppr)
DATA<-read.genepop("DATA.gen", ncode = 3)

grp <- find.clusters(DATA, max.n.clust=15)

dapc.DATA <- dapc(DATA, grp$grp)

scatter.dapc(dapc.DATA, posi.da="bottomright", bg="white", cell=0, cstar=0, solid=.4, cex=2,clab=0, leg=TRUE, scree.pca=TRUE, posi.pca="bottomleft")

compoplot(dapc.DATA)

write.csv(dapc.DATA$tab[,1:2], "dapc_pca_DATA.csv")

write.csv(dapc.DATA$grp, "dapc_member_DATA.csv")

###############AMOVA with poppr###############################
DATA<-read.PLINK("DATA.raw", map.file="DATA.map", parallel = FALSE)

##Asign individuals to population (Using T. tridentatus as example)
pop <- c("TS", "TS", "JP", "JP", "TS", "TS", "JP", "JP", "JP", "TS", "TS", "TS", "JP", "TS", "TS", "TS", "TS", "JP", "TS", "TS", "TS", "TS", "JP", "TS", "TS", "TS", "JP", "JP", "TS", "TS", "TS", "JP", "JP", "TS", "TS", "TS", "JP", "JP", "TS", "TS", "TS", "JP", "JP", "TS", "TS", "TS", "JP", "JP", "TS", "TS", "TS", "JP", "JP", "TS", "TS", "JP", "JP", "TS", "TS", "JP", "JP", "TS") 
pop(DATA) <- pop

amo<-poppr.amova(DATA, ~pop, within = TRUE)

amov <- cbind(amo$results,amo$componentsofcovariance)

rownames(amov) <- c("Between pop", "Within pop", "Within samples", "Total")

write.csv(amov, "amova.csv")

##################Population-based statistics######################################
##Load the library
library(PopGenome)

DATA <- readData("DIRECTORY/DATA", format="VCF", SNP.DATA = TRUE)

##Asign individuals to population (Using T. tridentatus as example)
pop.A <- c("DJ15","DJF01002","DJF01003","DJF01005","DJF01006","DJN02001","DJN02002","DJN02003","DJN02004","DJN02005","DJN02006","DJN03002","DJN03003","DJN03004","DJO01001","DJO01002","DJO02001","DJS01002","DJS01004","DJY01002","DJY01003","DJY01004","DJY01005","DJY01006")

pop.B <- c("DEK01001","DEM02","DEM03","DEM05","DPP01002","DPP01003","DPP01004","DPP01005","DPP01006","TTW02","TTW03","TTW06","TTW12","TTW15","DBH03","DBH04","DBH05","DBH08","DBH09","DGD01","DGD02","DGD03","DGD04","DGD05","DJA01001","DJA01002","DJA01003","DSY01","DSY02","DSY03","DVD01","DVD02","DVD03","DVF01","DVF02","DVG08","DZJ01","DZJ02")

DATA <- set.populations(DATA,list(pop.A, pop.B))

DATA <- neutrality.stats(DATA)

DATA <- F_ST.stats(DATA)

##Get population pairwise Fst
DATA@nuc.F_ST.pairwise

##Get nucleotide diversity
DATA@nuc.diversity.within

##Get population pairwise Dxy
DATA@nuc.diversity.between

##Get Tajima's D
DATA@Tajima.D
