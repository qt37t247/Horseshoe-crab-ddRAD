# Analytical framework adopting Chen et al., 2022

The analyses are for three Asian horseshoe crab species 

The Mangrove horseshoe crab, CR

The Coastal horseshoe crab, TG

The Tri-spine horseshoe crab, TT


## Prepare allele frequency 

Allele frequency is calculated with PLINK v1.9 for individuals of each sampling site using "--keep" and "--freq"

We then collect all the ".frq" files to generate a allele frequency table using the R script ("freq_cal.R") attached.

The output is "*.alfreq", which is to be used as input for gradientForest


## Prepare environmental variable 

We collect environment data from BIO-ORACLE and extract site variables using the R script ("geo_prep.R") attached.

The output is "*_env.csv", which is to be used as input for gradientForest


## Compute genotype-climate association with gradientForest 

We calculate genotype-climate association with a modified R script ("*_gradientforest.R") attached.

R script can output important vairants, three pdf files for PCA, transformed genotype-climate association across the distribution ranges and the genomic offset.


## Ecological niche modelling with ENMwizard and biomod2


## Landscape genetic analyses for the connectivity


## Reference 

Chen, Y., Jiang, Z., Fan, P., Ericson, P. G., Song, G., Luo, X., ... & Qu, Y. (2022). The combination of genomic offset and niche modelling provides insights into climate change-driven vulnerability. Nature Communications, 13(1), 4821.
