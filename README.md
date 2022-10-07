# Scripts used for the ddRADSeq data of three horseshoe crab species

Manuscript title: The Quaternary fate of horseshoe crabs – a living fossil – appears tightly linked to sea level  dynamics in Sundaland.

Additional information:


## Step 1. NGS reads alignment (Align.txt)

Programs used: 

STACKS (process_radtags): https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php 

BWA (MEM): https://github.com/lh3/bwa

SAMTOOLS: http://www.htslib.org/


## Step 2a. SNP calling and filtering (SNPs.txt)

Programs used: 

STACKS (ref_map.pl): https://catchenlab.life.illinois.edu/stacks/comp/ref_map.php 

PLINK: https://www.cog-genomics.org/plink/


## Step 2b. Site frequency spectrum calculation (SFS.txt)

Program used: 

ANGSD (angsd & realSFS): http://www.popgen.dk/angsd/index.php/SFS_Estimation


## Step 3a. R-based population genetic analyses (R-based_analyses.R)

Packages used: 

SNPRelate: https://bioconductor.org/packages/release/bioc/html/SNPRelate.html

poppr: https://cran.r-project.org/web/packages/poppr/index.html

PopGenome: https://cran.r-project.org/web/packages/PopGenome/index.html


## Step 3b. Reconstruction of demographical history

Scripts for estimating effective population size over time (Folder "Stairwayplot")
Stairway Plot2: https://github.com/xiaoming-liu/stairway-plot-v2

Scripts for simulating demographical events (Folder "FastSimCoal")
FastSimCoal2: http://cmpg.unibe.ch/software/fastsimcoal27/
