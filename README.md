# Scripts used for the ddRADSeq data of three Asian horseshoe crab species

In this project, we collected over 250 horseshoe crabs from 10 counties around Asia. We produced first comprehesive population genetic structure for the Asian horseshoe crab populations, reconstructed the evolutionary history, and projected their viability under impending anthropogenic climate change scenarios.

Manuscript title: Evolution and viability of Asian horseshoe crabs appear tightly linked to geo-climatic dynamics in Sunda Shelf.


## Data availability

Sequencing reads are archived in NCBI under Bioproject PRJNA1127623. 

To download the data, you may use the command below (may need to install relevant NCBI tools accordingly):

```bash
project='PRJNA1127623'
esearch -db sra -query $project | efetch -format runinfo > runinfo.csv
cat runinfo.csv | cut -d "," -f 1 > SRR.numbers
sed '1d' SRR.numbers | parallel -j 12 fastq-dump --split-files --origfmt --gzip {}
```


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
