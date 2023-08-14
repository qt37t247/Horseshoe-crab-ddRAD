## Collect allele frequency files for all sampling sites
files <- list.files(pattern = ".frq$", recursive = T)
frqs <- do.call("cbind", lapply(files, function(fn){data.frame(Filename=fn, read.table(fn, sep = "" , header = T , na.strings ="", stringsAsFactors= F))}))
write.csv(frqs, "raw_frequncy.csv")

## Retain allele frequency only, other information convert to row and column names
frqx <- frqs[,(1:length(files)*7-1]
colnames(frqx) <- files
rownames(frqx) <- gsub(" ", "", paste(frqs$CHR, "_", frqs$SNP))
write.csv(t(frqx), "*.alfreq") ## make it universal
