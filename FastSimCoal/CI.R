library(DescTools)
files <- list.files(pattern = ".bestlhoods$", recursive = T)
bests <- do.call("rbind", lapply(files, function(fn){data.frame(Filename=fn, read.table(fn, header = T))}))
row.names(bests) <- bests[,1]
bests <- bests[,-1]
write.csv(bests, "bests.csv")
ConfI <- data.frame(matrix(ncol = 3, nrow = ncol(bests)-2), row.names = colnames(bests[,1:(ncol(bests)-2)]))
colnames(ConfI) <- c("Median", "Lower", "Upper") 
for (i in 1:(ncol(bests)-2)){
  ConfI[i,] <- MedianCI(bests[,i])
}
write.csv(ConfI, "Confidence_interval.csv")
