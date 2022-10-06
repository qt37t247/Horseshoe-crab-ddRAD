# This script (by Joana Meier, modified by Qian Tang) calculates AIC from fsc modeling results
# Run in the folder with the highest likelihood
# Run it with R CMD BATCH AIC.R

# Read model name
args=sub('\\.tpl$','',list.files(pattern = "\\.tpl$"))

# Checks if model name was given
if(length(args)<1){
  stop("ERROR: No input / model name given\nUsage: fsc-calculateAIC.R modelname")
}

# Check if model.bestlhoods file exists
if(file.exists(paste(args[1],".bestlhoods",sep=""))){
  bestlhoods<-read.delim(paste(args[1],".bestlhoods",sep=""))
}else{
  stop(paste("ERROR: Aborted. No file ",args[1],".bestlhoods file exists",sep=""))
}

# Check if model.est file exists
if(file.exists(paste(args[1],".est",sep=""))){
  est<-readLines(paste(args[1],".est",sep=""))
}else{
  stop(paste("ERROR: Aborted. No file ",args[1],".est file exists in this directory!\nUsage: fsc-calculateAIC.R modelname",sep=""))
}

# Count number of parameters
k<-(grep("RULES",est))-(grep("//all N are",est)+1)

# Calculate AIC
AIC<-2*k-2*(bestlhoods$MaxEstLhood/log10(exp(1)))

# Calculate delta-likelihood
deltaL<-bestlhoods$MaxObsLhood-bestlhoods$MaxEstLhood

# Output model.AIC file in simulation folder
write.table(cbind(deltaL,AIC),paste(args[1],".AIC",sep=""),row.names = F,col.names = T,sep = "\t",quote = F)

