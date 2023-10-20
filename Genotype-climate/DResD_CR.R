#load genetic distance matrix (lower-left triangle)
M.gen.diss <- read.table(file = 'CR0L_gdist.txt', sep = '\t', header = T, row.names = 1, fill = T)
M.gen.diss <- as.dist(M.gen.diss, diag = TRUE, upper = TRUE)
#load geographical coordinates and create geo distance matrix
sample.points <- read.table(file = 'CR.coordinates', sep = '\t', header = T, row.names = 1, fill = T)
names(sample.points) <- c('x', 'y')
M.geo.dist <- read.table(file = 'CR_dist.txt', sep = '\t', header = T, row.names = 1, fill = T)
M.geo.dist <- as.dist(M.geo.dist, diag = TRUE, upper = TRUE)
#Data preparation
IBD.data <- data.frame(cbind(as.vector(M.geo.dist), as.vector(M.gen.diss)))
names(IBD.data) <- c('geo.dist', 'gen.diss')
#Estimate of the IBD 
IBD.model <- lm(IBD.data$gen.diss~IBD.data$geo.dist)
###Plot IBD###
pdf("CR_IBD.pdf")
plot(IBD.data, pch=19, col="orange")
abline(IBD.model, col="orange")
points(IBD.data, pch=15, col="blue")
abline(IBD.model, col="blue")
points(IBD.data, pch=17, col="green")
abline(IBD.model, col="green")

dev.off()
#estimate IBD residual matrix
M.IBD.res <- M.gen.diss - (coef(IBD.model)[1] + coef(IBD.model)[2] * M.geo.dist)
#Define grid step and margin width 
step <- 0.1
margin <- 3
#Define grid area
x.min <- min(sample.points$x) - margin
x.max <- max(sample.points$x) + margin
y.min <- min(sample.points$y) - margin
y.max <- max(sample.points$y) + margin
#Generate grid table
grid <- expand.grid(x = seq(from = x.min, to = x.max, by = step), y = seq(from = y.min, to = y.max, by = step))
names(grid) <- c('x', 'y')
#Calculate midpoints
n <- dim(sample.points)[1]
M.x1 <- as.dist(matrix(data = sample.points$x, nrow = n, ncol = n, byrow = F))
M.y1 <- as.dist(matrix(data = sample.points$y, nrow = n, ncol = n, byrow = F))
M.x2 <- as.dist(matrix(data = sample.points$x, nrow = n, ncol = n, byrow = T))
M.y2 <- as.dist(matrix(data = sample.points$y, nrow = n, ncol = n, byrow = T))
M.mid.x <- (M.x1 + M.x2) / 2
M.mid.y <- (M.y1 + M.y2) / 2
min.dist <- 0
max.dist <- 5000
midpoints <- data.frame(as.vector(M.mid.x[(M.geo.dist > min.dist) & (M.geo.dist < max.dist)]),
	as.vector(M.mid.y[(M.geo.dist > min.dist) & (M.geo.dist < max.dist)]),
	as.vector(M.IBD.res[(M.geo.dist > min.dist) & (M.geo.dist < max.dist)]))
names(midpoints) <- c('x', 'y', 'z')
#Calculate weighted variance
weighted.var <- function(x, w, na.rm = FALSE){
	if (na.rm){
		w <- w[i <- !is.na(x)]
		x <- x[i]
		}
	sum.w <- sum(w)
	sum.w2 <- sum(w ** 2)
	mean.w <- sum(x * w) / sum(w)
	(sum.w / (sum.w ** 2 - sum.w2)) * sum(w * (x - mean.w) ** 2, na.rm = na.rm)
	}
shift <- 1
calculation <- function(gridx, gridy, midpoints, shift){
	w <- (shift + sqrt((gridx - midpoints[,1]) ** 2 + (gridy - midpoints[,2]) ** 2)) ** -2
	return(c(weighted.mean(midpoints[,3], w), weighted.var(midpoints[,3], w) ** 0.5))
	}
grid <- cbind(grid, t(mapply(calculation, grid[,1], grid[,2], MoreArgs = list(midpoints, shift))))
names(grid) <- c('x', 'y', 'z', 'SD')
n.resamp <- 1000
grid <- cbind(grid, 0, 0, 0)
names(grid) <- c(names(grid)[1:4], 'prob.z>H0', 'H0.cl1', 'H0.cl2')
H0.lower.tail <- H0.upper.tail <- matrix(nrow = dim(grid)[1], ncol = floor(n.resamp * 0.025))
perm.rowscols = function (m) {
	d <- dim(m)
	s <- sample(1:d[1])
	m[s,s]
	}
#Randomization
for(i in 1:n.resamp) {
	M.IBD.res.resamp <- as.dist(perm.rowscols(as.matrix(M.IBD.res)))
	midpoints.resamp <- data.frame(midpoints[,1:2], as.vector(M.IBD.res.resamp[(M.geo.dist > min.dist) & (M.geo.dist < max.dist)]))
	resamp <- t(mapply(calculation, grid[,1], grid[,2], MoreArgs = list(midpoints.resamp, step)))
	grid[,5] <- grid[,5] + (grid[,3] > resamp[,1]) / n.resamp

	if(i < n.resamp * 0.025) {
		H0.lower.tail[,i] <- resamp[,1]
		H0.upper.tail[,i] <- resamp[,1]
		}
	else {

		H0.lower.tail <- t(apply(cbind(resamp, H0.lower.tail), MARGIN = 1 , FUN = sort))[,1:dim(H0.lower.tail)[2]]
		H0.upper.tail <- t(apply(cbind(resamp, H0.upper.tail), MARGIN = 1 , FUN =sort))[,2:(dim(H0.lower.tail)[2]+1)]
		}
	next
	}

grid[,6] <- apply(H0.lower.tail, MARGIN = 1, FUN = max)
grid[,7] <- apply(H0.upper.tail, MARGIN = 1, FUN = min)
#Bootstrapping
n.bts <- 250
grid <- cbind(grid, 1)
names(grid) <- c(names(grid)[1:7], 'power')
for(i in 1:n.bts) {
	s <- sample(1:dim(midpoints)[1], replace = T)
	midpoints.bts <- midpoints[s,]
	bts <- t(mapply(calculation, grid[,1], grid[,2], MoreArgs = list(midpoints.bts,step)))
	grid[,8] <- grid[,8] - (grid[,6] < bts[,1] & grid[,7] > bts[,1]) / n.bts
	next
	}
#Output results
write.csv(midpoints, file ='CR0_MidPoint.csv')
write.csv(grid, file = 'CR0_DresD.csv')

CR <- read.csv("CR_DresD.csv")


##########Visualization##############
CR <- read.csv("CR_DresD.csv")
grid <- CR[, -1]
colnames(grid) <- c("x", "y", "z", "SD", "prob", "cl1", "cl2", "power")

library(ggplot2)
ggplot(grid, aes(x = x, y = y, colour=z, fill = z)) + 
  geom_point(size=2, shape=15) +
  scale_color_gradient2(low="darkgreen", high="red", mid="white", midpoint=0)+
  geom_point(data = sample.points, shape=21, color="black", fill="yellow",stroke=1)+
  stat_contour(aes(z = power),bins=4,size=0.5,col="red",breaks=c(0.8))+
  stat_contour(aes(z = prob),bins=4,size=0.5,col="purple",breaks=c(0.025))+
  stat_contour(aes(z = prob),bins=4,size=0.5,col="purple",breaks=c(0.975))

library(dplyr)
gridf <- grid %>% filter(prob > 0.975|prob<0.025)
gridf <- gridf %>% filter(power > 0.8)

write.csv(gridf, "CR_DRD.csv")

