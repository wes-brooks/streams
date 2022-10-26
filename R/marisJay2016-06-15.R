library(SSN)
library(MASS)
library(viridis)

SSNpath <- "data/maris_ssn/maris4.ssn"
maris <- importSSN(SSNpath, predpts='preds')

DF = getSSNdata.frame(maris)
str(DF)
SLDF = as.SpatialLinesDataFrame(maris)
plot(SLDF)
SPDFp = as.SpatialPointsDataFrame(maris, data = 'preds')
plot(SPDFp, add = TRUE, pch = 19, col = 'red', cex = .4)
SPDF = as.SpatialPointsDataFrame(maris)
plot(SPDF, add = TRUE, col = 'blue')

xyUnique = DF[!duplicated(paste0(DF$NEAR_X,DF$NEAR_Y)),
	c('NEAR_X','NEAR_Y')]
plot(xyUnique, pch = 19, col = rgb(0,0,0,.05), cex = 4)

#kernel density estimator using MASS in a 200 x 200 grid
# using default bandwidth
kernelOut = kde2d(xyUnique$NEAR_X,
    xyUnique$NEAR_Y, n = 200)

#breaks based on quantiles
ncat = 50 # number of break categories
brks = c(min(kernelOut$z), quantile(kernelOut$z,(1:(ncat - 1))/ncat),
		max(kernelOut$z))

#quantile map
image(kernelOut, col = viridis(ncat), breaks = brks)

#breaks based on even spacing
brks = c(min(kernelOut$z), min(kernelOut$z) + 
		(range(kernelOut$z)[2] - range(kernelOut$z)[1])*(1:(ncat - 1))/ncat,
		max(kernelOut$z))

#quantile map
image(kernelOut, col = viridis(ncat), breaks = brks)

#Offset surface will be based on log of effort (count per unit effort)
logKernOut = kernelOut
logKernOut$z = log(logKernOut$z)

#breaks based on even spacing
brks = c(min(logKernOut$z), min(logKernOut$z) + 
		(range(logKernOut$z)[2] - range(logKernOut$z)[1])*(1:(ncat - 1))/ncat,
		max(logKernOut$z))

#quantile map
image(logKernOut, col = viridis(ncat), breaks = brks,
  main='Offset (Effort) Surface',
	xlab = 'NEAR_X', ylab = 'NEAR_Y')
points(xyUnique, pch = 19, col = rgb(0,0,0,.1), cex = 1)

#export it
png('/mnt/Hitachi2GB/00NMML/activePapers/fishSppDist/logOffset.png')
  image(logKernOut, col = viridis(ncat), breaks = brks,
    main='Offset (Effort) Surface',
	  xlab = 'NEAR_X', ylab = 'NEAR_Y')
  points(xyUnique, pch = 19, col = rgb(0,0,0,.1), cex = 1)
dev.off()
