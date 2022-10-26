# import and clean the MARIS and SSN data:
library("SSN")
library('dplyr')
library('sp')
library('rgdal')
library('stringi')
library('Matrix')
library("TMB")
library('MASS')
library('spatstat')

ssnpath <- "data/maris_ssn/maris4.ssn"

# import the data
maris <- importSSN(ssnpath, predpts='preds')

# import the data supplement that gives full coverage for the quadrature knots
supplement <- read.csv("data/quadpts.csv")
supplement$OBSPRED_ID <- 1:nrow(supplement) + max(maris@predpoints@SSNPoints[[1]]@point.data$OBSPRED_ID)
colnames(supplement)[23:24] <- c('NEAR_X', 'NEAR_Y')
supplement$NEAR_ANGLE <- 0
supplement$NEAR_DIST <- 0
supplement$pid <- 1:nrow(supplement) + max(maris@predpoints@SSNPoints[[1]]@point.data$pid)
supplement$locID <- supplement$siteID <- as.factor(1:nrow(supplement) + max(as.integer(levels(maris@predpoints@SSNPoints[[1]]@point.data$locID))))
colnames(supplement)[colnames(supplement) == 'ELEV_NSI'] <- 'Elev_NSI'
supplement$ratio <- 0.5
supplement$NEAR_FID <- -1
supplement$netID <- as.factor(supplement$netID)

supplement$FID <- NULL
supplement$ORIG_FID <- NULL
supplement$quad <- NULL
supplement$COMID <- NULL
supplement$GNIS_NAME <- NULL
supplement$FTYPE <- NULL
supplement$areaPI <- NULL
supplement$Length <-  NULL

s2 <- rbind(maris@predpoints@SSNPoints[[1]]@point.data, supplement)

maris@predpoints@SSNPoints[[1]]@point.data <- s2
maris@predpoints@SSNPoints[[1]]@point.coords <- as.matrix(s2[,c('NEAR_X', 'NEAR_Y')])
colnames(maris@predpoints@SSNPoints[[1]]@point.coords) <- c('coords.x1', 'coords.x2')
maris@predpoints@SSNPoints[[1]]@network.point.coords <- s2[,c('netID', 'rid', 'upDist')]
colnames(maris@predpoints@SSNPoints[[1]]@network.point.coords) <- c('NetworkID', 'SegmentID', 'DistanceUpstream')


# Get the IDs of distinct networks in the prediction dataset
netids.pred <- unique(maris@predpoints@SSNPoints[[1]]@point.data$netID)
netids.pred <- as.integer(levels(netids.pred)[netids.pred])

# Get the IDs of distinct networks in the observation data
netids.obs <- unique(maris@obspoints@SSNPoints[[1]]@point.data$netID)
netids.obs <- as.integer(levels(netids.obs)[netids.obs])

# import stream shapefiles
edges <- readOGR(paste(ssnpath, 'edges.shp', sep='/'), layer='edges')
netids <- unique(edges@data$netID)

# correct some small errors in the data
edges@data$upDist[edges@data$netID == 30] <- edges@data$upDist[edges@data$netID == 30] + edges@data$upDist[edges@data$rid == 12473]
edges@data$upDist[edges@data$netID == 84] <- edges@data$upDist[edges@data$netID == 84] + edges@data$upDist[edges@data$rid == 20199]
edges@data$upDist[edges@data$netID == 87] <- edges@data$upDist[edges@data$netID == 87] + edges@data$upDist[edges@data$rid == 20212]
edges@data$upDist[edges@data$netID == 91] <- edges@data$upDist[edges@data$netID == 91] + edges@data$upDist[edges@data$rid == 22498]
edges@data$netID[edges@data$netID == 30] <- 82
edges@data$netID[edges@data$netID == 84] <- 82
edges@data$netID[edges@data$netID == 87] <- 86
edges@data$netID[edges@data$netID == 91] <- 90

# drop the reach with rid 39142 because it connects two points that are also connected by 39143
indx <- edges@data$rid != 39142
edges@lines <- edges@lines[indx]
edges@data <- edges@data[indx,]

# import the binary codes that indicate flow connections in a stream network, and remove whitespace
connect <- dbConnect(RSQLite:::SQLite(), paste(ssnpath, "binaryID.db", sep='/'))
bincodes <- list()
meta <- matrix(NA, 0, 2)
run <- 0
for (netID in netids) {
  bincodes[[netID]] <- dbReadTable(connect, paste('net', netID, sep=''))
  bincodes[[netID]]$binaryID <- sapply(bincodes[[netID]]$binaryID, function(z) stri_replace_all_charclass(z, "\\p{WHITE_SPACE}", ""))
  run <- run + sum(edges@data$Length[edges@data$netID == netID])
  
  meta <- rbind(meta, cbind(bincodes[[netID]], netID))
}
meta <- as.data.frame(meta)

# Fix some minor data errors:

# network 30 is actually an offshoot of network 82
bincodes[[30]]$binaryID <- paste0(bincodes[[82]]$binaryID[bincodes[[82]]$rid == 12473], bincodes[[30]]$binaryID)
bincodes[[82]] <- rbind(bincodes[[82]], bincodes[[30]])
bincodes[[30]] <- data.frame(rid=vector(), binaryID=vector())
meta$netID[meta$netID == 30] <- 82
meta$binaryID[meta$netID == 82] <- bincodes[[82]]$binaryID

# network 84 is actually an offshoot of network 82
bincodes[[84]]$binaryID <- paste0(bincodes[[82]]$binaryID[bincodes[[82]]$rid == 20199], bincodes[[84]]$binaryID)
bincodes[[82]] <- rbind(bincodes[[82]], bincodes[[84]])
bincodes[[84]] <- data.frame(rid=vector(), binaryID=vector())
meta$netID[meta$netID == 84] <- 82
meta$binaryID[meta$netID == 82] <- bincodes[[82]]$binaryID

# network 87 is actually an offshoot of network 86
bincodes[[87]]$binaryID <- paste0(bincodes[[86]]$binaryID[bincodes[[86]]$rid == 20212], c('2', '20'))
bincodes[[86]] <- rbind(bincodes[[86]], bincodes[[87]])
bincodes[[87]] <- data.frame(rid=vector(), binaryID=vector())
meta$netID[meta$netID == 87] <- 86
meta$binaryID[meta$netID == 86] <- bincodes[[86]]$binaryID

# network 91 is actually an offshoot of network 90
bincodes[[91]]$binaryID <- paste0(bincodes[[90]]$binaryID[bincodes[[90]]$rid == 22498], bincodes[[91]]$binaryID)
bincodes[[90]] <- rbind(bincodes[[90]], bincodes[[91]])
bincodes[[91]] <- data.frame(rid=vector(), binaryID=vector())
meta$netID[meta$netID == 91] <- 90
meta$binaryID[meta$netID == 90] <- bincodes[[90]]$binaryID

# add the binaryID column to the edges data.
indx <- match(edges@data$rid, meta$rid)
edges@data$binaryID <- meta$binaryID[indx]

# identify the downstream reach for each reach
downstr <- vector()
for (i in 1:nrow(edges@data)) {
  cat(i, '\n')
  b.ref <- edges@data$binaryID[i]
  len.ref <- nchar(b.ref)
  indx.net <- edges@data$netID == edges@data$netID[i]
  
  # The downstream reach is the one on this net that matches the first n-1 bits of the current reach.
  edges.net <- edges@data[indx.net,]
  downstr.indx <- edges.net$binaryID == substr(b.ref, 1, len.ref - 1)
  
  if (any(downstr.indx)) {
    downstr <- c(downstr, edges.net$rid[downstr.indx][1])
  } else downstr <- c(downstr, NA)
}
edges@data$downstr <- downstr

# Remove observation points on rids that have no quadrature knots (there are only a few)
obs.rids <- unique(maris@obspoints@SSNPoints[[1]]@point.data$rid)
rid.noquad <- obs.rids[!(obs.rids %in% maris@predpoints@SSNPoints[[1]]@point.data$rid)]
indx.noquad <-  maris@obspoints@SSNPoints[[1]]@point.data$rid %in% rid.noquad
maris@obspoints@SSNPoints[[1]]@point.data <- maris@obspoints@SSNPoints[[1]]@point.data[!indx.noquad,]

save(maris, edges, file='data/intermediate-products/maris-ssn.rdata')