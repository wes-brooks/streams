# import and clean the observations, quadrature points, and stream network data:
library("SSN")
library('dplyr')
library('sp')
library('rgdal')
library('stringi')
library('Matrix')
library("TMB")
library('MASS')
library('spatstat')



load('data/intermediate-products/maris-ssn.rdata')
source('code/spatstat-override/linnet.R')

# create the rid.lookup object:
rid <- edges@data$rid
rid.lookup <- data.frame(rid=rid, seg=1:length(rid))
colnames(rid.lookup) <- c('rid', 'seg')

# extract the endpoints of the stream segments
ends <- sapply(edges@lines, function(x) c(head(x@Lines[[1]]@coords, 1), tail(x@Lines[[1]]@coords, 1))) %>% t
ends <- cbind(ends, edges@data$rid, edges@data$netID)
ends <- as.data.frame(rbind(ends[,c(1, 2, 5, 6)], ends[,c(3, 4, 5, 6)]))
colnames(ends) <- c('x', 'y', 'rid', 'netID')

# fix some errors in the data
ends$netID[ends$netID == 30] <- 82
ends$netID[ends$netID == 84] <- 82
ends$netID[ends$netID == 87] <- 86
ends$netID[ends$netID == 91] <- 90



# create a ppp object holding the vertices
win <- owin(xrange = range(ends$x), yrange = range(ends$y))
pts <- ends %>% group_by(x, y) %>% distinct(.keep_all = TRUE)
pts$pid <- 1:nrow(pts)
pts$rid <- NULL
ends <- right_join(pts, ends, by=c('x', 'y', 'netID'))
verts <- ppp(x=pts$x, y=pts$y, window=win)


# create the linnet object defining the network
m <- Matrix(FALSE, nrow(pts), nrow(pts)) # matrix indicating which points are connected
m.rid <- Matrix(0, nrow(pts), nrow(pts)) # matrix indicating the reach id of connected points

# identify which reaches touch each other
for (pp in unique(ends$pid)) {
  cat(paste(pp, '\n'))
  tangent <- ends$rid[ends$pid == pp]
  connected <- ends$pid[ends$rid %in% tangent & ends$pid != pp]
  m[pp, connected] <- m[connected, pp] <- TRUE
  
  m.rid[pp, ends$pid[ends$rid %in% tangent & ends$pid != pp]] <- ends$rid[ends$rid %in% tangent & ends$pid != pp]
}




# create a spatstat linear network object using sparse matrices
# create a linear network object using sparse matrices
net <- MyLinnet(vertices=verts, m=m, m.rid=m.rid, sparse=TRUE, warn=FALSE)

# this is the key that matches rid to seg
indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$netID <- edges@data$netID[indx]
rid.lookup$Length <- edges@data$Length[indx]
rid.lookup$upDist <- edges@data$upDist[indx]

# match rid to seg on the linear network
indx <- match(net$lines$rid, rid.lookup$rid)
net$lines$seg <- rid.lookup$seg[indx]
net$lines$Length <- rid.lookup$Length[indx]






# New GMRF precision matrix based on individual stream segments:
# First, let's consider only rids on networks with at least one point.
obs.nets <- as.numeric(as.character(unique(maris@obspoints@SSNPoints[[1]]@point.data$netID)))
indx.valid.net <- edges@data$netID %in% obs.nets
gmrf.lookup <- edges@data[indx.valid.net, c('netID', 'rid')]
gmrf.lookup$gmrfID <- 1:nrow(gmrf.lookup)

# define the random effects precision matrix
Sigma.mask <- Matrix(0, sum(indx.valid.net), sum(indx.valid.net))

# populate the random effects precison matrix of the GMRF
for (i in 1:nrow(edges@data)) {
  colID <- gmrf.lookup$rid == edges@data$rid[i]

  if (edges@data$netID[i] %in% obs.nets && !is.na(edges@data$downstr[i])) {
    downstrID <- gmrf.lookup$rid == edges@data$downstr[i]
    Sigma.mask[downstrID, colID] <- Sigma.mask[colID, downstrID] <- -1
  }
}
diag(Sigma.mask) <- -colSums(Sigma.mask)

save(m, m.rid, net, Sigma.mask, file='data/intermediate-products/stream-network-structure.rdata')