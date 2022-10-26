pairdist2.lpp <- function(X, net, ssnobj, lookup, predpts=NULL, binTable=NULL) {
    
    #mylpp <- lpp(X$coxdata, Q$domain)
    
    # create a data.frame of data at the points
    point.data <- data.frame(
        ratio = 1 - X$data$tp,
        NEAR_X = X$data$x,
        NEAR_Y = X$data$y
    )
    point.data$pid <- 1:nrow(point.data)
    point.data$locID <- point.data$pid
    
    indx <- match(X$data$seg, lookup$seg)
    point.data$rid <- as.factor(lookup$rid[indx])
    point.data$netID <- as.factor(lookup$netID[indx])
    
    indx <- match(point.data$rid, lookup$rid)
    point.data$upDist <- lookup$upDist[indx] - X$data$tp * lookup$Length[indx]

    # create a data frame holding the coordinates of each point on the plane
    pcoords <- data.frame(coords.x1 = point.data$x, coords.x2 = point.data$y)
    pcoords <- as.matrix(pcoords)
        
    # create a data.frame giving the coordinates of each point on the network system
    net.coords <- data.frame(DistanceUpstream = point.data$upDist)
    net.coords$NetworkID <- point.data$netID
    net.coords$SegmentID <- point.data$rid
    attr(net.coords, 'locID') <- point.data$locID
    
    # create the SpatialStreamNetwork with its component parts\
    ssnobj@obspoints@SSNPoints[[1]] <- new('SSNPoint', network.point.coords=net.coords, point.data=point.data, points.bbox=ssnobj@bbox, point.coords=pcoords)
    
    # matrices of distance and connection
    D <- createDistMatInMemory(ssnobj, predpts, binTable=binTable)
    C <- list()
    for (j in names(D[['obs']])) if (!is.null(D[['obs']][[j]])) C[[j]] <- Matrix(TRUE, nrow(D[['obs']][[j]]), ncol(D[['obs']][[j]])) 
    
    dist <- bdiag(D[['obs']])
    dist <- dist + t(dist)
    
    flag <- bdiag(C)
    dist[!flag] <- Inf
    
    # reorder the output to match the ordering of the inputs
    reorder <- unlist(sapply(names(D$obs), function(nm) {which(ssnobj@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID == nm)}))
    dist <- dist[order(reorder), order(reorder)]
    dist
}
