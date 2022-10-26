# getting the number of active edges at a given distance from the observation point:
# lpp <- bt.lpp


nedges <- function(netnum, lpp, edges, distances) {
  i <- 1
  ndist <- nrow(distances)
  
  cnt.ends <- matrix(NA, ndist, 0)
  
  # for (netnum in unique(lpp$data$netID)) {
  # for (netnum in c(23)) {
    
    # retrieve the serialized network of distances
    f <- file(paste0("~/Dropbox/streams/data/maris_ssn/maris4.ssn/distance/obs/dist.net.new.", netnum, ".RData"), open='rb')
    D <- unserialize(f)
    close(f)

    D <- D$obs[[1]]
    
    # subset the data objects to just focus on this network
    bb <- edges[edges$netID == netnum,]
    ppdata <- as.data.frame(lpp$data[lpp$data$netID == netnum,])
    
    indx <- match(ppdata$rid, bb$rid)
    ppdata$upDist <- bb$upDist[indx] - ppdata$tp * bb$Length[indx]
    
    # calculate the distances along the stream network
    for (rid in unique(ppdata$rid)) {
      # consider simultaneously all the points that are on the same stream reach
      indx.rid <- which(bb$rid == rid)[1]
      vertex.dist <- D[indx.rid,] + D[,indx.rid]
      # vertex.dist <- c(vertex.dist, bb$upDist[bb$rid == rid])
      
      # subset the matrix of distances between points to focus on those on this reach.
      indx <- lpp$data$rid == rid
      empiricDistanceLocal <- distances[, indx, drop=FALSE]
      
      first.pt <- TRUE
      indx <- ppdata$rid == rid
      for (j in 1:sum(indx)) {
        id <- which(indx)[j]
        cat(i, "\n")
        i <- i+1 
        
        if (first.pt) {
          d1 <- bb$upDist[indx.rid] - ppdata$upDist[[id]]
        } else {
          d1 <- uD.prev - ppdata$upDist[[id]]
        }
        uD.prev <- ppdata$upDist[[id]]
        first.pt <- FALSE
        
        # upstr.indx <- substr(c(bb$binaryID, '0'), 1, nchar(bb$binaryID[indx.rid])) == bb$binaryID[indx.rid]
        upstr.indx <- substr(bb$binaryID, 1, nchar(bb$binaryID[indx.rid])) == bb$binaryID[indx.rid]
        
        downstr.binary <- substr(bb$binaryID[indx.rid], 1, nchar(bb$binaryID[indx.rid]) - 1)
        # if (nchar(bb$binaryID[indx.rid])==1) {
        #   downstr.indx <- c(rep(FALSE, nrow(bb)), TRUE)
        # } else {
        #   downstr.indx <- c(bb$binaryID == downstr.binary, FALSE)
        # }
        
        # adjust the in and out distances for reaches in a direct downstream path from here
        out1 <- vertex.dist
        in1 <- vertex.dist - bb$Length
        
        # discharge.seq <- D[indx.rid,, drop=TRUE] == 0
        # 
        # in1 <- pmax(vertex.dist - edgedata$Length, 0)
        # in1[discharge.seq] <- D[, indx.rid, drop=TRUE][discharge.seq]
        
        if (nchar(bb$binaryID[bb$rid == rid]) > 1) {
          for (k in 1:(nchar(bb$binaryID[bb$rid == rid])-1)) {
            in1[bb$binaryID == substr(bb$binaryID[bb$rid == rid], 1, k)] <- vertex.dist[bb$binaryID == substr(bb$binaryID[bb$rid == rid], 1, k+1)]
          }
        }
        
        in1[upstr.indx] <- in1[upstr.indx] + d1
        in1[!upstr.indx] <- in1[!upstr.indx] - d1
        in1 <- c(in1, ppdata$upDist[[id]] - bb$Length[bb$binaryID == '1'])
        in1 <- pmax(in1, 0)
        
        out1[upstr.indx] <- out1[upstr.indx] + d1
        out1[!upstr.indx] <- out1[!upstr.indx] - d1
        out1 <- c(out1, ppdata$upDist[[id]])
        
        
        # vertex.dist[upstr.indx] <- vertex.dist[upstr.indx] + d1
        # vertex.dist[!upstr.indx] <- vertex.dist[!upstr.indx] - d1
        # in1 <- pmax(vertex.dist - c(bb$Length, bb$Length[bb$binaryID=='1']), 0)
        # in1 <- pmax(vertex.dist - bb$Length, 0)
        
        # in1[downstr.indx] <- 0
        # out1 <- vertex.dist
        cbind(in1, out1) -> inout1
        
        # r <- seq(0, 100000, length.out = 2000)
        cnt.ends <- cbind(cnt.ends, sapply(as.vector(empiricDistanceLocal[, j]), function(rr) sum(inout1[,1] < rr & inout1[,2]>rr)))
      }
    }
  
  # }
  
  cnt.ends
}
# The first row of cnt.ends should be all zeroes.
# cnt.ends[1,] <- rep(2, ncol(cnt.ends))
