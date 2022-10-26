# getting the number of active edges at a given distance from the observation point:
# lpp <- bt.lpp


kern <- function(netnum, obspts, quadpts, edges, bw, nzw.rids=NULL, distpath="~/Dropbox/streams/data/maris_ssn/maris4.ssn/distance/obs/") {
  i <- 1

  # retrieve the serialized network of distances
  f <- file(paste0(distpath, "dist.net.new.", netnum, ".RData"), open='rb')
  D <- unserialize(f)
  close(f)
  
  D <- D$obs[[1]]
  
  # subset the data objects to just focus on this network
  edgedata <- edges[edges$netID == netnum,]
  quaddata <- as.data.frame(quadpts[quadpts$netID == netnum,])
  obsdata <- obspts[obspts$netID == netnum,]
  
  indx <- match(quaddata$rid, edgedata$rid)
  quaddata$upDist <- edgedata$upDist[indx] - (1 - quaddata$ratio) * edgedata$Length[indx]
  
  rid.obs.match <- match(obsdata$rid, edgedata$rid)
  obsdata$Length <- edgedata$Length[rid.obs.match]
  
  kernel.wt <- matrix(0, nrow(quaddata), length(bw))
  
  # calculate the distances along the stream network
  for (rid in unique(quaddata$rid)) {
    if (any(edgedata$rid == rid) && is.null(nzw.rids) || rid %in% nzw.rids) {
      # consider simultaneously all the points that are on the same stream reach
      indx.rid <- which(edgedata$rid == rid)[1]
      vertex.dist <- D[indx.rid,] + D[,indx.rid]
      
      # 
      upstr.indx <- substr(c(edgedata$binaryID, '0'), 1, nchar(edgedata$binaryID[indx.rid])) == edgedata$binaryID[indx.rid]
      in1 <- pmax(vertex.dist - edgedata$Length, 0)
      
      # mark which reaches are close enough to potentially have observations with nonzero weight
      indx.nonzero.weight <- in1 < max(bw) + edgedata$Length[edgedata$rid == rid]
      nzw.edgedata <- edgedata[indx.nonzero.weight,]
      nzw.obsdata <- obsdata[obsdata$rid %in% nzw.edgedata$rid,]
      indx.rid.quad <- quaddata$rid == rid
      rid.quaddata <- quaddata[indx.rid.quad,]
      
      # quit early if no observation points can have weight
      if (nrow(nzw.obsdata) == 0) {
        rid.kernel.wt <- matrix(0, sum(rid.quaddata$rid == rid), length(bw))
      } else {
        # find the distance from each potentially nonzero observation to the vertex at the top of this reach.
        nzw.vertex.dist <- vertex.dist[indx.nonzero.weight]
        indx.obs.upstr <- nzw.obsdata$rid %in% edgedata$rid[upstr.indx]
        indx.obs.downstr <- (!indx.obs.upstr) & (nzw.obsdata$rid != rid)
        
        ridbin <- edgedata$binaryID[indx.rid]
        indx.discharge <- as.vector(sapply(edgedata$binaryID, function(bin) nchar(bin) <= nchar(ridbin) && substr(ridbin, 1, nchar(bin)) == bin))
        indx.obs.discharge <- nzw.obsdata$rid %in% edgedata$rid[indx.discharge]
        indx.obs.rid <- nzw.obsdata$rid == rid
        
        adjustment <- rep(0, nrow(nzw.obsdata))
        adjustment[!indx.obs.discharge] <- -(1 - nzw.obsdata$ratio[!indx.obs.discharge]) * nzw.obsdata$Length[!indx.obs.discharge]
        adjustment[indx.obs.discharge] <- (1 - nzw.obsdata$ratio[indx.obs.discharge]) * nzw.obsdata$Length[indx.obs.discharge]
        
        # get the distance of each observation point from the vertex of this rid
        pt.match <- match(nzw.obsdata$rid, nzw.edgedata$rid)
        obs.dist <- nzw.vertex.dist[pt.match] + adjustment
        
        first.pt <- TRUE
        nzw.indx.rid <- which(nzw.edgedata$rid == rid)[1]
        rid.quads <- which(rid.quaddata$rid == rid)
        rid.kernel.wt <- matrix(0, 0, length(bw))
        for (j in 1:length(rid.quads)) {
          id <- rid.quads[j]
          # cat(i, "\n")
          i <- i+1 
          
          if (first.pt) {
            d1 <- nzw.edgedata$upDist[nzw.indx.rid] - rid.quaddata$upDist[[id]]
          } else {
            d1 <- uD.prev - rid.quaddata$upDist[[id]]
          }
          uD.prev <- rid.quaddata$upDist[[id]]
          first.pt <- FALSE
          
          obs.dist[indx.obs.upstr] <- obs.dist[indx.obs.upstr] + d1
          obs.dist[indx.obs.downstr] <- obs.dist[indx.obs.downstr] - d1
          obs.dist[indx.obs.rid] <- abs(rid.quaddata$upDist[[id]] - nzw.obsdata$upDist[indx.obs.rid])
          
          rid.kernel.wt <- rbind(rid.kernel.wt, colSums(bisquare(obs.dist, bw)) / bw)
        }
        # cat(rid.kernel.wt, '\n')
      }
    
      kernel.wt[indx.rid.quad,] <- rid.kernel.wt
      cat(sum(kernel.wt != 0), rid, '\n')
    }
  }
  
  kernel.wt
}



kernel.wrap <- function(kernel.function, lookup, obspts, edges, bw, nzw.rids=NULL, distpath="~/Dropbox/streams/data/maris_ssn/maris4.ssn/distance/obs/", logresult=FALSE) {
  
  fun <- function(x, y, seg, tp) {
    # prepare a result object
    out <- matrix(NA, nrow=length(seg), ncol=length(bw))
    seg.indx <- match(seg, lookup$seg)
    nets <- unique(lookup$netID[seg.indx])
    
    # convert the (x, y, seg, tp) to centerpoints as expected by the kernel function:
    centerpoints <- data.frame(rid=lookup$rid[seg.indx], ratio=1-tp, netID=lookup$netID[seg.indx])
    
    
    # loop through the networks
    for (nn in nets) {
      indx.net <- lookup$netID[seg.indx] == nn
      
      # calculate the kernel density on netID nn
      out[indx.net,] <- kernel.function(netnum=nn, obspts=obspts, quadpts=centerpoints, edges=edges, bw=bw, nzw.rids=nzw.rids, distpath=distpath)
    }
    
    # if there's only one bandwidth, then return a vector instead of a matrix.
    if (length(bw) == 1)
      out <- as.vector(out)
    
    if (logresult)
      out <- log(out)
    
    # return the result
    out
  }
}
