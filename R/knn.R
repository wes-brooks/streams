# getting the number of active edges at a given distance from the observation point:
# lpp <- bt.lpp


knn <- function(netnum, obspts, quadpts, edges, k=1) {
  i <- 1
  
  # retrieve the serialized network of distances
  f <- file(paste0("~/Dropbox/streams/data/maris_ssn/maris4.ssn/distance/obs/dist.net.new.", netnum, ".RData"), open='rb')
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
  
  nn.dist <- matrix(Inf, nrow(quaddata), length(k))
  
  # calculate the distances along the stream network
  for (rid in unique(quaddata$rid)) {
    # consider simultaneously all the points that are on the same stream reach
    if (any(edgedata$rid == rid)) {
      indx.rid <- which(edgedata$rid == rid)[1]
      vertex.dist <- D[indx.rid,] + D[,indx.rid]
      
      # 
      upstr.indx <- substr(c(edgedata$binaryID, '0'), 1, nchar(edgedata$binaryID[indx.rid])) == edgedata$binaryID[indx.rid]
      discharge.seq <- rev(sapply(1:nchar(edgedata$binaryID[indx.rid]), function(n) which(edgedata$binaryID == substr(edgedata$binaryID[indx.rid], 1, n))[1]))
      indx.discharge <- (1:nrow(edgedata@data)) %in% discharge.seq
      n.discharge <- length(discharge.seq)
      
      in1 <- vertex.dist
      in1[!indx.discharge] <- pmax(vertex.dist[!indx.discharge] - edgedata$Length[!indx.discharge], 0)
      out1 <- vertex.dist
      if (n.discharge > 1) {
        out1[discharge.seq] <- D[, indx.rid, drop=TRUE][discharge.seq] + edgedata$Length[discharge.seq]
      }
      
      # this is the table that we'll use to search for the nearest neighbor(s)
      nn.table <- data.frame(rid=edgedata@data$rid, d.in=in1, d.out=out1, npts=sapply(edgedata@data$rid, function(r) {sum(obsdata$rid==r)}))
      
      # mark which reaches are close enough to potentially have observations with nonzero weight
      indx.nonzero.weight <- nn.table$d.in < min(nn.table$d.out[nn.table$npts > 0]) + edgedata$Length[edgedata$rid == rid]
      nzw.edgedata <- edgedata[indx.nonzero.weight,]
      nzw.obsdata <- obsdata[obsdata$rid %in% nzw.edgedata$rid,]
      indx.rid.quad <- quaddata$rid == rid
      rid.quaddata <- quaddata[indx.rid.quad,]
      # nzw.nn.table <- nn.table[indx.nonzero.weight,]
      
      # quit early if no observation points are reachable 
      if (sum(nn.table$npts) < min(k)) {
        rid.nn.dist <- matrix(Inf, sum(rid.quaddata$rid == rid), length(k))
      } else {
        # find the distance from each potentially nonzero observation to the vertex at the top of this reach.
        nzw.vertex.dist <- vertex.dist[indx.nonzero.weight]
        indx.obs.upstr <- (nzw.obsdata$rid %in% edgedata$rid[upstr.indx]) & (nzw.obsdata$rid != rid)
        indx.obs.downstr <- (!indx.obs.upstr) & (nzw.obsdata$rid != rid)
        indx.obs.discharge <- nzw.obsdata$rid %in% edgedata$rid[indx.discharge]
        indx.obs.rid <- nzw.obsdata$rid == rid
        
        #
        adjustment <- rep(0, nrow(nzw.obsdata))
        adjustment[!indx.obs.discharge] <- -(1 - nzw.obsdata$ratio[!indx.obs.discharge]) * nzw.obsdata$Length[!indx.obs.discharge]
        adjustment[indx.obs.discharge] <- (1 - nzw.obsdata$ratio[indx.obs.discharge]) * nzw.obsdata$Length[indx.obs.discharge]
        
        # 
        # indx.edge.upstr <- (nzw.edgedata$rid %in% edgedata$rid[upstr.indx]) & (nzw.edgedata$rid != rid)
        # indx.edge.downstr <- (!indx.edge.upstr) & (nzw.edgedata$rid != rid)
        # indx.edge.rid <- nzw.edgedata$rid == rid
        
        # get the distance of each observation point from the vertex of this rid
        pt.match <- match(nzw.obsdata$rid, nzw.edgedata$rid)
        obs.dist <- nzw.vertex.dist[pt.match] + adjustment
        # nzw.nn.table$nearest <- sapply(1:nrow(nzw.edgedata), function(r) ifelse(r %in% pt.match, min(obs.dist[pt.match == r]), Inf))
        
        first.pt <- TRUE
        nzw.indx.rid <- which(nzw.edgedata$rid == rid)[1]
        rid.quads <- which(rid.quaddata$rid == rid)
        rid.nn.dist <- matrix(Inf, 0, length(k))
        for (j in 1:length(rid.quads)) {
          id <- rid.quads[j]
          
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
          
          
          rid.nn.dist <- rbind(rid.nn.dist, sort(obs.dist)[k])
        }
      }
      
      nn.dist[indx.rid.quad] <- rid.nn.dist
      cat(sum(nn.dist != Inf, na.rm=TRUE), rid, '\n')
    }
  }
  
  nn.dist
}
