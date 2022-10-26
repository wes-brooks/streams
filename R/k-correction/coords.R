# Need network.point.coords for every junction
nn <- 20
# for (nn in sort(unique(basis$netID))) {
for (nn in missednets) {
  predpts <- maris@predpoints
  
  pred.indx <- edges@data$netID == nn
  
  
  predpt2 <- list()
  
  point.coords <- t(sapply(which(pred.indx), function(i) {tmp <- edges@lines[i][[1]]@Lines[[1]]@coords; tmp[1,]}))
  colnames(point.coords) <- c('coords.x1', 'coords.x2')
  
  point.data <- edges@data[pred.indx,]
  indx <- match(point.data$rid, edges@data$rid)
  network.point.coords <- data.frame(NetworkID=nn, SegmentID=as.factor(edges@data$rid[indx]), DistanceUpstream=edges@data$upDist[indx])
  attr(network.point.coords, 'locID') <- 1:nrow(network.point.coords)
  
  bounds <- t(apply(point.coords, 2, range))
  colnames(bounds) <- c('min', 'max')
  vvv <- new('SSNPoint',
             network.point.coords=network.point.coords,
             point.coords=point.coords,
             point.data=point.data,
             points.bbox=bounds,
             proj4string=predpts@SSNPoints[[1]]@proj4string)
  
  
  
  
  m2@network.line.coords <- network.point.coords
  m2@obspoints@SSNPoints[[1]] <- vvv

  
  time <- system.time(myDmat <- createDistMatInMemory(m2, binTable=bincodes))
  cat(nn, ":", time, "\n")
  
  filehandle <- file(file.path(maris@path, "distance", "obs", paste0('dist.net.new.', nn, '.RData')), open='wb')
  serialize(myDmat, filehandle, ascii=FALSE)
  close(filehandle)
}
