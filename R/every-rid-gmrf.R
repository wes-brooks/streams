# kill any presence points that are on a rid with no quadrature points
obs.rids <- unique(maris@obspoints@SSNPoints[[1]]@point.data$rid)
rid.noquad <- obs.rids[!(obs.rids %in% maris@predpoints@SSNPoints[[1]]@point.data$rid)]
indx.noquad <-  maris@obspoints@SSNPoints[[1]]@point.data$rid %in% rid.noquad
maris@obspoints@SSNPoints[[1]]@point.data <- maris@obspoints@SSNPoints[[1]]@point.data[!indx.noquad,]
obs.rids <- unique(maris@obspoints@SSNPoints[[1]]@point.data$rid)

# New GMRF precision matrix based on individual stream segments:
# First, let's consider only rids on networks with at least one point.
obs.nets <- as.numeric(as.character(unique(maris@obspoints@SSNPoints[[1]]@point.data$netID)))
indx.valid.net <- edges@data$netID %in% obs.nets
gmrf.lookup <- edges@data[indx.valid.net, c('netID', 'rid')]
gmrf.lookup$gmrfID <- 1:nrow(gmrf.lookup)





# identify the downstream reach for each reach
downstr <- vector()
for (i in 1:nrow(edges@data)) {
  cat(i, '\n')
  b.ref <- edges@data$binaryID[i]
  len.ref <- nchar(b.ref)
  indx.net <- edges@data$netID == edges@data$netID[i]
  
  edges.net <- edges@data[indx.net,]
  downstr.indx <- edges.net$binaryID == substr(b.ref, 1, len.ref - 1)
  
  if (any(downstr.indx)) {
    downstr <- c(downstr, edges.net$rid[downstr.indx][1])
  } else downstr <- c(downstr, NA)
}
edges@data$downstr <- downstr





# new method for figuring sampling intensity:
# observation sites / length
# indx <- match(XX$rid, edges@data$rid)
indx.q <- match(maris@predpoints@SSNPoints[[1]]@point.data$rid, edges@data$rid)
indx.cox <- match(maris@obspoints@SSNPoints[[1]]@point.data$rid, edges@data$rid)


obslocs <- maris@obspoints@SSNPoints[[1]]@point.data %>% group_by(rid) %>% summarize(cnt=n())
indx.obs <- match(edges@data$rid, obslocs$rid)
obslocs$cnt[indx.obs] -> edges@data$cnt
edges@data$cnt[is.na(edges@data$cnt)] <- 0
edges@data$effort <- edges@data$cnt / edges@data$Length

# XX$offset <- log(edges@data$effort[indx])


offset.cox <- log(edges@data$effort[indx.cox])
offset.q <- log(edges@data$effort[indx.q])





# define the random effects design matrix and the precision matrix
Z <- Matrix(0, nrow(XX), sum(indx.valid.net))
Sigma.mask <- Matrix(0, sum(indx.valid.net), sum(indx.valid.net))

# populate the random effects design matrix and the precison matrix of the GMRF
for (i in 1:nrow(edges@data)) {
  colID <- gmrf.lookup$rid == edges@data$rid[i]
  rowID <- XX$rid == edges@data$rid[i]

  if (edges@data$netID[i] %in% obs.nets && !is.na(edges@data$downstr[i])) {
    downstrID <- gmrf.lookup$rid == edges@data$downstr[i]
    Sigma.mask[downstrID, colID] <- Sigma.mask[colID, downstrID] <- -1
    # indx.match <- edges@data$rid == edges@data$downstr[downstrID]
    # if (any(indx.match)) {
    #   Sigma.mask[i, edges@data$rid[indx.match] + 1] <- Sigma.mask[edges@data$rid[indx.match] + 1, i] <- -1
    # }
  }
}
diag(Sigma.mask) <- -colSums(Sigma.mask)



effort.fun.kernel <- function(x, y, seg, tp, data, lookup, kernel.function, coef=1, raw=FALSE) {
  
  out <- kernel.function(x, y, seg, tp)
  
  if (!raw)
    out <- exp(coef * log(out))
  
  out
}


effort.wrapper.kernel <- function(data, lookup, coef, kernel.function) {
  wrapped <- function(x, y, seg, tp, raw=FALSE) {
    effort.fun.kernel(x, y, seg, tp, data, lookup, kernel.function, coef, raw)
  }
  
  wrapped
}


effort.fun <- function(x, y, seg, tp, data, lookup, coef=1, raw=FALSE) {
  indx <- match(seg, lookup$seg)
  indx <- match(lookup$rid[indx], data$rid)
  
  if (raw) {
    out <- ifelse(!is.na(indx), data$effort[indx], 0)
  } else {
    out <- ifelse(!is.na(indx), exp(coef * log(data$effort[indx])), 0)
  }
  
  out
}



effort.wrapper <- function(func, data, lookup, coef) {
  wrapped <- function(x, y, seg, tp, raw=FALSE) {
    func(x, y, seg, tp, data, lookup, coef, raw)
  }
  
  return(wrapped)
}