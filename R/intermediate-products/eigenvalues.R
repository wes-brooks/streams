# compute the eigenvalues of the precision matrix for each network:

load('data/intermediate-products/maris-ssn.rdata')
load('data/intermediate-products/stream-network-structure.rdata')

obs.nets <- as.numeric(as.character(unique(maris@obspoints@SSNPoints[[1]]@point.data$netID)))


#########################################################################################
# establish the order of reaches in the GMRF specification:
indx.valid.net <- edges@data$netID %in% obs.nets
gmrf.lookup <- edges@data[indx.valid.net, c('netID', 'rid')]
gmrf.lookup$gmrfID <- 1:nrow(gmrf.lookup)


########################################################################################
# recall which eigenvalues are connected with which networks (eigenvalues are sorted by network, not in the same order as gmrf.lookup)
eigenvalues <- vector()
for (nn in obs.nets) {
  indx <- gmrf.lookup$netID == nn
  submat <- Sigma.mask[indx, indx]
  eigenvalues <- c(eigenvalues, eigen(submat)$values)
}

save(eigenvalues, file='data/intermediate-products/eigenvalues.rdata')