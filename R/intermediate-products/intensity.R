source('code/kernel.R')
source('code/bisquare.R')

# estimate sampling intensity via kernel weights:
# make a copy of the maris object to hold the sampling-intensity object
ssn <- maris
bws <- c(1e3, 2e3, 5e3, 10e3, 20e3, 50e3) # kernel bandwidths in meters


################################
sampling.point.data <- ssn@obspoints@SSNPoints[[1]]@point.data %>% group_by(siteID) %>% distinct(.keep_all=TRUE) %>% as.data.frame


# kernel weight obspoints
# calculate the sampling intensity at each observation point:
cat('Calculating kernel weights for observation locations:\n')
kernel.wt.pt <- matrix(0, nrow(maris@obspoints@SSNPoints[[1]]@point.data), length(bws))
for (netID in unique(maris@obspoints@SSNPoints[[1]]@point.data$netID)[-c(268, 314, 321, 323)]) {
  # for (netID in nnn) {
  cat("net:", netID, "\n")
  indx.net <- maris@obspoints@SSNPoints[[1]]@point.data$netID == netID
  kernel.wt.pt[indx.net,] <- kern(netID, sampling.point.data, maris@obspoints@SSNPoints[[1]]@point.data, edges, bws, nzw.rids=obs.rids)
}
save(kernel.wt.pt, file='data/intermediate-products/kernel.wt.pt.rdata')
cat('Done!\n')



# calculate the sampling intensity at each quadrature knot:
cat('Calculating kernel weights for quadrature points:\n')
kernel.wt.q <- matrix(0, nrow(maris@predpoints@SSNPoints[[1]]@point.data), length(bws))
for (netID in unique(maris@predpoints@SSNPoints[[1]]@point.data$netID)[-c(268, 314, 321, 323)]) {
  cat("net:", netID, "\n")
  indx.net <- maris@predpoints@SSNPoints[[1]]@point.data$netID == netID
  kernel.wt.q[indx.net,] <- kern(netID, sampling.point.data, maris@predpoints@SSNPoints[[1]]@point.data, edges, bws, nzw.rids=obs.rids)
}
save(kernel.wt.q, file='data/intermediate-products/kernel.wt.q.rdata')
cat('Done!\n') 



