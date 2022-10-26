npts <- vector()

for (filename in dir(paste0(maris@path, "/distance/obs/"))) {
  pth <- paste0(maris@path, "/distance/obs/", filename)
  ff <- file(pth, open='rb')
  D <- unserialize(ff)
  close(ff)
  
  npts <- c(npts, nrow(D))
}


nets <- vector()
for (filename in dir(paste0(maris@path, "/distance/obs/"))) {
  indx <- gregexpr("[0-9+]", filename)[[1]]
  if (length(indx > 1)) {
    nets <- c(nets, substr(filename, indx[1], tail(indx, 1)))
  } else {
    nets <- c(nets, substr(filename, indx[1], indx[1]))
  }
}



nvert <- vector()
npt <- vector()
for (n in nets) {
  nvert <- c(nvert, sum(maris@network.line.coords$NetworkID == n))
  npt <- c(npt, sum(maris@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID == n))
}