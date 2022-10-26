# depict the basis function algorithm
load('data/basis-function-figure-data.rdata')

# work on just the shoshone and green rivers
grn <- net
grn$lines$marks[!(grn$lines %in% c())] <- NA
indx <- !is.na(grn$lines$marks)

# strip out unnecessary entries in the network data
grn$lines$ends <- grn$lines$ends[indx,]
grn$lines$n <- sum(indx)
grn$lines$marks <- grn$lines$marks[indx]
grn$lines$rid <- grn$lines$rid[indx]
grn$lines$Length <- grn$lines$Length[indx]
grn$lines$known <- grn$lines$known[indx]

# adjust the window extents
grn$lines$window$xrange <- range(grn$lines$ends[,c('x0', 'x1')])
grn$lines$window$yrange <- range(grn$lines$ends[,c('y0', 'y1')])


# plotting setup
layout(matrix(1:6, 2, 3, byrow=TRUE))
par(oma=c(0,0,0,0), mar=c(0,0,0,0))

# plot the decomposition
for (i in 1:6) {
    grn$lines$marks <- basefuns[, i, drop=TRUE]
    plot(grn, main='', ribn=7)
    
    for (br in sort(unique(basis$br))[1:(i+1)]) {
        indx <- basis$br == br
        indx <- rid.lookup$seg[rid.lookup$rid == basis$rid[indx][which.min(basis$binlen[indx])]]
        
        points(net$lines$ends[indx, c(1,2)], pch=19, cex=1.5)
    }
}

# Add a legend
legend(x='topleft', title="subnetwork", legend=1:7, bty='n', col=rainbow(7)[c(1,7,2,3,4,5,6)], pch=15, cex=1.5, pt.cex=3, title.adj=0)

