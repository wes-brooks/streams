
# get the intensity into the same data as the paths
indx <- match(edges@data$rid, XX$rid)
edges@data$lambda <- lambda[indx]


X.pred <- model.matrix(~ poly(logMS_Hist, 1) + poly(S1_93_11, 2) + poly(SLOPE, 1), data=XX)
eta.pred <- as.vector(X.pred %*% obj$cox$beta[1:5] + Z %*% obj$cox$M)

XX2 <- XX[!is.na(eta.pred),]

indx <- match(edges@data$rid, XX2$rid)
edges@data$eta.pred <- eta.pred[indx]

# extract the path data from the structure of SSN
paths <- sapply(edges@lines, function(x) as.data.frame(x@Lines[[1]]@coords), simplify=FALSE)
for (i in seq_len(length(paths))) {
  paths[[i]]$seg <- i
  paths[[i]]$lambda <- edges@data$lambda[i]
  paths[[i]]$eta.pred <- edges@data$eta.pred[i]
  paths[[i]]$rid <- edges@data$rid[i]
}

# convert paths from a list of data.frames into a data.frame
do.call(rbind.data.frame, paths) -> tmp
colnames(tmp) <- c('x', 'y', 'seg', 'lambda', 'eta.pred', 'rid')
colnames(tmp) <- c('x', 'y', 'seg', 'eta.pred', 'rid')

# filter out the smallest streams
indx <- match(tmp$rid, edges@data$rid)
tmp$major <- (edges@data$TotDASqKM > 15)[indx]

# plot the river network with coloration due to intensity
#ggplot(tmp) + aes(x=x, y=y, group=seg, color=eta.pred) + geom_path()
tmp2 <- tmp[tmp$major,]
tmp2$eta.pred[tmp2$eta.pred < -25] <- -25

#ggplot(tmp2) + aes(x=x, y=y, group=seg, color=eta.pred) + theme_minimal() + scale_color_gradient(na.value=NA, low='grey90', high='black') + geom_path()
ggplot(tmp2) + aes(x=x, y=y, group=seg, color=eta.pred) + theme_minimal()  + theme(axis.ticks=element_blank(), axis.text=element_blank()) +
    labs(colour=' brown trout\nlog-intensity', x='', y='') +
    scale_color_gradientn(na.value=NA, colours=c('grey60', 'tan', 'red')) + geom_path() 
