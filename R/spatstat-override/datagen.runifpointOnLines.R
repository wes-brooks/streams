datagen.runifpointOnLines <- function (n, L, ..., nzw.rids=NULL, correction.factor=NULL) {
    m <- length(n)
    ismarked <- (m > 1)
    if (m == 0 || (m == 1 && n == 0)) 
        return(data.frame(x = numeric(0), y = numeric(0), seg = integer(0), 
                          tp = numeric(0)))
    
    # if we have specified the only some lines have non-zero weight then use only those
    if (!is.null(nzw.rids)) {
      indx <- L$lines$rid %in% nzw.rids
    } else indx <- rep(TRUE, length(L$lines$Length))
    
    # specify a correction factor because some lines have nonzero intensity only on part of the lines.
    if (is.null(correction.factor)) {
      correction.factor <- cbind(rid=L$lines$rid, scale=rep(1, length(L$lines$Length)))
    }
    scale.indx <- match(L$lines$rid, correction.factor$rid)
    
    #
    len <- L$lines$Length[indx] * correction.factor$scale[scale.indx][indx]
    sumlen <- sum(len)
    cumlen <- cumsum(len)
    cum0len <- c(0, cumlen)
    Ldf <- as.data.frame(as.psp(L))[indx,]
    x0 <- with(Ldf, x0)
    y0 <- with(Ldf, y0)
    dx <- with(Ldf, x1 - x0)
    dy <- with(Ldf, y1 - y0)
    if (ismarked) {
        markvalues <- names(n)
        if (sum(nzchar(markvalues)) < m) 
            markvalues <- paste(1:m)
    }
    out <- data.frame(x = numeric(0), y = numeric(0), seg = integer(0), 
                      tp = numeric(0))
    if (ismarked) 
        out <- cbind(out, data.frame(marks = character(0)))
    for (j in 1:m) {
        if (n[[j]] > 0) {
            uu <- runif(n[[j]], min = 0, max = sumlen)
            kk <- findInterval(uu, cum0len, rightmost.closed = TRUE, 
                               all.inside = TRUE)
            tt <- (uu - cum0len[kk])/len[kk]
            tt[!is.finite(tt)] <- 0
            x <- x0[kk] + tt * dx[kk]
            y <- y0[kk] + tt * dy[kk]
            if (!ismarked) {
                out <- data.frame(x = x, y = y, seg = L$lines$seg[indx][kk], tp = tt)
            }
            else {
                outj <- data.frame(x = x, y = y, seg = L$lines$seg[indx][kk], tp = tt, 
                                   marks = markvalues[j])
                out <- rbind(out, outj)
            }
        }
    }
    if (ismarked) 
        out$marks <- factor(out$marks, levels = markvalues)
    return(out)
}