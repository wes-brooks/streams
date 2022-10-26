pointsOnLines2 <- function (L, eps = NULL, np = 1000, shortok = TRUE) {
    # stopifnot(is.psp(X))
    # X <- as.psp(L)
    # len <- lengths.psp(X)
    len <- L$lines$Length
    nseg <- length(len)
    if (is.null(eps)) {
        stopifnot(is.numeric(np) && length(np) == 1)
        stopifnot(is.finite(np) && np > 0)
        eps <- sum(len)/np
    }
    else {
        stopifnot(is.numeric(eps) && length(eps) == 1)
        stopifnot(is.finite(eps) && eps > 0)
    }
    Xdf <- as.data.frame(L$lines$ends)
    xmid <- with(Xdf, (x0 + x1)/2)
    ymid <- with(Xdf, (y0 + y1)/2)
    if (any(short <- (len <= eps)) && shortok) {
        Z <- data.frame(x = xmid[short], y = ymid[short])
    }
    else Z <- data.frame(x = numeric(0), y = numeric(0))
    for (i in (1:nseg)[!short]) {
        leni <- len[i]
        nwhole <- floor(leni/eps)
        if (leni/eps - nwhole < 0.5 && nwhole > 2) 
            nwhole <- nwhole - 1
        rump <- (leni - nwhole * eps)/2
        brks <- c(0, rump + (0:nwhole) * eps, leni)
        nbrks <- length(brks)
        ss <- (brks[-1] + brks[-nbrks])/2
        x <- with(Xdf, x0[i] + (ss/leni) * (x1[i] - x0[i]))
        y <- with(Xdf, y0[i] + (ss/leni) * (y1[i] - y0[i]))
        Z <- rbind(Z, data.frame(x = x, y = y))
    }
    Z <- lpp(as.ppp(Z, W = L$window), L=L)
    return(Z)
}