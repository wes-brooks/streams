dist2dpath2 <- function (dist, method = "C") 
{
    stopifnot(is.matrix(dist) && isSymmetric(dist))
    stopifnot(all(diag(dist) == 0))
    findist <- dist[is.finite(dist)]
    if (any(findist < 0)) 
        stop("Some distances are negative")
    n <- nrow(dist)
    if (n <= 1) 
        return(dist)
    cols <- col(dist)
    tol <- .Machine$double.eps
    posdist <- findist[findist > 0]
    if (length(posdist) > 0) {
        shortest <- min(posdist)
        tol2 <- shortest/max(n, 1024)
        tol <- max(tol, tol2)
    }
    switch(method, interpreted = {
        dpathnew <- dpath <- dist
        changed <- TRUE
        while (changed) {
            for (j in 1:n) dpathnew[, j] <- apply(dpath + dist[j, 
                                                               ][cols], 1, min)
            unequal <- (dpathnew != dpath)
            changed <- any(unequal) & any(abs(dpathnew - dpath)[unequal] > 
                                              tol)
            dpath <- dpathnew
        }
    }, C = {
        adj <- is.finite(dist)
        diag(adj) <- TRUE
        d <- dist
        d[!adj] <- -1
        z <- .C("Ddist2dpath", nv = as.integer(n), d = as.double(d), 
                adj = as.integer(adj), dpath = as.double(numeric(n * 
                                                                     n)), tol = as.double(tol), niter = as.integer(integer(1)), 
                status = as.integer(integer(1)))
        if (z$status == -1) warning(paste("C algorithm did not converge to tolerance", 
                                          tol, "after", z$niter, "iterations", "on", n, "vertices and", 
                                          sum(adj) - n, "edges"))
        dpath <- matrix(z$dpath, n, n)
        dpath[dpath < 0] <- Inf
    }, stop(paste("Unrecognised method", sQuote(method))))
    return(dpath)
}