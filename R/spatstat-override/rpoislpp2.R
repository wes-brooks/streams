rpoislpp2 <- function (lambda, L, ..., nsim = 1, nzw.rids=NULL, correction.factor=NULL, loglambda=FALSE) {
    if (missing(L) || is.null(L)) {
        if (!inherits(lambda, c("linim", "linfun"))) 
            stop("L is missing", call. = FALSE)
        L <- as.linnet(lambda)
    }
    else verifyclass(L, "linnet")
    result <- vector(mode = "list", length = nsim)
    for (i in 1:nsim) {
        X <- datagen.rpoisppOnLines2(lambda, L, nzw.rids=nzw.rids,correction.factor=correction.factor, loglambda=loglambda, ...)
        Y <- lpp(X, L)
        if (nsim == 1) 
            return(Y)
        result[[i]] <- Y
    }
    Y <- as.solist(Y)
    names(Y) <- paste("Simulation", 1:nsim)
    return(Y)
}