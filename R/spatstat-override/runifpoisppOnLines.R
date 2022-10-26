runifpoisppOnLines <- function (lambda, L, nsim = 1) {
    if (!is.numeric(lambda) || !all(is.finite(lambda) && (lambda >= 0))) 
        stop("lambda should be a finite, nonnegative number or numbers")
    if (!is.psp(L)) 
        L <- as.psp(L)
    W <- as.owin(L)
    result <- vector(mode = "list", length = nsim)
    for (i in 1:nsim) {
        X <- datagen.runifpoisppOnLines(lambda, L)
        Y <- ppp(X$x, X$y, marks = X$marks, window = W, check = FALSE)
        result[[i]] <- Y
    }
    if (nsim == 1) 
        return(result[[1]])
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
}