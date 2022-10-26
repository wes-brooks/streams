linearKinhom2 <- function (X, net, lambda, r = NULL, ..., correction = "Ang", normalise = TRUE, lookup=NULL, ssn=NULL)  {

    # X <- obj$dwpr$X
    # stopifnot(inherits(X, "lpp"))
    
    correction <- pickoption("correction", correction, c(none = "none", Ang = "Ang", best = "Ang"), multi = FALSE)
    
    if (is.null(lambda)) 
        linearK(X, r = r, ..., correction = correction)

    lengthL <- sum(lookup$Length)
    lambdaX <- head(lambda, npoints(X))
    # lambdaX <- getlambda2.lpp(lambda, X, ...)
    invlam <- 1/lambdaX
    invlam2 <- outer(invlam, invlam, "*")
    
    denom <- if (!normalise) lengthL
    else sum(invlam)
    
    K <- linearKengine2(X, net = net, lookup=lookup, ssn=ssn, ..., r = r, reweight = invlam2, denom = denom, correction = correction)
    
    switch(correction, Ang = {
        ylab <- quote(K[L, inhom](r))
        fname <- c("K", "list(L, inhom)")
    }, none = {
        ylab <- quote(K[net, inhom](r))
        fname <- c("K", "list(net, inhom)")
    })
    
    K <- rebadge.fv(K, new.fname = fname, new.ylab = ylab)
    
    attr(K, "dangerous") <- attr(lambdaX, "dangerous")
    return(K)
}
