linearK2 <- function (X, net, r = NULL, ..., D = NULL, correction = "Ang",  ssn=NULL, lookup=NULL)  {
    #X <- obj$dwpr$X
    
    #stopifnot(inherits(X, "lpp"))
    correction <- pickoption("correction", correction, c(none = "none", 
                                                         Ang = "Ang", best = "Ang"), multi = FALSE)
    sX <- summary(X)
    np <- sX$npoints
    lengthL <- sX$totlength
    denom <- np * (np - 1)/lengthL
    K <- linearKengine2(X, net = net, r = r, denom = denom, correction = correction,  ssn=ssn, lookup=lookup, ...)
    switch(correction, Ang = {
        ylab <- quote(K[L](r))
        fname <- c("K", "L")
    }, none = {
        ylab <- quote(K[net](r))
        fname <- c("K", "net")
    })
    K <- rebadge.fv(K, new.ylab = ylab, new.fname = fname)
    return(K)
}
