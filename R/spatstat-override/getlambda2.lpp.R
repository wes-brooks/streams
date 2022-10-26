getlambda2.lpp <- function (lambda, X, ..., update = TRUE) {
    lambdaname <- deparse(substitute(lambda))
    #XX <- as.ppp(X)
    XX <- X
    missup <- missing(update)
    if (update && (is.lppm(lambda) || is.ppm(lambda))) {
        danger <- FALSE
        lambda <- if (is.lppm(lambda)) 
            update(lambda, X)
        else update(lambda, XX)
        if (missup) 
            warn.once("lin.inhom.update", "The behaviour of linearKinhom and similar functions", 
                      "when lambda is an lppm object", "has changed (in spatstat 1.41-0 and later)", 
                      "See help(linearKinhom)")
    }
    else danger <- TRUE
    lambdaX <- if (is.vector(lambda)) 
        lambda
    else if (is.function(lambda)) 
        lambda(x=XX$data$x, y=XX$data$y, seg=XX$data$seg, tp=XX$data$tp, ...)
    else if (is.im(lambda)) 
        safelookup(lambda, XX)
    else if (inherits(lambda, "linim")) 
        safelookup(as.im(lambda), XX)
    else if (is.ppm(lambda) || inherits(lambda, "lppm")) 
        predict(lambda, locations = as.data.frame(XX))
    else stop(paste(lambdaname, "should be", "a numeric vector, function, pixel image, or fitted model"))
    if (!is.numeric(lambdaX)) 
        stop(paste("Values of", lambdaname, "are not numeric"))
    if ((nv <- length(lambdaX)) != (np <- npoints(X))) 
        stop(paste("Obtained", nv, "values of", lambdaname, "but point pattern contains", 
                   np, "points"))
    if (any(lambdaX < 0)) 
        stop(paste("Negative values of", lambdaname, "obtained"))
    if (any(lambdaX == 0)) 
        stop(paste("Zero values of", lambdaname, "obtained"))
    if (danger) 
        attr(lambdaX, "dangerous") <- lambdaname
    return(lambdaX)
}