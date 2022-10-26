ppllengine2 <- function (X, Y, action = "project", check = FALSE)  {
    # stopifnot(is.ppp(X))
    stopifnot(is.psp(Y))
    stopifnot(action %in% c("distance", "identify", "project"))
    if (Y$n == 0) 
        stop("Segment pattern Y contains 0 segments; projection undefined")
    if (npoints(X) == 0) {
        nowt <- numeric(0)
        none <- integer(0)
        switch(action, identify = return(none), distance = return(list(dist = nowt, 
                                               which = none)), project = return(list(Xproj = X, 
                                                                                     mapXY = none, d = nowt, tp = nowt)))
    }
    XX <- as.matrix(as.data.frame(unmark(X)))
    YY <- as.matrix(as.data.frame(unmark(Y)))
    huge <- max(diameter(as.rectangle(X$domain$window)), diameter(as.rectangle(Y$window)))
    d <- distppllmin(XX, YY, huge^2)
    mapXY <- d$min.which
    if (action == "identify") 
        return(mapXY)
    else if (action == "distance") 
        return(data.frame(dist = d$min.d, which = mapXY))
    alldata <- as.data.frame(cbind(XX, YY[mapXY, , drop = FALSE]))
    colnames(alldata) <- c("x", "y", "x0", "y0", "x1", "y1")
    dx <- with(alldata, x1 - x0)
    dy <- with(alldata, y1 - y0)
    leng <- sqrt(dx^2 + dy^2)
    co <- dx/leng
    si <- dy/leng
    xv <- with(alldata, x - x0)
    yv <- with(alldata, y - y0)
    xpr <- xv * co + yv * si
    ok <- is.finite(xpr)
    left <- !ok | (xpr <= 0)
    right <- ok & (xpr >= leng)
    xr <- with(alldata, ifelseAX(left, 0, ifelseXY(right, leng, 
                                                   xpr)))
    xproj <- with(alldata, x0 + ifelseXB(ok, xr * co, 0))
    yproj <- with(alldata, y0 + ifelseXB(ok, xr * si, 0))
    Xproj <- ppp(xproj, yproj, window = X$domain$window, marks = X$marks, 
                 check = check)
    tp <- xr/leng
    tp[!is.finite(tp)] <- 0
    return(list(Xproj = Xproj, mapXY = mapXY, d = d$min.d, tp = tp))
}