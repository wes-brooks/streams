handle.r.b.args2 <- function (r = NULL, breaks = NULL, window, pixeps = NULL, rmaxdefault = NULL) 
{
    if (!is.null(r) && !is.null(breaks)) 
        stop(paste("Do not specify both", sQuote("r"), "and", 
                   sQuote("breaks")))
    if (!is.null(breaks)) {
        breaks <- as.breakpts(breaks)
    }
    else if (!is.null(r)) {
        breaks <- breakpts.from.r(r)
    }
    else {
        rmax <- rmaxdefault
        if (is.null(rmax)) rmax <- diameter(Frame(window))
        if (is.null(pixeps)) {
            pixeps <- if (is.mask(window)) 
                min(window$xstep, window$ystep)
            else rmax/128
        }
        rstep <- pixeps/4
        breaks <- make.even.breaks(rmax, bstep = rstep)
    }
    return(breaks)
}