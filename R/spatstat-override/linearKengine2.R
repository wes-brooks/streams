linearKengine2 <- function (X, ..., net = NULL, r = NULL, reweight = NULL, denom = 1, correction = "Ang", showworking = FALSE, ssn=NULL, lookup=NULL) {
    # X <- obj$dwpr$X
    
    np <- npoints(X)
    L <- net
    #Y <- as.ppp(X)
    Y <- X
    W <- net$window
    rmaxdefault <- 0.98 * circumradius(L$window)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max
    
    type <- if (correction == "Ang")  "L"
    else "net"
    
    fname <- c("K", type)
    ylab <- substitute(K[type](r), list(type = type))
    
    if (np < 2) {
        zeroes <- numeric(length(r))
        df <- data.frame(r = r, est = zeroes)
        K <- fv(df, "r", ylab, "est", . ~ r, c(0, rmax), c("r", 
                                                           makefvlabel(NULL, "hat", fname)), c("distance argument r", 
                                                                                               "estimated %s"), fname = fname)
        return(K)
    }
    D <- as.matrix(pairdist2.lpp(X, net, ssn, lookup))
    if (correction == "none" && is.null(reweight)) {
        K <- compileK(D, r, denom = denom, fname = fname)
        K <- rebadge.fv(K, ylab, fname)
        unitname(K) <- unitname(X)
        return(K)
    }
    if (correction == "none") 
        edgewt <- 1
    else {
        m <- matrix(1, np, np)
        
        # if ("ssnlpp" %in% class(Y)) {
        #   m <- countends.ssnlpp(L, Y, D)
        # } else {
        #   for (j in 1:np) m[-j, j] <- countends2(L, Y[-j], D[-j, j])
        # }
        
        raw <- matrix(NA, nrow(D), 0)
        reorder <- vector()
        for (nn in unique(Y$data$netID)) {
          cat(nn, "\n")
          raw <- cbind(raw, nedges(nn, Y, edges, D))
          
          # for (rr in unique(Y$data$rid[Y$data$netID == nn])) {
          #   reorder <- c(reorder, which(Y$data$rid == rr))
          # }
          reorder <- c(reorder, which(Y$data$netID == nn))
        }
        
        m <- raw[,order(reorder)]
        
        if (any(uhoh <- (m == 0))) {
            warning("Internal error: disc boundary count equal to zero")
            m[uhoh] <- 1
        }
        edgewt <- 1/m
    }
    
    wt <- edgewt
    if (!is.null(reweight)) 
        wt <- wt * reweight

    K <- compileK(D, r, weights = wt, denom = denom, fname = fname)
    K <- bind.fv(K, data.frame(theo = r), makefvlabel(NULL, NULL, 
                                                      fname, "theo"), "theoretical Poisson %s")
    K <- rebadge.fv(K, ylab, fname)
    unitname(K) <- unitname(X)
    fvnames(K, ".") <- rev(fvnames(K, "."))
    if (showworking) 
        attr(K, "working") <- list(D = D, wt = wt)
    attr(K, "correction") <- correction
    return(K)
}
