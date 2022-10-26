simulate2.lppm <- function (object, nsim = 1, lambda=NULL, ..., correction.factor=NULL, nzw.rids=NULL, new.coef = NULL, progress = (nsim > 1), drop = FALSE, lmax=NULL, loglambda=FALSE) {
    starttime <- proc.time()
    #if (!is.poisson(object$fit)) 
    #    stop("Simulation of non-Poisson models is not yet implemented")

    #if (is.null(lmax)) {
    #    lambda <- predict(object, ..., new.coef = new.coef)
    #    lmax <- if (is.im(lambda)) 
    #        max(lambda)
    #    else unlist(lapply(lambda, max))
    #}    
    L <- as.linnet(object$X)
    result <- vector(mode = "list", length = nsim)
    pstate <- list()
    for (i in seq_len(nsim)) {
        if (progress) 
            pstate <- progressreport(i, nsim, state = pstate)
        result[[i]] <- rpoislpp2(lambda, L, lmax = lmax, nzw.rids=nzw.rids, correction.factor=correction.factor, loglambda=loglambda, ...)
        result[[i]]$data <- as.data.frame(result[[i]]$data)
        
        covars <- list()
        for (cc in names(object$Q$covars)) {
            covarfun <- object$Q$covars[[cc]]
            covars[[cc]] <- covarfun(x=result[[i]]$data$x, y=result[[i]]$data$y, seg=result[[i]]$data$seg, tp=result[[i]]$data$tp)
        }
        covars <- as.data.frame(covars)
        
        result[[i]]$data <- cbind(result[[i]]$data, covars)
    }
    
    if (nsim == 1 && drop) {
        result <- result[[1]]
    }
    else {
        result <- as.solist(result)
        if (nsim > 0) 
            names(result) <- paste("Simulation", 1:nsim)
    }
    result <- timed(result, starttime = starttime)
    return(result)
}
