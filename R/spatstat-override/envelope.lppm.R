envelope.lppm <- function (Y, fun = linearK, nsim = 99, nrank = 1, ..., funargs = list(), 
          funYargs = funargs, simulate = NULL, verbose = TRUE, transform = NULL, 
          global = FALSE, ginterval = NULL, use.theory = NULL, alternative = c("two.sided", 
                                                                               "less", "greater"), scale = NULL, clamp = FALSE, savefuns = FALSE, 
          savepatterns = FALSE, nsim2 = nsim, VARIANCE = FALSE, nSD = 2, 
          Yname = NULL, do.pwrong = FALSE, envir.simul = NULL) {
    cl <- short.deparse(sys.call())
    if (is.null(Yname)) 
        Yname <- short.deparse(substitute(Y))
    if (is.null(fun)) 
        fun <- linearK
    if ("clipdata" %in% names(list(...))) 
        stop(paste("The argument", sQuote("clipdata"), "is not available for envelope.pp3"))
    envir.user <- if (!is.null(envir.simul)) 
        envir.simul
    else parent.frame()
    envir.here <- sys.frame(sys.nframe())
    if (is.null(simulate)) {
        if (!is.poisson.ppm(Y$fit)) 
            stop("Simulation of non-Poisson models is not yet implemented")
        X <- Y$X
        MODEL <- Y
        NETWORK <- X$domain 
        lambdaFit <- predict(MODEL)
        LMAX <- if (is.im(lambdaFit)) 
            max(lambdaFit)
        else unlist(lapply(lambdaFit, max))
        simexpr <- expression(rpoislpp(lambdaFit, NETWORK, lmax = LMAX))
        dont.complain.about(NETWORK, LMAX)
        simrecipe <- simulrecipe(type = "lppm", expr = simexpr, 
                                 envir = envir.here, csr = FALSE)
    }
    else {
        simrecipe <- simulate
        X <- Y
    }
    envelopeEngine(X = X, fun = fun, simul = simrecipe, nsim = nsim, 
                   nrank = nrank, ..., funargs = funargs, funYargs = funYargs, 
                   verbose = verbose, clipdata = FALSE, transform = transform, 
                   global = global, ginterval = ginterval, use.theory = use.theory, 
                   alternative = alternative, scale = scale, clamp = clamp, 
                   savefuns = savefuns, savepatterns = savepatterns, nsim2 = nsim2, 
                   VARIANCE = VARIANCE, nSD = nSD, Yname = Yname, cl = cl, 
                   envir.user = envir.user, do.pwrong = do.pwrong)
}