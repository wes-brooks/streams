envelopeEngine <- function (X, fun, simul, nsim = 99, nrank = 1, ..., funargs = list(), 
          funYargs = funargs, verbose = TRUE, clipdata = TRUE, transform = NULL, 
          global = FALSE, ginterval = NULL, use.theory = NULL, alternative = c("two.sided", 
                                                                               "less", "greater"), scale = NULL, clamp = FALSE, savefuns = FALSE, 
          savepatterns = FALSE, saveresultof = NULL, weights = NULL, 
          nsim2 = nsim, VARIANCE = FALSE, nSD = 2, Yname = NULL, maxnerr = nsim, 
          internal = NULL, cl = NULL, envir.user = envir.user, expected.arg = "r", 
          do.pwrong = FALSE, foreignclass = NULL, collectrubbish = FALSE) {
    envir.here <- sys.frame(sys.nframe())
    alternative <- match.arg(alternative)
    foreignclass <- as.character(foreignclass)
    if (length(foreignclass) != 0 && clipdata) {
        warning(paste("Ignoring clipdata=TRUE:", "I don't know how to clip objects of class", 
                      sQuote(paste(foreignclass, collapse = ","))))
        clipdata <- FALSE
    }
    Xclass <- if (is.ppp(X)) 
        "ppp"
    else if (is.pp3(X)) 
        "pp3"
    else if (is.ppx(X)) 
        "ppx"
    else if (inherits(X, foreignclass)) 
        foreignclass
    else stop("Unrecognised class of point pattern")
    Xobjectname <- paste("point pattern of class", sQuote(Xclass))
    if (use.weights <- !is.null(weights)) {
        if (is.numeric(weights)) {
            compute.weights <- FALSE
            weightfun <- NULL
        }
        else if (is.function(weights)) {
            compute.weights <- TRUE
            weightfun <- weights
            weights <- NULL
        }
        else stop("weights should be either a function or a numeric vector")
    }
    else compute.weights <- FALSE
    patterns.only <- identical(internal$eject, "patterns")
    if (savevalues <- !is.null(saveresultof)) {
        stopifnot(is.function(saveresultof))
        SavedValues <- list()
    }
    if (inherits(simul, "simulrecipe")) {
        simtype <- simul$type
        simexpr <- simul$expr
        envir <- simul$envir
        csr <- simul$csr
        pois <- simul$pois
        constraints <- simul$constraints
    }
    else {
        simulate <- simul
        csr <- FALSE
        if (!is.null(icsr <- internal$csr)) 
            csr <- icsr
        pois <- csr
        constraints <- ""
        if (inherits(simulate, "envelope")) {
            simpat <- attr(simulate, "simpatterns")
            if (!is.null(simpat)) 
                simulate <- simpat
            else stop(paste("The argument", sQuote("simulate"), 
                            "is an envelope object but does not contain", 
                            "any saved point patterns."))
        }
        if (is.expression(simulate)) {
            simtype <- "expr"
            simexpr <- simulate
            envir <- envir.user
        }
        else if (is.list(simulate) && all(sapply(simulate, inherits, 
                                                 what = Xclass))) {
            simtype <- "list"
            SimDataList <- simulate
            simexpr <- expression(SimDataList[[i]])
            dont.complain.about(SimDataList)
            envir <- envir.here
            i <- 1
            if (!is.null(mess <- attr(simulate, "internal"))) {
                csr <- !is.null(mc <- mess$csr) && mc
            }
        }
        else stop(paste(sQuote("simulate"), "should be an expression, or a list of point patterns"))
    }
    if (clipdata) {
        Xsim <- eval(simexpr, envir = envir)
        if (!inherits(Xsim, Xclass)) 
            switch(simtype, csr = stop(paste("Internal error:", 
                                             Xobjectname, "not generated")), rmh = stop(paste("Internal error: rmh did not return an", 
                                                                                              Xobjectname)), kppm = stop(paste("Internal error: simulate.kppm did not return an", 
                                                                                                                               Xobjectname)), expr = stop(paste("Evaluating the expression", 
                                                                                                                                                                sQuote("simulate"), "did not yield an", Xobjectname)), 
                   list = stop(paste("Internal error: list entry was not an", 
                                     Xobjectname)), stop(paste("Internal error:", 
                                                               Xobjectname, "not generated")))
        clipwin <- Xsim$window
        if (!is.subset.owin(clipwin, X$window)) 
            warning("Window containing simulated patterns is not a subset of data window")
    }
    if (is.null(fun)) 
        stop("Internal error: fun is NULL")
    fname <- if (is.name(substitute(fun))) 
        short.deparse(substitute(fun))
    else if (is.character(fun)) 
        fun
    else "fun"
    fname <- sQuote(fname)
    if (is.character(fun)) {
        gotfun <- try(get(fun, mode = "function"))
        if (inherits(gotfun, "try-error")) 
            stop(paste("Could not find a function named", sQuote(fun)))
        fun <- gotfun
    }
    else if (!is.function(fun)) 
        stop(paste("unrecognised format for function", fname))
    fargs <- names(formals(fun))
    if (!any(c(expected.arg, "...") %in% fargs)) 
        stop(paste(fname, "should have", ngettext(length(expected.arg), 
                                                  "an argument", "arguments"), "named", commasep(sQuote(expected.arg)), 
                   "or a", sQuote("..."), "argument"))
    usecorrection <- any(c("correction", "...") %in% fargs)
    if ((nrank%%1) != 0) 
        stop("nrank must be an integer")
    if ((nsim%%1) != 0) 
        stop("nsim must be an integer")
    stopifnot(nrank > 0 && nrank < nsim/2)
    rgiven <- any(expected.arg %in% names(list(...)))
    if (tran <- !is.null(transform)) {
        stopifnot(is.expression(transform))
        transform.funX <- inject.expr("with(funX,.)", transform)
        transform.funXsim <- inject.expr("with(funXsim,.)", transform)
    }
    if (!is.null(ginterval)) 
        stopifnot(is.numeric(ginterval) && length(ginterval) == 
                      2)
    Xarg <- if (!clipdata) 
        X
    else X[clipwin]
    corrx <- if (usecorrection) 
        list(correction = "best")
    else NULL
    funX <- do.call(fun, resolve.defaults(list(Xarg), list(...), 
                                          funYargs, corrx))
    if (!inherits(funX, "fv")) 
        stop(paste("The function", fname, "must return an object of class", 
                   sQuote("fv")))
    if (!is.null(dang <- attr(funX, "dangerous")) && any(uhoh <- dang %in% 
                                                         names(list(...)))) {
        nuh <- sum(uhoh)
        warning(paste("Envelope may be invalid;", ngettext(nuh, 
                                                           "argument", "arguments"), commasep(sQuote(dang[uhoh])), 
                      ngettext(nuh, "appears", "appear"), "to have been fixed."), 
                call. = FALSE)
    }
    argname <- fvnames(funX, ".x")
    valname <- fvnames(funX, ".y")
    has.theo <- "theo" %in% fvnames(funX, "*")
    csr.theo <- csr && has.theo
    use.theory <- if (is.null(use.theory)) 
        csr.theo
    else (use.theory && has.theo)
    if (tran) {
        if (use.theory) 
            funX <- funX[, c(argname, valname, "theo")]
        else funX <- funX[, c(argname, valname)]
        funX <- eval(transform.funX)
    }
    rvals <- funX[[argname]]
    alim <- attr(funX, "alim")
    if (global && is.null(ginterval)) 
        ginterval <- if (rgiven || is.null(alim)) 
            range(rvals)
    else alim
    dual <- (global && !use.theory && !VARIANCE)
    Nsim <- if (!dual) 
        nsim
    else (nsim + nsim2)
    if (simtype == "list" && Nsim > length(SimDataList)) 
        stop(paste("Number of simulations", paren(if (!dual) 
            paste(nsim)
            else paste(nsim, "+", nsim2, "=", Nsim)), "exceeds number of point pattern datasets supplied"))
    if (patterns.only) {
        if (verbose) {
            action <- if (simtype == "list") 
                "Extracting"
            else "Generating"
            descrip <- switch(simtype, csr = "simulations of CSR", 
                              rmh = paste("simulated realisations of fitted", 
                                          if (pois) "Poisson" else "Gibbs", "model"), 
                              kppm = "simulated realisations of fitted cluster model", 
                              expr = "simulations by evaluating expression", 
                              list = "point patterns from list", "simulated realisations")
            if (!is.null(constraints) && nzchar(constraints)) 
                descrip <- paste(descrip, constraints)
            explan <- if (dual) 
                paren(paste(nsim2, "to estimate the mean and", 
                            nsim, "to calculate envelopes"))
            else ""
            splat(action, Nsim, descrip, explan, "...")
        }
        XsimList <- list()
        sstate <- list()
        for (i in 1:Nsim) {
            if (verbose) 
                sstate <- progressreport(i, Nsim, state = sstate)
            Xsim <- eval(simexpr, envir = envir)
            if (!inherits(Xsim, Xclass)) 
                switch(simtype, csr = {
                    stop(paste("Internal error:", Xobjectname, 
                               "not generated"))
                }, rmh = {
                    stop(paste("Internal error: rmh did not return an", 
                               Xobjectname))
                }, kppm = {
                    stop(paste("Internal error: simulate.kppm did not return an", 
                               Xobjectname))
                }, expr = {
                    stop(paste("Evaluating the expression", sQuote("simulate"), 
                               "did not yield an", Xobjectname))
                }, list = {
                    stop(paste("Internal error: list entry was not an", 
                               Xobjectname))
                }, stop(paste("Internal error:", Xobjectname, 
                              "not generated")))
            XsimList[[i]] <- Xsim
        }
        if (verbose) {
            cat(paste("Done.\n"))
            flush.console()
        }
        attr(XsimList, "internal") <- list(csr = csr)
        return(XsimList)
    }
    envelopeInfo <- list(call = cl, Yname = Yname, valname = valname, 
                         csr = csr, csr.theo = csr.theo, use.theory = use.theory, 
                         pois = pois, simtype = simtype, constraints = constraints, 
                         nrank = nrank, nsim = nsim, Nsim = Nsim, global = global, 
                         ginterval = ginterval, dual = dual, nsim2 = nsim2, VARIANCE = VARIANCE, 
                         nSD = nSD, alternative = alternative, scale = scale, 
                         clamp = clamp, use.weights = use.weights, do.pwrong = do.pwrong)
    if (verbose) {
        action <- if (simtype == "list") 
            "Extracting"
        else "Generating"
        descrip <- switch(simtype, csr = "simulations of CSR", 
                          rmh = paste("simulated realisations of fitted", if (pois) "Poisson" else "Gibbs", 
                                      "model"), kppm = "simulated realisations of fitted cluster model", 
                          expr = "simulations by evaluating expression", list = "point patterns from list", 
                          "simulated patterns")
        if (!is.null(constraints) && nzchar(constraints)) 
            descrip <- paste(descrip, constraints)
        explan <- if (dual) 
            paren(paste(nsim2, "to estimate the mean and", nsim, 
                        "to calculate envelopes"))
        else ""
        splat(action, Nsim, descrip, explan, "...")
    }
    catchpatterns <- savepatterns && simtype != "list"
    Caughtpatterns <- list()
    nrvals <- length(rvals)
    simvals <- matrix(, nrow = nrvals, ncol = Nsim)
    if (compute.weights) 
        weights <- numeric(Nsim)
    if (identical(expected.arg, "r")) {
        inferred.r.args <- list(r = rvals)
    }
    else if (identical(expected.arg, c("rmax", "nrval"))) {
        inferred.r.args <- list(rmax = max(rvals), nrval = length(rvals))
    }
    else stop(paste("Don't know how to infer values of", commasep(expected.arg)))
    funargs <- resolve.defaults(funargs, inferred.r.args, list(...), 
                                if (usecorrection) 
                                    list(correction = "best")
                                else NULL)
    nerr <- 0
    if (verbose) 
        pstate <- list()
    for (i in 1:Nsim) {
        ok <- FALSE
        while (!ok) {
            Xsim <- eval(simexpr, envir = envir)
            if (!inherits(Xsim, Xclass)) 
                switch(simtype, csr = stop(paste("Internal error:", 
                                                 Xobjectname, "not generated")), rmh = stop(paste("Internal error: rmh did not return an", 
                                                                                                  Xobjectname)), kppm = stop(paste("Internal error:", 
                                                                                                                                   "simulate.kppm did not return an", Xobjectname)), 
                       expr = stop(paste("Evaluating the expression", 
                                         sQuote("simulate"), "did not yield an", Xobjectname)), 
                       list = stop(paste("Internal error: list entry was not an", 
                                         Xobjectname)), stop(paste("Internal error:", 
                                                                   Xobjectname, "not generated")))
            if (catchpatterns) 
                Caughtpatterns[[i]] <- Xsim
            if (savevalues) 
                SavedValues[[i]] <- saveresultof(Xsim)
            if (compute.weights) {
                wti <- weightfun(Xsim)
                if (!is.numeric(wti)) 
                    stop("weightfun did not return a numeric value")
                if (length(wti) != 1) 
                    stop("weightfun should return a single numeric value")
                weights[i] <- wti
            }
            funXsim <- try(do.call(fun, append(list(Xsim), funargs)))
            ok <- !inherits(funXsim, "try-error")
            if (!ok) {
                nerr <- nerr + 1
                if (nerr > maxnerr) 
                    stop("Exceeded maximum number of errors")
                cat("[retrying]\n")
            }
        }
        if (i == 1) {
            if (!inherits(funXsim, "fv")) 
                stop(paste("When applied to a simulated pattern, the function", 
                           fname, "did not return an object of class", 
                           sQuote("fv")))
            argname.sim <- fvnames(funXsim, ".x")
            valname.sim <- fvnames(funXsim, ".y")
            if (argname.sim != argname) 
                stop(paste("The objects returned by", fname, 
                           "when applied to a simulated pattern", "and to the data pattern", 
                           "are incompatible. They have different argument names", 
                           sQuote(argname.sim), "and", sQuote(argname), 
                           "respectively"))
            if (valname.sim != valname) 
                stop(paste("When", fname, "is applied to a simulated pattern", 
                           "it provides an estimate named", sQuote(valname.sim), 
                           "whereas the estimate for the data pattern is named", 
                           sQuote(valname), ". Try using the argument", 
                           sQuote("correction"), "to make them compatible"))
            rfunX <- with(funX, ".x")
            rfunXsim <- with(funXsim, ".x")
            if (!identical(rfunX, rfunXsim)) 
                stop(paste("When", fname, "is applied to a simulated pattern,", 
                           "the values of the argument", sQuote(argname.sim), 
                           "are different from those used for the data."))
        }
        if (tran) {
            if (use.theory) 
                funXsim <- funXsim[, c(argname, valname, "theo")]
            else funXsim <- funXsim[, c(argname, valname)]
            funXsim <- eval(transform.funXsim)
        }
        simvals.i <- funXsim[[valname]]
        if (length(simvals.i) != nrvals) 
            stop("Vectors of function values have incompatible lengths")
        simvals[, i] <- funXsim[[valname]]
        if (verbose) 
            pstate <- progressreport(i, Nsim, state = pstate)
        if (collectrubbish) {
            rm(Xsim)
            rm(funXsim)
            gc()
        }
    }
    if (verbose) {
        cat("\nDone.\n")
        flush.console()
    }
    if (savefuns) {
        alldata <- cbind(rvals, simvals)
        simnames <- paste("sim", 1:Nsim, sep = "")
        colnames(alldata) <- c("r", simnames)
        alldata <- as.data.frame(alldata)
        SimFuns <- fv(alldata, argu = "r", ylab = attr(funX, 
                                                       "ylab"), valu = "sim1", fmla = deparse(. ~ r), alim = attr(funX, 
                                                                                                                  "alim"), labl = names(alldata), desc = c("distance argument r", 
                                                                                                                                                           paste("Simulation ", 1:Nsim, sep = "")), fname = attr(funX, 
                                                                                                                                                                                                                 "fname"), yexp = attr(funX, "yexp"), unitname = unitname(funX))
        fvnames(SimFuns, ".") <- simnames
    }
    if (savepatterns) 
        SimPats <- if (simtype == "list") 
            SimDataList
    else Caughtpatterns
    etype <- if (global) 
        "global"
    else if (VARIANCE) 
        "variance"
    else "pointwise"
    if (dual) {
        jsim <- 1:nsim
        jsim.mean <- nsim + 1:nsim2
    }
    else {
        jsim <- jsim.mean <- NULL
    }
    result <- envelope.matrix(simvals, funX = funX, jsim = jsim, 
                              jsim.mean = jsim.mean, type = etype, alternative = alternative, 
                              scale = scale, clamp = clamp, csr = csr, use.theory = use.theory, 
                              nrank = nrank, ginterval = ginterval, nSD = nSD, Yname = Yname, 
                              do.pwrong = do.pwrong, weights = weights)
    attr(result, "einfo") <- envelopeInfo
    if (savefuns) 
        attr(result, "simfuns") <- SimFuns
    if (savepatterns) {
        attr(result, "simpatterns") <- SimPats
        attr(result, "datapattern") <- X
    }
    if (use.weights) 
        attr(result, "weights") <- weights
    if (savevalues) {
        attr(result, "simvalues") <- SavedValues
        attr(result, "datavalue") <- saveresultof(X)
    }
    return(result)
}