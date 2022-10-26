#' Estimate the regression parameters of a Cox process model using the variational approximation and assuming independence of the random effects.
#'
#' \code{cox.variational.indep} uses a variational approximation to estimate the parameters of a Cox process regression model with spatial random effects.
#' For this function, the variational approximation to the posterior distribution of the spatial random effects is a multivariate normal with
#' diagonal covariance matrix.
#'
#' @param y vector of response data
#' @param X design matrix for fized effects
#' @param S design matrix for the spatial random effects
#' @param wt vector of observation weights
#' @param beta.start starting values for iteration to estimate the fixed effect coefficients
#' @param tau.start initial value of the precision of the random effects
#' @param tol tolerance for judging convergence of the algorithm
#' @param verbose logical indicating whether to write detailed progress reports to standard output
#' @param hess logical indicating whether to compute the hessian after convergence (slow)
#' @param sd logical indicating whether to report the estimated standard errors of the coefficients
#'
#' @return list composed of these elements:
#'
#' \code{beta}: estimated vector of fixed effect regression coefficients
#'
#' \code{M}: estimated mean vector for the posterior of the spatial random effects at the converged value of the variational approximation
#'
#' \code{diagV}: vector of diagonal entries of the estimated covariance matrix for the posterior of the spatial random effects at convergence of the variational approximation
#' 
#' \code{ltau}: estimated precision component for the spatial random effect
#'
#' \code{hessian}: estimated hessian matrix for \code{beta}, \code{M}, \code{log(diagV)} and \code{ltau} at convergence (estimated by call to \code{optim})
#'
#' \code{sd}: estimated standard errors of \code{beta}, \code{M}, and \code{ltau} at convergence (estimated by call to \code{sdreport})
#'
#' \code{object}: converged model object, returned by \code{nlminb}
#'
#' \code{neg.log.lik}: negative of the variational lower bound on the marginal log-likelihood at convergence


# quick QR factorization for matrices where no row has more than one nonzero entry.
QR.special <- function(A) {
    r <- sqrt(colSums(A))
    q <- t(t(A) / r)
    
    list(q=q, r=r)  
}


cox.variational.ICAR <- function(y, X, S, Q=NULL, sindx, wt, beta.start, means.matrix, u.start=NULL, logV.start=rep(0.01, ncol(S)), ltau.start=1, tol=sqrt(.Machine$double.eps), verbose=TRUE, hess, sd, logDetQ=NULL, maxit=10000) {
    # Start by estimating an optimal log(tau), assuming the given beta.start and u=rep(0,p)
    beta <- beta.start
    logV <- logV.start
    ltau <- ltau.start

    r <- ncol(S)
    p <- ncol(X)
    n <- nrow(X)
    
    if (is.null(logDetQ)) {
      eQ <- eigen(Q)
      logDetQ <- sum(log(eQ$values[eQ$values > sqrt(.Machine$double.eps)]))
    }
    
    if (is.null(u.start)) {
      u <- (1:ncol(S)) / ncol(S) / 100
    } else u <- u.start
    
    # Estimate variance of the variational approximation
    system.time(obj <- MakeADFun(data=list(
        y=y,
        X=X,
        S=S,
        Q=Q,
        logdetQ=logDetQ,
        wt=wt,
        s_indx=sindx,
        means_mat=means.matrix
    ),
    parameters=list(
        ltau = ltau,
        M=u,
        logV = logV,
        beta=beta),
    DLL="tmb_ICAR"))
    
    if (length(obj$par) > 1000) {
      system.time(res <- optim(obj$par, obj$fn, obj$gr, method='CG', control=list(maxit=maxit, reltol=tol)))
    } else {
      system.time(res <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=maxit, eval.max=maxit, rel.tol=tol)))
    }
    # run the variational approximation
    #res <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(iter.max=10000, eval.max=20000))
    
    # extract parameter estimates from the result
    beta <- tail(res$par, p)
    logV <- res$par[1 + r + 1:r]
    u <- res$par[1 + 1:r]
    ltau <- res$par[1]

    # bundle up the return object
    out <- list(beta=beta, M=u, logV=logV, ltau=ltau, neg.loglik=res$value, res=res, obj=obj)
    if (hess) try(out$hessian <- Matrix(optimHess(res$par, fn=obj$fn, gr=obj$gr, y=y, X=X, S=S, wt=wt, control=list(reltol=tol))), silent=TRUE)
    if (sd) try(out$sd <- sdreport(obj), silent=TRUE)
    
    # return the results
    out
}

