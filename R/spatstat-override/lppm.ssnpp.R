lppm.ssnpp <- function (X, Q, epsilon=1e-6, ...)  {
    
    
   
  
    varnames <- Q$varnames
    
    
    Xname <- short.deparse(substitute(X))
    callstring <- paste(short.deparse(sys.call()), collapse = "")
    cl <- match.call()
    nama <- names(list(...))
    resv <- c("method", "forcefit")
    if (any(clash <- resv %in% nama)) 
        warning(paste(ngettext(sum(clash), "Argument", "Arguments"), 
                      commasep(sQuote(resv[clash])), "must not be used"))
    # stopifnot(inherits(X, "ssnpp"))
    
    data.points <- list()
    for (i in 1:length(Q$covars)) 
        data.points[[varnames[i]]] <- Q$covars[[i]](x=X$data$x, y=X$data$y, seg=X$data$seg, tp=X$data$tp)
    
    data.points <- as.data.frame(data.points)
    data <- rbind(data.points, Q$data.quad[,varnames])
    
    # trend <- as.formula(Q$trend)
    # mf <- model.frame(trend, data=data)
    # mm <- model.matrix(trend, data=mf)
    
    points <- lpp(X$data[,c('seg','tp')], L=X$domain)
    
        quadrature <- quad(data = points, dummy = Q$quad.points, w = c(rep(epsilon, nrow(X$data)), Q$wt), param = list(dummy = list(method = "Spaced about every 1km along the stream network"), weight = list(method = "Counting weights based on segment length")))
    #Q <- linequad(X, eps = eps, nd = nd)
    fit <- ppm(quadrature, ..., method = "mpl", forcefit = TRUE, trend=Q$trend, data=data)
    if (!is.poisson.ppm(fit))  
        warning("Non-Poisson models currently use Euclidean distance")
    out <- list(X = X, fit = fit, Xname = Xname, call = cl, callstring = callstring, Q=Q)
    class(out) <- "lppm"
    return(out)
    
    
    
    

    # 
    # ux = unique(quad.x)
    # uy = unique(quad.y)
    # ux = sort(ux)
    # uy = sort(uy)
    # nx = length(ux)
    # ny = length(uy)
    # 
    # quad.mask = matrix(NA, ny, nx, dimnames = list(uy, ux))
    # col.ref   = match(quad.x, ux)
    # row.ref   = match(quad.y, uy)
    # 
    # all.vec = rep(NA, max(row.ref)*max(col.ref))
    # vec.ref = (col.ref - 1)*max(row.ref) + row.ref
    # all.vec[vec.ref] = 1
    # num.vec = all.vec
    # num.vec[is.na(all.vec)] = 0
    
    
    call.1 = quote(lppm(Q, trend = trend, data=model.frame(dwpr), interaction = Poisson(), correction = "none"))
    if (is.ai)
    {
        call.1 = bquote(ppm(Q = Q, trend = trend, covariates = cov.list, interaction = AreaInter(.(fit$r)), correction = "none"))
    }
    
    class(fit)       = "lppm"
    fit$fitter       = "glm"
    fit$coef         = fit$coefficients
    fit$method       = "mpl"
    fit$projected    = FALSE
    fit$trend        = trend
    fit$interaction  = NULL
    if (is.ai)
    {
        fit$interaction = AreaInter(fit$r)
    }
    fit.int          = list(name = "Poisson process", creator = "Poisson", family = NULL, pot = NULL, par = NULL, parnames = NULL)
    class(fit.int)   = "interact"
    fit$fitin        = list(interaction = Poisson(), coefs = fit$beta, Vnames = character(0), IsOffset = logical(0))
    if (is.ai)
    {
        fit$fitin        = list(interaction = AreaInter(fit$r), coefs = fit$beta, Vnames = "Interaction", IsOffset = FALSE)
    }
    class(fit$fitin) = c("fii", "list")
    fit$Q            = Q
    fit$maxlogpl     = fit$loglik
    fit$covariates   = cov.list
    fit$covfunargs   = list()
    fit$correction   = "none"
    fit$rbord        = 0
    
    glmdata          = data.frame(fit$prior.weights, fit$data$resp, fit$data[,-1], TRUE)
    if (is.ai == FALSE)
    {
        names(glmdata)   = c(".mpl.W", ".mpl.Y", names(cov.list), ".mpl.SUBSET")
    }
    if (is.ai)
    {
        names(glmdata)   = c(".mpl.W", ".mpl.Y", names(cov.list), "Interaction", ".mpl.SUBSET")
    }
    fit$version     = list(major = 1, minor = 31, release = 0)
    fit$problems    = list()
    fit$call        = call.1
    fit$callstring  = "character"
    fit$callframe   = environment()
    
    terms.int                 = terms(glmfit.form)
    terms.ppm                 = terms(trend)
    fit$terms                 = terms.ppm
    fit$internal$glmfit$terms = terms.int
    
    # glm.1 = glm(fit$pres/fit$prior.weights ~ as.matrix(fit$data[,-1]), weights = fit$prior.weights, family = poisson())
    # 
    # glm.1$coefficients      = fit$coefficients
    # glm.1$fitted.values     = fit$fitted.values
    # glm.1$residuals         = (fit$pres/fit$prior.weights - fit$fitted.values) / fit$fitted.values
    # glm.1$linear.predictors = log(fit$fitted.values)
    
    fam   = poisson()
    vari  = fam$variance(fit$fitted.values)
    deriv = 1/vari
    wt    = fit$prior.weights
    weii  = wt*1/(deriv^2*vari)
    w.1   = weii
    Xw    = t(as.vector(sqrt(weii)) * t(t(fit$data)))
    q.1   = qr(t(Xw))
    
    glm.1$qr      = q.1
    glm.1$weights = w.1
    glm.1$family  = quasi(log)
    glm.1$data    = glmdata
    glm.1$formula = glmfit.form
    
    fit$internal  = list(glmfit = dwpr, glmdata = glmdata)
    if (is.ai)
    {
        fit$internal  = list(glmfit = glm.1, glmdata = glmdata, Vnames = "Interaction", IsOffset = FALSE)
    }
    
    fit
}