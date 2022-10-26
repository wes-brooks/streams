quad.dwpr <- function(obj) {
    
    is.ai  = FALSE #is.numeric(fit$pt.interactions)
    quad.x = obj$loc.q$x #fit$x[fit$pres == 0]
    quad.y = obj$loc.q$y #fit$y[fit$pres == 0]
    
    
    
    
    qwt <- tail(obj$wt, length(quad.x))
    
    
    
    quad.win  = obj$window
    quad.dat = lpp(ppp(quad.x, quad.y, window = quad.win, check = FALSE), L=net)

    
    
    # find the names of the covariates
    num.var = length(obj$covar.names)
    cov.list = vector('list', num.var)
    varnames <- obj$covar.names
    
    # create a linfun object to represent the value of each covariate on the linear network
    for(var in 1:length(varnames)) {
        f <- linfun.factory(obj$data, varnames[[var]])
        cov.list[[var]] <- linfun(f, net)
    }
    
    # trend is the right hand side of the model formula
    names(cov.list) <- varnames
    trend <- paste("~ ", tail(as.character(obj$formula), 1), sep='')
    trend <- formula(trend)

    list(quad.points = quad.dat, covars = cov.list, trend = trend, varnames=obj$covar.names, data.quad=obj$quaddata, wt=qwt)

}