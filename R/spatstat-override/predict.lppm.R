predict.lppm <- function (object, ..., type = "trend", locations = NULL, new.coef = NULL) {
    type <- pickoption("type", type, c(trend = "trend", cif = "cif", 
                                       lambda = "cif"))
    X <- object$X
    fit <- object$fit
    L <- as.linnet(X)
    if (!is.null(locations)) {
        values <- predict(fit, locations = locations, type = type, 
                          new.coef = new.coef)
        return(values)
    }
    Llines <- as.psp(L)
    linemask <- as.mask.psp(Llines, ...)
    lineimage <- as.im(linemask)
    xx <- rasterx.mask(linemask)
    yy <- rastery.mask(linemask)
    mm <- linemask$m
    xx <- as.vector(xx[mm])
    yy <- as.vector(yy[mm])
    pixelcentres <- ppp(xx, yy, window = as.rectangle(linemask), 
                        check = FALSE)
    pixdf <- data.frame(xc = xx, yc = yy)
    p2s <- project2segment(pixelcentres, Llines)
    projloc <- as.data.frame(p2s$Xproj)
    projmap <- as.data.frame(p2s[c("mapXY", "tp")])
    projdata <- cbind(pixdf, projloc, projmap)
    if (!is.multitype(fit)) {
        values <- predict(fit, locations = projloc, type = type, 
                          new.coef = new.coef)
        Z <- lineimage
        Z[pixelcentres] <- values
        df <- cbind(projdata, values)
        out <- linim(L, Z, df = df)
    }
    else {
        lev <- levels(marks(data.ppm(fit)))
        out <- list()
        for (k in seq(length(lev))) {
            markk <- factor(lev[k], levels = lev)
            locnk <- cbind(projloc, data.frame(marks = markk))
            values <- predict(fit, locations = locnk, type = type, 
                              new.coef = new.coef)
            Z <- lineimage
            Z[pixelcentres] <- values
            df <- cbind(projdata, values)
            out[[k]] <- linim(L, Z, df = df)
        }
        names(out) <- as.character(lev)
        class(out) <- as.solist(out)
    }
    return(out)
}