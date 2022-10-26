# find the hessian of the fixed-effect coefficients by wrapping the objective and 
# gradient functions so that only the fixed effets can vary.
# training: allow the regression coefficients to vary. parameters masked TRUE can vary.

subset.wrapper <- function(fun, mask, masked.pars, vector.result=FALSE) {
  paramvec <- masked.pars
  
  out <- function(x=NULL) {
    if (!is.null(x))
      paramvec[mask] <- x
    
    if (vector.result)
      fun(paramvec)[,mask, drop=FALSE]
    else
      fun(paramvec)
  }
}