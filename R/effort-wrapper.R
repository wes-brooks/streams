effort.fun.kernel <- function(x, y, seg, tp, data, lookup, kernel.function, coef=1, raw=FALSE, logresult=FALSE) {
  
  # straight out of the kernel, this is raw and unlogged.
  out <- kernel.function(x, y, seg, tp)
  
  if (!raw && logresult) {
    out <- coef * log(out)
  } else if (!raw && !logresult) {
    out <- exp(coef * log(out))
  } else if (raw && logresult) {
    out <- log(out)
  } #else if (raw && !logresult) {
    #out <- out
  #}
  
  out
}


effort.wrapper.kernel <- function(data, lookup, coef, kernel.function) {
  wrapped <- function(x, y, seg, tp, raw=FALSE, logresult=FALSE) {
    effort.fun.kernel(x, y, seg, tp, data, lookup, kernel.function, coef, raw, logresult)
  }
  
  wrapped
}


effort.fun <- function(x, y, seg, tp, data, lookup, coef=1, raw=FALSE, logresult=FALSE) {
  indx <- match(seg, lookup$seg)
  indx <- match(lookup$rid[indx], data$rid)
  
  if (raw && logresult) {
    out <- ifelse(!is.na(indx), log(data$effort[indx]), 0)
  } else if (!raw && !logresult) {
    out <- ifelse(!is.na(indx), exp(coef * log(data$effort[indx])), 0)
  } else if (!raw && logresult) {
    out <- ifelse(!is.na(indx), coef * log(data$effort[indx]), 0)
  } else { # if (raw && !logresult) {
    out <- ifelse(!is.na(indx), data$effort[indx], 0)
  }
  
  out
}



effort.wrapper <- function(func, data, lookup, coef) {
  wrapped <- function(x, y, seg, tp, raw=FALSE, logresult=FALSE) {
    func(x, y, seg, tp, data, lookup, coef, raw, logresult)
  }
  
  return(wrapped)
}