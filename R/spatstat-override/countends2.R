countends2 <- function(L, x=locator(1), r, toler=NULL) {
  # L is the linear network (object of class "linnet")
  # x is the centre point of the disc
  # r is the radius of the disc
  #
  stopifnot(inherits(L, "linnet"))
  lines <- L$lines
  vertices <- L$vertices
  lengths <- lengths.psp(lines)
  dpath <- L$dpath
  win <- L$window
  nv <- vertices$n
  ns <- lines$n
  # get x
  # if(missing(x))
  #   x <- clickppp(1, win, add=TRUE)
  # else
  #   x <- as.ppp(x, win)
  #
  np <- npoints(x)
  if(length(r) != np)
    stop("Length of vector r does not match number of points in x")
  # project x to nearest segment
  pro <- project2segment(x, lines)
  # which segment?
  startsegment <- pro$mapXY
  # parametric position of x along this segment
  startfraction <- pro$tp
  
  # convert indices to C 
  seg0 <- startsegment - 1L
  from0 <- L$from - 1L
  to0   <- L$to - 1L
  if(is.null(toler)) {
    toler <- 0.001 * min(lengths[lengths > 0])
    toler <- max(.Machine$double.xmin, toler, finite=TRUE)
  } else {
    check.1.real(toler)
    stopifnot(toler > 0)
  }
  zz <- .C("Ccountends",
           np = as.integer(np),
           f = as.double(startfraction),
           seg = as.integer(seg0),
           r = as.double(r), 
           nv = as.integer(nv), 
           xv = as.double(vertices$x),
           yv = as.double(vertices$y),  
           ns = as.integer(ns),
           from = as.integer(from0),
           to = as.integer(to0), 
           dpath = as.double(dpath),
           lengths = as.double(lengths),
           toler=as.double(toler),
           nendpoints = as.integer(integer(np)))
  zz$nendpoints
}


# countends.ssnlpp <- function(L, x=locator(1), r, toler=NULL) {
#   # L is the linear network (object of class "linnet")
#   # x is the centre point of the disc
#   # r is the radius of the disc
#   #
#   stopifnot(inherits(L, "linnet"))
#   lines <- L$lines
#   vertices <- L$vertices
#   lengths <- lengths.psp(lines)
#   dpath <- L$dpath
#   win <- L$window
#   nv <- vertices$n
#   ns <- lines$n
#   # get x
#   # if(missing(x))
#   #   x <- clickppp(1, win, add=TRUE)
#   # else
#   #   x <- as.ppp(x, win)
#   #
#   np <- npoints(x)
#   if(length(r) != np)
#     stop("Length of vector r does not match number of points in x")
#   # project x to nearest segment
#   pro <- project2segment(x, lines)
#   # which segment?
#   startsegment <- pro$mapXY
#   # parametric position of x along this segment
#   startfraction <- pro$tp
#   
#   # convert indices to C 
#   seg0 <- startsegment - 1L
#   from0 <- L$from - 1L
#   to0   <- L$to - 1L
#   if(is.null(toler)) {
#     toler <- 0.001 * min(lengths[lengths > 0])
#     toler <- max(.Machine$double.xmin, toler, finite=TRUE)
#   } else {
#     check.1.real(toler)
#     stopifnot(toler > 0)
#   }
#   zz <- .C("Ccountends",
#            np = as.integer(np),
#            f = as.double(startfraction),
#            seg = as.integer(seg0),
#            r = as.double(r), 
#            nv = as.integer(nv), 
#            xv = as.double(vertices$x),
#            yv = as.double(vertices$y),  
#            ns = as.integer(ns),
#            from = as.integer(from0),
#            to = as.integer(to0), 
#            dpath = as.double(dpath),
#            lengths = as.double(lengths),
#            toler=as.double(toler),
#            nendpoints = as.integer(integer(np)))
#   zz$nendpoints
# }
