
##extent(range(c(stxy[,1], enxy[,1])), range(c(stxy[,2], enxy[,2]))) + diam * 2
buf <- function(d, b, n = 0) {
  ((d %/% b) + n) * b
}

#' Build large parent raster
#'
#' Create a large sparse raster from an extent and resolution. This grid is never meant to be populated with data
#' but used as a specification for smaller child windows that are aligned exactly within this larger, high-resolution parent. 
#' 
#' Currently the resolution is required to be square. 
#' @param x extent
#' @param rr resolution
#'
#' @return RasterLayer
#' @export
#'
buildparent <- function(x, rr) {
  ## x the extent of all data
  ext <- extent(buf(xmin(x), rr),
                buf(xmax(x), rr, 1),
                buf(ymin(x), rr),
                buf(ymax(x), rr, 1))

  ##xl <- ((range(x[,1] + diam * c(-1, 1)) %/% rr) + c(-1, 1) ) * rr
  ##  yl <- ((range(x[,2] + diam * c(-1, 1)) %/% rr) + c(-1, 1) ) * rr

  x <- raster(ext)
  res(x) <- rr
  x
}
extent.psp <- function(x) {
  extent(as.owin(x)$xrange, as.owin(x)$yrange)
}


#' Build a point or line segment pattern object within a raster window. 
#' 
#' This point or line segment pattern is designed to be a child window from a large raster parent, sharing exact alignment and resolution. 
#'
#' Function \code{bppp} is used for point patterns, \code{bpsp} for line segment patterns
#' @param x1 coordinate, x,y two values
#' @param x2 coordinates, x,y two values
#' @param pare sparse parent raster, from buildparent
#' @param diam diameter used to buffer the grids
#'
#' @return spatstat point or line segment pattern with built-in raster window
#' @export
#'
#' @importFrom raster extent res crop 
#' @importFrom spatstat owin psp ppp
bppp <- function(x1, pare, diam = 0) {
  x <- x1
  e <- extent(x[1] + c(-1, 1) *  diam, x[2] + c(-1, 1) *  diam)
  rr <- res(pare)[1]
  ## xl <- ((range(x[,1] + diam * c(-1, 1)) %/% rr) + c(-1, 1) ) * rr
  ##  yl <- ((range(x[,2] + diam * c(-1, 1)) %/% rr) + c(-1, 1) ) * rr
  ext <- extent(buf(xmin(e), rr),
                buf(xmax(e), rr, 1),
                buf(ymin(e), rr),
                buf(ymax(e), rr, 1))

  r <- crop(pare, ext + diam * 2, snap = "out")
  r[] <- 1L
  msk <- as.matrix(r)
  msk <- matrix(TRUE, nrow(msk), ncol(msk))
  w <- owin(c(xmin(r), xmax(r)), c(ymin(r), ymax(r)), mask = msk)
  ppp(x1[1], x1[2], window = w)
}


#' @export
#' @rdname bppp
bpsp <- function(x1, x2, pare, diam = 0) {
  x <- rbind(x1, x2)
  e <- extent(t(x))  ## ouch, two rows!
  rr <- res(pare)[1]
  ## xl <- ((range(x[,1] + diam * c(-1, 1)) %/% rr) + c(-1, 1) ) * rr
  ##  yl <- ((range(x[,2] + diam * c(-1, 1)) %/% rr) + c(-1, 1) ) * rr
  ext <- extent(buf(xmin(e), rr),
                buf(xmax(e), rr, 1),
                buf(ymin(e), rr),
                buf(ymax(e), rr, 1))

  r <- crop(pare, ext + diam * 2, snap = "out")
  r[] <- 1L
  #msk <- as.matrix(r)
  msk <- matrix(TRUE, nrow(r), ncol(r))
  w <- owin(c(xmin(r), xmax(r)), c(ymin(r), ymax(r)), mask = msk)
  psp(x1[1], x1[2], x2[1], x2[2], window = w)
}
