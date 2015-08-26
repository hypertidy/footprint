
##extent(range(c(stxy[,1], enxy[,1])), range(c(stxy[,2], enxy[,2]))) + diam * 2
buf <- function(d, b, n = 0) {
  ((d %/% b) + n) * b
}

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
  msk <- as.matrix(r)
  msk <- matrix(TRUE, nrow(msk), ncol(msk))
  w <- owin(c(xmin(r), xmax(r)), c(ymin(r), ymax(r)), mask = msk)
  psp(x1[1], x1[2], x2[1], x2[2], window = w)
}
