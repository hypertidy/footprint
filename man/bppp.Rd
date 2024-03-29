% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/footprint.R
\name{bppp}
\alias{bppp}
\alias{bpsp}
\title{Build a point or line segment pattern object within a raster window.}
\usage{
bppp(x1, pare, diam = 0)

bpsp(x1, x2, pare, diam = 0)
}
\arguments{
\item{x1}{coordinate, x,y two values}

\item{pare}{sparse parent raster, from buildparent}

\item{diam}{diameter used to buffer the grids}

\item{x2}{coordinates, x,y two values}
}
\value{
spatstat point or line segment pattern with built-in raster window
}
\description{
This point or line segment pattern is designed to be a child window from a large raster parent, sharing exact alignment and resolution.
}
\details{
Function \code{bppp} is used for point patterns, \code{bpsp} for line segment patterns
}
