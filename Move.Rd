\name{Move}
\alias{Move}
\alias{Move.MCBS}
\title{Moves Single Photon in Media}
\description{
    This function calculates the move of single photon in media.
}
\usage{
Move(myObject)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
}
\details{
    The \code{Move} function, alias \code{Move.MCBS}, calculates updated position coordinates \code{c(x,y,z)} according to the motion vector \code{c(u,v,w)} and the travel length until next interaction event.
    The movement of single photon is calculated in 3 dimensional space. The length of the straight section a photon can travel until the next interaction event is estimated as \code{-log( rnd )} where \code{rnd} is a uniformly distributed random variable from range 0-1.
}
\value{
    The \code{Move} function returns an MCBS object as \code{\link{class}}. Typically the same object is used in argument and result.
    As a result, the 3D position of the photon is updated:
    \item{x}{position along x-axis}
    \item{y}{position along y-axis}
    \item{z}{position along z-axis, the depth}
}
\note{
    Please note that this function is designed to be used in \code{\link{Simulation}} function. You may use it separately in case you design your own procedure. In such a case, please note that function \code{\link{Randomize}} prepares trajectory data for each photon and steps are done using move index \code{myObject$midx}.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{Chart}} for plot of calculated photon flux profile.
    \code{\link{Export}} for export of photon flux with corresponding radii.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Launch}} for start position of single photon and initial direction in media.
    \code{\link{Bounce}} for interaction with surface and computation of photon flux leaving media.
    \code{\link{Absorb}} for absorption of single photon energy in interaction with media.
    \code{\link{Scatter}} for single photon scattering interaction with media.
    \code{\link{Randomize}} for adjustment of random trajectory vectors.
}
\examples{
## Apple simulation data according to Qin and Lu (2006).
## DOI: 10.13031/2013.20862

# create object
cfgMedia <- c(0.63,30,0,1.4)
cfgSimulation <- c(1e6,0.05,1e-9,3,0.012)
apple <- MCBS(cfgMedia,cfgSimulation)
apple <- Setup(apple)
apple <- Randomize(apple)

# launch with 10° incident angle and move first photon
apple$idx <- 1
apple <- Launch(apple,10)
# first move
apple$midx <- 1
cat("Position",apple$x,":",apple$y,":",apple$z,"\n")
apple <- Move(apple)
cat("Position",apple$x,":",apple$y,":",apple$z,"\n")
# second move
apple$midx <- 2
apple <- Move(apple)
cat("Position",apple$x,":",apple$y,":",apple$z,"\n")

}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
