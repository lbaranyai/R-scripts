\name{Bounce}
\alias{Bounce}
\alias{Bounce.MCBS}
\title{Bounce Photon Interaction with Surface}
\description{
    This function is applied when photon trajectory goes out of the media or reaches the surface (z <= 0). The photon flux is summarized upon leaving semi-infinite media through the surface.
}
\usage{
Bounce(myObject)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
}
\details{
    The \code{Bounce} function, alias \code{Bounce.MCBS}, calculates internal reflection if incident angle to surface boundary was higher than critical angle. In all other cases, Fresnel reflection is calculated according to the direction and refractive index of media. Remaining photon weight that leaves the media is collected and summarized by the position. Position is calculated as distance from incident point. This way, radial averaging is performed with rigns of width of the spatial resolution. The internally reflected part of photon energy still follows trajectory until reaches limiting energy or expected trajectory length.
}
\value{
    The \code{Bounce} function returns an MCBS object as \code{\link{class}}. Typically the same object is used in argument and result.
    The function collects photon flux in an array:
    \item{heat}{vector of photon flux scattered back from media. One element represent a ring of width of the spatial resolution. The center of the ring is the incident point.}
}
\note{
    Please note thet this is an internal function of \code{\link{Simulation}}. Do not use it separately, unless you build your own simulation procedure. The function uses the coordinates \code{c(x,y,z)} of the current photon and the \code{c(u,v,w)} direction of its last move.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{Chart}} for plot of calculated photon flux profile.
    \code{\link{Export}} for export of photon flux with corresponding radii.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Launch}} for start position of single photon and initial direction in media.
    \code{\link{Move}} for moving single photon forward.
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

# launch and move first photon
apple$idx <- 1
apple <- Launch(apple)
# first move
apple$midx <- 1
apple <- Move(apple)
if (apple$z <= 0) {
 apple <- Bounce(apple)
}

}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
