\name{Absorb}
\alias{Absorb}
\alias{Absorb.MCBS}
\title{Absorbtion of Single Photon Energy in Media}
\description{
    This function calculates how photon energy decreases during interaction with media.
}
\usage{
Absorb(myObject)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
}
\details{
    The \code{Absorb} function, alias \code{Absorb.MCBS}, calculates the new energy level after interaction event in the media. If photon energy falls below 0.001, a roulette algorithm makes decision when photon has 10\% chance to survive with higher energy and 90\% chance to loose all energy. Please keep in mind, that photons start with \code{weight = 1} but the first interaction happens at boundary, where photons enter the media.
}
\value{
    The \code{Absorb} function returns an MCBS object as \code{\link{class}}. Typically the same object is used in argument and result.
    One parameter of a single photon is updated in this procedure:
    \item{weight}{photon energy. Its initial value is 1 and decrease according to the albedo.}
}
\note{
    Please note that this function is designed to be used in \code{\link{Simulation}} function. You may use it separately in case you design your own procedure. In such a case, please note that function \code{\link{Randomize}} prepares trajectory data for each photon and iterations are done using move index \code{myObject$midx}. Function \code{Absorb} uses random data in roulette algorithm to make decision about photon survival at low energy level.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{Chart}} for plot of calculated photon flux profile.
    \code{\link{Export}} for export of photon flux with corresponding radii.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Launch}} for start position of single photon and initial direction in media.
    \code{\link{Bounce}} for interaction with surface and computation of photon flux leaving media.
    \code{\link{Move}} for moving single photon forward.
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
# initial weight in media
cat("Photon energy level",apple$weight,"\n")
# first interaction event
apple$midx <- 1
apple <- Absorb(apple)
cat("Photon energy level",apple$weight,"\n")

}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
