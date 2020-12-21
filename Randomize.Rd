\name{Randomize}
\alias{Randomize}
\alias{Randomize.MCBS}
\title{Random Numbers Selection for Single Photon Trajectory}
\description{
    This function prepares vectors of random numbers to accelerate trajectory computation. Random numbers are selected from uniform distribution.
}
\usage{
Randomize(myObject)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
}
\details{
    The \code{Randomize} function, alias \code{Randomize.MCBS}, puts random numbers into the container vectors created with the object. During trajectory computation, each interaction event occurs randomly. Prepared container vectors accelerate computation, since values for each event are accessible with \code{myObject$midx} move index. Vector length \code{myObject$MAXLEN} is adjusted by \code{\link{Setup}} function based on transport albedo and limiting energy.
}
\value{
    The \code{Randomize} function returns an MCBS object as \code{\link{class}}. Typically the same object is used in argument and result.
    There are five containers of random numbers. Values are in the range of 0 to +1 unless indicated in description.
    \item{rmv}{length of linear section between interaction events. Used by \code{\link{Move}} function.}
    \item{rabs}{decision values for roulette algorithm. Used by \code{\link{Absorb}} function.}
    \item{rx1}{new scattering direction in the range of -1 to +1. Used by \code{\link{Scatter}} function.}
    \item{rx1}{new scattering direction in the range of -1 to +1. Used by \code{\link{Scatter}} function.}
    \item{rmu}{variable for Heyney-Greenstein phase function. Used by \code{\link{Scatter}} function.}
}
\note{
    Please note that this function is designed to be used in \code{\link{Simulation}} function. You may use it separately in case you design your own procedure. In such a case, please note that function \code{\link{Randomize}} prepares random values for one photon and has to be run always before calculating single photon trajectory.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{Chart}} for plot of calculated photon flux profile.
    \code{\link{Export}} for export of photon flux with corresponding radii.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Move}} for moving single photon forward.
    \code{\link{Absorb}} for absorption of single photon energy in interaction with media.
    \code{\link{Scatter}} for single photon scattering interaction with media.
}
\examples{
## Apple simulation data according to Qin and Lu (2006).
## DOI: 10.13031/2013.20862

# create object
cfgMedia <- c(0.63,30,0,1.4)
cfgSimulation <- c(1e6,0.05,1e-9,3,0.012)
apple <- MCBS(cfgMedia,cfgSimulation)
apple <- Setup(apple)

# first five random values of move length
apple <- Randomize(apple)
apple$rmv[1:5]
}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
