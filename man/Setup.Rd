\name{Setup}
\alias{Setup}
\alias{Setup.MCBS}
\title{Initial Configuration of MCBS Class}
\description{
    This function performs the basic configuration of Monte Carlo Backscattering (\link{MCBS}) class. It calculates specular reflection, transport albedo, critical angle and similar parameters.
}
\usage{
Setup(myObject)
}
\arguments{
    \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
}
\details{
    The function \code{Setup}, alias \code{Setup.MCBS}, performs the initial calculation of essential parameters for computation. Additionally, computes the maximum expected trajectory length and prepares the launch position of all photons with random polar coordinates within circular beam.
}
\value{
    The \code{Setup} function returns an MCBS object as \code{\link{class}}. Typically the same object is used in argument and result.
    The new parameters calculated in \code{Setup} are:
    \item{albedo}{transport albedo}
    \item{rs}{specular reflection}
    \item{cangle}{cosine value of critical angle}
    \item{MAXLEN}{maximum expected length of photon trajectory. This length depends on the albedo and limiting energy.}
    There are additional containers for random launch position of photons. Launch position is given with polar coordinates (radius, angle).
}
\note{
    Please note that this function is designed to be used in \code{\link{Simulation}} function. You may use it separately in case you design your own procedure.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{Chart}} for plot of calculated photon flux profile.
    \code{\link{Export}} for export of photon flux with corresponding radii.
    \code{\link{Launch}} for start position of single photon and initial direction in media.
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
# transport albedo
apple$albedo
# specular reflection
apple$rs
}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
