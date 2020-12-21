\name{print}
\alias{print}
\alias{print.MCBS}
\title{Print text report on console}
\description{
    This function prints text report about simulation parameters and calculated coeeficients (albedo, critical angle) and total backscattered reflection.
}
\usage{
print(myObject)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
}
\details{
    The \code{print} function, alias \code{print.MCBS}, provides text description of simulation. Results of Monte Carlo simulation computed by the \code{\link{Simulation}} function are given as \code{data.frame}. Text report includes: absorption coefficient, scattering coefficient, anisotropy factor, refractive index, transport albedo, specular reflection, critical angle, backscattered total reflection.
}
\value{
    Results of simulation are provided as \code{data.frame} in table format. Print function lists adjusted parameters on console, does not return any value.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{Chart}} for plot of calculated photon flux profile.
    \code{\link{Export}} for export of photon flux with corresponding radii.
}
\examples{
## Apple simulation data according to Qin and Lu (2006).
## DOI: 10.13031/2013.20862

# create object
cfgMedia <- c(0.63,30,0,1.4)
# the 100 low number of photons is only for demonstration
cfgSimulation <- c(100,0.05,1e-9,3,0.1)
apple <- MCBS(cfgMedia,cfgSimulation)
apple <- Simulation(apple)

# show simulation report
print(apple)
}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
