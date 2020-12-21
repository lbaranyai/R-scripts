\name{Chart}
\alias{Chart}
\alias{Chart.MCBS}
\title{Plot Backscattering Photon Flux Profile}
\description{
    This function creates a plot of diffuse reflectance (backscattering) intensity measured on the surface. The photon flux is shown by radius (distance from incident point).
}
\usage{
Chart(myObject,isLog=FALSE)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
  \item{isLog}{Logical variable. Default value is FALSE. If set to TRUE, photon flux will be shown in logarithmic scale.}
}
\details{
    The \code{Chart} function, alias \code{Chart.MCBS}, plots simulation results. Photon flux is computed by \code{\link{Simulation}} function. The chart shows flux according to the spatial resolution of the simulation. In case you wish to make your customized chart, please use the \code{\link{Export}} function to receive data.
}
\value{
    This function does not have any result value.
}
\note{
    In case you wish to customize your chart, please use \code{\link{Export}} function to receive data.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{print}} show adjusted parameters and total backscattered reflection.
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

# show simulation result
Chart(apple)
# show simulation result with logarithmic scale
Chart(apple,TRUE)
}
\keyword{classes}
\keyword{methods}
\keyword{aplot}
