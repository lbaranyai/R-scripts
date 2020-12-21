\name{Export}
\alias{Export}
\alias{Export.MCBS}
\title{Export Simulation Result}
\description{
    This function exports simulation results in table format. It can be used to save results into file or make customized charts.
}
\usage{
Export(myObject)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
}
\details{
    The \code{Export} function, alias \code{Export.MCBS}, provides a table with column names. Results of Monte Carlo simulation computed by the \code{\link{Simulation}} function are given as \code{data.frame}. The radii are automatically calculated according to the spatial resolution.
}
\value{
    Results are provided as \code{data.frame} in table format. Two columns represent radii and photon flux.
    \item{Radius}{radius of concentric rings, where photon flux is calculated. Starts with incident point at radius of zero (r=0). The unit is cm.}
    \item{Flux}{calculated photon flux, 1/cm^2}
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{print}} show adjusted parameters and total backscattered reflection.
    \code{\link{Chart}} for plot of calculated photon flux profile.
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

# show photon flux near incident point
tbl <- Export(apple)
tbl[1:3,]

# save results into file
write.table(Export(apple),tempfile("apple",fileext=".dat"))
}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
