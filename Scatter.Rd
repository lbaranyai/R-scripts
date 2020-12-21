\name{Scatter}
\alias{Scatter}
\alias{Scatter.MCBS}
\title{Scattering of Single Photon in Media}
\description{
    Scattering is the interaction event in the media when photon hits structural elements and changes direction as a result of collision.
}
\usage{
Scatter(myObject)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
}
\details{
     The \code{Scatter} function, alias \code{Scatter.MCBS}, calculates the new direction of photon travel. In isotropic media, any new 3D direction can occur with equal probability. In anisotropic media, the new direction is decided according to the Heyney-Greenstein phase function. The isotropic and anisotropic behavior of media is shown by anisotropy factor (g). The anisotropy factor is zero (g = 0) for isotropic media, g = -1 means backward reflection (backscattering) is the most likely, while g = +1 means the most likely result is that photon continues forward (forward scattering). Biological tissues were reported to have high positive anisotropy factor g > 0.8.
}
\value{
    The \code{Scatter} function returns an MCBS object as \code{\link{class}}. Typically the same object is used in argument and result.
    As a result of \code{Scatter} function, new direction is given by \code{c(u,v,w)}. These values are used in function \code{\link{Move}} to calculate next position.
    \item{u}{3D vector coordinate along x-axis}
    \item{v}{3D vector coordinate along y-axis}
    \item{w}{3D vector coordinate along z-axis}
}
\references{
    Henyey, L.G., Greenstein, J.L. (1941)
    Diffuse Radiation in the Galaxy.
    \emph{Astrophysical Journal}, \bold{93}, 70--83.
    \doi{10.1086/144246}.

    Jacques, S.L. (1998)
    Light Distributions from Point, Line and Plane Sources for Photochemical Reactions and Fluorescence in Turbid Biological Tissues.
    \emph{Photochemistry and Photobiology}, \bold{67(1)}: 23-–32.
    \doi{10.1111/j.1751-1097.1998.tb05161.x}.
}
\note{
    Please note that this function is designed to be used in \code{\link{Simulation}} function. You may use it separately in case you design your own procedure. In such a case, please note that function \code{\link{Randomize}} prepares trajectory data for each photon and iterations are done using move index \code{myObject$midx}. Results of \code{Scatter} function are used in \code{\link{Move}} function.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{Chart}} for plot of calculated photon flux profile.
    \code{\link{Export}} for export of photon flux with corresponding radii.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Randomize}} for adjustment of random trajectory vectors.
    \code{\link{Launch}} for start position of single photon and initial direction in media.
    \code{\link{Bounce}} for interaction with surface and computation of photon flux leaving media.
    \code{\link{Move}} for moving single photon forward.
    \code{\link{Absorb}} for absorption of single photon energy in interaction with media.
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
apple <- Launch(apple)

# first interaction event
cat("Direction",apple$u,":",apple$v,":",apple$w,"\n")
apple$midx <- 1
apple <- Scatter(apple)
cat("Direction",apple$u,":",apple$v,":",apple$w,"\n")
}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
