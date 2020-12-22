\name{Launch}
\alias{Launch}
\alias{Launch.MCBS}
\title{Launch Photons into Semi-infinite Media}
\description{
    This function starts single photon into the semi-infinite media from surface (z = 0). The 3D launch direction \code{c(u,v,w)} points into the media, according to the incident angle (relative to normal) and refractive index. If incident angle was zero, the incidence is normal (perpendicular to surface). If incidence angle was not zero, elliptical distortion of circular light beam is also considered.
}
\usage{
Launch(myObject,iAngle=0)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
  \item{iAngle}{This optional parameter means the incident angle relative to normal. The default value is 0 (the incidence is perpendicular to surface).}
}
\details{
    The \code{Launch} function, alias \code{Launch.MCBS}, calculates the start position of the single photon on the surface. The coordinates \code{c(x,y,z)} are computed as z = 0 and \code{c(x,y)} represent random position in circular light beam. The values z > 0 are inside the media. Beam intensity distribution is considered to be flat, what means that positions have equal probability to occur. Starting polar coordinates are adjusted randomly using uniform distribution by function \code{\link{Setup}}. If incident angle was not zero, elliptical distortion is corrected along the y-axis.
}
\value{
    The \code{Launch} function returns an MCBS object as \code{\link{class}}. Typically the same object is used in argument and result.
    Position and starting 3D vector are calculated:
    \item{x}{position along x-axis, calculated from random polar coordinates.}
    \item{y}{position along y-axis, calculated from random polar coordinates. If incident angle was not zero, this value is distorted.}
    \item{z}{initial depth, z=0 on the surface.}
    \item{u}{launch vector coordinate along x-axis.}
    \item{v}{launch vector coordinate along y-axis. If incident angle was not zero, this value is calculated according to the incident angle and refractive index.}
    \item{w}{launch vector coordinate along z-axis. If incident angle was not zero, this value is calculated according to the incident angle and refractive index.}
}
\references{
    Jacques, S.L. (1998)
    Light Distributions from Point, Line and Plane Sources for Photochemical Reactions and Fluorescence in Turbid Biological Tissues.
    \emph{Photochemistry and Photobiology}, \bold{67(1)}: 23–-32.
    \doi{10.1111/j.1751-1097.1998.tb05161.x}.

    Wang, L., Jacques, S.L., Zheng, L. (1997)
    CONV-convolution for responses to a finite diameter photon beam incident on multi-layered tissues.
    \emph{Computer Methods and Programs in Biomedicine}, \bold{54}: 141–-150.
    \doi{10.1016/S0169-2607(97)00021-7}.
}
\note{
    Theoretical models assume the incident angle is zero. During measurement, camera and lightsource cannot occupy the same direction. According to the Fresnel reflection equations, incident angle below 20° is recommended. The non-zero incident angle requires correction in computation by means of refraction of the initial direction inside media and elliptical distortion of circular light beam on the surface.
    Please note that this function is not intended to use independently. It calculates start values for photon of index \code{myObject$idx}. If you want to use it to build your custom simulation, do not forget to assign index, such as \code{myObject$idx <- 1}. Index number must be between 1 and \code{myObject$MAXLEN}.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{Chart}} for plot of calculated photon flux profile.
    \code{\link{Export}} for export of photon flux with corresponding radii.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Bounce}} for interaction with surface and computation of photon flux leaving media.
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

# calculate first photon
apple$idx <- 1
apple <- Launch(apple)
# see vector coordinates
cat("Start position",apple$x,":",apple$y,"\n")
cat("Direction",apple$u,":",apple$v,":",apple$w,"\n")

# calculate first photon with 10° incident angle
apple$idx <- 1
apple <- Launch(apple,10)
# see vector coordinates
cat("Start position",apple$x,":",apple$y,"\n")
cat("Direction",apple$u,":",apple$v,":",apple$w,"\n")
}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
