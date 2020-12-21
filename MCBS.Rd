\name{MCBS}
\alias{MCBS}
\title{Monte Carlo Backscattering class}
\description{
  This function is the constructor of Monte Carlo Backscattering (MCBS) class. Diffuse reflection, also called as backscattering, is the phenomena when light enters into the media, such as biological tissue, and after traveling some distance it leaves the media or tissue through the surface.
}
\usage{
MCBS(myTissue,myLight)
}
\arguments{
  \item{myTissue}{collection of optical parameters of the semi-infinite medium. It is a list of numbers: \code{c(mua,mus,g,n)} as absorption coefficient (1/cm), scattering coefficient (1/cm), anisortopy factor, refractive index.}
  \item{myLight}{collection of computational parameters. It is a list of numbers: \code{c(photons,beamradius,limitingenergy,haloradius,resolution)} as number of photons used in calculation, radius of incident light beam (cm), photon energy threshold to stop tracing if energy decreased too much, radius of surface area where photon flux is collected (cm), spatial resolution (cm/pixel). Due to the simulation of optical measurement, spatial resolution is expressed on the basis of pixel (px).}
}
\details{
    The \code{MCBS} function creates and initializes a class using the given parameters. The two arguments are separated according to their role. The first one describes media in which photons travel. The second defines simulation parameters. See arguments for details.
    
    The optical parameters describing observed media typically has the relationship of \code{mua << mus}, which means absorption coefficient is much smaller that scattering coefficient. Many research papers publish the reduced scattering coefficient \code{rmus = mus*(1-g)}. The general problem of comparison and validation is that \code{list(mus=2,g=0)} results the same reduced scattering coefficient like \code{list(mus=20,g=0.9)}. The value of g shows the most likely direction of light after interaction event in media: g = 1 means continue forward, g = 0 means there is no expected direction and media is isotropic, g = -1 means reflected backward. In biological tissues, anisotropy factor (g) is reported to be in the range of 0.8 - 0.99.
    
    Since argument \code{myLight} is reflecting measurement conditions and hardware setup, it can be considered constant for one system. Please note that quality of results depends on the number of photons. The more is used the better quality you get. The minimum suggested number is 100,000 photons. If you calculate the amount of photons with sensor or camera integration time, this number will increase above billion.
    The argument \code{myTissue} includes optical parameters of the sample, this is expected to change during experiment such as storage, grading, etc.
}
\value{
    The \code{MCBS} function returns an object of \code{\link{class}}. Package functions are used to access their parameters for reporting, adjustment and computation. Direct manipulation of object parameters is highly not recommended. Please use the provided functions to reach what is necessary or use them read-only.
    An object of class \code{"MCBS"} is a list containing at least the following components:
    \item{mu_a}{absorption coeffieint, 1/cm}
    \item{mu_s}{scattering coefficient, 1/cm}
    \item{g}{anisotropy factor}
    \item{n}{relative refractive index of media (the other media is assumed to be air)}
    \item{albedo}{transport albedo, the descrease of photon energy after one absorption event in media. The value is initialized by \code{\link{Setup}} function.}
    \item{cangle}{cosine of critical angle. The value is initialized by \code{\link{Setup}} function.}
    \item{rs}{specular reflection. The value is calculated on the basis of refractive index by \code{\link{Setup}} function.}
    Other parameters include containers for calculation. Please do not use them directly. See examples how to extract data.
}
\seealso{
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Randomize}} for adjustment of random trajectory vectors.
    \code{\link{Simulation}} for running the simulation with adjusted parameters.
    \code{\link{print}} show adjusted parameters and total backscattered reflection.
    \code{\link{Chart}} for plot of calculated photon flux profile.
    \code{\link{Export}} for export of photon flux with corresponding radii.
}
\references{
    Jacques, S.L. (1998)
    Light Distributions from Point, Line and Plane Sources for Photochemical Reactions and Fluorescence in Turbid Biological Tissues.
    \emph{Photochemistry and Photobiology}, \bold{67(1)}: 23-–32.
    \doi{10.1111/j.1751-1097.1998.tb05161.x}.

    Wang, L., Jacques, S.L., Zheng, L. (1995)
    MCML – Monte Carlo modeling of light transport in multi-layered tissues.
    \emph{Computer Methods and Programs in Biomedicine}, \bold{47}: 131-–146.
    \doi{10.1016/0169-2607(95)01640-F}.

    Wang, L., Jacques, S.L., Zheng, L. (1997)
    CONV-convolution for responses to a finite diameter photon beam incident on multi-layered tissues.
    \emph{Computer Methods and Programs in Biomedicine}, \bold{54}: 141–-150.
    \doi{10.1016/S0169-2607(97)00021-7}.

    Zołek, N.S., Liebert, A., Maniewski, R. (2006)
    Optimization of the Monte Carlo code for modeling of photon migration in tissue.
    \emph{Computer Methods and Programs in Biomedicine}, \bold{84}: 50–-57.
    \doi{10.1016/j.cmpb.2006.07.007}.

    Baranyai, L., Zude, M. (2009)
    Analysis of laser light propagation in kiwifruit using backscattering imaging and Monte Carlo simulation.
    \emph{Computers and Electronics in Agriculture}, \bold{69}: 33-–39.
    \doi{10.1016/j.compag.2009.06.011}.
    
    Baranyai, L. (2020)
    Laser induced diffuse reflectance imaging – Monte Carlo simulation of backscattering measured on the surface.
    \emph{MethodsX}, \bold{7}: 100958.
    \doi{10.1016/j.mex.2020.100958}.

}
\author{
    The object design was inspired by the C code "Small Monte Carlo" written by Scott Prahl (http://omlc.org).
}
\examples{
## Apple simulation data according to Qin and Lu (2006).
## DOI: 10.13031/2013.20862

cfgMedia <- c(
0.63,  # absorption, 670 nm
30,    # reduced scattering, 670 nm
0,     # isotropic tissue assumed
1.4)   # refractive index

cfgSimulation <- c(
1e6,   # 1 million photons
0.05,  # 1 mm diameter (0.05 cm radius) laser light beam
1e-9,  # limiting energy level
3,     # 3 cm radius is computed
0.012) # 0.012 cm/pixel resolution is assumed

apple <- MCBS(cfgMedia,cfgSimulation)

cfgMedia <- c(
0.63,  # absorption, 670 nm
300,   # scattering, 670 nm
0.9,   # anisotropic tissue with forward scattering is assumed
1.4)   # refractive index
apple <- MCBS(cfgMedia,cfgSimulation)
}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
