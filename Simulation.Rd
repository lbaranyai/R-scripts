\name{Simulation}
\alias{Simulation}
\alias{Simulation.MCBS}
\title{Runs the Monte Carlo Simulation}
\description{
    This function performs all actions to run the simulation. It calls other functions that prepare computation and calculate photon trajectory. Please keep in mind, that this technique can be slow depending on the number of photons and computational power of your hardware.
}
\usage{
Simulation(myObject,iAngle=0)
}
\arguments{
  \item{myObject}{The mandatory parameter of this function is the MCBS class created by \code{\link{MCBS}} constructor function.}
  \item{iAngle}{The optional parameter \code{iAngle} defines the incident angle of light beam. It is measured from normal position, perpendicular to the surface. The default value is zero (iAngle=0). This parameter is passed to \code{\link{Launch}} function.}
}
\details{
    The function \code{Simulation}, alias \code{Simulation.MCBS}, was designed to take every steps such as \code{Setup} initial conditions and compute each photon trajectory individually. It was written keeping clarity in mind. In case you want to make your own procedure, this function can be used as a template.
    
    The result of simulation is photon flux, which is expressed as relative fluence rate (F/P, 1/cm^2) where F is the fluence rate and P is the power of light source. The fluence rate (F) is computed in concentric rings around the incident point.
}
\value{
    The \code{Simulation} function returns an MCBS object as \code{\link{class}}. Typically the same object is used in argument and result.
    The main objective of the simulation is to calculate spatial photon flux. The flux is computed with radial averaging relative to the incident point. Each value represent the average flux of the ring with the width of resolution (1 pixel).
    \item{heat}{calculated photon flux, 1/cm^2}
}
\references{
    Farrell, T.J., Patterson, M.S., Wilson, B. (1992)
    A diffusion theory model of spatially resolved, steady-state diffuse reflectance for the noninvasive determination of tissue optical properties in vivo.
    \emph{Medical Physics}, \bold{19(4)}: 879--888.
    \doi{10.1118/1.596777}.

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
\note{
    Please note that quality of results depends on the number of photons. The more is used the better quality you get. The minimum suggested number is 100,000 photons. If you calculate the amount of photons with sensor or camera integration time, this number will increase above billion.
    Validation of computations is difficult, because other Monte Carlo programs or the diffusion theory model can be used as reference. The diffusion theory model does not consider many parameters included in Monte Carlo simulation, such as beam radius, incident angle, anisotropy. Please note that models and reports usually use the reduced scattering coefficient \code{rms = (1-g)*ms} what means that different scattering and anisotropy parameters, like \code{list(ms=20,g=0)} and \code{list(ms=200,g=0.9)}, result the same reduced scattering coefficient.
    Do not access and manipulate the result data vector \code{myObject$heat} directly. You can receive results by using the \code{\link{Export}} function.
}
\seealso{
    \code{\link{MCBS}} for construction of object with initial input parameters.
    \code{\link{Setup}} for initial computation of specular reflection and transport albedo.
    \code{\link{Randomize}} for adjustment of random trajectory vectors.
    \code{\link{print}} show adjusted parameters and total backscattered reflection.
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

# run simulation with default incident angle
apple <- Simulation(apple)
print(apple)

# run simulation with 10 deg incident angle
apple <- Simulation(apple,10)
Chart(apple)
}
\keyword{classes}
\keyword{methods}
\keyword{datagen}
\keyword{models}
