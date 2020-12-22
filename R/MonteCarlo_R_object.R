## Monte Carlo simulation of photon trajectory in biological tissue
## by Laszlo Baranyai, lbaranyai@atb-potsdam.de

#
# Monte Carlo backscattering simulation object
#
# Initial parameters:
# myTissue <- c(mua,mus,g,n)
# myLight <- c(photons,beamradius,limitingenergy,haloradius,resolution)

## Constructor
MCBS <- function(myTissue,myLight)
{
  me <- list(
    # initial parameters
    mu_a    = myTissue[1], # absorption, 1/cm
    mu_s    = myTissue[2], # scattering, 1/cm
    g       = myTissue[3], # anisotropy
    n       = myTissue[4], # refractive index
    photons = myLight[1],  # number of photons
    beamr   = myLight[2],  # beam radius, cm
    limit   = myLight[3],  # limiting energy
    radius  = myLight[4],  # halo area radius, cm
    rpx     = myLight[5],  # resolution, cm/bin
    # containers and variables for calculation
    heat    = rep(0,1000),
    idx     = 1,
    x       = 0,
    y       = 0,
    z       = 0,
    u       = 0,
    v       = 0,
    w       = 0,
    weight  = 0,
    rs      = 0,
    rd      = 0,
    bit     = 0,
    albedo  = 0,
    cangle  = 0,
    MAXLEN  = 0,
    midx    = 0,
    rbr     = rep(0,1000),
    rba     = rep(0,1000),
    rmv     = rep(0,1000),
    rabs    = rep(0,1000),
    rx1     = rep(0,1000),
    rx2     = rep(0,1000),
    rmu     = rep(0,1000)
  )

  # set the name of the class
  class(me) <- append(class(me),"MCBS")
  return(me)
}

## Initialize computation
Setup <- function(myObject) { UseMethod("Setup",myObject) }
Setup.MCBS <- function(myObject)
{
  # container for intesity profile
  myObject$heat <- rep(0,1 + round(myObject$radius / myObject$rpx))
  # transport albedo = intensity decrease after interaction event
  myObject$albedo <- myObject$mu_s / (myObject$mu_s + myObject$mu_a)
  # specular reflection
  myObject$rs <- ( (myObject$n - 1) / (myObject$n + 1) )^2
  # critical angle
  myObject$cangle <- sqrt(1.0 - 1.0/(myObject$n^2))
  # maximum trajectory length in events
  myObject$MAXLEN <- round( log(myObject$limit)/log(myObject$albedo) )+1
  # beam start position, radius and angle
  myObject$rbr <- myObject$beamr*sqrt(runif(myObject$photons))
  myObject$rba <- runif(myObject$photons,min=0,max=2*pi)
  # reset summary values
  myObject$rd <- 0
  myObject$bit <- 0
  
  return(myObject)
}

## Launch single photon within beam
# Parameter:
# iAngle = incident angle, default = 0 deg (relative to normal)
Launch <- function(myObject,iAngle=0) { UseMethod("Launch",myObject) }
Launch.MCBS <- function(myObject,iAngle=0)
{
  # initial photon weight
  myObject$weight <- 1.0 - myObject$rs
  # internal angle after refraction
  myAngle <- asin(sin(iAngle*pi/180)/myObject$n)
  # incident direction
  myObject$u <- 0
  myObject$v <- sin(myAngle)
  myObject$w <- cos(myAngle)
  # start position
  myObject$x <- myObject$rbr[myObject$idx] * cos(myObject$rba[myObject$idx])
  myObject$y <- myObject$rbr[myObject$idx] * sin(myObject$rba[myObject$idx])/cos(myAngle)
  myObject$z <- 0

  return(myObject)
}

## Bounce interaction with surface
Bounce <- function(myObject) { UseMethod("Bounce",myObject) }
Bounce.MCBS <- function(myObject)
{
  myObject$w <- -1*myObject$w
  myObject$z <- -1*myObject$z
  # check for internal reflection, then nothing to do
  if (myObject$w > myObject$cangle) {
    t <- sqrt(1.0-(1.0-myObject$w^2)*myObject$n^2)
    temp1 <- (myObject$w - myObject$n*t)/(myObject$w + myObject$n*t)
    temp <- (t - myObject$n*myObject$w)/(t + myObject$n*myObject$w)
    # Fresnel reflection
    rf <- (temp1*temp1+temp*temp)/2.0
    myObject$rd <- myObject$rd + (1.0-rf) * myObject$weight
    # collect leaving photons by radius
    # Lambertian correction to normal direction
    lcc <- abs(myObject$w) / sqrt(myObject$u^2 + myObject$v^2 + myObject$w^2)
    lcc <- myObject$n * sqrt(1.0-lcc^2)
    if (lcc^2 > 1) {
      # failsafe check
      lcc <- 0
    } else {
      lcc <- sqrt(1.0 - lcc^2)
    }  
    # compute radius
    r <- sqrt(myObject$x^2 + myObject$y^2)
    r <- round(r / myObject$rpx) + 1
    if (r <= length(myObject$heat)) {
      myObject$heat[r] <- myObject$heat[r] + lcc * (1.0-rf) * myObject$weight
    }
    # continue travel inside
    myObject$weight <- myObject$weight - (1.0-rf) * myObject$weight;
  }
  
  return(myObject)
}

## Move photon within tissue
Move <- function(myObject) { UseMethod("Move",myObject) }
Move.MCBS <- function(myObject)
{
  # travel length
  d <- -1.0*log(myObject$rmv[myObject$midx])
  # new position
  myObject$x <- myObject$x + d * myObject$u
  myObject$y <- myObject$y + d * myObject$v
  myObject$z <- myObject$z + d * myObject$w

  return(myObject)
}

## Absorb energy in tissue
Absorb <- function(myObject) { UseMethod("Absorb",myObject) }
Absorb.MCBS <- function(myObject)
{
  # decrease energy due to absorption
  myObject$weight <- myObject$weight * myObject$albedo
  # roulette
  if (myObject$weight < 0.001) {
    myObject$bit <- myObject$bit - myObject$weight
    if (myObject$rabs[myObject$midx] > 0.1) {
      myObject$weight <- 0
    } else {
      myObject$weight <- 10*myObject$weight
    }
    myObject$bit <- myObject$bit + myObject$weight
  }

  return(myObject)
}

## Scatter light in tissue
Scatter <- function(myObject) { UseMethod("Scatter",myObject) }
Scatter.MCBS <- function(myObject)
{
  x1 <- myObject$rx1[myObject$midx]
  x2 <- myObject$rx2[myObject$midx]
  x3 <- x1^2 + x2^2
  if (myObject$g == 0) {
    # isotropic tissue
    if (x3<=0 || x3>1) {
      # computational problem, do nothing
      myObject$u <- 0
      myObject$v <- 0
      myObject$w <- 0
    } else {
     myObject$u <- 2.0 * x3 - 1.0
     myObject$v <- x1 * sqrt((1-myObject$u^2)/x3)
     myObject$w <- x2 * sqrt((1-myObject$u^2)/x3)
    }
  } else {
    # anisotropic tissue
    mu <- (1-myObject$g^2)/(1-myObject$g+2.0*myObject$g*myObject$rmu[myObject$midx])
    mu <- (1 + myObject$g^2-mu*mu)/2.0/myObject$g
    # failsafe check
    if (mu^2 > 1) { mu <- 1 }
    if ( myObject$w^2 < 0.81 ) {
      t <- mu * myObject$u + sqrt((1-mu*mu)/(1-myObject$w^2)/x3) * (x1*myObject$u*myObject$w-x2*myObject$v)
      myObject$v <- mu * myObject$v + sqrt((1-mu*mu)/(1-myObject$w^2)/x3) * (x1*myObject$v*myObject$w+x2*myObject$u)
      myObject$w <- mu * myObject$w - sqrt((1-mu*mu)*(1-myObject$w^2)/x3) * x1
    } else {
      t <- mu * myObject$u + sqrt((1-mu*mu)/(1-myObject$v^2)/x3) * (x1*myObject$u*myObject$v + x2*myObject$w)
      myObject$w <- mu * myObject$w + sqrt((1-mu*mu)/(1-myObject$v^2)/x3) * (x1*myObject$v*myObject$w - x2*myObject$u)
      myObject$v <- mu * myObject$v - sqrt((1-mu*mu)*(1-myObject$v^2)/x3) * x1;
    }
    myObject$u = t;
  }

  return(myObject)
}

## Fill random data for single photon trajectory
Randomize <- function(myObject) { UseMethod("Randomize",myObject) }
Randomize.MCBS <- function(myObject)
{
  # maximum trajectory length is computed to MAXLEN
  # random numbers are selected from uniform distribution
  # move length, 0-1
  myObject$rmv <- runif(myObject$MAXLEN)
  # absorption roulette, 0-1
  myObject$rabs <- runif(myObject$MAXLEN)
  # new direction after scattering, from -1 to +1
  myObject$rx1 <- runif(myObject$MAXLEN,min=-1,max=1)
  myObject$rx2 <- runif(myObject$MAXLEN,min=-1,max=1)
  # Heyney-Greenstein phase function random variable, 0-1
  myObject$rmu <- runif(myObject$MAXLEN)

  return(myObject)
}

## Run Monte Carlo simulation
# Parameter:
# iAngle = incident angle, default = 0 deg (relative to normal)
Simulation <- function(myObject,iAngle=0) { UseMethod("Simulation",myObject) }
Simulation.MCBS <- function(myObject,iAngle=0)
{
  # initialize computation
  myObject <- Setup(myObject)
  # compute all photons
  for (i in 1:myObject$photons) {
    # get random values
    # ACCELERATE myObject <- Randomize(myObject)
    myObject$rmv <- runif(myObject$MAXLEN)
    myObject$rabs <- runif(myObject$MAXLEN)
    myObject$rx1 <- runif(myObject$MAXLEN,min=-1,max=1)
    myObject$rx2 <- runif(myObject$MAXLEN,min=-1,max=1)
    myObject$rmu <- runif(myObject$MAXLEN)
    # EOF Randomize()
    # setup photon index
    myObject$idx <- i
    # follow trajectory of single photon
    # ACCELERATE myObject <- Launch(myObject,iAngle)
    myObject$weight <- 1.0 - myObject$rs
    myAngle <- asin(sin(iAngle*pi/180)/myObject$n)
    myObject$u <- 0
    myObject$v <- sin(myAngle)
    myObject$w <- cos(myAngle)
    myObject$x <- myObject$rbr[myObject$idx] * cos(myObject$rba[myObject$idx])
    myObject$y <- myObject$rbr[myObject$idx] * sin(myObject$rba[myObject$idx])/cos(myAngle)
    myObject$z <- 0
    # EOF Launch()
    myObject$midx <- 1
    while (myObject$weight > myObject$limit) {
      # ACCELERATE myObject <- Move(myObject)
      d <- -1.0*log(myObject$rmv[myObject$midx])
      myObject$x <- myObject$x + d * myObject$u
      myObject$y <- myObject$y + d * myObject$v
      myObject$z <- myObject$z + d * myObject$w
      # EOF Move()
      if (myObject$z <= 0) {
        # ACCELERATE myObject <- Bounce(myObject)
        myObject$w <- -1*myObject$w
        myObject$z <- -1*myObject$z
        if (myObject$w > myObject$cangle) {
          t <- sqrt(1.0-(1.0-myObject$w^2)*myObject$n^2)
          temp1 <- (myObject$w - myObject$n*t)/(myObject$w + myObject$n*t)
          temp <- (t - myObject$n*myObject$w)/(t + myObject$n*myObject$w)
          rf <- (temp1*temp1+temp*temp)/2.0
          myObject$rd <- myObject$rd + (1.0-rf) * myObject$weight
          lcc <- abs(myObject$w) / sqrt(myObject$u^2 + myObject$v^2 + myObject$w^2)
          lcc <- myObject$n * sqrt(1.0-lcc^2)
          if (lcc^2 > 1) {
            lcc <- 0
          } else {
            lcc <- sqrt(1.0 - lcc^2)
          }  
          r <- sqrt(myObject$x^2 + myObject$y^2)
          r <- round(r / myObject$rpx) + 1
          if (r <= length(myObject$heat)) {
            myObject$heat[r] <- myObject$heat[r] + lcc * (1.0-rf) * myObject$weight
          }
          myObject$weight <- myObject$weight - (1.0-rf) * myObject$weight;
        }
        # EOF Bounce()
      }
      # ACCELERATE myObject <- Absorb(myObject)
      myObject$weight <- myObject$weight * myObject$albedo
      if (myObject$weight < 0.001) {
        myObject$bit <- myObject$bit - myObject$weight
        if (myObject$rabs[myObject$midx] > 0.1) {
          myObject$weight <- 0
        } else {
          myObject$weight <- 10*myObject$weight
        }
        myObject$bit <- myObject$bit + myObject$weight
      }
      # EOF Absorb()
      # ACCELERATE myObject <- Scatter(myObject)
      x1 <- myObject$rx1[myObject$midx]
      x2 <- myObject$rx2[myObject$midx]
      x3 <- x1^2 + x2^2
      if (myObject$g == 0) {
        if (x3<=0 || x3>1) {
          myObject$u <- 0
          myObject$v <- 0
          myObject$w <- 0
        } else {
          myObject$u <- 2.0 * x3 - 1.0
          myObject$v <- x1 * sqrt((1-myObject$u^2)/x3)
          myObject$w <- x2 * sqrt((1-myObject$u^2)/x3)
        }
      } else {
        mu <- (1-myObject$g^2)/(1-myObject$g+2.0*myObject$g*myObject$rmu[myObject$midx])
        mu <- (1 + myObject$g^2-mu*mu)/2.0/myObject$g
        if (mu^2 > 1) { mu <- 1 }
        if ( myObject$w^2 < 0.81 ) {
          t <- mu * myObject$u + sqrt((1-mu*mu)/(1-myObject$w^2)/x3) * (x1*myObject$u*myObject$w-x2*myObject$v)
          myObject$v <- mu * myObject$v + sqrt((1-mu*mu)/(1-myObject$w^2)/x3) * (x1*myObject$v*myObject$w+x2*myObject$u)
          myObject$w <- mu * myObject$w - sqrt((1-mu*mu)*(1-myObject$w^2)/x3) * x1
        } else {
          t <- mu * myObject$u + sqrt((1-mu*mu)/(1-myObject$v^2)/x3) * (x1*myObject$u*myObject$v + x2*myObject$w)
          myObject$w <- mu * myObject$w + sqrt((1-mu*mu)/(1-myObject$v^2)/x3) * (x1*myObject$v*myObject$w - x2*myObject$u)
          myObject$v <- mu * myObject$v - sqrt((1-mu*mu)*(1-myObject$v^2)/x3) * x1;
        }
        myObject$u = t;
      }
      # EOF Scatter()
      myObject$midx <- myObject$midx + 1
      # if trajectory is too long, force quit
      if (myObject$midx > myObject$MAXLEN) myObject$weight <- 0;
    }
  }
  # calculate photon flux
  tmp <- seq(1,length(myObject$heat),by=1) - 1
  tmp <- pi*(2*tmp*myObject$rpx^2 + myObject$rpx^2) # ring area A = pi*((r+dr)^2 - r^2)
  myObject$heat <- myObject$heat/myObject$photons/tmp

  return(myObject)
}

## Print information about simulation
print <- function(myObject) { UseMethod("print",myObject) }
print.MCBS <- function(myObject)
{
  # show general message
  cat("Absorption =",myObject$mu_a,"1/cm\n")
  cat("Scattering =",myObject$mu_s,"1/cm\n")
  cat("Anisotropy =",myObject$g,"\n")
  cat("Refractive index =",myObject$n,"\n")
  if (myObject$rd == 0) {
    warning("No simulation was performed yet!")
  } else {
    # show derived attributes
    cat("Albedo =",myObject$albedo,"\n")
    cat("Specular reflection =",myObject$rs,"\n")
    cat("Critical angle =",acos(myObject$cangle)*180/pi,"deg\n")
    cat("Backscattered reflection =",myObject$rd/(myObject$bit + myObject$photons),"\n")
  }
}

## Draw chart of simulation result
Chart <- function(myObject,isLog=FALSE) { UseMethod("Chart",myObject) }
Chart.MCBS <- function(myObject,isLog=FALSE)
{
  # x axis data in cm
  x <- seq(1,length(myObject$heat),by=1)
  x <- (x - 1)*myObject$rpx
  if (isLog == TRUE) {
    plot(x,myObject$heat,log="y",col="blue",xlab="Radius, cm",ylab="Photon flux, 1/cm2")
  } else {
    plot(x,myObject$heat,col="blue",xlab="Radius, cm",ylab="Photon flux, 1/cm2")
  }
}

## Export intensity profile data
Export <- function(myObject) { UseMethod("Export",myObject) }
Export.MCBS <- function(myObject)
{
  # x axis data in cm
  x <- seq(1,length(myObject$heat),by=1)
  Radius <- (x - 1)*myObject$rpx
  Flux <- myObject$heat

  return( data.frame(Radius,Flux) )
}
