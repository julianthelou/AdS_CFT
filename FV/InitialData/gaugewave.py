from pylab import *

# Allows to just dump fortran stuff here and plot it 
# to understand the gaugewave

SIN=sin
COS=cos
SQRT=sqrt
ICA=0.1
Pi=pi
xGP=lambda d: x
x=linspace(0,1)
tGP=0

HH     = 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP))
dxH    = -2.0*Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP)) 
alp    = sqrt(HH)
Kxx    = - Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP))/SQRT( 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP))  )
traceK = Kxx/HH
phi    = ( 1.0 / HH)**(1.0/6.0)
GtildeX = 2.0/(3.0*HH**(5.0/3.0))*dxH
gammatXX = phi**2*HH

ion()
plot(x,HH)


