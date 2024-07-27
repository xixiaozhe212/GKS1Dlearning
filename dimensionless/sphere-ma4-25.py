#======================================================#
# Description:                                         #
# Initial condition for sod shock tube problem         #
# Solver:                                              #
#                                                      #
#======================================================#

# ==================================================== #

from math import pi, sqrt, exp
import numpy as np

# ==================================================== #

def Tv2Kv(Tv):
    "Calculate Kv"
    ThetaV = 3371
    Ratio = ThetaV/Tv
    Kv = 2 * Ratio / (exp(Ratio) - 1)
    return Kv

# Boltzmann constant m^2 kg s^-2 K^-1
k = 1.38065e-23

# variable soft sphere (VSS) model
alpha = 1.0  # hard sphere
omega = 0.74

# Nitrogen molecular mass, kg
m = 46.5e-27
# Nitrogen molecular variable hard sphere (VHS) diameter, m
d = 4.17e-10
# Nitrogen gas specific gas constant
R = k/m

## Temperature, K
T = 65.0
Tref = 273

# Internal Degrees Of Freedom (iDoF)
Kv = Tv2Kv(T)
Kr = 0
iDoF = Kv + Kr
# Adiabatic constant
gamma = (5.0 + iDoF)/(3.0 + iDoF)

# characteristic length, m
L = 0.002

# Mach number
Ma = 4.25

## Velocity, m/2
U = Ma*sqrt(gamma*R*T)
# Most probable speed of molecular, m/2
C = sqrt(2.0*R*T)
# Knudsen number
Kn = 0.031
# Mean free path
mfp = L*Kn

# Reference Viscosity of VSS model
muRef = 5.0*(alpha+1.0)*(alpha+2.0)*sqrt(m*k*Tref/pi) \
      /(4.0*alpha*(5.0-2.0*omega)*(7.0-2.0*omega)*d*d)

# Inflow viscosity
mu = muRef*(T/Tref)**omega

# initial density field, Kg/m^3

# Inflow density by VSS model
rho = 4.0*alpha*( 5.0-2.0*omega )*( 7.0-2.0*omega ) / ( 5.0*( alpha+1 )*( alpha+2 ) ) \
    * sqrt( m/( 2.0*pi*kB*T) ) * mu / mfp

# Most probable speed of molecular, m/s
C = sqrt(2*R*T)

print("Density             (rho) = ", rho)
print("Temperature         (T)   = ", T)
print("Velocity            (U)   = ", U)
print("Viscosity           (mu)  = ", muRef)
print("Spefic gas constant (R)   = ", R)
print("Knudsen number      (Kn)  = ", Kn)
print("Mach number         (Ma)  = ", Ma)
print("Reynolds number     (Re)  = ", rho*U*L/mu)
print("Most probable speed (C)   = ", C)
print("Stop time           (t)   = ", 0.2*L/C)
