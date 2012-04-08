#! /usr/bin/env python

import scipy as sp
from math import sqrt, exp

# ==============================================================================
# Project Constants

# Drag Params
Cd = 2.0 # unitless
A = 3.0 # m**2
m = 970.0 # kg
rho_0 = 3.614E-13 # kg / m**3
ref_r = 7078136.3 # m
h_step = 88667.0 # m
drag_C = 0.5 * Cd * (A / m)

# Other constants
dTheta = 7.2921158553E-5 # rad / s
J2 = 1.082626925638815E-3
mu = 3.986004415E+14  # m**3 / s**2
R = 6378136.3 # m

# Tracking station coordinates (ECF)
stn101 = [-5127510.0,-3794160.0, 0.0]
stn337 = [3860910.0, 3238490.0, 3898094.0]
stn394 = [549505.0,-1380872.0, 6182197.0]

# Initial State (units in meters and secons)
state0 = [757700.0,  # X
          5222607.0, # Y
          4851500.0, # Z
          2213.21,   # dX
          4678.34,   # dY 
         -5371.30,   # dZ
          mu,        # mu
          J2,        # J2
          Cd,        # Cd
          stn101[0], # X_1
          stn101[1], # Y_1
          stn101[2], # Z_1
          stn337[0], # X_1
          stn337[1], # Y_1
          stn337[2], # Z_1
          stn394[0], # X_1
          stn394[1], # Y_1
          stn394[2], # Z_1
          ]

# Initial STM is just identity, arranged into a list
stm0 = sp.identity(18, float).reshape(1,324).tolist()[0]

# ==============================================================================
# Project Functions

def r(vec):
    return sqrt( float( vec[0]**2 ) +
                 float( vec[1]**2 ) +
                 float( vec[2]**2 ) )


def acc_j2(U,comp):
    """
    This function augments the two-body EOMs with a J2 term.
    """
    if comp == 'x' or comp == 'y':
        return (1 - 1.5*J2*(R/r(U))**2*(5*(U[2]/r(U))**2-1))
    elif comp =='z':
        return (1 - 1.5*J2*(R/r(U))**2*(5*(U[2]/r(U))**2-3))


def rho_A(U):
    return rho_0*exp(-(r(U)-ref_r)/h_step)

def V_A(U):
    return sqrt( (U[3]+U[1]*dTheta)**2 +
                 (U[4]-U[0]*dTheta)**2 + U[5]**2 )

def acc_drag(U,comp):
    """
    This function augments the two-body EOMs with a drag term.

    """
    def drag_prefix(U):
        return -drag_C*rho_A(U)*V_A(U)
    if comp == 'x':
        return drag_prefix(U)*(U[3]+U[1]*dTheta)
    elif comp == 'y':
        return drag_prefix(U)*(U[4]-U[0]*dTheta)
    elif comp =='z':
        return drag_prefix(U)*(U[5])

