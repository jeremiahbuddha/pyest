#! /usr/bin/env python

# ==============================================================================
# Project Constants

# DRAG PARAMS
Cd = 2.0 # unitless
A = 3.6E-6 # km**2
m = 1350.0 # kg
rho_0 = 4.0E-4 # kg/km**3
ref_r = 7298.1450 # km
h_step = 200.0 # km
dTheta = 7.29211585530066E-5 # rad/s
drag_C = 0.5*Cd*(A/m)

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
        return (1 - 1.5*j2*(R/r(U))**2*(5*(U[2]/r(U))**2-1))
    elif comp =='z':
        return (1 - 1.5*j2*(R/r(U))**2*(5*(U[2]/r(U))**2-3))


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

