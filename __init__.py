#! /usr/bin/env python

import scipy as sp
from math import sqrt, exp
from pyest.data import *

# ==============================================================================
# Project Functions

def r(vec):
    return sqrt( float(vec[0])**2  +
                 float(vec[1])**2  +
                 float(vec[2])**2  )

def split_state(U):
    """
    Extracts X,Y,Z,dX,dY,dZ and magnitude r from the project
    state vector U.

    """
    return U[0], U[1], U[2], U[3], U[4], U[5], r(U)

def split_stn_coords(U, stn):
    """
    Evaluates drag parameters for a particular value of the
    state vector U.

    """
    if stn == 101:
        return U[9], U[10], U[11], 9
    elif stn == 337:
        return U[12], U[13], U[14], 12
    elif stn == 394:
        return U[15], U[16], U[17], 15

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

