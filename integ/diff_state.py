#! /usr/bin/env python

from pyest import *
from pyest.integ import A_matrix
from scipy import matrix, array

# ==============================================================================
# ==============================================================================
# Define the differential state vector and stm

# NOTE: The equations for J2 and Drag accelerations are given in the solutions
# to ASEN5070 HW#2

def Udot(t,U):

    # Evaluate the differential state vector
    dState = [
        # State elements
        U[3],    # X_dot
        U[4],    # Y_dot
        U[5],    # Z_dot
       -mu*U[0]/r(U)**3*acc_j2(U,'x') + acc_drag(U,'x'), # DX_dot
       -mu*U[1]/r(U)**3*acc_j2(U,'y') + acc_drag(U,'y'), # DY_dot
       -mu*U[2]/r(U)**3*acc_j2(U,'z') + acc_drag(U,'z'), # DY_dot
        0.0,     # mu_dot
        0.0,     # J2_dot
        0.0,     # Cd_dot
        0.0,     # X1
        0.0,     # Y1
        0.0,     # Z1
        0.0,     # X2
        0.0,     # Y2
        0.0,     # Z2
        0.0,     # X3
        0.0,     # Y3
        0.0,     # Z3
        ]

    # Evaluate A-matrix at current value of U
    A = A_matrix(U)
    # Get the STM
    STM = matrix(U[18:342]).reshape((18,18))

    # Evaluate the differential STM
    dSTM = A*STM
    dSTM = dSTM.reshape((1,342)).tolist()[0]

    # Return the differential state + stm vector
    return dState + dSTM
