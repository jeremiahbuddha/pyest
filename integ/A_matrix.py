#! /usr/bin/env python

from math import pi, sqrt, sin, cos, tan, acos, asin, atan, atan2
from scipy import matrix
from pyest.data import *
from pyest import *

def eval_drag(U):
    """
    Evaluates drag parameters for a particular value of the
    state vector U.

    """
    return rho_A(U), V_A(U)

# ==============================================================================
# Define the partial derivative A-Matrix for the projects state

def A_matrix(U):
    """
    Evaluate the projects A-matrix for a particular value of the
    state vector U.

    """
    # Get the state elements in a more readable form
    X, Y, Z, dX, dY, dZ, r = split_state(U)
    rho_A, V_A = eval_drag(U)

    # Most of the A-matrix are zero's. Define a matrix of all zero's here, and
    # then fill in the non-zero terms below.
    Atrix = sp.zeros((18, 18), float)

    # ==========================================================================
    # Partial of dX wrt STATE
    # dX/dX
    Atrix[0, 3] = 1.0

    # ==========================================================================
    # Partial of dY wrt STATE
    # dY/dY
    Atrix[1, 4] = 1.0

    # ==========================================================================
    # Partial of dZ wrt STATE
    # dZ/dZ
    Atrix[2, 5] = 1.0

    # ==========================================================================
    # Partial of ddX wrt STATE

    # ddX/X
    Atrix[3, 0] = -mu / r**3 * (1 - (3.0/2.0) * J2 * (R / r)**2 * (5 * (Z / r)**2 - 1)) + \
              3 * mu * X**2 / r**5 * (1 - (5.0/2.0) * J2 * (R / r)**2 * (7 * (Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * X * (dX + dTheta * Y) / (r * h_step) + \
             -drag_C * rho_A * (-dTheta * dY + dTheta**2 * X) * (dX + dTheta * Y) / V_A

    # ddX/Y
    Atrix[3, 1] = 3 * mu * X * Y / r**5 * (1 - (5.0/2.0) * J2 *(R / r)**2 * (7 * (Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * Y * (dX + dTheta * Y) / (r * h_step) + \
             -drag_C * rho_A * (dTheta * dX + dTheta**2 * Y) * (dX + dTheta * Y) / V_A + \
             -drag_C * rho_A * V_A * dTheta

    # ddX/Z
    Atrix[3, 2] = 3 * mu * X * Z / r**5 * (1 - (5.0/2.0) * J2 *(R / r)**2 * ( 7 * (Z / r)**2 - 3)) + \
              drag_C * rho_A * V_A * Z * (dX + dTheta * Y) / (r * h_step)

    # ddX/dX
    Atrix[3, 3] = -drag_C * rho_A * (dX + dTheta * Y)**2 / V_A + \
             -drag_C * rho_A * V_A

    # ddX/dY
    Atrix[3, 4] = -drag_C * rho_A * (dY - dTheta * X) * (dX + dTheta * Y) / V_A 

    # ddX/dZ
    Atrix[3, 5] = -drag_C * rho_A * dZ * (dX + dTheta * Y) / V_A 

    # ddX/mu
    Atrix[3, 6] = -X / r**3 * (1 - (3.0/2.0) * J2 * (R / r)**2 * (5 * (Z / r)**2 - 1))

    # ddX/J2
    Atrix[3, 7] = (3.0/2.0) * mu * X / r**3 * (R / r)**2 * (5 * (Z / r)**2 - 1)

    # ddX/Cd
    Atrix[3, 8] = -(1.0/2.0) * (A / m) * rho_A * V_A * (dX + dTheta * Y)

    # ==========================================================================
    # Partial of ddY wrt STATE

    # ddY/X
    Atrix[4, 0] = 3 * mu * X * Y / r**5 * (1 - (5.0/2.0) * J2 * (R / r)**2 * (7 * (Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * X * (dY - dTheta * X) / (r * h_step) + \
             -drag_C * rho_A * (dTheta**2 * X - dTheta * dY) * (dY - dTheta * X) / V_A + \
              drag_C * rho_A * V_A * dTheta

    # ddY/Y
    Atrix[4, 1] = -mu / r**3 * (1 - (3.0/2.0) * J2 * (R / r)**2 * (5 * (Z / r)**2 - 1)) + \
              3 * mu * Y**2 / r**5 * (1 - (5.0/2.0) * J2 * (R / r)**2 * (7 * (Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * Y * (dY - dTheta * X) / (r * h_step) + \
             -drag_C * rho_A * (dTheta * dX + dTheta**2 * Y) * (dY - dTheta * X) / V_A

    # ddY/Z
    Atrix[4, 2] = 3 * mu * Y * Z / r**5 * (1 - (5.0/2.0) * J2 * (R / r)**2 * (7 * (Z / r)**2 - 3)) + \
              drag_C * rho_A * V_A * Z * (dY - dTheta * X) / (r * h_step)

    # ddY/dX
    Atrix[4, 3] = -drag_C * rho_A * (dY - dTheta * X) * (dX + dTheta * Y) / V_A

    # ddY/dY
    Atrix[4, 4] = -drag_C * rho_A * (dY - dTheta * X)**2 / V_A + \
             -drag_C * rho_A * V_A

    # ddY/dZ
    Atrix[4, 5] = -drag_C * rho_A * dZ * (dY - dTheta * X) / V_A

    # ddY/mu
    Atrix[4, 6] = -Y / r**3 * (1 - (3.0/2.0) * J2 * (R / r)**2 * (5 * (Z / r)**2 - 1))

    # ddY/J2
    Atrix[4, 7] = (3.0/2.0) * mu * Y / r**3 * (R / r)**2 * (5 * (Z / r)**2 - 1)

    # ddY/Cd
    Atrix[4, 8] = -(1.0/2.0) * (A / m) * rho_A * V_A * (dY - dTheta * X)

    # ==========================================================================
    # Partial of ddZ wrt STATE

    # ddZ/X
    Atrix[5, 0] = 3 * mu * X * Z / r**5 * (1 - (5.0/2.0) * J2 * (R / r)**2 * (7 * (Z / r)**2 - 3)) + \
              drag_C * rho_A * V_A * dZ * X / (r * h_step) + \
             -drag_C * rho_A * dZ * (dTheta**2 * X - dTheta * dY) / V_A
              

    # ddZ/Y
    Atrix[5, 1] = 3 * mu * Y * Z / r**5 * (1 - (5.0/2.0) * J2 * (R / r)**2 * (7 * (Z / r)**2 - 3)) + \
              drag_C * rho_A * V_A * dZ * Y / (r * h_step) + \
             -drag_C * rho_A * dZ * (dTheta * dX + dTheta**2 * Y) / V_A


    # ddZ/Z
    Atrix[5, 2] = -mu / r**3 * (1 - (3.0/2.0) * J2 * (R / r)**2 * (5 * (Z / r)**2 - 3)) + \
              3 * mu * Z**2 / r**5 * (1 - (5.0/2.0) * J2 * (R / r)**2 * (7 * (Z / r)**2 - 5)) + \
              drag_C * rho_A * V_A * Z * dZ / (r * h_step)

    # ddZ/dX
    Atrix[5, 3] = -drag_C * rho_A * dZ * (dX + dTheta * Y) / V_A 

    # ddZ/dY
    Atrix[5, 4] = -drag_C * rho_A * dZ * (dY - dTheta * X) / V_A 

    # ddZ/dZ
    Atrix[5, 5] = -drag_C * rho_A * dZ**2 / V_A + \
             -drag_C * rho_A * V_A

    # ddZ/mu
    Atrix[5, 6] = -Z / r**3 * (1 - (3.0/2.0) * J2 * (R / r)**2 * (5 * (Z / r)**2 - 3))

    # ddZ/J2
    Atrix[5, 7] = (3.0/2.0) * mu * Z / r**3 * (R / r)**2 * (5 * (Z / r)**2 - 3)

    # ddZ/Cd
    Atrix[5, 8] = -(1.0/2.0) * (A / m) * rho_A * V_A * dZ

    return Atrix

