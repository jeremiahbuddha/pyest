#! /usr/bin/env python

from math import pi, sqrt, sin, cos, tan, acos, asin, atan, atan2
from scipy import matrix

def split_state(U):
    """
    Extracts X,Y,Z,dX,dY,dZ and magnitude r from the project
    state vector U.

    """
    return U[0], U[1], U[2], U[3], U[4], U[5], r(U)

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
    A = sp.zeros((18, 18), float)

    # ==========================================================================
    # Partial of dX wrt STATE
    # dX/dX
    A[0, 3] = 1.0

    # ==========================================================================
    # Partial of dY wrt STATE
    # dY/dY
    A[1, 4] = 1.0

    # ==========================================================================
    # Partial of dZ wrt STATE
    # dZ/dZ
    A[2, 5] = 1.0

    # ==========================================================================
    # Partial of ddX wrt STATE

    # ddX/X
    A[3, 0] = -mu / r**3 * (1 - 1.5 * J2 * (R / r)**2 * (5 * (Z / r)**2 - 1)) + \
              3 * mu * X**2 / r**5 * (1 - 2.5 * J2 * (R / r)**2 * (7 * (Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * X * (dX + dTheta * Y) / (r * H) + \
             -drag_C * rho_A * (-dTheta * dY + dTheta**2 * X) * (dX + dTheta * Y) / V_A

    # ddX/Y
    A[3, 1] = -3 * mu * X * Y / r**5 * (1 - 2.5 * J2 *(R / r)**2 * (7 *(Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * Y * (dX + dTheta * Y) / (r * H) + \
             -drag_C * rho_A * (dTheta * dX + dTheta**2 * Y) * (dX + dTheta * Y) / V_A + \
             -drag_C * rho_A * V_A * dTheta

    # ddX/Z
    A[3, 2] = 3 * mu * X * Z / r**5 * (1 - 2.5 * J2 *(R / r)**2*( 7 * (Z / r)**2 - 3)) + \
              drag_C * rho_A * V_A * Z * (dX + dTheta * Y) / (r * H)

    # ddX/dX
    A[3, 3] = -drag_C * rho_A * (dX + dTheta * Y)**2 / V_A + \
             -drag_C * rho_A * V_A

    # ddX/dY
    A[3, 4] = -drag_C * rho_A * (dY - dTheta * X) * (dX + dTheta * Y) / V_A 

    # ddX/dZ
    A[3, 5] = -drag_C * rho_A * dZ * (dX + dTheta * Y) / V_A 

    # ddX/mu
    A[3, 6] = -X / r**3 * (1 - 1.5 * J2 * (R / r)**2 * (5 * (Z / r)**2 - 1))

    # ddX/J2
    A[3, 7] = 1.5 * mu * X / r**3 * (R / r)**2 * (5 * (Z / r)**2 - 1)

    # ddX/Cd
    A[3, 8] = -0.5 * (A / m) * rho_A * V_A * (dX + dTheta * Y)

    # ==========================================================================
    # Partial of ddY wrt STATE

    # ddY/X
    A[4, 0] = 3 * mu * X * Y / r**5 * (1 - 2.5 * J2 *(R / r)**2 * (7 * (Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * X * (dY - dTheta * X) / (r * H) + \
             -drag_C * rho_A * (dTheta**2 * X - dTheta * dY) * (dY - dTheta * X) / V_A + \
              drag_C * rho_A * V_A * dTheta

    # ddY/Y
    A[4, 1] = -mu / r**3 * (1 - 1.5 * J2 * (R / r)**2 * (5 * (Z / r)**2 - 1)) + \
              3 * mu * Y**2 / r**5 * (1 - 2.5 * J2 * (R / r)**2 * (7 * (Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * Y * (dY - dTheta * X) / (r * H) + \
             -drag_C * rho_A * (dTheta * dX + dTheta**2 * Y) * (dY - dTheta * X) / V_A

    # ddY/Z
    A[4, 2] = 3 * mu * Y * Z / r**5 * (1 - 2.5 * J2 * (R / r)**2 * (7 * (Z / r)**2 - 3)) + \
              drag_C * rho_A * V_A * Z * (dY - dTheta * X) / (r * H)

    # ddY/dX
    A[4, 3] = -drag_C * rho_A * (dY - dTheta * X) * (dX + dTheta * Y) / V_A # FISHY same as for ddX/dY

    # ddY/dY
    A[4, 4] = -drag_C * rho_A * (dY - dTheta * X)**2 / V_A + \
             -drag_C * rho_A * V_A

    # ddY/dZ
    A[4, 5] = -drag_C * rho_A * dZ * (dY - dTheta * X) / V_A

    # ddY/mu
    A[4, 6] = -Y / r**3 * (1 - 1.5 * J2 * (R / r)**2 * (5 * (Z / r)**2 - 1))

    # ddY/J2
    A[4, 7] = 1.5 * mu * Y / r**3 * (R / r)**2 * (5 * (Z / r)**2 - 1)

    # ddY/Cd
    A[4, 8] = -0.5 * (A / m) * rho_A * V_A * (dY - dTheta * X)

    # ==========================================================================
    # Partial of ddZ wrt STATE

    # ddZ/X
    A[5, 0] = 3 * mu * X * Z / r**5 * (1 - 2.5 * J2 * (R / r)**2 * (7 * (Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * dZ * X * / (r * H) + \
             -drag_C * rho_A * dZ * (dTheta**2 * X - dTheta * dY) / V_A
              

    # ddZ/Y
    A[5, 1] = 3 * mu * Y * Z / r**5 * (1 - 2.5 * J2 * (R / r)**2 * (7 * (Z / r)**2 - 1)) + \
              drag_C * rho_A * V_A * dZ * Y * / (r * H) + \
             -drag_C * rho_A * dZ * (dTheta * dX + dTheta**2 * Y) / V_A


    # ddZ/Z
    A[5, 2] = -mu / r**3 * (1 - 1.5 * J2 * (R / r)**2 * (5 * (Z / r)**2 - 3)) + \
              3 * mu * Z**2 / r**5 * (1 - 2.5 * J2 * (R / r)**2 * (7 * (Z / r)**2 - 5)) + \
              drag_C * rho_A * V_A * Z * dZ / (r * H)

    # ddZ/dX
    A[5, 3] = -drag_C * rho_A * dZ * (dX + dTheta * Y) / V_A 

    # ddZ/dY
    A[5, 4] = -drag_C * rho_A * dZ * (dY - dTheta * X) / V_A 

    # ddZ/dZ
    A[5, 5] = -drag_C * rho_A * dZ**2 / V_A + \
             -drag_C * rho_A * V_A

    # ddZ/mu
    A[5, 6] = -Z / r**3 * (1 - 1.5 * J2 * (R / r)**2 * (5 * (Z / r)**2 - 3))

    # ddZ/J2
    A[5, 7] = 1.5 * mu * Z / r**3 * (R / r)**2 * (5 * (Z / r)**2 - 3)

    # ddZ/Cd
    A[5, 8] = -0.5 * (A / m) * rho_A * V_A * dZ

    return A

