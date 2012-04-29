#! /usr/bin/env python

from math import pi, sqrt, sin, cos, tan, acos, asin, atan, atan2
from scipy import matrix
from pyest.data import *
from pyest import *

# ==============================================================================
# Define the partial derivative H~-Matrix for the projects state

def Htilda_matrix(U,t,stn):
    """
    Evaluate the projects H~-matrix for a particular value of the
    state vector U.

    """
    # Calculate theta at current time
    theta = dTheta * t

    # Get the state elements in a more readable form
    X, Y, Z, dX, dY, dZ, r = split_state(U)

    # Get station coords 
    Xs, Ys, Zs, stn_indx = split_stn_coords(U, stn)

    # Calulate range
    rng = sqrt(X**2 + Y**2 + Z**2 + Xs**2 + Ys**2 + Zs**2 \
          - 2 * (X * Xs + Y * Ys) * cos(theta) + 2*(X * Ys - Y * Xs) * sin(theta) \
          - 2 * Z * Zs)

    # Calulate range-rate
    dRng = (X * dX + Y * dY + Z * dZ - (dX * Xs + dY * Ys) * cos(theta) + \
           dTheta * (X * Xs + Y * Ys) * sin(theta) + \
           (dX * Ys - dY * Xs) * sin(theta) + \
           dTheta * (X * Ys - Y * Xs) * cos(theta) - dZ * Zs) / rng

    # Most of the Htilda-matrix are zero's. Define a matrix of all zero's here, and
    # then fill in the non-zero terms below.
    Htilda = sp.zeros((2, 18), float)

    # ==========================================================================
    # Partial of range wrt STATE

    # rng/X
    Htilda[0, 0] = (X - Xs * cos(theta) + Ys * sin(theta)) / rng

    # rng/Y
    Htilda[0, 1] = (Y - Ys * cos(theta) - Xs * sin(theta)) / rng

    # rng/Z
    Htilda[0, 2] = (Z - Zs) / rng

    # rng/Xs
    Htilda[0, stn_indx] = (Xs - X * cos(theta) - Y * sin(theta)) / rng

    # rng/Ys
    Htilda[0, stn_indx + 1] = (Ys - Y * cos(theta) + X * sin(theta)) / rng

    # rng/Zs
    Htilda[0, stn_indx + 2] = (Zs - Z) / rng


    # ==========================================================================
    # Partial of Range-rate wrt STATE

    # dRng/X
    Htilda[1, 0] = (dX + dTheta * Xs * sin(theta) + dTheta * Ys * cos(theta)) / rng \
                   - dRng * (X - Xs * cos(theta) + Ys * sin(theta)) / rng**2                    

    # dRng/Y
    Htilda[1, 1] = (dY + dTheta * Ys * sin(theta) - dTheta * Xs * cos(theta)) / rng \
                   - dRng * (Y - Ys * cos(theta) - Xs * sin(theta)) / rng**2

    # dRng/Z
    Htilda[1, 2] = dZ / rng - dRng * (Z - Zs) / rng**2

    # dRng/dX
    Htilda[1, 3] = (X - Xs * cos(theta) + Ys * sin(theta)) / rng

    # dRng/dY
    Htilda[1, 4] = (Y - Ys * cos(theta) - Xs * sin(theta)) / rng

    # dRng/dZ
    Htilda[1, 5] = (Z - Zs) / rng

    # dRng/Xs
    Htilda[1, stn_indx] = (-dX * cos(theta) + dTheta * X * sin(theta) \
                          - dY * sin(theta) - dTheta * Y * cos(theta)) / rng \
                          - dRng * (Xs - X * cos(theta) - Y * sin(theta)) / rng**2

    # dRng/Ys
    Htilda[1, stn_indx + 1] = (-dY * cos(theta) + dTheta * Y * sin(theta) \
                          + dX * sin(theta) + dTheta * X * cos(theta)) / rng \
                          - dRng * (Ys - Y * cos(theta) + X * sin(theta)) / rng**2

    # dRng/Zs
    Htilda[1, stn_indx + 2] = -dZ / rng - dRng * (Zs - Z) / rng**2


    return Htilda

