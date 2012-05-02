#! /usr/bin/env python

import os
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from pyest.integ import Udot, Htilda_matrix
#from pyest.data import TRACKING_DICT, INITIAL_X0, INITIAL_x0, INITIAL_P0, W
from pyest.data import *
from pyest import *
from pyest.test import *

# =============================================================================
# =============================================================================
# MAIN EXECUTABLE SCRIPT FOR PYEST PROJECT 

T_END = 18340.0
T_DELTA = 20.0
OBS = TRACKING_DICT
I = matrix(eye(18))

# ==============================================================================
# FUNCTION DEFINITIONS

def initialize_batch(X0, P0, x0):
    """
    Generate t=0 values for a new iteration from an initial state, covariance
    and a-priori estimate.
    """
    # Get initial state and STM
    X0_list = X0.T.tolist()[0]
    STM0 = sp.matrix(sp.eye(18))
    STM0_list = sp.eye(18).reshape(1,324).tolist()[0]
    
    # Accumulate measurement at t=0
    obs0 = OBS[0]
    stn0 = obs0[0]
    comp0, Htilda0 = Htilda_matrix(X0_list, 0, stn0)
    ###comp0 = Htilda0 * X0
    resid0 = [ obs0[1] - float(comp0[0]),
               obs0[2] - float(comp0[1]) ]
    y0 = sp.matrix([resid0]).T
    H0 = Htilda0 * STM0

    # Sum matrices without apriori's (to check against solutions) 
    ### L0 = H0.T * W * H0
    ### N0 = H0.T * W * y0

    L0 = P0.I + H0.T * W * H0
    N0 = P0.I * x0 + H0.T * W * y0

    # Initialize the integrator
    eom = ode(Udot).set_integrator('dop853', atol=1.0E-10, rtol=1.0E-9)
    eom.set_initial_value(X0_list + STM0_list, 0)

    return [STM0, comp0, resid0, Htilda0, H0, L0, N0, eom]

def iterate_batch(X0, P0, x0):

    # Get t=0 parameters for this iteration
    IC = initialize_batch(X0, P0, x0)

    # Initialize container arrays - I am probably storing more values
    # than I really need to
    X      = { 0 : X0 }
    stm    = { 0 : IC[0] }
    comp   = { 0 : IC[1] }
    resids = { 0 : IC[2] }
    Htilda = { 0 : IC[3] }
    H      = { 0 : IC[4] }
    L      = IC[5]
    N      = IC[6]
    eom    = IC[7]

    # Integrate from t=0 to t=T_END
    while eom.successful() and eom.t < T_END:

        # Step state vector to current time
        this_t = eom.t + T_DELTA
        eom.integrate(this_t)

        # Split out state vector and STM for this time
        this_all = eom.y.tolist()

        this_X_list = this_all[0:18]
        this_X = sp.matrix(this_X_list).T
        X[this_t] = this_X
        
        this_stm_list = this_all[18:]
        this_stm = sp.matrix(this_stm_list).reshape(18,18)
        stm[this_t] = this_stm

        # If there is a measurement at the current time, then process it.
        if this_t in OBS.keys():

            this_obs = OBS[this_t]

            # Calculate predicted observations
            this_stn = OBS[this_t][0]
            this_comp, this_Htilda = Htilda_matrix(this_X, this_t, this_stn)
            this_H = this_Htilda * this_stm

            #this_comp = [float(this_Htilda[0] * this_X), # Range 
            #             float(this_Htilda[1] * this_X)] # Range-rate

            this_resid = [ this_obs[1] - this_comp[0],
                           this_obs[2] - this_comp[1]]

            this_y = sp.matrix([this_resid]).T

            # Accumulate matrices                       
            L = L + this_H.T * W * this_H
            N = N + this_H.T * W * this_y

            # Save values for bookkeeping
            comp[this_t]   = this_comp
            resids[this_t]  = this_resid
            Htilda[this_t] = this_Htilda
            H[this_t]      = this_H
                 

    # Invert the information matrix using cholesky decomposition
    L_ = np.linalg.cholesky(L)
    P = (L_.I).T * L_.I
    # Solve the normal equations for best estimate of X
    x_hat = P * N
 
#    print "### SUM: H.T * W * H"  
#    print L
#    print "### SUM: H.T * W * y" 
#    print N
#    exit()
 
    new_X = X0 + x_hat
    new_P = P
    new_x = x0 - x_hat

    return new_X, new_P, new_x, resids, [x_hat, X, stm, comp, resids, Htilda, H]

def initialize_sequential(X_bar0, P_bar0, x_bar0):
    """
    Generate t=0 values for a new iteration from an initial state, covariance
    and a-priori estimate.
    """
    # Get initial state and STM and initialize integrator
    X_bar0_list = X_bar0.T.tolist()[0]
    STM0 = sp.matrix(sp.eye(18))
    STM0_list = sp.eye(18).reshape(1,324).tolist()[0]

    eom = ode(Udot).set_integrator('dop853', atol=1.0E-10, rtol=1.0E-9)
    eom.set_initial_value(X0_list + STM0_list, 0)

    # Perform measurement update for t=0 observation
    obs0 = OBS[0]
    stn0 = obs0[0]
    y0, Htilda0 = Htilda_matrix(X_bar0_list, 0, stn0)
    comp0 = Htilda0 * X_bar0
    resid0 = [ obs0[1] - float(comp0[0]),
               obs0[2] - float(comp0[1]) ]
    y0 = sp.matrix([resid0]).T
    K0 = P_bar0 * Htilda0.T * (Htilda0 * P_bar0 * Htilda0.T + W.I).I
    x_hat0 = x_bar0 + K0 * (y0 - Htilda0 * x_bar0)
    P0 = (I - K0 * Htilda0) * P_bar0 

    return [STM0, comp0, resid0, Htilda0, x_hat0, P0, eom]


def iterate_sequential(X0, P0, x0):

    # Get t=0 parameters for this iteration
    IC = get_initial_params(X0, P0, x0)

    # Initialize container arrays
    X      = { 0 : X0 }
    stm    = { 0 : IC[0] }
    comp   = { 0 : IC[1] }
    resids = { 0 : IC[2] }
    Htilda = { 0 : IC[3] }
    x_hats = { 0 : IC[4] }
    Ps     = { 0 : IC[5] }
    eom    = IC[6]

    # Integrate from t=0 to t=T_END
    while eom.successful() and eom.t < T_END:

        # Step state vector to current time
        this_t = eom.t + T_DELTA
        eom.integrate(this_t)

        # Split out state vector and STM for this time
        this_all = eom.y.tolist()

        this_X_list = this_all[0:18]
        this_X = sp.matrix(this_X_list).T
        X[this_t] = this_X

        this_stm_list = this_all[18:]
        this_stm = sp.matrix(this_stm_list).reshape(18,18)
        stm[this_t] = this_stm

        # If there is a measurement at the current time, then process it.
        if this_t in OBS.keys():
 
            prev_t = this_t - T_DELTA

            this_obs = OBS[this_t]
          
            # Perform time update
            stm_n_nm1 = this_stm * stm[prev_t].I
            this_x_bar = stm_n_nm1 * x_hats[prev_t]
            this_P_bar = stm_n_nm1 * Ps[prev_t] * stm_n_nm1.T

            

            # Calculate predicted observations
            this_stn = OBS[this_t][0]
            this_Htilda = Htilda_matrix(this_X, this_t, this_stn)
            this_H = this_Htilda * this_stm

            this_comp = [float(this_Htilda[0] * this_X), # Range 
                         float(this_Htilda[1] * this_X)] # Range-rate

            this_resid = [ this_obs[1] - this_comp[0],
                           this_obs[2] - this_comp[1]]

            this_y = sp.matrix([this_resid]).T

            # Accumulate matrices                       
            L = L + this_H.T * W * this_H
            N = N + this_H.T * W * this_y



def plot_resids(resids, title = "Residuals (obs - com)"):

    # Compute the RMS for range and range-rate
    rng_rms = 0.0
    dRng_rms = 0.0
    for resid in resids.values():
        rng_rms += resid[0]**2
        dRng_rms += resid[1]**2
    rng_rms = sqrt(rng_rms / len(resids))
    dRng_rms = sqrt(dRng_rms / len(resids))

    # Create a figure
    resid_fig = plt.figure()

    # add sma plot
    rng_resids = [ float(res[0]) for res in resids.values() ]
    rng_ax  = resid_fig.add_subplot('211', title=title)
    rng_ax.plot(resids.keys(),rng_resids,'o')
    rng_ax.set_ylabel("Range (m)")
    rng_ax.annotate("Range RMS: {0}".format(rng_rms),
               xy=(10,-20), xycoords='axes points')

    # add ecc plot
    dRng_resids = [ float(res[1]) for res in resids.values() ]
    dRng_ax  = resid_fig.add_subplot('212')
    dRng_ax.plot(resids.keys(),dRng_resids,'o')
    dRng_ax.set_ylabel("Range-Rate (m/s)")
    dRng_ax.set_xlabel("Seconds from Epoch")
    dRng_ax.annotate("Range-rate RMS: {0}".format(dRng_rms), 
                xy=(10,-20), xycoords='axes points')

    # Show figure
    resid_fig.show()

    return rng_rms, dRng_rms

def print_state(X):
    """
    Pretty-print a matrix() 18x1 state.

    """
    out = ''
    for coord in range(18):
        out += "{0}".format(STATE_VARS[coord])
        val = float(X[coord])
        out += " {0: 2.4e}\n".format(val)

    print out

def print_covariance(P):
    """
    Pretty-print a matrix() 18x18 covariance.

    """
    def b(string):
        """
        Turns a given string blue.

        """
        return "\033[94m{0}\033[0m".format(string)

    out = "     "
    # Print out header with state variables
    for var in STATE_VARS:
        out += "{0:9s}  ".format(var)

    # Print out correlation / covariance matrix    
    for row in range(18):
        out += "\n{0:3s} ".format(STATE_VARS[row])
        for col in range(18):
            # Print correlations on lower diagnal
            if col < row:
                out += "{0: 2.2e}, ".format(float(P[row,col]/(sqrt(P[row,row]) * sqrt(P[col,col]) )))
            # Highlight variances in blue
            elif row == col:
                out += b("{0: 2.2e}, ".format(float(P[row,col])))
            else:
                out += "{0: 2.2e}, ".format(float(P[row,col]))
    
    out += "\n"

    print out


if __name__ == '__main__':

    # ==============================================================================
    # RUN THE BATCH FILTER

    # The first fit uses the apriori vals / covariance
    X1, P1, x1, resids1, TEST_DATA = iterate(INITIAL_X0, INITIAL_P0, INITIAL_x0)
    # Plot the residuals
    rng_rms, dRng_rms = plot_resids(resids1, title="Residuals (obs - com) for Iter 1")

    print "\n### x_hat for iteration 1:"
    print_state(TEST_DATA[0])
    print "\n### Cov for iteration 1:"
    print_covariance(P1)
    
    ## Second fit
    X2, P2, x2, resids2, TEST_DATA = iterate(X1, INITIAL_P0, x1)
    # Plot the residuals
    rng_rms, dRng_rms = plot_resids(resids2, title="Residuals (obs - com) for Iter 2")

    print "x_hat for iteration 2:"
    print_state(TEST_DATA[0])
    print "Cov for iteration 2:"
    print_covariance(P2)

    ## Third fit
    X3, P3, x3, resids3, TEST_DATA = iterate(X2, INITIAL_P0, x2)
    ## Plot the residuals
    rng_rms, dRng_rms = plot_resids(resids3, title="Residuals (obs - com) for Iter 3")

    print "x_hat for iteration 3:"
    print_state(TEST_DATA[0])
    print "Cov for iteration 3:"
    print_covariance(P3)

    raw_input("Press enter to exit")
