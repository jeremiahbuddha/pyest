
import scipy as sp
import numpy as np
from scipy import matrix
from scipy.integrate import ode
from pyest.filter import Udot, Htilda_matrix
from pyest.data import *

def initialize_batch(X_bar0, P_bar0, x_bar0):
    """
    Generate t=0 values for a new iteration from an initial state, covariance
    and a-priori estimate.
    """
    # Get initial state and STM and initialize integrator
    X_bar0_list = X_bar0.T.tolist()[0]
    stm0 = sp.matrix(sp.eye(18))
    stm0_list = sp.eye(18).reshape(1,324).tolist()[0]

    eom = ode(Udot).set_integrator('dop853', atol=1.0E-10, rtol=1.0E-9)
    eom.set_initial_value(X_bar0_list + stm0_list, 0)

    # Accumulate measurement at t=0
    obs0 = OBS[0]
    stn0 = obs0[0]
    comp0, Htilda0 = Htilda_matrix(X_bar0_list, 0, stn0)
    resid0 = [ obs0[1] - float(comp0[0]),
               obs0[2] - float(comp0[1]) ]
    y0 = sp.matrix([resid0]).T
    H0 = Htilda0 * stm0

    L0 = P_bar0.I + H0.T * W * H0
    N0 = P_bar0.I * x_bar0 + H0.T * W * y0

    return [stm0, comp0, resid0, Htilda0, H0, L0, N0, eom]

def iterate_batch(X_bar0, P_bar0, x_bar0):

    # Get t=0 parameters for this iteration
    IC = initialize_batch(X_bar0, P_bar0, x_bar0)

    # Initialize container arrays - I am probably storing more values
    # than I really need to
    Xs      = { 0 : X_bar0 }
    stms    = { 0 : IC[0] }
    comps   = { 0 : IC[1] }
    resids  = { 0 : IC[2] }
    Htildas = { 0 : IC[3] }
    Hs      = { 0 : IC[4] }
    L       = IC[5]
    N       = IC[6]
    eom     = IC[7]
    # Integrate from t=0 to t=T_END
    while eom.successful() and eom.t < T_END:

        # Step state vector to current time
        this_t = eom.t + T_DELTA
        eom.integrate(this_t)

        # Split out state vector and STM for this time
        this_all = eom.y.tolist()

        this_X_list = this_all[0:18]
        this_X = sp.matrix(this_X_list).T
        Xs[this_t] = this_X

        this_stm_list = this_all[18:]
        this_stm = sp.matrix(this_stm_list).reshape(18,18)
        stms[this_t] = this_stm

        # If there is a measurement at the current time, then process it.
        if this_t in OBS.keys():

            this_obs = OBS[this_t]

            # Calculate predicted observations
            this_stn = OBS[this_t][0]
            this_comp, this_Htilda = Htilda_matrix(this_X, this_t, this_stn)
            this_H = this_Htilda * this_stm

            this_resid = [ this_obs[1] - float(this_comp[0]),
                           this_obs[2] - float(this_comp[1])]

            this_y = sp.matrix([this_resid]).T

            # Accumulate matrices                       
            L = L + this_H.T * W * this_H
            N = N + this_H.T * W * this_y

            # Save values for bookkeeping
            comps[this_t]   = this_comp
            resids[this_t]  = this_resid
            Htildas[this_t] = this_Htilda
            Hs[this_t]      = this_H

    # Invert the information matrix using cholesky decomposition
    L_ = np.linalg.cholesky(L)
    P = (L_.I).T * L_.I
    # Solve the normal equations for best estimate of X
    x_hat = P * N

    # Compute postfit residuals
    presids = {}
    for time,resid in resids.items():
        offset= Hs[time] * x_hat
        presid = [ resid[0] - float(offset[0]),
                   resid[1] - float(offset[1])]
        presids[time] = presid

    # Update state and state deviation vectors
    new_X = X_bar0 + x_hat
    new_P = P
    new_x = x_bar0 - x_hat

    return [new_X, new_P, new_x, x_hat, resids, presids, Xs, stms, comps, Htildas, Hs]
