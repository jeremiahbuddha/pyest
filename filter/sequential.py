
import scipy as sp
import numpy as np
from scipy import matrix
from scipy.integrate import ode
from pyest.filter import Udot, Htilda_matrix
from pyest.data import *

def initialize_sequential(X_bar0, P_bar0, x_bar0):
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

    # Perform measurement update for t=0 observation
    obs0 = OBS[0]
    stn0 = obs0[0]
    comp0, Htilda0 = Htilda_matrix(X_bar0_list, 0, stn0)
    resid0 = [ obs0[1] - float(comp0[0]),
               obs0[2] - float(comp0[1]) ]
    y0 = sp.matrix([resid0]).T
    K0 = P_bar0 * Htilda0.T * (Htilda0 * P_bar0 * Htilda0.T + W.I).I
    x_hat0 = x_bar0 + K0 * (y0 - Htilda0 * x_bar0)
    P0 = (I - K0 * Htilda0) * P_bar0
    #P0 = (I - K0 * Htilda0) * P_bar0 * (I - K0 * Htilda0).T + K0 * W.I * K0.T

    return [stm0, comp0, resid0, Htilda0, x_hat0, P0, eom]

def iterate_sequential(X_bar0, P_bar0, x_bar0):

    # Get t=0 parameters for this iteration
    IC = initialize_sequential(X_bar0, P_bar0, x_bar0)

    # Initialize container arrays
    Xs      = { 0 : X_bar0 }
    stms    = { 0 : IC[0] }
    comps   = { 0 : IC[1] }
    resids  = { 0 : IC[2] }
    Htildas = { 0 : IC[3] }
    x_hats  = { 0 : IC[4] }
    Ps      = { 0 : IC[5] }
    eom     = IC[6]

    prev_stm = IC[0]
    prev_x_hat = IC[4]
    prev_P = IC[5]
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
            this_stn = this_obs[0]

            # Perform time update
            stm_n_nm1 = this_stm * prev_stm.I
            this_x_bar = stm_n_nm1 * prev_x_hat
            this_P_bar = stm_n_nm1 * prev_P * stm_n_nm1.T

            this_comp, this_Htilda = Htilda_matrix(this_X, this_t, this_stn)
            this_resid = [ this_obs[1] - float(this_comp[0]),
                           this_obs[2] - float(this_comp[1])]

            this_y = sp.matrix([this_resid]).T

            this_K = this_P_bar * this_Htilda.T * (this_Htilda * this_P_bar * this_Htilda.T + W.I).I

            # Perform measurement update
            this_x_hat = this_x_bar + this_K * (this_y - this_Htilda * this_x_bar)
            this_P = (I - this_K * this_Htilda) * this_P_bar
            # Joseph update
            #this_P = (I - this_K * this_Htilda) * this_P_bar * (I - this_K * this_Htilda).T +\
            #         this_K * W.I * this_K.T

            prev_stm = this_stm
            prev_x_hat = this_x_hat
            prev_P = this_P

            # Save parameters
            comps[this_t] = this_comp
            resids[this_t] = this_resid
            Htildas[this_t] = this_Htilda
            x_hats[this_t] = this_x_hat
            Ps[this_t] = this_P

    x_hat = this_stm.I * this_x_hat
    P = this_stm.I * this_P * this_stm.I.T

    # Compute postfit residuals
    presids = {}
    for time,resid in resids.items():
        offset= Htildas[time] * x_hats[time]
        presid = [ resid[0] - float(offset[0]),
                   resid[1] - float(offset[1])]
        presids[time] = presid

    new_X = X_bar0 + x_hat
    new_P = P
    new_x = x_bar0 - x_hat

    return [new_X, new_P, new_x, x_hat, resids, presids, Xs, stms, comps, Htildas, Ps]

