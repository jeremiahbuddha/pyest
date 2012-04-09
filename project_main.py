#! /usr/bin/env python

import scipy as sp
from scipy.integrate import ode
from pyest.integ import Udot, Htilda_matrix
from pyest.data import TRACKING_DICT
from pyest import *

# Set up the integrator (Runga-Kutta 8/5)
eom = ode(Udot).set_integrator('dop853', atol=1.0E-12, rtol=1.0E-12)

# Initalize the integrator
ref_state0 = state0 + stm0
eom.set_initial_value(ref_state0, 0)

 # Integrate over the time interval, saving a state at every time step
t_max = 18340.0
dt = 20.0
observed = TRACKING_DICT
computed = {}
resids = {}
while eom.successful() and eom.t < t_max:
    # Step state vector to current time
    this_time = eom.t+dt
    eom.integrate(this_time)
    this_all = eom.y.tolist()

    if this_time in observed.keys():
        # Calculate predicted observations
        this_stn = observed[this_time][0]
        Htilda = Htilda_matrix(this_all, this_time, this_stn)

        state = this_all[0:18] 
        state_mtrx = sp.matrix(state).T

        computed[this_time] = [Htilda[0] * state_mtrx, # Range 
                               Htilda[1] * state_mtrx] # Range-rate

        resids[this_time] = [observed[this_time][1] - computed[this_time][0],
                             observed[this_time][2] - computed[this_time][1]]

STM_matrix = sp.matrix(this_all[18:]).reshape(18,18)
