#! /usr/bin/env python

import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import ode
from pyest.integ import Udot, Htilda_matrix
from pyest.data import TRACKING_DICT
from pyest import *

# =============================================================================
# =============================================================================

# Set up the integrator (Runga-Kutta 8/5)
eom = ode(Udot).set_integrator('dop853', atol=1.0E-12, rtol=1.0E-12)

# Initalize the integrator
ref_state0 = state0 + stm0
eom.set_initial_value(ref_state0, 0)

 # Integrate over the time interval, saving a state at every time step
t_max = 18340.0
dt = 20.0
states = {}
observed = TRACKING_DICT
computed = {}
resids = {}
while eom.successful() and eom.t < t_max:
    # Step state vector to current time
    this_t = eom.t+dt
    eom.integrate(this_t)
    this_all = eom.y.tolist()
    
    # Save the current state
    states[this_t] = this_all

    if this_t in observed.keys():
        # Calculate predicted observations
        this_stn = observed[this_t][0]
        Htilda = Htilda_matrix(this_all, this_t, this_stn)

        state = this_all[0:18] 
        state_mtrx = sp.matrix(state).T

        computed[this_t] = [Htilda[0] * state_mtrx, # Range 
                            Htilda[1] * state_mtrx] # Range-rate

        resids[this_t] = [observed[this_t][1] - computed[this_t][0],
                          observed[this_t][2] - computed[this_t][1]]


# Print the State and STM at 18340 seconds
end_state = states[18340][:18]
end_STM = states[18340][18:]
#end_STM = sp.matrix(states[18340][18:]).reshape(18,18)

state_vars = [ 'X  ', 'Y  ', 'Z  ', 'dX ', 'dY ', 'dZ ',
               'mu ', 'J_2', 'C_d', 'X_1', 'Y_1', 'Z_1',
               'X_2', 'Y_2', 'Z_2', 'X_3', 'Y_3', 'Z_3' ]

print "\n### The state at t=18340"
for var,val in zip(state_vars,end_state):
    print "   {0}: {1}".format(var,val)

print "\n### The STM at t=18340"
row_start = range(0,162,18)
for indx in row_start:
    print "   {0:+0.4e},{1:+0.4e},{2:+0.4e},{3:+0.4e},{4:+0.4e},"+\
         "{5:+0.4e},{6:+0.4e},{7:+0.4e},{8:+0.4e}".format(
        end_STM[indx],end_STM[indx+1],end_STM[indx+2],
        end_STM[indx+3],end_STM[indx+4],end_STM[indx+5],
        end_STM[indx+6],end_STM[indx+7],end_STM[indx+8])


# Compute the RMS for range and range-rate
rng_rms = 0.0
dRng_rms = 0.0
for resid in resids.values():
    rng_rms += resid[0]**2
    dRng_rms += resid[1]**2
rng_rms = sqrt(rng_rms / len(resids))
dRng_rms = sqrt(dRng_rms / len(resids))

print "\n### The RMS values"
print "   Range: {0}".format(rng_rms)
print "   Range-rate: {0}".format(dRng_rms)

# Create a figure
resid_fig = plt.figure()

# add sma plot
rng_resids = [ float(res[0]) for res in resids.values() ]
rng_ax  = resid_fig.add_subplot('211', title="Residuals (obs - com)")
rng_ax.plot(resids.keys(),rng_resids,'o')
rng_ax.set_ylabel("Range (m)")

# add ecc plot
dRng_resids = [ float(res[1]) for res in resids.values() ]
dRng_ax  = resid_fig.add_subplot('212')
dRng_ax.plot(resids.keys(),dRng_resids,'o')
dRng_ax.set_ylabel("Range-Rate (m/s)")
dRng_ax.set_xlabel("Seconds from Epoch")

# Show figure
resid_fig.show()
raw_input("Press 'Enter' to exit")

