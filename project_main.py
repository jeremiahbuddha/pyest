#! /usr/bin/env python

import scipy as sp
from scipy.integrate import ode
from pyest.integ import Udot, A_matrix, Htilda_matrix
from pyest import state0, stm0

# Set up the integrator (Runga-Kutta 8/5)
eom = ode(Udot).set_integrator('dop853', atol=1.0E-12, rtol=1.0E-12)

ref_state0 = state0 + stm0
# Initalize the integrator
eom.set_initial_value(ref_state0, 0)

 # Integrate over the time interval, saving a state at every time step
t_max = 18340.0
dt = 20.0
#while eom.successful() and eom.t < t_max:
eom.integrate(eom.t+t_max)
this_time = eom.t
this_all = eom.y.tolist()

#print 'State Vector:'
#print this_all[0:9]

STM_matrix = sp.matrix(this_all[18:]).reshape(18,18)

print 'STM:'
print STM_matrix

