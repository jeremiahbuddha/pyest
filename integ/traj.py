#! /usr/bin/env python

from math import pi, sqrt, sin, cos, tan, acos, asin, atan, atan2,exp
from numpy import array, cross, dot
from scipy import matrix, array
from scipy.integrate import odeint, ode
from ASEN5070.integ import State, planets, deg2rad, rad2deg, mag

# ============================================================================
# ============================================================================
# Traj class

class Traj(object):

    def __init__(self,initial_state,integ_len=None,integ_type="two_body", 
            step_size=20):
        self._state0 = initial_state
        # If integ_len isn't specified, set to two orbital periods
        if not integ_len:
            self._integ_len = 2*self._state0.get_period()
        else:
            self._integ_len = integ_len
        self._integ_type = integ_type.lower()
        self._step_size = step_size
        self._center = self._state0.get_center()
#        self._mu = planets[self._center][0]
        self._mu = 1.0
        self._j2 = planets[self._center][3]
        self._R  = planets[self._center][1]
        self._all_states = []
        self._integrate()

    def _get_EOM(self):
        """
        This function returns the two-body equations of motion, optionally
        including the effects of j2 and drag.

        """
        # Define contants used in the eom
        mu = self._mu
        j2 = self._j2
        R  = self._R
        
        def j2_suffix(U,comp):
            """
            This function augments the two-body EOMs with a J2 term.

            """
            if comp == 'x' or comp == 'y':
                return (1-1.5*j2*(R/mag(U))**2*(5*(U[2]/mag(U))**2-1))
            elif comp =='z':
                return (1-1.5*j2*(R/mag(U))**2*(5*(U[2]/mag(U))**2-3))

        def drag_suffix(U,comp):
            """
            This function augments the two-body EOMs with a J2 term.

            """
            Cd = 2.0 # unitless
            A = 3.6E-6 # km**2
            m = 1350.0 # kg
            rho_0 = 4.0E-4 # kg/km**3
            ref_r = 7298.1450 # km
            h_step = 200.0 # km
            theta_dot = 7.29211585530066E-5 # rad/s
            def rho_A(U):
                return rho_0*exp(-(mag(U)-ref_r)/h_step)
            def V_A(U):
                return sqrt( (U[3]+U[1]*theta_dot)**2 +
                             (U[4]-U[0]*theta_dot)**2 + U[5]**2 )
            def prefix(U):
                return -0.5*Cd*(A/m)*rho_A(U)*V_A(U)

            if comp == 'x':
                return prefix(U)*(U[3]+U[1]*theta_dot)
            elif comp == 'y':
                return prefix(U)*(U[4]-U[0]*theta_dot)
            elif comp =='z':
                return prefix(U)*(U[5])

        if self._integ_type == 'two_body':
            def two_body_EOM(t,U):
                return [ 
                    U[3], U[4], U[5], 
                    -mu*U[0]/mag(U)**3, 
                    -mu*U[1]/mag(U)**3, 
                    -mu*U[2]/mag(U)**3 
                    ]
            return two_body_EOM
           
        if self._integ_type == 'two_body_j2':

            def two_body_j2_EOM(t,U):
                return [
                    U[3], U[4], U[5],
                    -mu*U[0]/mag(U)**3*j2_suffix(U,'x'),
                    -mu*U[1]/mag(U)**3*j2_suffix(U,'y'),
                    -mu*U[2]/mag(U)**3*j2_suffix(U,'z'),
                    ]
            return two_body_j2_EOM
 
        if self._integ_type == 'two_body_drag':
            def two_body_drag_EOM(t,U):
                return [
                    U[3], U[4], U[5],
                    -mu*U[0]/mag(U)**3+drag_suffix(U,'x'),
                    -mu*U[1]/mag(U)**3+drag_suffix(U,'y'),
                    -mu*U[2]/mag(U)**3+drag_suffix(U,'z')
                    ]
            return two_body_drag_EOM

        if self._integ_type == 'all':
            def all_EOM(t,U):
                return [
                    U[3], U[4], U[5],
                    -mu*U[0]/mag(U)**3*j2_suffix(U,'x')+drag_suffix(U,'x'),
                    -mu*U[1]/mag(U)**3*j2_suffix(U,'y')+drag_suffix(U,'y'),
                    -mu*U[2]/mag(U)**3*j2_suffix(U,'z')+drag_suffix(U,'z')
                    ]
            return all_EOM



    def _integrate(self):

        # Reset state list
        self._all_states = []

        # Define our equations of motion
        Udot = self._get_EOM()
 
        # Set up the integrator (Runga-Kutta 8/5)
        eom = ode(Udot).set_integrator('dop853',
                                      atol=1.0E-12,
                                      rtol=1.0E-12)
        eom.set_initial_value(self._state0.get_cartesian(),0)

        # Integrate over the time interval, saving a state at every time step
        t_max = self._integ_len
        dt = self._step_size
        while eom.successful() and eom.t < t_max:
            eom.integrate(eom.t+dt)
            this_time = eom.t
            this_cartesian = eom.y.tolist()
            this_state = State(time=this_time, cartesian=this_cartesian)
            self._all_states.append(this_state)

    def get_states(self):
        return self._all_states

    def get_times(self):
        return [ state.get_time() for state in self._all_states ]

    def get_positions(self):
        return [ state.get_cartesian()[0:3] for state in self._all_states ]  

    def get_vels(self):
        return [ state.get_cartesian()[3:6] for state in self._all_states ]

    def get_accels(self):
        return [ [ -( self._mu/mag(state.get_cartesian()[0:3])**3 )*pos 
            for pos in state.get_cartesian()[0:3] ] 
            for state in self._all_states ]


