#! /usr/bin/env python

from math import pi, sqrt, sin, cos, tan, acos, asin, atan, atan2
from numpy import cross, dot
from scipy import matrix, array
from scipy.integrate import odeint, ode
import matplotlib.pyplot as plt
from ASEN5070.integ import State, Traj, planets, rad2deg, deg2rad


def split_state(U):
    return U[0],U[1],U[2],U[3],U[4],U[5], r(U)

def eval_drag(U):
    return rho_A(U), V_A(U)

def A_matrix(state):

    X,Y,Z,dX,dY,dZ,r = split_state(state)
    rho_A, V_A = eval_drag(state)

    # Create a matrix of zeros
    A = sp.zeros( (18,18),float )

### partial of dX wrt STATE
    # dX/dX
    A[0,3] = 1.0

### partial of dY wrt STATE

    # dY/dY
    A[1,4] = 1.0

### partial of dZ wrt STATE
    # dZ/dZ
    A[2,5] = 1.0

## partial of ddX wrt STATE
    # ddX/X
    A[3,0] = -mu/r**3*(1 - 1.5*j2*(R/r)**2*(5*(Z/r)**2-1)) + \
              3*mu*X**2/r**5*(1 - 2.5*j2*(R/r)**2*(7*(Z/r)**2-1)) +\
              drag_C*rho_A*V_A*X*(dX+dTheta*Y)/(r*H) +\
             -drag_C*rho_A*(-dTheta*dY + dTheta**2*X)*(dX+dTheta*Y)/V_A
    # ddX/Y
    A[3,1] = -3*mu*X*Y/r**5*(1 - 2.5*j2*(R/r)**2*(7*(Z/r)**2-1)) +\
              drag_C*rho_A*V_A*Y*(dX+dTheta*Y)/(r*H) +\
             -drag_C*rho_A*(dTheta*dX + dTheta**2*Y)*(dX+dTheta*Y)/V_A +\
             -drag_C*rho_A*V_A*dTheta
    # ddX/Z
    A[3,2] = 3*mu*X*Z/r**5*(1 - 2.5*j2*(R/r)**2*(7*(Z/r)**2-3)) +\
              drag_C*rho_A*V_A*Z*(dX+dTheta*Y)/(r*H)
    # ddX/dX
    A[3,3] = -drag_C*rho_A*(dX + dTheta*Y)**2/V_A +\
             -drag_C*rho_A*V_A
    # ddX/dY
    A[3,4] = -drag_C*rho_A*(dY - dTheta*X)*(dX + dTheta*Y)/V_A 
    # ddX/dZ
    A[3,5] = -drag_C*rho_A*dZ*(dX + dTheta*Y)/V_A 
    # ddX/mu
    A[3,6] = -X/r**3*(1 - 1.5*J2(R/r)**2*(5*(Z/r)**2 - 1))
    # ddX/J2
    A[3,7] = 1.5*mu*X/r**3*(R/r)**2*(5*(Z/r)**2 - 1)
    # ddX/Cd
    A[3,8] = -0.5*(A/m)*rho_A*V_A*(dX + dTheta*Y)

## partial of ddY wrt STATE
    # ddY/X
    A[4,0] = 3*mu*X*Y/r**5*(1 - 2.5*j2*(R/r)**2*(7*(Z/r)**2-1)) +\
              drag_C*rho_A*V_A*X*(dY - dTheta*X)/(r*H) +\
             -drag_C*rho_A*(dTheta**2*X - dTheta*dY)*(dY - dTheta*X)/V_A+\
              drag_C*rho_A*V_A*dTheta
    # ddY/Y
    A[4,1] = -mu/r**3*(1 - 1.5*j2*(R/r)**2*(5*(Z/r)**2-1)) + \
              3*mu*Y**2/r**5*(1 - 2.5*j2*(R/r)**2*(7*(Z/r)**2-1)) +\
              drag_C*rho_A*V_A*Y*(dY-dTheta*X)/(r*H) +\
             -drag_C*rho_A*(dTheta*dX + dTheta**2*Y)*(dY-dTheta*X)/V_A
    # ddY/Z
    A[4,2] = 3*mu*Y*Z/r**5*(1 - 2.5*j2*(R/r)**2*(7*(Z/r)**2-3)) +\
              drag_C*rho_A*V_A*Z*(dY - dTheta*X)/(r*H)
    # ddY/dX
    A[4,3] = -drag_C*rho_A*(dY - dTheta*X)*(dX + dTheta*Y)/V_A # FISHY same as for ddX/dY
    # ddY/dY
    A[4,4] = -drag_C*rho_A*(dY - dTheta*X)**2/V_A +\
             -drag_C*rho_A*V_A
    # ddY/dZ
    A[4,5] = -drag_C*rho_A*dZ*(dY - dTheta*X)/V_A
    # ddY/mu
    A[4,6] = -Y/r**3*(1 - 1.5*J2(R/r)**2*(5*(Z/r)**2 - 1))
    # ddY/J2
    A[4,7] = 1.5*mu*Y/r**3*(R/r)**2*(5*(Z/r)**2 - 1)
    # ddY/Cd
    A[4,8] = -0.5*(A/m)*rho_A*V_A*(dY - dTheta*X)


## partial of ddZ wrt STATE
    # ddZ/X
    A[5,0] = 3*mu*X*Z/r**5*(1 - 2.5*j2*(R/r)**2*(7*(Z/r)**2-1)) +\
              drag_C*rho_A*V_A*dZ*X*/(r*H) +\
             -drag_C*rho_A*dZ*(dTheta**2*X - dTheta*dY)/V_A
              
    # ddZ/Y
    A[5,1] = 3*mu*Y*Z/r**5*(1 - 2.5*j2*(R/r)**2*(7*(Z/r)**2-1)) +\
              drag_C*rho_A*V_A*dZ*Y*/(r*H) +\
             -drag_C*rho_A*dZ*(dTheta*dX + dTheta**2*Y)/V_A

    # ddZ/Z
    A[5,2] = -mu/r**3*(1 - 1.5*j2*(R/r)**2*(5*(Z/r)**2-3)) + \
              3*mu*Z**2/r**5*(1 - 2.5*j2*(R/r)**2*(7*(Z/r)**2-5)) +\
              drag_C*rho_A*V_A*Z*dZ/(r*H)
    # ddZ/dX
    A[5,3] = -drag_C*rho_A*dZ*(dX + dTheta*Y)/V_A 
    # ddZ/dY
    A[5,4] = -drag_C*rho_A*dZ*(dY - dTheta*X)/V_A 
    # ddZ/dZ
    A[5,5] = -drag_C*rho_A*dZ**2/V_A +\
             -drag_C*rho_A*V_A
    # ddZ/mu
    A[5,6] = -Z/r**3*(1 - 1.5*J2(R/r)**2*(5*(Z/r)**2 - 3))
    # ddZ/J2
    A[5,7] = 1.5*mu*Z/r**3*(R/r)**2*(5*(Z/r)**2 - 3)
    # ddZ/Cd
    A[5,8] = -0.5*(A/m)*rho_A*V_A*dZ

# ==============================================================================
# ==============================================================================
# Define the differential state vector and stm
# NOTE: The equations for J2 and Drag accelerations are given in the solutions
# to ASEN5070 HW#2

def Udot(t,U):
    dState = [
        # State elements
        U[3],    # X_dot
        U[4],    # Y_dot
        U[5],    # Z_dot
       -mu*U[0]/r(U)**3*acc_j2(U,'x') + acc_drag(U,'x'), # DX_dot
       -mu*U[1]/r(U)**3*acc_j2(U,'y') + acc_drag(U,'y'), # DY_dot
       -mu*U[2]/r(U)**3*acc_j2(U,'z') + acc_drag(U,'z'), # DY_dot
        0.0,     # mu_dot
        0.0,     # J2_dot
        0.0,     # Cd_dot
        0.0,     # X1
        0.0,     # Y1
        0.0,     # Z1
        0.0,     # X2
        0.0,     # Y2
        0.0,     # Z2
        0.0,     # X3
        0.0,     # Y3
        0.0,     # Z3
        ]

    A = A_matrix(U)
    STM = matrix(U[18:342]).reshape((18,18))

    dSTM = A*STM
    dSTM = dSTM.reshape((1,342)).tolist()[0]

    return dState + dSTM

# Integrate the EOM's
# Set up the integrator (Runga-Kutta 8/5)
eom = ode(Udot).set_integrator('dop853',
                               atol=1.0E-12,
                               rtol=1.0E-12)

eom.set_initial_value(ref_state0,0)

 # Integrate over the time interval, saving a state at every time step
t_max = 100.0
dt = 10.0
while eom.successful() and eom.t < t_max:
    eom.integrate(eom.t+dt)
    this_time = eom.t
    this_all = eom.y.tolist()
    this_state = this_all[0:4]
    this_stm = matrix( [ [ this_all[4],this_all[8], this_all[12],this_all[16] ],
                         [ this_all[5],this_all[9], this_all[13],this_all[17] ],
                         [ this_all[6],this_all[10],this_all[14],this_all[18] ],
                         [ this_all[7],this_all[11],this_all[15],this_all[19] ] ] )
    ref_states.append([this_time,this_state,this_stm])

for state in ref_states:
    print "\nt={0}, pos={1}, vel={2}, \nstm={3}".format(state[0],
                                           state[1][0:2],
                                           state[1][2:],
                                           state[2])
    
print """
# ==============================================================================
# c.) Demonstrate that the STM is symplectic by multiplying it by it's inverse
#     at every time stamp and confirming the result is Identity

"""
print ref_states[-1][2]*ref_states[-1][2].getI()


print """
# ==============================================================================
# d.) Demonstrate that the STM is symplectic by multiplying it by it's inverse
#     at every time stamp and confirming the result is Identity

"""
for i,state in enumerate(ref_states):
    print "\nX({0}) - X*({0}) = {1}".format(state[0],array(true_states[i][1])-array(state[1]))
    print "STM({0},0)*dX(0) = {1}".format(state[0],(state[2]*matrix([ [1.0E-6,-1.0E-6,1.0E-6,1.0E-6]] ).getT()).getT())


