
import matplotlib.pyplot as plt
from math import sqrt, exp
from pyest.data import *

STATE_VARS= [ 'X  ', 'Y  ', 'Z  ', 'dX ', 'dY ', 'dZ ',
              'mu ', 'J_2', 'C_d', 'X_1', 'Y_1', 'Z_1',
              'X_2', 'Y_2', 'Z_2', 'X_3', 'Y_3', 'Z_3' ]

# ==============================================================================
# FUNCTION DEFINITIONS

def r(vec):
    return sqrt( float(vec[0])**2  +
                 float(vec[1])**2  +
                 float(vec[2])**2  )

def split_state(U):
    """
    Extracts X,Y,Z,dX,dY,dZ and magnitude r from the project
    state vector U.

    """
    return U[0], U[1], U[2], U[3], U[4], U[5], r(U)

def split_stn_coords(U, stn):
    """
    Evaluates drag parameters for a particular value of the
    state vector U.

    """
    if stn == 101:
        return U[9], U[10], U[11], 9
    elif stn == 337:
        return U[12], U[13], U[14], 12
    elif stn == 394:
        return U[15], U[16], U[17], 15

def acc_j2(U,comp):
    """
    This function augments the two-body EOMs with a J2 term.
    """
    if comp == 'x' or comp == 'y':
        return (1 - 1.5*J2*(R/r(U))**2*(5*(U[2]/r(U))**2-1))
    elif comp =='z':
        return (1 - 1.5*J2*(R/r(U))**2*(5*(U[2]/r(U))**2-3))

def rho_A(U):
    return rho_0*exp(-(r(U)-ref_r)/h_step)

def V_A(U):
    return sqrt( (U[3]+U[1]*dTheta)**2 +
                 (U[4]-U[0]*dTheta)**2 + U[5]**2 )

def acc_drag(U,comp):
    """
    This function augments the two-body EOMs with a drag term.

    """
    def drag_prefix(U):
        return -drag_C*rho_A(U)*V_A(U)
    if comp == 'x':
        return drag_prefix(U)*(U[3]+U[1]*dTheta)
    elif comp == 'y':
        return drag_prefix(U)*(U[4]-U[0]*dTheta)
    elif comp =='z':
        return drag_prefix(U)*(U[5])

def get_postfit_resids(run_data):
    # Generate postfit residuals
    x_hat = run_data[0]
    resids = run_data[4]
    Hs = run_data[6]
    presids = {}
    for time in range(0, int(T_END), int(T_DELTA)):
        if time in OBS.keys():
            obs = OBS[time]
            resid = resids[time]
            offset= Hs[time] * x_hat
            presid = [ resid[0] - float(offset[0]),
                       resid[1] - float(offset[1])]
            presids[time] = presid
    return presids

def plot_resids(iter_num, resids, presids, filter_type):

    # Create a figure
    resid_fig = plt.figure()

    # ======================================================================
    # Plot prefit resids
    # Compute the RMS for range and range-rate for prefit
    rng_rms = 0.0
    dRng_rms = 0.0
    for res in resids.values():
        rng_rms += res[0]**2
        dRng_rms += res[1]**2
    rng_rms = sqrt(rng_rms / len(resids))
    dRng_rms = sqrt(dRng_rms / len(resids))

    # add prefit rng
    rng_resids = [ float(res[0]) for res in resids.values() ]
    rng_ax  = resid_fig.add_subplot('221', 
                  title="{0} Prefit Residuals for Iter {1}".format(filter_type,iter_num) )
    rng_ax.plot(resids.keys(),rng_resids,'o')
    rng_ax.set_ylabel("Range (m)")
    rng_ax.annotate("Range RMS: {0}".format(rng_rms),
               xy=(10,-20), xycoords='axes points')

    # add prefit rng_rate
    dRng_resids = [ float(res[1]) for res in resids.values() ]
    dRng_ax  = resid_fig.add_subplot('223')
    dRng_ax.plot(resids.keys(),dRng_resids,'o')
    dRng_ax.set_ylabel("Range-Rate (m/s)")
    dRng_ax.set_xlabel("Seconds from Epoch")
    dRng_ax.annotate("Range-rate RMS: {0}".format(dRng_rms),
                xy=(10,-20), xycoords='axes points')

    # ======================================================================
    # Plot postfit resids
    # Compute the RMS for range and range-rate for postfit
    rng_rms = 0.0
    dRng_rms = 0.0
    for res in presids.values():
        rng_rms += res[0]**2
        dRng_rms += res[1]**2
    rng_rms = sqrt(rng_rms / len(presids))
    dRng_rms = sqrt(dRng_rms / len(presids))

    # add postfit rng
    rng_presids = [ float(res[0]) for res in presids.values() ]
    rng_ax  = resid_fig.add_subplot('222',
                  title="{0} Postfit Residuals for Iter {1}".format(filter_type,iter_num) )
    rng_ax.plot(presids.keys(),rng_presids,'o')
    rng_ax.set_ylabel("Range (m)")
    rng_ax.annotate("Range RMS: {0}".format(rng_rms),
               xy=(10,-20), xycoords='axes points')

    # add postfit rng_rate
    dRng_presids = [ float(res[1]) for res in presids.values() ]
    dRng_ax  = resid_fig.add_subplot('224')
    dRng_ax.plot(presids.keys(),dRng_presids,'o')
    dRng_ax.set_ylabel("Range-Rate (m/s)")
    dRng_ax.set_xlabel("Seconds from Epoch")
    dRng_ax.annotate("Range-rate RMS: {0}".format(dRng_rms),
                xy=(10,-20), xycoords='axes points')

    # Show figure
    resid_fig.show()

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

def print_latex_state(Xs):
    """
    Pretty-print a matrix() 18x1 state.

    """

    out = ''
    iters = Xs.keys()
    iters.sort()

    for num in iters:
        out += ' Iter {0} & '.format(num)

    out += '\\ \n'
 
    for coord in range(18):
        out += "${0}$ & ".format(STATE_VARS[coord])
        for num in iters:
            val = float(Xs[num][coord])
            out += " {0: 2.4e} &".format(val)
        out += '\\ \n'

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
                out += "{0: 2.2f},     ".format(float(P[row,col]/(sqrt(P[row,row]) * sqrt(P[col,col]) )))
            # Highlight variances in blue
            elif row == col:
                out += b("{0: 2.2e}, ".format(float(P[row,col])))
            else:
                out += "{0: 2.2e}, ".format(float(P[row,col]))
    
    out += "\n"

    print out

def plot_covariance_trace(Ps):

    # Get covariance times
    times = Ps.keys()
    pos_times = []
    pos_traces = []
    neg_times = []
    neg_traces = []
    for time in times:
        full_P =  Ps[time]
        pos_vel_P = full_P[0:6,0:6]

        trace = float(pos_vel_P.trace())
        if trace > 0:
            pos_times.append(time)
            pos_traces.append(trace)
        else:
            neg_times.append(time)
            neg_traces.append(abs(trace))

    # Create a figure
    fig = plt.figure()

    # add trace plot
    ax  = fig.add_subplot('111', title="Trace of covariance for Joseph")
    ax.plot(pos_times,pos_traces,'bo', label = 'Positive Vals')
    ax.plot(neg_times,neg_traces,'ro', label = 'Negative Vals')
    ax.set_xlabel("Time")
    ax.set_ylabel("Trace")
    ax.set_yscale('log')
    ax.legend()
 
    # Show figure
    fig.show()

# ==============================================================================
# Make submodules importable from pyest.filter namespace
from A_matrix import A_matrix
from Htilda_matrix import Htilda_matrix
from diff_state import Udot
from batch import initialize_batch, iterate_batch
from sequential import initialize_sequential, iterate_sequential

