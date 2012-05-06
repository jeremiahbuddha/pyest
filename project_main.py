#! /usr/bin/env python

import sys
from pyest.filter import *
from pyest.data import *

# =============================================================================
# =============================================================================
# MAIN EXECUTABLE SCRIPT FOR PYEST PROJECT 

def run_batch(num_iters):

    # Initialize apriori state and state deviation vectors
    last_X = INITIAL_X0
    last_x = INITIAL_x0

    x_hats = {}
    Ps = {}
    # ==============================================================================
    # RUN THE BATCH FILTER
    for num in range(num_iters):
        
        # The first fit uses the apriori vals / covariance
        run_data = iterate_batch(last_X, INITIAL_P0, last_x)

        # NOTE:
        # run_data = [
        # 0   new_X,   Best estimate for state vector X
        # 1   new_P,   Covariance associated with new_X
        # 2   new_x,   Updated state deviation apriori for next iter
        # 3   x_hat,   Best estimate of state deviation 
        # 4   resids,  Dictionary of prefit residual values, indexed by time
        # 5   presids, Dictionary of postfit residual values, indexed by time
        # 6   Xs,      Dictionary of integrated states, indexed by time
        # 7   stms,    Dictionary of state transition matrices, indexed by time
        # 8   comps,   Dictionary of computed observation values, indexed by time
        # 9   Htildas, Dictionary of Htilda matrices, indexed by time
        # 10  Hs       Dictionary of H matrices, indexed by time
        #     ]

        plot_resids(num, run_data[4], run_data[5], 'Batch')

        #print "\n### Batch x_hat for iteration {0}:".format(num)
        #print_state(run_data[0])
        #print "\n### Batch Cov for iteration {0}:".format(num)
        #print_covariance(this_P)

        x_hats[num] = run_data[3]
        Ps[num] = run_data[1]       
 
        last_X = run_data[0]
        last_x = run_data[2]

    print_latex_state(x_hats)
    print_latex_cov(run_data[1])

def run_sequential(num_iters):

    # Initialize apriori state and state deviation vectors
    last_X = INITIAL_X0
    last_x = INITIAL_x0

    x_hats = {}
    Ps = {}
    # ==============================================================================
    # RUN THE SEQUENTIAL FILTER
    for num in range(num_iters):

        # The first fit uses the apriori vals / covariance
        run_data = iterate_sequential(last_X, INITIAL_P0, last_x)

        # NOTE:
        # run_data = [
        # 0   new_X,   Best estimate for state vector X
        # 1   new_P,   Covariance associated with new_X
        # 2   new_x,   Updated state deviation apriori for next iter
        # 3   x_hat,   Best estimate of state deviation 
        # 4   resids,  Dictionary of prefit residual values, indexed by time
        # 5   presids, Dictionary of postfit residual values, indexed by time
        # 6   Xs,      Dictionary of integrated states, indexed by time
        # 7   stms,    Dictionary of state transition matrices, indexed by time
        # 8   comps,   Dictionary of computed observation values, indexed by time
        # 9   Htildas, Dictionary of Htilda matrices, indexed by time
        # 10  Ps       Dictionary of covariance matrices, indexed by time
        #     ]

        plot_resids(num, run_data[4], run_data[5], 'Sequential')

        #print "\n### Sequential x_hat for iteration {0}:".format(num)
        #print_state(run_data[0])
        #print "\n### Sequential Cov for iteration {0}:".format(num)
        #print_covariance(this_P)

        x_hats[num] = run_data[3]
        Ps[num] = run_data[1]

        last_X = run_data[0]
        last_x = run_data[2]

    print_latex_state(x_hats)
    print_latex_cov(run_data[1])

if __name__ == '__main__':

    filter_name = sys.argv[1]
    num_iters = int(sys.argv[2])

    if filter_name.lower() == 'batch':
        run_batch(num_iters)
    elif filter_name.lower() == 'sequential':
        run_sequential(num_iters)
    else:
        print "Please select 'batch' or 'sequential' for filter type"
        exit()

    raw_input("Press enter to exit")
