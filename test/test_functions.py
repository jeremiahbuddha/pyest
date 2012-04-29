#! /usr/bin/env python

import unittest as ut
import scipy as sp
from pyest.integ import A_matrix, Htilda_matrix
from pyest.data import *
from pyest.test import *
from pyest.project_main import iterate

# =============================================================================
# Test Integ functions
# =============================================================================
class test_integ(ut.TestCase):

    def test_A_matrix(self):
        # Get initial state
        X0 = INITIAL_X0.T.tolist()[0]
        # Make the call to A_matrix function
        MY_A = A_matrix(X0)
        REF_A = sp.matrix(FIRST_CALL_TO_A)
        for row in range(18):
            for col in range(18):
                param = "{0}_dot wrt {1}".format(STATE_VARS[row],STATE_VARS[col])
                my_val = float(MY_A[row,col])
                ref_val = float(REF_A[row,col])
                ###print param+": {0}".format(abs(my_val - ref_val))
                self.assertAlmostEqual(my_val, ref_val, delta = min(1.E-7,abs(ref_val*1.E-7)))

    def test_Htilda_matrix(self):
        # Get initial state
        X0 = INITIAL_X0.T.tolist()[0]
        # Make the call to Htilda function
        MY_Htilda  = Htilda_matrix(X0,0,337)
        REF_Htilda = sp.matrix(FIRST_CALL_TO_HTILDA)
        for row in range(2):
            for col in range(18):
                param = "{0} wrt {1}".format(MEASUREMENT_VARS[row],STATE_VARS[col])
                my_val = float(MY_Htilda[row,col])
                ref_val = float(REF_Htilda[row,col])
                ###print abs(my_val - ref_val)
                self.assertAlmostEqual(my_val, ref_val, delta = min(1.E-7,abs(ref_val*1.E-7)))

    def test_iterate(self):

        # Use the iterate function to make first pass-through of data
        X1, P1, x1, resids1, TEST_VARS = iterate(INITIAL_X0, INITIAL_P0, INITIAL_x0)

        # Test state position and velocity at t=20 to 120
        MY_Xs = TEST_VARS[1]
        REF_Xs = P1_POS_VEL_T_20_120
        for time in [0,20,40,60,80,100,120]:
            my_pos_vel  = MY_Xs[time].T.tolist()[0]
            ref_pos_vel = REF_Xs[time]
            for coord in range(6):
                param = "t={0}, {1}".format(time,STATE_VARS[coord])
                my_val = my_pos_vel[coord]
                ref_val = ref_pos_vel[coord]
                ###print param+": {0}".format(abs(my_val - ref_val))
                self.assertAlmostEqual(my_val, ref_val, delta = min(1.E-3,abs(ref_val*1.E-3)))

        # Test state transition matrix at t=20 to 120
        MY_STMs = TEST_VARS[2]
        REF_STMs = P1_PHI_T_20_18340
        for time in [20, 18340]:
            my_stm  = MY_STMs[time]
            ref_stm = sp.matrix(REF_STMs[time])
            for row in range(18):
                for col in range(18):
                    param = "{0} wrt {1}0".format(STATE_VARS[row],STATE_VARS[col])
                    my_val = float(my_stm[row,col])
                    ref_val = float(ref_stm[row,col])
                    print param+": {0}".format(abs(my_val - ref_val))
                    #self.assertAlmostEqual(my_val, ref_val, delta = min(1.E-3,abs(ref_val*1.E-3)))
 

        ## TEST
        #MY_x_hat = TEST_VARS[0]
#
#        MY_X = TEST_VARS[1]
#        MY_stm = TEST_VARS[2]
#        MY_comp = TEST_VARS[3]
#        MY_resids = TEST_VARS[4]
#        MY_Htilda = TEST_VARS[5]
#        MY_H = TEST_VARS[6]




if __name__ == '__main__':
    ut.main()


