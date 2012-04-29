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
                #if ref_val != 0.0:
                #    print param+": {0} %".format(abs(my_val - ref_val)/abs(ref_val) * 100.0)
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
                #if ref_val != 0.0:
                #    print param+": {0} %".format(abs(my_val - ref_val)/abs(ref_val) * 100.0)
                self.assertAlmostEqual(my_val, ref_val, delta = min(1.E-7,abs(ref_val*1.E-7)))

    def test_iterate(self):

        # Use the iterate function to make first pass-through of data
        X1, P1, x1, resids1, TEST_VARS = iterate(INITIAL_X0, INITIAL_P0, INITIAL_x0)
 
        # NOTE: TEST_VARS = [x_hat, X, stm, comp, resids, Htilda, H]

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
                #if ref_val != 0.0:
                #    print param+": {0} %".format(abs(my_val - ref_val)/abs(ref_val) * 100.0)
                self.assertAlmostEqual(my_val, ref_val, delta = min(1.E-3,abs(ref_val*1.E-3)))

        # Test state transition matrix at t=20 to 120
        print "\nNOT TESTING STM"
        MY_STMs = TEST_VARS[2]
        REF_STMs = P1_PHI_T_20_18340
        for time in [20, 18340]:
            my_stm  = MY_STMs[time]
            ref_stm = sp.matrix(REF_STMs[time])
            for row in range(18):
                for col in range(18):
                    param = "t={0}, {1} wrt {2}0".format(time, STATE_VARS[row], STATE_VARS[col])
                    my_val = float(my_stm[row,col])
                    ref_val = float(ref_stm[row,col])
                    #if ref_val != 0.0:
                    #    print param+": {0} %".format(abs(my_val - ref_val)/abs(ref_val) * 100.0)
                    #self.assertAlmostEqual(my_val, ref_val, delta = min(1.E-3,abs(ref_val*1.E-3)))

        # Test H matrices
        print "\nNOT TESTING H"
        MY_Hs  = TEST_VARS[6]
        REF_Hs = P1_H_T_0_20_18340
        for time in [0, 20, 18340]:
            my_H  = MY_Hs[time]
            ref_H = sp.matrix(REF_Hs[time])
            for row in range(2):
                for col in range(18):
                    param = "t={0}, {1} wrt {2}".format(time, MEASUREMENT_VARS[row], STATE_VARS[col])
                    my_val = float(my_H[row,col])
                    ref_val = float(ref_H[row,col])
                    #if ref_val != 0.0:
                    #    print param+": {0} %".format(abs(my_val - ref_val)/abs(ref_val)*100.0)
                    #self.assertAlmostEqual(my_val, ref_val, delta = min(1.E-7,abs(ref_val*1.E-7)))


        # Test x_hat after first pass through
        print "\nNOT TESTING x_hat"
        MY_xhat = TEST_VARS[0]
        REF_xhat = P1_XHAT
        for coord in range(18):
            param = "{0}".format(STATE_VARS[coord])
            my_val = float(MY_xhat[coord])
            ref_val = REF_xhat[coord]
            #print param+": {0:2.4f} %".format(abs(my_val - ref_val)/abs(ref_val)*100.0)
            #self.assertAlmostEqual(my_val, ref_val, delta = min(1.E-7,abs(ref_val*1.E-7)))


if __name__ == '__main__':
    ut.main()
