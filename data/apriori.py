
from scipy import matrix, eye


APVALS = matrix( [      
    0.0E+00, # x
    0.0E+00, # y
    0.0E+00, # z 
    0.0E+00, # dx
    0.0E+00, # dy
    0.0E+00, # dz
    0.0E+00, # mu
    0.0E+00, # J2
    0.0E+00, # Cd
    0.0E+00, # x1
    0.0E+00, # y1
    0.0E+00, # z1
    0.0E+00, # x2
    0.0E+00, # y2
    0.0E+00, # z2 
    0.0E+00, # x3
    0.0E+00, # y3
    0.0E+00, # z3 
    ]).T

x0 = APVALS

APCOV = matrix( [ 
    1.0E+06, # x
    1.0E+06, # y
    1.0E+06, # z 
    1.0E+06, # dx
    1.0E+06, # dy
    1.0E+06, # dz
    1.0E+20, # mu
    1.0E+06, # J2
    1.0E+06, # Cd
    1.0E-10, # x1
    1.0E-10, # y1
    1.0E-10, # z1
    1.0E+06, # x2
    1.0E+06, # y2
    1.0E+06, # z2 
    1.0E+06, # x3
    1.0E+06, # y3
    1.0E+06, # z3 
    ]).T

P0 = eye(18) * APCOV_VALS
 
