
from scipy import matrix, eye

range_sigma = 0.01 # meters
range_rate_sigma = 0.001 # meters/sec

W = matrix( [ [ 1.0 / range_sigma**2 , 0.0                       ],
              [ 0.0                  , 1.0 / range_rate_sigma**2 ] ] )

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

scl = 1.0

APCOV = matrix( [ 
    scl*1.0E+06, # x
    scl*1.0E+06, # y
    scl*1.0E+06, # z 
    scl*1.0E+06, # dx
    scl*1.0E+06, # dy
    scl*1.0E+06, # dz
    scl*1.0E+20, # mu
    scl*1.0E+06, # J2
    scl*1.0E+06, # Cd
    scl*1.0E-10, # x1
    scl*1.0E-10, # y1
    scl*1.0E-10, # z1
    scl*1.0E+06, # x2
    scl*1.0E+06, # y2
    scl*1.0E+06, # z2 
    scl*1.0E+06, # x3
    scl*1.0E+06, # y3
    scl*1.0E+06, # z3 
    ]).T

P0 = matrix(eye(18))
for i in range(18):
    P0[i,i] = APCOV[i]
