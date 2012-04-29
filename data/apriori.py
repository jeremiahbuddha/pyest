
from scipy import matrix, eye

# ==============================================================================
# Project Constants

# Earth constants
dTheta = 7.29211585530066E-5
J2 = 1.082626925638815E-3
mu = 3.986004415E+14  # m**3 / s**2
R = 6378136.3 # m

# Drag Params
Cd = 2.0 # unitless
A = 3.0 # m**2
m = 970.0 # kg
rho_0 = 3.614E-13 # kg / m**3
ref_r = 7078136.3 # m
h_step = 88667.0 # m
drag_C = (1.0 / 2.0) * Cd * (A / m)

# Tracking station coordinates (ECF)
stn101 = [-5127510.0,-3794160.0, 0.0]
stn337 = [3860910.0, 3238490.0, 3898094.0]
stn394 = [549505.0,-1380872.0, 6182197.0]

# ==============================================================================
# PROJECT MEASUREMENT COVARIANCE
range_sigma = 0.01 # meters
range_rate_sigma = 0.001 # meters/sec

W = matrix( [ [ 1.0 / range_sigma**2 , 0.0                       ],
              [ 0.0                  , 1.0 / range_rate_sigma**2 ] ] )

# ==============================================================================
# PROJECT STATE AND STATE COVARIANCE

# Initial State (units in meters and secons)
INITIAL_X0 = matrix([
      757700.0,  # X
      5222607.0, # Y
      4851500.0, # Z
      2213.21,   # dX
      4678.34,   # dY 
     -5371.30,   # dZ
      mu,        # mu
      J2,        # J2
      Cd,        # Cd
      stn101[0], # X_1
      stn101[1], # Y_1
      stn101[2], # Z_1
      stn337[0], # X_1
      stn337[1], # Y_1
      stn337[2], # Z_1
      stn394[0], # X_1
      stn394[1], # Y_1
      stn394[2], # Z_1
      ]).T

INITIAL_x0 = matrix( [      
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

INITIAL_P0 = matrix(eye(18))
for i in range(18):
    INITIAL_P0[i,i] = APCOV[i]
