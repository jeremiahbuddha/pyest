#! /usr/bin/env python

from math import pi, sqrt 

# ============================================================================
# ============================================================================
# Constants and functions

deg2rad = pi/180.0
rad2deg = 1/deg2rad

# Planetary dictionary
planets = {     #         GM Eq.      Eq.Rad  Mean Orb Rad.     J2
                #   (km**3/sec**2)     (km)       (km)      (unitless)
    'mercury' : [        22032.080,   2439.7,   57909225.7,  0          ],
    'venus'   : [       324858.599,   6051.8,  108209019.0,  0          ],
    'earth'   : [       398600.400,   6378.145,  149598016.0,  0.00108248 ],
#    'earth'   : [       398600.433,   6378.1,  149598016.0,  0          ],
    'mars'    : [        42828.314,   3397.0,  227936834.0,  0          ],
    'jupiter' : [    126712767.858,  71492.0,  778412700.0,  0          ],
    'saturn'  : [     37940626.061,  60268.0, 1426725400.0,  0          ],
    'uranus'  : [      5794549.007,  25559.0, 2870972200.0,  0          ],
    'neptune' : [      6836534.064,  24764.0, 4498252900.0,  0          ],
    'sun'     : [ 132712440017.987, 696000.0,     'N/A'   ,  0          ],
    'titan'   : [         8978.100,   2575.5,     'N/A'   ,  0          ],
    }

def mag( vec ):
    return sqrt( float( vec[0]**2 ) + 
                   float( vec[1]**2 ) + 
                   float( vec[2]**2 ) )

# Import module classes for easier access
from state import State
from traj import Traj
