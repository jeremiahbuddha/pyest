#! /usr/bin/env python

from math import pi, sqrt, sin, cos, tan, acos, asin, atan, atan2
from numpy import array, cross, dot
from scipy import matrix, array
from ASEN5070.integ import planets, deg2rad, rad2deg, mag

def angle_0_to_2pi(ang):
    if ang < 0.0:
        ang = 2*pi + ang
    return ang

# =============================================================================
# =============================================================================
# State Class
class State(object):

    def __init__(self, time = 0, cartesian = [], classical = [], center='earth'):
        self._center = center
        self._mu = planets[center][0]
        self._time = time
        if cartesian and not classical:
            self._cartesian = cartesian
            self._classical = self._cartesian_to_classical(cartesian)
        elif classical and not cartesian:
            self._classical = classical
            self._cartesian = self._classical_to_cartesian(classical)
        elif not classical and not cartesian:
            raise UserError(
                "Must supply either a cartesian or classical state.")
        elif classical and cartesian:
            raise UserError(
                "Must supply either a cartesian or classical state, not both.")

    def _cartesian_to_classical(self,cartesian):
        # Split input cartesian list into position and velocity arrays
        r = array([cartesian[0], cartesian[1], cartesian[2]])
        v = array([cartesian[3], cartesian[4], cartesian[5]])
        r_mag = mag(r)
        v_mag = mag(v)
        # Compute sepcific angular momentum
        h = cross(r,v)
        h_mag = mag(h)        
        # Compute specific energy
        specific_energy = v_mag**2/2.0 - self._mu/r_mag 
        # Compute SMA
        sma = -self._mu/(2*specific_energy)
        # Compute eccentricity
        ecc = sqrt(1 - h_mag**2/(sma*self._mu))
        # Compute inclination
        inc = acos(h[2]/h_mag)
        # Compute Right Ascension of Ascending Node
        raan = atan2(h[0],-h[1])
        raan = angle_0_to_2pi(raan) 
        # Compute the argument of latitude
        arg_lat = atan2(r[2]/sin(inc),(r[0]*cos(raan)+r[1]*sin(raan)))
        arg_lat = angle_0_to_2pi(arg_lat)
        # Compute the true anomaly
        TA = acos( (sma*(1-ecc**2) - r_mag)/(ecc*r_mag) ) 
        if dot(r,v) < 0:
            TA = 2*pi - TA
        TA = angle_0_to_2pi(TA)
        # Compute the Argument of periapse
        arg_p = arg_lat - TA
        arg_p = angle_0_to_2pi(arg_p)
        # Compute the Eccentric Anomaly
        EA = 2*atan( sqrt( (1-ecc)/(1+ecc) )*tan(TA/2) )  
        EA = angle_0_to_2pi(EA)
        # Compute the Mean Anomaly
        MA = EA - ecc*sin(EA)
        MA = angle_0_to_2pi(MA)
        # Compute the orbital period
        T = ( 2*pi/self._mu**2 )*( h_mag/sqrt(1-ecc**2) )**3 
        # Compute time since periapsis
        time_since_peri = MA*T/( 2*pi ) 

        return [sma, ecc, inc, raan, arg_p, MA] 

    def _classical_to_cartesian(self,classical):
        # Split the input classical list into elements
        sma, ecc, inc, raan, arg_p, MA = classical
        # Solve keplers equation for Eccentric Anomaly
        EA_iter = MA
        for i in xrange(0,10000):
            EA = MA + ecc*sin(EA_iter)
            if ( EA - EA_iter < 0.0001):
                break;
            else:
                EA_iter = EA
        # Compute True anomaly
        TA = 2*atan( sqrt( (1+ecc)/(1-ecc) )*tan( EA/2 ) )
        if EA > pi:
            TA = 2*pi - TA
        # Compute the argument of latitude
        arg_lat = TA + arg_p
        # Compute position magnitude
        r_mag = sma*( 1-ecc*cos(EA) )
        # Compute specific angular momentum magnitude
        h_mag = sqrt( self._mu*sma*( 1-ecc**2 ) )
        # Compute cartesian position
        X = r_mag*( cos(raan)*cos(arg_lat) - sin(raan)*sin(arg_lat)*cos(inc) )
        Y = r_mag*( sin(raan)*cos(arg_lat) + cos(raan)*sin(arg_lat)*cos(inc) )
        Z = r_mag*( sin(inc)*sin(arg_lat) )
        # Compute cartesian velocity
        preamble = h_mag*ecc*sin(TA)/( r_mag*sma*(1-ecc**2) ) 
        dX =  X*preamble - ( h_mag/r_mag )*( cos(raan)*sin(arg_lat) + sin(raan)*cos(arg_lat)*cos(inc) )
        dY =  Y*preamble - ( h_mag/r_mag )*( sin(raan)*sin(arg_lat) - cos(raan)*cos(arg_lat)*cos(inc) )
        dZ =  Z*preamble + ( h_mag/r_mag )*( sin(inc)*cos(arg_lat ) )
        
        return [ X, Y, Z, dX, dY, dZ ]

    def get_time(self):
        return self._time

    def get_mu(self):
        return self._mu

    def get_cartesian(self):
        return self._cartesian

    def get_classical(self):
        return self._classical

    def get_period(self):
        return ( 2.0*pi/sqrt(self._mu) )*( self._classical[0]**(1.5) )

    def get_center(self):
        return self._center

