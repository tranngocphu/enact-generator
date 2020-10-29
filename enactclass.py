''' 
Implementation of the paper:
EnAcT: Generating Aircraft Encounters using a Spherical Earth Model
http://www.atmseminar.org/seminarContent/seminar13/papers/ATM_Seminar_2019_paper_63.pdf
'''

from math import cos, sin, asin, acos, sqrt
import numpy as np
import pandas as pd
from scipy import optimize

from constants import *
from tools import to360, toPI


'''
ConflictInput Class
Construct objects for conflict inputs params
Input arguments from users are in aviation units,
the params are converted to SI units.
This serves as input interface for ConflictGenerator Class
'''
class ConflictInput: 

    def __init__(
        self,
        lat0,
        lon0,
        hdg0, 
        encounter_level, 
        encounter_angle, 
        h_distance, 
        v_distance, 
        ownship_phase, 
        intruder_phase,
        ownship_type, 
        intruder_type, 
        ownship_cs, 
        intruder_cs,
        look_ahead,
        duration_after,
        bada_data_path
    ):        
        self.look_ahead = look_ahead * 60         # look ahead time, convert mins to secs
        self.duration_after = duration_after * 60 # duration after the conflict for trajectory generation
        self.ownship_cs = ownship_cs              # CallSign of the ownship
        self.intruder_cs = intruder_cs            # CallSign of the intruder
        self.V0 = v_distance * ft2m               # Vertical separation at CPA, convert from nm to m
        self.H0 = h_distance * nm2m               # Horizontal separation at CPA, convert feet to m
        self.cpa_level = encounter_level          # Ownship FL at CPA
        self.phi00     = lat0 * deg2rad           # Latitude of ownship at CPA, convert from degree to radian
        self.lambda00  = lon0 * deg2rad           # Longitude of ownship at CPA, convert from degree to radian
        self.h00 = encounter_level * 100 * ft2m   # Ownship altitude at CPA
        self.theta00 = hdg0 * deg2rad             # Heading of the ownship at CPA
        self.dtheta0 = encounter_angle * deg2rad  # Encounter angle at CPA, convert from degree to radian
        self.ownship_phase  = ownship_phase       # Phase of the ownship        
        self.intruder_phase = intruder_phase      # Phase of the intruder
        self.ownship_type = ownship_type          # Ownship aircraft type
        self.intruder_type = intruder_type        # Intruder aircraft type
        self.bada_data_path = bada_data_path      # path/to/bada/data/file



class ConflictGenerator:
    
    def __init__( self, params ):
        self.look_ahead     = params.look_ahead
        self.duration_after = params.duration_after
        self.ownship_cs     = params.ownship_cs
        self.intruder_cs    = params.intruder_cs
        self.V0             = params.V0
        self.H0             = params.H0
        self.cpa_level      = params.cpa_level  
        self.phi00          = params.phi00    
        self.lambda00       = params.lambda00 
        self.h00            = params.h00
        self.h10            = None
        self.theta00        = params.theta00
        self.dtheta0        = params.dtheta0
        self.theta10        = None   
        self.delta0         = None        
        self.ownship_type   = params.ownship_type        
        self.ownship_phase  = params.ownship_phase
        self.intruder_type  = params.intruder_type
        self.intruder_phase = params.intruder_phase
        self.beta0 = None
        self.tas0  = None
        self.tas1  = None
        self.vs0   = None
        self.vs1   = None
        self.dvs   = None        
        self.N00   = None
        self.D00   = None
        self.N10   = None
        self.D10   = None      
        self.track = None
        self.bada_data_path = params.bada_data_path
        self.pre_process()
        self.compute_cpa()
        
        
    
    def pre_process( self ):                     
        self.get_speed(self.cpa_level, self.cpa_level)  # get aircraft speed from BADA
        self.delta0 = self.H0 / R                       # central angle of 2 aircraft at CPA
        # self.theta00 = toPI(self.theta00)             # bring theta00 to (-pi, pi] if neccessary
        self.theta10 = self.theta00 + self.dtheta0      # find heading of intruder  
        # self.theta10 = toPI(self.theta10)             # bring theta10 to (-pi, pi] if neccessary
        self.h10 = self.h00 + self.V0    
            
    
    def get_speed( self, fl0, fl1 ):        
        ptf = pd.read_csv(self.bada_data_path)  # read performance table from BADA
        
        # OWNSHIP        
        ownship     = ptf[(ptf.ac_type == self.ownship_type) & (ptf.fl == str(fl0))]  # select row according to ac_type and FL
        self.tas0   = ownship.iloc[0][self.ownship_phase + "_TAS"] * knot2ms  # read true air speed
        self.omega0 = self.tas0 / (2*R) # convert to angular speed
        
        # get vertical speed of ownship:
        if self.ownship_phase != "CR" :
            self.vs0 = ownship.iloc[0][self.ownship_phase + "_ROCD_nom_m"] * (-1)**(int(self.phases[0]=="DES")) * ftmin2ms           
        else :
            self.vs0 = 0  # ownship is leveling        
        
        # INTRUDER        
        intruder    = ptf[(ptf.ac_type == self.intruder_type) & (ptf.fl == str(fl1))]
        self.tas1   = intruder.iloc[0][self.intruder_phase + "_TAS"] * knot2ms        
        self.omega1 = self.tas1 / (2*R)      
        
        # get vertical speed of intruder:
        if self.intruder_phase != "CR" :
            self.vs1 = intruder.iloc[0][self.intruder_phase + "_ROCD_nom_m"] * (-1)**(int(self.phases[1]=="DES")) * ftmin2ms         
        else :
            self.vs1 = 0  # intruder is leveling
        self.dvs = self.vs1 - self.vs0


    def f_beta0( self, beta0 ):
        K = self.V0 * self.dvs / self.H0
        g = cos(self.phi00) * sin(self.delta0) * cos(beta0) + sin(self.phi00) * cos(self.delta0)
        dg = - cos(self.phi00) * sin(self.delta0) * sin(beta0)        
        f = cos(self.phi00) * cos(self.delta0) * cos(self.theta10) *  cos(beta0) - \
            sin(self.phi00) * sin(self.delta0) * cos(self.theta10) + \
            cos(self.phi00) * sin(self.theta10) * sin(beta0)
        df = cos(self.phi00) * sin(self.theta10) * cos(beta0) - cos(self.phi00) * cos(self.theta10) * cos(self.delta0) * sin(beta0)        
        fbeta0 = self.tas1 * ( f / sqrt(1 - g**2) ) - self.tas0 * cos(beta0 - self.theta00) + K
        dfbeta0 = self.tas0 * sin(beta0 - self.theta00) + self.tas1 * df / sqrt(1 - g**2) + self.tas1 * f * g * dg / (1 - g**2)**(3/2)
        return fbeta0, dfbeta0


    def find_beta0( self ):
        guess = 0
        done = False
        while not done :            
            solution = optimize.root_scalar(self.f_beta0, x0=guess, fprime=True, method='newton', xtol=1e-6, maxiter=20)
            done = solution.converged 
            guess += PI / 6
            if guess > 2 * PI :
                solution = False
                break
        return solution


    def compute_cpa( self ) :
        beta0 = self.find_beta0() 

        if beta0 is False :
            self.beta0 = False
            return 
        else :    
            self.beta0 = toPI(beta0.root)
            self.beta0 = beta0.root

        # compute N00, north00, east00, D00
        self.N00 = np.array([ \
            [ cos(self.phi00) * cos(self.lambda00) ], \
            [ cos(self.phi00) * sin(self.lambda00) ], \
            [ sin(self.phi00) ]
        ])        
        
        self.north00 = np.array([
            [-sin(self.phi00) * cos(self.lambda00)],
            [-sin(self.phi00) * sin(self.lambda00)],
            [ cos(self.phi00) ]
        ])

        self.east00 = np.array([
            [-sin(self.lambda00)],
            [ cos(self.lambda00)],
            [ 0 ]
        ])

        self.D00 = cos(self.theta00) * self.north00 + sin(self.theta00) * self.east00

        # compute B0
        self.B0 = cos(self.beta0) * self.north00 + sin(self.beta0) * self.east00

        #compute N10
        self.N10 = cos(self.delta0) * self.N00 + sin(self.delta0) * self.B0

        # compute Sine and Cosine of phi10, lambda10
        self.sin_phi10    = self.N10[2][0]  # Sin(Phi10) = Z-component of N10, equation 7
        self.cos_phi10    = sqrt(1-self.sin_phi10**2)
        self.sin_lambda10 = self.N10[1][0] / self.cos_phi10 # equation 7, Y-component of N10
        self.cos_lambda10 = self.N10[0][0] / self.cos_phi10 # equation 7, X-component of N10

        # compute north10, east10, and D10
        self.north10 = np.array([
            [-self.sin_phi10 * self.cos_lambda10],
            [-self.sin_phi10 * self.sin_lambda10],
            [ self.cos_phi10 ]
        ])

        self.east10 = np.array([
            [-self.sin_lambda10],
            [ self.cos_lambda10],
            [ 0 ]
        ])

        self.D10 = cos(self.theta10) * self.north10 + sin(self.theta10) * self.east10  
    
    
    def ownship_track( self, t ):
        """Compute the location of the ownship at time t.
        Time t=0 gives the location at CPA. Negative time 
        indicates location before CPA.

        Arguments:
            t (integer) : time
        
        Returns:
           [lat, lon, alt] in [degree, degree, meter]
        
        """
        self.N0 = cos(self.omega0*t) * self.N00 + sin(self.omega0*t) * self.D00
        sin_phi = self.N0[2][0]
        cos_phi = sqrt(1-sin_phi**2)
        cos_lambda = self.N0[0][0] / cos_phi
        phi0 = asin(sin_phi)
        lambda0 = acos(cos_lambda)         
        h0 = self.h00 + self.vs0 * t
        return phi0*rad2deg, lambda0*rad2deg, h0
        
    
    def intruder_track( self, t ):
        """Compute the location of the intruder at time t.
        Time t=0 gives the location at CPA. Negative time 
        indicates location before CPA.

        Arguments:
            t (integer) : time
        
        Returns:
           [lat, lon, alt] in [degree, degree, meter]
        
        """
        self.N1 = cos(self.omega1*t) * self.N10 + sin(self.omega1*t) * self.D10
        sin_phi = self.N1[2][0]
        cos_phi = sqrt(1-sin_phi**2)
        cos_lambda = self.N1[0][0] / cos_phi
        phi1 = asin(sin_phi)
        lambda1 = acos(cos_lambda)     
        h1 = self.h10 + self.vs1 * t    
        return phi1*rad2deg, lambda1*rad2deg, h1


    def compute_track( self, save_as='' ):

        if self.beta0 is False :
            self.track = "No beta0 found."
            return False

        if save_as == 'df':
            df = pd.DataFrame(columns=[
                "CallSign",
                "Destination",
                "ExerciseTime",
                "FightLevel",
                "Heading",
                "Latitude",
                "Longitude",
                "NavigationStatus",
                "SSRCode",
                "Speed",
                "TrackNumber",
                "Type",
                "X",
                "Y"
            ])
        else:
            time = []
            ownship_latlon = []        
            ownship_alt = []
            ownship_fl = []
            intruder_latlon = []
            intruder_alt = []
            intruder_fl = []


        for t in range( -self.look_ahead, self.duration_after ) :
            lat0, lon0, h0 = self.ownship_track(t)
            lat1, lon1, h1 = self.intruder_track(t)    
            
            if save_as == 'df':
                df = df.append({
                    "CallSign": self.ownship_cs,
                    "Destination": "",
                    "ExerciseTime": t + self.look_ahead,
                    "FightLevel": h0,
                    "Heading": self.theta00,
                    "Latitude": lat0,
                    "Longitude": lon0,
                    "NavigationStatus": "",
                    "SSRCode": "",
                    "Speed": self.tas0,
                    "TrackNumber": "",
                    "Type": self.ownship_type,
                    "X": "",
                    "Y": ""                
                }, ignore_index=True)
            
                df = df.append({
                    "CallSign": self.intruder_cs,
                    "Destination": "",
                    "ExerciseTime": t + self.look_ahead,
                    "FightLevel": h1,
                    "Heading": self.theta10,
                    "Latitude": lat1,
                    "Longitude": lon1,
                    "NavigationStatus": "",
                    "SSRCode": "",
                    "Speed": self.tas1,
                    "TrackNumber": "",
                    "Type": self.intruder_type,
                    "X": "",
                    "Y": ""                
                }, ignore_index=True)

            else :
                time.append(t + self.look_ahead)            
            
                ownship_latlon.append([round(lat0,3), round(lon0,3)])
                ownship_alt.append(round(h0))
                ownship_fl.append(round(h0*m2ft/100))
            
                intruder_latlon.append([round(lat1,3), round(lon1,3)])
                intruder_alt.append(round(h1))
                intruder_fl.append(round(h1*m2ft/100))

        
        if save_as is 'df':
            self.track = df
        else:        
            self.track = {                
                self.ownship_cs: {
                    'time': time,
                    'latlon': ownship_latlon,
                    'alt': ownship_alt,
                    'fl': ownship_fl
                },
                self.intruder_cs: {
                    'time': time,
                    'latlon': intruder_latlon,
                    'alt': intruder_alt,
                    'fl': intruder_fl
                },
            }
        
        return True
