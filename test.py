from enactclass import *

# conflict inputs are in aviation units

lat0 = 2.000  # lat of ownship at CPA
lon0 = 101.000  # lon of ownship at CPA
counter_level = 350  # hundred of ft
v_distance = 0 # ft
h_distance = 3 # nm
hdg0 = 90  # degree
counter_angle = 30  # degree
ownship_phase = "CR"
intruder_phase = "CR"
ownship_type = "A320"
intruder_type = "A320"
look_ahead = 5  # minutes
duration_after = 10 # miutes after the conflict
ownship_cs  = "AUT001"
intruder_cs = "AUT002"
bada_data_path = '/home/phu/projects/enact-generator/bada36_PTF.csv'

conflict_params = ConflictInput(
    lat0, lon0, hdg0, \
    counter_level, counter_angle, h_distance, v_distance, \
    ownship_phase, intruder_phase, ownship_type, intruder_type, \
    ownship_cs, intruder_cs, look_ahead, duration_after, bada_data_path
)

conflict = ConflictGenerator(conflict_params)

conflict.compute_track(save_as='df')  # conflict.track will be a dataframe

print(conflict.beta0)

print(conflict.track)