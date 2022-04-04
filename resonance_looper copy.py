# Resonance Looper
#
# Author: Joseph A'Hearn
# Created 03/20/2017
#
# This program loops through mean longitude values and semimajor axis values 
#   near the corotation eccentricity resonance and runs a simulation at each
#   set of values in the resonance space in order to create a resonance map.

from subprocess import call
import sys
import os
# WARNING: the following 4 lines are read from the resonance map scripts.
mean_longitude_min = 180.5
mean_longitude_max = 239.0
semimajor_axis_min = 167506.5
semimajor_axis_max = 167570.5

mean_longitude_step_size = 0.5
semimajor_axis_step_size = 0.5

mean_longitude = mean_longitude_min
semimajor_axis = semimajor_axis_min

#mean_longitude = 210
#semimajor_axis = 167537.0

call('python reset_big.py', shell=True)

while(mean_longitude < mean_longitude_max):
	call('python geo2xyz.py ' + str(mean_longitude) + ' ' + str(semimajor_axis), shell=True) # this modifies big.in
	while(semimajor_axis < semimajor_axis_max):
		call('python geo2xyz.py ' + str(mean_longitude) + ' ' + str(semimajor_axis), shell=True) # this modifies big.in
		call('sh ./simulate_and_plot_with_fit.sh', shell=True)
		call('mv Res_Args_w_sin_fit.png ' + str(mean_longitude) + '_' + str(semimajor_axis) + '.png', shell=True)
		call('rm *mp *out *aei *geo', shell=True)
		call('python reset_big.py', shell=True)
		semimajor_axis += semimajor_axis_step_size
	mean_longitude += mean_longitude_step_size
	semimajor_axis = semimajor_axis_min

call('python resonance_map.py', shell=True)
call('python resonance_map_LER.py', shell=True)