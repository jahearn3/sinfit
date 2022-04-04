# Sliding Box Method
# Author: Joseph A'Hearn
# Created 04/11/2017
#
# This program is plotting a Fourier spectrum of the resonant argument
#   in a dynamic timeframe in order to track the change in amplitude over time

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.fftpack import fft

def list_of_bodies(): # excludes Mimas
	body_list = []
	filename = 'big.in'
	with open(filename) as f:
		content = f.readlines()
	content = [x.strip() for x in content] 
	for i in range(10, len(content)):
		if((i - 10) % 4 == 0): # singles out the lines containing the names of the bodies (except Mimas)
			line_with_body = content[i]
			body = line_with_body[:8]
			body = body.strip(' ')
			body_list.append(body)
	return body_list

def import_mimas_data():
	# Getting the data from Mimas
	filename = 'MIMAS.geo'
	mimdata = np.loadtxt(filename, skiprows=4)
	t = mimdata[:,0]
	t_yr = t / 365.25
	lambda_mim = mimdata[:,4]
	lngper_mim = mimdata[:,5]
	return t_yr, lambda_mim, lngper_mim



body_list = list_of_bodies()
t_yr, lambda_mim, lngper_mim = import_mimas_data()
simulation_time = np.amax(t_yr)
with open("fourier_output.txt", "a") as myfile:
	for body in body_list:
		boddata = np.loadtxt(body + '.geo', skiprows=4)
		a_bod = boddata[:,1]
		lambda_bod = boddata[:,4]
		phi_CER = ((7 * lambda_mim) - (6 * lambda_bod) - lngper_mim) % 360
		data_points = len(phi_CER) 

		# split into a chunk of time that will slide until the end
		end_of_first_interval_idx = int(round(3 * data_points / 5))
		t_i = t_yr[:end_of_first_interval_idx]
		phi_i = phi_CER[:end_of_first_interval_idx]

		j = 0
		slide = 10
		while(end_of_first_interval_idx + (slide * j) < data_points):

			t = np.linspace(0, t_yr[end_of_first_interval_idx] * 365.25, int(round(3 * data_points / 5)))
			# taking the fourier transform by brute force
			n_samples = 2**10
			sample_period = np.linspace(2**9, 2**11, n_samples)
			y_1_omega = []
			y_2_omega = []
			amplitude_fourier = []
			period_by_fourier = []
			for i in range(n_samples):
				y_1_omega.append(np.sum(phi_i * np.cos(2 * np.pi * t / sample_period[i])))
				y_2_omega.append(np.sum(phi_i * np.sin(2 * np.pi * t / sample_period[i])))
				amplitude_fourier.append(2 * np.sqrt(y_1_omega[i]**2 + y_2_omega[i]**2) / data_points)
	
			print(np.amax(amplitude_fourier), sample_period[np.argmax(amplitude_fourier)])
			myfile.write('\t' + str(np.amax(amplitude_fourier)) + '\t' + str(sample_period[np.argmax(amplitude_fourier)]))
	
			fig = plt.figure()
			ax1 = fig.add_subplot(211)
			ax1.set_xlim([t_i[0], t_i[-1]])
			ax1.plot(t_i, phi_i)
			ax1.set_xlabel('Time (years)')
			ax1.set_ylabel(r'$\phi_{CER}$')
			ax2 = fig.add_subplot(212)
			ax2.plot(sample_period, amplitude_fourier)
			ax2.set_xlabel('Period (days)')
			ax2.set_ylabel('Amplitude (deg)')
			
			plt.tight_layout()
			ax1.figure.savefig("Fourier" + str(j) + ".png")
			plt.clf()
			# slide the box
			j += 1	
			t_i = t_yr[(slide * j) : (end_of_first_interval_idx + (slide * j))]
			phi_i = phi_CER[(slide * j) : (end_of_first_interval_idx + (slide * j))]
				
