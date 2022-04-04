# Fourier Spectrum of the CER resonant argument
# Author: Joseph A'Hearn
# Created 04/11/2017
#
# This program is plotting a Fourier spectrum of the resonant argument
# 

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
	#with open(filename) as f:
	#	line_count = np.sum(1 for _ in f)
	#	data_points = line_count - 4
	mimdata = np.loadtxt(filename, skiprows=4)
	t = mimdata[:,0]
	t_yr = t / 365.25
	lambda_mim = mimdata[:,4]
	lngper_mim = mimdata[:,5]
	return t_yr, lambda_mim, lngper_mim


black = True
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
		t = np.linspace(0, simulation_time * 365.25, data_points)
		# taking the fourier transform by brute force
		#expected_period = 3
		#period_margin_of_error = 2
		#expected_freq = 2 * np.pi / expected_period
		#freq_margin_of_error = 2 * np.pi / (expected_period - period_margin_of_error)
		n_samples = 2**10
		#sample_freq = np.linspace(0.003, 0.03, n_samples)
		#sample_period = np.linspace(10, 3500, n_samples)
		sample_period = np.linspace(2**9, 2**11, n_samples)
		y_1_omega = []
		y_2_omega = []
		amplitude_fourier = []
		period_by_fourier = []
		phi = phi_CER - np.mean(phi_CER)
		for i in range(n_samples):
			#y_1_omega.append(np.sum(phi_CER * np.cos(sample_freq[i] * t)))
			#y_2_omega.append(np.sum(phi_CER * np.sin(sample_freq[i] * t)))
			y_1_omega.append(np.sum(phi * np.cos(2 * np.pi * t / sample_period[i])))
			y_2_omega.append(np.sum(phi * np.sin(2 * np.pi * t / sample_period[i])))
			amplitude_fourier.append(2 * np.sqrt(y_1_omega[i]**2 + y_2_omega[i]**2) / data_points)
			#period_by_fourier.append(2 * np.pi / sample_freq[i])

		#transform_of_phi_CER = fft(phi_CER)
		#xf = np.linspace(2 * np.pi / 3500, 2 * np.pi / 10, data_points)
		#Tf = 2 * np.pi / xf

		print(np.amax(amplitude_fourier))
		print(sample_period[np.argmax(amplitude_fourier)])
		#print(period_by_fourier[np.argmax(amplitude_fourier)])
		if (black == True):
			fig = plt.figure(facecolor='black')
		else:
			fig = plt.figure()

		ax1 = fig.add_subplot(211)
		ax1.set_xlim([0, simulation_time])
		ax1.set_xlabel('Time (years)')
		ax1.set_ylabel(r'$\phi_{CER}$ (deg)')
		ax1.set_title("Resonant argument vs. time (7:6 CER with Mimas)")

		if(black==True):
			ax1.spines['top'].set_visible(False)
			ax1.spines['right'].set_visible(False)
			ax1.spines['bottom'].set_linewidth(0.5)
			ax1.spines['left'].set_linewidth(0.5)
			ax1.spines['bottom'].set_color('white')
			ax1.spines['left'].set_color('white')
			ax1.title.set_color('white')
			ax1.yaxis.label.set_color('white')
			ax1.xaxis.label.set_color('white')
			ax1.tick_params(axis='x', colors='white')
			ax1.tick_params(axis='y', colors='white')
			ax1.tick_params(axis='both', direction='in')
			ax1.get_xaxis().tick_bottom()
			ax1.get_yaxis().tick_left()
			ax1.scatter(t_yr,phi_CER, c='y', s=0.05)
		else:
			ax1.scatter(t_yr,phi_CER, s=0.05)

		ax2 = fig.add_subplot(212)
		ax2.set_xlabel('Period (days)')
		ax2.set_ylabel('Amplitude (deg)')
		#ax2.set_xlim([0, 3000])
		ax2.set_title("Fourier spectrum of resonant argument (7:6 CER with Mimas)")

		if(black==True):
			ax2.spines['top'].set_visible(False)
			ax2.spines['right'].set_visible(False)
			ax2.spines['bottom'].set_linewidth(0.5)
			ax2.spines['left'].set_linewidth(0.5)
			ax2.spines['bottom'].set_color('white')
			ax2.spines['left'].set_color('white')
			ax2.title.set_color('white')
			ax2.yaxis.label.set_color('white')
			ax2.xaxis.label.set_color('white')
			ax2.tick_params(axis='x', colors='white')
			ax2.tick_params(axis='y', colors='white')
			ax2.tick_params(axis='both', direction='in')
			ax2.get_xaxis().tick_bottom()
			ax2.get_yaxis().tick_left()
			ax2.plot(sample_period, amplitude_fourier, c='y', label='Amplitude = ' + str("{:.2f}".format(np.amax(amplitude_fourier))) + ' deg\nPeriod = ' + str("{:.2f}".format(sample_period[np.argmax(amplitude_fourier)])) + ' days')
		else:
			ax2.plot(sample_period, amplitude_fourier, label='Amplitude = ' + str(np.amax(amplitude_fourier)) + ' deg\nPeriod = ' + str(sample_period[np.argmax(amplitude_fourier)]) + ' days')
		
		plt.legend(loc=0, prop={'size':7})
		plt.tight_layout()
		if(black==True):
			ax1.figure.savefig("Fourier_black.png", facecolor=fig.get_facecolor(), transparent=True)
		else:
			ax1.figure.savefig("Fourier.png")
		plt.clf()		