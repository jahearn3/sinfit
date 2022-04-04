# Sine Fit for both the CER and the LER
# Author: Joseph A'Hearn
# Created 03/28/2017
#
# This program is for fitting data points to a sine function of varying amplitude
# 

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import leastsq 
import random


# Parameters that often change from one simulation to another
body_list = ["AEGAEON"]
period_guess = 4.0
black = True
fits = True

# Functions

def interval_minimum(t, phi, a, b):
	phi_min = np.amin(phi[a:b])
	i_min = np.argmin(phi[a:b])
	t_cut = t[a:b]
	t_min = t_cut[i_min]
	return t_min, phi_min

def interval_maximum(t, phi, a, b):
	phi_max = np.amax(phi[a:b])
	i_max = np.argmax(phi[a:b])
	t_cut = t[a:b]
	t_max = t_cut[i_max]
	return t_max, phi_max

def sine_fit(simulation_time, A_0, growth, shift, freq, phase):
	t_fit = np.linspace(0, simulation_time, simulation_time * 10000)
	phi_fit = A_0 * (1 + (growth * (t_fit / simulation_time))) * np.sin((freq * t_fit) + phase) + shift
	return t_fit, phi_fit

def sine_lstsq(t, phi, A_guess, growth_guess, shift_guess, freq_guess, phase_guess, data_points): 
	t_lstsq = np.linspace(0, np.amax(t), data_points)
	optimize_func = lambda x: x[0] * (1 + (x[1] * (t_lstsq / np.amax(t)))) * np.sin((x[2] * t_lstsq) + x[3]) + x[4] - phi
	A_est, growth_est, freq_est, phase_est, shift_est = leastsq(optimize_func, [A_guess, growth_guess, freq_guess, phase_guess, shift_guess])[0]
	phi_lstsq = A_est * (1 + (growth_est * (t_lstsq / np.amax(t)))) * np.sin((freq_est * t_lstsq) + phase_est) + shift_est 

	# average A(t) over the first 4 years (period cycle) should be at the 2 years mark
	# A(t) = A_est * (1 + (growth_est * (t_lstsq / np.amax(t))))
	#period = 2 * np.pi / freq_est

	A_t = A_est * (1 + (0.5 * growth_est ))


	return t_lstsq, phi_lstsq, A_t, growth_est

#def sine_lstsq_constant_A(t, phi, A_guess, shift_guess, freq_guess, phase_guess, data_points):
#	t_lstsq = np.linspace(0, np.amax(t), data_points)
#	optimize_func = lambda x: x[0] * np.sin((x[1] * t_lstsq) + x[2]) + x[3] - phi
#	A_est, freq_est, phase_est, shift_est = leastsq(optimize_func, [A_guess, freq_guess, phase_guess, shift_guess])[0]
#	phi_lstsq = A_est * np.sin((freq_est * t_lstsq) + phase_est) + shift_est 
#	return t_lstsq, phi_lstsq, A_est
#
#def sine_lstsq_linear_A(t, phi, A_guess, growth_guess, shift_guess, freq_guess, phase_guess, data_points):
#	t_lstsq = np.linspace(0, np.amax(t), data_points)
#	optimize_func = lambda x: x[0] * (1 + (x[1] * (t_lstsq / np.amax(t)))) * np.sin((x[2] * t_lstsq) + x[3]) + x[4] - phi
#	A_est, growth_est, freq_est, phase_est, shift_est = leastsq(optimize_func, [A_guess, growth_guess, freq_guess, phase_guess, shift_guess])[0]
#	phi_lstsq = A_est * (1 + (growth_est * (t_lstsq / np.amax(t)))) * np.sin((freq_est * t_lstsq) + phase_est) + shift_est 
#	return t_lstsq, phi_lstsq, A_est, growth_est
#
#def sine_lstsq_sinusoidal_A():
#	t_lstsq = np.linspace(0, np.amax(t), data_points)
#	#optimize_func = lambda x: x[0] * (1 + (x[1] * (t_lstsq / np.amax(t)))) * np.sin((x[2] * t_lstsq) + x[3]) + x[4] - phi
#
#
#def sine_lstsq_exponential_A():
#	t_lstsq = np.linspace(0, np.amax(t), data_points)


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

def initial_estimates(t_yr, period_guess):
	data_points = len(t_yr)
	simulation_time = np.amax(t_yr)
	time_step = simulation_time / data_points # in terms of years
	i_period = int(np.rint(period_guess / time_step)) # estimate of the index at which one libration cycle is completed
	n_periods = np.floor(simulation_time / period_guess)
	return data_points, simulation_time, time_step, i_period, n_periods

#def initialize_plot(simulation_time, resonance, phi_CER, phi_LER):
#	# Getting the plot started
#	fig = plt.figure()
#	ax1 = plt.subplot2grid((6,6),(0, 0), rowspan=6, colspan=5)
#	ax1.set_title("7:6 CER")
#	ax1.set_xlabel('t (years)')
#	ax1.set_ylabel(r'$\phi_{CER}$ (deg)')
#	ax1.set_xlim([0, simulation_time])
#	ax1.scatter(t_yr, phi_CER, s=0.25, label=str(body))
#	ax2 = plt.subplot2grid((6,6),(3, 0), rowspan=3, colspan=5)
#	ax2.set_title("7:6 LER")
#	ax2.set_xlabel('t (years)')
#	ax2.set_ylabel(r'$\phi_{LER}$ (deg)')
#	ax2.set_xlim([0, simulation_time])
#	ax2.scatter(t_yr, phi_LER, s=0.25, label=str(body))

def CER_fit(phi, n_periods, t_yr, i_period, data_points, fits):
	try_no = 0
	if(np.amax(phi) - np.amin(phi) > 354):
		print(" Outside of resonance. ")
		A_est = 360 
		growth_est = 0
		rms_deviation = 0
	else:
		t_minima = []
		phi_minima = []
		t_maxima = []
		phi_maxima = []
		k = 0
		while(k < n_periods):
			minima_t, minima_phi = interval_minimum(t_yr, phi, k * i_period, i_period * (1 + k))
			t_minima.append(minima_t)
			phi_minima.append(minima_phi)
			maxima_t, maxima_phi = interval_maximum(t_yr, phi, k * i_period, i_period * (1 + k))
			t_maxima.append(maxima_t)
			phi_maxima.append(maxima_phi)
			k += 1
		# Seeing if there are any extrema in the incomplete period at the end of the simulation
		minima_t, minima_phi = interval_minimum(t_yr, phi, k * i_period, len(t_yr) - 1)
		t_minima.append(minima_t)
		phi_minima.append(minima_phi)
		maxima_t, maxima_phi = interval_maximum(t_yr, phi, k * i_period, len(t_yr) - 1)
		t_maxima.append(maxima_t)
		phi_maxima.append(maxima_phi)
		# initial amplitude: use the first max and the first min
		A_0 = (phi_maxima[0] - phi_minima[0]) / 2
		# final amplitude: use the last max and the last min
		A_f = (phi_maxima[len(phi_maxima) - 1] - phi_minima[len(phi_minima) - 1]) / 2
		#A_avg = (A_0 + A_f) / 2
		growth = (A_f - A_0) / A_0
		# shift: guess at the midpoint position of the sine function
		shift = phi_minima[0] + A_0
		# period and frequency: x difference between extrema (should be about 4 years); assume it stays constant
		maxima_period = t_maxima[1] - t_maxima[0]
		minima_period = t_minima[1] - t_minima[0]
		period = (maxima_period + minima_period) / 2
		freq = 2 * np.pi / period
		# phase: using phi/amplitude
		if (np.absolute(np.arcsin((shift - phi[0]) / A_0)) <= 1):
			phase = np.arcsin((shift - phi[0]) / A_0) + np.pi # not sure why adding pi makes this work better
		else:
			# we need another way to compute the phase if an invalid value is encountered in arcsin (e.g. slightly greater than 1)
			phase = np.pi
		t_fit, phi_fit = sine_fit(simulation_time, A_0, growth, shift, freq, phase)
		#ax.plot(t_fit, phi_fit, c='g', label=(str(body) + ' initial guess')) 
		
		t_lstsq, phi_lstsq, A_est, growth_est = sine_lstsq(t_yr, phi, A_0, growth, shift, freq, phase, data_points)

		# 10 attempts to fix when the magnitude of the growth constant > 1
		attempts = 0
		while(np.absolute(growth_est) > 1.0 and attempts < 10):
			delta0 = random.uniform(3.0, 10.0)
			delta1 = random.uniform(1.0, 2.0)
			t_lstsq, phi_lstsq, A_est, growth_est = sine_lstsq(t_yr, phi, (A_0 + A_est) / 2, growth_est / delta0, shift, freq * delta1, phase, data_points)
			attempts += 1
			try_no += 1

		if(np.absolute(growth_est) > 1.0):
			A_est = np.nan

		variance = 0
		for j in range(data_points):
			variance += (phi[j] - phi_lstsq[j])**2
		rms_deviation = np.sqrt(variance / (j + 1))
		# 10 attempts to fix when it's way off
		attempts = 0
		while(rms_deviation > 15 and attempts < 10):
			delta1 = random.uniform(0.7, 1.3)
			delta2 = random.uniform(0.7, 1.3)
			delta3 = random.uniform(0.7, 1.3)
			t_lstsq, phi_lstsq, A_est, growth_est = sine_lstsq(t_yr, phi, delta1 * A_0, growth, shift, delta2 * freq, delta3 * phase, data_points)
			variance = 0
			for i in range(data_points):
				variance += (phi[i] - phi_lstsq[i])**2
			rms_deviation = np.sqrt(variance / (i + 1))
			attempts += 1
			try_no += 100

		if(rms_deviation > 20):
			A_est = (A_0 + A_f) / 2 
		
		if(np.absolute(growth_est) > 1.0):
			A_est = np.nan
		if(fits == True):
			if(rms_deviation > 20 or np.absolute(growth_est) > 1.0):
				ax1.plot(t_lstsq, phi_lstsq, c='y') # yellow means the recorded amplitude does not correspond to the fit; the fit is probably off
			else:
				ax1.plot(t_lstsq, phi_lstsq, c='r', label='Amplitude = ' + str("{:.2f}".format(A_est)) + ' deg') 
		#ax1.plot(t_lstsq, phi_lstsq, c='r', label=(str(body) + ' least squares fit')) # should include parameters in label or put them on the plot in some other way
		
	return A_est, growth_est, rms_deviation, try_no

def LER_analysis(phi_LER, data_points, simulation_time, fits):
	variance = 0
	varimean = 0
	drops = 0
	rises = 0
	phi_mean = np.mean(phi_LER)
	for i in range(data_points):
		variance += (phi_LER[i])**2
		varimean += (phi_LER[i] - phi_mean)**2
		if(i > 0): # to count the number of times it goes over the +- 180 edge
			if((phi_LER[i])**2 > 100 and (phi_LER[i-1])**2 > 100):
				if(phi_LER[i] > 0 and phi_LER[i - 1] < 0): # from - to +
					drops += 1
				if(phi_LER[i] < 0 and phi_LER[i - 1] > 0): # from + to -
					rises += 1
	deviation_from_zero = np.sqrt(variance / data_points)
	deviation_from_mean = np.sqrt(varimean / data_points)
	if(fits == True):
		ax2.plot((0.0, simulation_time), (deviation_from_zero, deviation_from_zero), 'b-', label='deviation from zero = ' + str(deviation_from_zero))
		ax2.plot((0.0, simulation_time), (-deviation_from_zero, -deviation_from_zero), 'b-')
		ax2.plot((0.0, simulation_time), (phi_mean, phi_mean), 'r-', label=r'$\phi_{mean}$ = ' + str(phi_mean))
		ax2.plot((0.0, simulation_time), (phi_mean + deviation_from_mean, phi_mean + deviation_from_mean), 'g-')
		ax2.plot((0.0, simulation_time), (phi_mean - deviation_from_mean, phi_mean - deviation_from_mean), 'g-', label='deviation from ' + r'$\phi_{mean}$ = ' + str(deviation_from_mean))

	return deviation_from_zero, phi_mean, deviation_from_mean, drops, rises

# Main ------------------------------------------------------------------------------------

t_yr, lambda_mim, lngper_mim = import_mimas_data()
data_points, simulation_time, time_step, i_period, n_periods = initial_estimates(t_yr, period_guess)

with open("amplitude_in_res_space.txt", "a") as myfile:
	for body in body_list:
		boddata = np.loadtxt(body + '.geo', skiprows=4)
		a_bod = boddata[:,1]
		lambda_bod = boddata[:,4]
		lngper_bod = boddata[:,5]
		phi_CER = ((7 * lambda_mim) - (6 * lambda_bod) - lngper_mim) % 360
		phi_LER = ((((7 * lambda_mim) - (6 * lambda_bod) - lngper_bod) + 180 ) % 360 ) - 180

		# Getting the plot started
		if(black == True):
			fig = plt.figure(facecolor='black')
		else:
			fig = plt.figure()
		ax1 = fig.add_subplot(211)
		ax1.set_title("7:6 CER")
		#ax1.set_xlabel('t (years)')
		ax1.set_ylabel(r'$\phi_{CER}$ (deg)')
		ax1.set_xlim([0, simulation_time])
		if(black == True):
			ax1.scatter(t_yr, phi_CER, c='chartreuse', s=0.25)
		else: 
			ax1.scatter(t_yr, phi_CER, s=0.25)
		ax2 = fig.add_subplot(212)
		ax2.set_title("7:6 LER")
		ax2.set_xlabel('t (years)')
		ax2.set_ylabel(r'$\phi_{LER}$ (deg)')
		ax2.set_xlim([0, simulation_time])
		ax2.set_ylim([-180, 180])
		if(black == True):
			ax2.scatter(t_yr, phi_LER, c='chartreuse', s=0.25)
		else: 
			ax2.scatter(t_yr, phi_LER, s=0.25)

		A_CER, growth_CER, rms_dev_CER, try_no = CER_fit(phi_CER, n_periods, t_yr, i_period, data_points, fits)
		dev_zero_LER, phi_mean_LER, dev_mean_LER, drops, rises = LER_analysis(phi_LER, data_points, simulation_time, fits) 
		
		# output amplitude, initial semimajor axis, and initial mean longitude into a txt file
		# later I can make a resonance map with these data 
		myfile.write(str(a_bod[0]) + '\t' + str(lambda_bod[0]) + '\t' + str(np.absolute(A_CER)) + '\t' + str(growth_CER) + '\t' + str(rms_dev_CER) + '\t' + str(try_no) + '\t' + str(dev_zero_LER) + '\t' + str(phi_mean_LER) + '\t' + str(dev_mean_LER) + '\t' + str(drops) + '\t' + str(rises) + '\n')
		print(str(a_bod[0]) + '\t' + str(lambda_bod[0]) + '\t' + str(np.absolute(A_CER)) + '\t' + str(growth_CER) + '\t' + str(rms_dev_CER) + '\t' + str(dev_zero_LER) + '\t' + str(phi_mean_LER) + '\t' + str(dev_mean_LER) + '\n')

# Finishing the plot
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':6})
ax1.legend(loc=0, prop={'size':6})
ax2.legend(loc=0, prop={'size':6})
plt.tight_layout()
if(black == True):
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
	ax1.figure.savefig("Res_Args_w_sin_fit_black.png", facecolor=fig.get_facecolor(), transparent=True)
else:
	ax1.figure.savefig("Res_Args_w_sin_fit.png")
plt.clf()


