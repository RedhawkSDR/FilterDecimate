#!/usr/bin/env python
#
# This file is protected by Copyright. Please refer to the COPYRIGHT file
# distributed with this source distribution.
#
# This file is part of REDHAWK.
#
# REDHAWK is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# REDHAWK is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.
#

import unittest
import ossie.utils.testing
import os
from omniORB import any, CORBA
from ossie.cf import CF
from ossie.utils.bulkio import bulkio_data_helpers
from ossie.utils import sb
from ossie import properties

from numpy import sin, arange, pi
from scipy.signal import lfilter, firwin, butter
from pylab import figure, plot, grid, show, title
import scipy.fftpack
import numpy as np
import sys

import time, struct
from ctypes import *

def display(f, string):
	f.write(string)
	sys.stdout.write(string)

f = open('unit_test.log','w')


#########################################################################################
#							Define Test Cases											#
#########################################################################################
case1 = {'Filter_Type':'lowpass', 'Number_of_Taps':100, 'Center_Frequency':0, 'Bandwidth':10000, 'Gain':1, 'Output_Rate':256000}
case2 = {'Filter_Type':'lowpass', 'Number_of_Taps':50, 'Center_Frequency':0, 'Bandwidth':10000, 'Gain':1, 'Output_Rate':256000}
case3 = {'Filter_Type':'lowpass', 'Number_of_Taps':37, 'Center_Frequency':0, 'Bandwidth':10000, 'Gain':1, 'Output_Rate':256000}
case4 = {'Filter_Type':'lowpass', 'Number_of_Taps':100, 'Center_Frequency':1000, 'Bandwidth':10000, 'Gain':1, 'Output_Rate':256000}
case5 = {'Filter_Type':'lowpass', 'Number_of_Taps':47, 'Center_Frequency':0, 'Bandwidth':6000, 'Gain':1, 'Output_Rate':256000}
case6 = {'Filter_Type':'lowpass', 'Number_of_Taps':100, 'Center_Frequency':0, 'Bandwidth':10000, 'Gain':100, 'Output_Rate':256000}
case7 = {'Filter_Type':'lowpass', 'Number_of_Taps':100, 'Center_Frequency':0, 'Bandwidth':10000, 'Gain':1, 'Output_Rate':32000}
case8 = {'Filter_Type':'lowpass', 'Number_of_Taps':100, 'Center_Frequency':0, 'Bandwidth':8000, 'Gain':50, 'Output_Rate':128000}
case9 = {'Filter_Type':'lowpass', 'Number_of_Taps':73, 'Center_Frequency':20, 'Bandwidth':11000, 'Gain':36, 'Output_Rate':64000}
case10 = {'Filter_Type':'lowpass', 'Number_of_Taps':100, 'Center_Frequency':5000, 'Bandwidth':10000, 'Gain':1, 'Output_Rate':256000}

lowpass_cases = [case1, case2, case3, case4, case5, case6, case7, case8, case9, case10]

case11 = {'Filter_Type':'bandpass', 'Number_of_Taps':101, 'Center_Frequency':15000, 'Bandwidth':1000, 'Gain':1, 'Output_Rate':256000}
case12 = {'Filter_Type':'bandpass', 'Number_of_Taps':50, 'Center_Frequency':15000, 'Bandwidth':1000, 'Gain':1, 'Output_Rate':256000}
case13 = {'Filter_Type':'bandpass', 'Number_of_Taps':37, 'Center_Frequency':15000, 'Bandwidth':1000, 'Gain':1, 'Output_Rate':256000}
case14 = {'Filter_Type':'bandpass', 'Number_of_Taps':100, 'Center_Frequency':5000, 'Bandwidth':1000, 'Gain':1, 'Output_Rate':128000}
case15 = {'Filter_Type':'bandpass', 'Number_of_Taps':47, 'Center_Frequency':5000, 'Bandwidth':1000, 'Gain':1, 'Output_Rate':64000}
case16 = {'Filter_Type':'bandpass', 'Number_of_Taps':100, 'Center_Frequency':5000, 'Bandwidth':1000, 'Gain':100, 'Output_Rate':32000}
case17 = {'Filter_Type':'bandpass', 'Number_of_Taps':100, 'Center_Frequency':30000, 'Bandwidth':6000, 'Gain':1, 'Output_Rate':256000}
case18 = {'Filter_Type':'bandpass', 'Number_of_Taps':100, 'Center_Frequency':30000, 'Bandwidth':8000, 'Gain':50, 'Output_Rate':256000}
case19 = {'Filter_Type':'bandpass', 'Number_of_Taps':73, 'Center_Frequency':30000, 'Bandwidth':1000, 'Gain':36, 'Output_Rate':256000}
case20 = {'Filter_Type':'bandpass', 'Number_of_Taps':100, 'Center_Frequency':30000, 'Bandwidth':2000, 'Gain':25, 'Output_Rate':256000}

bandpass_cases = [case11, case12, case13, case14, case15, case16, case17, case18, case19, case20]

passed = {1:True, 2:True, 3:True, 4:True, 5:True, 6:True, 7:True, 8:True, 9:True, 10:True, 11:True, 12:True, 13:True, 14:True, 15:True, 16:True, 17:True, 18:True, 19:True, 20:True, 21:True, 22:True, 23:True, 24:True, 25:True, 26:True, 27:True, 28:True, 29:True, 30:True, 31:True, 32:True, 33:True, 34:True, 35:True, 36:True, 37:True, 38:True, 39:True, 40:True} 

passed_count = 0


#########################################################################################
#							UNIT TEST 1 - Lowpass - Real Data							#
#########################################################################################

display(f,'***************************************************************************\n')
display(f,'******************** Unit Test 1 - Lowpass - Real Data ********************\n')
display(f,'***************************************************************************\n\n')

#------------------------------------------------
# Create components and connections
#------------------------------------------------
display(f,'* Creating Components and Connections\n')
filt = sb.launch('../FilterDecimate.spd.xml',execparams={"DEBUG_LEVEL":5})

inputS=sb.DataSource(bytesPerPush=64);
outputS=sb.DataSink();
inputS.connect(filt,providesPortName='dataFloat_in')
filt.connect(outputS,providesPortName='floatIn')
	
#------------------------------------------------
# Start Components
#------------------------------------------------
display(f,'* Starting Component\n\n') 
sb.start()

case = 0
for ii in lowpass_cases:
	case += 1
	passed_count += 1
	display(f, '*** TEST CASE ' + str(passed_count) + ' ***\n')

	#------------------------------------------------
	# Create a signal for demonstration.
	#------------------------------------------------
	# 16384 samples of (5000Hz + 15000Hz + 30000 Hz) real signal at 256 kHz
	display(f,'* Creating Test Signal with Parameters\n')
	sample_rate = 256000.
	display(f,'sample_rate = 256000\n')
	nsamples = 16384.
	display(f,'nsamples = 10000\n')
	
	F_5KHz = 5000.
	display(f,'F_1KHz = 5000.\n')
	A_5KHz = 10.0
	display(f,'A_1KHz = 10.0\n')
	
	F_15KHz = 15000.
	display(f,'F_5KHz = 15000.\n')
	A_15KHz = 7.0
	display(f,'A_5KHz = 10.0\n')
	
	F_30KHz = 30000.
	display(f,'F_15KHz = 30000.\n')
	A_30KHz = 5.
	display(f,'A_15KHz = 5.0\n\n')
	
	t = arange(nsamples) / sample_rate
	t_dec = arange(nsamples/(sample_rate / ii['Output_Rate'])) / sample_rate
	signal = A_5KHz * sin(2*pi*F_5KHz*t) + A_15KHz * sin(2*pi*F_15KHz*t) + A_30KHz*sin(2*pi*F_30KHz*t)
	
	#------------------------------------------------
	# Configure component filter properties
	#------------------------------------------------
	display(f,'* FilterDecimate Component Configuration\n')
	filt.Filter_Type = ii['Filter_Type']
	display(f,'filt.Filter_Type = ' + ii['Filter_Type'] + '\n')
	filt.Number_of_Taps = ii['Number_of_Taps']
	display(f,'filt.Number_of_Taps = ' + str(ii['Number_of_Taps']) + '\n')
	filt.Center_Frequency = ii['Center_Frequency']
	display(f,'filt.Center_Frequency = ' + str(ii['Center_Frequency']) + '\n')
	filt.Bandwidth = ii['Bandwidth']
	display(f,'filt.Bandwidth = ' + str(ii['Bandwidth']) + '\n')
	filt.Gain = ii['Gain']
	display(f,'filt.Gain = ' + str(ii['Gain']) + '\n')
	filt.Output_Rate = ii['Output_Rate']
	display(f,'filt.Output_Rate = ' + str(ii['Output_Rate']) + '\n\n')

	#------------------------------------------------
	# Push Data to Component and get data
	#------------------------------------------------
	data = []
	for i in signal:
		data.append(float(i))
	
	inputS.push(data, False, "low_real", 256000.0,False)
	time.sleep(.5)
	received_data = outputS.getData()
	display(f,'* Filtered with FilterDecimates Component\n')

	#------------------------------------------------
	# Create a FIR filter and apply it to signal.
	#------------------------------------------------
	display(f,'* NumPy Lowpass Filter Configuration\n')
	
	nyq_rate = ii['Output_Rate']/2
	display(f,'nyq_rate = ' + str(ii['Output_Rate']) + '/2 = ' + str(nyq_rate) + '\n')
	cutoff_hz = ii['Bandwidth']
	display(f,'cutoff_hz = ' + str(ii['Bandwidth']) + '\n')
	numtaps = ii['Number_of_Taps']
	display(f,'numtaps = ' + str(ii['Number_of_Taps']) + '\n\n')
	
	# Use firwin to create lowpass filter coefficients
	fir_coeff = firwin(numtaps, cutoff_hz/nyq_rate, window=('kaiser',6))
	
	# Filter data using firwin generated coefficients
	py_data = lfilter(fir_coeff, 1, signal)

	# Decimate Data	
	py_dec = arange(nsamples/(sample_rate / ii['Output_Rate']))
	for i,j in enumerate(py_data):
		if i % int(sample_rate / ii['Output_Rate']) == 0:
			py_dec[i/int(sample_rate / ii['Output_Rate'])] = j

	display(f,'* Filtered & Decimated Signal with NumPy \n\n')

	#------------------------------------------------
	# Adjust filter gains and compare filtered data
	#------------------------------------------------
	total = 0;
	true_max = 0;
	length = len(received_data)

	comp_max = max(received_data[int(length*.25):int(length*.75)])
	py_max = max(py_dec[int(length*.25):int(length*.75)])
	gain_adjust = comp_max/py_max
	
	if comp_max > (py_max*gain_adjust):
		true_max = comp_max
	else:
		true_max = py_max*gain_adjust

	py_dec *= gain_adjust

	for j in range(int(length*.25), int(length*.75)):
		total += abs(py_dec[j]-received_data[j])
	
	avg = 100*(total/length)/true_max

	#------------------------------------------------
	# Check if variance in data is small enough
	#------------------------------------------------
	display(f,'* Comparing Output Data\n')
	passed[passed_count] = True
	if len(received_data) == 0 :
	    passed[passed_count] = False
	if avg > 20:
		passed[passed_count] = False
	
	#------------------------------------------------
	# Output Results
	#------------------------------------------------
	if passed[passed_count]:
		print "-------------------------------------------------------"
		print "LOWPASS TEST w/ Real Data ", case, ".......................",u'\u2714'	
		print "-------------------------------------------------------"
		f.write("-------------------------------------------------------\n")
		f.write("LOWPASS TEST w/ Real Data " + str(case) + "......................."+u'\u2714'.encode('utf8'))
		f.write("\n-------------------------------------------------------\n")
	else:
		print "-------------------------------------------------------"
		print "LOWPASS TEST w/ Real Data ", case, ".......................",u'\u2718'
		print "-------------------------------------------------------"
		f.write("-------------------------------------------------------\n")
		f.write("LOWPASS TEST w/ Real Data " + str(case) + "......................."+u'\u2718'.encode('utf8'))
		f.write("\n-------------------------------------------------------\n")
	
	display(f,'\n\n')

	#------------------------------------------------
	# Uncomment this section to see plots of data
	#------------------------------------------------
	#
	#display(f, 'Average variation between outputs is: ' + str(avg) + '% of max value\n\n' )
	#
	#FFT_in = abs(scipy.fft(signal))
	#FFT_npout = abs(scipy.fft(py_dec))
	#FFT_filt = abs(scipy.fft(received_data))
	#freqs = scipy.fftpack.fftfreq(signal.size, t[1]-t[0])
	#freqs_dec = scipy.fftpack.fftfreq(length, t[1]-t[0])
	#
	#figure(0)
	#title('Input Signal FFT')
	#plot(freqs,20*scipy.log10(FFT_in),'o')
	#
	#figure(1)
	#title('Component & NumPy Output FFT')
	#plot(freqs_dec,20*scipy.log10(FFT_npout),'x')
	#plot(freqs_dec,20*scipy.log10(FFT_filt),'+')
	#
	#figure(2)
	#title('Input Signal Time Plot')
	#plot(t, signal, 'y')
	#
	#figure(3)
	#title('Component & NumPy Output Plot')
	#plot(t_dec, py_dec, 'b')
	#plot(t_dec, received_data, 'g')
	#		
	#grid(True)
	#show()

#------------------------------------------------
# Stop Component
#------------------------------------------------
sb.stop()


#########################################################################################
#							UNIT TEST 2 - Bandpass - Real Data							#
#########################################################################################

display(f,'****************************************************************************\n')
display(f,'******************** Unit Test 2 - Bandpass - Real Data ********************\n')
display(f,'****************************************************************************\n\n')

#------------------------------------------------
# Create components and connections
#------------------------------------------------
display(f,'* Creating Components and Connections\n')
filt = sb.launch('../FilterDecimate.spd.xml',execparams={"DEBUG_LEVEL":5})

inputS=sb.DataSource(bytesPerPush=64);
outputS=sb.DataSink();
inputS.connect(filt,providesPortName='dataFloat_in')
filt.connect(outputS,providesPortName='floatIn')

#------------------------------------------------
# Start Components
#------------------------------------------------
display(f,'* Starting Component\n\n') 
sb.start()

case = 0
for ii in bandpass_cases:
	case += 1
	passed_count += 1
	display(f, '*** TEST CASE ' + str(passed_count) + ' ***\n')

	#------------------------------------------------
	# Create a signal for demonstration.
	#------------------------------------------------
	# 16384 samples of (5000Hz + 15000Hz + 30000 Hz) real signal at 256 kHz
	display(f,'* Creating Test Signal with Parameters\n')
	sample_rate = 256000.
	display(f,'sample_rate = 256000\n')
	nsamples = 16384.
	display(f,'nsamples = 10000\n')
	
	F_5KHz = 5000.
	display(f,'F_1KHz = 5000.\n')
	A_5KHz = 10.0
	display(f,'A_1KHz = 10.0\n')
	
	F_15KHz = 15000.
	display(f,'F_5KHz = 15000.\n')
	A_15KHz = 7.0
	display(f,'A_5KHz = 10.0\n')
	
	F_30KHz = 30000.
	display(f,'F_15KHz = 30000.\n')
	A_30KHz = 5.
	display(f,'A_15KHz = 5.0\n\n')
	
	t = arange(nsamples) / sample_rate
	t_dec = arange(nsamples/(sample_rate / ii['Output_Rate'])) / sample_rate
	signal = A_5KHz * sin(2*pi*F_5KHz*t) + A_15KHz * sin(2*pi*F_15KHz*t) + A_30KHz*sin(2*pi*F_30KHz*t)
	
	#------------------------------------------------
	# Configure component filter properties
	#------------------------------------------------
	display(f,'* FilterDecimate Component Configuration\n')
	filt.Filter_Type = ii['Filter_Type']
	display(f,'filt.Filter_Type = ' + ii['Filter_Type'] + '\n')
	filt.Number_of_Taps = ii['Number_of_Taps']
	display(f,'filt.Number_of_Taps = ' + str(ii['Number_of_Taps']) + '\n')
	filt.Center_Frequency = ii['Center_Frequency']
	display(f,'filt.Center_Frequency = ' + str(ii['Center_Frequency']) + '\n')
	filt.Bandwidth = ii['Bandwidth']
	display(f,'filt.Bandwidth = ' + str(ii['Bandwidth']) + '\n')
	filt.Gain = ii['Gain']
	display(f,'filt.Gain = ' + str(ii['Gain']) + '\n')
	filt.Output_Rate = ii['Output_Rate']
	display(f,'filt.Output_Rate = ' + str(ii['Output_Rate']) + '\n\n')
	
	#------------------------------------------------
	# Push Data to Component and get data
	#------------------------------------------------
	data = []
	
	for i in signal:
		data.append(float(i))
	
	inputS.push(data, False, "bandpass_real", 256000.0,False)
	time.sleep(.5)
	received_data = outputS.getData()
	display(f,'* Filtered with FilterDecimate Component\n')
	
	#------------------------------------------------
	# Create a bandpass filter and apply it to signal.
	#------------------------------------------------
	display(f,'* NumPy Bandpass Filter Configuration\n')
	
	nyq_rate = sample_rate/2
	display(f,'nyq_rate = sample_rate/2\n')
	cutoff1_hz = ii['Center_Frequency'] - (ii['Bandwidth']/2)
	display(f,'cutoff1_hz = ' + str(ii['Center_Frequency'] - (ii['Bandwidth']/2)) + '\n')
	cutoff2_hz = ii['Center_Frequency'] + (ii['Bandwidth']/2)
	display(f,'cutoff2_hz = ' + str(ii['Center_Frequency'] + (ii['Bandwidth']/2)) + '\n')
	
	#Use firwin to create bandpass filter coefficients
	b, a = butter(1, [cutoff1_hz/nyq_rate, cutoff2_hz/nyq_rate], btype='band')
	
	#Filter data using firwin generated coefficients
	py_data = lfilter(b, a, signal)
	
	#Decimate Data	
	py_dec = arange(nsamples/(sample_rate/ii['Output_Rate']))
	for i,j in enumerate(py_data):
		if i % int(sample_rate / ii['Output_Rate']) == 0:
			py_dec[i/int(sample_rate/ii['Output_Rate'])] = j
	
	display(f,'* Filtered & Decimated Signal with NumPy \n\n')
	
	#------------------------------------------------
	# Shift outputs to align properly for comparison
	#------------------------------------------------
	last_sign = 0
	sign = 0
	count = 0
	py_i = 0
	if py_dec[0] > 0:
		last_sign = 1
	else:
		last_sign = -1
	
	for i, j in enumerate(py_dec):
		if j > 0:
			sign = 1
		else:
			sign = -1
		
		if last_sign != sign:
			count += 1
			last_sign = sign
	
			if count == 10:
				py_i = i
	
	last_sign = 0
	sign = 0
	count = 0
	comp_i = 0
	if py_dec[0] > 0:
		last_sign = 1
	else:
		last_sign = -1
	
	for i, j in enumerate(received_data):
		if j > 0:
			sign = 1
		else:
			sign = -1
		
		if last_sign != sign:
			count += 1
			last_sign = sign
	
			if count == 10:
				comp_i = i
	
	difference = py_i - comp_i
	py_dec = np.roll(py_dec,-difference)

	#------------------------------------------------
	# Adjust filter gains and compare filtered data
	#------------------------------------------------
	total = 0;
	true_max = 0;
	length = len(py_dec)

	comp_max = max(received_data[int(length*.25):int(length*.75)])
	py_max = max(py_dec[int(length*.25):int(length*.75)])
	gain_adjust = comp_max/py_max
	
	if comp_max > (py_max*gain_adjust):
		true_max = comp_max
	else:
		true_max = py_max*gain_adjust

	py_dec *= gain_adjust

	for j in range(int(length*.25), int(length*.75)):
		total += abs(py_dec[j]-received_data[j])
	
	avg = 100*(total/length)/true_max
	
	#------------------------------------------------
	# Check if variance in data is small enough
	#------------------------------------------------
	display(f,'* Comparing Output Data\n')
	passed[passed_count] = True
	if len(received_data) == 0 :
	    passed[passed_count] = False
	if avg >20:
		passed[passed_count] = False
	
	#------------------------------------------------
	# Output Results
	#------------------------------------------------
	if passed[passed_count]:
		print "--------------------------------------------------------"
		print "BANDPASS TEST w/ Real Data ", case, ".......................",u'\u2714'	
		print "--------------------------------------------------------"
		f.write("--------------------------------------------------------\n")
		f.write("BANDPASS TEST w/ Real Data " + str(case) + "......................."+u'\u2714'.encode('utf8'))
		f.write("\n--------------------------------------------------------\n")
	else:
		print "--------------------------------------------------------"
		print "BANDPASS TEST w/ Real Data ", case, ".......................",u'\u2718'
		print "--------------------------------------------------------"
		f.write("--------------------------------------------------------\n")
		f.write("BANDPASS TEST w/ Real Data " + str(case) + "......................."+u'\u2718'.encode('utf8'))
		f.write("\n--------------------------------------------------------\n")
	
	display(f,'\n\n')

	#------------------------------------------------
	# Uncomment this section to see plots of data
	#------------------------------------------------
	#
	#display(f, 'Average variation between outputs is: ' + str(avg) + '% of max value\n\n' )
	#
	#FFT_in = abs(scipy.fft(signal))
	#FFT_npout = abs(scipy.fft(py_dec))
	#FFT_filt = abs(scipy.fft(received_data))
	#freqs = scipy.fftpack.fftfreq(signal.size, t[1]-t[0])
	#freqs_dec = scipy.fftpack.fftfreq(length, t[1]-t[0])
	#
	#figure(0)
	#title('Input Signal FFT')
	#plot(freqs,20*scipy.log10(FFT_in),'o')
	#
	#figure(1)
	#title('Component & NumPy Output FFT')
	#plot(freqs_dec,20*scipy.log10(FFT_npout),'x')
	#plot(freqs_dec,20*scipy.log10(FFT_filt),'+')
	#
	#figure(2)
	#title('Input Signal Time Plot')
	#plot(t, signal, 'y')
	#
	#
	#figure(3)
	#title('Component & NumPy Output Plot')
	#plot(t_dec, py_dec, 'b')
	#plot(t_dec, received_data, 'g')
	#		
	#grid(True)
	#show()

#------------------------------------------------
# Stop Component
#------------------------------------------------
sb.stop()


#########################################################################################
#							UNIT TEST 3 - Lowpass - Complex Data						#
#########################################################################################

display(f,'******************************************************************************\n')
display(f,'******************** Unit Test 3 - Lowpass - Complex Data ********************\n')
display(f,'******************************************************************************\n\n')

#------------------------------------------------
# Create components and connections
#------------------------------------------------
display(f,'* Creating Components and Connections\n')
filt = sb.launch('../FilterDecimate.spd.xml',execparams={"DEBUG_LEVEL":5})

inputS=sb.DataSource(bytesPerPush=64);
outputS=sb.DataSink();
inputS.connect(filt,providesPortName='dataFloat_in')
filt.connect(outputS,providesPortName='floatIn')

#------------------------------------------------
# Start Components
#------------------------------------------------
display(f,'* Starting Component\n\n') 
sb.start()

case = 0
for ii in lowpass_cases:
	case += 1
	passed_count += 1
	display(f, '*** TEST CASE ' + str(passed_count) + ' ***\n')

	#------------------------------------------------
	# Create a signal for demonstration.
	#------------------------------------------------
	# 16384 samples of (5000Hz + 15000Hz + 30000 Hz) complex signal at 256 kHz
	display(f,'* Creating Test Signal with Parameters\n')
	sample_rate = 256000.
	display(f,'sample_rate = 256000\n')
	nsamples = 16384.
	display(f,'nsamples = 10000\n')
	
	F_5KHz = 5000.
	display(f,'F_1KHz = 5000.\n')
	A_5KHz = 10.0
	display(f,'A_1KHz = 10.0\n')
	
	F_15KHz = 15000.
	display(f,'F_5KHz = 15000.\n')
	A_15KHz = 7.0
	display(f,'A_5KHz = 10.0\n')
	
	F_30KHz = 30000.
	display(f,'F_15KHz = 30000.\n')
	A_30KHz = 5.
	display(f,'A_15KHz = 5.0\n\n')
	
	t = arange(nsamples) / sample_rate
	t_dec = arange(nsamples/(sample_rate / ii['Output_Rate'])) / sample_rate
	signal = A_5KHz * np.exp(1j*2*pi*F_5KHz*t) + A_15KHz * np.exp(1j*2*pi*F_15KHz*t) + A_30KHz * np.exp(1j*2*pi*F_30KHz*t)

	signal_redhawk = np.zeros(len(signal)*2)
	for jj in np.arange(0, len(signal)*2, 2):
		signal_redhawk[jj] = signal[jj/2].real
		signal_redhawk[jj+1] = signal[jj/2].imag

	#------------------------------------------------
	# Configure component filter properties
	#------------------------------------------------
	display(f,'*  FilterDecimate Component Configuration\n')
	filt.Filter_Type = ii['Filter_Type']
	display(f,'filt.Filter_Type = ' + ii['Filter_Type'] + '\n')
	filt.Number_of_Taps = ii['Number_of_Taps']
	display(f,'filt.Number_of_Taps = ' + str(ii['Number_of_Taps']) + '\n')
	filt.Center_Frequency = ii['Center_Frequency']
	display(f,'filt.Center_Frequency = ' + str(ii['Center_Frequency']) + '\n')
	filt.Bandwidth = ii['Bandwidth']
	display(f,'filt.Bandwidth = ' + str(ii['Bandwidth']) + '\n')
	filt.Gain = ii['Gain']
	display(f,'filt.Gain = ' + str(ii['Gain']) + '\n')
	filt.Output_Rate = ii['Output_Rate']
	display(f,'filt.Output_Rate = ' + str(ii['Output_Rate']) + '\n\n')

	#------------------------------------------------
	# Push Data to Component and get data
	#------------------------------------------------
	data = []
	for i in signal_redhawk:
		data.append(float(i))
	
	inputS.push(data, False, "low_complex", 256000.0, True)
	time.sleep(1.5)

	received_data_redhawk = arange((nsamples*2) / (sample_rate / ii['Output_Rate']),dtype=np.complex128)
	received_data_redhawk = outputS.getData()
	received_data = np.ndarray(len(received_data_redhawk)/2, dtype=complex)

	for jj in np.arange(0, len(received_data_redhawk), 2):
		received_data[jj/2] = received_data_redhawk[jj] + 1j*received_data_redhawk[jj+1]
	
	display(f,'* Filtered with FilterDecimate Component\n')

	#------------------------------------------------
	# Create a FIR filter and apply it to signal.
	#------------------------------------------------
	display(f,'* NumPy Lowpass Filter Configuration\n')
	
	nyq_rate = ii['Output_Rate']/2
	display(f,'nyq_rate = ' + str(ii['Output_Rate']) + '/2 = ' + str(nyq_rate) + '\n')
	cutoff_hz = ii['Bandwidth']
	display(f,'cutoff_hz = ' + str(ii['Bandwidth']) + '\n')
	numtaps = ii['Number_of_Taps']
	display(f,'numtaps = ' + str(ii['Number_of_Taps']) + '\n\n')
	
	# Use firwin to create lowpass filter coefficients
	fir_coeff = firwin(numtaps, cutoff_hz/nyq_rate, window=('kaiser',6))
	
	# Filter data using firwin generated coefficients
	py_data = lfilter(fir_coeff, 1, signal)

	# Decimate Data	
	py_dec = arange(nsamples/(sample_rate / ii['Output_Rate']),dtype=np.complex128)
	for i,j in enumerate(py_data):
		if i % int(sample_rate / ii['Output_Rate']) == 0:
			py_dec[i/int(sample_rate / ii['Output_Rate'])] = j

	display(f,'* Filtered & Decimated Signal with NumPy \n\n')

	#------------------------------------------------
	# Adjust filter gains and compare filtered data
	#------------------------------------------------
	total = 0;
	true_max = 0;
	length = len(received_data)

	comp_max = max(np.real(received_data[int(length*.25):int(length*.75)]))
	py_max = max(np.real(py_dec[int(length*.25):int(length*.75)]))
	gain_adjust = comp_max/py_max
	
	if comp_max > (py_max*gain_adjust):
		true_max = comp_max
	else:
		true_max = py_max*gain_adjust

	py_dec *= gain_adjust

	for j in range(int(length*.25), int(length*.75)):
		total += abs(py_dec[j]-received_data[j])

	avg = 100*(total/length)/true_max

	#------------------------------------------------
	# Check if variance in data is small enough
	#------------------------------------------------
	display(f,'* Comparing Output Data\n')
	passed[passed_count] = True
	if len(received_data) == 0 :
	    passed[passed_count] = False
	if avg > 20:
		passed[passed_count] = False

	#------------------------------------------------
	# Output Results
	#------------------------------------------------
	if passed[passed_count]:
		print "----------------------------------------------------------"
		print "LOWPASS TEST w/ Complex Data ", case, ".......................",u'\u2714'	
		print "----------------------------------------------------------"
		f.write("----------------------------------------------------------\n")
		f.write("LOWPASS TEST w/ Complex Data " + str(case) + "......................."+u'\u2714'.encode('utf8'))
		f.write("\n----------------------------------------------------------\n")
	else:
		print "----------------------------------------------------------"
		print "LOWPASS TEST w/ Complex Data ", case, ".......................",u'\u2718'
		print "----------------------------------------------------------"
		f.write("----------------------------------------------------------\n")
		f.write("LOWPASS TEST w/ Complex Data " + str(case) + "......................."+u'\u2718'.encode('utf8'))
		f.write("\n----------------------------------------------------------\n")
	
	display(f,'\n\n')

	#------------------------------------------------
	# Uncomment this section to see plots of data
	#------------------------------------------------
	#
	#display(f, 'Average variation between outputs is: ' + str(avg) + '% of max value\n\n' )
	#
	#FFT_in = abs(scipy.fft(signal))
	#FFT_npout = abs(scipy.fft(py_dec))
	#FFT_filt = abs(scipy.fft(received_data))
	#freqs = scipy.fftpack.fftfreq(signal.size, t[1]-t[0])
	#freqs_dec = scipy.fftpack.fftfreq(length, t[1]-t[0])
	#
	#figure(0)
	#title('Input Signal FFT')
	#plot(freqs,20*scipy.log10(FFT_in),'o')
	#
	#figure(1)
	#title('Component & NumPy Output FFT')
	#plot(freqs_dec,20*scipy.log10(FFT_npout),'x')
	#plot(freqs_dec,20*scipy.log10(FFT_filt),'+')
	#
	#figure(2)
	#title('Input Signal Time Plot')
	#plot(t, signal, 'y')
	#
	#
	#figure(3)
	#title('Component & NumPy Output Plot')
	#plot(t_dec, np.real(py_dec), 'b')
	#plot(t_dec, np.real(received_data), 'g')
	#		
	#figure(4)
	#title('Component & NumPy Output Plot')
	#plot(t_dec, np.imag(py_dec), 'b')
	#plot(t_dec, np.imag(received_data), 'g')
	#
	#grid(True)
	#show()

#------------------------------------------------
# Stop Component
#------------------------------------------------
sb.stop()


#########################################################################################
#							UNIT TEST 4 - Bandpass - Complex Data						#
#########################################################################################

display(f,'*******************************************************************************\n')
display(f,'******************** Unit Test 4 - Bandpass - Complex Data ********************\n')
display(f,'*******************************************************************************\n\n')

#------------------------------------------------
# Create components and connections
#------------------------------------------------
display(f,'* Creating Components and Connections\n')
filt = sb.launch('../FilterDecimate.spd.xml',execparams={"DEBUG_LEVEL":5})

inputS=sb.DataSource(bytesPerPush=64);
outputS=sb.DataSink();
inputS.connect(filt,providesPortName='dataFloat_in')
filt.connect(outputS,providesPortName='floatIn')

#------------------------------------------------
# Start Components
#------------------------------------------------
display(f,'* Starting Component\n\n') 
sb.start()

case = 0
for ii in bandpass_cases:
	case += 1
	passed_count += 1
	display(f, '*** TEST CASE ' + str(passed_count) + ' ***\n')

	#------------------------------------------------
	# Create a signal for demonstration.
	#------------------------------------------------
	# 16384 samples of (5000Hz + 15000Hz + 30000 Hz) complex signal at 256 kHz
	display(f,'* Creating Test Signal with Parameters\n')
	sample_rate = 256000.
	display(f,'sample_rate = 256000\n')
	nsamples = 16384.
	display(f,'nsamples = 10000\n')
	
	F_5KHz = 5000.
	display(f,'F_1KHz = 5000.\n')
	A_5KHz = 10.0
	display(f,'A_1KHz = 10.0\n')
	
	F_15KHz = 15000.
	display(f,'F_5KHz = 15000.\n')
	A_15KHz = 7.0
	display(f,'A_5KHz = 10.0\n')
	
	F_30KHz = 30000.
	display(f,'F_15KHz = 30000.\n')
	A_30KHz = 5.
	display(f,'A_15KHz = 5.0\n\n')
	
	t = arange(nsamples) / sample_rate
	t_dec = arange(nsamples/(sample_rate / ii['Output_Rate'])) / sample_rate
	signal = A_5KHz * np.exp(1j*2*pi*F_5KHz*t) + A_15KHz * np.exp(1j*2*pi*F_15KHz*t) + A_30KHz * np.exp(1j*2*pi*F_30KHz*t)

	signal_redhawk = np.zeros(len(signal)*2)
	for jj in np.arange(0, len(signal)*2, 2):
		signal_redhawk[jj] = signal[jj/2].real
		signal_redhawk[jj+1] = signal[jj/2].imag

	#------------------------------------------------
	# Configure component filter properties
	#------------------------------------------------
	display(f,'*  FilterDecimate Component Configuration\n')
	filt.Filter_Type = ii['Filter_Type']
	display(f,'filt.Filter_Type = ' + ii['Filter_Type'] + '\n')
	filt.Number_of_Taps = ii['Number_of_Taps']
	display(f,'filt.Number_of_Taps = ' + str(ii['Number_of_Taps']) + '\n')
	filt.Center_Frequency = ii['Center_Frequency']
	display(f,'filt.Center_Frequency = ' + str(ii['Center_Frequency']) + '\n')
	filt.Bandwidth = ii['Bandwidth']
	display(f,'filt.Bandwidth = ' + str(ii['Bandwidth']) + '\n')
	filt.Gain = ii['Gain']
	display(f,'filt.Gain = ' + str(ii['Gain']) + '\n')
	filt.Output_Rate = ii['Output_Rate']
	display(f,'filt.Output_Rate = ' + str(ii['Output_Rate']) + '\n\n')

	#------------------------------------------------
	# Push Data to Component and get data
	#------------------------------------------------
	data = []
	for i in signal_redhawk:
		data.append(float(i))
	
	inputS.push(data, False, "low_complex", 256000.0, True)
	time.sleep(1.5)

	received_data_redhawk = arange((nsamples*2) / (sample_rate / ii['Output_Rate']),dtype=np.complex128)
	received_data_redhawk = outputS.getData()
	received_data = np.ndarray(len(received_data_redhawk)/2, dtype=complex)

	for jj in np.arange(0, len(received_data_redhawk), 2):
		received_data[jj/2] = received_data_redhawk[jj] + 1j*received_data_redhawk[jj+1]
	
	display(f,'* Filtered with FilterDecimate Component\n')

	#-----------------------------------------------------
	# Create a FIR bandpass filter and apply it to signal.
	#-----------------------------------------------------
	display(f,'* NumPy Lowpass Filter Configuration\n')
	
	nyq_rate = ii['Output_Rate']/2
	display(f,'nyq_rate = ' + str(ii['Output_Rate']) + '/2 = ' + str(nyq_rate) + '\n')
	cutoff_hz = ii['Bandwidth']
	display(f,'cutoff_hz = ' + str(ii['Bandwidth']) + '\n')
	numtaps = ii['Number_of_Taps']
	display(f,'numtaps = ' + str(ii['Number_of_Taps']) + '\n\n')
	
	# Use firwin to create lowpass filter coefficients
	fir_coeff = firwin(numtaps, cutoff_hz/nyq_rate, window=('kaiser',6))

	t = arange(numtaps)
	shift = np.exp(1j * 2.0 * np.pi * ii['Center_Frequency'] * t / sample_rate)
	fir_coeff = np.real(fir_coeff*shift)

	# Filter data using firwin generated coefficients
	py_data = lfilter(fir_coeff, 1, signal)

	# Decimate Data	
	py_dec = arange(nsamples/(sample_rate / ii['Output_Rate']),dtype=np.complex128)
	for i,j in enumerate(py_data):
		if i % int(sample_rate / ii['Output_Rate']) == 0:
			py_dec[i/int(sample_rate / ii['Output_Rate'])] = j
	
	#------------------------------------------------
	# Shift outputs to align properly for comparison
	#------------------------------------------------
	last_sign = 0
	sign = 0
	count = 0
	py_i = 0
	if py_dec[0] > 0:
		last_sign = 1
	else:
		last_sign = -1
	
	for i, j in enumerate(py_dec):
		if j > 0:
			sign = 1
		else:
			sign = -1
		
		if last_sign != sign:
			count += 1
			last_sign = sign
	
			if count == 10:
				py_i = i
	
	last_sign = 0
	sign = 0
	count = 0
	comp_i = 0
	if py_dec[0] > 0:
		last_sign = 1
	else:
		last_sign = -1
	
	for i, j in enumerate(received_data):
		if j > 0:
			sign = 1
		else:
			sign = -1
		
		if last_sign != sign:
			count += 1
			last_sign = sign
	
			if count == 10:
				comp_i = i
	
	difference = py_i - comp_i
	py_dec = np.roll(py_dec,-difference)

	#------------------------------------------------
	# Adjust filter gains and compare filtered data
	#------------------------------------------------
	total = 0;
	true_max = 0;
	length = len(received_data)

	comp_max = max(np.real(received_data[int(length*.25):int(length*.75)]))
	py_max = max(np.real(py_dec[int(length*.25):int(length*.75)]))
	gain_adjust = comp_max/py_max
	
	if comp_max > (py_max*gain_adjust):
		true_max = comp_max
	else:
		true_max = py_max*gain_adjust

	py_dec *= gain_adjust

	for j in range(int(length*.25), int(length*.75)):
		total += abs(py_dec[j]-received_data[j])

	avg = 100*(total/length)/true_max
	
	#------------------------------------------------
	# Check if variance in data is small enough
	#------------------------------------------------
	display(f,'* Comparing Output Data\n')
	passed[passed_count] = True
	if len(received_data) == 0 :
	    passed[passed_count] = False
	if avg >20:
		passed[passed_count] = False
	
	#------------------------------------------------
	# Output Results
	#------------------------------------------------
	if passed[passed_count]:
		print "-----------------------------------------------------------"
		print "BANDPASS TEST w/ Complex Data ", case, ".......................",u'\u2714'	
		print "-----------------------------------------------------------"
		f.write("-----------------------------------------------------------\n")
		f.write("BANDPASS TEST w/ Complex Data " + str(case) + "......................."+u'\u2714'.encode('utf8'))
		f.write("\n-----------------------------------------------------------\n")
	else:
		print "-----------------------------------------------------------"
		print "BANDPASS TEST w/ Complex Data ", case, ".......................",u'\u2718'
		print "-----------------------------------------------------------"
		f.write("-----------------------------------------------------------\n")
		f.write("BANDPASS TEST w/ Complex Data " + str(case) + "......................."+u'\u2718'.encode('utf8'))
		f.write("\n-----------------------------------------------------------\n")
	
	display(f,'\n\n')

	#------------------------------------------------
	# Uncomment this section to see plots of data
	#------------------------------------------------
	#
	#display(f, 'Average variation between outputs is: ' + str(avg) + '% of max value\n\n' )
	#
	#FFT_in = abs(scipy.fft(signal))
	#FFT_npout = abs(scipy.fft(py_dec))
	#FFT_filt = abs(scipy.fft(received_data))
	#freqs = scipy.fftpack.fftfreq(signal.size, t[1]-t[0])
	#freqs_dec = scipy.fftpack.fftfreq(length, t[1]-t[0])
	#
	#figure(0)
	#title('Input Signal FFT')
	#plot(freqs,20*scipy.log10(FFT_in),'o')
	#
	#figure(1)
	#title('Component & NumPy Output FFT')
	#plot(freqs_dec,20*scipy.log10(FFT_npout),'x')
	#plot(freqs_dec,20*scipy.log10(FFT_filt),'+')
	#
	#figure(2)
	#title('Input Signal Time Plot')
	#plot(t, signal, 'y')
	#
	#
	#figure(3)
	#title('Component & NumPy Output Plot')
	#plot(t_dec, np.real(py_dec), 'b')
	#plot(t_dec, np.real(received_data), 'g')
	#		
	#figure(4)
	#title('Component & NumPy Output Plot')
	#plot(t_dec, np.imag(py_dec), 'b')
	#plot(t_dec, np.imag(received_data), 'g')
	#
	#grid(True)
	#show()

#------------------------------------------------
# Stop Component
#------------------------------------------------
sb.stop()


#########################################################################################
#							Print Final Results											#
#########################################################################################
display(f, '\n\n###########################################################\n')
display(f, '#                   FINAL RESULTS                         #\n')
display(f, '###########################################################\n\n')

for i in range(len(passed)+1):
	if i > 0 and i <= 10:
		if passed[i]:
			print "LOWPASS TEST w/ Real Data ", i, ".......................",u'\u2714'
			f.write("\nLOWPASS TEST w/ Real Data " + str(i) + "......................."+u'\u2714'.encode('utf8'))
		else:
			print "LOWPASS TEST w/ Real Data ", i, ".......................",u'\u2718'
			f.write("\nLOWPASS TEST w/ Real Data " + str(i) + "......................."+u'\u2718'.encode('utf8'))
	if i > 10 and i <= 20:
		if passed[i]:
			print "BANDPASS TEST w/ Real Data ", i-10, ".......................",u'\u2714'
			f.write("\nBANDPASS TEST w/ Real Data " + str(i-10) + "......................."+u'\u2714'.encode('utf8'))
		else:
			print "BANDPASS TEST w/ Real Data ", i-10, ".......................",u'\u2718'
			f.write("\nBANDPASS TEST w/ Real Data " + str(i-10) + "......................."+u'\u2718'.encode('utf8'))
	if i > 20 and i <= 30:
		if passed[i]:
			print "LOWPASS TEST w/ Complex Data ", i-20, ".......................",u'\u2714'
			f.write("\nLOWPASS TEST w/ Complex Data " + str(i-20) + "......................."+u'\u2714'.encode('utf8'))		
		else:
			print "LOWPASS TEST w/ Complex Data ", i-20, ".......................",u'\u2718'
			f.write("\nLOWPASS TEST w/ Complex Data " + str(i-20) + "......................."+u'\u2718'.encode('utf8'))
	if i > 30 and i <= 40:
		if passed[i]:
			print "BANDPASS TEST w/ Complex Data ", i-30, ".......................",u'\u2714'	
			f.write("\nBANDPASS TEST w/ Complex Data " + str(i-30) + "......................."+u'\u2714'.encode('utf8'))
		else:
			print "BANDPASS TEST w/ Complex Data ", i-30, ".......................",u'\u2718'	
			f.write("\nBANDPASS TEST w/ Complex Data " + str(i-30) + "......................."+u'\u2718'.encode('utf8'))
	
display(f, '\n###########################################################\n\n')
display(f, 'Results written to "unit_test.log"\n\n')

f.close()
	
