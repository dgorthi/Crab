#! /usr/bin/python
import struct
import time
import pylab
import time
import numpy as np
import os
from pprint import pprint
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from filterbank_utils import *
import sys
from numpy.fft import fft, fftshift
from random import randint
from numpy import absolute


power = np.empty([0])
def convert_direct(filename,bottom=218, top=1242):
	""" Read data in Deepthi's file format 

	This uses np.fromfile and numpy arrays instead of struct,
	which should make it faster than the original version

	"""
	# FILENAME EDITED 9-27 BY SAM, REPLACE BEFORE RUNNING
	# filename= 'data_2017-09-27_15-30-41_723342346'
	n_data = 512*4*1024*18  # num_channels_per_packet*4_pols_per_channels*num_packets_per_buffer = Number of data points per buffer
	n_flag = 1024*18      # Number of flags per buffer
	print("Reading in binary data...")
	with open(filename, 'r') as fp:
		file_size = os.stat(filename).st_size
		no_cpy_buf = int(file_size/(n_data*2)) # file_size/(n_data_points_per_buffer*2_bytes_per_data_point) = no_cpy_buf

		raw_data = np.zeros((no_cpy_buf, n_data), dtype='<H')
		for ii in range(no_cpy_buf):    
			#copy = struct.unpack(struct_fmt, fp.read(struct_len))
			data = np.fromfile(fp, dtype='<H', count=n_data)
			#print(len(data))
			flag = np.fromfile(fp, dtype='<B', count=n_flag)
			zeros = np.zeros(2048, dtype='<H')
			print(data)
			"""
			for y in range(len(flag)):
				if flag[y] == 1:
					data[2048*y:2048*(y+1)] = zeros
			"""
			raw_data[ii] = data
	#noise = data
	print("Reading binary data done")
	total_data = raw_data.flatten()
	xx      = total_data[0::4]
	yy      = total_data[1::4]
	global power
	power = np.sqrt((np.asarray(xx)/3.5)**2+np.asarray(yy)**2)
	print(power[0:1000])
	power = power.reshape(len(power)/1536,1536)
	power = np.asarray([l[bottom:top] for l in power])
	power = power.flatten()
	power = power/36
	idx=np.where(power>255)
	power[idx] = 255
	power = np.uint8(power)

	return power

#def get_data_2(filename='data_2017-09-27_15-30-41_723342346',bottom=2200, top=2800):
def get_data_2(filename,bottom=218, top=1242):
	""" Read data in Deepthi's file format 

	This uses np.fromfile and numpy arrays instead of struct,
	which should make it faster than the original version

	"""
	# FILENAME EDITED 9-27 BY SAM, REPLACE BEFORE RUNNING
	# filename= 'data_2017-09-27_15-30-41_723342346'
	n_data = 512*4*1024*18  # num_channels_per_packet*4_pols_per_channels*num_packets_per_buffer = Number of data points per buffer
	n_flag = 1024*18      # Number of flags per buffer
	print("Reading in binary data...")
	with open(filename, 'r') as fp:
		file_size = os.stat(filename).st_size
		no_cpy_buf = int(file_size/(n_data*2)) # file_size/(n_data_points_per_buffer*2_bytes_per_data_point) = no_cpy_buf

		raw_data = np.zeros((no_cpy_buf, n_data), dtype='<H')
		for ii in range(no_cpy_buf):    
			#copy = struct.unpack(struct_fmt, fp.read(struct_len))
			data = np.fromfile(fp, dtype='<H', count=n_data)
			#print(len(data))
			flag = np.fromfile(fp, dtype='<B', count=n_flag)
			zeros = np.zeros(2048, dtype='<H')
			for y in range(len(flag)):
				if flag[y] == 1:
					data[2048*y:2048*(y+1)] = zeros
			raw_data[ii] = data
	#noise = data
	print("Reading binary data done")
	total_data = raw_data.flatten()
	xx      = total_data[0::4]
	yy      = total_data[1::4]
	global power
	print("Cutting unwanted channels...")
	power = np.sqrt((np.asarray(xx)/3.5)**2+np.asarray(yy)**2)
	#power = power[0:1536*1024*10]
	n_spectra = no_cpy_buf*n_data/4/1536
	power = power.reshape(len(power)/1536,1536)
	power = np.asarray([l[bottom:top] for l in power])
	power_fold = np.sum(power, axis = 1)
	power_fold_mean = 0
	"""
	for i in range(len(power_fold)):
		if i < len(power_fold)-1024:
			power_fold_mean = np.mean(power_fold[i:i+1024])	
		power[i,:]= power[i,:]-power_fold_mean/1024
		power_fold[i] = power_fold[i]-power_fold_mean
	"""
	#power = power.flatten()
	#power = power.reshape(len(power)/1024,1024)
	#power = np.transpose(power)
	print("Cutting unwanted channels done")
	#power = radar_removal(power, n_spectra)
	"""
	f, axarr = plt.subplots(2, sharex=True)
	g, xarr = plt.subplots(2, sharex= True)
	h, arr = plt.subplots(2, sharex= True)
	axarr[0].plot(power[1017],'r')
	xarr[0].plot(np.array(power_fold).flatten(),'r')
	arr[0].pcolormesh(power,cmap='coolwarm')
	"""
	n_channels = 1024
	thresh_arg = 3
	n_sample_thresh = 1024

	print('Broadband Start')
	start = time.time()
	for i in range(1024/n_channels):
		RFI_removal(thresh_1=3,thresh_2=1.5, thresh_3=.75, thresh_arg=3, n_sample_thresh_2=80, n_sample_thresh_3=160, start_channel=n_channels*i, n_channel=n_channels, n_sample_thresh=1024)
	end = time.time()
	print('Broadband Done')
	print('Time for Broadband: ' + str(end-start))

	n_channels = 1
	thresh_arg = 1
	n_sample_thresh = 1024
	print('Narrowband Start')
	start = time.time()
	for i in range(1024/n_channels):
		#print('n_channels*i'+str(n_channels*i))
		RFI_removal(thresh_1=2.5,thresh_2=1.5, thresh_3=.75, thresh_arg=thresh_arg, n_sample_thresh_2=80, n_sample_thresh_3=160, start_channel=n_channels*i, n_channel=n_channels, n_sample_thresh=n_sample_thresh)
	end = time.time()
	print('Time for Narrowband: ' + str(end-start))
	"""axarr[0].set_title('BEFORE red, AFTER blue')
	axarr[1].plot(power[1017], 'b')	
	arr[1].pcolormesh(power,cmap='coolwarm')
	#power = np.transpose(power)
	power_fold = [sum(l) for l in power]
	xarr[1].plot(np.array(power_fold).flatten(),'b')
	plt.show()
"""
	power = power/36
	idx=np.where(power>255)
	power[idx] = 255
	power = np.uint8(power)

	return power

def RFI_removal(thresh_1=3, thresh_2=2, thresh_3=1, thresh_arg=1, n_sample_thresh_2=80, n_sample_thresh_3=160, start_channel=0, n_channel=1024, n_sample_thresh=1024):
        """Removes RFI by simple threshold test
        data: 2D Data Array to be passed in with rows representing channels and columns representing time
        thresh_1: Standard Deviation is multiplied by this value to set first threshold
        thresh_2: Standard Deviation is multiplied by this value to set second threshold
        thresh_3: Standard is multiplied by this value to set third threshold
        thresh_arg: If 1, use threshold 1
                    If 2, use threshold 1 and 2
                    If 3, use threshold 1, 2 and 3
        n_sample_thresh_2: Number of samples to apply threshold 2 to around sample picked out by threshold 1
        n_sample_thresh_3: Number of samples to apply threshold 3 to around sample picked out by threshold 1
        channel_start: Start channel for RFI removal
        n_channel: Number of channels to analyze at a time (1 for narrowband, 1024 for broadband, 500 for midband)
        n_sample_thresh: Number of samples to use for calculation of median and standard deviation for threshold

        For each time sample, takes N channels (N between 1 and 1024) and computes total energy for those N channels (sums the 	       energies up) and if value exceeds threshold replace with average value for channel across time +/- noise for
	all individual channels within the N channels.
        If thresh_arg is set to 2 or 3, then for each time sample that exceeds the threshold value, compare the surrounding	       n_sample_thresh_2 and n_sample_thresh_3 samples to the theshold values set by thresh_2 and thresh_3 respectively,
	and if exceeds the theshold value, replace each channel with average value across time for the channel +/- noise

	Threshold value is calculated by taking mean and standard deviation of total power for each time sample across all
	time samples and calculating mean +/- thresh_#*SD.

	Noise is calculated by taking gaussian noise about average value of channel across all of time with standard 
	deviation being standard deviation of power for channel across time
        """
	#start = time.time()
	global power
        rows = power.shape[0]                            #number of channels
        columns = power.shape[1]                         #number of time samples
        for i in range(rows/n_sample_thresh):        #break data into n_sample_thresh chunks to analyze at a time
		power_sample = power[n_sample_thresh*i:n_sample_thresh*(i+1),start_channel:start_channel+n_channel]
        	rows = power_sample.shape[0]                            #number of channels
        	columns = power_sample.shape[1]                         #number of time samples

		#axis 0 is columns, axis 1 is rows
        	avg_power_across_time_per_channel = []
		avg_total_power_across_time = []
        	total_power_of_channels_per_time = np.sum(power_sample,axis=1)
		avg_power_across_time_per_channel = np.mean(power_sample, axis =0)
		SD_avg_power_across_time_per_channel = np.std(power_sample, axis =0) 
		
		avg_total_power_across_time = np.mean(total_power_of_channels_per_time)
		SD_total_power_across_time = np.std(total_power_of_channels_per_time)
		
		thresh1_top = avg_total_power_across_time + thresh_1*SD_total_power_across_time
		thresh1_bottom = avg_total_power_across_time - thresh_1*SD_total_power_across_time
		thresh2_top = avg_total_power_across_time + thresh_2*SD_total_power_across_time
		thresh2_bottom = avg_total_power_across_time - thresh_2*SD_total_power_across_time
		thresh3_top = avg_total_power_across_time + thresh_3*SD_total_power_across_time
		thresh3_bottom = avg_total_power_across_time - thresh_3*SD_total_power_across_time
	
		#noise = np.random.normal(0,SD_total_power_across_time/1024,n_channel)
		for j in xrange(n_sample_thresh):
			if total_power_of_channels_per_time[j] < thresh1_bottom or total_power_of_channels_per_time[j] > thresh1_top:		
				for m in xrange(n_channel):
					power_sample[j,m] = avg_power_across_time_per_channel[m] #+noise[m] #np.random.normal(avg_power_across_time_per_channel[m],SD_avg_power_across_time_per_channel[m], 1)[0]#len(power[start_channel+i:start_channel+i+1,n_sample_thresh*i+j]))

				if thresh_arg == 2 or thresh_arg == 3:
					for k in xrange(n_sample_thresh_2):
						if j-n_sample_thresh_2/2+k >= 0 and j-n_sample_thresh_2/2+k < n_sample_thresh-1:
							if total_power_of_channels_per_time[j-n_sample_thresh_2/2+k] < thresh2_bottom or total_power_of_channels_per_time[j-n_sample_thresh_2/2+k] > thresh2_top:
								for l in xrange(n_channel):
									power_sample[(j-n_sample_thresh_2/2+k),l] =  avg_power_across_time_per_channel[l]#+noise[l]#np.random.normal(avg_power_across_time_per_channel[l],SD_avg_power_across_time_per_channel[l],1)[0]

				if thresh_arg == 3:		
					for k in xrange(n_sample_thresh_3):
						if j-n_sample_thresh_3/2+k >= 0 and j-n_sample_thresh_3/2+k < n_sample_thresh-1:
							if total_power_of_channels_per_time[j-n_sample_thresh_3/2+k] < thresh3_bottom or total_power_of_channels_per_time[j-n_sample_thresh_3/2+k] > thresh3_top:
								for l in xrange(n_channel):
									power_sample[(j-n_sample_thresh_3/2+k),l] =  avg_power_across_time_per_channel[l] # +noise[l]#np.random.normal(avg_power_across_time_per_channel[l],SD_avg_power_across_time_per_channel[l],1)[0]
	power[n_sample_thresh*i:n_sample_thresh*(i+1),start_channel:start_channel+n_channel] = power_sample	


def radar_removal(power, n_spectra):
	power = power.reshape((len(power)/1536,1536))
	power = [l[bottom:top] for l in power]
	power_fold = [sum(l) for l in power]
	power = np.array(power)
	power = power.flatten()
	#remove 12 second periodic signal
	indices = [0]
	index = np.argmax(power_fold)
	i = 0
	#12016 is number of spectra between peaks of pulses
	while (index-12003) >= 0:
		indices[0] = index
		index = index-12003
	indices[0] = index
	while (indices[-1]<n_spectra):
		i+=1
		indices.append(indices[0]+12003*i)
	indices.pop()
	for i in range(len(indices)):
		indices1 = np.argmax(power_fold[(indices[i]-100):(indices[i]+100)])
		indices[i]= indices[i]-100+indices1
	f, axarr = plt.subplots(2, sharex= True)
	g, xarr = plt.subplots(2, sharex= True)
	xarr[0].plot(power)
	axarr[0].plot(power_fold)
	for i in range(61440/1024):
		print np.std(power_fold[i*1024:(i+1)*1024])/1536
		print np.mean(power_fold[i*1024:(i+1)*1024])/1536
	for i in range(len(indices)):
		bottom_index = 1024*indices[i]-100*1024
		top_index = 1024*indices[i]+100*1024
		bottom_replace = 1024*indices[i]-251*1024
		top_replace = 1024*indices[i]-51*1024
		noise = np.random.normal(np.mean(power[0:1000]),np.std(power[0:1000]),len(power[bottom_replace:top_replace]))
		power[bottom_index:top_index] = noise 
	xarr[1].plot(power, 'r')
	power = power.reshape((len(power)/1024,1024))
	power_fold = [sum(l) for l in power]
	axarr[1].plot(power_fold,'r')
	plt.show()
	return power


def RFI_removal1(data):
	"""Simple threshold RFI removal for 1000 time samples and 1024 channels...
 
	FOR BROADBAND:
	For each time sample, calculate the average power across the channels. Calculate the median and standard deviation of these values. 
	Create threshold of median +/- 3SD. For each time sample, compare the average power to this threshold value. If lies outside of this 
	threshold, replace the value of each channel in the time sample with the average value of that channel across the 1000 time samples + some noise.
 
	FOR NARROWBAND:
	For each channel, calculate the average power across time. Calculate the median and standard deviation of the values. Create threshold 
	of median +/- 3SD. For each channel, compare the average power for each time for that channel to the threshold value. If lies outside
	of threshold, replace just the sample that exceeds the threshold with the median value of the average power for that channel + some noise. 
 
	FOR 12 SECOND PULSE:
	Create threshold same way as broadband removal. For each time sample compare value to threshold. If outside threshold, replace the time 
	sample with gaussian noise about the median value of the averages of the time samples. Lower threshold to 2SD and then compare surrounding 
	50 time samples. If time sample outside threshold, replace with gaussian noise. Lower threshold to 1SD, and compare surround 100 time samples. 
	If time sample outside threshold, replace with gaussian noise.
	""" 
 
	rows = data.shape[0]
	columns = data.shape[1]
	power_fold = [sum(l) for l in data]
	power_mean = np.ones(rows*columns)*np.mean(power_fold)/1024
	mean_channel = np.mean(data, axis = 0)
	median_across_channels = np.median(mean_channel)	#calculate median of these mean values
	std_across_channels = np.std(mean_channel)		#calculate median of these mean values
	threshold_top = median_across_channels+std_across_channels*3	#set thresholds
	threshold_bottom = median_across_channels-std_across_channels*3
	for i in range(columns):
		noise = np.random.normal(median_across_channels,std_across_channels,len(data[:,i]))
		if mean_channel[i] < threshold_bottom or mean_channel[i] > threshold_top: #if mean of all channels for a given time sample is outside threshold range	
			data[:,i] = noise	#replace that time sample with noise
			#data[:,i] = 0

	mean_time = np.mean(data, axis = 1)
	mean_time = []
	for i in range(rows):
		mean_time.append(np.mean(data[i,:]))
	median_across_time = np.median(mean_time)	#calculate median of these mean values
	std_across_time = np.std(mean_time)		#calculate median of these mean values
	threshold_top = median_across_time+std_across_time*3	#set thresholds
	threshold_bottom = median_across_time-std_across_time*3
	for i in range(columns):
		noise = np.random.normal(median_across_channels,std_across_channels,len(data[:,i]))
		if mean_channel[i] < threshold_bottom or mean_channel[i] > threshold_top: #if mean of all channels for a given time sample is outside threshold range	
			data[:,i] = noise	#replace that time sample with noise
			#data[:,i] = 0

	mean_time = np.mean(data, axis = 1)
	mean_time = []
	for i in range(rows):
		mean_time.append(np.mean(data[i,:]))
	median_across_time = np.median(mean_time)	#calculate median of these mean values
	std_across_time = np.std(mean_time)		#calculate median of these mean values
	threshold_top = median_across_time+std_across_time*3	#set thresholds
	threshold_bottom = median_across_time-std_across_time*3
	for i in range(rows):
		noise = np.random.normal(median_across_time,std_across_time,len(data[i,:]))
		if mean_time[i] < threshold_bottom or mean_time[i] > threshold_top: #if mean of all channels for a given time sample is outside threshold range	
			#data[i,:] = 0
			data[i,:] = noise	#replace that time sample with noise
		data[i,:] = data[i,:]-np.mean(data[i,:])
	#data = data.flatten()-power_mean
	return data

def d():
	[l[bottom:top] for l in power]
	global power
        rows = power.shape[0]                            #number of channels
        columns = power.shape[1]                         #number of time samples
        for i in range(columns/n_sample_thresh):        #break data into n_sample_thresh chunks to analyze at a time
		power_sample = power[start_channel:start_channel+n_channel, n_sample_thresh*i:n_sample_thresh*(i+1)]
        	rows = power_sample.shape[0]                            #number of channels
        	columns = power.shape[1]                         #number of time samples
        	avg_power_across_time_per_channel = []
		avg_total_power_across_time = []
        	total_power_of_channels_per_time = []
		SD_avg_power_across_time_per_channel = []
                for j in range(n_sample_thresh):
                        total_power_of_channels_per_time.append(np.mean(power[start_channel:start_channel+n_channel,n_sample_thresh*i+j])) 
                for j in range(n_channel):
			avg_power_across_time_per_channel.append(np.mean(power[start_channel+j,n_sample_thresh*i:n_sample_thresh*(i+1)]))
			SD_avg_power_across_time_per_channel.append(np.std(power[start_channel+j,n_sample_thresh*i:n_sample_thresh*(i+1)]))
		avg_total_power_across_time = np.mean(total_power_of_channels_per_time)
		SD_total_power_across_time = np.std(total_power_of_channels_per_time)
		
		thresh1_top = avg_total_power_across_time + thresh_1*SD_total_power_across_time
		thresh1_bottom = avg_total_power_across_time - thresh_1*SD_total_power_across_time
		thresh2_top = avg_total_power_across_time + thresh_2*SD_total_power_across_time
		thresh2_bottom = avg_total_power_across_time - thresh_2*SD_total_power_across_time
		thresh3_top = avg_total_power_across_time + thresh_3*SD_total_power_across_time
		thresh3_bottom = avg_total_power_across_time - thresh_3*SD_total_power_across_time
		#print("len(avg_power_across_time_per_channel)")
		#print(len(avg_power_across_time_per_channel))
		#print("\ntotal_power_of_channels_per_time")
		#print(len(total_power_of_channels_per_time))
		#print(start_channel)
		for j in range(n_sample_thresh):
			if total_power_of_channels_per_time[j] < thresh1_bottom or total_power_of_channels_per_time[j] > thresh1_top:		
				start = time.time()
				for m in range(n_channel):
					power[start_channel+m,n_sample_thresh*i+j] = np.random.normal(avg_power_across_time_per_channel[m],SD_avg_power_across_time_per_channel[m], 1)[0]#len(power[start_channel+i:start_channel+i+1,n_sample_thresh*i+j]))
				end = time.time()
				print("first block: " + str(end-start))
				start = time.time()
				if thresh_arg == 2 or thresh_arg == 3:
					for k in range(n_sample_thresh_2):
						if j-n_sample_thresh_2/2+k >= 0 and j-n_sample_thresh_2/2+k < n_sample_thresh-1:
							if total_power_of_channels_per_time[j-n_sample_thresh_2/2+k] < thresh2_bottom or total_power_of_channels_per_time[j-n_sample_thresh_2/2+k] > thresh2_top:
								for l in range(n_channel):
									power[start_channel+l,n_sample_thresh*i+(j-n_sample_thresh_2/2+k)] =  np.random.normal(avg_power_across_time_per_channel[l],SD_avg_power_across_time_per_channel[l],1)[0]

				if thresh_arg == 3:		
					for k in range(n_sample_thresh_3):
						if j-n_sample_thresh_3/2+k >= 0 and j-n_sample_thresh_3/2+k < n_sample_thresh-1:
							if total_power_of_channels_per_time[j-n_sample_thresh_3/2+k] < thresh3_bottom or total_power_of_channels_per_time[j-n_sample_thresh_3/2+k] > thresh3_top:
								for l in range(n_channel):
									power[start_channel+l,n_sample_thresh*i+(j-n_sample_thresh_3/2+k)] = np.random.normal(avg_power_across_time_per_channel[l],SD_avg_power_across_time_per_channel[l],1)[0]

				end = time.time()
				print("second block: " + str(end-start))
def get_data_1(filename,bottom=240, top=1264):
    """ Read data in Deepthi's file format 
    
    This uses np.fromfile and numpy arrays instead of struct,
    which should make it faster than the original version
    """
    n_data = 2097152*9  # Number of data points per buffer
    n_flag = 1024*9      # Number of flags per buffer
    
    with open(filename, 'r') as fp:
        file_size = os.stat(filename).st_size
        no_cpy_buf = int(file_size/(2097152*18))
        noise = np.random.normal(0,1000,1024)
        #print no_cpy_buf
        
        raw_data = np.zeros((no_cpy_buf, n_data), dtype='<H')
        for ii in range(no_cpy_buf):    
            #copy = struct.unpack(struct_fmt, fp.read(struct_len))
            data = np.fromfile(fp, dtype='<H', count=n_data)
	    #print(len(data))
	    flag = np.fromfile(fp, dtype='<B', count=n_flag)
	    zeros = np.zeros(2048, dtype='<H')
            for y in range(len(flag)):
                if flag[y] == 1:
                    data[2048*y:2048*(y+1)] = zeros
	    raw_data[ii] = data
		#noise = data
        total_data = raw_data.flatten()
	xx      = total_data[0::4]
        yy      = total_data[1::4]

	xx = xx.reshape((len(xx)/1536,1536))
	xx = [(l[bottom:top]) for l in xx]
	xx = np.array(xx)
	xx = xx.flatten()

	yy = yy.reshape((len(yy)/1536,1536))
	yy = [(l[bottom:top]) for l in yy]
	yy = np.array(yy)
	yy = yy.flatten()
	
	power = np.sqrt((np.array(xx)/3.5)**2+np.array(yy)**2)

	return power


if __name__ == '__main__':
	#header_filename = 'HEADER.txt'
	filename = ''
	filename_out = sys.argv[1]
	#bottom = int(sys.argv[2])
	#top = int(sys.argv[3])

	os.remove('filename8.txt')
	os.system('ls /data0/data_sep-27-2017/data_2017-09-27_15-20-12_509768885 > filename8.txt')
	fp_fn=open("filename8.txt","r")
	filename_r=fp_fn.readline()
        filename=filename_r[0:-1]
        #if filename_r=='': break
        print ("open raw file: %s ..."%filename)
	#filename = raw_input("Please enter filename to read:\n")

	# Calculate the MJD time from filename
	TIME = fp_fn.readline()[5:] 
	date = TIME[0:10]
	hour = TIME[11:13]
	minute = TIME[14:16]
	second = TIME[17:19]
	nanosecond = TIME[20:]
	hour = str((int(hour)+7)%24)
	TIME = date+"T"+ hour + ":" + minute + ":" + second + "." + nanosecond
	TIME = Time(TIME, format='isot', scale='utc') 
	TIME = TIME.mjd
	#print(float(time)+240000.5)
	fp_fn.seek(0) 
	#dec_dms = (str(int(c.dec.dms[0]))+":" + str(int(c.dec.dms[1]))+ ":" +str(c.dec.dms[2]))

	# Generate a filterbank header
	fil_header = {}
	fil_header['nbits'] = 8
	fil_header['tsamp'] = float(2**8*8192/2.1/10**9)  #vacc_len*8192/2.1G .00099865*2.**5
	fil_header['nchans'] = 1024 #top-bottom
	fil_header['nifs'] = 1
	fil_header['data_type'] = 1
	fil_header['fch1'] = 1644.72656-.25634765/2.
	fil_header['foff'] = -.25634765
	fil_header['tstart'] = TIME
	fil_header['telescope_id']  = 1      # This is pretty much just made up
	fil_header['source_name']   = 'B0329+54' 
	fil_header['data_type']     = 1       # Filterbank data
	fil_header['rawdatafile']   = filename
	#fil_header['src_raj']       = Angle('03:32:59.368', unit='hour')
	#fil_header['src_raj']       = Angle('05:34:31.973', unit='hour')
	#fil_header['src_dej']	= Angle(str(dec_dms) + ' degrees')
	#fil_header['az_start']      = tel_az
	#fil_header['za_start']      = tel_alt

	# Now, serialize the header into sigproc format
	fil_header_str = generate_sigproc_header(fil_header)

	# Save to disk
	#filename_out = raw_input('Please enter file name for output (ex: filterbank_test.fil): ')
	#bottom = int(raw_input('Bottom'))
	#top = int(raw_input('Top'))

	#filename_out = 'Test_June_10_2017.fil'
	print("creating %s ..." % filename_out)

	#xx_file = open('xx_file', 'w')
	#yy_file = open('yy_file', 'w')
	with open(filename_out, 'w') as fbfile:
		fil_header_str = generate_sigproc_header(fil_header)    
		fbfile.write(fil_header_str)
		while True:
			filename_r=fp_fn.readline()
			filename=filename_r[0:-1]
			if filename_r=='':
				break
			print ("open raw file: %s ..."%filename)
			data = get_data_2(filename)#,bottom,top)
			data.tofile(fbfile)
			#power = get_data_2(filename,bottom,top)
			#xx.tofile(xx_file, sep = "\n")
			#yy.tofile(yy_file, sep = "\n ")
			print ("%s conversion done"%filename)
	fp_fn.close()
	fbfile.close()
	#fp_pos.close()
	#xx_file.close()
	#yy_file.close()

	#plot_data('xx_file','yy_file')
