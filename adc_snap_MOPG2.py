#! /usr/bin/python
import corr,time,numpy,struct,sys,logging,pylab
from optparse import OptionParser
from numpy import int32,uint32,array,zeros,arange, log10
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.gridspec import GridSpec
from numpy.fft import fft, fftshift
from numpy import absolute

HOST = '192.168.1.51'
#HOST = 'rpi-crab.casper.pvt'
#
def exit_fail():
	print 'FAILURE DETECTED. Log entries:\n',lh.printMessages()
	try:
		fpga.stop()
	except: pass
	raise
	exit()

def exit_clean():
	try:
		fpga.stop()
	except: pass
	exit()

def myplot():
	plt.clf()
	fpga.write_int('adc_scope1_snapshot_ctrl',0)
	chl_a=fpga.snapshot_get('adc_scope0_snapshot', man_trig=True, man_valid=True)
	fpga.write_int('adc_scope1_snapshot_ctrl',1)
	chl_b=fpga.snapshot_get('adc_scope1_snapshot', man_trig=True, man_valid=True)
	#a=fpga.get_adc_snapshots('snapshot0')
	#file=open('li0604.bin','wb')
	#value_a=struct.unpack('<8192B', chl_a[0:])
	value_a=numpy.fromstring(chl_a['data'], dtype = numpy.int8)
	value_b=numpy.fromstring(chl_b['data'], dtype = numpy.int8)

	adc1 = value_a[0::2]
	adc2 = value_a[1::2]
	adc3 = value_b[0::2]
	adc4 = value_b[1::2]
 	x1 = arange(0,8191,2)
	x2 = arange(1,8192,2)	
	
	adc0_fft = absolute(fftshift(fft(value_a)))
	adc1_fft = absolute(fftshift(fft(value_b)))
	x3 = arange(-len(adc0_fft)/2,len(adc0_fft)/2,1)

	#ax1 = plt.subplot2grid((3,2), (0,0))
	#ax1 = plt.subplot(gs1[0, 0])
	ax1 = plt.subplot(gs1[0])
	#plt.subplot(3,2,1)
	plt.title('ADC0')
	plt.ylabel('Level')
	plt.xlim(0,1024) 
	plt.xlabel('ADC0 Zoom in 8x')
	plt.plot(value_a,'g') 
	plt.plot(x1,adc1,'b') 
	plt.plot(x2,adc2,'r')


	#ax2 = plt.subplot2grid((3,2), (0,1))
	ax2 = plt.subplot(gs1[0, 1])	
	plt.title('ADC1')
	#plt.ylabel('Level')
	plt.xlim(0,400) 
	plt.xlabel('ADC0 Zoom in 20x')
	plt.plot(value_b,'g')
	plt.plot(x1,adc3,'b') 
	plt.plot(x2,adc4,'r')

	"""
	#ax3 = plt.subplot2grid((3,2), (1,0))
	ax3 = plt.subplot(gs1[1, 0])
	plt.xlabel('ADC0 Zoom in 20x')
	plt.ylabel('Level')
	plt.xlim(0,400)
	plt.plot(value_a,'g') 
	plt.plot(x1,adc1,'b'+'o') 
	plt.plot(x2,adc2,'r'+'*')


	#ax4 = plt.subplot2grid((3,2), (1,1))
	ax4 = plt.subplot(gs1[1, 1])
	plt.xlabel('ADC1 Zoom in 20x')
	plt.xlim(0,400)
	plt.plot(value_b,'g')
	plt.plot(x1,adc3,'b'+'o') 
	plt.plot(x2,adc4,'r'+'*')

	#ax3 = plt.subplot2grid((3,2), (2,0))
	ax5 = plt.subplot(gs1[2, 0])
	plt.xlabel('ADC0 Zoom in 80x')
	plt.ylabel('Level')
	plt.xlim(0,100)
	ADC0,=plt.plot(value_a,'g',label='ADC0') 
	coreA,=plt.plot(x1,adc1,'b'+'o',label='coreA') 
	coreB,=plt.plot(x2,adc2,'r'+'*',label='coreB')
	plt.legend([ADC0,coreA,coreB], ['ADC0', 'coreA','coreB'])


	#ax4 = plt.subplot2grid((3,2), (2,1))
	ax6 = plt.subplot(gs1[2, 1])
	plt.xlabel('ADC1 Zoom in 80x')
	plt.xlim(0,100)
	ADC1,=plt.plot(value_b,'g',label='ADC1')
	coreC,=plt.plot(x1,adc3,'b'+'o',label='coreC') 
	coreD,=plt.plot(x2,adc4,'r'+'*',label='coreD')
	plt.legend([ADC1,coreC,coreD], ['ADC1', 'coreC','coreD'])
	"""	
	ax6 = plt.subplot(gs1[1,0])
	plt.xlabel('ADC0 Voltage Histogram')
	plt.hist(adc1,256, color = 'r')
	plt.hist(adc2,256)

	ax7 = plt.subplot(gs1[1,1])
	plt.xlabel('ADC1 Voltage Histogram')
	plt.hist(adc3,256)
	plt.hist(adc4,256)
	
	ax8 = plt.subplot(gs1[2,0])
	plt.xlabel('ADC0 FFT')
	plt.xlim(0,4096)
	plt.ylim(0,5)
	plt.plot(x3,log10(adc0_fft))

	ax9 = plt.subplot(gs1[2,1])
	plt.xlabel('ADC1 FFT')
	plt.xlim(0,4096)
	plt.ylim(0,5)
	plt.plot(x3,log10(adc1_fft))
	fig.canvas.draw()
	fig.canvas.manager.window.after(100, myplot)

if __name__ == '__main__':

	print 'Connecting to board:', HOST
	fpga = corr.katcp_wrapper.FpgaClient(HOST)
	time.sleep(0.1)

	if fpga.is_connected():
		print 'DONE'
	else:
		print 'ERROR: Failed to connect to server %s!' %(HOST)
		sys.exit(0);

	fig = plt.figure(figsize=(12,6))
	'''ax1 = plt.subplot2grid((3,2), (0,0))
	ax2 = plt.subplot2grid((3,2), (0,1))
	ax3 = plt.subplot2grid((3,2), (1,0))
	ax4 = plt.subplot2grid((3,2), (1,1))
	ax5 = plt.subplot2grid((3,2), (2,0))
	ax6 = plt.subplot2grid((3,2), (2,1))'''
	plt.suptitle("ADC Snapshot")

	gs1 = GridSpec(3,2)
	#gs1.update(left=0.04, right=0.05, hspace=0.25)
	gs1.update(hspace=0.25)

	#ax2 = plt.subplot(gs1[0, 1])
	#ax3 = plt.subplot(gs1[1, 0])
	#ax4 = plt.subplot(gs1[1, 1])
	#ax5 = plt.subplot(gs1[-1, 1])
	#ax6 = plt.subplot(gs1[-1, -1])
	

	fig.canvas.manager.window.after(200, myplot)
	plt.show()
