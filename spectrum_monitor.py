#! /usr/bin/python
import socket,pylab,matplotlib,math,corr,array
import struct
import time

paa=4096*[0.0]
pbb=4096*[0.0]
ab_real=4096*[0.0]
ab_img=4096*[0.0]
xx = 4096*[0]
yy = 4096*[0]
xy_real = 4096*[0]
xy_img =  4096*[0]
acclen=0
gain=0
xy_img_tmp=512*[0]
xy_real_tmp=512*[0]
yy_tmp=512*[0]
xx_tmp=512*[0]
frame=0
def get_data():
	global acclen,gain,frame
	gain=snap.read_int('vacc_shift')

	acclen=snap.read_int('vacc_acc_len')/256

	for i in range(0,8):
		print "grab number",i,"packet done!"
		data, addr = sock.recvfrom(4096+8)
		header = struct.unpack('<Q', data[0:8])[0]
		payload_len = len(data) - 8 # subtract 8 bytes of header
		for k in range(0,512):
			xx_tmp[k]=struct.unpack('<H', data[(8+8*(k+0)):(8+8*(k+0) + 2)])	
			yy_tmp[k]=struct.unpack('<H', data[(8+8*k+2):(2+8+8*k+2)])
			xy_real_tmp[k]=struct.unpack('<h', data[(8+8*k+4):(2+8+8*k+4)])
			xy_img_tmp[k]=struct.unpack('<h', data[(8+8*k+6):(2+8+8*k+6)])
		#data1 = struct.unpack('<2048h', data[8:])

		if ((header >> 9) % 1  == 0):
			print time.time(), 
			print 'received %d bytes' % payload_len,
			print 'from', addr,	  
			print 'HEADER: %d, TIMESTAMP: %d, CHANNEL %4d, SEQ: %d' % (header, header >> 12, header & 0xfff, header >> 9)
		SEQ = header >> 9
		CHANNEL = header & 0xfff
		#print 'data 1', data1[:]


		#xx_tmp = struct.unpack('<512H',data[3::4]))

		'''xy_img_tmp = data1[0::4]
		xy_real_tmp = data1[1::4]
		yy_tmp = data1[2::4]
		xx_tmp = data1[3::4]'''

		print int(CHANNEL)
		for j in range(0,512):
			xx[j+int(CHANNEL)] = xx_tmp[j]
			yy[j+int(CHANNEL)] = yy_tmp[j]
			xy_real[j+int(CHANNEL)] = xy_real_tmp[j]
			xy_img[j+int(CHANNEL)] = xy_img_tmp[j]		

		file1.write(str(yy_tmp))
		frame=frame+1
	print "grab number",frame,"frame done!"
	return SEQ,xx,yy,xy_real,xy_img


def plot_spectrum():
	pp=2
	p2=1 #subplot parameter

	matplotlib.pyplot.clf()
	seq,paa,pbb,ab_real,ab_img = get_data()
	#print seq
	#print acclen
	#print gain
	'''for i in range(0,4096):
      		paa[i]=paa[i]/acclen
      		paa[i]=12.7-6*log_g0+10*math.log10(paa[i]+0.0000001) 
		pbb[i]=pbb[i]/acclen*8
      		pbb[i]=12.7-6*log_g0+10*math.log10(pbb[i]+0.0000001) '''

	#print paa
	pylab.subplot(411)
	#pylab.title('SEQ is '+str(seq),bbox=dict(facecolor='red', alpha=0.5))
	pylab.title('SEQ is '+str(seq))
    	#pylab.title('xx')	
	pylab.semilogy(paa,color="g")
    	pylab.xlim(0,4096)
	pylab.ylabel('xx')
	#pylab.ylabel('Power(dBm)')

	#print pbb
	pylab.subplot(412)
    	#pylab.title('yy')	
	pylab.semilogy(pbb,color="b")
    	pylab.xlim(0,4096)
    	#pylab.ylim(-120,0)
    	pylab.ylabel('yy')
	#pylab.ylabel('Power(dBm)')

	#print pab_real
	pylab.subplot(413)	
	pylab.semilogy(ab_real,color="y")
    	pylab.xlim(0,4096)
	pylab.ylabel('xy_real')
	#pylab.ylabel('Power(dBm)')

	#print pab_img
	pylab.subplot(414)	
	pylab.semilogy(ab_img,color="r")
    	pylab.xlim(0,4096)
	pylab.ylabel('xy_img')
	#pylab.ylabel('Power(dBm)')

	fig.canvas.draw()	
	fig.canvas.manager.window.after(200,plot_spectrum)     
	return True

if __name__ == '__main__':

	snap=corr.katcp_wrapper.FpgaClient('192.168.1.51')
	time.sleep(0.1)
	if (snap.is_connected()==True):
		print "snap connected, done!"

	IP = "10.10.12.35" #bind on all IP addresses
	PORT = 10000

	sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
	sock.bind((IP, PORT))
	file1 = open("yy.txt","w")
	if PORT != -1:
		print "10GbE port connect done!"
		fig = matplotlib.pyplot.figure()
		#gobject.timeout_add(1000,plot_spectrum)
		fig.canvas.manager.window.after(200,plot_spectrum)
		matplotlib.pyplot.show()
	file1.close()
