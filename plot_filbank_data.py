#! /usr/bin/env python

import numpy as np
from sigpyproc.Readers import FilReader
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import argparse

parser = argparse.ArgumentParser(description='Open a filterbank file and plot data',\
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('files',type=str,nargs='+',
                    help='List of filterbank files you want to read')
parser.add_argument('-s','--start',type=int, default=0,
                    help='Start sample you want to plot from')
parser.add_argument('-nsam','--samples',type=int,default=1000,
                    help='Number of samples you want to plot')
parser.add_argument('-t','--time',action='store_true',default=False,
                    help='Plot a time series (average over freq axis)')
parser.add_argument('-f','--frequency',action='store_true',default=False,
                    help='Plot a spectrum averaged over the time samples provided')
parser.add_argument('--save', type=str, default=None,
                    help='Save the plot with the provided filename as a pdf file')
args = parser.parse_args()

# Open filterbank file and 
# load data from sample s to e:

filbank = {}

for filename in args.files:
    fp = FilReader(filename)
    filbank[filename] = fp.readBlock(args.start,args.samples)

fp = FilReader(args.files[0])
f_top = fp.header['ftop']
f_bottom = fp.header['fbottom']
nchans = fp.header['nchans']
freq = np.linspace(f_top,f_bottom,num=nchans)

if args.frequency:
    fig,ax = plt.subplots(1,1)
    for k,v in filbank.items():
	stats = []
	stats.append("%0.3e"%((np.mean(np.sum(v, axis = 1)[300:1200]))))
	stats.append("%0.3e"%((np.median(np.sum(v, axis = 1)[300:1200]))))
	stats.append("%0.3e"%((np.std(np.sum(v, axis = 1)[300:1200]))))
        ax.semilogy(freq,np.sum(v,axis=1),label=k)
	#plt.figtext(1, 0.8, 'Mean: ' + str(stats[0]) + '\nMedian: ' + str(stats[1]) + '\nSTD: ' + str(stats[2]))
	print(k)
	print('Mean: ' + str(stats[0]) + '\nMedian: ' + str(stats[1]) + '\nSTD: ' + str(stats[2])+ '\n\n')
    ax.set_ylabel('Amplitude',fontsize=22)
    ax.set_xlabel('Frequency [MHz]',fontsize=22)
    ax.legend()

elif args.time:
    fig,ax = plt.subplots(1,1)
    for k,v in filbank.items():
        ts = np.sum(v,axis=0)
        md = np.median(ts)
        ax.plot(ts-md,label=k, alpha = 0.5)
    ax.set_ylabel('Amplitude',fontsize=22)
    ax.set_xlabel('Time',fontsize=22)
    ax.legend()

else:
    for i,d in enumerate(filbank.items()):
        plt.figure(i) 
        idx = d[1].nonzero()[0]
        data = d[1][np.nonzero(d[1])].reshape(-1,args.samples)
        plt.imshow(data,interpolation='none',aspect='auto',\
                   extent=(0,len(data[0]),freq[idx.min()],\
                   freq[idx.max()]), cmap=plt.cm.brg_r,\
                   norm=LogNorm(vmin=data.min(),vmax=data.max()))
        plt.title(d[0])
        plt.colorbar()

if args.save:
    fig.savefig(args.save+'.pdf')

plt.show()
