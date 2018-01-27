import numpy as np
from sigpyproc.Readers import FilReader
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Open a filterbank file and plot data',\
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('files',type=str,nargs='+',
                    help='List of filterbank files you want to read')
parser.add_argument('-s','--start',type=int, default=0,
                    help='Start sample you want to plot from')
parser.add_argument('-e','--end',type=int,default=30000,
                    help='End sample you want to plot to')
parser.add_argument('-t','--time',action='store_true',default=False,
                    help='Plot a time series (average over freq axis)')
parser.add_argument('-f','--frequency',action='store_true',default=False,
                    help='Plot a spectrum averaged over the time samples provided')

args = parser.parse_args()

# Open filterbank file and 
# load data from sample s to e:

filbank = {}

for filename in args.files:
    fp = FilReader(filename)
    filbank[filename] = fp.readBlock(args.start,args.end)

if args.frequency:
    fp = FilReader(args.files[0])
    f_top = fp.header['ftop']
    f_bottom = fp.header['fbottom']
    nchans = fp.header['nchans']
    freq = np.linspace(f_top,f_bottom,num=nchans)

    fig,ax = plt.subplots(1,1)
    for k,v in filbank.items():
        ax.semilogy(freq,np.sum(v,axis=1),label=k)
    ax.legend()

elif args.time:
    fig,ax = plt.subplots(1,1)
    for k,v in filbank.items():
        ts = np.sum(v,axis=0)
        md = np.median(ts)
        ax.plot(ts-md,label=k)
    ax.legend()

else:
    fig,ax = plt.subplots(1,len(filbank.keys()))
    for i,d in enumerate(filbank.items()):
        ax[i].imshow(d[1])
        ax[i].set_title(d[0])

plt.show()
