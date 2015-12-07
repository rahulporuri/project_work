######################################################################
# Python code to produce plots
######################################################################
# Importing libraries
import matplotlib
from numpy import *
from pylab import *
import pylab
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from matplotlib import rc, rcParams, pyplot
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show
from matplotlib.ticker import MaxNLocator
#from pyx import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
######################################################################
figure(1,figsize=(8,6))
######################################################################
rc('font', family='TimesRoman', weight = 'extra bold', size = 18.0)
rc('text', usetex=True)
rc('axes', linewidth = 2, labelsize = 'large')  
rc('xtick', labelsize= 'medium')
rcParams['xtick.major.size'] = 8.0 
rcParams['xtick.minor.size'] = 4.0
rcParams['xtick.major.pad'] = 8.0 
rcParams['xtick.minor.pad'] = 8.0
rc('ytick', labelsize= 'medium')  
rcParams['ytick.major.size'] = 8.0 
rcParams['ytick.minor.size'] = 0.0
rcParams['ytick.major.pad'] = 8.0 
rcParams['ytick.minor.pad'] = 8.0
rc('lines', linewidth = 1.5, markeredgewidth=1.5)
rc('savefig', dpi=300)
######################################################################
# Fixing space around the plot
pylab.axes([0.125,0.125,0.825,0.825])

#ax.xaxis.set_major_locator( xmajorLocator )
#ax.xaxis.set_minor_locator( xminorLocator )
#ax.yaxis.set_major_locator( ymajorLocator )
#ax.yaxis.set_minor_locator( yminorLocator )
######################################################################
#Reading data from files and plotting them
#read from textfile to x and y
#data = genfromtxt('amm-hnla_131.txt')
data = genfromtxt('final_data')
#data = loadtxt('amm-hnla_131.txt')

k1 = data[:,0]
k2 = data[:,1]
k3 = data[:,2]

x = k3/k1
y = k2/k1

z = data[:,10]

#x = data[:,0]
#y = data[:,1]
#z = data[:,2]

#print z

xi = np.linspace(0,1,400)
yi = np.linspace(0.5,1,400)
zi = griddata(x,y,z,xi,yi,interp='linear')
f = pyplot.figure()
ax = f.gca()
#CS = plt.contour(xi,yi,zi)#,15,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,zi,50,cmap=plt.cm.jet)#,vmax=abs(zi).max(), vmin=-abs(zi).max())
#plt.hexbin(x,y,z) 
#plt.contour(xi,yi,zi)
#plt.figure()
#im = plt.imshow(z, interpolation='bilinear', origin='lower',cmap=cm.gray, extent=(-3,3,-2,2))
#levels = np.arange(-1.2, 1.6, 0.2)
#CS = plt.contour(z, levels,origin='lower',linewidths=2,extent=(-3,3,-2,2))
plt.axis([0, 1, 0.5, 1])
pylab.xlabel(r'$k_3/k_1$')
pylab.ylabel(r'$k_2/k_1$')
#mn=int(z.min())                  # colorbar min value
#mx=int(z.max())                  # colorbar max value
#md=(mx-mn)/2                     # colorbar midpoint value
divider = make_axes_locatable(ax)
cax = divider.append_axes("right",size="3%", pad=0.25)
m1=round(zi.min(),2)
m2=round(zi.max(),2)
md=round(m1+(m2-m1)/2,2)
print m1,m2,md
cbar=plt.colorbar(ticks=[m1,md,m2], cax=cax)
#cbar.set_ticklabels([m1,md,m2])
cbar.set_ticklabels([md,m2,1]) 
#cbar.set_position(right,0.05,1)
ax.set_aspect(1)
#cb = plt.colorbar()
#plt.axis().set_aspect(0.5, adjustable='box')
#Make a contour plot of hnl against k31 and k21
######################################################################
#pylab.savefig('amm-hnla-b1e-2.eps',bbox_inches='tight')
pylab.show()
######################################################################
