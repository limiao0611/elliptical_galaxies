import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy import *

f= open('e_t1.dat')
t,etot,eth= loadtxt(f,usecols=(0,2,4),unpack=True)

print (etot)
delta_e = []
for i in range(len(etot)-1):
    delta_e.append( (etot[i+1]-etot[i])/1e51 )

print (delta_e)
print (mean(delta_e),std(delta_e))
