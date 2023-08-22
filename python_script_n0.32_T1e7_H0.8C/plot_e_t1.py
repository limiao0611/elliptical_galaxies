import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy import *

f= open('e_t1.dat')
t,etot,eth= loadtxt(f,usecols=(0,2,4),unpack=True)
plt.plot(t,etot,label='Etot')
plt.show()
