import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy import *

f= open('f_unshocked_t.dat')
t,f1e2,f1e3= loadtxt(f,usecols=(0,2,3),unpack=True)
t1=array(range(1,300,5))
plt.plot(t1,exp(-t1*1.0/60),label='t_decay=60 Myr')
plt.plot(t1,exp(-t1*1.0/700),label='t_decay=t_cool')
plt.scatter(t,f1e2,label='SN_Colour<1e2')
plt.scatter(t,f1e3,label='SN_Colour<1e3')
f.close()


#plt.plot(t,etot,label='Etot')
plt.xlabel('t [Myr]')
plt.ylabel('volume fraction')
plt.yscale('log')
plt.ylim(1e-6,1)
plt.legend(loc=3)
plt.grid()
plt.savefig('f_unshocked_t.png')
plt.show()

