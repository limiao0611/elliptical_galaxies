from numpy import *
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import romberg


l_unit = 5e3 #length unit: 5kpc
Ia_frac = 0.15
II_low_frac = 0.6
h_Ia = 325./l_unit#pc
h_II_low = 120./l_unit  #pc
h_II_high = 360./l_unit  #pc


def Ia(z):
    return Ia_frac * exp(-abs(z)/h_Ia) /h_Ia *0.5

def II_low(z):
    return (1-Ia_frac)*II_low_frac / ( (2*3.14)**0.5 *  h_II_low) * exp( -z**2/(2*h_II_low**2) )


def II_high(z):
    return (1-Ia_frac)*(1-II_low_frac) / ( (2*3.14)**0.5 *  h_II_high) * exp( -z**2/(2*h_II_high**2) )

fn = ('SN_pos.dat')
f=open(fn,'r')
ypos,zpos = loadtxt(f,usecols=(2,3),unpack=True)

z = arange(-1,1,0.01)
plt.plot(z, Ia(z),label='Ia')
plt.plot(z, II_low(z),label='II_low')
plt.plot(z, II_high(z),label='II_high')
plt.plot(z, Ia(z)+II_low(z)+II_high(z),label='All')
#plt.hist(zpos,bins=50,normed=True,label='from sim')

plt.ylim(0.01,15)
plt.yscale('log')
plt.legend(loc=2)
plt.savefig('SN_dist_ana_only.png')

print 'cumulative:Ia=', quad(Ia,-1,1)
print 'cumulative:II_low=', quad(II_low,-1,1)
print 'cumulative:II_high=', quad(II_high,-1,1)
