import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
yt.enable_parallelism()
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import glob
from numpy import *


def get_fn(i):
    a0=''
    if (i<10):
      a0='000'
    if (10<=i<100):
      a0='00'
    if (100<=i<999):
      a0='0'
    filen='DD'+a0+str(i)+'/sb_'+a0+str(i)
    return filen

def _z_den_flux(field,data):
    return (data['density']*data['z-velocity'])


def see(i):
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()
   hot = ad.cut_region("obj['temperature']>3e4")
   prof = yt.create_profile(ad,'z','z_den_flux',units = {'z':'kpc'}, weight_field='cell_volume') 
   z=prof.x.value
   z_flux = prof['z_den_flux'].value
   prof_hot = yt.create_profile(hot,'z','z_den_flux',units={'z':'kpc'},weight_field='cell_volume')
   z_hot = prof_hot.x.value
   z_hot_flux = prof_hot['z_den_flux'].value
   time=round(ds.current_time.in_units('Myr'),1)
   cc=['g','r','c','m','y','k','orange','crimson','pink']
   plt.plot(z,z_flux,label='t='+str(time)+'Myr',lw=2,c=cc[n])
   plt.plot(z_hot,z_hot_flux,'--',lw=2.,c=cc[n])
   global n
   n=n+1

n=0
yt.add_field('z_den_flux',function=_z_den_flux, units  = 'Msun/yr/kpc**2')
num = [40,58]
num=[80,100,130,160,190,210,250,276]
num = [130,180,230,260,295]
num=[20,25,30,35,39]
num=[80,100,120,140,160,166]
#num = [1000]
fig=plt.figure()
for i in num:
    see(i)

plt.plot([0,2.5],[10**-0.8, 10**-0.8],':',label='SF rate',lw=3.5,c='blue')
plt.xlim(1e-2,3)
plt.xlabel('z [kpc]')
plt.ylabel('z-density flux  [Msun/kpc^2/yr]')
plt.yscale('log')
plt.legend(loc=1)
    
plt.savefig('z_den_flux_evolution'+str(num)+'.png')


