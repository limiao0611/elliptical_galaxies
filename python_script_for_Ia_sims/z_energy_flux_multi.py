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

def _z_ek_flux(field,data):
    return (0.5*data['density']*data['z-velocity']*   (data['x-velocity']*data['x-velocity'] + data['y-velocity']*data['y-velocity'] +data['z-velocity']*data['z-velocity'] ))

def _z_eth_flux(field,data):
    return data['density']*data["GasEnergy"]*data['z-velocity']


def see(i):
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()

   z_flux_ek=[]
   z_flux_eth=[]
   z=[]
   dz= 0.02
   z1=-0.5
   z2 = z1+dz
   energy_inject=0.00349 #MW case, in unit of erg/s/cm**2
   while (z2<0.5):
     bx=ds.box([0,0,z1],[0.03,0.03,z2])
     z_flux_ek.append(mean(bx['z_ek_flux'])/energy_inject)
     z_flux_eth.append(mean(bx['z_eth_flux'])/energy_inject)
     z.append((z1+z2)/2.*5)
     z1=z1+dz
     z2=z1+dz
   time=round(ds.current_time.in_units('Myr'),1)
   cc=['g','r','c','m','y','k','orange','crimson','pink']
   plt.plot(z,z_flux_ek,label='Ek:t='+str(time)+'Myr',lw=2,c=cc[n])
   plt.plot(z,z_flux_eth,'--',lw=2,c=cc[n])
#   plt.plot(z_hot,z_hot_flux,'--',lw=2.,c=cc[n])
   global n
   n=n+1
n=0
yt.add_field('z_ek_flux',function=_z_ek_flux, units  = 'erg/s/cm**2')
yt.add_field('z_eth_flux',function=_z_eth_flux, units  = 'erg/s/cm**2')
num=[30,35,40,45,50,53]
fig=plt.figure()
for i in num:
    see(i)
#plt.plot([-2.5,2.5],[0.00349,0.00349],':',label='SN injection rate',lw=3.5,c='blue')
plt.xlim(-3,3)
plt.xlabel('z [kpc]')
plt.ylabel('energy loading')
#plt.yscale('log')
plt.legend(loc='best', fancybox=True, framealpha=0.4)
plt.savefig('z_energy_loading_evolution'+str(num)+'.png')


