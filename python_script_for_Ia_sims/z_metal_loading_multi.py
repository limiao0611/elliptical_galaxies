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

def _z_metal_flux(field,data):
    return data['SN_Colour']*data['z-velocity']


def see(i):
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()

   z_metal_loading=[]
   z=[]
   dz= 0.02
   z1=-0.5
   z2 = z1+dz
   metal_inject=10**-0.8*204.4 # 204.4 * SFR density, in unit of Msun/yr/kpc^2
   while (z2<0.5):
     bx=ds.box([0,0,z1],[0.03,0.03,z2])
     z_metal_loading.append(mean(bx['z_metal_flux'])/metal_inject)
     z.append((z1+z2)/2.*5)
     z1=z1+dz
     z2=z1+dz
   time=round(ds.current_time.in_units('Myr'),1)
   cc=['g','r','c','m','y','k','orange','crimson','pink']
   plt.plot(z,z_metal_loading,label='t='+str(time)+'Myr',lw=2,c=cc[n])
#   plt.plot(z_hot,z_hot_flux,'--',lw=2.,c=cc[n])
   global n
   n=n+1

n=0
yt.add_field('z_metal_flux',function = _z_metal_flux,units = 'Msun/kpc**2/yr')
num=[25,30,35,40,45,50,53]
fig=plt.figure()
for i in num:
    see(i)

plt.plot([-2.5,2.5],[ 0,0],':',lw=3.5,c='blue')
plt.xlim(-3,3)
plt.xlabel('z [kpc]')
plt.ylabel('Metal loading')
#plt.yscale('log')
plt.legend(loc='best', fancybox=True, framealpha=0.4)
plt.savefig('z_metal_loading_evolution'+str(num)+'.png')


