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
   warm = ad.cut_region("obj['temperature']<3e4")
   cold =ad.cut_region("obj['temperature']<5e2")   
   prof_hot = yt.create_profile(hot,'z','z-velocity',units={'z':'kpc', 'z-velocity':'km/s'},weight_field='cell_mass')
   z_hot = prof_hot.x.value
   z_hot_vel = prof_hot['z-velocity'].value
   prof_warm = yt.create_profile(warm,'z','z-velocity',units={'z':'kpc', 'z-velocity':'km/s'},weight_field='cell_mass')
   z_warm = prof_warm.x.value
   z_warm_vel = prof_warm['z-velocity'].value
   prof_cold = yt.create_profile(cold,'z','z-velocity',units={'z':'kpc', 'z-velocity':'km/s'},weight_field='cell_mass')
   z_cold = prof_cold.x.value
   z_cold_vel = prof_cold['z-velocity'].value
   time=round(ds.current_time.in_units('Myr'),1)
   cc=['g','r','c','m','y','k','orange','crimson','pink']
   plt.plot(z_hot,z_hot_vel,'-',lw=2.,c=cc[n], label='t='+str(time)+'Myr')
   plt.plot(z_warm,z_warm_vel,'--',lw=2.,c=cc[n])
   plt.plot(z_cold,z_cold_vel,':',lw=2.,c=cc[n])
   global n
   n=n+1

n=0
num = [200,300,350,400,460,500,527]
fig=plt.figure()
for i in num:
    see(i)

#plt.plot([0,2.5],[0.015,0.015],':',label='SF rate',lw=3.5,c='blue')
#plt.xlim(0,3)
plt.xlabel('z [kpc]')
plt.ylabel('z-velocity  [km/s]')
plt.ylim(10.,1e3)
plt.yscale('log')
plt.legend(loc=2)
    
plt.savefig('z-velocity_z_evolution_mass_weighted'+str(num)+'.png')


