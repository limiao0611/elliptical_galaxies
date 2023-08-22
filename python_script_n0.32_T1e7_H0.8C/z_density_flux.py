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
    return data['density']*data['z-velocity']


def see(i):
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()
#   print ad['density']
#   print ad['z-velocity'].in_cgs()
#   print ad['z_den_flux']
   sl = ds.slice('z',0.25)
   print sl['z_den_flux']
   print mean(sl['z_den_flux'])
   plot= yt.ProfilePlot(ad,'z','z_den_flux',weight_field='cell_volume', accumulation = False)
#   plot.set_ylim('z_den_flux',1e-24,1e-13)
   plot.set_log('z',log=False)
   plot.set_unit('z','kpc')
#   plt.plot([0,2.5],[0.015,0.015],label='SF rate')
#   plt.legend(loc=1)
   plot.save() 

#yt.add_field('z_den_flux',function=_z_den_flux, units  = 'g/s/cm**2')
yt.add_field('z_den_flux',function=_z_den_flux, units  = 'Msun/yr/kpc**2')
num = [500,900]
for i in num:
    see(i)
    
