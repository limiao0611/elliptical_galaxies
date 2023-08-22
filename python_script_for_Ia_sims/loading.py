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
    return abs(0.5*data['density']*data['z-velocity']*   (data['x-velocity']*data['x-velocity'] + data['y-velocity']*data['y-velocity'] +data['z-velocity']*data['z-velocity'] ))

def _z_eth_flux(field,data):
    return abs(data['density']*data["GasEnergy"]*data['z-velocity'])

def _z_etot_flux(field,data):
    return abs( data['density']*data["GasEnergy"]*data['z-velocity']) + abs(0.5*data['density']*data['z-velocity']*   (data['x-velocity']*data['x-velocity'] + data['y-velocity']*data['y-velocity'] +data['z-velocity']*data['z-velocity'] ))

def _z_metal_flux(field,data):
    return abs(data['SN_Colour']*data['z-velocity'])

def _z_den_flux(field,data):
    return abs(data['density']*data['z-velocity'])


def see(i):
   f=open('loading.dat','a')
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()
   halo = ad.cut_region("obj['z']>0.2") and ad.cut_region("obj['z']<-0.2")


   metal_inject=10**-0.8*204.4 # 204.4 * SFR density, in unit of Msun/yr/kpc^2
   energy_inject=0.00349 #case, in unit of erg/s/cm**2
   SF_density = 10**-0.8 

   global z_ml,z_el,z_mel
   
   ml=mean(halo["z_den_flux"])/SF_density
   el=mean(halo["z_etot_flux"])/energy_inject
   mel= mean(halo["z_metal_flux"]/metal_inject)
   z_ml.append(ml)
   z_el.append(el)
   z_mel.append(mel)

   time  = round(ds.current_time.in_units('Myr'),2 )
   print >>f, i, time, ml,el,mel
   print z_ml, mean(z_ml),std(z_ml)
   print z_el, mean(z_el),std(z_el)
   print z_mel, mean(z_mel), std(z_mel)
   f.close()

n=0
yt.add_field('z_metal_flux',function = _z_metal_flux,units = 'Msun/kpc**2/yr')
yt.add_field('z_ek_flux',function=_z_ek_flux, units  = 'erg/s/cm**2')
yt.add_field('z_eth_flux',function=_z_eth_flux, units  = 'erg/s/cm**2')
yt.add_field('z_etot_flux',function=_z_etot_flux, units  = 'erg/s/cm**2')
yt.add_field('z_den_flux',function=_z_den_flux, units  = 'Msun/yr/kpc**2')

global f
z_ml=[] #mass loading
z_el=[]  # energy loading
z_mel=[]  # metal loading

num=range(50,170,4)
for i in num:
    see(i)
f=open('loading.dat','a')
print >>f, mean(z_ml), mean(z_el), mean(z_mel)
print >>f, std(z_ml),std(z_el), std(z_mel)
