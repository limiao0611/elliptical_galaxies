#Input: 
# dt_SN 
# Edot_PEH
# lx
#

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
yt.enable_parallelism()
from yt import YTQuantity
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
    return (0.5*data['density']*data['z-velocity']*data['z-velocity']*data['z-velocity'])

def _z_eth_flux(field,data):
    return data['density']*data["GasEnergy"]*data['z-velocity']

def _z_e_flux(field,data):
    return 0.5*data['density']*data['z-velocity']*data['z-velocity']*data['z-velocity']+ data['density']*data["GasEnergy"]*data['z-velocity'] 

def _e_dot_cool(field,data):
    return data['cell_mass']*data["GasEnergy"]/data['Cooling_Time']

def _e_tot(field,data):
    global specific_energy_unit
    return data['TotalEnergy']*data['cell_mass']*specific_energy_unit

def _e_th(field,data):
    return data['GasEnergy']*data['cell_mass']

def _e_k(field,data):
    return data['kinetic_energy']*data['cell_volume']

def see(i):
   f=open('Edot_t.dat','a')
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()
   print ds.length_unit
#   z_flux_ek=[]
#   z_flux_eth=[]
#   z=[]
   dz= 0.02
   z1=-0.5
   z2 = z1+dz
#   energy_inject=1e-4 #MW case, in unit of erg/s/cm**2
   lx =0.03
   bx=ds.box([0,0,z1],[lx,lx,z2])
   z_flux1=  - mean(bx['z_e_flux'])

   z2=0.5
   z1=z2-dz
   bx=ds.box([0,0,z1],[lx,lx,z2])
   z_flux2= mean(bx['z_e_flux'])

   area = (lx*ds.length_unit)**2
   Edot_wind = (z_flux1 + z_flux2)* area 

   Edot_rad = sum(ad['e_dot_cool'])
   Etot=sum(ad['e_tot'])
   Eth = sum(ad['e_th']) 
   Ek = sum(ad['e_k'])
   Mass = sum(ad['cell_mass'])
   print 'mass=', Mass
   Edot_PEH = Mass* YTQuantity(3.5e-25 / 1.67e-24 ,'erg/s/g')

   time=ds.current_time.in_units('Myr')
   global n
   n=n+1
   global time_last, Ek_last,Etot_last,Eth_last, Edot_SN
   del_t = time-time_last
   print 'del_t=',del_t
   Edot_tot_allgas = ((Etot-Etot_last)/del_t ).in_units('erg/s')
         
   Ek_last=Ek
   Eth_last=Eth
   Etot_last=Etot
   time_last=time
   Edot_all_minus= Edot_wind+Edot_rad + Edot_tot_allgas
   print >> f, i,time, Edot_wind, Edot_rad, Edot_tot_allgas, Edot_all_minus, Edot_SN, Edot_PEH , Etot,Eth,Ek,Mass

   f.close()

n=0

fn1=get_fn(46)
ds=yt.load(fn1)
specific_energy_unit=ds.length_unit**2/ds.time_unit**2
energy_unit = ds.mass_unit* specific_energy_unit
time_last=0.* ds.time_unit
Etot_last=0.* energy_unit
Ek_last=0.* energy_unit
Eth_last=0.* energy_unit
dt_SN=4.19e4  #yr
Edot_SN = YTQuantity(1e51/dt_SN/3.15e7,'erg/s') 


yt.add_field('z_ek_flux',function=_z_ek_flux, units  = 'erg/s/cm**2')
yt.add_field('z_eth_flux',function=_z_eth_flux, units  = 'erg/s/cm**2')
yt.add_field('z_e_flux',function=_z_e_flux, units  = 'erg/s/cm**2')
yt.add_field('e_dot_cool',function=_e_dot_cool, units  = 'erg/s')
yt.add_field('e_tot',function=_e_tot, units  = 'erg')
yt.add_field('e_th',function=_e_th, units  = 'erg')
yt.add_field('e_k',function=_e_k, units  = 'erg')

num=range(50,170,4)
for i in num:
    see(i)


