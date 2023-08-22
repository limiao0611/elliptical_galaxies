import yt
import matplotlib
import matplotlib.pyplot as plt
from yt.units import second, gram, parsec,centimeter, erg
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

def _velocity(field,data):
    return  (data['x-velocity']*data['x-velocity'] + data['y-velocity']*data['y-velocity'] +data['z-velocity']*data['z-velocity'] )**0.5


yt.add_field('z_ek_flux',function=_z_ek_flux, units  = 'erg/s/cm**2')
yt.add_field('velocity',function=_velocity, units  = 'cm/s')


f1=open("yt_field_list",'w')
f2=open("yt_derived_field_list",'w')
num = range(56,60)
num=[300,301]
num=[298,299]
num=range(217,218)


for i  in num:
    filen = get_fn(i)
    print (filen)
    ds = yt.load(filen)
    ad = ds.all_data()
    specific_energy_unit = ds.length_unit**2/ds.time_unit**2
    a = ds.current_time
    b= ds.time_unit
    time = a*b/second/3.15e7 
    total1 = sum(ad["enzo","TotalEnergy"]*ad['cell_mass']*specific_energy_unit).in_units('erg')
    thermal1 = sum(ad["enzo","GasEnergy"]*ad['cell_mass']).in_units('erg')
    T_min = min(ad["enzo","Temperature"])
    field='velocity'
    v_dispersion_mass = (ad.quantities.weighted_average_quantity(field, 'cell_mass')).in_units('km/s')
    v_dispersion_volume = (ad.quantities.weighted_average_quantity(field, 'cell_volume')).in_units('km/s')
    t_cool = (ad.mean("cooling_time", weight='cell_volume')).in_units('Myr')

#    print>>f3, (time,total1,thermal1)
    f3=open("e_t1.dat",'a')
    print (time,total1,thermal1,T_min,v_dispersion_mass,v_dispersion_volume,t_cool)
    print  (time,total1,thermal1,T_min,v_dispersion_mass,v_dispersion_volume,t_cool,file=f3)
    f3.close()

f1.close()
f2.close()
    
