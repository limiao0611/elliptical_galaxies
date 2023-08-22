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


num = range(56,60)
num=[300,301]
num=[298,299]
num=range(6,10,1)


for i  in num:
    filen = get_fn(i)
    print (filen)
    ds = yt.load(filen)
    ad = ds.all_data()
    time= (ds.current_time).in_units('Myr')

    unshock1 =  ad["SN_Colour"]<1e2
    unshock2 =  ad["SN_Colour"]<1e3
    num_cells = ds.domain_dimensions[0]*ds.domain_dimensions[1]*ds.domain_dimensions[2]
    num_1 = size(ad['SN_Colour'][unshock1])
    num_2 = size(ad['SN_Colour'][unshock2])
    f_1 =num_1*1.0/num_cells
    f_2 =num_2*1.0/num_cells


    f3=open("f_unshocked_t.dat",'a')
    print (time,f_1,f_2)
    print  (time,f_1,f_2,file=f3)
    f3.close()

    
