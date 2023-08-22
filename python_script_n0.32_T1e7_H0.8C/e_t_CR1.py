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

f1=open("yt_field_list",'w')
f2=open("yt_derived_field_list",'w')
f3=open("e_t_CR.dat",'a')
num = range(900,1000)
#num = range(1349,1352)
for i  in num:
    filen = get_fn(i)
    print filen
    ds = yt.load(filen)
    ad = ds.all_data()
    print ds.current_time, ds.time_unit,ds.mass_unit,ds.length_unit
#    print >>f1, ds.field_list
#    print >>f2, ds.derived_field_list
#    print ad["gas","density"]
#    print ad["enzo","TotalEnergy"] #total energy density
#    print ad["enzo","GasEnergy"]   # thermal energy density
#    print ad["gas",'pressure']
#    print ad["enzo",'Temperature']
#    print ad["enzo",'x-velocity']
#    print ad["enzo",'CREnergyDensity']
    a = ds.current_time
    energy_density_unit = ds.mass_unit/ds.length_unit/ds.time_unit**2
    specific_energy_unit = ds.length_unit**2/ds.time_unit**2
    print '############'

    total_energy = sum( (ad["enzo","TotalEnergy"].in_cgs())* ad['gas','cell_mass'])    
    CR_energy = sum( ad["enzo",'CREnergyDensity'] *ad['gas','cell_mass']/ad['gas','density']) *energy_density_unit
    print >>f3,ds.current_time, total_energy,CR_energy
