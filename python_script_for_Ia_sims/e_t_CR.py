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
num = range(5,1000)
#num = range(1349,1352)
for i  in num:
    filen = get_fn(i)
    print filen
    ds = yt.load(filen)
    ad = ds.all_data()
    print ds.current_time, ds.time_unit,ds.mass_unit,ds.length_unit
    print >>f1, ds.field_list
    print >>f2, ds.derived_field_list
    print ad["gas","density"]
    print ad["enzo","TotalEnergy"] #total energy density
#    print ad["enzo","GasEnergy"]   # thermal energy density
    print ad["gas",'pressure']
    print ad["enzo",'Temperature']
    print ad["enzo",'x-velocity']
    print ad["enzo",'CREnergyDensity']
    a = ds.current_time
    b= ds.time_unit
    print 'a*b=',a*b
    print 'b/second=',b/second
    print size(ad["enzo","TotalEnergy"])
    print ds.domain_dimensions[0]
    num_cells = ds.domain_dimensions[0]*ds.domain_dimensions[1]*ds.domain_dimensions[2]
    print num_cells
    c =sum(ad["enzo","TotalEnergy"][0:num_cells]) 
    print 'c=',c
    energy_density_unit = ds.mass_unit/ds.length_unit/ds.time_unit**2
    print energy_density_unit
    total_energy= energy_density_unit * c* (ds.length_unit/ds.domain_dimensions[0])**3
#    total_energy = 0. 
#    thermal_energy = 0. 
    specific_energy_unit = ds.length_unit**2/ds.time_unit**2
    print '############'
    print  ad["enzo","TotalEnergy"][0:2] * ad["gas","density"][0:2]

#    for i in range(num_cells):
    total_energy=sum(ad["enzo","TotalEnergy"][0:num_cells]*specific_energy_unit * ad["gas","density"][0:num_cells] * (ds.length_unit/ds.domain_dimensions[0])**3 )
#    thermal_energy= sum(ad["enzo","GasEnergy"][0:num_cells]*specific_energy_unit * ad["gas","density"][0:num_cells] * (ds.length_unit/ds.domain_dimensions[0])**3 )

#    d = sum(ad["enzo","GasEnergy"][0:num_cells])
#    thermal_energy= energy_density_unit * d* (ds.length_unit/ds.domain_dimensions[0])**3
    time = a*b/second/3.15e7 
    pressure_total = sum( ad["gas",'pressure'][0:num_cells])*(ds.length_unit/ds.domain_dimensions[0])**3*1.5
    f_vh=0.
    hot = ad["temperature"].in_units('K') > 2e5
    num_hot = size(ad['temperature'][hot])
    f_vh =num_hot*1.0/num_cells
    print 'f_vh',f_vh
#    kinetic_energy_total = sum (ad['gas','kinetic_energy']*ad['index','cell_volume']) #.in_units['erg'] )
#    thermal_energy = total_energy-kinetic_energy_total
#    print>>f3, time, total_energy, thermal_energy,kinetic_energy_total, pressure_total,f_vh
    total_energy = sum( (ad["enzo","TotalEnergy"].in_cgs())* ad['gas','cell_mass'])    
    CR_energy = sum( ad["enzo",'CREnergyDensity'] *ad['gas','cell_mass']/ad['gas','density']) *energy_density_unit
    print >>f3,time, total_energy,CR_energy
