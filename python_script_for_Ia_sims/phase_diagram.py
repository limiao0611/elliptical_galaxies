import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
from yt import *
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

def _metallicity(field,data):
    return data["SN_Colour"]/data["density"]/3.07e4 # normalized to Ia ejecta

def _cool_rate(field,data):
    return data["cell_mass"]*data["GasEnergy"]/data["cooling_time"]

def _cool_time_inv(field,data):
    return 1./data["cooling_time"]

def _SN_Colour_fraction(field,data):
    return data["SN_Colour"]/data["density"]


def see(i):
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()

   a= "metallicity"
   b=  "cooling_time"
   c = "cell_mass"
   d= None

#   hot = ad.cut_region("obj['temperature']>3e4")
#   midplane = ad.cut_region("obj['z']<0.1 and obj['z']>-0.1")
#   midplane = ad.cut_region("obj['z']<0.1") and ad.cut_region("obj['z']>-0.1")
#   midplane = ds.box([0,0,-0.1],[0.24,0.24,0.1])
#   halo = ad.cut_region("obj['z']>0.2") and ad.cut_region("obj['z']<-0.2")

   plot = PhasePlot(ad,a,b ,[c] ,weight_field=d)
   plot.set_unit("cooling_time","Myr")
   plot.save(a+'_'+b+'_'+c+'_all_'+str(i)+'.png')
#   plot.save(a+'_'+b+'_'+c+'_all_'+str(i)+d+'-weighted.png')

#   plot = PhasePlot(midplane,a,b,c,weight_field=d)
#   plot.save(a+'_'+b+'_'+c+'_midplane1_'+str(i)+'.png')
#   plot.save(a+'_'+b+'_'+c+'_midplane1_'+str(i)+d+'-weighted.png')

#   plot = PhasePlot(halo,a,b,c,weight_field=d)
#   plot.save(a+'_'+b+'_'+c+'halo'+str(i)+'.png')
#   plot.save(a+'_'+b+'_'+c+'halo'+str(i)+d+'-weighted.png')

yt.add_field('metallicity',function=_metallicity)
yt.add_field('cool_rate',function=_cool_rate,units="erg/s")
yt.add_field('SN_Colour_fraction',function=_SN_Colour_fraction,units="")
yt.add_field('cool_time_inv',function=_cool_time_inv,units="1/s")
num=range(1,100,20)
#num=[1,100,200]
for i in num:
    see(i)

