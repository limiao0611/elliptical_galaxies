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

def _z_den_flux(field,data):
    return (data['density']*data['z-velocity'])


def see(i):
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()
#   hot = ad.cut_region("obj['temperature']>3e4")
#   cr =  ad.cut_region("obj['z']>-0.0001") and ad.cut_region("obj['z']<0.0001")
   cr =  ds.box([0,0,-0.0004],[0.05,0.05,0.])
   print cr['temperature']
   print len(cr['temperature'])
   print max(cr['temperature'])   
   print min(cr['temperature'])   
n=0
num=[1]
#num = [1000]
for i in num:
    see(i)

