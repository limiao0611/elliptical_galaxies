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
   cr =  ad.cut_region("obj['z']<0.")
#   cr = ad.cut_region("obj['z-velocity']<0.") and ad_cutregion("obj['z']>-0.")
   plot = yt.ProfilePlot(cr,"z-velocity",["cell_mass"],weight_field=None,x_log=False)
   plot.set_unit('z-velocity','km/s')
#   plot.set_xlim(-1e3,1e3)
   plot.save()

num=[1]
num = [91]
for i in num:
    see(i)

