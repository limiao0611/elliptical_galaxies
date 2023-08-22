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
#   cr =  ad.cut_region("obj['z']<0.")
#   cr = ad.cut_region("obj['z-velocity']<0.") and ad_cutregion("obj['z']>-0.")
   x='Density'
   x='Temperature'
   time=ds.current_time.in_units('Myr')
   plot = yt.ProfilePlot(ad, x ,["cell_volume"],weight_field=None,x_log=True,y_log=None,accumulation=None,fractional=True,label=str(time))
#   plot.set_unit(x,'g/cm**3')
#   plot.set_ylim('Cell_volume',1e-3,1)

#   plot.set_xlim(-1e3,1e3)

#   plot.annotate_timestamp(corner='upper_left',time_unit='Myr',text_args={'color':'black'})
#   plot.set_title("time")
   plot.save()

num=[1]
num = range(50,400,50)
for i in num:
    see(i)

