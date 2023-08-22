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

font = {'family' : 'serif',
#        'sans-serif':'Helvetica',
#        'weight' : 'bold',
        'size'   : 17}

matplotlib.rc('font', **font)


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
    return abs(data['density']*data['z-velocity'])

def _CR_Energy_Density(field,data):
    return data['CREnergyDensity']* data.ds.mass_unit/data.ds.length_unit/data.ds.time_unit**2

def ReadCoolingCurve():
    FileIn = open("cool_rates.in_300K")
    data = loadtxt(FileIn)
    FileIn.close()
    return (data[:,0],data[:,2])

(LogT, LogCoolRate) = ReadCoolingCurve()

def _cooling_rate(field,data):
    return 10.0**interp(log10(data["temperature"]), LogT, LogCoolRate)* YTQuantity(1,"erg*cm**3/s")

m_H = YTQuantity(1.67e-24,'g')
mu = 0.6*m_H
k_b =YTQuantity(1.38e-16, 'erg/K')

def _cooling_rate_per_volume(field,data):
    return 1.76*0.4* data["cooling_rate"] *data["density"]**2/m_H/mu

def _t_cool_isob(field,data):
    return 2.5 * data["density"]/mu * k_b * data["temperature"]/data["cooling_rate_per_volume"]

def _t_cool_isoc(field,data):
    return 3./5* data["t_cool_isob"]

def _cooling_rate_per_cell(field,data):
    return data["cooling_rate_per_volume"]*data["cell_volume"]

def _SN_Colour_fraction(field,data):
    norm_factor = 30675.  # 1 mass in code unit/1msun
    return data["SN_Colour"]/data["density"]/norm_factor


def _number_density(field,data):
    return abs(data['density']/YTQuantity(1.67e-24,'g'))

def see(i):
  fn = get_fn(i)
  ds = yt.load(fn)
  fig = plt.figure()
  grid = AxesGrid(fig, (0.093,0.038,0.82,0.95),
                nrows_ncols = (2,2),
                axes_pad = 1.05,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="10%",
                cbar_pad="0%")
#  fields = ['density','temperature','pressure','z-velocity']#,'CR_Energy_Density','z_den_flux']
#  fields = ['density','temperature','pressure','SN_Colour']#,'CR_Energy_Density','z_den_flux']
  fields = ['number_density','temperature','pressure','SN_Colour_fraction']#,'CR_Energy_Density','z_den_flux']
  c1 = [0.0588,0.0588,0.0]
  p=yt.SlicePlot(ds,'x',fields,center=c1,axes_unit='kpc',fontsize=17)

  c = ['magma_r', 'plasma', 'arbre',  'YlGn', 'arbre', 'magma', 'magma', 'arbre']
  for m in range(len(fields)):
    p.set_cmap(fields[m], c[m] )

  p.set_zlim('number_density',1e-3,1.0)
  p.set_zlim('temperature',1e4,1e7)
  p.set_zlim('SN_Colour_fraction',1e-5,0.1)
  p.set_zlim('pressure',1e-12,1e-10)
  p.set_colorbar_label('SN_Colour_fraction', r'f$_{\rm{color}}$')
  p.annotate_timestamp(corner='upper_left',time_unit='Myr',text_args={'color':'white', 'size': 28})
  for j, field in enumerate(fields):
    plot = p.plots[field]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]
  p._setup_plots()
  if yt.is_root():
     plt.savefig('multi_x_slice_'+str(i)+'_2.png')

yt.add_field('z_den_flux',function=_z_den_flux, units  = 'Msun/yr/kpc**2')
yt.add_field('cooling_rate',function=_cooling_rate,units="erg*cm**3/s")
yt.add_field('cooling_rate_per_volume',function=_cooling_rate_per_volume,units="erg/cm**3/s")
yt.add_field('cooling_rate_per_cell',function=_cooling_rate_per_cell,units="erg/s")
yt.add_field('t_cool_isob',function=_t_cool_isob, units="Myr")
yt.add_field('t_cool_isoc',function=_t_cool_isoc, units="Myr")
yt.add_field('SN_Colour_fraction',function=_SN_Colour_fraction,units="")
yt.add_field('number_density',function=_number_density, units  = '1/cm**3')
num=[1,10,30,50,70,80,90]
num=[70,90]
for i in num:
   see(i)


