import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
yt.enable_parallelism()
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import glob

def _z_den_flux(field,data):
    return abs(data['density']*data['z-velocity'])
yt.add_field('z_den_flux',function=_z_den_flux, units  = 'Msun/yr/kpc**2')


ts=yt.load(glob.glob('DD*/sb_????'))
storage ={}
for sto, ds in ts.piter(storage=storage):

  fig = plt.figure()

  grid = AxesGrid(fig, (0.075,0.075,0.70,0.90), 
                nrows_ncols = (3,2),
                axes_pad = 1,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="20%",
                cbar_pad="0%")

# fields = ['density','temperature','pressure','CREnergyDensity']
  fields = ['density','temperature','pressure','z-velocity','CREnergyDensity','z_den_flux']

  c1 = [0.05,0.05,0.0]
  p=yt.SlicePlot(ds,'x',fields,center=c1,axes_unit='kpc',fontsize=12)

#  p.set_unit('z-velocity','km/s')
  p.set_zlim('density',1e-28,1e-20)
  p.set_zlim('temperature',5,1e8)
  p.set_zlim('pressure',1e-15,1e-10)
  p.set_zlim('CREnergyDensity',1e-9,1e-2)
  p.set_zlim('z_den_flux',1e-3,1.)
  p.set_zlim('z-velocity',-2e3,2e3)
  p.annotate_timestamp(corner='uppler_left',time_unit='Myr',text_args={'color':'black'})

#  p.zoom(5)

  for j, field in enumerate(fields):
    plot = p.plots[field]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]

  p._setup_plots()
  p.annotate_timestamp(corner='uppler_left',time_unit='Myr',text_args={'color':'black'})
#  plt.title('t='+str(ds.current_time*ds.time_unit))
  plt.savefig('./x_slice/multi_x_slice_'+str(ds)+'.png')
  sto.result_id=str(ds)

