import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
yt.enable_parallelism()
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import glob

#num=range(10,1200)

ts=yt.load(glob.glob('DD*/sb_????'))
storage ={}
for sto, ds in ts.piter(storage=storage):

  fig = plt.figure()

  grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 2),
                axes_pad = 1.0,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%")

  fields = ['density','temperature','pressure','CREnergyDensity']
  c1 = [0.05,0.05,0.0]
  p=yt.SlicePlot(ds,'z',fields,center=c1)

#  p.set_unit('z-velocity','km/s')
  p.set_zlim('density',1e-28,1e-20)
  p.set_zlim('temperature',5,1e8)
  p.set_zlim('pressure',1e-15,1e-10)
#  p.set_zlim('z-velocity',-2e3,2e3)



  for j, field in enumerate(fields):
    plot = p.plots[field]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]


  p._setup_plots()
  p.annotate_timestamp(corner='uppler_left',time_unit='Myr',text_args={'color':'black'})
  plt.savefig('./z_slice/multi_z_slice_'+str(ds)+'.png')
  sto.result_id=str(ds)

