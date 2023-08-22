import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
#yt.enable_parallelism()
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

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

def see(i):
  fn = get_fn(i)

  ds = yt.load(fn)

  c1 = [0.,0.,0.0]
  #c1=[0.545,0.599,0.602]
  #c1 =[0.891576, 0.130623, 0.186343]
#  sp=ds.sphere(c1, (1, 'kpc'))

#  slc = yt.ProjectionPlot(ds, 'x', ['z-velocity'],center=c1,weight_field='density').annotate_grids()
#  slc.save()
#  slc = yt.ProjectionPlot(ds, 'x', ['TotalEnergy'],center=c1,weight_field='density').annotate_grids()
#  slc.save()
  slc = yt.ProjectionPlot(ds, 'x', ['density'],center=c1,weight_field='density', width=(2,'kpc'))
  slc.save()

num = range(260,340,10)
for i in num:
   see(i)

