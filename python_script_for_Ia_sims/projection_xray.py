import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
#yt.enable_parallelism()
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

font = {'family' : 'serif',
#        'sans-serif':'Helvetica',
#        'weight' : 'bold',
        'size'   : 23}

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

def see(i):
  fn = get_fn(i)

  ds = yt.load(fn)

  c1 = [0.0208,0.0208,0.0]
  #c1=[0.545,0.599,0.602]
  #c1 =[0.891576, 0.130623, 0.186343]
#  sp=ds.sphere(c1, (1, 'kpc'))

#  slc = yt.ProjectionPlot(ds, 'x', ['z-velocity'],center=c1,weight_field='density').annotate_grids()
#  slc.save()
#  slc = yt.ProjectionPlot(ds, 'x', ['TotalEnergy'],center=c1,weight_field='density').annotate_grids()
#  slc.save()
  yt.add_xray_emissivity_field(ds, 0.5, 3.5 ,with_metals=True, coonstant_metallicity=0.5)

  slc = yt.ProjectionPlot(ds, 'x', ['xray_emissivity_0.5_3.5_keV'],center=c1,weight_field=None, fontsize=25) #.annotate_grids()
  slc.set_zlim('xray_emissivity_0.5_3.5_keV', 1e-5, 1e-3)
  slc.set_cmap('xray_emissivity_0.5_3.5_keV', 'plasma')
  slc.set_colorbar_label('xray_emissivity_0.5_3.5_keV', 'Projected X-ray Emission 0.5-3.5 keV  [erg/cm$^2$/s] ')
  slc.annotate_timestamp(corner='upper_left',time_unit='Myr',text_args={'color':'white', 'size':45})
  slc.save()

num = [202,300,330]
num=[1,50,100,150,200,250,300,350,390,440]
num=[50,250,440]
for i in num:
   see(i)

