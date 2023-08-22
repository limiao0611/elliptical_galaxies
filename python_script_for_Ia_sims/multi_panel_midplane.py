import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
#yt.enable_parallelism()
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def _CR_Energy_Density(field,data):
    return data['CREnergyDensity']* data.ds.mass_unit/data.ds.length_unit/data.ds.time_unit**2


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
  a0=''
  if (i<10):
      a0='000'
  if (10<=i<100):
      a0='00'
  if (100<=i<999):
      a0='0'

  fn = get_fn(i)

  ds = yt.load(fn)

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

  fields = ['density','temperature','pressure','SN_Colour']
  c1 = [0.0588,0.0588, 0.001 ]
  p=yt.SlicePlot(ds,'z',fields,center=c1)

#  p.set_unit('z-velocity','km/s')
  p.set_zlim('density',1e-27,1e-25)
  p.set_zlim('temperature',1e6,1e8)
  p.set_zlim('pressure',1e-11,1e-9)
#  p.set_zlim('CR_Energy_Density',1e-15,1e-10)
#  p.set_zlim('z-velocity',-2e3,2e3)



  for j, field in enumerate(fields):
    plot = p.plots[field]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]


  p._setup_plots()
  p.annotate_timestamp(corner='upper_left',time_unit='Myr',text_args={'color':'black'})
  if yt.is_root():
    plt.savefig('multi_z_slice_'+str(ds)+'_z'+str(c1[2])+'.png')
  #c1=[0.545,0.599,0.602]
  #c1 =[0.891576, 0.130623, 0.186343]
#  sp=ds.sphere(c1, (1, 'kpc'))

#  slc = yt.SlicePlot(ds, 'y', ['density','temperature','pressure','z-velocity'],center=c1).save()
#  slc = yt.SlicePlot(ds, 'z', ['density','temperature','pressure','z-velocity'],center=c1).save()


#  c2 = [0.087,0.087,0.17]

#  c2= [0.05,0.05,0.0]
#  slc = yt.SlicePlot(ds, 'z', ['density','temperature','pressure','z-velocity'],center=c2).save()

#yt.add_field('CR_Energy_Density',function=_CR_Energy_Density, units  = 'erg/cm**3')


num =[0,1,10,100,200,261]
num = [0,1,10,100,161]
num = [400,547]
num = [575]
num = [0,1,10,36]
num = [0,1,10,30]
num=[380,500,625]
num=[1050,1100]
#num = [1,2,3,4,5]
#num = [181]
num=range(10,1200)
num=[1,10,100,300]
num=[500,700,869]
num=[120,125,130,135,140,145,150,155,160]
num=[1,10,20,40,80,100,130,160]
num=[50,70,89]
num=[100,140]
num=[1,10,20,30,39]
num=[202,300,400,500,600,700]
num=[10,50,100,150]
num=[10,50,100,150,202,300,400,500,600,700,800]
num=[900,1000,1100]
num=range(1200,1800,100)
num=range(2100,2800,100)
num=range(10,230,10)
for i in num:
   see(i)


