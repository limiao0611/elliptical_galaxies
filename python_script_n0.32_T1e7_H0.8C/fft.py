import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
from yt import *
#yt.enable_parallelism()
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy import *
from numpy.fft import *

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

def _number_density(field,data):
    return abs(data['density']/YTQuantity(1.67e-24,'g'))



def see(i):
  fn = get_fn(i)

  ds = yt.load(fn)
  ad = ds.all_data()

  c1 = [0.5,0.5,0.5]

  print (ad['density'].size)
  kk = fftn(ad['density'])
  print (ad['density'])
  print (kk)


  max_level = ds.index.max_level
  print('max_level=',max_level)
  dims1 = ds.domain_dimensions
  nx, ny, nz= dims1
  print(nx,ny,nz)
  low = ds.domain_left_edge
  irho='density'
  iu='x-velocity'
  cube = ds.covering_grid(max_level, left_edge=low,
                            dims=dims1,
                            fields=[irho, iu])

  rho = cube[irho].d
  ux=cube[iu].d
  print (rho.size)
  print (rho.shape)

  ru=fftn(ux)
  print ("ux_k=",ru)
  
  print (ux[2][5][10])
  print ('-----------------------')
  print (ifftn(ru)[2][5][10])

yt.add_field('number_density',function=_number_density, units  = '1/cm**3')

num = [202,300,330]
num=range(600,601)
for i in num:
   see(i)

