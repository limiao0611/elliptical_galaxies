import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
from yt import *
yt.enable_parallelism()
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

def _v_mag_compress(field,data):
  nx = ny = nz = len(data['density'])
  if (nx > 16.1):
    nx, ny, nz = data.ds.domain_dimensions

  L = ((data.ds.domain_right_edge - data.ds.domain_left_edge).in_units('kpc')).d
  print ('L=',L)
  kx = fftfreq(nx,d=L[0]*1.0/nx)* YTQuantity(1.,'1/kpc')
  ky = fftfreq(nx,d=L[1]*1.0/ny)* YTQuantity(1.,'1/kpc')
  kz = fftfreq(nx,d=L[2]*1.0/nz)* YTQuantity(1.,'1/kpc')

  print ('kx.shape=',kx.shape)
  kx3d, ky3d, kz3d = meshgrid(kx, ky, kz, indexing="ij")
  k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)
#  print ('k=',k[3][3][3])
  q1 = data['velocity_divergence'].in_units('1/s')
  print ('q1.shape=',q1.shape)
  print ('nx=',nx)
  q1_3d=reshape(q1,(nx,ny,nz))
  q1_3d_k = fftn(q1_3d)*YTQuantity(1.0,'kpc**3/s')
#  print ('q1_3d=',q1_3d[3][3][3])
#  print ('q1_3d_k=',q1_3d_k[3][3][3])
  phi_3d_k = q1_3d_k/k**2
  phi_3d_k[0][0][0] = 0.0 * YTQuantity(1.0,'kpc**5/s')
  phi_3d = ifftn(phi_3d_k)  * YTQuantity(1.0,'kpc**2/s')
  phi_3d_real = phi_3d.real
  print ('phi in v_mag_compress=', phi_3d_real.ravel())

  norm_factor = -0.027
  res = (data.ds.domain_right_edge[0] - data.ds.domain_left_edge[0]).in_units('kpc')/nx
  print ('res in v_mag_compress =', res)
  vx, vy, vz = gradient(phi_3d_real,edge_order=2)/res * norm_factor # * YTQuantity(1.0,'kpc**2/s')

  print ('vy in v_mag_compress =', vy.ravel())
  print ('vz in v_mag_compress =', vz.ravel())
  v_mag_3d = sqrt(vx**2 + vy**2 + vz**2) #.in_units('km/s') # compressional component

  v_mag=v_mag_3d.ravel()
  return v_mag

def _velocity(field,data):
    return  (data['x-velocity']*data['x-velocity'] + data['y-velocity']*data['y-velocity'] +data['z-velocity']*data['z-velocity'] )**0.5



def _z_den_flux(field,data):
    return abs(data['density']*data['z-velocity'])

def _CR_Energy_Density(field,data):
    return data['CREnergyDensity']* data.ds.mass_unit/data.ds.length_unit/data.ds.time_unit**2

def see(i):
  fn = get_fn(i)
  ds = yt.load(fn)
  ad=ds.all_data()
  q1 = ad["v_mag_compress"]
  fig = plt.figure()
  grid = AxesGrid(fig, (0.075,0.075,0.70,0.90), 
                nrows_ncols = (2,2),
                axes_pad = 1,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="40%",
                cbar_pad="0%")
#  fields = ['density','temperature','pressure','z-velocity']#,'CR_Energy_Density','z_den_flux']
#  fields = ['density','temperature','pressure','SN_Colour']#,'CR_Energy_Density','z_den_flux']
  fields = ['velocity','temperature','v_mag_compress','SN_Colour']#,'CR_Energy_Density','z_den_flux']
  c1 = [0.1176,0.1176,0.0]
#  c1=[0.092026 ,       0.146422,        -0.024736]
  p=yt.SlicePlot(ds,'x',fields,center=c1,axes_unit='kpc',fontsize=12)
#  p.set_unit('z-velocity','km/s')
#  p.set_zlim('density',1e-27,1e-25)
  p.set_zlim('temperature',1e4,1e8)
#  p.set_zlim('pressure',1e-13,1e-9)
#  p.set_zlim('pressure',1e-11,1e-9)
#  p.set_zlim('z-velocity',-2e3,2e3)
#  p.set_zlim('entropy',5,200)
  p.set_zlim('SN_Colour',1e2,1e6)
  p.annotate_timestamp(corner='upper_left',time_unit='Myr',text_args={'color':'black'})
  for j, field in enumerate(fields):
    plot = p.plots[field]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]
  p._setup_plots()
  p.annotate_timestamp(corner='upper_left',time_unit='Myr',text_args={'color':'black'})
  if yt.is_root():
     plt.savefig('multi_x_slice_'+str(i)+'_v_compress.png')

yt.add_field('z_den_flux',function=_z_den_flux, units  = 'Msun/yr/kpc**2')
yt.add_field('velocity',function=_velocity, units  = 'km/s')
yt.add_field('v_mag_compress',function=_v_mag_compress, units  = 'km/s')
num=[1,15,26]
num=[202,300,400,500,600,700]
num=[10,50,100,150,202,300,400,500,600,700,800]
num=[900,1000,1100]
num=range(1200,1800,100)
num=range(10,170,10)
num=range(0,2200,200)
num=range(2800,4000,100)
num=range(50,500,50)
num=[460,470,480,490]
num=[1,2,3,4]
num=range(140,260,10)
num=[5,105,255]
num=range(50,550,50)
num=[1,50,100,300,500,700]
num=range(300,900,200)
for i in num:
   see(i)


