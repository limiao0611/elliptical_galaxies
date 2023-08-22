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
from scipy import stats

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

def _velocity(field,data):
    return  (data['x-velocity']*data['x-velocity'] + data['y-velocity']*data['y-velocity'] +data['z-velocity']*data['z-velocity'] )**0.5


def _number_density(field,data):
#    print ('data1.shape=',data['density'].shape)
    return abs(data['density']/YTQuantity(1.67e-24,'g'))

def _num_density(field,data):
    nx = ny = nz = len(data['density'])
    if (nx > 16.1):
      nx, ny, nz = data.ds.domain_dimensions
    num=range(1, nx*ny*nz)*YTQuantity(1.0,'1/cm**3')
    return num


def _phi_v(field,data):
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
  q1_3d=reshape(q1,(nx,ny,nz))
  q1_3d_k = fftn(q1_3d)*YTQuantity(1.0,'kpc**3/s')
#  print ('q1_3d=',q1_3d[3][3][3])
#  print ('q1_3d_k=',q1_3d_k[3][3][3])
  phi_3d_k = q1_3d_k/k**2
  phi_3d_k[0][0][0] = 0.0* YTQuantity(1.0,'kpc**5/s')
  phi_3d = ifftn(phi_3d_k)* YTQuantity(1.0,'kpc**2/s')
  phi = phi_3d.ravel()
  phi_real = phi.real

  return phi_real


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
   
 
def see(i):

#  global kx, ky, kz, nx, ny, nz    
  fn = get_fn(i)
  ds = yt.load(fn)

  nx, ny, nz = ds.domain_dimensions

  f=open('derived_field_list','w')
  print (ds.derived_field_list, file=f)
  ad = ds.all_data()

  c1 = [0.5,0.5,0.5]


#  cr =  ds.box([0,0,-0.0004],[0.05,0.05,0.])
#  a = ds.retrieve_ghost_zones(3,'density')

#  grad_vx= ds.add_gradient_fields(("gas","velocity_x"))
#  print (ad[grad_vx[0]])
#  grad_n=ds.add_gradient_fields(('gas','number_density'))
#  print ('grad_n=',ad[grad_n[0]])

  grad_v_div = ds.add_gradient_fields(('gas','velocity_divergence'))
  print ('grad_v_div=',ad[grad_v_div[3]])

  phi_v1=ad['phi_v']
  phi_v1_3d = reshape(phi_v1,(nx, ny, nz))
  res = (ds.domain_right_edge[0] - ds.domain_left_edge[0]).in_units('kpc')/nx
  vx, vy, vz = gradient(phi_v1_3d,edge_order=2)/res
  
  print ('vx=',vx.in_units('km/s'))
  print ('vy=',vy.in_units('km/s'))
#  print ('vz=',vz.in_units('km/s'))
  print ('vx.shape=',vx.shape)

  print ('x-velocity=',ad['x-velocity'].in_units('km/s'))
  print ('y-velocity=',ad['y-velocity'].in_units('km/s'))
  print ('z-velocity=',ad['z-velocity'].in_units('km/s'))

  print ('x-velocity=',ad['velocity_x'].in_units('km/s'))
  print ('y-velocity=',ad['velocity_y'].in_units('km/s'))
  print ('z-velocity=',ad['velocity_z'].in_units('km/s'))


  print ('vx.ravel=',(vx.in_units('km/s')).ravel() )
  print ('res=',res)


  ux = reshape(ad['x-velocity'],(nx, ny,nz))
  uy = reshape(ad['y-velocity'],(nx, ny,nz))
  uz = reshape(ad['z-velocity'],(nx, ny,nz))

 
#   norm_fac=0.000001
  f=open('norm_fac.dat','a')
  ra = arange(-1.7, -1.3, 0.1)
  for norm_fac in arange(-0.04,-0.02, 0.001):
    wx= (ux -  vx*norm_fac).in_units('km/s')
    wy= (uy -  vy*norm_fac).in_units('km/s')
    wz= (uz -  vz*norm_fac).in_units('km/s')

    divergence_w_3d = ((gradient(wx[0],edge_order=2)+gradient(wy[1], edge_order=2)+ gradient(wz[2], edge_order=2))/res).in_units('1/s')
    divergence_w=divergence_w_3d.ravel()
    print ('norm_fac, divergence_w=',norm_fac,mean(abs(divergence_w)),  mean(divergence_w), std(divergence_w) ,file=f)
    print ('norm_fac, divergence_w=',norm_fac,mean(abs(divergence_w)),  mean(divergence_w), std(divergence_w) )

#    print ('wx=',wx)
#    print ('ux=',ux.in_units('km/s'))

  f.close()

  norm_fac= -0.027
  print ('normed compressional vx.ravel=',(norm_fac*vx.in_units('km/s')).ravel() )
  print ('x-velocity=',ad['x-velocity'].in_units('km/s'))

  wx= (ux -  vx*norm_fac).in_units('km/s')
  wy= (uy -  vy*norm_fac).in_units('km/s')
  wz= (uz -  vz*norm_fac).in_units('km/s')

  wx =wx.ravel()
  wy =wy.ravel()
  wz =wz.ravel()

  vx =vx.ravel()*norm_fac
  vy =vy.ravel()*norm_fac
  vz =vz.ravel()*norm_fac


  u_mag = sqrt(ad['x-velocity']**2 + ad['y-velocity']**2 + ad['z-velocity']**2).in_units('km/s')
  w_mag = sqrt(wx**2 + wy**2 + wz**2).in_units('km/s') # solenoidal component
  v_mag = sqrt(vx**2 + vy**2 + vz**2).in_units('km/s') # compressional component
 
  print ('total u_mag=',u_mag)
  print ('total velocity=',ad['velocity'])
  print ('solenoidal w_mag=',w_mag)
  print ('compressional v_mag=',v_mag)

  print ('mean, std total velocity=',mean(ad['velocity']), std(ad['velocity']))
  print ('mean,std u_mag',mean(u_mag), std(u_mag))
  print ('mean,std compressional v_mag=',mean(v_mag), std(v_mag))
  print ('mean,std solenoidal w_mag=',mean(w_mag) , std(w_mag))


  print ('v_mag_compress=',ad["v_mag_compress"])
  print ('mean, std v_mag_compress=',mean(ad["v_mag_compress"]),  std(ad["v_mag_compress"]))

  print ('phi in main = ', phi_v1)
  print ('vy in main =', vy)
  print ('vz in main =', vz)


yt.add_field(('gas','phi_v'),function=_phi_v,units='kpc**2/s')
yt.add_field('number_density',function=_number_density, units  = '1/cm**3')
yt.add_field('num_density',function=_num_density, units  = '1/cm**3')
yt.add_field('velocity',function=_velocity, units  = 'km/s')
yt.add_field('v_mag_compress',function=_v_mag_compress, units  = 'km/s')

num=[600]
for i in num:
   see(i)

