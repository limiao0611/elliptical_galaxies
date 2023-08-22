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
  iq='velocity_divergence'
  cube = ds.covering_grid(max_level, left_edge=low,
                            dims=dims1,
                            fields=[irho, iu,iq],num_ghost_zones=3)

  rho = cube[irho].d
  ux=cube[iu].d
  q = cube[iq].d
  print (rho.size)
  print (rho.shape)

  ru=fftn(ux)
  print ("ux_k=",ru)
  
  print (ux[2][5][10])
  print ('-----------------------')
  print (ifftn(ru)[2][5][10])

# wavenumbers
  L = (ds.domain_right_edge - ds.domain_left_edge).d

  kx = np.fft.fftfreq(nx)*nx/L[0]
  ky = np.fft.fftfreq(ny)*ny/L[1]
  kz = np.fft.fftfreq(nz)*nz/L[2]

  print ('kx.shape=',kx.shape)


# bin the Fourier KE into radial kbins
  kx3d, ky3d, kz3d = meshgrid(kx, ky, kz, indexing="ij")
  k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

  print ('kx3d=',kx3d)
  print ('kx3d.shape=',kx3d.shape)

  print ('k.shape=',k.shape)

  q1 = ad['velocity_divergence']
  print ('q1=',q1)
  print ('q1.shape=',q1.shape)

  q1_3d=reshape(q1,(nx,ny,nz))
  print ('q1_3d=',q1_3d)
  print ('kx=',kx)
  print ('kx.shape=',kx.shape)

  q1_3d_k = fftn(q1_3d)
  print ('q1_3d_k=',q1_3d_k)
  print ('k.shape=',k.shape)
  print ('q1_3d.shape=',q1_3d.shape)
#  print ('k=',k)

  phi_3d_k = q1_3d_k/k**2


  print ('kx=',kx)
  print ('ky=',ky)
  print ('kz=',kz)

  phi_3d_k[0][0][0] = 0.0

  phi_3d = ifftn(phi_3d_k)
  print ('phi_3d=',phi_3d)  

  phi = phi_3d.ravel()
  print (phi)
  print (phi.real)
  print ('real mode=', dot(phi.real,phi.real))
  print ('imag mode=', dot(phi.imag,phi.imag))
  print (mean(abs( phi.real**2/phi.imag**2) ))
  print (phi.shape)
  phi_real = phi.real
  print (phi_real.shape)
  print (phi_real)
 
#  grad_fields = ds.add_gradient_fields(("gas","temperature")) 
  temp_grad_field = ds.add_gradient_fields(('gas','temperature'))
  print (temp_grad_field)
  bb = ad[temp_grad_field[0]]
  print ('bb=',bb)

#  a1= (ad['density']).in_units('g/cm**3')
#  a1_3d=reshape(a1,(nx, ny, nz))
#  aa = a1_3d - rho*YTQuantity(1.0,'g/cm**3')
#  print ('aa=',aa)




yt.add_field('number_density',function=_number_density, units  = '1/cm**3')

num = [202,300,330]
num=range(600,601)
num=[140]
for i in num:
   see(i)

