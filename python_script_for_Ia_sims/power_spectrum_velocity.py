import numpy as np
# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
from numpy import *
from matplotlib.pyplot import cm



import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 2.0 #set the value globally


from matplotlib import rc

font = {'family' : 'serif',
        'size'   : 22}


matplotlib.rc('font', **font)



n=0.08
T=1e7
tc0 = 133.

fig=plt.figure(figsize=(10,8))

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


def doit(ds, data_num):

    # a FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this.

    max_level = ds.index.max_level

    print ("max_level=",max_level)
#    ref = int(np.product(ds.ref_factors[0:max_level])) # refine factors

    ref = 1
    low = ds.domain_left_edge
    dims = ds.domain_dimensions*ref

    nx, ny, nz = dims

#    nindex_rho = 1./3.
    nindex_rho = 0.

    Kk = np.zeros( (nx//2+1, ny//2+1, nz//2+1))

    for vel in [("gas", "velocity_x"), ("gas", "velocity_y"),
                ("gas", "velocity_z")]:

        Kk += 0.5*fft_comp(ds, ("gas", "density"), vel,
                           nindex_rho, max_level, low, dims)

    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).in_units('kpc').d

    kx = np.fft.rfftfreq(nx)*nx/L[0]
    ky = np.fft.rfftfreq(ny)*ny/L[1]
    kz = np.fft.rfftfreq(nz)*nz/L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
    kmax = np.min(0.5*dims/L)

    kbins = np.arange(kmin, kmax, kmin)
    N = len(kbins)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    whichbin = np.digitize(k.flat, kbins)
    ncount = np.bincount(whichbin)

    E_spectrum = np.zeros(len(ncount)-1)

    for n in range(1,len(ncount)):
        E_spectrum[n-1] = np.sum(Kk.flat[whichbin==n])
    E_spectrum = E_spectrum[1:N]


    k = 0.5*(kbins[0:N-1] + kbins[1:N])

    index = np.argmax(E_spectrum)
    kmax = k[index]
    Emax = E_spectrum[index] * k[index]**2
    Emin=min(E_spectrum) 
    t = round(  np.log(np.exp(ds.current_time.in_units("Myr")))    ,1)
    print ("hey", len(k))
    if (data_num==1):
      plt.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5",label=r"$k^{-5/3}$", lw=2.0)
      plt.loglog([10*kmin, 10*kmin], [Emin*1e2,2*Emax* kmax**2],label=r'k$=1/(2R_{\rm{fade}}$)',ls='--',color='grey',lw=2 )

    t=round(t/tc0,1)
#    plt.loglog(k, rho_spectrum, label= r"$t/t_{c,0}=$"+str(t), color=color[m] , lw=2.5)


    plt.loglog(k[0:len(k)-1], E_spectrum[0:len(k)-1]*k[0:len(k)-1]**2 , label= r"$t/t_{c,0}=$"+str(t), color=color[m] , lw=2.5) #, label= str(t)+"Myr")
    plt.xlabel(r"$k (1/kpc) $")
#    plt.ylabel(r"$0.5 ( |vx(k)|^2+ |vy(k)|^2 + |vz(k)|^2)*k^2 $")
  
    plt.ylabel(r"$\widetilde{E}(k)$")
    plt.ylim(1e12, 1e15)
    plt.legend(loc=3,fontsize=17,ncol=2)
    plt.tick_params(axis='both', which='major',
                direction='in', length=8, right=True, top=True)
    plt.tick_params(axis='both', which='minor',
                direction='in', length=5, right=True, top=True)

    plt.savefig("velocity_power_spectrum_square"+str(data_num)+"_1.png")

def fft_comp(ds, irho, iu, nindex_rho, level, low, delta ):

    cube = ds.covering_grid(level, left_edge=low,
                            dims=delta,
                            fields=[irho, iu])

    rho = cube[irho].d
    u = cube[iu].d

    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be
    # too much information here.  fftn puts the positive freq terms in
    # the first half of the axes -- that's what we keep.  Our
    # normalization has an '8' to account for this clipping to one
    # octant.
    ru = np.fft.fftn(rho**nindex_rho * u)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
    ru = 8.0*ru/(nx*ny*nz)

    return np.abs(ru)**2


num=[1,100,250,350,370,390]
num=[1,50,130,210,290,370]
color=cm.magma(linspace(0,0.9, len(num) ))
color=cm.hot_r(np.linspace(0.3,1,len(num)))
m=0
#num=[390]
for i in num:
  fn = get_fn(i)
  ds = yt.load(fn)
  doit(ds,i)
  m=m+1

