import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
from matplotlib.pyplot import *
from numpy import * 
import numpy as np

pc= 3.086e18

T=3e6
r_fade =73. *pc

l_box = 20*r_fade

k=1.38e-16
gamma=5./3.
mu = 0.6*1.67e-24
cs = sqrt(gamma*k*T/mu)
Myr = 3.15e13
t_sc = 0.5* l_box/cs/Myr # sound crossing time of the box
print ('t_sc/Myr=',t_sc)

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



def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))

def see(i):
  fn = get_fn(i)

  ds = yt.load(fn)

  t=ds.current_time.in_units('Myr') 
  ad = ds.all_data()

  hot = ad.cut_region('obj["temperature"].in_units("K")>3e5')
#  mach_rms = sqrt(mean(ad['mach_number']**2))
#  rho_mean_V = mean(ad['density'])
#  rho_mean_V_log = log10(rho_mean_V)
#  rho_std=std(ad['density'])

#  print ('den=', ad['density'])
#  print ('log den=',  log(ad['density']) )
#  mu=mean( log(ad['density']) )
#  a = std(  log(ad['density']) )


  mach_rms_h = sqrt(mean(hot['mach_number']**2))
  print (mach_rms_h)
  rho_mean_V_h = mean(hot['density'])
  rho_std_h =std(hot['density']) 
  rho_fluc_h =  rho_std_h/rho_mean_V_h
  s_h =  log(hot['density']/rho_mean_V_h)
  s_mean_h =  mean(s_h)
  s_std_h =  std(s_h)

  s_mean_h_m, s_std_h_m = weighted_avg_and_std(s_h, hot['density'])

  print (s_mean_h, s_std_h, - s_std_h**2/2.)
  print (s_mean_h_m, s_std_h_m)


  b_rho_h = rho_std_h/rho_mean_V_h/mach_rms_h 
  b_s_h = sqrt(  exp(s_std_h**2)-1.   )/mach_rms_h


 
#  rho_std = std(log10(ad['density']))
#  print ('mean=', rho_mean_V)
#  print ('std=',rho_std)  
#  print (i, t, rho_mean_V, rho_std, mach_rms, rho_std/rho_mean_V)
  f=open('b2.dat','a')
# print (i, t, rho_mean_V, rho_std, mach_rms, rho_std/rho_mean_V, rho_std/rho_mean_V/mach_rms, mu,a  , file=f)
  print (i, t, t_sc, mach_rms_h, rho_mean_V_h, rho_std_h, rho_fluc_h,  s_mean_h,  s_std_h, s_mean_h_m, s_std_h_m,  b_rho_h ,  b_s_h   , file=f)
  f.close()
num=range(3,300,20)
#num=range(1,15,2)
#num=[50]
for i in num:
   see(i)


