import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
from matplotlib.pyplot import *
from numpy import * 

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

  t=ds.current_time.in_units('Myr') 
  ad = ds.all_data()

  print (len(ad['density']) , len(ad['mach_number']**2) )
  print (ad['mach_number'])
  mach_rms = sqrt(mean(ad['mach_number']**2))
  rho_mean_V = mean(ad['density'])
#  rho_mean_V_log = log10(rho_mean_V)
  rho_std=std(ad['density'])
  
#  rho_std = std(log10(ad['density']))
  print ('mean=', rho_mean_V)
  print ('std=',rho_std)  
  print (i, t, rho_mean_V, rho_std, mach_rms, rho_std/rho_mean_V)
  f=open('b.dat','a')
  print (i, t, rho_mean_V, rho_std, mach_rms, rho_std/rho_mean_V, rho_std/rho_mean_V/mach_rms,   file=f)
  f.close()
num=[50,100,200,300,400]
num=range(10,280,20)
for i in num:
   see(i)


