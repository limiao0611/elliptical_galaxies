# change parameters for different runs: t_over_tcool;  title name

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
from yt import *
yt.enable_parallelism()
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import glob
from numpy import *

from matplotlib.pyplot import cm 


run='n0.02-T3e6'
tc0= 90. # in Myr

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



#print get_fn(5),get_fn(50),get_fn(512),get_fn(5000),
profiles=[]
labels=[]
plot_specs=[]
num=[1,3,6,19,73,173,280]
num = [10,30,50,70,100,140,180]
num=[30,70,140]
m=0
color=cm.cool(np.linspace(0,1,len(num)))
mu=[]
a=[]
b_s=[]
b=[]
for i in num:
    filen = get_fn(i)
    ds = yt.load(filen)
    ad=ds.all_data() 
#    hot= ad.cut_region("obj['density']<3e-25")
    hot = ad["temperature"].in_units('K') > 1e4

#    print ('hot den=', ad['density'][hot])
    print (ad['density'])
    mu1= mean( log(ad['density'][hot]) ) 
    mu.append( mu1)
    a1 = std(log(ad['density'][hot]))
    a.append(a1)
#   mach=mean(ad['mach_number'])
    mach_rms = sqrt(mean(ad['mach_number']**2))
    b_s1 = sqrt((exp(a1**2)-1.)/mach_rms**2)
    b_s.append(b_s1)

    rho_mean_V = mean(ad['density'][hot])
#  rho_mean_V_log = log10(rho_mean_V)
    rho_std=std(ad['density'][hot])

    b1=  rho_std/rho_mean_V/mach_rms
    b.append(b1)

    profiles.append(yt.create_profile(ad,'density','cell_volume', weight_field=None))
    time=ds.current_time.in_units('Myr')
    t_over_tcool = (time/tc0).round(2)
    labels.append("$t/t_{c,0}=$"+ str(t_over_tcool)  )
    m=m+1

def normal(x, mu, a): # mu is the mean , a is std dev.
    return 1./sqrt(2*3.1415926)/a *exp(- (x-mu)**2/2./a**2  )


for n in range(len(num)):
    plt.plot(profiles[n].x, profiles[n]['cell_volume'], label = labels[n] ,color=color[n])
    
    x1=profiles[len(num)-1].x

    
    y1=normal(log(x1), mu[n], a[n])
    norm_factor = max(y1)/max(profiles[n]['cell_volume'])
    y1 = y1/norm_factor
    lab='log-normal,mu='+str(round(mu[n],2))+';a='+str(round(a[n],2))+';b_s='+str(b_s[n].round(2))+ ';b_rho,h='+str(b[n].round(2))
    lab='log-normal;b_s='+str(b_s[n].round(2))+ ';b_rho,h='+str(b[n].round(2))
    plt.plot(x1, y1, label=lab, color=color[n], ls='--')


plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('density [g/cm^3]')
plt.ylabel('cell volume [cm^3]')
plt.ylim(1e59,1e65)
plt.title(run)
plt.savefig('density_dist_compare_'+str(num)+'.png')

