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
from matplotlib import rc

font = {'family' : 'serif',
        'sans-serif':'Helvetica',
#        'weight' : 'bold',
        'size'   : 21}

matplotlib.rc('font', **font)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)




run='n0.08-T1e7'
tc0= 133. # in Myr

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
num=[100,300,500]
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
    mu1= mean( log(ad['number_density'][hot]) ) 
    mu.append( mu1)
    a1 = std(log(ad['number_density'][hot]))
    a.append(a1)
#   mach=mean(ad['mach_number'])
    mach_rms = sqrt(mean(ad['mach_number']**2))
    b_s1 = sqrt((exp(a1**2)-1.)/mach_rms**2)
    b_s.append(b_s1)

    rho_mean_V = mean(ad['number_density'][hot])
#  rho_mean_V_log = log10(rho_mean_V)
    rho_std=std(ad['number_density'][hot])

    b1=  rho_std/rho_mean_V/mach_rms
    b.append(b1)

    print ('a1, rho_std=',a1,  rho_std)
    print ('sigma_s^2, ln(1+sigma_rho^2)=', a1**2, log(1+ (rho_std/rho_mean_V) **2)  )


    profiles.append(yt.create_profile(ad,'number_density','cell_volume', weight_field=None,fractional=True, extrema={'number_density': (6e-3, 2.4e2)}) )
    time=ds.current_time.in_units('Myr')
    t_over_tcool = (time/tc0).round(1)
    labels.append("$t/t_{c,0}=$"+ str(t_over_tcool)+'; $b_s=$'+str(b_s[m].round(1))+ ';$b_{n,h}=$'+str(b[m].round(1)) )
    m=m+1

def normal(x, mu, a): # mu is the mean , a is std dev.
    return 1./sqrt(2*3.1415926)/a *exp(- (x-mu)**2/2./a**2  )

fig=plt.figure(figsize=(10,8))
for n in range(len(num)):
    plt.plot(profiles[n].x, profiles[n]['cell_volume'], label = labels[n] ,color=color[n],lw=2.5)
    
    x1=profiles[len(num)-1].x

    
    y1=normal(log(x1), mu[n], a[n])
    norm_factor = max(y1)/max(profiles[n]['cell_volume'])
    y1 = y1/norm_factor
#    lab='log-normal,mu='+str(round(mu[n],2))+';a='+str(round(a[n],2))+';b_s='+str(b_s[n].round(2))+ ';b_rho,h='+str(b[n].round(2))
#    lab='log-normal;b_s='+str(b_s[n].round(2))+ ';b_rho,h='+str(b[n].round(2))
#    lab='log-normal'
#    lab=''
    plt.plot(x1, y1,color=color[n], ls='--')


plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number Density [1/cm$^3$]')
plt.ylabel('Volume Fraction')
plt.ylim(1e-7,1)
plt.legend(loc=1)
plt.title(run)
plt.savefig('density_dist_compare_'+str(num)+'_2.png')

