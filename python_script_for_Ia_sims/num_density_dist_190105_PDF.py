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
#variable n should be number of curves to plot (I skipped this earlier thinking that it is obvious when looking at picture - srry my bad mistake xD): n=len(array_of_curves_to_plot)
#version 1:
#for i,c in zip(range(n),color):
#   ax1.plot(x, y,c=c)

import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 2.0 #set the value globally

from matplotlib import rc

font = {'family' : 'serif',
        'sans-serif':'Helvetica',
#        'weight' : 'bold',
        'size'   : 23}

matplotlib.rc('font', **font)

run='n0.02-T3e6'


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

def ReadCoolingCurve():
    FileIn = open("cool_rates.in_300K")
    data = loadtxt(FileIn)
    FileIn.close()
    return (data[:,0],data[:,2])

(LogT, LogCoolRate) = ReadCoolingCurve()

def _cooling_rate(field,data):
    return 10.0**interp(log10(data["temperature"]), LogT, LogCoolRate)* YTQuantity(1,"erg*cm**3/s")



m_H = YTQuantity(1.67e-24,'g')
mu = 0.6*m_H
k_b =YTQuantity(1.38e-16, 'erg/K')

def _cooling_rate_per_volume(field,data):
    return 1.76*0.4* data["cooling_rate"] *data["density"]**2/m_H/mu

def _t_cool_isob(field,data):
    return 2.5 * data["density"]/mu * k_b * data["temperature"]/data["cooling_rate_per_volume"]

def _t_cool_isoc(field,data):
    return 3./5* data["t_cool_isob"]

def _cooling_rate_per_cell(field,data):
    return data["cooling_rate_per_volume"]*data["cell_volume"]

def _SN_Colour_fraction(field,data):
    norm_factor = 30675.  # 1 mass in code unit/1msun
    return data["SN_Colour"]/data["density"]/norm_factor

yt.add_field('cooling_rate',function=_cooling_rate,units="erg*cm**3/s")
yt.add_field('cooling_rate_per_volume',function=_cooling_rate_per_volume,units="erg/cm**3/s")
yt.add_field('cooling_rate_per_cell',function=_cooling_rate_per_cell,units="erg/s")
yt.add_field('t_cool_isob',function=_t_cool_isob, units="Myr")
yt.add_field('t_cool_isoc',function=_t_cool_isoc, units="Myr")
yt.add_field('SN_Colour_fraction',function=_SN_Colour_fraction,units="")


#print get_fn(5),get_fn(50),get_fn(512),get_fn(5000),
profiles=[]
labels=[]
plot_specs=[]
num=[1,3,6,19,73,173,280]
num = [10,30,50,70,100,140,180]
n=0
color=cm.cool(np.linspace(0,1,len(num)))
color=cm.hot_r(np.linspace(0.3,1,len(num)))
for i in num:
    print ("color=",color)
    print ("color[1]=",color[1])
    filen = get_fn(i)
    ds = yt.load(filen)
    ad=ds.all_data()
    xmin = 6e-4
    xmax=2.4e2
    nbin=64 
    profiles.append(yt.create_profile(ad,'number_density','cell_volume', weight_field=None, fractional=True, extrema={'number_density': (xmin,xmax)},  n_bins = nbin) )
    t_over_tcool = round(i/90.,1)
    labels.append("$t/t_{c,0}=$"+ str(t_over_tcool)  )
#    labels.append(str(i)+('Myr'))
#    plot_specs.append(dict(linewidth=2, alpha=0.7,c=color[n]))
#    plot = yt.ProfilePlot(my_sphere,'radius','pressure')
    n=n+1


log10_bin_size = log10(xmax/xmin)/nbin
fig=plt.figure(figsize=(11,8))
for n in range(len(num)):
    plt.plot(profiles[n].x, profiles[n]['cell_volume']/log10_bin_size, label = labels[n] ,color=color[n])



plt.legend(fontsize=22)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number Density n [cm$^{-3}$]')
#plt.ylabel('Volume Fraction')
plt.ylabel(r'$\frac{dV}{V\ dlog_{10} n}$ ')
ymin=5e-7/log10_bin_size
ymax=1/log10_bin_size
plt.ylim(ymin,ymax)
plt.title(run)
plt.plot([0.02, 0.02],[ymin,ymax], c='navy', ls='--',lw=2)
plt.tick_params(axis='both', which='major',
                direction='in', length=8, right=True, top=True)
plt.tick_params(axis='both', which='minor',
                direction='in', length=5, right=True, top=True)

plt.savefig('num_density_dist'+str(num)+'_fractional_2_PDF.png',bbox='tight')

#plot = yt.ProfilePlot.from_profiles(profiles,labels=labels,plot_specs=plot_specs)    
#plot.set_font(({'family':'sans-serif', 'style':'italic',
#               'weight':'bold', 'size':24, 'color':'blue'}))
#plot.set_font_size(24)
#plot.set_unit('radius', 'pc')
#plot.set_xscale('linear')
#plot.save('num_density_dist'+str(num)+'_fractional_2.png')

