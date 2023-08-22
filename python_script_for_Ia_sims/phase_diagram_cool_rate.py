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

def _metallicity(field,data):
    return data["SN_Colour"]/data["density"]

def _cool_rate(field,data):
    return data["cell_mass"]*data["GasEnergy"]/data["cooling_time"]

def _cool_time_inv(field,data):
    return 1./data["cooling_time"]

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


def see(i):
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()

   a ='SN_Colour_fraction'
   b='t_cool_isoc'
   c = "cell_mass"
   d=None

   plot = PhasePlot(ad, a,b, [c],weight_field=d)
   plot.save(a+ '_' + b+ '_' + c +'_' +str(i)+'_all.png')
   time = ds.current_time.in_units("Myr")
   cool_rate_hot_gas =sum(ad["cooling_rate_per_cell"])
   f1=open("cool_rate_total.dat",'a')
   print (time, cool_rate_hot_gas,file=f1)
   f1.close()


yt.add_field('metallicity',function=_metallicity)
#yt.add_field('cool_rate',function=_cool_rate,units="erg/s")
yt.add_field('cool_time_inv',function=_cool_time_inv,units="1/s")
yt.add_field('cooling_rate',function=_cooling_rate,units="erg*cm**3/s")
yt.add_field('cooling_rate_per_volume',function=_cooling_rate_per_volume,units="erg/cm**3/s")
yt.add_field('cooling_rate_per_cell',function=_cooling_rate_per_cell,units="erg/s")
yt.add_field('t_cool_isob',function=_t_cool_isob, units="Myr")
yt.add_field('t_cool_isoc',function=_t_cool_isoc, units="Myr")
yt.add_field('SN_Colour_fraction',function=_SN_Colour_fraction,units="")

num=[80,100,130]
num=range(300,5000,1000)
num=[300,500,1000,2000,3000,4000]
num=[1,30,50,70,80,90]
for i in num:
    see(i)

