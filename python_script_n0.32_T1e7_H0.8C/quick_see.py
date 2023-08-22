import yt
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

  c1 = [0.5,0.5,0.5]
  #c1=[0.545,0.599,0.602]
  #c1 =[0.891576, 0.130623, 0.186343]
#  sp=ds.sphere(c1, (20, 'pc'))

  ad=ds.all_data()
  print ad['gas','density']
#  print 'min_den=', min(ad['gas','density'])
  print ad['enzo','Temperature']
  print ad['enzo','z-velocity']
  low = abs(ad['index','z'])<0.003
  print 'density at low altitude=',ad['gas','density'][low]
  print 'z-vel at low altitude=',ad['enzo','z-velocity'][low]
  print 'temp at low altitude=',ad['enzo','Temperature'][low]
  high = ad['index','z']>0.48
  print 'density at high altitude=',ad['gas','density'][high]
  print 'temp at high altitude=',ad['enzo','Temperature'][high]
  print 'sum den=',sum(ad['gas','density']) 
num =[1,10,100,580]
num = [1,4]
for i in num:
   see(i)


