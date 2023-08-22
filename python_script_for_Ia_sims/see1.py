import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import yt
from matplotlib.pyplot import *

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
 

  c1 = [0.015,0.015,0.0]
  #c1=[0.545,0.599,0.602]
  #c1 =[0.891576, 0.130623, 0.186343]
#  sp=ds.sphere(c1, (1, 'kpc'))

  slc = yt.SlicePlot(ds, 'x', ['density','temperature','pressure','z-velocity','Cooling_Time'],center=c1,width=((150.,'pc'),(750.,'pc')))
#  slc.zoom(10)
#  colorbar=slc.cb
#  colorbar(fraction=1.5)
  slc.set_zlim('density',5e-29,2e-21)
  slc.set_zlim('temperature',1e2,1e8)
  slc.set_font_size(50)
  slc.annotate_line([0.015,0.0,0.03],[0.015,0.03,0.03],plot_args={'color':'white','linestyle':'--','linewidth':10.} )
  slc.annotate_line([0.015,0.0,-0.03],[0.015,0.03,-0.03],plot_args={'color':'white','linestyle':'--','linewidth':10.} )
  slc.save()
#  slc = yt.SlicePlot(ds, 'y', ['density','temperature','pressure','z-velocity'],center=c1).save()
#  slc = yt.SlicePlot(ds, 'z', ['density','temperature','pressure','z-velocity'],center=c1).save()


#  c2 = [0.087,0.087,0.17]
  c2= [0.015,0.015,0.0]
#  slc = yt.SlicePlot(ds, 'z', ['density','temperature','pressure','z-velocity'],center=c2).save()
#  slc = yt.SlicePlot(ds, 'z', ['density','temperature','pressure','z-velocity'],center=c2).save(str(c2)+'.png')
#  slc = yt.SlicePlot(ds, 'y', ['density','temperature','pressure','z-velocity'],center=c2).save()
#  slc = yt.SlicePlot(ds, 'x', ['density','temperature','pressure','z-velocity'],center=c2).save()

num =[0,1,10,100,200,261]
num = [0,1,10,100,161]
num = [400,547]
num = [575]
num = [0,1,10,36]
num = [0,1,10,30]
num=[380,500,625]
num=[90,100,110,120,130,140,150,160,166]
num=[150,160,165]
for i in num:
   see(i)


