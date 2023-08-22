import yt
import matplotlib
import matplotlib.pyplot as plt

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
num = [1,3,6,19,73]
for i in num:
    filen = get_fn(i)
    ds = yt.load(filen)
    my_sphere = ds.sphere("c",(100., 'pc'))    
    
# Calculate and store the bulk velocity for the sphere.
#    bulk_velocity = my_sphere.quantities.bulk_velocity()
#    my_sphere.set_field_parameter('bulk_velocity', bulk_velocity)
    pn = 'density'
    pn = 'pressure'
    pn = ['density','temperature','pressure']
    profiles.append(yt.create_profile(my_sphere,'radius',pn))

    labels.append(str(i)+('*1.5e4yr'))
    plot_specs.append(dict(linewidth=2, alpha=0.7))
#    plot = yt.ProfilePlot(my_sphere,'radius','pressure')
plot = yt.ProfilePlot.from_profiles(profiles,labels=labels,plot_specs=plot_specs)    
plot.set_unit('radius', 'pc')
#plot.set_xscale('linear')
plot.save()

