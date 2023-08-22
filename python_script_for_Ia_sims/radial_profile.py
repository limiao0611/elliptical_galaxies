import yt
import matplotlib
import matplotlib.pyplot as plt
# Load the dataset.
ds = yt.load("DD0070/sb_0070")

# Create a sphere of radius 100 kpc in the center of the box.
my_sphere = ds.sphere("c", (100.0, "pc"))

# Create a profile of the average density vs. radius.
plot = yt.ProfilePlot(my_sphere, "radius", "pressure")
#                      ,weight_field="cell_mass")

# Change the units of the radius into kpc (and not the default in cgs)

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

num=[40,60,90]
for i in num:
    filen = get_fn(i)
    ds = yt.load(filen)
    my_sphere = ds.sphere("c",(100., 'pc'))    
    plot = yt.ProfilePlot(my_sphere,'radius','pressure')
    plot.set_unit('radius', 'pc')
# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
    plot.save()
