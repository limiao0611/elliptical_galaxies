matplotlib.use('Agg')
import yt
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


def see(i):
   fn = get_fn(i)
   ds = yt.load(fn)
   ad = ds.all_data()
   # Get the minimum and maximum densities.
   quan="density"
   mi, ma = ad.quantities.extrema(quan)
   mi=1e-27
   ma=1e-19

# Create a transfer function to map field values to colors.
# We bump up our minimum to cut out some of the background fluid
   tf = yt.ColorTransferFunction((np.log10(mi)+1, np.log10(ma)))

# Add five Gaussians, evenly spaced between the min and
# max specified above with widths of 0.02 and using the
# spectral colormap.
   tf.add_layers(5, w=0.002, colormap="spectral")

# Choose a center for the render.
   c = [0.015, 0.015, 0.0]

# Choose a vector representing the viewing direction.
   L = [0.0, 0.1, 0.0]

# Set the width of the image.
# Decreasing or increasing this value
# results in a zoom in or out.
   W = 0.06

# The number of pixels along one side of the image.
# The final image will have Npixel^2 pixels.
   Npixels = 512

# Create a camera object.
# This object creates the images and
# can be moved and rotated.
   cam = ds.camera(c, L, W, Npixels, tf)

# Create a snapshot.
# The return value of this function could also be accepted, modified (or saved
# for later manipulation) and then put written out using write_bitmap.
# clip_ratio applies a maximum to the function, which is set to that value
# times the .std() of the array.
   cam.snapshot("volume_rendered_"+quan+str(i)+".png", clip_ratio=8.0)

num=range(20,70,5)
for i in num:
    see(i)


