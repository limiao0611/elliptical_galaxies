import yt



ds = yt.load('DD0152/sb_0152')
sp=ds.sphere([0.5,0.5,0.5], (10, 'pc'))
print sp['density']
print sp['pressure']

slc = yt.SlicePlot(ds, 'z', 'density').save()
slc = yt.SlicePlot(ds, 'x', 'density').save()
slc = yt.SlicePlot(ds, 'y', 'density').save()
slc = yt.SlicePlot(ds, 'z', 'pressure').save()
slc = yt.SlicePlot(ds, 'x', 'pressure').save()
slc = yt.SlicePlot(ds, 'y', 'pressure').save()
slc = yt.SlicePlot(ds, 'z', 'temperature').save()
slc = yt.SlicePlot(ds, 'x', 'temperature').save()
slc = yt.SlicePlot(ds, 'y', 'temperature').save()

