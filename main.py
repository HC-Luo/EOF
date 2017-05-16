import numpy as np
import matplotlib.pyplot as plt
import struct
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from matplotlib.ticker import NullFormatter
import matplotlib.ticker as mticker
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
from pylab import *


month = 1
NX = 144
NY = 73
NT = 732
MNH = M = year = 61
N = 41*21
KS = 1
KV = 8
KVT = 8
pi = np.pi
filename = '/Users/leonidas/AnacondaProjects/sx/Climatesx/copy-2/hgt500.grd'
fileread = file(filename,'rb')
z0 = []
for i in range(NX*NY*NT):
    zz, = struct.unpack('f',fileread.read(4))
    z0.append(zz)
z = np.zeros((NX,NY,NT),float)

for i in range(NX):
        for j in range(NY):
            for t in range(NT):
                z[i,j,t] = z0[i+NX*j+NX*NY*t]
hh = np.zeros((NX,NY,12,year))
h = np.zeros((41,21,61))
F = np.zeros((41,21,61))
for t in range(year):
    for k in range(12):
        for j in range(NY):
            for i in range(NX):
                hh[i,j,k,t] = z[i,j,k+t*12]
for i in range(41):
    for j in range(21):
        for t in range(61):
            h[i,j,t] = hh[i+16,j+44,month-1,t]
for i in range(41):
    for j in range(21):
        for t in range(61):
            F[i,j,t] = (h[i,j,t]-np.mean(h[i,j,0:61]))/np.std(h[i,j,0:61])
f = F.reshape((41*21,61))
A = np.dot(f,f.T)
E,EOF = np.linalg.eig(A)


m,n = np.shape(EOF)
data = np.row_stack((E.T,EOF))
index = np.argsort(data, axis=-1)
sorted_data = data[:,index[0,:]]
E,EOF = sorted_data[0,:],sorted_data[1:m+1,:]

tempt = -1
PC = np.dot(EOF.T,f)
EOF_JAN = EOF[:,tempt]
eof = EOF_JAN.reshape((41,21))

#==================================plotting=====================================
left, width = 0.1, 0.65
bottom, height = 0, 0.8
bottom_h = bottom + height - 0.23
plt.figure(1, figsize=(9, 10))


ax1 = plt.axes([left, bottom_h, width, 0.3])
ax1.plot(PC[tempt,:])
ax1.plot(np.zeros(61))
plt.xlim(0,60)
plt.ylim(-40,40)
plt.xticks(np.arange(13)*5,1943+np.arange(14)*5)
ax = gca()
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(2))


scale = '50m'
lim=[40,140,20,70]
xstep = 20
ystep = 10
ax2 = plt.axes([left, bottom, width, height],projection=ccrs.PlateCarree())
land = cfeature.NaturalEarthFeature('physical', 'land', scale,edgecolor='face',facecolor=cfeature.COLORS['land'])
#ax2.add_feature(land, facecolor='0.75',alpha = 0.5)
ax2.coastlines(scale)
lon_formatter = LongitudeFormatter(zero_direction_label=False)
lat_formatter = LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.set_xticks(np.arange(lim[0],lim[1]+xstep,xstep), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(lim[2],lim[3]+ystep,ystep), crs=ccrs.PlateCarree())

dd = 0.02
lev = np.arange(-0.03,0.17+dd,dd)
tic = np.arange(-0.03,0.17+dd*5,5*dd)
x = np.arange(40,140+2.5,2.5)
y = np.arange(20,70+2.5,2.5)
lons, lats = np.meshgrid(x, y)
#pic = ax2.contourf(lons, lats, eof.T,levels = lev,linewidths = 2,cmap = 'jet')
#plt.colorbar(pic,shrink = 1,ticks = tic,orientation = 'horizontal',anchor = (left,bottom+2))
#plt.clabel(pic, inline=1, fontsize=10,fmt = '%1.1f')
pic = ax2.contourf(lons, lats, eof.T,linewidths = 2,cmap = 'jet')
plt.colorbar(pic,shrink = 1,orientation = 'horizontal',anchor = (left,bottom+2))
#plt.savefig('/Users/leonidas/AnacondaProjects/sx/Climatesx/copy-2/Q1.pdf')
plt.show()
