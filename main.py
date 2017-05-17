import numpy as np
import matplotlib.pyplot as plt
import struct
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from matplotlib.ticker import MultipleLocator
from pylab import gca

#================================parameters=====================================
month          = 1    #January
Eigenvector    = 1    #The largest eigenvalue
longtitude_min = 0    #data minimun longtitude
longtitude_max = 360  #data maximun longtitude
latitude_min   = -90  #data minimun latitude
latitude_max   = 90   #data maximun latitude
D              = 2.5  #resolution
lonleft        = 40   #needed minimun longtitude
lonright       = 140  #needed maximun longtitude
latup          = 70   #needed minimun latitude
latdown        = 20   #needed maximun latitude
beginyear      = 1943 #the first year
year           = 61   #total year
NX             = int((longtitude_max-longtitude_min)/D)
NY             = int((latitude_max-latitude_min)/D)+1
NT             = 12*year
NXlim          = int((lonright-lonleft)/D+1)
NYlim          = int((latup-latdown)/D+1)
N              = NXlim*NYlim
#================================read data======================================
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
h = np.zeros((NXlim,NYlim,year))
F = np.zeros((NXlim,NYlim,year))
for t in range(year):
    for k in range(12):
        for j in range(NY):
            for i in range(NX):
                hh[i,j,k,t] = z[i,j,k+t*12]
for i in range(NXlim):
    for j in range(NYlim):
        for t in range(year):
            h[i,j,t] = hh[i+16,j+44,month-1,t]

#================================initializing===================================
#normalizing
for i in range(NXlim):
    for j in range(NYlim):
        for t in range(year):
            F[i,j,t] = (h[i,j,t]-np.mean(h[i,j,0:year]))/np.std(h[i,j,0:year])
f = F.reshape((N,year))

#==================calculating eigenvalues and eigenvectors=====================
A = np.dot(f,f.T)
E,EOF = np.linalg.eig(A)
m,n = np.shape(EOF)
data = np.row_stack((E.T,EOF))
index = np.argsort(data, axis=-1)
sorted_data = data[:,index[0,:]]
E,EOF = sorted_data[0,:],sorted_data[1:m+1,:]

#======================time coefficients and eof matrix=========================
tempt = -1*Eigenvector
PC = np.dot(EOF.T,f) #time coefficients
EOF_MONTH = EOF[:,tempt]
eof = EOF_MONTH.reshape((NXlim,NYlim))

#====variance contribution ratio and cumulative variance contribution ratio=====
PH = E[tempt]/sum(E) #variance contribution ratio
SS = 0
for i in np.arange(1,Eigenvector+1):
    SS = SS+E[-i]
SPH = SS/sum(E)      #cumulative variance contribution ratio
print PH,SPH

#==================================plotting=====================================
left, width = 0.1, 0.65
bottom, height = 0, 0.8
bottom_h = bottom + height - 0.23
plt.figure(1, figsize=(9, 10))

ax1 = plt.axes([left, bottom_h, width, 0.3])
ax1.plot(PC[tempt,:])
ax1.plot(np.zeros(year))
plt.xlim(0,year-1)
plt.ylim(-40,40)
plt.xticks(np.arange(13)*5,beginyear+np.arange(14)*5)
ax = gca()
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(2))

scale = '50m'
lim = [lonleft,lonright,latdown,latup]
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

#dd = 0.02
#lev = np.arange(-0.03,0.17+dd,dd)
#tic = np.arange(-0.03,0.17+dd*5,5*dd)
x = np.arange(lonleft,lonright+D,D)
y = np.arange(latdown,latup+D,D)
lons, lats = np.meshgrid(x, y)
#pic = ax2.contourf(lons, lats, eof.T,levels = lev,linewidths = 2,cmap = 'jet')
#plt.colorbar(pic,shrink = 1,ticks = tic,orientation = 'horizontal',anchor = (left,bottom+2))
#plt.clabel(pic, inline=1, fontsize=10,fmt = '%1.1f')
pic = ax2.contourf(lons, lats, eof.T,linewidths = 2,cmap = 'jet')
plt.colorbar(pic,shrink = 1,orientation = 'horizontal',anchor = (left,bottom+2))
#plt.savefig('/Users/leonidas/AnacondaProjects/sx/Climatesx/copy-2/Q1.pdf')
plt.show()
