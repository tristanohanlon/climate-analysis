from pyhdf import SD
import datetime as dt  # Python standard library datetime  module
import numpy as np
import matplotlib.pyplot as plt
f=SD.SD('2010/CERES_2010_tcc.hdf')
#view datasets
f.datasets()

tcc=f.select('tcc')
lat=f.select('latitude (dimension)')
lon=f.select('longitude (dimension)')
tcc=tcc.get()
lon=lon.get()
lat=lat.get()

a = np.reshape(tcc,180)

plt.figure()
fig, ax = plt.subplots()
ax.plot(lat,a,'-b', label='Total Cloud Cover')
#ax.axis('equal')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=3, fancybox=True, shadow=True);

#plt.plot(lat[140:161],tcc[140:161])
plt.xlabel('Latitude')
plt.ylabel('Cloud Fraction %')
plt.title('CERES Total Cloud Fraction 2010')
plt.grid(True)
plt.show()