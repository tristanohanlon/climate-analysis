from pyhdf import SD
import datetime as dt  # Python standard library datetime  module
import numpy as np
import matplotlib.pyplot as plt
f=SD.SD('2016/CER_SSF1deg-Month_Terra-MODIS_Edition4A_401405.201601')
#view datasets
f.datasets()


tc=f.select('cld_amount_zon') 
tlc=f.select('cld_amount_liq_zon') 
tic=f.select('cld_amount_ice_zon') 
lat=f.select('latitude')
lon=f.select('longitude')
tc=tc.get()
tlc=tlc.get()
tic=tic.get()
lon=lon.get()
lat=lat.get()
tcc=tc[-1,:] #selects the last index of the array (total cloud)
tlcc=tlc[-1,:]
ticc=tic[-1,:]



plt.figure()
fig, ax = plt.subplots()
ax.plot(lat[140:161],tcc[140:161], '-b', label='Total Cloud Cover')
ax.plot(lat[140:161],tlcc[140:161], '--r', label='Liquid Cloud Cover')
ax.plot(lat[140:161],ticc[140:161], '--g', label='Ice Cloud Cover')
#ax.axis('equal')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=3, fancybox=True, shadow=True);

#plt.plot(lat[140:161],tcc[140:161])
plt.xlabel('Latitude')
plt.ylabel('Cloud Fraction %')
plt.title('CERES Total Cloud Fraction over the Southern Ocean for 01/2016')
plt.grid(True)
plt.show()