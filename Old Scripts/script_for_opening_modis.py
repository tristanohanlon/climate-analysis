from pyhdf import SD
import datetime as dt  # Python standard library datetime  module
import numpy as np
import matplotlib.pyplot as plt
f=SD.SD('D:/MOD06_L2/MOD06_L2.A2010001.1425.061.2017308135133.hdf')
#view datasets
f.datasets()
'Cloud_Phase_Infrared_1km'
'Cloud_Phase_Optical_Properties'
'Cloud_Fraction_Day_Nadir'
subdataset_name='Cloud_Phase_Optical_Properties'#'Cloud_Fraction_Day_Nadir'#'Cloud_Phase_Optical_Properties'
sds=f.select(subdataset_name)
data=sds.get()*sds.attributes()['scale_factor']
lat=f.select('Latitude')
lon=f.select('Longitude')
lon=lon.get()
lat=lat.get()

plt.figure()
plt.pcolormesh(lon,lat, data[:2030,:1350][::5,::5])
plt.show()
f=SD.SD('D:/MOD021KM/MOD021KM.A2010110.1855.061.2017254184547.hdf')

'EV_500_Aggr1km_RefSB'
'EV_250_Aggr1km_RefSB'
'EV_1KM_RefSB'
sds=f.select('EV_250_Aggr1km_RefSB')
data=sds.get()[0]*sds.attributes()['reflectance_scales'][0]
plt.figure(dpi=300)
plt.imshow(data,cmap='gray',vmin=0,vmax=0.3)
plt.show()