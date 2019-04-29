from pyhdf import SD
import datetime as dt  # Python standard library datetime  module
import numpy as np
import matplotlib.pyplot as plt
f=SD.SD('2016/CER_SSF1deg-Month_Terra-MODIS_Edition4A_401405.201601')
#view datasets
f.datasets()

subdataset_name='cld_amount_zon'#'Cloud_Fraction_Day_Nadir'#'Cloud_Phase_Optical_Properties'
sds=f.select(subdataset_name)
#data=sds.get()*sds.attributes()['scale_factor']
lat=f.select('latitude')
lon=f.select('longitude')
data=sds.get()
lon=lon.get()
lat=lat.get()
tcc=data[-1,:] #selects the last index of the array (total cloud)
#data_av=np.mean(data, axis=4)

plt.figure()
plt.plot(lat[140:161],tcc[140:161])
plt.show()
