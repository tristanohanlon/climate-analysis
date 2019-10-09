import datetime as dt  # Python standard library datetime  module
import numpy as np
from pyhdf import SD
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

f = SD.SD('Data/modis_test.hdf')
#view datasets
datasets_dic = f.datasets()


lat=f.select('Latitude')
lon=f.select('Longitude')
lon=lon.get()
lat=lat.get()


sds_obj = f.select('Cloud_Fraction_Nadir_Day') # select sds

data = sds_obj.get()*0.01 # get sds data

#for key, value in sds_obj.attributes().items():
   # print (key, value)
   # if key == 'add_offset':
   #     add_offset = value  
  #  if key == 'scale_factor':
   #     scale_factor = value

#data = (data - add_offset) * scale_factor



#print (data)
   
#data_av = np.mean(data, axis=0) 
plt.figure()
plt.contourf(lon, lat, data)

lat_bins=np.arange(-25,-5,0.5)
lat_hist=np.digitize(lat,bins=lat_bins)
list1=[]
for i in range(len(lat_bins)):
    location=lat_hist==i
    cf_lat=data[location]
    nw_data=np.nanmean(cf_lat[(cf_lat<1.0)&(cf_lat>-0.01)])
    list1.append(nw_data)
plt.figure()
plt.plot(lat_bins,list1)
plt.xlabel('Latitude')
plt.ylabel('Total Cloud Fraction')
plt.title('Total Cloud Fraction vs Latitude')
plt.grid(True)


#plt.contourf(lon,lat,data_av)


