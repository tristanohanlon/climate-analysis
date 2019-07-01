# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This will create a dataset of time averaged global cloud cover, liquid water cloud fraction and ice water cloud fraction with altitude.
The code can select either global or southern ocean data.
Data is stored in the 2D arrays: 

gfdl_tcc_alt
gfdl_tclw_alt
gfdl_tciw_alt

[:,0] = alt
[:,1] = cloud fraction

The data has already been scaled and offsetted.
"""
import time
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

#Home PC
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/GFDL/', 'r') # Uni Laptop

dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/cl_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
cf = dataset.variables['cl'][:] #Extract cloud_area_fraction_in_atmosphere_layer, keyed to time, lev, lon and lat (420, 33, 180, 288)

dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/clw_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
lw = dataset.variables['clw'][:] #Extract Mass Fraction of Cloud Liquid Water (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/cli_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
iw = dataset.variables['cli'][:] #Extract Mass Fraction of Cloud Ice (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/ta_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
temp = dataset.variables['ta'][:] #Extract air temperature, keyed to time, plev, lon and lat (420, 19, 180, 288)
plev = dataset.variables['plev'][:] #Extract air pressure corresponding to temperature (19,)

dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/pfull_AERmon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
pressure = dataset.variables['pfull'][:] #Extract air pressure, keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/pfull_AERmon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
ps = dataset.variables['ps'][:] #Extract surface air pressure, keyed to time, lon and lat (420, 180, 288)

end = time.time()
print('Importing data took:', end - start, 's')

#---------------tcc-------------#
start = time.time()

cf = np.array(tcc)
cf= np.mean(cf, axis=0) #Average fraction of cloud cover over time
cf = np.mean(cf, axis=-1) #Average fraction of cloud cover over longitude
"""
#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
tcc = np.vstack((lat, tcc)).T #creates a (721,38) array

#Select latitudes over the southern ocean
tcc = tcc[tcc[:,0]>=-70]
tcc = tcc[tcc[:,0]<=-50]

#Split the combined array into just the tcc data, eliminating the first coloumn of latitude
tcc = tcc[:,1:38]
"""
tcc = np.mean(tcc, axis=-1) #Average fraction of cloud cover over latitude
ecmwf_tcc_plevel = np.vstack((alt, cf)).T 

end = time.time()
print('Averaging tcc data took:', end - start, 's')

#---------------tciw-------------#

start = time.time()

tciw = np.array(tciw)
tciw = np.mean(tciw, axis=0) #Average specific cloud ice water content (kg/kg) over time
tciw = np.mean(tciw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude
"""
#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
tciw = np.vstack((lat, tciw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
tciw = tciw[tciw[:,0]>=-70]
tciw = tciw[tciw[:,0]<=-50]

#Split the combined array into just the tciw data, eliminating the first coloumn of latitude
tciw = tciw[:,1:38]
"""
tciw = np.mean(tciw, axis=-1) * 10000 #Average specific cloud ice water content (kg/kg) over latitude, scale up by 10000

end = time.time()
print('Averaging tciw data took:', end - start, 's')

#---------------tclw-------------#

start = time.time()

tclw = np.array(tclw)
tclw = np.mean(tclw, axis=0) #Average specific cloud liquid water content (kg/kg) over time
tclw = np.mean(tclw, axis=-1) #Average specific cloud liquid water content (kg/kg) over longitude
"""
#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
tclw = np.vstack((lat, tclw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
tclw = tclw[tclw[:,0]>=-70]
tclw = tclw[tclw[:,0]<=-50]

#Split the combined array into just the tclw data, eliminating the first coloumn of latitude
tclw = tclw[:,1:38]
"""
tclw = np.mean(tclw, axis=-1) * 10000 #Average specific cloud liquid water content (kg/kg) over latitude, scale up by 10000

end = time.time()
print('Averaging tclw data took:', end - start, 's')

#---------------temperature-------------#

start = time.time()

temp = np.array(temp)
temp = np.mean(temp, axis=0) #Average air temperature over time
temp = np.mean(temp, axis=-1) #Average air temperature over longitude
"""
#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
temp = np.vstack((lat, temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
temp = temp[temp[:,0]>=-70]
temp = temp[temp[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
temp = temp[:,1:38]
"""
temp = np.mean(temp, axis=-1)  #Average air temperature over latitude

end = time.time()
print('Averaging temp data took:', end - start, 's')

#---------------combine-------------#

start = time.time()

# Join the two lists as if they were two columns side by side, into a list of two elements each
gfdl_tcc_alt = np.vstack((plevel, tcc)).T 
gfdl_tclw_alt = np.vstack((plevel, tclw)).T
gfdl_tciw_alt = np.vstack((plevel, tciw)).T
gfdl_temp_alt = np.vstack((plevel, temp)).T

end = time.time()
print('Creating the combined arrays took:', end - start, 's')

#----------------------------#

plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(gfdl_tcc_alt[:,1],gfdl_tcc_alt[:,0], '-r', label='Fraction Cloud Cover')
ax2.plot(gfdl_tclw_alt[:,1],gfdl_tclw_alt[:,0], '-b', label='Specific Cloud Liquid Water Content')
ax2.plot(gfdl_tciw_alt[:,1],gfdl_tciw_alt[:,0], '--b', label='Specific Cloud Ice Water Content')

#ax.axis('equal')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Cloud Liquid and Ice Water Content (kg/kg) x $10^{-4}$')
ax1.set_ylabel('Pressure Level (hPa)')
plt.title('Cloud Fraction and Phase vs Pressure Level GFDL.AM4 2010')
plt.gca().invert_yaxis()

plt.grid(True)
plt.show()
