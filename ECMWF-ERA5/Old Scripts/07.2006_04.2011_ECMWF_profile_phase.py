# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This will create a dataset of time averaged global cloud cover, liquid water cloud fraction and ice water cloud fraction with altitude.
The code can select either global or southern ocean data.

Data is stored in the 2D arrays: 



[:,0] = alt
[:,1] = cloud fraction

The data has already been scaled and offsetted.
"""
import time
import numpy as np
import math
from netCDF4 import Dataset
import matplotlib.pyplot as plt

############################################################################### 2006

# Uni Laptop
dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2006_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
#dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2006_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

lat = dataset.variables['latitude'][:] #Extract latitude data (721)
plevel = dataset.variables['level'][:] #Extract pressure (millibars) level data (37 levels)
a_cf = dataset.variables['cc'][:] #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
a_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
a_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
a_temp = dataset.variables['t'][:] #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')

#---------------cf-------------#
start = time.time()

a_cf = np.array(a_cf)
a_cf = np.take(a_cf, [6,7,8,9,10,11], axis=0) #select months July to December 2006
a_cf = np.mean(a_cf, axis=0) #Average fraction of cloud cover over time
a_cf = np.mean(a_cf, axis=-1) #Average fraction of cloud cover over longitude
a_cf_alt_lat = a_cf

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
a_cf_so = np.vstack((lat, a_cf)).T #creates a (721,38) array

#Select latitudes over the southern ocean
a_cf_so = a_cf_so[a_cf_so[:,0]>=-70]
a_cf_so = a_cf_so[a_cf_so[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
a_cf_so = a_cf_so[:,1:38]

a_cf = np.mean(a_cf, axis=-1) #Average fraction of cloud cover over latitude
a_cf_so = np.mean(a_cf_so, axis=0) #Average southern ocean fraction of cloud cover over latitude

end = time.time()
print('Averaging cf data took:', end - start, 's')

#---------------iw-------------#

start = time.time()

a_iw = np.array(a_iw)
a_iw = np.take(a_iw, [6,7,8,9,10,11], axis=0) #select months July to December 2006
a_iw = np.mean(a_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
a_iw = np.mean(a_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude
a_iw_alt_lat = a_iw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
a_iw_so = np.vstack((lat, a_iw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
a_iw_so = a_iw_so[a_iw_so[:,0]>=-70]
a_iw_so = a_iw_so[a_iw_so[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
a_iw_so = a_iw_so[:,1:38]

a_iw = np.mean(a_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over latitude
a_iw_so = np.mean(a_iw_so, axis=0) #Average southern ocean specific cloud ice water content (kg/kg) over latitude

end = time.time()
print('Averaging iw data took:', end - start, 's')

#---------------lw-------------#

start = time.time()

a_lw = np.array(a_lw)
a_lw = np.take(a_lw, [6,7,8,9,10,11], axis=0) #select months July to December 2006
a_lw = np.mean(a_lw, axis=0) #Average specific cloud liquid water content (kg/kg) over time
a_lw = np.mean(a_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over longitude
a_lw_alt_lat = a_lw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
a_lw_so = np.vstack((lat, a_lw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
a_lw_so = a_lw_so[a_lw_so[:,0]>=-70]
a_lw_so = a_lw_so[a_lw_so[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
a_lw_so = a_lw_so[:,1:38]

a_lw = np.mean(a_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over latitude
a_lw_so = np.mean(a_lw_so, axis=0) #Average southern ocean specific cloud liquid water content (kg/kg) over latitude

end = time.time()
print('Averaging lw data took:', end - start, 's')

#---------------temperature-------------#

start = time.time()

a_temp = np.array(a_temp)
a_temp = np.take(a_temp, [6,7,8,9,10,11], axis=0) #select months July to December 2006
a_temp = np.mean(a_temp, axis=0) #Average air temperature over time
a_temp = np.mean(a_temp, axis=-1) #Average air temperature over longitude
a_temp_alt_lat = a_temp

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
a_temp_so = np.vstack((lat, a_temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
a_temp_so = a_temp_so[a_temp_so[:,0]>=-70]
a_temp_so = a_temp_so[a_temp_so[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
a_temp_so = a_temp_so[:,1:38]

a_temp = np.mean(a_temp, axis=-1) #Average air temperature over latitude
a_temp_so = np.mean(a_temp_so, axis=0)  #Average southern ocean air temperature over latitude

end = time.time()
print('Averaging temp data took:', end - start, 's')

############################################################################### 2007

# Uni Laptop
dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2007_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
#dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2007_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

b_cf = dataset.variables['cc'][:] #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
b_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
b_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
b_temp = dataset.variables['t'][:] #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')

#---------------cf-------------#
start = time.time()

b_cf = np.array(b_cf)
b_cf = np.mean(b_cf, axis=0) #Average fraction of cloud cover over time
b_cf = np.mean(b_cf, axis=-1) #Average fraction of cloud cover over longitude
b_cf_alt_lat = b_cf

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
b_cf_so = np.vstack((lat, b_cf)).T #creates a (721,38) array

#Select latitudes over the southern ocean
b_cf_so = b_cf_so[b_cf_so[:,0]>=-70]
b_cf_so = b_cf_so[b_cf_so[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
b_cf_so = b_cf_so[:,1:38]

b_cf = np.mean(b_cf, axis=-1) #Average fraction of cloud cover over latitude
b_cf_so = np.mean(b_cf_so, axis=0) #Average southern ocean fraction of cloud cover over latitude

end = time.time()
print('Averaging cf data took:', end - start, 's')

#---------------iw-------------#

start = time.time()

b_iw = np.array(b_iw)
b_iw = np.mean(b_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
b_iw = np.mean(b_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude
b_iw_alt_lat = b_iw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
b_iw_so = np.vstack((lat, b_iw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
b_iw_so = b_iw_so[b_iw_so[:,0]>=-70]
b_iw_so = b_iw_so[b_iw_so[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
b_iw_so = b_iw_so[:,1:38]

b_iw = np.mean(b_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over latitude
b_iw_so = np.mean(b_iw_so, axis=0) #Average southern ocean specific cloud ice water content (kg/kg) over latitude

end = time.time()
print('Averaging iw data took:', end - start, 's')

#---------------lw-------------#

start = time.time()

b_lw = np.array(b_lw)
b_lw = np.mean(b_lw, axis=0) #Average specific cloud liquid water content (kg/kg) over time
b_lw = np.mean(b_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over longitude
b_lw_alt_lat = b_lw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
b_lw_so = np.vstack((lat, b_lw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
b_lw_so = b_lw_so[b_lw_so[:,0]>=-70]
b_lw_so = b_lw_so[b_lw_so[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
b_lw_so = b_lw_so[:,1:38]

b_lw = np.mean(b_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over latitude
b_lw_so = np.mean(b_lw_so, axis=0) #Average southern ocean specific cloud liquid water content (kg/kg) over latitude

end = time.time()
print('Averaging lw data took:', end - start, 's')

#---------------temperature-------------#

start = time.time()

b_temp = np.array(b_temp)
b_temp = np.mean(b_temp, axis=0) #Average air temperature over time
b_temp = np.mean(b_temp, axis=-1) #Average air temperature over longitude
b_temp_alt_lat = b_temp

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
b_temp_so = np.vstack((lat, b_temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
b_temp_so = b_temp_so[b_temp_so[:,0]>=-70]
b_temp_so = b_temp_so[b_temp_so[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
b_temp_so = b_temp_so[:,1:38]

b_temp = np.mean(b_temp, axis=-1) #Average air temperature over latitude
b_temp_so = np.mean(b_temp_so, axis=0)  #Average southern ocean air temperature over latitude

end = time.time()
print('Averaging temp data took:', end - start, 's')

############################################################################### 2008

# Uni Laptop
dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2008_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
#dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2008_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

c_cf = dataset.variables['cc'][:] #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
c_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
c_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
c_temp = dataset.variables['t'][:] #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')

#---------------cf-------------#
start = time.time()

c_cf = np.array(c_cf)
c_cf = np.mean(c_cf, axis=0) #Average fraction of cloud cover over time
c_cf = np.mean(c_cf, axis=-1) #Average fraction of cloud cover over longitude
c_cf_alt_lat = c_cf

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
c_cf_so = np.vstack((lat, c_cf)).T #creates a (721,38) array

#Select latitudes over the southern ocean
c_cf_so = c_cf_so[c_cf_so[:,0]>=-70]
c_cf_so = c_cf_so[c_cf_so[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
c_cf_so = c_cf_so[:,1:38]

c_cf = np.mean(c_cf, axis=-1) #Average fraction of cloud cover over latitude
c_cf_so = np.mean(c_cf_so, axis=0) #Average southern ocean fraction of cloud cover over latitude

end = time.time()
print('Averaging cf data took:', end - start, 's')

#---------------iw-------------#

start = time.time()

c_iw = np.array(c_iw)
c_iw = np.mean(c_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
c_iw = np.mean(c_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude
c_iw_alt_lat = c_iw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
c_iw_so = np.vstack((lat, c_iw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
c_iw_so = c_iw_so[c_iw_so[:,0]>=-70]
c_iw_so = c_iw_so[c_iw_so[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
c_iw_so = c_iw_so[:,1:38]

c_iw = np.mean(c_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over latitude
c_iw_so = np.mean(c_iw_so, axis=0) #Average southern ocean specific cloud ice water content (kg/kg) over latitude

end = time.time()
print('Averaging iw data took:', end - start, 's')

#---------------lw-------------#

start = time.time()

c_lw = np.array(c_lw)
c_lw = np.mean(c_lw, axis=0) #Average specific cloud liquid water content (kg/kg) over time
c_lw = np.mean(c_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over longitude
c_lw_alt_lat = c_lw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
c_lw_so = np.vstack((lat, c_lw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
c_lw_so = c_lw_so[c_lw_so[:,0]>=-70]
c_lw_so = c_lw_so[c_lw_so[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
c_lw_so = c_lw_so[:,1:38]

c_lw = np.mean(c_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over latitude
c_lw_so = np.mean(c_lw_so, axis=0) #Average southern ocean specific cloud liquid water content (kg/kg) over latitude

end = time.time()
print('Averaging lw data took:', end - start, 's')

#---------------temperature-------------#

start = time.time()

c_temp = np.array(c_temp)
c_temp = np.mean(c_temp, axis=0) #Average air temperature over time
c_temp = np.mean(c_temp, axis=-1) #Average air temperature over longitude
c_temp_alt_lat = c_temp

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
c_temp_so = np.vstack((lat, c_temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
c_temp_so = c_temp_so[c_temp_so[:,0]>=-70]
c_temp_so = c_temp_so[c_temp_so[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
c_temp_so = c_temp_so[:,1:38]

c_temp = np.mean(c_temp, axis=-1) #Average air temperature over latitude
c_temp_so = np.mean(c_temp_so, axis=0)  #Average southern ocean air temperature over latitude

end = time.time()
print('Averaging temp data took:', end - start, 's')

############################################################################### 2009

# Uni Laptop
dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2009_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
#dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2009_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

d_cf = dataset.variables['cc'][:] #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
d_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
d_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
d_temp = dataset.variables['t'][:] #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')

#---------------cf-------------#
start = time.time()

d_cf = np.array(d_cf)
d_cf = np.mean(d_cf, axis=0) #Average fraction of cloud cover over time
d_cf = np.mean(d_cf, axis=-1) #Average fraction of cloud cover over longitude
d_cf_alt_lat = d_cf

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
d_cf_so = np.vstack((lat, d_cf)).T #creates a (721,38) array

#Select latitudes over the southern ocean
d_cf_so = d_cf_so[d_cf_so[:,0]>=-70]
d_cf_so = d_cf_so[d_cf_so[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
d_cf_so = d_cf_so[:,1:38]

d_cf = np.mean(d_cf, axis=-1) #Average fraction of cloud cover over latitude
d_cf_so = np.mean(d_cf_so, axis=0) #Average southern ocean fraction of cloud cover over latitude

end = time.time()
print('Averaging cf data took:', end - start, 's')

#---------------iw-------------#

start = time.time()

d_iw = np.array(d_iw)
d_iw = np.mean(d_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
d_iw = np.mean(d_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude
d_iw_alt_lat = d_iw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
d_iw_so = np.vstack((lat, d_iw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
d_iw_so = d_iw_so[d_iw_so[:,0]>=-70]
d_iw_so = d_iw_so[d_iw_so[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
d_iw_so = d_iw_so[:,1:38]

d_iw = np.mean(d_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over latitude
d_iw_so = np.mean(d_iw_so, axis=0) #Average southern ocean specific cloud ice water content (kg/kg) over latitude

end = time.time()
print('Averaging iw data took:', end - start, 's')

#---------------lw-------------#

start = time.time()

d_lw = np.array(d_lw)
d_lw = np.mean(d_lw, axis=0) #Average specific cloud liquid water content (kg/kg) over time
d_lw = np.mean(d_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over longitude
d_lw_alt_lat = d_lw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
d_lw_so = np.vstack((lat, d_lw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
d_lw_so = d_lw_so[d_lw_so[:,0]>=-70]
d_lw_so = d_lw_so[d_lw_so[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
d_lw_so = d_lw_so[:,1:38]

d_lw = np.mean(d_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over latitude
d_lw_so = np.mean(d_lw_so, axis=0) #Average southern ocean specific cloud liquid water content (kg/kg) over latitude

end = time.time()
print('Averaging lw data took:', end - start, 's')

#---------------temperature-------------#

start = time.time()

d_temp = np.array(d_temp)
d_temp = np.mean(d_temp, axis=0) #Average air temperature over time
d_temp = np.mean(d_temp, axis=-1) #Average air temperature over longitude
d_temp_alt_lat = d_temp

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
d_temp_so = np.vstack((lat, d_temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
d_temp_so = d_temp_so[d_temp_so[:,0]>=-70]
d_temp_so = d_temp_so[d_temp_so[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
d_temp_so = d_temp_so[:,1:38]

d_temp = np.mean(d_temp, axis=-1) #Average air temperature over latitude
d_temp_so = np.mean(d_temp_so, axis=0)  #Average southern ocean air temperature over latitude

end = time.time()
print('Averaging temp data took:', end - start, 's')

############################################################################### 2010

# Uni Laptop
dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2010_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
#dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2010_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

e_cf = dataset.variables['cc'][:] #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
e_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
e_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
e_temp = dataset.variables['t'][:] #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')

#---------------cf-------------#
start = time.time()

e_cf = np.array(e_cf)
e_cf = np.mean(e_cf, axis=0) #Average fraction of cloud cover over time
e_cf = np.mean(e_cf, axis=-1) #Average fraction of cloud cover over longitude
e_cf_alt_lat = e_cf

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
e_cf_so = np.vstack((lat, e_cf)).T #creates a (721,38) array

#Select latitudes over the southern ocean
e_cf_so = e_cf_so[e_cf_so[:,0]>=-70]
e_cf_so = e_cf_so[e_cf_so[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
e_cf_so = e_cf_so[:,1:38]

e_cf = np.mean(e_cf, axis=-1) #Average fraction of cloud cover over latitude
e_cf_so = np.mean(e_cf_so, axis=0) #Average southern ocean fraction of cloud cover over latitude

end = time.time()
print('Averaging cf data took:', end - start, 's')

#---------------iw-------------#

start = time.time()

e_iw = np.array(e_iw)
e_iw = np.mean(e_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
e_iw = np.mean(e_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude
e_iw_alt_lat = e_iw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
e_iw_so = np.vstack((lat, e_iw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
e_iw_so = e_iw_so[e_iw_so[:,0]>=-70]
e_iw_so = e_iw_so[e_iw_so[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
e_iw_so = e_iw_so[:,1:38]

e_iw = np.mean(e_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over latitude
e_iw_so = np.mean(e_iw_so, axis=0) #Average southern ocean specific cloud ice water content (kg/kg) over latitude

end = time.time()
print('Averaging iw data took:', end - start, 's')

#---------------lw-------------#

start = time.time()

e_lw = np.array(e_lw)
e_lw = np.mean(e_lw, axis=0) #Average specific cloud liquid water content (kg/kg) over time
e_lw = np.mean(e_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over longitude
e_lw_alt_lat = e_lw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
e_lw_so = np.vstack((lat, e_lw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
e_lw_so = e_lw_so[e_lw_so[:,0]>=-70]
e_lw_so = e_lw_so[e_lw_so[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
e_lw_so = e_lw_so[:,1:38]

e_lw = np.mean(e_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over latitude
e_lw_so = np.mean(e_lw_so, axis=0) #Average southern ocean specific cloud liquid water content (kg/kg) over latitude

end = time.time()
print('Averaging lw data took:', end - start, 's')

#---------------temperature-------------#

start = time.time()

e_temp = np.array(e_temp)
e_temp = np.mean(e_temp, axis=0) #Average air temperature over time
e_temp = np.mean(e_temp, axis=-1) #Average air temperature over longitude
e_temp_alt_lat = e_temp

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
e_temp_so = np.vstack((lat, e_temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
e_temp_so = e_temp_so[e_temp_so[:,0]>=-70]
e_temp_so = e_temp_so[e_temp_so[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
e_temp_so = e_temp_so[:,1:38]

e_temp = np.mean(e_temp, axis=-1) #Average air temperature over latitude
e_temp_so = np.mean(e_temp_so, axis=0)  #Average southern ocean air temperature over latitude

end = time.time()
print('Averaging temp data took:', end - start, 's')

############################################################################### 2011

# Uni Laptop
dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2011_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
#dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2011_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

f_cf = dataset.variables['cc'][:] #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
f_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
f_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
f_temp = dataset.variables['t'][:] #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')

#---------------cf-------------#
start = time.time()

f_cf = np.array(f_cf)
f_cf = np.take(f_cf, [1,2,3], axis=0) #select months July to December 2006
f_cf = np.mean(f_cf, axis=0) #Average fraction of cloud cover over time
f_cf = np.mean(f_cf, axis=-1) #Average fraction of cloud cover over longitude
f_cf_alt_lat = f_cf

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
f_cf_so = np.vstack((lat, f_cf)).T #creates a (721,38) array

#Select latitudes over the southern ocean
f_cf_so = f_cf_so[f_cf_so[:,0]>=-70]
f_cf_so = f_cf_so[f_cf_so[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
f_cf_so = f_cf_so[:,1:38]

f_cf = np.mean(f_cf, axis=-1) #Average fraction of cloud cover over latitude
f_cf_so = np.mean(f_cf_so, axis=0) #Average southern ocean fraction of cloud cover over latitude

end = time.time()
print('Averaging cf data took:', end - start, 's')

#---------------iw-------------#

start = time.time()

f_iw = np.array(f_iw)
f_iw = np.take(f_iw, [1,2,3], axis=0) #select months July to December 2006
f_iw = np.mean(f_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
f_iw = np.mean(f_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude
f_iw_alt_lat = f_iw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
f_iw_so = np.vstack((lat, f_iw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
f_iw_so = f_iw_so[f_iw_so[:,0]>=-70]
f_iw_so = f_iw_so[f_iw_so[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
f_iw_so = f_iw_so[:,1:38]

f_iw = np.mean(f_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over latitude
f_iw_so = np.mean(f_iw_so, axis=0) #Average southern ocean specific cloud ice water content (kg/kg) over latitude

end = time.time()
print('Averaging iw data took:', end - start, 's')

#---------------lw-------------#

start = time.time()

f_lw = np.array(f_lw)
f_lw = np.take(f_lw, [1,2,3], axis=0) #select months July to December 2006
f_lw = np.mean(f_lw, axis=0) #Average specific cloud liquid water content (kg/kg) over time
f_lw = np.mean(f_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over longitude
f_lw_alt_lat = f_lw

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
f_lw_so = np.vstack((lat, f_lw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
f_lw_so = f_lw_so[f_lw_so[:,0]>=-70]
f_lw_so = f_lw_so[f_lw_so[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
f_lw_so = f_lw_so[:,1:38]

f_lw = np.mean(f_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over latitude
f_lw_so = np.mean(f_lw_so, axis=0) #Average southern ocean specific cloud liquid water content (kg/kg) over latitude

end = time.time()
print('Averaging lw data took:', end - start, 's')

#---------------temperature-------------#

start = time.time()

f_temp = np.array(f_temp)
f_temp = np.take(f_temp, [1,2,3], axis=0) #select months Feb to April 2011
f_temp = np.mean(f_temp, axis=0) #Average air temperature over time
f_temp = np.mean(f_temp, axis=-1) #Average air temperature over longitude
f_temp_alt_lat = f_temp

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
f_temp_so = np.vstack((lat, f_temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
f_temp_so = f_temp_so[f_temp_so[:,0]>=-70]
f_temp_so = f_temp_so[f_temp_so[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
f_temp_so = f_temp_so[:,1:38]

f_temp = np.mean(f_temp, axis=-1) #Average air temperature over latitude
f_temp_so = np.mean(f_temp_so, axis=0)  #Average southern ocean air temperature over latitude

end = time.time()
print('Averaging temp data took:', end - start, 's')



###############################################################################
###############################################################################
###############################################################################

#---combine each dataset with common plevel---#

a_cf = np.vstack((plevel, a_cf)).T 
a_lw = np.vstack((plevel, a_lw)).T 
a_iw = np.vstack((plevel, a_iw)).T 
a_cf_so = np.vstack((plevel, a_cf_so)).T 
a_lw_so = np.vstack((plevel, a_lw_so)).T 
a_iw_so = np.vstack((plevel, a_iw_so)).T 
a_temp = np.vstack((plevel, a_temp)).T 
a_temp_so = np.vstack((plevel, a_temp_so)).T 

b_cf = np.vstack((plevel, b_cf)).T 
b_lw = np.vstack((plevel, b_lw)).T 
b_iw = np.vstack((plevel, b_iw)).T 
b_cf_so = np.vstack((plevel, b_cf_so)).T 
b_lw_so = np.vstack((plevel, b_lw_so)).T 
b_iw_so = np.vstack((plevel, b_iw_so)).T 
b_temp = np.vstack((plevel, b_temp)).T 
b_temp_so = np.vstack((plevel, b_temp_so)).T 

c_cf = np.vstack((plevel, c_cf)).T 
c_lw = np.vstack((plevel, c_lw)).T 
c_iw = np.vstack((plevel, c_iw)).T 
c_cf_so = np.vstack((plevel, c_cf_so)).T 
c_lw_so = np.vstack((plevel, c_lw_so)).T 
c_iw_so = np.vstack((plevel, c_iw_so)).T 
c_temp = np.vstack((plevel, c_temp)).T 
c_temp_so = np.vstack((plevel, c_temp_so)).T 

d_cf = np.vstack((plevel, d_cf)).T 
d_lw = np.vstack((plevel, d_lw)).T 
d_iw = np.vstack((plevel, d_iw)).T 
d_cf_so = np.vstack((plevel, d_cf_so)).T 
d_lw_so = np.vstack((plevel, d_lw_so)).T 
d_iw_so = np.vstack((plevel, d_iw_so)).T 
d_temp = np.vstack((plevel, d_temp)).T 
d_temp_so = np.vstack((plevel, d_temp_so)).T 

e_cf = np.vstack((plevel, e_cf)).T 
e_lw = np.vstack((plevel, e_lw)).T 
e_iw = np.vstack((plevel, e_iw)).T 
e_cf_so = np.vstack((plevel, e_cf_so)).T 
e_lw_so = np.vstack((plevel, e_lw_so)).T 
e_iw_so = np.vstack((plevel, e_iw_so)).T 
e_temp = np.vstack((plevel, e_temp)).T 
e_temp_so = np.vstack((plevel, e_temp_so)).T 

f_cf = np.vstack((plevel, f_cf)).T 
f_lw = np.vstack((plevel, f_lw)).T 
f_iw = np.vstack((plevel, f_iw)).T 
f_cf_so = np.vstack((plevel, f_cf_so)).T 
f_lw_so = np.vstack((plevel, f_lw_so)).T 
f_iw_so = np.vstack((plevel, f_iw_so)).T 
f_temp = np.vstack((plevel, f_temp)).T 
f_temp_so = np.vstack((plevel, f_temp_so)).T 

###############################################################################
###############################################################################
###############################################################################

#---reduce each dataset based on common plevel---#
############################################################################### cf

combined = np.vstack((a_cf, b_cf, c_cf, d_cf, e_cf, f_cf))

#print("get unique alt")
unique = np.unique(a_cf[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
cf = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        cf[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
cf[i] = [current_alt, average]
cf = cf[:,1]

############################################################################### cf_so

combined = np.vstack((a_cf_so, b_cf_so, c_cf_so, d_cf_so, e_cf_so, f_cf_so))

#print("get unique alt")
unique = np.unique(a_cf_so[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
cf_so = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        cf_so[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
cf_so[i] = [current_alt, average]
cf_so = cf_so[:,1]

############################################################################### lw

combined = np.vstack((a_lw, b_lw, c_lw, d_lw, e_lw, f_lw))

#print("get unique alt")
unique = np.unique(a_lw[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
lw = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        lw[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
lw[i] = [current_alt, average]
lw = lw[:,1]


############################################################################### lw_so

combined = np.vstack((a_lw_so, b_lw_so, c_lw_so, d_lw_so, e_lw_so, f_lw_so))

#print("get unique alt")
unique = np.unique(a_lw_so[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
lw_so = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        lw_so[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
lw_so[i] = [current_alt, average]
lw_so = lw_so[:,1]


############################################################################### iw

combined = np.vstack((a_iw, b_iw, c_iw, d_iw, e_iw, f_iw))

#print("get unique alt")
unique = np.unique(a_iw[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
iw = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        iw[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
iw[i] = [current_alt, average]
iw = iw[:,1]

############################################################################### iw_so

combined = np.vstack((a_iw_so, b_iw_so, c_iw_so, d_iw_so, e_iw_so, f_iw_so))

#print("get unique alt")
unique = np.unique(a_iw_so[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
iw_so = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        iw_so[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
iw_so[i] = [current_alt, average]
iw_so = iw_so[:,1]

############################################################################### temp

combined = np.vstack((a_temp, b_temp, c_temp, d_temp, e_temp, f_temp))

#print("get unique alt")
unique = np.unique(a_temp[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
temp = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        temp[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
temp[i] = [current_alt, average]
temp = temp[:,1]


############################################################################### temp_so

combined = np.vstack((a_temp_so, b_temp_so, c_temp_so, d_temp_so, e_temp_so, f_temp_so))

#print("get unique alt")
unique = np.unique(a_temp_so[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
temp_so = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        temp_so[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
temp_so[i] = [current_alt, average]
temp_so = temp_so[:,1]

###############################################################################
###############################################################################
###############################################################################

#---Calculate Altitude---#

#https://www.mide.com/pages/air-pressure-at-altitude-calculator
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
alt_t = np.empty((plevel.size,1),dtype=float)
alt_t_so = np.empty((plevel.size,1),dtype=float)
alt_p = np.empty((plevel.size,1),dtype=float)
alt_ts = np.empty((plevel.size,1),dtype=float)
alt_ts_so = np.empty((plevel.size,1),dtype=float)


# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in temp:
    newalt = (288.15 - item)/6.49
    alt_t[i] = [newalt]
    i+=1

#alt_t[i] = [newalt]

# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in temp_so:
    newalt = (288.15 - item)/6.49
    alt_t_so[i] = [newalt]
    i+=1

#alt_t_so[i] = [newalt]

# Iterate through all of the pressure elements (lower stratosphere 11km < h <25km)
i = 0
for item in plevel:
    newalt = (1.73 - math.log(item/226.5))/0.157
    alt_p[i] = [newalt]
    i+=1

#alt_p[i] = [newalt]


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in temp:
    newalt = (item-141.94)/2.99
    alt_ts[i] = [newalt]
    i+=1

#alt_ts[i] = [newalt]

# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in temp_so:
    newalt = (item-141.94)/2.99
    alt_ts_so[i] = [newalt]
    i+=1

#alt_ts_so[i] = [newalt]

alt = alt_t
alt_so = alt_t_so

"""
#---------------combine-------------#

alt_so = np.transpose(alt_so)
alt = np.transpose(alt)
# Join the two lists as if they were two columns side by side, into a list of two elements each
ecmwf_tcc_alt = np.vstack((alt, cf)).T
ecmwf_tclw_alt = np.vstack((alt, lw)).T
ecmwf_tciw_alt = np.vstack((alt, iw)).T
ecmwf_temp_alt = np.vstack((alt, temp)).T
ecmwf_plevel_alt = np.vstack((alt, plevel)).T

ecmwf_tcc_temp = np.vstack((temp, cf)).T 
ecmwf_tclw_temp = np.vstack((temp, lw)).T
ecmwf_tciw_temp = np.vstack((temp, iw)).T

ecmwf_tcc_alt_so = np.vstack((alt_so, cf_so)).T 
ecmwf_tclw_alt_so = np.vstack((alt_so, lw_so)).T
ecmwf_tciw_alt_so = np.vstack((alt_so, iw_so)).T
ecmwf_temp_alt_so = np.vstack((alt_so, temp_so)).T
ecmwf_plevel_alt_so = np.vstack((alt_so, plevel)).T

ecmwf_tcc_temp_so = np.vstack((temp_so, cf_so)).T 
ecmwf_tclw_temp_so = np.vstack((temp_so, lw_so)).T
ecmwf_tciw_temp_so = np.vstack((temp_so, iw_so)).T

#----------------------------#


###############################################################################
import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') #Uni Laptop

with h5py.File('07.2006_04.2011_ECMWF.h5', 'w') as p:
    p.create_dataset('lat', data=lat)
    p.create_dataset('alt', data=alt)
    
    p.create_dataset('tcc', data=ecmwf_tcc_lat)
    p.create_dataset('tclw', data=ecmwf_tclw_lat)
    p.create_dataset('tciw', data=ecmwf_tciw_lat)  
    
    p.create_dataset('cf', data=ecmwf_tcc_alt)
    p.create_dataset('cf_so', data=ecmwf_tcc_alt_so)
    p.create_dataset('lw', data=ecmwf_tclw_alt)
    p.create_dataset('lw_so', data=ecmwf_tclw_alt_so)
    p.create_dataset('iw', data=ecmwf_tciw_alt)
    p.create_dataset('iw_so', data=ecmwf_tciw_alt_so)
    
    p.create_dataset('temp_alt_lat', data=temp_alt_lat)  
    p.create_dataset('cf_alt_lat', data=cf_alt_lat)  
    p.create_dataset('lw_alt_lat', data=lw_alt_lat)
    p.create_dataset('iw_alt_lat', data=iw_alt_lat)


    p.create_dataset('temp', data=ecmwf_temp_alt)
    p.create_dataset('temp_so', data=ecmwf_temp_alt_so)
    p.create_dataset('pressure', data=ecmwf_plevel_alt)
    p.create_dataset('pressure_so', data=ecmwf_plevel_alt_so)

    p.create_dataset('cf_t', data=ecmwf_tcc_temp)
    p.create_dataset('cf_t_so', data=ecmwf_tcc_temp_so)
    p.create_dataset('lw_t', data=ecmwf_tclw_temp)
    p.create_dataset('lw_t_so', data=ecmwf_tclw_temp_so)
    p.create_dataset('iw_t', data=ecmwf_tciw_temp)
    p.create_dataset('iw_t_so', data=ecmwf_tciw_temp_so)

    p.close()
"""