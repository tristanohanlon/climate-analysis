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
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2006_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2006_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

lat = dataset.variables['latitude'][:] #Extract latitude data (721)
plevel = dataset.variables['level'][:] #Extract pressure (millibars) level data (37 levels)
a_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
a_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')


#---------------iw-------------#

start = time.time()

a_iw = np.array(a_iw)
a_iw = np.take(a_iw, [6,7,8,9,10,11], axis=0) #select months July to December 2006
a_iw = np.mean(a_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
a_iw = np.mean(a_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude

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


############################################################################### 2007

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2007_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2007_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

b_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
b_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')


#---------------iw-------------#

start = time.time()

b_iw = np.array(b_iw)
b_iw = np.mean(b_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
b_iw = np.mean(b_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude

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


############################################################################### 2008

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2008_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2008_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

c_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
c_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')


#---------------iw-------------#

start = time.time()

c_iw = np.array(c_iw)
c_iw = np.mean(c_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
c_iw = np.mean(c_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude

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


############################################################################### 2009

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2009_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2009_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

d_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
d_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')


#---------------iw-------------#

start = time.time()

d_iw = np.array(d_iw)
d_iw = np.mean(d_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
d_iw = np.mean(d_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude

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



############################################################################### 2010

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2010_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2010_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

e_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
e_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')


#---------------iw-------------#

start = time.time()

e_iw = np.array(e_iw)
e_iw = np.mean(e_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
e_iw = np.mean(e_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude

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


############################################################################### 2011

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2011_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2011_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

f_iw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
f_lw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')


#---------------iw-------------#

start = time.time()

f_iw = np.array(f_iw)
f_iw = np.take(f_iw, [1,2,3], axis=0) #select months July to December 2006
f_iw = np.mean(f_iw, axis=0) #Average specific cloud ice water content (kg/kg) over time
f_iw = np.mean(f_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude

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



###############################################################################




###############################################################################
import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') #Uni Laptop

with h5py.File('07.2006_04.2011_ECMWF.h5', 'w') as p:
    
    p.create_dataset('tcc', data=ecmwf_tcc_lat)
    p.create_dataset('tclw', data=ecmwf_tclw_lat)
    p.create_dataset('tciw', data=ecmwf_tciw_lat)  
    
    p.create_dataset('cf', data=ecmwf_tcc_alt)
    p.create_dataset('cf_so', data=ecmwf_tcc_alt_so)
    p.create_dataset('lw', data=ecmwf_tclw_alt)
    p.create_dataset('lw_so', data=ecmwf_tclw_alt_so)
    p.create_dataset('iw', data=ecmwf_tciw_alt)
    p.create_dataset('iw_so', data=ecmwf_tciw_alt_so)

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