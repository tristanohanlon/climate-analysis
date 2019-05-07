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
import math as math
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import h5py
import os
from scipy import integrate

start = time.time()
os.chdir('E:/University/University/MSc/Models/Data/GFDL/') #Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/GFDL/') #Uni Laptop

dataset = Dataset('cl_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
cf = dataset.variables['cl'][:] #Extract cloud_area_fraction_in_atmosphere_layer, keyed to time, lev, lon and lat (420, 33, 180, 288)
lat = dataset.variables['lat'][:] #Extract latitude, (180,)

dataset = Dataset('clw_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
lw = dataset.variables['clw'][:] #Extract Mass Fraction of Cloud Liquid Water (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('cli_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
iw = dataset.variables['cli'][:] #Extract Mass Fraction of Cloud Ice (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('ta_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
temp = dataset.variables['ta'][:] #Extract air temperature, keyed to time, plev, lon and lat (420, 19, 180, 288)
plev = dataset.variables['plev'][:] #Extract air pressure corresponding to temperature (19,)

dataset = Dataset('pfull_AERmon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
pressure = dataset.variables['pfull'][:] #Extract air pressure, keyed to time, level, lon and lat (420, 33, 180, 288)

end = time.time()
print('Importing data took:', end - start, 's')

#----------------------------#

#Select the months from 2010
cf = np.take(cf, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
lw = np.take(lw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
iw = np.take(iw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
temp = np.take(temp, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
pressure = np.take(pressure, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
#Average over time
cf = np.mean(cf, axis=0)
lw = np.mean(lw, axis=0)
iw = np.mean(iw, axis=0)
temp = np.mean(temp, axis=0)
pressure = np.mean(pressure, axis=0)
#Average over longitude
cf = np.mean(cf, axis=-1)
lw = np.mean(lw, axis=-1)
iw = np.mean(iw, axis=-1)
temp = np.mean(temp, axis=-1)
pressure = np.mean(pressure, axis=-1)
#Select Southern ocean Latitudes

tcf = np.vstack((lat, cf)).T #creates a (180,34) array
tcf = tcf[tcf[:,0]>=-70]
tcf = tcf[tcf[:,0]<=-50]
cf = tcf[:,1:34] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude

tclw = np.vstack((lat, lw)).T #creates a (180,34) array
tclw = tclw[tclw[:,0]>=-70]
tclw = tclw[tclw[:,0]<=-50]
lw = tclw[:,1:34] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude

tciw = np.vstack((lat, iw)).T #creates a (180,34) array
tciw = tciw[tciw[:,0]>=-70]
tciw = tciw[tciw[:,0]<=-50]
iw = tciw[:,1:34] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude

tcp = np.vstack((lat, pressure)).T #creates a (180,34) array
tcp = tcp[tcp[:,0]>=-70]
tcp = tcp[tcp[:,0]<=-50]
pressure = tcp[:,1:34] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude

#Average over latitude
cf = np.mean(cf, axis=0) / 100
lw = np.mean(lw, axis=0)
iw = np.mean(iw, axis=0)
temp = np.mean(temp, axis=-1)
pressure = np.mean(pressure, axis=0) / 100

#---convert pressure levels to altitude---#
i = 0
alt_trop = np.empty(pressure.size)

for value in pressure:
    alt = (288.19 - 288.08*math.pow((value/1012.90), (1/5.256))) / 0.00649  
    alt_trop[i] = alt
    i+=1
alt_trop = alt_trop / 1000    
 

i = 0
alt_strat = np.empty(pressure.size)    

for value in pressure:
    alt = ((1.73 - math.log(value/226.50)) / 0.000157)
    alt_strat[i] = alt
    i+=1
alt_strat = alt_strat / 1000   
    
alt_trop = alt_trop[0:20]       
alt_strat = alt_strat[20:33]  
alt = np.concatenate((alt_trop, alt_strat), axis = 0)
"""
#---get comparison altitude from smaller temp and plev arrays---#
i = 0
alt_comp_trop = np.empty(gplev.size)

for value in gplev:
    altc = (288.19 - 288.08*math.pow((value/101290), (1/5.256))) / 0.00649  
    alt_comp_trop[i] = altc
    i+=1
alt_comp_trop = alt_comp_trop / 1000    
alt_comp_trop = alt_comp_trop[0:9]  

i = 0
alt_comp_strat = np.empty(gplev.size)    

for value in gplev:
    altc = ((1.73 - math.log(value/22650)) / 0.000157) / 1000  
    alt_comp_strat[i] = altc
    i+=1
alt_comp_strat=alt_comp_strat[9:15]  
alt_comp = np.concatenate((alt_comp_trop, alt_comp_strat), axis = 0)

#---get mass of air path (kg/m^2)---#

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = (plev) / (286.9 * temp)

ap = integrate.trapz(air_density, alt_comp) #mass of air path (kg/m^2)
"""
#---create temperature profile---#
#troposphere from pressure data
#stratosphere assumed from Earth Atmosphere Model 
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html

temp_profile = 288.19 - 6.49*alt

#manually add on 216.69K values to altitudes above 11km to 25km 

###############################################################################

#---create datasets---#

cf = np.vstack((alt, cf)).T
lw = np.vstack((alt, lw)).T
iw = np.vstack((alt, iw)).T
temp = np.vstack((alt, temp_profile)).T
pressure = np.vstack((alt, pressure)).T
    
cf_t = np.vstack((temp_profile, cf[:,1])).T
lw_t = np.vstack((temp_profile, lw[:,1])).T
iw_t = np.vstack((temp_profile, iw[:,1])).T
    
cf_p = np.vstack((pressure[:,1], cf[:,1])).T
lw_p = np.vstack((pressure[:,1], lw[:,1])).T
iw_p = np.vstack((pressure[:,1], iw[:,1])).T


plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(cf[:,1],cf[:,0], '-r', label='Fraction Cloud Cover')
ax2.plot(lw[:,1],lw[:,0], '-b', label='Specific Cloud Liquid Water Content')
ax2.plot(iw[:,1],iw[:,0], '--b', label='Specific Cloud Ice Water Content')

#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
          ncol=4, fancybox=True, shadow=True);
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Cloud Liquid and Ice Water Content (kg/kg) x $10^{-4}$')
ax1.set_ylabel('Altitude (km)')
plt.title('Cloud Fraction and Phase vs Altitude GFDL.AM4 2010')
#plt.gca().invert_yaxis()

plt.grid(True)
plt.show()



















