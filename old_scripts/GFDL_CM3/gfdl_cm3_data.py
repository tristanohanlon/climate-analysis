# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon

"""
import time
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import h5py
import math
from scipy import integrate



########################################---get variables---########################################

os.chdir('E:/University/University/MSc/Models/Data/CMIP5/gfdl_cm3_rcp4.5') #Home PC
#os.chdir('D:/MSc/Models/Data/CMIP5/gfdl_cm3_rcp4.5') #ext HDD

f = Dataset('cl_Amon_GFDL-CM3_rcp45_r1i1p1_200601-201012.nc', 'r')

################---get latitude---################
lat = np.array(f.variables['lat'][:])

################---get pressure levels---################
plev = np.array(f.variables['lev'][:]) #in hPa
a = np.array(f.variables['a'][:]) #in hPa
b = np.array(f.variables['b'][:]) #in hPa
p0 = np.array(f.variables['p0'][:]) #in hPa
ps = np.array(f.variables['ps'][:]) #in hPa
ps = np.mean(ps, axis = 0)
ps = np.mean(ps, axis = 0)
ps = np.mean(ps, axis = 0)

p = a*p0 + b*ps

p = np.array(p)

################---get cloud fraction---################
cf1 = np.array(f.variables['cl'][:])
cf1 = cf1[6:] # get values from 07.2006 to 12.2010

f = Dataset('cl_Amon_GFDL-CM3_rcp45_r1i1p1_201101-201512.nc', 'r')
cf2 = np.array(f.variables['cl'][:])
cf2 = cf2[1:4] # get values from 02.2011 to 04.2011

cf = np.concatenate((cf1,cf2))
cf = np.mean(cf, axis = 0) / 100


################---get total cloud cover---################
f = Dataset('clt_Amon_GFDL-CM3_rcp45_r1i1p1_200601-201012.nc', 'r')
tcc1 = np.array(f.variables['clt'][:])
tcc1 = tcc1[6:] # get values from 07.2006 to 12.2010

f = Dataset('clt_Amon_GFDL-CM3_rcp45_r1i1p1_201101-201512.nc', 'r')
tcc2 = np.array(f.variables['clt'][:])
tcc2 = tcc2[1:4] # get values from 02.2011 to 04.2011

tcc = np.concatenate((tcc1,tcc2))
tcc = np.mean(tcc, axis = 0)
tcc = np.mean(tcc, axis = -1) / 100


################---get cloud liquid content################
f = Dataset('clw_Amon_GFDL-CM3_rcp45_r1i1p1_200601-201012.nc', 'r')
lw1 = np.array(f.variables['clw'][:])
lw1 = lw1[6:] # get values from 07.2006 to 12.2010

f = Dataset('clw_Amon_GFDL-CM3_rcp45_r1i1p1_201101-201512.nc', 'r')
lw2 = np.array(f.variables['clw'][:])
lw2 = lw2[1:4] # get values from 02.2011 to 04.2011

lw = np.concatenate((lw1,lw2))
lw = np.mean(lw, axis = 0)



################---get cloud ice content################
f = Dataset('cli_Amon_GFDL-CM3_rcp45_r1i1p1_200601-201012.nc', 'r')
iw1 = np.array(f.variables['cli'][:])
iw1 = iw1[6:] # get values from 07.2006 to 12.2010

f = Dataset('cli_Amon_GFDL-CM3_rcp45_r1i1p1_201101-201512.nc', 'r')
iw2 = np.array(f.variables['cli'][:])
iw2 = iw2[1:4] # get values from 02.2011 to 04.2011

iw = np.concatenate((iw1,iw2))
iw = np.mean(iw, axis = 0)


################---get temperature################
f = Dataset('ta_Amon_GFDL-CM3_rcp45_r1i1p1_200601-201012.nc', 'r')
T1 = np.array(f.variables['ta'][:])
T1 = T1[6:] # get values from 07.2006 to 12.2010
p_temp = np.array(f.variables['plev'][:])

f = Dataset('ta_Amon_GFDL-CM3_rcp45_r1i1p1_201101-201512.nc', 'r')
T2 = np.array(f.variables['ta'][:])
T2 = T2[1:4] # get values from 02.2011 to 04.2011

T = np.concatenate((T1,T2))
T = np.mean(T, axis = 0)


###############################################################################

#---convert pressure levels to altitude---#

#https://www.mide.com/pages/air-pressure-at-altitude-calculator
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
alt_t = np.empty((plev.size,1),dtype=float)
alt_p = np.empty((plev.size,1),dtype=float)
alt_ts = np.empty((plev.size,1),dtype=float)

# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in p:
    newalt = (288.19 - 288.08*((item/101290)**(1/5.256)))/6.49
    alt_t[i] = [newalt]
    i+=1


# Iterate through all of the pressure elements (lower stratosphere 11km < h <25km)
i = 0
for item in p:
    newalt = (1.73 - math.log(item/22650))/0.157
    alt_p[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in p:
    newalt = (216.6*((item/2488)**(1/-11.388)) - 141.94)/2.99
    alt_ts[i] = [newalt]
    i+=1


#manually adjust alt and alt_so arrays usinf alt_p and alt_ts
alt = alt_t

#Need to get extra temperature levels based on pressure values

#---create temperature profile---#
#troposphere from pressure data
#stratosphere assumed from Earth Atmosphere Model 
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html

temp_t = 288.19 - 6.49*alt
temp_s = 141.94 + 2.99*alt
#manually add on 216.69K values to altitudes above 11km to 25km 

temp = temp_t
temp_so = temp_t
"""
os.chdir('E:/University/University/MSc/Models/climate-analysis/CAM5/reduced_datasets') #Home PC
with h5py.File('CAM5_alt.h5', 'w') as p:
    
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('lat', data=lat)
    p.create_dataset('plev', data=plev)
    p.create_dataset('lw', data=lw)
    p.create_dataset('iw', data=iw)
    p.create_dataset('temp', data=temp)
    p.create_dataset('alt', data=alt)
    
 
    p.close()
"""
############################################################################### Get Temperature altitudes

#---convert pressure levels to altitude---#

#https://www.mide.com/pages/air-pressure-at-altitude-calculator
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
temp_alt_t = np.empty((p_temp.size,1),dtype=float)
temp_alt_p = np.empty((p_temp.size,1),dtype=float)
temp_alt_ts = np.empty((p_temp.size,1),dtype=float)

# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in p_temp:
    newalt = (288.19 - 288.08*((item/101290)**(1/5.256)))/6.49
    temp_alt_t[i] = [newalt]
    i+=1


# Iterate through all of the pressure elements (lower stratosphere 11km < h <25km)
i = 0
for item in p_temp:
    newalt = (1.73 - math.log(item/22650))/0.157
    temp_alt_p[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in p_temp:
    newalt = (216.6*((item/2488)**(1/-11.388)) - 141.94)/2.99
    temp_alt_ts[i] = [newalt]
    i+=1


#manually adjust alt and alt_so arrays usinf alt_p and alt_ts
alt_temp = temp_alt_t

#Need to get extra temperature levels based on pressure values

#---create temperature profile---#
#troposphere from pressure data
#stratosphere assumed from Earth Atmosphere Model 
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html

temp_t = 288.19 - 6.49*alt
temp_s = 141.94 + 2.99*alt
#manually add on 216.69K values to altitudes above 11km to 25km 

temp = temp_t
###############################################################################

#---combine arrays---#
# since lw and iw are in kg/kg - need to convert to LWP and IWC in kgm^-2
# Get density levels
# Integrate density with altitude to get air path AP
# multiply lw and iw by AP

alt_m = np.hstack(alt*1000)
p = np.vstack(p/100)
pressure = np.hstack((alt, p))
temp_g = np.hstack((alt, temp))

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = (pressure[:,1] * 100) / (286.9 * temp_g[:,1])

ap = integrate.trapz(air_density, alt_m)

tclw = np.mean(lw , axis = 0)
tclw = np.mean(tclw  , axis = -1) * ap

tciw = np.mean(iw , axis = 0)
tciw = np.mean(tciw , axis = -1) * ap

tcc = np.vstack((lat, tcc)).T # Join the two lists as if they were two columns side by side, into a list of two elements each
tclw = np.vstack((lat, tclw)).T
tciw = np.vstack((lat, tciw)).T

#----------------------------#
alt = np.hstack(alt)
alt_temp = np.hstack(alt_temp)

cf_g = np.mean(cf, axis = -1)
cf_g = np.mean(cf_g, axis = -1)
cf_g = np.vstack((alt, cf_g)).T

lw_g = np.mean(lw, axis = -1)
lw_g = np.mean(lw_g, axis = -1)
lw_g = np.vstack((alt, lw_g)).T

iw_g = np.mean(iw, axis = -1)
iw_g = np.mean(iw_g, axis = -1)
iw_g = np.vstack((alt, iw_g)).T

temp_alt_lat = np.mean(T, axis = -1) # altitude levels correspond with alt_temp
cf_alt_lat = np.mean(cf, axis = -1)
lw_alt_lat = np.mean(lw, axis = -1)
iw_alt_lat = np.mean(iw, axis = -1)

#Select Southern ocean Latitudes
cf_so = np.mean(cf, axis = -1)
cf_so = np.transpose(cf_so)

cf_so = np.hstack((np.vstack(lat), cf_so)) #creates a (180,34) array
cf_so = cf_so[cf_so[:,0]>=-70]
cf_so = cf_so[cf_so[:,0]<=-50]
cf_so = cf_so[:,1:] #Split the combined array into just the tccf data, eliminating the first coloumn of latitude
cf_so = np.mean(cf_so, axis = 0)
cf_so = np.vstack((alt, cf_so)).T

lw_so = np.mean(lw, axis = -1)
lw_so = np.transpose(lw_so)

lw_so = np.hstack((np.vstack(lat), lw_so)) #creates a (180,34) array
lw_so = lw_so[lw_so[:,0]>=-70]
lw_so = lw_so[lw_so[:,0]<=-50]
lw_so = lw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
lw_so = np.mean(lw_so, axis = 0)
lw_so = np.vstack((alt, lw_so)).T


iw_so = np.mean(iw, axis = -1)
iw_so = np.transpose(iw_so)

iw_so = np.hstack((np.vstack(lat), iw_so)) #creates a (180,34) array
iw_so = iw_so[iw_so[:,0]>=-70]
iw_so = iw_so[iw_so[:,0]<=-50]
iw_so = iw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
iw_so = np.mean(iw_so, axis = 0)
iw_so = np.vstack((alt, iw_so)).T


temp_so = np.mean(T, axis = -1)
temp_so = np.transpose(temp_so)

temp_so = np.hstack((np.vstack(lat), temp_so)) #creates a (180,34) array
temp_so = temp_so[temp_so[:,0]>=-70]
temp_so = temp_so[temp_so[:,0]<=-50]
temp_so = temp_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
temp_so = np.mean(temp_so, axis = 0)
temp_so = np.vstack((alt_temp, temp_so)).T


cf_t = np.vstack((temp_g[:,1], cf_g[:,1])).T
cf_t_so = np.vstack((temp_g[:,1], cf_so[:,1])).T
lw_t = np.vstack((temp_g[:,1], lw_g[:,1])).T
lw_t_so = np.vstack((temp_g[:,1], lw_so[:,1])).T
iw_t = np.vstack((temp_g[:,1], iw_g[:,1])).T
iw_t_so = np.vstack((temp_g[:,1], iw_so[:,1])).T




os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL_CM3/reduced_datasets') #Home PC
with h5py.File('07.2006_04.2011_GFDL_CM3.h5', 'w') as p:
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_temp', data=alt_temp)  
    p.create_dataset('lat', data=lat)  
    p.create_dataset('air_density', data=air_density)  
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    p.create_dataset('cf', data=cf_g)
    p.create_dataset('lw', data=lw_g)
    p.create_dataset('iw', data=iw_g)
    p.create_dataset('temp', data=temp_g)
    p.create_dataset('cf_so', data=cf_so)
    p.create_dataset('lw_so', data=lw_so)
    p.create_dataset('iw_so', data=iw_so)
    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure)
    p.create_dataset('temp_alt_lat', data=temp_alt_lat)  
    p.create_dataset('cf_alt_lat', data=cf_alt_lat)  
    p.create_dataset('lw_alt_lat', data=lw_alt_lat)
    p.create_dataset('iw_alt_lat', data=iw_alt_lat)
    p.create_dataset('cf_t', data=cf_t)
    p.create_dataset('cf_t_so', data=cf_t_so)
    p.create_dataset('lw_t', data=lw_t)
    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('iw_t', data=iw_t)
    p.create_dataset('iw_t_so', data=iw_t_so)
    p.close()


plt.figure()
fig, ax1 = plt.subplots()

#ax2 = ax1.twinx()

ax1.contourf(lat, alt ,temp_alt_lat)
#ax2.plot(tclw[:,0],tclw[:,1], '-b', label='Liquid Water Content')
#ax2.plot(tciw[:,0],tciw[:,1], '--b', label='Ice Water Content')

#ax.axis('equal')
#ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
#          ncol=4, fancybox=True, shadow=True);
#ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
 #         ncol=4, fancybox=True, shadow=True);
           
ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
#ax2.set_ylabel('Liquid and Ice Water Content ($kgm^{-2}$)')

plt.title('Cloud Fraction and Phase Content vs Latitude - GFDL.AM4 - July 2006 to April 2011')

plt.grid(True)
plt.show()
