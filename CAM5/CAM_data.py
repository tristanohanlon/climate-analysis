# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon

This will create a dataset of time averaged global total cloud cover with latitude.
Time period is from 01.1980 - 12.2014

"""
import time
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import h5py
import math
from scipy import integrate


start = time.time()

#---get latitude and cf---#
os.chdir('E:/University/University/MSc/Models/Data/CAM5/tcc/') #Home PC
#os.chdir('D:/MSc/Models/Data/CESM1(CAM5)/tcc') #ext HDD

lat = []
cf = np.empty((900, 192, 288),dtype=float)
 # create a blank array to add cloud amount data

# Load every file in the directory
for filename in os.listdir(): 
    # Load the file
    f = Dataset(filename, 'r')
    # Get the cloud cover data as a list. Since this is the 'cloud free' data, need to invert later
    cf = np.concatenate((cf, f.variables['CLDTOT'][:]))
cf = cf[900:3600]
cf1 = cf[182:1826] # remove initial zero values and get values from 07.2006 to 12.2010
cf2 = cf[1857:1946] # get values from 02.2011 to 04.2011
cf = np.concatenate((cf1,cf2))
lat = f.variables['lat'][:]
lat = np.array(lat)

plev = f.variables['lev'][:] #in hPa

tcc = np.mean(cf, axis = 0)
tcc = np.mean(tcc, axis = -1)

#---get surface pressure---#

os.chdir('E:/University/University/MSc/Models/Data/CAM5/ps') #Home PC
#os.chdir('D:/MSc/Models/Data/CESM1(CAM5)/ps') #ext HDD
ps = np.empty((900, 192, 288),dtype=float)
 # create a blank array to add cloud amount data

# Load every file in the directory
for filename in os.listdir(): 
    # Load the file
    f = Dataset(filename, 'r')
    # Get the cloud cover data as a list. Since this is the 'cloud free' data, need to invert later
    ps = np.concatenate((ps, f.variables['PS'][:]))
ps = ps[900:3600]
ps1 = ps[182:1826] # remove initial zero values and get values from 07.2006 to 12.2010
ps2 = ps[1857:1946] # get values from 02.2011 to 04.2011
ps = np.concatenate((ps1,ps2))

ps = np.mean(ps, axis = 0)


end = time.time()
print('Importing data from files to lists took:', end - start, 's')

###############################################################################

#---convert pressure levels to altitude---#

#https://www.mide.com/pages/air-pressure-at-altitude-calculator
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
alt_t = np.empty((plev.size,1),dtype=float)
alt_p = np.empty((plev.size,1),dtype=float)
alt_ts = np.empty((plev.size,1),dtype=float)

# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in plev:
    newalt = (288.19 - 288.08*((item/1012.9)**(1/5.256)))/6.49
    alt_t[i] = [newalt]
    i+=1


# Iterate through all of the pressure elements (lower stratosphere 11km < h <25km)
i = 0
for item in plev:
    newalt = (1.73 - math.log(item/226.50))/0.157
    alt_p[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in plev:
    newalt = (216.6*((item/24.88)**(1/-11.388)) - 141.94)/2.99
    alt_ts[i] = [newalt]
    i+=1


#manually adjust alt and alt_so arrays usinf alt_p and alt_ts
alt = alt_t


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


###############################################################################


#---get lw---#

os.chdir('E:/University/University/MSc/Models/Data/CAM5/tclw') #Home PC
#os.chdir('D:/MSc/Models/Data/CESM1(CAM5)/tclw') #ext HDD
lw1 = np.empty((900, 30, 192, 288),dtype=float)
lw2 = np.empty((900, 30, 192, 288),dtype=float)
lw3 = np.empty((900, 30, 192, 288),dtype=float)


 # create a blank array to add cloud amount data

# Load every file in the directory
for filename in os.listdir(): 
    # Load the file
    f = Dataset('b.e11.BRCP45C5CNBDRD.f09_g16.003.cam.h0.CLDLIQ.200601-208012.nc', 'r')
    # Get the cloud cover data as a list. Since this is the 'cloud free' data, need to invert later
    lw3 = np.array(f.variables['CLDLIQ'][:])
    
lw = np.concatenate((lw1, lw2, lw3))
lw1 = lw[182:1826] # remove initial zero values and get values from 07.2006 to 12.2010
lw2 = lw[1857:1946] # get values from 02.2011 to 04.2011
lw = np.concatenate((lw1,lw2))

lw = np.mean(lw, axis = 0)

#---get iw---#

os.chdir('E:/University/University/MSc/Models/Data/CAM5/tciw') #Home PC
#os.chdir('D:/MSc/Models/Data/CESM1(CAM5)/tciw') #ext HDD
iw1 = np.empty((900, 30, 192, 288),dtype=float)
iw2 = np.empty((900, 30, 192, 288),dtype=float)
iw3 = np.empty((900, 30, 192, 288),dtype=float)

 # create a blank array to add cloud amount data

# Load every file in the directory
for filename in os.listdir(): 
    # Load the file
    f = Dataset('b.e11.BRCP45C5CNBDRD.f09_g16.003.cam.h0.CLDICE.200601-208012.nc', 'r')
    # Get the cloud cover data as a list. Since this is the 'cloud free' data, need to invert later
    iw3 = np.array(f.variables['CLDICE'][:])
    
iw = np.concatenate((iw1, iw2, iw3))
iw1 = iw[182:1826] # remove initial zero values and get values from 07.2006 to 12.2010
iw2 = iw[1857:1946] # get values from 02.2011 to 04.2011
iw = np.concatenate((iw1,iw2))
iw = np.mean(iw, axis = 0)


#---get temp---#

os.chdir('E:/University/University/MSc/Models/Data/CAM5/T') #Home PC
#os.chdir('D:/MSc/Models/Data/CESM1(CAM5)/ps') #ext HDD
temp1 = np.empty((900, 30, 192, 288),dtype=float)
temp2 = np.empty((900, 30, 192, 288),dtype=float)
temp3 = np.empty((900, 30, 192, 288),dtype=float)

 # create a blank array to add cloud amount data

# Load every file in the directory
for filename in os.listdir(): 
    # Load the file
    f = Dataset('b.e11.BRCP45C5CNBDRD.f09_g16.003.cam.h0.T.200601-208012.nc', 'r')
    # Get the cloud cover data as a list. Since this is the 'cloud free' data, need to invert later
    temp3 = np.array(f.variables['T'][:])
    
temp = np.concatenate((temp1, temp2, temp3))
temp1 = temp[182:1826] # remove initial zero values and get values from 07.2006 to 12.2010
temp2 = temp[1857:1946] # get values from 02.2011 to 04.2011
temp = np.concatenate((temp1,temp2))
temp = np.mean(temp, axis = 0)


       
###############################################################################

#---combine arrays---#
# since lw and iw are in kg/kg - need to convert to LWP and IWC in kgm^-2
# Get density levels
# Integrate density with altitude to get air path AP
# multiply lw and iw by AP

alt_m = alt*1000

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = (pressure[:,1] * 100) / (286.9 * temp_g)

ap = integrate.trapz(air_density, alt_m)

tclw = np.mean(lw , axis = 0)
tclw = np.mean(tclw  , axis = -1) * -ap

tciw = np.mean(iw , axis = 0)
tciw = np.mean(tciw , axis = -1) * -ap


tcc = np.vstack((lat, tcc)).T # Join the two lists as if they were two columns side by side, into a list of two elements each
tclw = np.vstack((lat, tclw)).T
tciw = np.vstack((lat, tciw)).T

pressure = np.vstack((alta, plev)).T
#----------------------------#

lw_g = np.mean(lw, axis = -1)
lw_g = np.mean(lw_g, axis = -1)
lw_g = np.vstack((alta, lw_g)).T

iw_g = np.mean(iw, axis = -1)
iw_g = np.mean(iw_g, axis = -1)
iw_g = np.vstack((alta, iw_g)).T

temp_g = np.mean(temp, axis = -1)
temp_g = np.mean(temp_g, axis = -1)
temp_g = np.vstack((alta, temp_g)).T

lw_alt_lat = np.mean(lw, axis = -1)
iw_alt_lat = np.mean(iw, axis = -1)

#Select Southern ocean Latitudes

lw_so = np.mean(lw, axis = -1)
lw_so = np.transpose(lw_so)

lw_so = np.hstack((lat, lw_so)) #creates a (180,34) array
lw_so = lw_so[lw_so[:,0]>=-70]
lw_so = lw_so[lw_so[:,0]<=-50]
lw_so = lw_so[:,1:31] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
lw_so = np.mean(lw_so, axis = 0)
lw_so = np.vstack((alta, lw_so)).T


iw_so = np.mean(iw, axis = -1)
iw_so = np.transpose(iw_so)

iw_so = np.hstack((lat, iw_so)) #creates a (180,34) array
iw_so = iw_so[iw_so[:,0]>=-70]
iw_so = iw_so[iw_so[:,0]<=-50]
iw_so = iw_so[:,1:31] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
iw_so = np.mean(iw_so, axis = 0)
iw_so = np.vstack((alta, iw_so)).T



temp_so = np.mean(temp, axis = -1)
temp_so = np.transpose(temp_so)

temp_so = np.hstack((lat, temp_so)) #creates a (180,34) array
temp_so = temp_so[temp_so[:,0]>=-70]
temp_so = temp_so[temp_so[:,0]<=-50]
temp_so = temp_so[:,1:31] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
temp_so = np.mean(temp_so, axis = 0)
temp_so = np.vstack((alta, temp_so)).T

lw_t = np.vstack((temp_g[:,1], lw_g[:,1])).T
lw_t_so = np.vstack((temp_so[:,1], lw_so[:,1])).T
iw_t = np.vstack((temp_g[:,1], iw_g[:,1])).T
iw_t_so = np.vstack((temp_so[:,1], iw_so[:,1])).T

os.chdir('E:/University/University/MSc/Models/climate-analysis/CAM5/reduced_datasets') #Home PC
with h5py.File('07.2006_04.2011_CAM5.h5', 'w') as p:
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('lat', data=lat)  
    p.create_dataset('air_density', data=air_density)  
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    p.create_dataset('lw', data=lw_g)
    p.create_dataset('iw', data=iw_g)
    p.create_dataset('temp', data=temp_g)
    p.create_dataset('lw_so', data=lw_so)
    p.create_dataset('iw_so', data=iw_so)
    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure)
    p.create_dataset('lw_alt_lat', data=lw_alt_lat)
    p.create_dataset('iw_alt_lat', data=iw_alt_lat)
    p.create_dataset('lw_t', data=lw_t)
    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('iw_t', data=iw_t)
    p.create_dataset('iw_t_so', data=iw_t_so)
    
 
    p.close()








end = time.time()
print('Averaging data and creating combined arrays took:', end - start, 's')

#Select latitudes over the southern ocean
#gfdl_tcc_lat = gfdl_tcc_lat[gfdl_tcc_lat[:,0]>=-70]
#gfdl_tcc_lat = gfdl_tcc_lat[gfdl_tcc_lat[:,0]<=-50]

plt.figure()
fig, ax1 = plt.subplots()

#ax2 = ax1.twinx()

ax1.plot(lat,tcc, '-r', label='Total Cloud Fraction')
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