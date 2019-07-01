# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will extract the cloud fraction profile against altitude. 
The global cloud fraction is already averaged over all longitude and latitude at each altitude layer.
[:,0] = altitude
[:,1] = cloud fraction
"""
import time
import numpy as np
import os
from pyhdf import SD
import matplotlib.pyplot as plt
import h5py
import math
from scipy import integrate


###############---get tcc data---###############
tcc=np.zeros((1, 180)) # create a blank array to add cloud amount data
counter=0

# The directory where your HDF files are stored
#os.chdir('E:/University/University/MSc/Models/Data/CERES') 
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CERES') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the total cloud fraction data
    tcc = tcc+f.select('cld_amount_zon').get()[-1:]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
lat = np.array(f.select('latitude').get()[:])
tcc = np.array((tcc/100)/counter)
tcc = np.vstack((lat, tcc)).T

###############---get cf data---###############
cf=np.zeros(4) # create a blank array to add cloud amount data
alt=np.zeros(4) # create a blank array to add altitude data
counter=0

# The directory where your HDF files are stored

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the altitude data as a list
    alt = alt+f.select('cld_eff_hgt_glob').get()[:-1]
    
    # Get the cloud fraction data as a list, excluding the last value which is the total for all altitudes
    cf = cf+f.select('cld_amount_glob').get()[:-1]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
alt = np.array(alt/counter)
cf = np.array(cf/counter)
cf = np.vstack((alt, (cf/100))).T


###############---get lwp---###############
lwp = np.zeros((1, 180)) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the liquid water path (gm^-2), excluding the last value which is the total for all altitudes
    lwp=lwp+f.select('cld_lwp_zon').get()[-1:]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
lwp = np.array((lwp/1000)/counter) #kgm^-2

lwp[lwp > 1] = None   


###############---get iwp---###############
iwp = np.zeros((1, 180)) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the liquid water path (gm^-2), excluding the last value which is the total for all altitudes
    iwp=iwp+f.select('cld_iwp_zon').get()[-1:]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
iwp = np.array((iwp/1000)/counter) #kgm^-2
iwp[iwp > 1] = None   


###############---get tclw_frac---###############
tclw_frac = np.zeros((1, 180)) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the liquid water path (gm^-2), excluding the last value which is the total for all altitudes
    tclw_frac=tclw_frac+f.select('cld_amount_liq_zon').get()[-1:]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
tclw_frac = np.mean(tclw_frac, axis = 0)    
tclw_frac = np.array((tclw_frac/100)/counter) 

#need to correct fill values

###############---get tciw_frac---###############
tciw_frac = np.zeros((1, 180)) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the liquid water path (gm^-2), excluding the last value which is the total for all altitudes
    tciw_frac=tciw_frac+f.select('cld_amount_ice_zon').get()[-1:]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
tciw_frac = np.mean(tciw_frac, axis = 0)    
tciw_frac = np.array((tciw_frac/100)/counter) 
#need to correct fill values

"""
###############---determine air path---###############
p = f.select('cld_eff_press_glob').get()[:-1]
temp = f.select('cld_eff_temp_glob').get()[:-1]

alt_m = np.hstack(alt*1000)
p = np.vstack(p)
pressure = np.hstack((np.vstack(alt), p))
temp_g = np.vstack((alt, temp)).T

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = (pressure[:,1] * 100) / (286.9 * temp_g[:,1])

ap = [18800.0,6686.0,4312.0,1925.0]
ap = np.array(ap)
#ap = - integrate.trapz(air_density, alt_m)

lw = np.divide((lw_frac, ap))

x = np.vstack(np.transpose(lwp / (iwp + lwp)))
y = np.vstack(tcc[:,1])
z = np.multiply(x,y)


###############---get iw_frac data---###############
iw_frac=np.zeros(4) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the ice water cloud fraction data as a list, excluding the last value which is the total for all altitudes
    iw_frac = iw_frac+f.select('cld_amount_ice_glob').get()[:-1]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
iw_frac = np.array(iw_frac/counter)
iw_frac = np.vstack((alt, (iw_frac/100))).T


###############---get lw_frac data---###############
lw_frac=np.zeros(4) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the liquid water cloud fraction data as a list, excluding the last value which is the total for all altitudes
    lw_frac = lw_frac+f.select('cld_amount_liq_glob').get()[:-1]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
lw_frac = np.array(lw_frac/counter)
lw_frac = np.vstack((alt, (lw_frac/100))).T

"""
###############---get cf_alt_lat data---###############
cf_alt_lat = np.zeros((4, 180)) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the cloud fraction data as a list, excluding the last value which is the total for all altitudes
    cf_alt_lat=cf_alt_lat+f.select('cld_amount_zon').get()[:-1]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
    
cf_alt_lat = np.array((cf_alt_lat/100)/counter)

###############---get lw_alt_lat data---###############
lw_alt_lat = np.zeros((4, 180)) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the liquid water cloud fraction data as a list, excluding the last value which is the total for all altitudes
    lw_alt_lat=lw_alt_lat+f.select('cld_amount_liq_zon').get()[:-1]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
    
lw_alt_lat = np.array((lw_alt_lat/100)/counter)
######################
lw = np.zeros((4, 180)) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the liquid water cloud fraction data as a list, excluding the last value which is the total for all altitudes
    lw_alt_lat=lw_alt_lat+f.select('cld_amount_liq_zon').get()[:-1]
    
    counter+=1
###############---get iw_alt_lat data---###############
iw_alt_lat = np.zeros((4, 180)) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the ice water cloud fraction data as a list, excluding the last value which is the total for all altitudes
    iw_alt_lat=iw_alt_lat+f.select('cld_amount_ice_zon').get()[:-1]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
    
iw_alt_lat = np.array((iw_alt_lat/100)/counter)

###############---get temp_alt_lat data---###############
temp_alt_lat = np.zeros((4, 180)) # create a blank array to add cloud amount data
counter=0

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the ice water cloud fraction data as a list, excluding the last value which is the total for all altitudes
    temp_alt_lat=temp_alt_lat+f.select('cld_eff_temp_zon').get()[:-1]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
    
temp_alt_lat = np.array(temp_alt_lat/counter)

temp_alt_lat[temp_alt_lat > 290] = None   

###############---create combined data---###############
temp = np.vstack((alt, temp)).T

tclw = np.vstack((lat, lwp)).T
tclw_frac = np.vstack((lat, tclw_frac)).T

tciw = np.vstack((lat, iwp)).T
tciw_frac = np.vstack((lat, tciw_frac)).T

cf_t = np.vstack((temp[:,1], cf[:,1])).T
lw_t = np.vstack((temp[:,1], lw_frac[:,1])).T
iw_t = np.vstack((temp[:,1], iw_frac[:,1])).T

#os.chdir('E:/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets') #Home PC
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets') #laptop

with h5py.File('07.2006_04.2011_CERES.h5', 'w') as p:
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('lat', data=lat)  
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    p.create_dataset('tclw_frac', data=tclw_frac)
    p.create_dataset('tciw_frac', data=tciw_frac)
    p.create_dataset('cf', data=cf)
    p.create_dataset('lw_frac', data=lw_frac)
    p.create_dataset('iw_frac', data=lw_frac)
    p.create_dataset('temp', data=temp)
#    p.create_dataset('cf_so', data=cf_so)
#    p.create_dataset('lw_so', data=lw_so)
#    p.create_dataset('iw_so', data=iw_so)
#    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure)
    p.create_dataset('temp_alt_lat', data=temp_alt_lat)  
    p.create_dataset('cf_alt_lat', data=cf_alt_lat)  
    p.create_dataset('lw_alt_lat', data=lw_alt_lat)
    p.create_dataset('iw_alt_lat', data=iw_alt_lat)
    p.create_dataset('cf_t', data=cf_t)
#    p.create_dataset('cf_t_so', data=cf_t_so)
    p.create_dataset('lw_t', data=lw_t)
#    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('iw_t', data=iw_t)
#    p.create_dataset('iw_t_so', data=iw_t_so)
    p.close()


"""
plt.figure()
plt.contourf(lat,alt,temp_alt_lat)
plt.show()
"""