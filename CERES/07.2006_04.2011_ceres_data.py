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
os.chdir('D:/MSc/Models/Data/CERES') 

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

######
#Test fraction

lw_frac = np.transpose(lwp / (lwp + iwp) * tcc[:,1])
iw_frac = np.transpose(iwp / (lwp + iwp) * tcc[:,1])

plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(lat,tcc[:,1], '-r', label='Cloud Fraction')
ax1.plot(lat,lw_frac, '-b', label='Liquid Fraction')
ax1.plot(lat,iw_frac, '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CERES Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.show()


"""
###############---create combined data---###############
tclw = np.vstack((lat, lwp)).T
tciw = np.vstack((lat, iwp)).T
tclw_frac = np.vstack((lat, np.transpose(lw_frac))).T
tciw_frac = np.vstack((lat, np.transpose(iw_frac))).T



#os.chdir('E:/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets') #Home PC
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets') #laptop
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets') #laptop

with h5py.File('07.2006_04.2011_CERES.h5', 'w') as p:
    
    p.create_dataset('lat', data=lat)  
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    p.create_dataset('tclw_frac', data=tclw_frac)
    p.create_dataset('tciw_frac', data=tciw_frac)

    p.close()


"""