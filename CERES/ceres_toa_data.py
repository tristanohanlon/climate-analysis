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


###############---get Clear Sky TOA SW Flux data---###############
clr_toa_sw_reg=np.zeros((180, 360)) # create a blank array to add cloud amount data
counter=0

# The directory where your HDF files are stored
#os.chdir('E:/University/University/MSc/Models/Data/CERES') 
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CERES') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the total cloud fraction data
    clr_toa_sw_reg = clr_toa_sw_reg+f.select('clr_toa_sw_reg').get()[:]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
lat = np.array(f.select('latitude').get()[:])
lon = np.array(f.select('longitude').get()[:])

clr_toa_sw_reg = np.array(clr_toa_sw_reg/counter)
clr_toa_sw_reg[clr_toa_sw_reg > 1400] = None   


###############---get Clear Sky TOA LW Flux data---###############
clr_toa_lw_reg=np.zeros((180, 360)) # create a blank array to add cloud amount data
counter=0

# The directory where your HDF files are stored
#os.chdir('E:/University/University/MSc/Models/Data/CERES') 
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CERES') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the total cloud fraction data
    clr_toa_lw_reg = clr_toa_lw_reg+f.select('clr_toa_lw_reg').get()[:]
    
    counter+=1

clr_toa_lw_reg = np.array(clr_toa_lw_reg/counter)
clr_toa_lw_reg[clr_toa_lw_reg > 500] = None   

###############---get Clear Sky TOA net Flux data---###############
clr_toa_net_reg=np.zeros((180, 360)) # create a blank array to add cloud amount data
counter=0

# The directory where your HDF files are stored
#os.chdir('E:/University/University/MSc/Models/Data/CERES') 
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CERES') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the total cloud fraction data
    clr_toa_net_reg = clr_toa_net_reg+f.select('clr_toa_net_reg').get()[:]
    
    counter+=1

clr_toa_net_reg = np.array(clr_toa_net_reg/counter)
clr_toa_net_reg[clr_toa_net_reg > 400] = None   


###############---get all Sky TOA SW Flux data---###############
all_toa_sw_reg=np.zeros((180, 360)) # create a blank array to add cloud amount data
counter=0

# The directory where your HDF files are stored
#os.chdir('E:/University/University/MSc/Models/Data/CERES') 
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CERES') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the total cloud fraction data
    all_toa_sw_reg = all_toa_sw_reg+f.select('all_toa_sw_reg').get()[:]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
lat = np.array(f.select('latitude').get()[:])
lon = np.array(f.select('longitude').get()[:])

all_toa_sw_reg = np.array(all_toa_sw_reg/counter)
all_toa_sw_reg[all_toa_sw_reg > 1400] = None   


###############---get all Sky TOA LW Flux data---###############
all_toa_lw_reg=np.zeros((180, 360)) # create a blank array to add cloud amount data
counter=0

# The directory where your HDF files are stored
#os.chdir('E:/University/University/MSc/Models/Data/CERES') 
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CERES') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the total cloud fraction data
    all_toa_lw_reg = all_toa_lw_reg+f.select('all_toa_lw_reg').get()[:]
    
    counter+=1

all_toa_lw_reg = np.array(all_toa_lw_reg/counter)
all_toa_lw_reg[all_toa_lw_reg > 500] = None   

###############---get all Sky TOA net Flux data---###############
all_toa_net_reg=np.zeros((180, 360)) # create a blank array to add cloud amount data
counter=0

# The directory where your HDF files are stored
#os.chdir('E:/University/University/MSc/Models/Data/CERES') 
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CERES') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the total cloud fraction data
    all_toa_net_reg = all_toa_net_reg+f.select('all_toa_net_reg').get()[:]
    
    counter+=1

all_toa_net_reg = np.array(all_toa_net_reg/counter)
all_toa_net_reg[all_toa_net_reg > 400] = None   


###############---get TOA SW insolation (outgoing) data---###############
toa_sw_insol_reg=np.zeros((180, 360)) # create a blank array to add cloud amount data
counter=0

# The directory where your HDF files are stored
#os.chdir('E:/University/University/MSc/Models/Data/CERES') 
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CERES') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the total cloud fraction data
    toa_sw_insol_reg = toa_sw_insol_reg+f.select('all_toa_sw_reg').get()[:]
    
    counter+=1

# Join the two lists as if they were two columns side by side, into a list of two elements each
toa_sw_insol_reg = np.array(toa_sw_insol_reg/counter)
toa_sw_insol_reg[toa_sw_insol_reg > 1400] = None   
"""
plt.figure()
plt.contourf(lon,lat,toa_sw_insol_reg)
plt.show()
"""