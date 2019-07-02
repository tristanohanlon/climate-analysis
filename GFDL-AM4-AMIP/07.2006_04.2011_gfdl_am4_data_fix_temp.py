# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon

This will create a dataset of time averaged global total cloud cover with latitude.
Time period is from 01.1980 - 12.2014

"""
import time
import sys
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import h5py
import math
from scipy import integrate
from sklearn.impute import SimpleImputer
from scipy import interpolate


#---get plev---#

#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/mri_esm2/amip/cf') #Home PC
os.chdir('D:/MSc/Models/Data/CMIP6/gfdl_am4/') #HDD
#os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CMIP6/gfdl_am4/') #Home PC

f = Dataset('cl_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r') # 420 months

lat = f.variables['lat'][:]
plev = f.variables['lev'][:] #in hPa
a = np.array(f.variables['ap'][:]) #in hPa
b = np.array(f.variables['b'][:]) #in hPa

#---get surface pressure---#

f = Dataset('ps_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r') # 420 months
ps1 = np.array(f.variables['ps'][318:372])
ps2 = np.array(f.variables['ps'][373:376])
ps = np.concatenate((ps1,ps2))

ps = np.mean(ps, axis = 0)

ps_so = np.mean(ps, axis = -1)
ps_so = np.transpose(ps_so)

ps_so = np.vstack((lat, ps_so)).T #creates a (180,34) array
ps_so = ps_so[ps_so[:,0]>=-70]
ps_so = ps_so[ps_so[:,0]<=-50]
ps_so = ps_so[:,1:] #Split the combined array into just the tccf data, eliminating the first coloumn of latitude
ps_so = np.mean(ps_so, axis = 0)


ps = np.mean(ps, axis = 0)
ps = np.mean(ps, axis = 0)

p = a + b*ps
p = np.array(p / 100) #hPa

p_so = a + b*ps_so
p_so = np.array(p_so / 100) #hPa

#---get temp---#

#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/mri_esm2/amip/') #Home PC
os.chdir('D:/MSc/Models/Data/CMIP6/gfdl_am4/') #HDD
#os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CMIP6/gfdl_am4/') #Home PC

f = Dataset('ta_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r') # 420 months
T1 = np.array(f.variables['ta'][318:372])
T2 = np.array(f.variables['ta'][373:376])
T = np.concatenate((T1,T2))
T[T > 400] = np.nan
T = np.nanmean(T, axis = 0)

plev_t = np.array(f.variables['plev'][:]) / 100



T_so = np.nanmean(T, axis = -1)
T_so = np.transpose(T_so)

T_so = np.hstack((np.vstack(lat), T_so)) #creates a (180,34) array
T_so = T_so[T_so[:,0]>=-70]
T_so = T_so[T_so[:,0]<=-50]
T_so = T_so[:,1:] #Split the combined array into just the tccf data, eliminating the first coloumn of latitude
T_so = np.nanmean(T_so, axis = 0)
#T_so = np.vstack((plev_t, T_so)).T

T = np.nanmean(T, axis = -1)
T = np.nanmean(T, axis = -1)


temp_plev = np.vstack((T,plev_t)).T


os.chdir('c:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets') #Home PC
b = h5py.File('2007_2008_gfdl_am4.h5', 'r')

lw_frac = b['lw_frac'][:]
lw_frac_so = b['lw_frac_so'][:]
b.close()

alt = lw_frac[:,0]
h = interpolate.interp1d(p, alt, fill_value="extrapolate", kind = 'cubic')
alt_so = h(p_so)



f = interpolate.interp1d(temp_plev[:,1], temp_plev[:,0], kind = 'cubic')
lw_frac = np.vstack((f(p),lw_frac)).T

g = interpolate.interp1d(plev_t, T_so, kind = 'cubic')
lw_frac_so = np.vstack((g(p_so),lw_frac_so)).T

fig, ax1 = plt.subplots()
ax1.plot(lw_frac_so[:,0],lw_frac_so[:,1], '-r', label='SO')
ax1.plot(lw_frac[:,0],lw_frac[:,1], '-b', label='G')


#----------------------------#




os.chdir('c:/Users/tristan/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets') #Home PC
with h5py.File('07.2006_04.2011_gfdl_am4.h5', 'w') as p:
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_temp', data=alt_t)  
    p.create_dataset('lat', data=lat)  
    p.create_dataset('air_density', data=air_density)  
    p.create_dataset('temp', data=temp_g)
    p.create_dataset('pressure', data=pressure)
        
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    p.create_dataset('tclw_frac', data=tclw_frac)
    p.create_dataset('tciw_frac', data=tciw_frac)
    
    p.create_dataset('cf', data=cf_g)
    p.create_dataset('lw', data=lw_g)
    p.create_dataset('iw', data=iw_g)
    p.create_dataset('lw_frac', data=lw_frac_g)
    p.create_dataset('iw_frac', data=iw_frac_g)
    
    p.create_dataset('cf_so', data=cf_so)
    p.create_dataset('lw_so', data=lw_so)
    p.create_dataset('iw_so', data=iw_so)
    p.create_dataset('lw_frac_so', data=lw_frac_so)
    p.create_dataset('iw_frac_so', data=iw_frac_so)  
    
    p.create_dataset('temp_alt_lat', data=temp_alt_lat)  
    p.create_dataset('cf_alt_lat', data=cf_alt_lat)  
    p.create_dataset('lw_alt_lat', data=lw_alt_lat)
    p.create_dataset('iw_alt_lat', data=iw_alt_lat)
    
    p.create_dataset('cf_t', data=cf_t)
    p.create_dataset('cf_t_so', data=cf_t_so)

    p.create_dataset('lw_t', data=lw_t)
    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('lw_frac_t', data=lw_frac_t)
    p.create_dataset('lw_frac_t_so', data=lw_frac_t_so)

    p.create_dataset('iw_t', data=iw_t)
    p.create_dataset('iw_t_so', data=iw_t_so)
    p.create_dataset('iw_frac_t', data=iw_frac_t)
    p.create_dataset('iw_frac_t_so', data=iw_frac_t_so)
    
 
    p.close()



end = time.time()
print('Averaging data and creating combined arrays took:', end - start, 's')


"""
fig, ax1 = plt.subplots()
ax1.plot(lat,tcc, '-r', label='Total Cloud Fraction')
"""