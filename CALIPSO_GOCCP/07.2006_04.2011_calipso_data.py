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




###########################---get latitude, longitude and cloud fraction data---###########################

#os.chdir('E:/University/University/MSc/Models/Data/CALIPSO-GOCCP') #Home PC
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CALIPSO-GOCCP') #Home PC
#os.chdir('D:/MSc/Models/Data/CALIPSO-GOCCP') #ext HDD

f = Dataset('Map_OPAQ330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
lat = np.array(f.variables['latitude'][:])
lon = np.array(f.variables['longitude'][:])



#opaque (thick) clouds
cc1 = np.array(f.variables['cltcalipso_opaque'][78:132]) # 7.2006 to 12.2020 - 54 months
cc2 = np.array(f.variables['cltcalipso_opaque'][133:136]) # 2.2011 to 4.2011 - 3 months
cc = np.concatenate((cc1,cc2))

cc[cc < 0] = None #set fill values to nan

cc = np.nanmean(cc, axis = 0) #average over time
tcc = np.nanmean(cc, axis = -1) #average over longitude (lat)



#thin clouds
tf1 = np.array(f.variables['cltcalipso_thin'][78:132]) # 7.2006 to 12.2020 - 54 months
tf2 = np.array(f.variables['cltcalipso_thin'][133:136]) # 2.2011 to 4.2011 - 3 months
tf = np.concatenate((tf1,tf2))

tf[tf < 0] = None #set fill values to nan

tf = np.nanmean(tf, axis = 0) #average over time



#total clouds
cc = cc + tf #(lat, lon)
tcc = np.nanmean(cc, axis = -1) #average over longitude (lat)



###########################---get lat - phase fractions---###########################

f = Dataset('3D_CloudFraction_Phase330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')



#ice water fraction
ice1 = np.array(f.variables['clcalipso_RPIC'][78:132]) # 2.2011 to 4.2011 - 3 months
ice2 = np.array(f.variables['clcalipso_RPIC'][133:136]) # 2.2011 to 4.2011 - 3 months
ice = np.concatenate((ice1,ice2))

ice[ice < 0] = None #set fill values to nan
ice = np.nanmean(ice, axis = 0) #average over time

ice_frac = np.nanmean(ice, axis = 0) #average over alt
ice_frac = np.nanmean(ice_frac, axis = -1) #average over lon

liq_frac = (1 - ice_frac) * tcc
ice_frac = ice_frac  * tcc


"""
plt.figure()
plt.contourf(lat, alt,ice_frac_alt_lat)
plt.show()


plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(lat,tcc[:,1], '-r', label='Cloud Fraction')
ax1.plot(lat,liq_frac, '-b', label='Liquid Fraction')
ax1.plot(lat,ice_frac, '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CALIPSO-GOCCP Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.show()
"""
###########################---get alt - cloud fraction---###########################


f = Dataset('3D_CloudFraction330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
alt = np.array(f.variables['alt_mid'][:]) # km

cff1 = np.array(f.variables['clcalipso'][78:132]) # 7.2006 to 12.2020 - 54 months
cff2 = np.array(f.variables['clcalipso'][133:136]) # 2.2011 to 4.2011 - 3 months
cff = np.concatenate((cff1,cff2))

cff[cff < 0] = None #set fill values to nan

cff = np.nanmean(cff, axis = 0) #average over time
cff = np.nanmean(cff, axis = -1) #average over longitude
cf = np.nanmean(cff, axis = -1) #average over latitude


###########################---get alt - phase fractions---###########################

#ice water fraction with alt and lat

ice_frac_alt_lat = np.nanmean(ice, axis = -1) #average over lon
liq_frac_alt_lat = (1 - ice_frac_alt_lat)

ice_frac_alt_lat = ice_frac_alt_lat *cff
liq_frac_alt_lat = liq_frac_alt_lat *cff

#ice water fraction
iw_frac = np.nanmean(ice, axis = -1) #average over lon
iw_frac = np.nanmean(iw_frac, axis = -1) #average over lat

lw_frac = (1 - iw_frac) * cf
iw_frac = iw_frac * cf

"""

plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(cf, alt, '-r', label='Cloud Fraction')
ax1.plot(lw_frac, alt, '-b', label='Liquid Fraction')
ax1.plot(iw_frac, alt, '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CALIPSO-GOCCP Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.show()
"""
###############---create southern ocean data---###############

cf_so = cff
cf_so = np.transpose(cf_so)

cf_so = np.hstack((np.vstack(lat), cf_so)) #creates a (180,34) array
cf_so = cf_so[cf_so[:,0]>=-70]
cf_so = cf_so[cf_so[:,0]<=-50]
cf_so = cf_so[:,1:] #Split the combined array into just the tccf data, eliminating the first coloumn of latitude
cf_so = np.mean(cf_so, axis = 0)
cf_so = np.vstack((alt, cf_so)).T

lw_so = liq_frac_alt_lat
lw_so = np.transpose(lw_so)

lw_so = np.hstack((np.vstack(lat), lw_so)) #creates a (180,34) array
lw_so = lw_so[lw_so[:,0]>=-70]
lw_so = lw_so[lw_so[:,0]<=-50]
lw_so = lw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
lw_so = np.mean(lw_so, axis = 0)
lw_so = np.vstack((alt, lw_so)).T


iw_so = ice_frac_alt_lat
iw_so = np.transpose(iw_so)

iw_so = np.hstack((np.vstack(lat), iw_so)) #creates a (180,34) array
iw_so = iw_so[iw_so[:,0]>=-70]
iw_so = iw_so[iw_so[:,0]<=-50]
iw_so = iw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
iw_so = np.mean(iw_so, axis = 0)
iw_so = np.vstack((alt, iw_so)).T


###########################---get alt - temp - phase fractions---###########################
#os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CALIPSO-GOCCP') #Home PC
os.chdir('D:/MSc/Models/Data/CALIPSO-GOCCP') #ext HDD

f = Dataset('3D_CloudFraction_Temp330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
alt_t = np.array(f.variables['temp_mid'][:]) + 273 # km (38)

cf_t1 = np.array(f.variables['cltemp'][78:132]) # 7.2006 to 12.2020 - 54 months
cf_t2 = np.array(f.variables['cltemp'][133:136]) # 2.2011 to 4.2011 - 3 months
cf_t = np.concatenate((cf_t1,cf_t2))

cf_t[cf_t < 0] = None #set fill values to nan

cf_t = np.nanmean(cf_t, axis = 0) #average over time
cf_t_lat =  np.nanmean(cf_t, axis = -1) #average over longitude

cf_t = np.nanmean(cf_t, axis = -1) #average over longitude
cf_t = np.nanmean(cf_t, axis = -1) #average over latitude
cf_t = np.vstack((alt_t, cf_t)).T


iw_t1 = np.array(f.variables['cltemp_phase'][78:132]) # 7.2006 to 12.2020 - 54 months
iw_t2 = np.array(f.variables['cltemp_phase'][133:136]) # 2.2011 to 4.2011 - 3 months
iw_t = np.concatenate((iw_t1,iw_t2))

iw_t[iw_t < 0] = None #set fill values to nan

iw_t = np.nanmean(iw_t, axis = 0) #average over time
iw_t_lat =  np.nanmean(iw_t, axis = -1) #average over longitude

lw_t_lat = (1 - iw_t_lat) * cf_t_lat 
iw_t_lat =  iw_t_lat * cf_t_lat 

iw_t = np.nanmean(iw_t_lat, axis = -1) #average over latitude
iw_t = np.vstack((alt_t, iw_t)).T

lw_t = np.nanmean(lw_t_lat, axis = -1) #average over latitude
lw_t = np.vstack((alt_t, lw_t)).T

cf_t_so = cf_t_lat
cf_t_so = np.transpose(cf_t_so)

cf_t_so = np.hstack((np.vstack(lat), cf_t_so)) #creates a (180,34) array
cf_t_so = cf_t_so[cf_t_so[:,0]>=-70]
cf_t_so = cf_t_so[cf_t_so[:,0]<=-50]
cf_t_so = cf_t_so[:,1:] #Split the combined array into just the tccf data, eliminating the first coloumn of latitude
cf_t_so = np.mean(cf_t_so, axis = 0)
cf_t_so = np.vstack((alt_t, cf_t_so)).T

lw_t_so = lw_t_lat
lw_t_so = np.transpose(lw_t_so)

lw_t_so = np.hstack((np.vstack(lat), lw_t_so)) #creates a (180,34) array
lw_t_so = lw_t_so[lw_t_so[:,0]>=-70]
lw_t_so = lw_t_so[lw_t_so[:,0]<=-50]
lw_t_so = lw_t_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
lw_t_so = np.mean(lw_t_so, axis = 0)
lw_t_so = np.vstack((alt_t, lw_t_so)).T


iw_t_so = iw_t_lat
iw_t_so = np.transpose(iw_t_so)

iw_t_so = np.hstack((np.vstack(lat), iw_t_so)) #creates a (180,34) array
iw_t_so = iw_t_so[iw_t_so[:,0]>=-70]
iw_t_so = iw_t_so[iw_t_so[:,0]<=-50]
iw_t_so = iw_t_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
iw_t_so = np.mean(iw_t_so, axis = 0)
iw_t_so = np.vstack((alt_t, iw_t_so)).T


###############---create combined data---###############
tcc = np.vstack((lat, tcc)).T
tclw = np.vstack((lat, liq_frac)).T
tciw = np.vstack((lat, ice_frac)).T

cf = np.vstack((alt, cf)).T
lw_frac = np.vstack((alt, lw_frac)).T
iw_frac = np.vstack((alt, iw_frac)).T
"""
#os.chdir('E:/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets') #Home PC
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets') #laptop
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets') #laptop

with h5py.File('07.2006_04.2011_CALIPSO.h5', 'w') as p:
    
    p.create_dataset('lat', data=lat)  
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_t', data=alt_t)
    
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw_frac', data=tclw)
    p.create_dataset('tciw_frac', data=tciw)
 
    p.create_dataset('cf_t_lat', data=cf_t_lat)
    p.create_dataset('lw_t_lat', data=lw_t_lat)
    p.create_dataset('iw_t_lat', data=iw_t_lat)
    
    p.create_dataset('cf', data=cf)
    p.create_dataset('lw_frac', data=lw_frac)
    p.create_dataset('iw_frac', data=iw_frac)
    
    p.create_dataset('cf_so', data=cf_so)
    p.create_dataset('lw_frac_so', data=lw_so)
    p.create_dataset('iw_frac_so', data=iw_so)   

    p.create_dataset('cf_t', data=cf_t)
    p.create_dataset('lw_t_frac', data=lw_t)
    p.create_dataset('iw_t_frac', data=iw_t)   

    p.create_dataset('cf_t_so', data=cf_t_so)
    p.create_dataset('lw_t_frac_so', data=lw_t_so)
    p.create_dataset('iw_t_frac_so', data=iw_t_so)   
    
    p.create_dataset('cf_alt_lat', data=cff)
    p.create_dataset('liq_frac_alt_lat', data=liq_frac_alt_lat)
    p.create_dataset('ice_frac_alt_lat', data=ice_frac_alt_lat)
    
    p.close()

"""
"""
 



plt.figure()
fig, ax1 = plt.subplots()

#ax1.plot(cf_t, alt_t, '-r', label='Cloud Fraction')
ax1.contourf(lat, alt_t[8:], cf_t_lat[8:])

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CALIPSO-GOCCP Global Cloud Fraction and Phase Fraction vs Latitude')
plt.gca().invert_yaxis()
plt.grid(True)
plt.show()

"""





