# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This import the reduced datasets from 2001 to 2005 CCCM, ECMWF and gfdl4. 
The code can select both global or southern ocean data.

"""
import time
import matplotlib.pyplot as plt
import h5py
import os
import numpy as np
import scipy.interpolate as interp

start = time.time()


#---Importing Data from Reduced Datasets---#
"""
# Uni Laptop
#ECMWF-ERA5 Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
a = h5py.File('2001_2005_ecmwf_era5.h5', 'r')

#GFDL-AM4-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('2001_2005_gfdl_am4.h5', 'r')

#GFDL-HIRAM-C360 Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL-HIRAM-C360/reduced_datasets')
h = h5py.File('2001_2005_gfdl_hiram.h5', 'r')

#MRI-ESM2-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
d = h5py.File('2001_2005_mri_esm2.h5', 'r')

#MRI-CGCM3-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/MRI-CGCM3-AMIP/reduced_datasets')
i = h5py.File('2001_2005_mri_cgcm3.h5', 'r')

#CESM1-CAM5-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CESM1-CAM5-AMIP/reduced_datasets')
c = h5py.File('2001_2005_cesm1_cam5.h5', 'r')

#CESM2-CAM6-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets')
e = h5py.File('2001_2005_cesm2_cam6.h5', 'r')

#GISS-E2R-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GISS-E2R-AMIP/reduced_datasets')
j = h5py.File('2001_2005_giss_e2r.h5', 'r')

#GISS-E21G-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GISS-E21G-AMIP/reduced_datasets')
k = h5py.File('2001_2005_giss_e21g.h5', 'r')

#IPSL-CM5A-LR-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/IPSL-CM5A-LR-AMIP/reduced_datasets')
l = h5py.File('2001_2005_ipsl_cm5a_lr.h5', 'r')

#IPSL-CM6A-LR-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/IPSL-CM6A-LR-AMIP/reduced_datasets')
m = h5py.File('2001_2005_ipsl_cm6a_lr.h5', 'r')

#MIROC5-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/MIROC5-AMIP/reduced_datasets')
n = h5py.File('2001_2005_miroc5.h5', 'r')

#MIROC6-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/MIROC6-AMIP/reduced_datasets')
o = h5py.File('2001_2005_miroc6.h5', 'r')
"""


# Home PC
#ECMWF Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
a = h5py.File('2001_2005_ecmwf_era5.h5', 'r')

#GFDL-AM4-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('2001_2005_gfdl_am4.h5', 'r')

#GFDL-HIRAM-C360 Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL-HIRAM-C360/reduced_datasets')
h = h5py.File('2001_2005_gfdl_hiram.h5', 'r')

#MRI-ESM2-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
d = h5py.File('2001_2005_mri_esm2.h5', 'r')

#MRI-CGCM3-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/MRI-CGCM3-AMIP/reduced_datasets')
i = h5py.File('2001_2005_mri_cgcm3.h5', 'r')

#CESM1-CAM5-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CESM1-CAM5-AMIP/reduced_datasets')
c = h5py.File('2001_2005_cesm1_cam5.h5', 'r')

#CESM2-CAM6-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets')
e = h5py.File('2001_2005_cesm2_cam6.h5', 'r')

#GISS-E2R-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GISS-E2R-AMIP/reduced_datasets')
j = h5py.File('2001_2005_giss_e2r.h5', 'r')

#GISS-E21G-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GISS-E21G-AMIP/reduced_datasets')
k = h5py.File('2001_2005_giss_e21g.h5', 'r')

#IPSL-CM5A-LR-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/IPSL-CM5A-LR-AMIP/reduced_datasets')
l = h5py.File('2001_2005_ipsl_cm5a_lr.h5', 'r')

#IPSL-CM6A-LR-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/IPSL-CM6A-LR-AMIP/reduced_datasets')
m = h5py.File('2001_2005_ipsl_cm6a_lr.h5', 'r')

#MIROC5-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/MIROC5-AMIP/reduced_datasets')
n = h5py.File('2001_2005_miroc5.h5', 'r')

#MIROC6-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/MIROC6-AMIP/reduced_datasets')
o = h5py.File('2001_2005_miroc6.h5', 'r')


"""
# Laptop
#ECMWF Data
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
#a = h5py.File('2001_2005_ecmwf.h5', 'r')

#GFDL-AM4-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('2001_2005_gfdl_am4.h5', 'r')

#GFDL-HIRAM-C360 Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/GFDL-HIRAM-C360/reduced_datasets')
h = h5py.File('2001_2005_gfdl_hiram.h5', 'r')

#MRI-ESM2-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
d = h5py.File('2001_2005_mri_esm2.h5', 'r')

#MRI-CGCM3-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/MRI-CGCM3-AMIP/reduced_datasets')
i = h5py.File('2001_2005_mri_cgcm3.h5', 'r')

#CESM1-CAM5-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CESM1-CAM5-AMIP/reduced_datasets')
c = h5py.File('2001_2005_cesm1_cam5.h5', 'r')

#CESM2-CAM6-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets')
e = h5py.File('2001_2005_cesm2_cam6.h5', 'r')

#GISS-E2R-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/GISS-E2R-AMIP/reduced_datasets')
j = h5py.File('2001_2005_giss_e2r.h5', 'r')

#GISS-E21G-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/GISS-E21G-AMIP/reduced_datasets')
k = h5py.File('2001_2005_giss_e21g.h5', 'r')

#IPSL-CM5A-LR-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/IPSL-CM5A-LR-AMIP/reduced_datasets')
l = h5py.File('2001_2005_ipsl_cm5a_lr.h5', 'r')

#IPSL-CM6A-LR-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/IPSL-CM6A-LR-AMIP/reduced_datasets')
m = h5py.File('2001_2005_ipsl_cm6a_lr.h5', 'r')

#MIROC5-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/MIROC5-AMIP/reduced_datasets')
n = h5py.File('2001_2005_miroc5.h5', 'r')

#MIROC6-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/MIROC6-AMIP/reduced_datasets')
o = h5py.File('2001_2005_miroc6.h5', 'r')

"""


############################################################################### ECMWF Data

#---ECMWF Global Latitude Data---#

ecmwf_tcc_lat_g = a['tcc'][:] # 0-1
ecmwf_tclw_lat_g = a['tclw'][:] #kgm^-2
ecmwf_tciw_lat_g = a['tciw'][:] #kgm^-2

#---ECMWF Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
ecmwf_tclw_frac_lat_g = a['tclw_frac'][:] 
ecmwf_tciw_frac_lat_g = a['tciw_frac'][:]


#---ECMWF Southern Ocean Latitude Data---#

ecmwf_tcc_lat_so = ecmwf_tcc_lat_g[ecmwf_tcc_lat_g[:,0]>=-70]
ecmwf_tcc_lat_so = ecmwf_tcc_lat_so[ecmwf_tcc_lat_so[:,0]<=-50] # 0-1

ecmwf_tclw_lat_so = ecmwf_tclw_lat_g[ecmwf_tclw_lat_g[:,0]>=-70]
ecmwf_tclw_lat_so = ecmwf_tclw_lat_so[ecmwf_tclw_lat_so[:,0]<=-50] #kgm^-2

ecmwf_tciw_lat_so = ecmwf_tciw_lat_g[ecmwf_tciw_lat_g[:,0]>=-70]
ecmwf_tciw_lat_so = ecmwf_tciw_lat_so[ecmwf_tciw_lat_so[:,0]<=-50] #kgm^-2


#---ECMWF lat-alt contour data---#

ecmwf_tcc_alt_lat = a['cf_alt_lat'][:] #kg/kg
ecmwf_tclw_alt_lat = a['lw_alt_lat'][:] #kg/kg
ecmwf_tciw_alt_lat = a['iw_alt_lat'][:] #kg/kg
ecmwf_temp_alt_lat = a['temp_alt_lat'][:] #kg/kg
ecmwf_lat = a['lat'][:]
ecmwf_alt = a['alt'][:]

ecmwf_tclw_frac_alt_lat = (ecmwf_tclw_alt_lat / (ecmwf_tclw_alt_lat + ecmwf_tciw_alt_lat)) * ecmwf_tcc_alt_lat
ecmwf_tciw_frac_alt_lat = (ecmwf_tciw_alt_lat / (ecmwf_tclw_alt_lat + ecmwf_tciw_alt_lat)) * ecmwf_tcc_alt_lat


#---ECMWF Global Profile---#

ecmwf_tcc_alt_g = a['cf'][:] # 0-1
ecmwf_tclw_alt_g = a['lw'][:] #kg/kg
ecmwf_tciw_alt_g = a['iw'][:] #kg/kg
ecmwf_temp_alt_g = a['temp'][:] #K
ecmwf_plevel_alt_g = a['pressure'][:] #hPa

ecmwf_tcc_temp_g = a['cf_t'][:] # 0-1
ecmwf_tclw_temp_g = a['lw_t'][:] #kg/kg
ecmwf_tciw_temp_g = a['iw_t'][:] #kg/kg

ecmwf_tclw_frac_temp_g = a['lw_frac_t'][:]
ecmwf_tciw_frac_temp_g = a['iw_frac_t'][:]

#---ECMWF Phase Profile Fractions---#

ecmwf_tclw_frac_alt_g = a['lw_frac'][:]
ecmwf_tciw_frac_alt_g = a['iw_frac'][:]
ecmwf_tclw_frac_alt_so = a['lw_frac_so'][:]
ecmwf_tciw_frac_alt_so = a['iw_frac_so'][:]


#---ECMWF Southern Ocean Profile---#

ecmwf_tcc_alt_so = a['cf_so'][:] # 0-1
ecmwf_tclw_alt_so = a['lw_so'][:] #kg/kg
ecmwf_tciw_alt_so = a['iw_so'][:] #kg/kg

ecmwf_tcc_temp_so = a['cf_t_so'][:] # 0-1
ecmwf_tclw_temp_so = a['lw_t_so'][:] #kg/kg
ecmwf_tciw_temp_so = a['iw_t_so'][:] #kg/kg

ecmwf_tclw_frac_temp_so = a['lw_frac_t_so'][:]
ecmwf_tciw_frac_temp_so = a['lw_frac_t_so'][:]



############################################################################### GFDL-AM4-AMIP Data

#---GFDL-AM4-AMIP Global Latitude Data---#

gfdl4_tcc_lat_g = b['tcc'][:] # 0-1
gfdl4_tclw_lat_g = b['tclw'][:] #kgm^-2
gfdl4_tciw_lat_g = b['tciw'][:] #kgm^-2

#---GFDL-AM4-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
gfdl4_tclw_frac_lat_g = b['tclw_frac'][:] 
gfdl4_tciw_frac_lat_g = b['tciw_frac'][:]


#---GFDL-AM4-AMIP Southern Ocean Latitude Data---#

gfdl4_tcc_lat_so = gfdl4_tcc_lat_g[gfdl4_tcc_lat_g[:,0]>=-70]
gfdl4_tcc_lat_so = gfdl4_tcc_lat_so[gfdl4_tcc_lat_so[:,0]<=-50] # 0-1

gfdl4_tclw_lat_so = gfdl4_tclw_lat_g[gfdl4_tclw_lat_g[:,0]>=-70]
gfdl4_tclw_lat_so = gfdl4_tclw_lat_so[gfdl4_tclw_lat_so[:,0]<=-50] #kgm^-2

gfdl4_tciw_lat_so = gfdl4_tciw_lat_g[gfdl4_tciw_lat_g[:,0]>=-70]
gfdl4_tciw_lat_so = gfdl4_tciw_lat_so[gfdl4_tciw_lat_so[:,0]<=-50] #kgm^-2


#---GFDL-AM4-AMIP lat-alt contour data---#

gfdl4_tcc_alt_lat = b['cf_alt_lat'][:] #kg/kg
gfdl4_tclw_alt_lat = b['lw_alt_lat'][:] #kg/kg
gfdl4_tciw_alt_lat = b['iw_alt_lat'][:] #kg/kg
gfdl4_temp_alt_lat = b['temp_alt_lat'][:] #kg/kg
gfdl4_lat = b['lat'][:]
gfdl4_alt = b['alt'][:]
gfdl4_alt_temp = b['alt_temp'][:]

gfdl4_tclw_frac_alt_lat = (gfdl4_tclw_alt_lat / (gfdl4_tclw_alt_lat + gfdl4_tciw_alt_lat)) * gfdl4_tcc_alt_lat
gfdl4_tciw_frac_alt_lat = (gfdl4_tciw_alt_lat / (gfdl4_tclw_alt_lat + gfdl4_tciw_alt_lat)) * gfdl4_tcc_alt_lat


#---GFDL-AM4-AMIP Global Profile---#

gfdl4_tcc_alt_g = b['cf'][:] # 0-1
gfdl4_tclw_alt_g = b['lw'][:] #kg/kg
gfdl4_tciw_alt_g = b['iw'][:] #kg/kg
gfdl4_temp_alt_g = b['temp'][:] #K
gfdl4_plevel_alt_g = b['pressure'][:] #hPa

gfdl4_tcc_temp_g = b['cf_t'][:] # 0-1
gfdl4_tclw_temp_g = b['lw_t'][:] #kg/kg
gfdl4_tciw_temp_g = b['iw_t'][:] #kg/kg

gfdl4_tclw_frac_temp_g = b['lw_frac_t'][:]
gfdl4_tciw_frac_temp_g = b['iw_frac_t'][:]

#---GFDL-AM4-AMIP Phase Profile Fractions---#

gfdl4_tclw_frac_alt_g = b['lw_frac'][:]
gfdl4_tciw_frac_alt_g = b['iw_frac'][:]
gfdl4_tclw_frac_alt_so = b['lw_frac_so'][:]
gfdl4_tciw_frac_alt_so = b['iw_frac_so'][:]


#---GFDL-AM4-AMIP Southern Ocean Profile---#

gfdl4_tcc_alt_so = b['cf_so'][:] # 0-1
gfdl4_tclw_alt_so = b['lw_so'][:] #kg/kg
gfdl4_tciw_alt_so = b['iw_so'][:] #kg/kg

gfdl4_tcc_temp_so = b['cf_t_so'][:] # 0-1
gfdl4_tclw_temp_so = b['lw_t_so'][:] #kg/kg
gfdl4_tciw_temp_so = b['iw_t_so'][:] #kg/kg

gfdl4_tclw_frac_temp_so = b['lw_frac_t_so'][:]
gfdl4_tciw_frac_temp_so = b['lw_frac_t_so'][:]


############################################################################### GFDL-HIRAM-C360 Data

#---GFDL-HIRAM-C360 Global Latitude Data---#

gfdl_hiram_tcc_lat_g = h['tcc'][:] # 0-1
gfdl_hiram_tclw_lat_g = h['tclw'][:] #kgm^-2
gfdl_hiram_tciw_lat_g = h['tciw'][:] #kgm^-2

#---GFDL-HIRAM-C360 Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
gfdl_hiram_tclw_frac_lat_g = h['tclw_frac'][:] 
gfdl_hiram_tciw_frac_lat_g = h['tciw_frac'][:]


#---GFDL-HIRAM-C360 Southern Ocean Latitude Data---#

gfdl_hiram_tcc_lat_so = gfdl_hiram_tcc_lat_g[gfdl_hiram_tcc_lat_g[:,0]>=-70]
gfdl_hiram_tcc_lat_so = gfdl_hiram_tcc_lat_so[gfdl_hiram_tcc_lat_so[:,0]<=-50] # 0-1

gfdl_hiram_tclw_lat_so = gfdl_hiram_tclw_lat_g[gfdl_hiram_tclw_lat_g[:,0]>=-70]
gfdl_hiram_tclw_lat_so = gfdl_hiram_tclw_lat_so[gfdl_hiram_tclw_lat_so[:,0]<=-50] #kgm^-2

gfdl_hiram_tciw_lat_so = gfdl_hiram_tciw_lat_g[gfdl_hiram_tciw_lat_g[:,0]>=-70]
gfdl_hiram_tciw_lat_so = gfdl_hiram_tciw_lat_so[gfdl_hiram_tciw_lat_so[:,0]<=-50] #kgm^-2


#---GFDL-HIRAM-C360 lat-alt contour data---#

gfdl_hiram_tcc_alt_lat = h['cf_alt_lat'][:] #kg/kg
gfdl_hiram_tclw_alt_lat = h['lw_alt_lat'][:] #kg/kg
gfdl_hiram_tciw_alt_lat = h['iw_alt_lat'][:] #kg/kg
gfdl_hiram_temp_alt_lat = h['temp_alt_lat'][:] #kg/kg
gfdl_hiram_lat = h['lat'][:]
gfdl_hiram_alt = h['alt'][:]
gfdl_hiram_alt_temp = h['alt_temp'][:]

gfdl_hiram_tclw_frac_alt_lat = (gfdl_hiram_tclw_alt_lat / (gfdl_hiram_tclw_alt_lat + gfdl_hiram_tciw_alt_lat)) * gfdl_hiram_tcc_alt_lat
gfdl_hiram_tciw_frac_alt_lat = (gfdl_hiram_tciw_alt_lat / (gfdl_hiram_tclw_alt_lat + gfdl_hiram_tciw_alt_lat)) * gfdl_hiram_tcc_alt_lat

#---GFDL-HIRAM-C360 Global Profile---#

gfdl_hiram_tcc_alt_g = h['cf'][:] # 0-1
gfdl_hiram_tclw_alt_g = h['lw'][:] #kg/kg
gfdl_hiram_tciw_alt_g = h['iw'][:] #kg/kg
gfdl_hiram_temp_alt_g = h['temp'][:] #K
gfdl_hiram_plevel_alt_g = h['pressure'][:] #hPa


gfdl_hiram_tcc_temp_g = h['cf_t'][:] # 0-1
gfdl_hiram_tclw_temp_g = h['lw_t'][:] #kg/kg
gfdl_hiram_tciw_temp_g = h['iw_t'][:] #kg/kg

gfdl_hiram_tclw_frac_temp_g = h['lw_frac_t'][:]
gfdl_hiram_tciw_frac_temp_g = h['iw_frac_t'][:]


#---GFDL-HIRAM-C360 Phase Profile Fractions---#

gfdl_hiram_tclw_frac_alt_g = h['lw_frac'][:]
gfdl_hiram_tciw_frac_alt_g = h['iw_frac'][:]
gfdl_hiram_tclw_frac_alt_so = h['lw_frac_so'][:]
gfdl_hiram_tciw_frac_alt_so = h['iw_frac_so'][:]

#---GFDL-HIRAM-C360 Southern Ocean Profile---#

gfdl_hiram_tcc_alt_so = h['cf_so'][:] # 0-1
gfdl_hiram_tclw_alt_so = h['lw_so'][:] #kg/kg
gfdl_hiram_tciw_alt_so = h['iw_so'][:] #kg/kg

gfdl_hiram_tcc_temp_so = h['cf_t_so'][:] # 0-1
gfdl_hiram_tclw_temp_so = h['lw_t_so'][:] #kg/kg
gfdl_hiram_tciw_temp_so = h['iw_t_so'][:] #kg/kg

gfdl_hiram_tclw_frac_temp_so = h['lw_frac_t_so'][:]
gfdl_hiram_tciw_frac_temp_so = h['lw_frac_t_so'][:]


############################################################################### MRI_ESM2-AMIP Data

#---MRI_ESM2-AMIP Global Latitude Data---#

mri_tcc_lat_g = d['tcc'][:] # 0-1
mri_tclw_lat_g = d['tclw'][:] #kgm^-2
mri_tciw_lat_g = d['tciw'][:] #kgm^-2


#---MRI_ESM2-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
mri_tclw_frac_lat_g = d['tclw_frac'][:] 
mri_tciw_frac_lat_g = d['tciw_frac'][:]



#---MRI_ESM2-AMIP Southern Ocean Latitude Data---#

mri_tcc_lat_so = mri_tcc_lat_g[mri_tcc_lat_g[:,0]>=-70]
mri_tcc_lat_so = mri_tcc_lat_so[mri_tcc_lat_so[:,0]<=-50] # 0-1

mri_tclw_lat_so = mri_tclw_lat_g[mri_tclw_lat_g[:,0]>=-70]
mri_tclw_lat_so = mri_tclw_lat_so[mri_tclw_lat_so[:,0]<=-50] #kgm^-2

mri_tciw_lat_so = mri_tciw_lat_g[mri_tciw_lat_g[:,0]>=-70]
mri_tciw_lat_so = mri_tciw_lat_so[mri_tciw_lat_so[:,0]<=-50] #kgm^-2


#---MRI_ESM2-AMIP lat-alt contour data---#

mri_tcc_alt_lat = d['cf_alt_lat'][:] #kg/kg
mri_tclw_alt_lat = d['lw_alt_lat'][:] #kg/kg
mri_tciw_alt_lat = d['iw_alt_lat'][:] #kg/kg
mri_temp_alt_lat = d['temp_alt_lat'][:] #kg/kg
mri_lat = d['lat'][:]
mri_alt = d['alt'][:]
mri_alt_temp = d['alt_temp'][:]


mri_tclw_frac_alt_lat = (mri_tclw_alt_lat / (mri_tclw_alt_lat + mri_tciw_alt_lat)) * mri_tcc_alt_lat
mri_tciw_frac_alt_lat = (mri_tciw_alt_lat / (mri_tclw_alt_lat + mri_tciw_alt_lat)) * mri_tcc_alt_lat


#---MRI_ESM2-AMIP Global Profile---#

mri_tcc_alt_g = d['cf'][:] # 0-1
mri_tclw_alt_g = d['lw'][:] #kg/kg
mri_tciw_alt_g = d['iw'][:] #kg/kg
mri_temp_alt_g = d['temp'][:] #K
mri_plevel_alt_g = d['pressure'][:] #hPa

mri_tcc_temp_g = d['cf_t'][:] # 0-1
mri_tclw_temp_g = d['lw_t'][:] #kg/kg
mri_tciw_temp_g = d['iw_t'][:] #kg/kg

mri_tclw_frac_temp_g = d['lw_frac_t'][:]
mri_tciw_frac_temp_g = d['iw_frac_t'][:]


#---MRI_ESM2-AMIP Phase Profile Fractions---#

mri_tclw_frac_alt_g = d['lw_frac'][:]
mri_tciw_frac_alt_g = d['iw_frac'][:]
mri_tclw_frac_alt_so = d['lw_frac_so'][:]
mri_tciw_frac_alt_so = d['iw_frac_so'][:]


#---MRI_ESM2-AMIP Southern Ocean Profile---#

mri_tcc_alt_so = d['cf_so'][:] # 0-1
mri_tclw_alt_so = d['lw_so'][:] #kg/kg
mri_tciw_alt_so = d['iw_so'][:] #kg/kg

mri_tcc_temp_so = d['cf_t_so'][:] # 0-1
mri_tclw_temp_so = d['lw_t_so'][:] #kg/kg
mri_tciw_temp_so = d['iw_t_so'][:] #kg/kg

mri_tclw_frac_temp_so = d['lw_frac_t_so'][:]
mri_tciw_frac_temp_so = d['lw_frac_t_so'][:]


############################################################################### MRI-CGCM3-AMIP Data

#---MRI-CGCM3-AMIP Global Latitude Data---#

mri_cgcm_tcc_lat_g = i['tcc'][:] # 0-1
mri_cgcm_tclw_lat_g = i['tclw'][:] #kgm^-2
mri_cgcm_tciw_lat_g = i['tciw'][:] #kgm^-2


#---MRI-CGCM3-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
mri_cgcm_tclw_frac_lat_g = i['tclw_frac'][:] 
mri_cgcm_tciw_frac_lat_g = i['tciw_frac'][:]



#---MRI-CGCM3-AMIP Southern Ocean Latitude Data---#

mri_cgcm_tcc_lat_so = mri_cgcm_tcc_lat_g[mri_cgcm_tcc_lat_g[:,0]>=-70]
mri_cgcm_tcc_lat_so = mri_cgcm_tcc_lat_so[mri_cgcm_tcc_lat_so[:,0]<=-50] # 0-1

mri_cgcm_tclw_lat_so = mri_cgcm_tclw_lat_g[mri_cgcm_tclw_lat_g[:,0]>=-70]
mri_cgcm_tclw_lat_so = mri_cgcm_tclw_lat_so[mri_cgcm_tclw_lat_so[:,0]<=-50] #kgm^-2

mri_cgcm_tciw_lat_so = mri_cgcm_tciw_lat_g[mri_cgcm_tciw_lat_g[:,0]>=-70]
mri_cgcm_tciw_lat_so = mri_cgcm_tciw_lat_so[mri_cgcm_tciw_lat_so[:,0]<=-50] #kgm^-2


#---MRI-CGCM3-AMIP lat-alt contour data---#

mri_cgcm_tcc_alt_lat = i['cf_alt_lat'][:] #kg/kg
mri_cgcm_tclw_alt_lat = i['lw_alt_lat'][:] #kg/kg
mri_cgcm_tciw_alt_lat = i['iw_alt_lat'][:] #kg/kg
mri_cgcm_temp_alt_lat = i['temp_alt_lat'][:] #kg/kg
mri_cgcm_lat = i['lat'][:]
mri_cgcm_alt = i['alt'][:]
mri_cgcm_alt_temp = i['alt_temp'][:]


mri_cgcm_tclw_frac_alt_lat = (mri_cgcm_tclw_alt_lat / (mri_cgcm_tclw_alt_lat + mri_cgcm_tciw_alt_lat)) * mri_cgcm_tcc_alt_lat
mri_cgcm_tciw_frac_alt_lat = (mri_cgcm_tciw_alt_lat / (mri_cgcm_tclw_alt_lat + mri_cgcm_tciw_alt_lat)) * mri_cgcm_tcc_alt_lat


#---MRI-CGCM3-AMIP Global Profile---#

mri_cgcm_tcc_alt_g = i['cf'][:] # 0-1
mri_cgcm_tclw_alt_g = i['lw'][:] #kg/kg
mri_cgcm_tciw_alt_g = i['iw'][:] #kg/kg
mri_cgcm_temp_alt_g = i['temp'][:] #K
mri_cgcm_plevel_alt_g = i['pressure'][:] #hPa


mri_cgcm_tcc_temp_g = i['cf_t'][:] # 0-1
mri_cgcm_tclw_temp_g = i['lw_t'][:] #kg/kg
mri_cgcm_tciw_temp_g = i['iw_t'][:] #kg/kg

mri_cgcm_tclw_frac_temp_g = i['lw_frac_t'][:]
mri_cgcm_tciw_frac_temp_g = i['iw_frac_t'][:]



#---MRI-CGCM3-AMIP Phase Profile Fractions---#

mri_cgcm_tclw_frac_alt_g = i['lw_frac'][:]
mri_cgcm_tciw_frac_alt_g = i['iw_frac'][:]
mri_cgcm_tclw_frac_alt_so = i['lw_frac_so'][:]
mri_cgcm_tciw_frac_alt_so = i['iw_frac_so'][:]


#---MRI-CGCM3-AMIP Southern Ocean Profile---#

mri_cgcm_tcc_alt_so = i['cf_so'][:] # 0-1
mri_cgcm_tclw_alt_so = i['lw_so'][:] #kg/kg
mri_cgcm_tciw_alt_so = i['iw_so'][:] #kg/kg


mri_cgcm_tcc_temp_so = i['cf_t_so'][:] # 0-1
mri_cgcm_tclw_temp_so = i['lw_t_so'][:] #kg/kg
mri_cgcm_tciw_temp_so = i['iw_t_so'][:] #kg/kg

mri_cgcm_tclw_frac_temp_so = i['lw_frac_t_so'][:]
mri_cgcm_tciw_frac_temp_so = i['lw_frac_t_so'][:]



############################################################################### CESM1-CAM5-AMIP Data

#---CESM1-CAM5-AMIP Global Latitude Data---#

cam5_tcc_lat_g = c['tcc'][:] # 0-1
cam5_tclw_lat_g = c['tclw'][:] #kgm^-2
cam5_tciw_lat_g = c['tciw'][:] #kgm^-2


#---CESM1-CAM5-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
cam5_tclw_frac_lat_g = c['tclw_frac'][:] 
cam5_tciw_frac_lat_g = c['tciw_frac'][:]


#---CESM1-CAM5-AMIP Southern Ocean Latitude Data---#

cam5_tcc_lat_so = cam5_tcc_lat_g[cam5_tcc_lat_g[:,0]>=-70]
cam5_tcc_lat_so = cam5_tcc_lat_so[cam5_tcc_lat_so[:,0]<=-50] # 0-1

cam5_tclw_lat_so = cam5_tclw_lat_g[cam5_tclw_lat_g[:,0]>=-70]
cam5_tclw_lat_so = cam5_tclw_lat_so[cam5_tclw_lat_so[:,0]<=-50] #kgm^-2

cam5_tciw_lat_so = cam5_tciw_lat_g[cam5_tciw_lat_g[:,0]>=-70]
cam5_tciw_lat_so = cam5_tciw_lat_so[cam5_tciw_lat_so[:,0]<=-50] #kgm^-2


#---CESM1-CAM5-AMIP lat-alt contour data---#

cam5_tcc_alt_lat = c['cf_alt_lat'][:] 
cam5_tclw_alt_lat = c['lw_alt_lat'][:] #kg/kg
cam5_tciw_alt_lat = c['iw_alt_lat'][:] #kg/kg
cam5_temp_alt_lat = c['temp_alt_lat'][:] #kg/kg
cam5_lat = c['lat'][:]
cam5_lat = np.hstack(cam5_lat)
cam5_alt = c['alt'][:]
cam5_alt = np.hstack(cam5_alt)
cam5_alt_temp = c['alt_temp'][:]
cam5_alt_temp = np.hstack(cam5_alt_temp)

cam5_tclw_frac_alt_lat = (cam5_tclw_alt_lat / (cam5_tclw_alt_lat + cam5_tciw_alt_lat)) * cam5_tcc_alt_lat
cam5_tciw_frac_alt_lat = (cam5_tciw_alt_lat / (cam5_tclw_alt_lat + cam5_tciw_alt_lat)) * cam5_tcc_alt_lat


#---CESM1-CAM5-AMIP Global Profile---#

cam5_tcc_alt_g = c['cf'][:] # 0-1
cam5_tclw_alt_g = c['lw'][:] #kg/kg
cam5_tciw_alt_g = c['iw'][:] #kg/kg
cam5_temp_alt_g = c['temp'][:] #K
cam5_plevel_alt_g = c['pressure'][:] #hPa

cam5_tcc_temp_g = c['cf_t'][:] # 0-1
cam5_tclw_temp_g = c['lw_t'][:] #kg/kg
cam5_tciw_temp_g = c['iw_t'][:] #kg/kg

cam5_tclw_frac_temp_g = c['lw_frac_t'][:]
cam5_tciw_frac_temp_g = c['iw_frac_t'][:]


#---CESM1-CAM5-AMIP Phase Profile Fractions---#

cam5_tclw_frac_alt_g = c['lw_frac'][:]
cam5_tciw_frac_alt_g = c['iw_frac'][:]
cam5_tclw_frac_alt_so = c['lw_frac_so'][:]
cam5_tciw_frac_alt_so = c['iw_frac_so'][:]


#---CESM1-CAM5-AMIP Southern Ocean Profile---#

cam5_tcc_alt_so = c['cf_so'][:] # 0-1
cam5_tclw_alt_so = c['lw_so'][:] #kg/kg
cam5_tciw_alt_so = c['iw_so'][:] #kg/kg

cam5_tcc_temp_so = c['cf_t_so'][:] # 0-1
cam5_tclw_temp_so = c['lw_t_so'][:] #kg/kg
cam5_tciw_temp_so = c['iw_t_so'][:] #kg/kg

cam5_tclw_frac_temp_so = c['lw_frac_t_so'][:]
cam5_tciw_frac_temp_so = c['lw_frac_t_so'][:]


############################################################################### CESM2-CAM6-AMIP Data

#---CESM2-CAM6-AMIP Global Latitude Data---#

cam6_tcc_lat_g = e['tcc'][:] # 0-1
cam6_tclw_lat_g = e['tclw'][:] #kgm^-2
cam6_tciw_lat_g = e['tciw'][:] #kgm^-2


#---CESM2-CAM6-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
cam6_tclw_frac_lat_g = e['tclw_frac'][:] 
cam6_tciw_frac_lat_g = e['tciw_frac'][:]


#---CESM2-CAM6-AMIP Southern Ocean Latitude Data---#

cam6_tcc_lat_so = cam6_tcc_lat_g[cam6_tcc_lat_g[:,0]>=-70]
cam6_tcc_lat_so = cam6_tcc_lat_so[cam6_tcc_lat_so[:,0]<=-50] # 0-1

cam6_tclw_lat_so = cam6_tclw_lat_g[cam6_tclw_lat_g[:,0]>=-70]
cam6_tclw_lat_so = cam6_tclw_lat_so[cam6_tclw_lat_so[:,0]<=-50] #kgm^-2

cam6_tciw_lat_so = cam6_tciw_lat_g[cam6_tciw_lat_g[:,0]>=-70]
cam6_tciw_lat_so = cam6_tciw_lat_so[cam6_tciw_lat_so[:,0]<=-50] #kgm^-2


#---CESM2-CAM6-AMIP lat-alt contour data---#

cam6_tcc_alt_lat = e['cf_alt_lat'][:] 
cam6_tclw_alt_lat = e['lw_alt_lat'][:] #kg/kg
cam6_tciw_alt_lat = e['iw_alt_lat'][:] #kg/kg
cam6_temp_alt_lat = e['temp_alt_lat'][:] #kg/kg
cam6_lat = e['lat'][:]
cam6_lat = np.hstack(cam6_lat)
cam6_alt = e['alt'][:]
cam6_alt = np.hstack(cam6_alt)

cam6_tclw_frac_alt_lat = (cam6_tclw_alt_lat / (cam6_tclw_alt_lat + cam6_tciw_alt_lat)) * cam6_tcc_alt_lat
cam6_tciw_frac_alt_lat = (cam6_tciw_alt_lat / (cam6_tclw_alt_lat + cam6_tciw_alt_lat)) * cam6_tcc_alt_lat


#---CESM2-CAM6-AMIP Global Profile---#

cam6_tcc_alt_g = e['cf'][:] # 0-1
cam6_tclw_alt_g = e['lw'][:] #kg/kg
cam6_tciw_alt_g = e['iw'][:] #kg/kg
cam6_temp_alt_g = e['temp'][:] #K
cam6_plevel_alt_g = e['pressure'][:] #hPa

cam6_tcc_temp_g = e['cf_t'][:] # 0-1
cam6_tclw_temp_g = e['lw_t'][:] #kg/kg
cam6_tciw_temp_g = e['iw_t'][:] #kg/kg

cam6_tclw_frac_temp_g = e['lw_frac_t'][:]
cam6_tciw_frac_temp_g = e['iw_frac_t'][:]


#---CESM2-CAM6-AMIP Phase Profile Fractions---#

cam6_tclw_frac_alt_g = e['lw_frac'][:]
cam6_tciw_frac_alt_g = e['iw_frac'][:]
cam6_tclw_frac_alt_so = e['lw_frac_so'][:]
cam6_tciw_frac_alt_so = e['iw_frac_so'][:]


#---CESM2-CAM6-AMIP Southern Ocean Profile---#

cam6_tcc_alt_so = e['cf_so'][:] # 0-1
cam6_tclw_alt_so = e['lw_so'][:] #kg/kg
cam6_tciw_alt_so = e['iw_so'][:] #kg/kg

cam6_tcc_temp_so = e['cf_t_so'][:] # 0-1
cam6_tclw_temp_so = e['lw_t_so'][:] #kg/kg
cam6_tciw_temp_so = e['iw_t_so'][:] #kg/kg

cam6_tclw_frac_temp_so = e['lw_frac_t_so'][:]
cam6_tciw_frac_temp_so = e['lw_frac_t_so'][:]

############################################################################### GISS-E2R-AMIP Data

#---GISS-E2R-AMIP Global Latitude Data---#

giss5_tcc_lat_g = j['tcc'][:] # 0-1
giss5_tclw_lat_g = j['tclw'][:] #kgm^-2
giss5_tciw_lat_g = j['tciw'][:] #kgm^-2


#---GISS-E2R-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
giss5_tclw_frac_lat_g = j['tclw_frac'][:] 
giss5_tciw_frac_lat_g = j['tciw_frac'][:]



#---GISS-E2R-AMIP Southern Ocean Latitude Data---#

giss5_tcc_lat_so = giss5_tcc_lat_g[giss5_tcc_lat_g[:,0]>=-70]
giss5_tcc_lat_so = giss5_tcc_lat_so[giss5_tcc_lat_so[:,0]<=-50] # 0-1

giss5_tclw_lat_so = giss5_tclw_lat_g[giss5_tclw_lat_g[:,0]>=-70]
giss5_tclw_lat_so = giss5_tclw_lat_so[giss5_tclw_lat_so[:,0]<=-50] #kgm^-2

giss5_tciw_lat_so = giss5_tciw_lat_g[giss5_tciw_lat_g[:,0]>=-70]
giss5_tciw_lat_so = giss5_tciw_lat_so[giss5_tciw_lat_so[:,0]<=-50] #kgm^-2


#---GISS-E2R-AMIP lat-alt contour data---#

giss5_tcc_alt_lat = j['cf_alt_lat'][:] #kg/kg
giss5_tclw_alt_lat = j['lw_alt_lat'][:] #kg/kg
giss5_tciw_alt_lat = j['iw_alt_lat'][:] #kg/kg
giss5_temp_alt_lat = j['temp_alt_lat'][:] #kg/kg
giss5_lat = j['lat'][:]
giss5_alt = j['alt'][:]
giss5_alt_temp = j['alt_temp'][:]


giss5_tclw_frac_alt_lat = (giss5_tclw_alt_lat / (giss5_tclw_alt_lat + giss5_tciw_alt_lat)) * giss5_tcc_alt_lat
giss5_tciw_frac_alt_lat = (giss5_tciw_alt_lat / (giss5_tclw_alt_lat + giss5_tciw_alt_lat)) * giss5_tcc_alt_lat


#---GISS-E2R-AMIP Global Profile---#

giss5_tcc_alt_g = j['cf'][:] # 0-1
giss5_tclw_alt_g = j['lw'][:] #kg/kg
giss5_tciw_alt_g = j['iw'][:] #kg/kg
giss5_temp_alt_g = j['temp'][:] #K
giss5_plevel_alt_g = j['pressure'][:] #hPa

giss5_tcc_temp_g = j['cf_t'][:] # 0-1
giss5_tclw_temp_g = j['lw_t'][:] #kg/kg
giss5_tciw_temp_g = j['iw_t'][:] #kg/kg

giss5_tclw_frac_temp_g = j['lw_frac_t'][:]
giss5_tciw_frac_temp_g = j['iw_frac_t'][:]


#---GISS-E2R-AMIP Phase Profile Fractions---#

giss5_tclw_frac_alt_g = j['lw_frac'][:]
giss5_tciw_frac_alt_g = j['iw_frac'][:]
giss5_tclw_frac_alt_so = j['lw_frac_so'][:]
giss5_tciw_frac_alt_so = j['iw_frac_so'][:]


#---GISS-E2R-AMIP Southern Ocean Profile---#

giss5_tcc_alt_so = j['cf_so'][:] # 0-1
giss5_tclw_alt_so = j['lw_so'][:] #kg/kg
giss5_tciw_alt_so = j['iw_so'][:] #kg/kg

giss5_tcc_temp_so = j['cf_t_so'][:] # 0-1
giss5_tclw_temp_so = j['lw_t_so'][:] #kg/kg
giss5_tciw_temp_so = j['iw_t_so'][:] #kg/kg

giss5_tclw_frac_temp_so = j['lw_frac_t_so'][:]
giss5_tciw_frac_temp_so = j['lw_frac_t_so'][:]


############################################################################### GISS-E21G-AMIP Data

#---GISS-E21G-AMIP Global Latitude Data---#

giss6_tcc_lat_g = k['tcc'][:] # 0-1
giss6_tclw_lat_g = k['tclw'][:] #kgm^-2
giss6_tciw_lat_g = k['tciw'][:] #kgm^-2


#---GISS-E21G-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
giss6_tclw_frac_lat_g = k['tclw_frac'][:] 
giss6_tciw_frac_lat_g = k['tciw_frac'][:]



#---GISS-E21G-AMIP Southern Ocean Latitude Data---#

giss6_tcc_lat_so = giss6_tcc_lat_g[giss6_tcc_lat_g[:,0]>=-70]
giss6_tcc_lat_so = giss6_tcc_lat_so[giss6_tcc_lat_so[:,0]<=-50] # 0-1

giss6_tclw_lat_so = giss6_tclw_lat_g[giss6_tclw_lat_g[:,0]>=-70]
giss6_tclw_lat_so = giss6_tclw_lat_so[giss6_tclw_lat_so[:,0]<=-50] #kgm^-2

giss6_tciw_lat_so = giss6_tciw_lat_g[giss6_tciw_lat_g[:,0]>=-70]
giss6_tciw_lat_so = giss6_tciw_lat_so[giss6_tciw_lat_so[:,0]<=-50] #kgm^-2


#---GISS-E21G-AMIP lat-alt contour data---#

giss6_tcc_alt_lat = k['cf_alt_lat'][:] #kg/kg
giss6_tclw_alt_lat = k['lw_alt_lat'][:] #kg/kg
giss6_tciw_alt_lat = k['iw_alt_lat'][:] #kg/kg
giss6_temp_alt_lat = k['temp_alt_lat'][:] #kg/kg
giss6_lat = k['lat'][:]
giss6_alt = k['alt'][:]
giss6_alt_temp = k['alt_temp'][:]


giss6_tclw_frac_alt_lat = (giss6_tclw_alt_lat / (giss6_tclw_alt_lat + giss6_tciw_alt_lat)) * giss6_tcc_alt_lat
giss6_tciw_frac_alt_lat = (giss6_tciw_alt_lat / (giss6_tclw_alt_lat + giss6_tciw_alt_lat)) * giss6_tcc_alt_lat


#---GISS-E21G-AMIP Global Profile---#

giss6_tcc_alt_g = k['cf'][:] # 0-1
giss6_tclw_alt_g = k['lw'][:] #kg/kg
giss6_tciw_alt_g = k['iw'][:] #kg/kg
giss6_temp_alt_g = k['temp'][:] #K
giss6_plevel_alt_g = k['pressure'][:] #hPa


giss6_tcc_temp_g = k['cf_t'][:] # 0-1
giss6_tclw_temp_g = k['lw_t'][:] #kg/kg
giss6_tciw_temp_g = k['iw_t'][:] #kg/kg

giss6_tclw_frac_temp_g = k['lw_frac_t'][:]
giss6_tciw_frac_temp_g = k['iw_frac_t'][:]



#---GISS-E21G-AMIP Phase Profile Fractions---#

giss6_tclw_frac_alt_g = k['lw_frac'][:]
giss6_tciw_frac_alt_g = k['iw_frac'][:]
giss6_tclw_frac_alt_so = k['lw_frac_so'][:]
giss6_tciw_frac_alt_so = k['iw_frac_so'][:]


#---GISS-E21G-AMIP Southern Ocean Profile---#

giss6_tcc_alt_so = k['cf_so'][:] # 0-1
giss6_tclw_alt_so = k['lw_so'][:] #kg/kg
giss6_tciw_alt_so = k['iw_so'][:] #kg/kg


giss6_tcc_temp_so = k['cf_t_so'][:] # 0-1
giss6_tclw_temp_so = k['lw_t_so'][:] #kg/kg
giss6_tciw_temp_so = k['iw_t_so'][:] #kg/kg

giss6_tclw_frac_temp_so = k['lw_frac_t_so'][:]
giss6_tciw_frac_temp_so = k['lw_frac_t_so'][:]


############################################################################### IPSL-CM5A-LR-AMIP Data

#---IPSL-CM5A-LR-AMIP Global Latitude Data---#

ipsl5_tcc_lat_g = l['tcc'][:] # 0-1
ipsl5_tclw_lat_g = l['tclw'][:] #kgm^-2
ipsl5_tciw_lat_g = l['tciw'][:] #kgm^-2


#---IPSL-CM5A-LR-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
ipsl5_tclw_frac_lat_g = l['tclw_frac'][:] 
ipsl5_tciw_frac_lat_g = l['tciw_frac'][:]



#---IPSL-CM5A-LR-AMIP Southern Ocean Latitude Data---#

ipsl5_tcc_lat_so = ipsl5_tcc_lat_g[ipsl5_tcc_lat_g[:,0]>=-70]
ipsl5_tcc_lat_so = ipsl5_tcc_lat_so[ipsl5_tcc_lat_so[:,0]<=-50] # 0-1

ipsl5_tclw_lat_so = ipsl5_tclw_lat_g[ipsl5_tclw_lat_g[:,0]>=-70]
ipsl5_tclw_lat_so = ipsl5_tclw_lat_so[ipsl5_tclw_lat_so[:,0]<=-50] #kgm^-2

ipsl5_tciw_lat_so = ipsl5_tciw_lat_g[ipsl5_tciw_lat_g[:,0]>=-70]
ipsl5_tciw_lat_so = ipsl5_tciw_lat_so[ipsl5_tciw_lat_so[:,0]<=-50] #kgm^-2


#---IPSL-CM5A-LR-AMIP lat-alt contour data---#

ipsl5_tcc_alt_lat = l['cf_alt_lat'][:] #kg/kg
ipsl5_tclw_alt_lat = l['lw_alt_lat'][:] #kg/kg
ipsl5_tciw_alt_lat = l['iw_alt_lat'][:] #kg/kg
ipsl5_temp_alt_lat = l['temp_alt_lat'][:] #kg/kg
ipsl5_lat = l['lat'][:]
ipsl5_alt = l['alt'][:]
ipsl5_alt_temp = l['alt_temp'][:]


ipsl5_tclw_frac_alt_lat = (ipsl5_tclw_alt_lat / (ipsl5_tclw_alt_lat + ipsl5_tciw_alt_lat)) * ipsl5_tcc_alt_lat
ipsl5_tciw_frac_alt_lat = (ipsl5_tciw_alt_lat / (ipsl5_tclw_alt_lat + ipsl5_tciw_alt_lat)) * ipsl5_tcc_alt_lat


#---IPSL-CM5A-LR-AMIP Global Profile---#

ipsl5_tcc_alt_g = l['cf'][:] # 0-1
ipsl5_tclw_alt_g = l['lw'][:] #kg/kg
ipsl5_tciw_alt_g = l['iw'][:] #kg/kg
ipsl5_temp_alt_g = l['temp'][:] #K
ipsl5_plevel_alt_g = l['pressure'][:] #hPa

ipsl5_tcc_temp_g = l['cf_t'][:] # 0-1
ipsl5_tclw_temp_g = l['lw_t'][:] #kg/kg
ipsl5_tciw_temp_g = l['iw_t'][:] #kg/kg

ipsl5_tclw_frac_temp_g = l['lw_frac_t'][:]
ipsl5_tciw_frac_temp_g = l['iw_frac_t'][:]


#---IPSL-CM5A-LR-AMIP Phase Profile Fractions---#

ipsl5_tclw_frac_alt_g = l['lw_frac'][:]
ipsl5_tciw_frac_alt_g = l['iw_frac'][:]
ipsl5_tclw_frac_alt_so = l['lw_frac_so'][:]
ipsl5_tciw_frac_alt_so = l['iw_frac_so'][:]


#---IPSL-CM5A-LR-AMIP Southern Ocean Profile---#

ipsl5_tcc_alt_so = l['cf_so'][:] # 0-1
ipsl5_tclw_alt_so = l['lw_so'][:] #kg/kg
ipsl5_tciw_alt_so = l['iw_so'][:] #kg/kg

ipsl5_tcc_temp_so = l['cf_t_so'][:] # 0-1
ipsl5_tclw_temp_so = l['lw_t_so'][:] #kg/kg
ipsl5_tciw_temp_so = l['iw_t_so'][:] #kg/kg

ipsl5_tclw_frac_temp_so = l['lw_frac_t_so'][:]
ipsl5_tciw_frac_temp_so = l['lw_frac_t_so'][:]


############################################################################### IPSL-CM6A-LR-AMIP Data

#---IPSL-CM6A-LR-AMIP Global Latitude Data---#

ipsl6_tcc_lat_g = m['tcc'][:] # 0-1
ipsl6_tclw_lat_g = m['tclw'][:] #kgm^-2
ipsl6_tciw_lat_g = m['tciw'][:] #kgm^-2


#---IPSL-CM6A-LR-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
ipsl6_tclw_frac_lat_g = m['tclw_frac'][:] 
ipsl6_tciw_frac_lat_g = m['tciw_frac'][:]



#---IPSL-CM6A-LR-AMIP Southern Ocean Latitude Data---#

ipsl6_tcc_lat_so = ipsl6_tcc_lat_g[ipsl6_tcc_lat_g[:,0]>=-70]
ipsl6_tcc_lat_so = ipsl6_tcc_lat_so[ipsl6_tcc_lat_so[:,0]<=-50] # 0-1

ipsl6_tclw_lat_so = ipsl6_tclw_lat_g[ipsl6_tclw_lat_g[:,0]>=-70]
ipsl6_tclw_lat_so = ipsl6_tclw_lat_so[ipsl6_tclw_lat_so[:,0]<=-50] #kgm^-2

ipsl6_tciw_lat_so = ipsl6_tciw_lat_g[ipsl6_tciw_lat_g[:,0]>=-70]
ipsl6_tciw_lat_so = ipsl6_tciw_lat_so[ipsl6_tciw_lat_so[:,0]<=-50] #kgm^-2


#---IPSL-CM6A-LR-AMIP lat-alt contour data---#

ipsl6_tcc_alt_lat = m['cf_alt_lat'][:] #kg/kg
ipsl6_tclw_alt_lat = m['lw_alt_lat'][:] #kg/kg
ipsl6_tciw_alt_lat = m['iw_alt_lat'][:] #kg/kg
ipsl6_temp_alt_lat = m['temp_alt_lat'][:] #kg/kg
ipsl6_lat = m['lat'][:]
ipsl6_alt = m['alt'][:]
ipsl6_alt_temp = m['alt_temp'][:]


ipsl6_tclw_frac_alt_lat = (ipsl6_tclw_alt_lat / (ipsl6_tclw_alt_lat + ipsl6_tciw_alt_lat)) * ipsl6_tcc_alt_lat
ipsl6_tciw_frac_alt_lat = (ipsl6_tciw_alt_lat / (ipsl6_tclw_alt_lat + ipsl6_tciw_alt_lat)) * ipsl6_tcc_alt_lat


#---IPSL-CM6A-LR-AMIP Global Profile---#

ipsl6_tcc_alt_g = m['cf'][:] # 0-1
ipsl6_tclw_alt_g = m['lw'][:] #kg/kg
ipsl6_tciw_alt_g = m['iw'][:] #kg/kg
ipsl6_temp_alt_g = m['temp'][:] #K
ipsl6_plevel_alt_g = m['pressure'][:] #hPa


ipsl6_tcc_temp_g = m['cf_t'][:] # 0-1
ipsl6_tclw_temp_g = m['lw_t'][:] #kg/kg
ipsl6_tciw_temp_g = m['iw_t'][:] #kg/kg

ipsl6_tclw_frac_temp_g = m['lw_frac_t'][:]
ipsl6_tciw_frac_temp_g = m['iw_frac_t'][:]



#---IPSL-CM6A-LR-AMIP Phase Profile Fractions---#

ipsl6_tclw_frac_alt_g = m['lw_frac'][:]
ipsl6_tciw_frac_alt_g = m['iw_frac'][:]
ipsl6_tclw_frac_alt_so = m['lw_frac_so'][:]
ipsl6_tciw_frac_alt_so = m['iw_frac_so'][:]


#---IPSL-CM6A-LR-AMIP Southern Ocean Profile---#

ipsl6_tcc_alt_so = m['cf_so'][:] # 0-1
ipsl6_tclw_alt_so = m['lw_so'][:] #kg/kg
ipsl6_tciw_alt_so = m['iw_so'][:] #kg/kg


ipsl6_tcc_temp_so = m['cf_t_so'][:] # 0-1
ipsl6_tclw_temp_so = m['lw_t_so'][:] #kg/kg
ipsl6_tciw_temp_so = m['iw_t_so'][:] #kg/kg

ipsl6_tclw_frac_temp_so = m['lw_frac_t_so'][:]
ipsl6_tciw_frac_temp_so = m['lw_frac_t_so'][:]

############################################################################### MIROC5-AMIP Data

#---MIROC5-AMIP Global Latitude Data---#

miroc5_tcc_lat_g = n['tcc'][:] # 0-1
miroc5_tclw_lat_g = n['tclw'][:] #kgm^-2
miroc5_tciw_lat_g = n['tciw'][:] #kgm^-2


#---MIROC5-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
miroc5_tclw_frac_lat_g = n['tclw_frac'][:] 
miroc5_tciw_frac_lat_g = n['tciw_frac'][:]



#---MIROC5-AMIP Southern Ocean Latitude Data---#

miroc5_tcc_lat_so = miroc5_tcc_lat_g[miroc5_tcc_lat_g[:,0]>=-70]
miroc5_tcc_lat_so = miroc5_tcc_lat_so[miroc5_tcc_lat_so[:,0]<=-50] # 0-1

miroc5_tclw_lat_so = miroc5_tclw_lat_g[miroc5_tclw_lat_g[:,0]>=-70]
miroc5_tclw_lat_so = miroc5_tclw_lat_so[miroc5_tclw_lat_so[:,0]<=-50] #kgm^-2

miroc5_tciw_lat_so = miroc5_tciw_lat_g[miroc5_tciw_lat_g[:,0]>=-70]
miroc5_tciw_lat_so = miroc5_tciw_lat_so[miroc5_tciw_lat_so[:,0]<=-50] #kgm^-2


#---MIROC5-AMIP lat-alt contour data---#

miroc5_tcc_alt_lat = n['cf_alt_lat'][:] #kg/kg
miroc5_tclw_alt_lat = n['lw_alt_lat'][:] #kg/kg
miroc5_tciw_alt_lat = n['iw_alt_lat'][:] #kg/kg
miroc5_temp_alt_lat = n['temp_alt_lat'][:] #kg/kg
miroc5_lat = n['lat'][:]
miroc5_alt = n['alt'][:]
miroc5_alt_temp = n['alt_temp'][:]


miroc5_tclw_frac_alt_lat = (miroc5_tclw_alt_lat / (miroc5_tclw_alt_lat + miroc5_tciw_alt_lat)) * miroc5_tcc_alt_lat
miroc5_tciw_frac_alt_lat = (miroc5_tciw_alt_lat / (miroc5_tclw_alt_lat + miroc5_tciw_alt_lat)) * miroc5_tcc_alt_lat


#---MIROC5-AMIP Global Profile---#

miroc5_tcc_alt_g = n['cf'][:] # 0-1
miroc5_tclw_alt_g = n['lw'][:] #kg/kg
miroc5_tciw_alt_g = n['iw'][:] #kg/kg
miroc5_temp_alt_g = n['temp'][:] #K
miroc5_plevel_alt_g = n['pressure'][:] #hPa

miroc5_tcc_temp_g = n['cf_t'][:] # 0-1
miroc5_tclw_temp_g = n['lw_t'][:] #kg/kg
miroc5_tciw_temp_g = n['iw_t'][:] #kg/kg

miroc5_tclw_frac_temp_g = n['lw_frac_t'][:]
miroc5_tciw_frac_temp_g = n['iw_frac_t'][:]


#---MIROC5-AMIP Phase Profile Fractions---#

miroc5_tclw_frac_alt_g = n['lw_frac'][:]
miroc5_tciw_frac_alt_g = n['iw_frac'][:]
miroc5_tclw_frac_alt_so = n['lw_frac_so'][:]
miroc5_tciw_frac_alt_so = n['iw_frac_so'][:]


#---MIROC5-AMIP Southern Ocean Profile---#

miroc5_tcc_alt_so = n['cf_so'][:] # 0-1
miroc5_tclw_alt_so = n['lw_so'][:] #kg/kg
miroc5_tciw_alt_so = n['iw_so'][:] #kg/kg

miroc5_tcc_temp_so = n['cf_t_so'][:] # 0-1
miroc5_tclw_temp_so = n['lw_t_so'][:] #kg/kg
miroc5_tciw_temp_so = n['iw_t_so'][:] #kg/kg

miroc5_tclw_frac_temp_so = n['lw_frac_t_so'][:]
miroc5_tciw_frac_temp_so = n['lw_frac_t_so'][:]


############################################################################### MIROC6-AMIP Data

#---MIROC6-AMIP Global Latitude Data---#

miroc6_tcc_lat_g = o['tcc'][:] # 0-1
miroc6_tclw_lat_g = o['tclw'][:] #kgm^-2
miroc6_tciw_lat_g = o['tciw'][:] #kgm^-2


#---MIROC6-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
miroc6_tclw_frac_lat_g = o['tclw_frac'][:] 
miroc6_tciw_frac_lat_g = o['tciw_frac'][:]



#---MIROC6-AMIP Southern Ocean Latitude Data---#

miroc6_tcc_lat_so = miroc6_tcc_lat_g[miroc6_tcc_lat_g[:,0]>=-70]
miroc6_tcc_lat_so = miroc6_tcc_lat_so[miroc6_tcc_lat_so[:,0]<=-50] # 0-1

miroc6_tclw_lat_so = miroc6_tclw_lat_g[miroc6_tclw_lat_g[:,0]>=-70]
miroc6_tclw_lat_so = miroc6_tclw_lat_so[miroc6_tclw_lat_so[:,0]<=-50] #kgm^-2

miroc6_tciw_lat_so = miroc6_tciw_lat_g[miroc6_tciw_lat_g[:,0]>=-70]
miroc6_tciw_lat_so = miroc6_tciw_lat_so[miroc6_tciw_lat_so[:,0]<=-50] #kgm^-2


#---MIROC6-AMIP lat-alt contour data---#

miroc6_tcc_alt_lat = o['cf_alt_lat'][:] #kg/kg
miroc6_tclw_alt_lat = o['lw_alt_lat'][:] #kg/kg
miroc6_tciw_alt_lat = o['iw_alt_lat'][:] #kg/kg
miroc6_temp_alt_lat = o['temp_alt_lat'][:] #kg/kg
miroc6_lat = o['lat'][:]
miroc6_alt = o['alt'][:]
miroc6_alt_temp = o['alt_temp'][:]


miroc6_tclw_frac_alt_lat = (miroc6_tclw_alt_lat / (miroc6_tclw_alt_lat + miroc6_tciw_alt_lat)) * miroc6_tcc_alt_lat
miroc6_tciw_frac_alt_lat = (miroc6_tciw_alt_lat / (miroc6_tclw_alt_lat + miroc6_tciw_alt_lat)) * miroc6_tcc_alt_lat


#---MIROC6-AMIP Global Profile---#

miroc6_tcc_alt_g = o['cf'][:] # 0-1
miroc6_tclw_alt_g = o['lw'][:] #kg/kg
miroc6_tciw_alt_g = o['iw'][:] #kg/kg
miroc6_temp_alt_g = o['temp'][:] #K
miroc6_plevel_alt_g = o['pressure'][:] #hPa


miroc6_tcc_temp_g = o['cf_t'][:] # 0-1
miroc6_tclw_temp_g = o['lw_t'][:] #kg/kg
miroc6_tciw_temp_g = o['iw_t'][:] #kg/kg

miroc6_tclw_frac_temp_g = o['lw_frac_t'][:]
miroc6_tciw_frac_temp_g = o['iw_frac_t'][:]



#---MIROC6-AMIP Phase Profile Fractions---#

miroc6_tclw_frac_alt_g = o['lw_frac'][:]
miroc6_tciw_frac_alt_g = o['iw_frac'][:]
miroc6_tclw_frac_alt_so = o['lw_frac_so'][:]
miroc6_tciw_frac_alt_so = o['iw_frac_so'][:]


#---MIROC6-AMIP Southern Ocean Profile---#

miroc6_tcc_alt_so = o['cf_so'][:] # 0-1
miroc6_tclw_alt_so = o['lw_so'][:] #kg/kg
miroc6_tciw_alt_so = o['iw_so'][:] #kg/kg


miroc6_tcc_temp_so = o['cf_t_so'][:] # 0-1
miroc6_tclw_temp_so = o['lw_t_so'][:] #kg/kg
miroc6_tciw_temp_so = o['iw_t_so'][:] #kg/kg

miroc6_tclw_frac_temp_so = o['lw_frac_t_so'][:]
miroc6_tciw_frac_temp_so = o['lw_frac_t_so'][:]

############################################################################### End importing Data

end = time.time()
print('Importing data took:', end - start, 's')


"""

os.chdir('c:/Users/toha006/University/University/MSc/Models/Images/')
"""
os.chdir('E:/University/University/MSc/Models/Images/Meeting 17.7 - 12 models')

############################################################################### Temperature Profiles


#---Plot Global Liquid Fraction with Temperature---#


"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tclw_frac_temp_g[17:,0],(ecmwf_tclw_frac_temp_g[17:,1]/ecmwf_tcc_temp_g[17:,1]), '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tclw_frac_temp_g[:19,0],(gfdl_hiram_tclw_frac_temp_g[:19,1]/gfdl_hiram_tcc_temp_g[:19,1]), '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_g[:22,0],(mri_cgcm_tclw_frac_temp_g[:22,1]/mri_cgcm_tcc_temp_g[:22,1]), '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_temp_g[:14,0],(cam5_tclw_frac_temp_g[:14,1]/cam5_tcc_temp_g[:14,1]), '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tclw_frac_temp_g[:22,0],(miroc5_tclw_frac_temp_g[:22,1]/miroc5_tcc_temp_g[:22,1]), '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tclw_frac_temp_g[:14,0],(ipsl5_tclw_frac_temp_g[:14,1]/ipsl5_tcc_temp_g[:14,1]), '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tclw_frac_temp_g[:18,0],(giss5_tclw_frac_temp_g[:18,1]/giss5_tcc_temp_g[:18,1]), '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(ecmwf_tclw_frac_temp_g[17:,0],(ecmwf_tclw_frac_temp_g[17:,1]/ecmwf_tcc_temp_g[17:,1]), '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tclw_frac_temp_g[:20,0],(gfdl4_tclw_frac_temp_g[:20,1]/gfdl4_tcc_temp_g[:20,1]), '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_g[:28,0],(mri_tclw_frac_temp_g[:28,1]/mri_tcc_temp_g[:28,1]), '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_g[16:,0],(cam6_tclw_frac_temp_g[16:,1]/cam6_tcc_temp_g[16:,1]), '-c', label='CMIP6-CESM2-CAM6-AMIP')
ax2.plot(miroc6_tclw_frac_temp_g[:24,0],(miroc6_tclw_frac_temp_g[:24,1]/miroc6_tcc_temp_g[:24,1]), '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tclw_frac_temp_g[:39,0],(ipsl6_tclw_frac_temp_g[:39,1]/ipsl6_tcc_temp_g[:39,1]), '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tclw_frac_temp_g[:19,0],(giss6_tclw_frac_temp_g[:19,1]/giss6_tcc_temp_g[:19,1]), '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')



ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Fraction')

#ax1.set_title('2001 to 2005 Global Liquid Fraction vs Temperature')

ax1.grid(True)
ax2.grid(True)
ax1.text(205, 1, 'a)')
ax2.text(205, 1, 'b)')
plt.savefig("2001_2005_global_liquid_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Global Liquid Cloud Fraction with Temperature---#


"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tclw_frac_temp_g[18:,0],ecmwf_tclw_frac_temp_g[18:,1], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tclw_frac_temp_g[:18,0],gfdl_hiram_tclw_frac_temp_g[:18,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_g[:22,0],mri_cgcm_tclw_frac_temp_g[:22,1], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_temp_g[:13,0],cam5_tclw_frac_temp_g[:13,1], '-c', label='CMIP5-CESM1-CAM5-AMIP')

ax2.plot(ecmwf_tclw_frac_temp_g[17:,0],ecmwf_tclw_frac_temp_g[17:,1], '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tclw_frac_temp_g[:20,0],gfdl4_tclw_frac_temp_g[:20,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_g[:28,0],mri_tclw_frac_temp_g[:28,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_g[16:,0],cam6_tclw_frac_temp_g[16:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')



ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Cloud Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Cloud Fraction')

#ax1.set_title('2001 to 2005 Global Liquid Cloud Fraction vs Temperature')

ax1.grid(True)
ax2.grid(True)
ax1.text(208, 0.3, 'a)')
ax2.text(208, 0.3, 'b)')
plt.savefig("2001_2005_global_liquid_frac_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Southern Ocean Liquid Fraction with Temperature---#


"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tclw_frac_temp_so[18:,0],(ecmwf_tclw_frac_temp_so[18:,1]/ecmwf_tcc_temp_so[18:,1]), '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tclw_frac_temp_so[:19,0],(gfdl_hiram_tclw_frac_temp_so[:19,1]/gfdl_hiram_tcc_temp_so[:19,1]), '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_so[:22,0],(mri_cgcm_tclw_frac_temp_so[:22,1]/mri_cgcm_tcc_temp_so[:22,1]), '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_temp_so[:14,0],(cam5_tclw_frac_temp_so[:14,1]/cam5_tcc_temp_so[:14,1]), '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tclw_frac_temp_so[:22,0],(miroc5_tclw_frac_temp_so[:22,1]/miroc5_tcc_temp_so[:22,1]), '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tclw_frac_temp_so[:14,0],(ipsl5_tclw_frac_temp_so[:14,1]/ipsl5_tcc_temp_so[:14,1]), '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tclw_frac_temp_so[:18,0],(giss5_tclw_frac_temp_so[:18,1]/giss5_tcc_temp_so[:18,1]), '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(ecmwf_tclw_frac_temp_so[18:,0],(ecmwf_tclw_frac_temp_so[18:,1]/ecmwf_tcc_temp_so[18:,1]), '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tclw_frac_temp_so[:20,0],(gfdl4_tclw_frac_temp_so[:20,1]/gfdl4_tcc_temp_so[:20,1]), '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_so[:28,0],(mri_tclw_frac_temp_so[:28,1]/mri_tcc_temp_so[:28,1]), '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_so[16:,0],(cam6_tclw_frac_temp_so[16:,1]/cam6_tcc_temp_so[16:,1]), '-c', label='CMIP6-CESM2-CAM6-AMIP')
ax2.plot(miroc6_tclw_frac_temp_so[:24,0],(miroc6_tclw_frac_temp_so[:24,1]/miroc6_tcc_temp_so[:24,1]), '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tclw_frac_temp_so[:39,0],(ipsl6_tclw_frac_temp_so[:39,1]/ipsl6_tcc_temp_so[:39,1]), '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tclw_frac_temp_so[:19,0],(giss6_tclw_frac_temp_so[:19,1]/giss6_tcc_temp_so[:19,1]), '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Cloud Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Cloud Fraction')

#ax1.set_title('2001 to 2005 Southern Ocean Liquid Cloud Fraction vs Temperature')

ax1.text(200, 1, 'a)')
ax2.text(200, 1, 'b)')
ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_so_liquid_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot Southern Ocean Liquid Cloud Fraction with Temperature---#


"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))


ax1.plot(ecmwf_tclw_frac_temp_so[18:,0],ecmwf_tclw_frac_temp_so[18:,1], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tclw_frac_temp_so[:18,0],gfdl_hiram_tclw_frac_temp_so[:18,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_so[:22,0],mri_cgcm_tclw_frac_temp_so[:22,1], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_temp_so[:13,0],cam5_tclw_frac_temp_so[:13,1], '-c', label='CMIP5-CESM1-CAM5-AMIP')

ax2.plot(ecmwf_tclw_frac_temp_so[17:,0],ecmwf_tclw_frac_temp_so[17:,1], '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tclw_frac_temp_so[:20,0],gfdl4_tclw_frac_temp_so[:20,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_so[:28,0],mri_tclw_frac_temp_so[:28,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_so[16:,0],cam6_tclw_frac_temp_so[16:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Cloud Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Cloud Fraction')

#ax1.set_title('2001 to 2005 Southern Ocean Liquid Cloud Fraction vs Temperature')

ax1.text(210, 0.3, 'a)')
ax2.text(210, 0.3, 'b)')
ax1.grid(True)
ax2.grid(True)
#plt.savefig("2001_2005_so_liquid_frac_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""

############################################################################### Global Latitude Plots

#---Plot Global Cloud Fraction with Latitude---#

"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tcc_lat_g[:,0],gfdl_hiram_tcc_lat_g[:,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tcc_lat_g[:,0],mri_cgcm_tcc_lat_g[:,1], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tcc_lat_g[:,0],cam5_tcc_lat_g[:,1], '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tcc_lat_g[:,0],miroc5_tcc_lat_g[:,1], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tcc_lat_g[:,0],ipsl5_tcc_lat_g[:,1], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tcc_lat_g[:,0],giss5_tcc_lat_g[:,1], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tcc_lat_g[:,0],gfdl4_tcc_lat_g[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_lat_g[:,0],mri_tcc_lat_g[:,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tcc_lat_g[:,0],cam6_tcc_lat_g[:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')
ax2.plot(miroc6_tcc_lat_g[:,0],miroc6_tcc_lat_g[:,1], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tcc_lat_g[:,0],ipsl6_tcc_lat_g[:,1], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tcc_lat_g[:,0],giss6_tcc_lat_g[:,1], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Cloud Fraction')
ax2.set_xlabel('Latitude')

#ax1.set_title ('2001 to 2005 Global Cloud Fraction vs Latitude')
ax1.text(-140, 1, 'a)')
ax2.text(-140, 1, 'b)')
ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_tcc_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot Global Comparison Liquid Cloud Fractions with Latitude Double Plot---#
"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(gfdl_hiram_tclw_frac_lat_g[:,0],gfdl_hiram_tclw_frac_lat_g[:,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_lat_g[:,0],mri_cgcm_tclw_frac_lat_g[:,1], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_lat_g[:,0],cam5_tclw_frac_lat_g[:,1], '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tclw_lat_g[:,0],miroc5_tclw_lat_g[:,1], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tclw_lat_g[:,0],ipsl5_tclw_lat_g[:,1], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tclw_lat_g[:,0],giss5_tclw_lat_g[:,1], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(gfdl4_tclw_frac_lat_g[:,0],gfdl4_tclw_frac_lat_g[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_lat_g[:,0],mri_tclw_frac_lat_g[:,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_lat_g[:,0],cam6_tclw_frac_lat_g[:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')
ax2.plot(miroc6_tclw_lat_g[:,0],miroc6_tclw_lat_g[:,1], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tclw_lat_g[:,0],ipsl6_tclw_lat_g[:,1], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tclw_lat_g[:,0],giss6_tclw_lat_g[:,1], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax2.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Liquid Water Fraction')
ax2.set_ylabel('Cloud Liquid Water Fraction')

ax1.set_title('2001 to 2005 Global Liquid Fractions of Clouds vs Latitude')

ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_tclw_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Global Comparison Ice Cloud Fractions with Latitude---#

"""
fig, ax = plt.subplots()

ax.plot(gfdl_hiram_tciw_frac_lat_g[:,0],gfdl_hiram_tciw_frac_lat_g[:,1], '--k', label='CMIP5-GFDL-HIRAM-AMIP')
ax.plot(gfdl4_tciw_frac_lat_g[:,0],gfdl4_tciw_frac_lat_g[:,1], '--g', label='CMIP6-GFDL-AM4-AMIP')
ax.plot(mri_cgcm_tciw_frac_lat_g[:,0],mri_cgcm_tciw_frac_lat_g[:,1], '--y', label='CMIP5-MRI_CGCM3-AMIP')
ax.plot(mri_tciw_frac_lat_g[:,0],mri_tciw_frac_lat_g[:,1], '--m', label='CMIP6-MRI_ESM2-AMIP')
ax.plot(cam5_tciw_frac_lat_g[:,0],cam5_tciw_frac_lat_g[:,1], '-b', label='CMIP5-CESM1-CAM5-AMIP')
ax.plot(cam6_tciw_frac_lat_g[:,0],cam6_tciw_frac_lat_g[:,1], '--c', label='CMIP6-CESM2-CAM6-AMIP')


ax.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Ice Water Fraction')

plt.title('2001 to 2005 Global Ice Fractions of Clouds vs Latitude')

plt.grid(True)
plt.savefig("2001_2005_tciw_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""




############################################################################### Global Altitude Plots

#---Plot Global Cloud Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tcc_alt_g[9:,1],ecmwf_tcc_alt_g[9:,0], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tcc_alt_g[:23,1],gfdl_hiram_tcc_alt_g[:23,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tcc_alt_g[:25,1],mri_cgcm_tcc_alt_g[:25,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tcc_alt_g[:16,1],cam5_tcc_alt_g[:16,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tcc_alt_g[:,1],miroc5_tcc_alt_g[:,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tcc_alt_g[:,1],ipsl5_tcc_alt_g[:,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tcc_alt_g[:,1],giss5_tcc_alt_g[:,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(ecmwf_tcc_alt_g[9:,1],ecmwf_tcc_alt_g[9:,0], '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tcc_alt_g[:23,1],gfdl4_tcc_alt_g[:23,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_alt_g[:42,1],mri_tcc_alt_g[:42,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tcc_alt_g[10:,1],cam6_tcc_alt_g[10:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tcc_alt_g[:,1],miroc6_tcc_alt_g[:,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tcc_alt_g[:,1],ipsl6_tcc_alt_g[:,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tcc_alt_g[:,1],giss6_tcc_alt_g[:,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Cloud Fraction')

#ax1.set_title('2001 to 2005 Global Cloud Fraction vs Altitude')

ax1.set_xlim(0, 0.35)
ax2.set_xlim(0, 0.35)

ax1.text(-0.1, 19, 'a)')
ax2.text(-0.05, 19, 'b)')

ax1.grid(True)
ax2.grid(True)

plt.savefig("2001_2005_tcc_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Cloud Liquid Water Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tclw_frac_alt_g[18:,1],ecmwf_tclw_frac_alt_g[18:,0], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tclw_frac_alt_g[:18,1],gfdl_hiram_tclw_frac_alt_g[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_alt_g[:19,1],mri_cgcm_tclw_frac_alt_g[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_alt_g[:13,1],cam5_tclw_frac_alt_g[:13,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tclw_frac_alt_g[:,1],miroc5_tclw_frac_alt_g[:,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tclw_frac_alt_g[:,1],ipsl5_tclw_frac_alt_g[:,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tclw_frac_alt_g[:,1],giss5_tclw_frac_alt_g[:,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(ecmwf_tclw_frac_alt_g[18:,1],ecmwf_tclw_frac_alt_g[18:,0], '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tclw_frac_alt_g[:19,1],gfdl4_tclw_frac_alt_g[:19,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_alt_g[:26,1],mri_tclw_frac_alt_g[:26,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_frac_alt_g[17:,1],cam6_tclw_frac_alt_g[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tclw_frac_alt_g[:,1],miroc6_tclw_frac_alt_g[:,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tclw_frac_alt_g[:,1],ipsl6_tclw_frac_alt_g[:,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tclw_frac_alt_g[:,1],giss6_tclw_frac_alt_g[:,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Liquid Water Fraction')
ax2.set_xlabel('Cloud Liquid Water Fraction')

#ax1.set_title('2001 to 2005 Global Cloud Liquid Water Fraction vs Altitude')

ax1.set_xlim(0, 0.35)
ax2.set_xlim(0, 0.35)

ax1.text(-0.1, 9, 'a)')
ax2.text(-0.05, 9, 'b)')
ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_tclw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Cloud Liquid Water Content Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tclw_alt_g[18:,1]*1000,ecmwf_tclw_alt_g[18:,0], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tclw_alt_g[:18,1]*1000,gfdl_hiram_tclw_alt_g[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_alt_g[:19,1]*1000,mri_cgcm_tclw_alt_g[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_alt_g[:13,1]*1000,cam5_tclw_alt_g[:13,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tclw_alt_g[:,1]*1000,miroc5_tclw_alt_g[:,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tclw_alt_g[:,1]*1000,ipsl5_tclw_alt_g[:,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tclw_alt_g[:,1]*1000,giss5_tclw_alt_g[:,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(ecmwf_tclw_alt_g[18:,1]*1000,ecmwf_tclw_alt_g[18:,0], '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tclw_alt_g[:19,1]*1000,gfdl4_tclw_alt_g[:19,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_alt_g[:26,1]*1000,mri_tclw_alt_g[:26,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_alt_g[17:,1]*1000,cam6_tclw_alt_g[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tclw_alt_g[:,1]*1000,miroc6_tclw_frac_alt_g[:,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tclw_alt_g[:,1]*1000,ipsl6_tclw_alt_g[:,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tclw_alt_g[:,1]*1000,giss6_tclw_alt_g[:,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Liquid Water Content (g/kg)')
ax2.set_xlabel('Liquid Water Content (g/kg)')

#ax1.set_title('2001 to 2005 Global Cloud Liquid Water Fraction vs Altitude')

ax1.set_xlim(0, 0.027)
ax2.set_xlim(0, 0.027)

ax1.text(-0.005, 9, 'a)')
ax2.text(-0.004, 9, 'b)')
ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_tclw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""


############################################################################### Southern Ocean Altitude Plots

#---Plot Southern Ocean Cloud Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tcc_alt_so[13:,1],ecmwf_tcc_alt_so[13:,0], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tcc_alt_so[:23,1],gfdl_hiram_tcc_alt_so[:23,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tcc_alt_so[:25,1],mri_cgcm_tcc_alt_so[:25,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tcc_alt_so[:16,1],cam5_tcc_alt_so[:16,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tcc_alt_so[:,1],miroc5_tcc_alt_so[:,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tcc_alt_so[:,1],ipsl5_tcc_alt_so[:,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tcc_alt_so[:,1],giss5_tcc_alt_so[:,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(ecmwf_tcc_alt_so[9:,1],ecmwf_tcc_alt_so[9:,0], '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tcc_alt_so[:23,1],gfdl4_tcc_alt_so[:23,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_alt_so[:42,1],mri_tcc_alt_so[:42,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tcc_alt_so[10:,1],cam6_tcc_alt_so[10:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tcc_alt_so[:,1],miroc6_tcc_alt_so[:,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tcc_alt_so[:,1],ipsl6_tcc_alt_so[:,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tcc_alt_so[:,1],giss6_tcc_alt_so[:,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Cloud Fraction')

#ax1.set_title('2001 to 2005 Southern Ocean Cloud Fraction vs Altitude')

ax1.set_xlim(0, 0.65)
ax2.set_xlim(0, 0.65)

ax1.text(-0.11, 19, 'a)')
ax2.text(-0.09, 19, 'b)')

ax1.grid(True)
ax2.grid(True)

plt.savefig("2001_2005_tcc_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Cloud Liquid Water Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tclw_frac_alt_so[18:,1],ecmwf_tclw_frac_alt_so[18:,0], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tclw_frac_alt_so[:18,1],gfdl_hiram_tclw_frac_alt_so[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_alt_so[:19,1],mri_cgcm_tclw_frac_alt_so[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_alt_so[:13,1],cam5_tclw_frac_alt_so[:13,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tclw_frac_alt_so[:,1],miroc5_tclw_frac_alt_so[:,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tclw_frac_alt_so[:,1],ipsl5_tclw_frac_alt_so[:,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tclw_frac_alt_so[:,1],giss5_tclw_frac_alt_so[:,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(ecmwf_tclw_frac_alt_so[18:,1],ecmwf_tclw_frac_alt_so[18:,0], '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tclw_frac_alt_so[:19,1],gfdl4_tclw_frac_alt_so[:19,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_alt_so[:26,1],mri_tclw_frac_alt_so[:26,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_frac_alt_so[17:,1],cam6_tclw_frac_alt_so[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tclw_frac_alt_so[:,1],miroc6_tclw_frac_alt_so[:,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tclw_frac_alt_so[:,1],ipsl6_tclw_frac_alt_so[:,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tclw_frac_alt_so[:,1],giss6_tclw_frac_alt_so[:,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Liquid Water Fraction')
ax2.set_xlabel('Cloud Liquid Water Fraction')

#ax1.set_title('2001 to 2005 Southern Ocean Cloud Liquid Water Fraction vs Altitude')

ax1.set_xlim(0, 0.35)
ax2.set_xlim(0, 0.35)

ax1.text(-0.10, 9, 'a)')
ax2.text(-0.05, 9, 'b)')

ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_tclw_frac_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot SO Cloud Liquid Water Content Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(ecmwf_tclw_alt_so[18:,1]*1000,ecmwf_tclw_alt_so[18:,0], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tclw_alt_so[:18,1]*1000,gfdl_hiram_tclw_alt_so[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_alt_so[:19,1]*1000,mri_cgcm_tclw_alt_so[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_alt_so[:13,1]*1000,cam5_tclw_alt_so[:13,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')
ax1.plot(miroc5_tclw_alt_so[:,1],miroc5_tclw_alt_so[:,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tclw_alt_so[:,1],ipsl5_tclw_alt_so[:,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tclw_alt_so[:,1],giss5_tclw_alt_so[:,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(ecmwf_tclw_alt_so[18:,1]*1000,ecmwf_tclw_alt_so[18:,0], '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tclw_alt_so[:19,1]*1000,gfdl4_tclw_alt_so[:19,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_alt_so[:26,1]*1000,mri_tclw_alt_so[:26,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_alt_so[17:,1]*1000,cam6_tclw_alt_so[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tclw_alt_so[:,1],miroc6_tclw_alt_so[:,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tclw_alt_so[:,1],ipsl6_tclw_alt_so[:,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tclw_alt_so[:,1],giss6_tclw_alt_so[:,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Liquid Water Content (g/kg)')
ax2.set_xlabel('Liquid Water Content (g/kg)')

#ax1.set_title('2001 to 2005 Global Cloud Liquid Water Fraction vs Altitude')

ax1.set_xlim(0, 0.055)
ax2.set_xlim(0, 0.055)

ax1.text(-0.006, 9, 'a)')
ax2.text(-0.006, 9, 'b)')
ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_tclw_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""
############################################################################### Contour plots

#---Combined Grid tclw_frac---#

fig, ax = plt.subplots(nrows=6, ncols=3, figsize=(10, 10))

ax[0, 2].contourf(ecmwf_lat, ecmwf_alt[19:], ecmwf_tclw_frac_alt_lat[19:], vmin=0, vmax=0.4)
ax[0, 2].set_xlabel('Latitude')
ax[0, 0].set_ylabel('Altitude (km)')
ax[0, 2].set_title('c) ECMWF-ERA5')
ecmwf_temp = ax[0, 2].contour(ecmwf_lat, ecmwf_alt[19:], (ecmwf_temp_alt_lat[19:] - 273.15), colors='grey')
ecmwf_temp.collections[5].set_linewidth(2)
ecmwf_temp.collections[5].set_color('white')
ax[0, 2].clabel(ecmwf_temp, inline=1, fontsize=10)

ax[0, 0].contourf(gfdl_hiram_lat, gfdl_hiram_alt[0:17], gfdl_hiram_tclw_frac_alt_lat[0:17], vmin=0, vmax=0.4)
ax[0, 0].set_title('a) CMIP5-GFDL-HIRAM')
gfdl_hiram_temp = ax[0, 0].contour(gfdl_hiram_lat, gfdl_hiram_alt_temp[1:7], (gfdl_hiram_temp_alt_lat[1:7] - 273.15), colors='grey')
gfdl_hiram_temp.collections[7].set_linewidth(2)
gfdl_hiram_temp.collections[7].set_color('white')
ax[0, 0].clabel(gfdl_hiram_temp, inline=1, fontsize=10)

ax[0, 1].contourf(gfdl4_lat, gfdl4_alt[0:18], gfdl4_tclw_frac_alt_lat[0:18], vmin=0, vmax=0.4)
ax[0, 1].set_title('b) CMIP6-GFDL-AM4')
gfdl4_temp = ax[0, 1].contour(gfdl4_lat, gfdl4_alt_temp[1:7], (gfdl4_temp_alt_lat[1:7] - 273.15), colors='grey')
gfdl4_temp.collections[6].set_linewidth(2)
gfdl4_temp.collections[6].set_color('white')
ax[0, 1].clabel(gfdl4_temp, inline=1, fontsize=10)

ax[1, 0].contourf(mri_cgcm_lat, mri_cgcm_alt[0:18], mri_cgcm_tclw_frac_alt_lat[0:18], vmin=0, vmax=0.4)
ax[1, 0].set_title('d) CMIP5-MRI-CGCM3')
ax[1, 0].set_ylabel('Altitude (km)')
mri_cgcm_temp = ax[1, 0].contour(mri_cgcm_lat, mri_cgcm_alt_temp[1:7], (mri_cgcm_temp_alt_lat[1:7] - 273.15), colors='grey')
mri_cgcm_temp.collections[5].set_linewidth(2)
mri_cgcm_temp.collections[5].set_color('white')
ax[1, 0].clabel(mri_cgcm_temp, inline=1, fontsize=10)

cont=ax[1, 1].contourf(mri_lat, mri_alt[0:25], mri_tclw_frac_alt_lat[0:25], vmin=0, vmax=0.4)
ax[1, 1].set_title('e) CMIP6-MRI_ESM2')
mri_temp = ax[1, 1].contour(mri_lat, mri_alt_temp[1:7], (mri_temp_alt_lat[1:7] - 273.15), colors='grey')
mri_temp.collections[5].set_linewidth(2)
mri_temp.collections[5].set_color('white')
ax[1, 1].clabel(mri_temp, inline=1, fontsize=10)

ax[2, 0].contourf(cam5_lat, cam5_alt[:13], cam5_tclw_frac_alt_lat[:13], vmin=0, vmax=0.4)
ax[2, 0].set_title('f) CMIP5-CESM1-CAM5')
ax[2, 0].set_ylabel('Altitude (km)')
cam5_temp = ax[2, 0].contour(cam5_lat, cam5_alt_temp[1:7], (cam5_temp_alt_lat[1:7] - 273.15), colors='grey')
cam5_temp.collections[6].set_linewidth(2)
cam5_temp.collections[6].set_color('white')
ax[2, 0].clabel(cam5_temp, inline=1, fontsize=10)

ax[2, 1].contourf(cam6_lat, cam6_alt[17:32], cam6_tclw_frac_alt_lat[17:32], vmin=0, vmax=0.4)
ax[2, 1].set_title('g) CMIP6-CESM2.1-CAM6')
cam6_temp = ax[2, 1].contour(cam6_lat, cam6_alt[17:32], (cam6_temp_alt_lat[17:32] - 273.15), colors='grey')
cam6_temp.collections[5].set_linewidth(2)
cam6_temp.collections[5].set_color('white')
ax[2, 1].clabel(cam6_temp, cam6_temp.levels[:5:], inline=1, fontsize=10)

ax[3, 0].contourf(miroc5_lat, miroc5_alt[:13], miroc5_tclw_frac_alt_lat[:13], vmin=0, vmax=0.4)
ax[3, 0].set_title('h) CMIP5-MIROC5')
ax[3, 0].set_ylabel('Altitude (km)')
miroc5_temp = ax[3, 0].contour(miroc5_lat, miroc5_alt_temp[1:7], (miroc5_temp_alt_lat[1:7] - 273.15), colors='grey')
miroc5_temp.collections[6].set_linewidth(2)
miroc5_temp.collections[6].set_color('white')
ax[3, 0].clabel(miroc5_temp, inline=1, fontsize=10)

ax[3, 1].contourf(miroc6_lat, miroc6_alt[17:32], miroc6_tclw_frac_alt_lat[17:32], vmin=0, vmax=0.4)
ax[3, 1].set_title('i) CMIP6-MIROC6')
miroc6_temp = ax[3, 1].contour(miroc6_lat, miroc6_alt[17:32], (miroc6_temp_alt_lat[17:32] - 273.15), colors='grey')
miroc6_temp.collections[5].set_linewidth(2)
miroc6_temp.collections[5].set_color('white')
ax[3, 1].clabel(miroc6_temp, miroc6_temp.levels[:5:], inline=1, fontsize=10)

ax[4, 0].contourf(giss5_lat, giss5_alt[:13], giss5_tclw_frac_alt_lat[:13], vmin=0, vmax=0.4)
ax[4, 0].set_title('j) CMIP5-GISS-E2R')
ax[4, 0].set_ylabel('Altitude (km)')
giss5_temp = ax[4, 0].contour(giss5_lat, giss5_alt_temp[1:7], (giss5_temp_alt_lat[1:7] - 273.15), colors='grey')
giss5_temp.collections[6].set_linewidth(2)
giss5_temp.collections[6].set_color('white')
ax[4, 0].clabel(giss5_temp, inline=1, fontsize=10)

ax[4, 1].contourf(giss6_lat, giss6_alt[17:32], giss6_tclw_frac_alt_lat[17:32], vmin=0, vmax=0.4)
ax[4, 1].set_title('k) CMIP6-GISS-E21G')
giss6_temp = ax[4, 1].contour(giss6_lat, giss6_alt[17:32], (giss6_temp_alt_lat[17:32] - 273.15), colors='grey')
giss6_temp.collections[5].set_linewidth(2)
giss6_temp.collections[5].set_color('white')
ax[4, 1].clabel(giss6_temp, giss6_temp.levels[:5:], inline=1, fontsize=10)

ax[5, 0].contourf(ipsl5_lat, ipsl5_alt[:13], ipsl5_tclw_frac_alt_lat[:13], vmin=0, vmax=0.4)
ax[5, 0].set_title('l) CMIP5-IPSL-CM5A-LR')
ax[5, 0].set_ylabel('Altitude (km)')
ax[5, 0].set_xlabel('Latitude')
ipsl5_temp = ax[5, 0].contour(ipsl5_lat, ipsl5_alt_temp[1:7], (ipsl5_temp_alt_lat[1:7] - 273.15), colors='grey')
ipsl5_temp.collections[6].set_linewidth(2)
ipsl5_temp.collections[6].set_color('white')
ax[5, 0].clabel(ipsl5_temp, inline=1, fontsize=10)

ax[5, 1].contourf(ipsl6_lat, ipsl6_alt[17:32], ipsl6_tclw_frac_alt_lat[17:32], vmin=0, vmax=0.4)
ax[5, 1].set_title('m) CMIP6-IPSL-CM6A-LR')
ax[5, 1].set_xlabel('Latitude')
ipsl6_temp = ax[5, 1].contour(ipsl6_lat, ipsl6_alt[17:32], (ipsl6_temp_alt_lat[17:32] - 273.15), colors='grey')
ipsl6_temp.collections[5].set_linewidth(2)
ipsl6_temp.collections[5].set_color('white')
ax[5, 1].clabel(ipsl6_temp, ipsl6_temp.levels[:5:], inline=1, fontsize=10)



ax[1, 2].remove()  # don't display empty ax
ax[2, 2].remove()  # don't display empty ax
ax[3, 2].remove()  # don't display empty ax
ax[4, 2].remove()  # don't display empty ax
ax[5, 2].remove()  # don't display empty ax




cbaxes = fig.add_axes([0.8, 0.2, 0.03, 0.3]) 
cbar = fig.colorbar(cont, cax=cbaxes)
cbar.set_clim(0, 0.42)
cbar.set_label('Cloud Liquid Water Fraction')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)
plt.savefig("2001_2005_contour_tclw.svg", format="svg", bbox_inches='tight')
plt.show()



#---Combined Grid tciw_frac---#
"""

fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(10, 10))

cont=ax[0, 2].contourf(ecmwf_lat, ecmwf_alt[9:], ecmwf_tciw_frac_alt_lat[9:], vmin=0, vmax=0.5)
ax[0, 2].set_xlabel('Latitude')
ax[0, 0].set_ylabel('Altitude (km)')
ax[0, 2].set_title('c) ECMWF-ERA5')
ecmwf_temp = ax[0, 2].contour(ecmwf_lat, ecmwf_alt[9:], (ecmwf_temp_alt_lat[9:] - 273.15), colors='grey')
ecmwf_temp.collections[6].set_linewidth(2)
ecmwf_temp.collections[6].set_color('white')
ax[0, 2].clabel(ecmwf_temp, inline=1, fontsize=10)

ax[0, 0].contourf(gfdl_hiram_lat, gfdl_hiram_alt[:25], gfdl_hiram_tciw_frac_alt_lat[:25], vmin=0, vmax=0.5)
ax[0, 0].set_title('a) CMIP5-GFDL-HIRAM')
gfdl_hiram_temp = ax[0, 0].contour(gfdl_hiram_lat, gfdl_hiram_alt_temp[1:13], (gfdl_hiram_temp_alt_lat[1:13] - 273.15), colors='grey')
gfdl_hiram_temp.collections[6].set_linewidth(2)
gfdl_hiram_temp.collections[6].set_color('white')
ax[0, 0].clabel(gfdl_hiram_temp, inline=1, fontsize=10)

ax[0, 1].contourf(gfdl4_lat, gfdl4_alt[:26], gfdl4_tciw_frac_alt_lat[:26], vmin=0, vmax=0.5)
ax[0, 1].set_title('b) CMIP6-GFDL-AM4')
gfdl4_temp = ax[0, 1].contour(gfdl4_lat, gfdl4_alt_temp[1:13], (gfdl4_temp_alt_lat[1:13] - 273.15), colors='grey')
gfdl4_temp.collections[6].set_linewidth(2)
gfdl4_temp.collections[6].set_color('white')
ax[0, 1].clabel(gfdl4_temp, inline=1, fontsize=10)

ax[1, 0].contourf(mri_cgcm_lat, mri_cgcm_alt[:30], mri_cgcm_tciw_frac_alt_lat[:30], vmin=0, vmax=0.5)
ax[1, 0].set_title('d) CMIP5-MRI-CGCM3')
ax[1, 0].set_ylabel('Altitude (km)')
mri_cgcm_temp = ax[1, 0].contour(mri_cgcm_lat, mri_cgcm_alt_temp[1:13], (mri_cgcm_temp_alt_lat[1:13] - 273.15), colors='grey')
mri_cgcm_temp.collections[6].set_linewidth(2)
mri_cgcm_temp.collections[6].set_color('white')
ax[1, 0].clabel(mri_cgcm_temp, inline=1, fontsize=10)

ax[1, 1].contourf(mri_lat, mri_alt[:44], mri_tciw_frac_alt_lat[:44], vmin=0, vmax=0.5)
ax[1, 1].set_title('e) CMIP6-MRI_ESM2')
mri_temp = ax[1, 1].contour(mri_lat, mri_alt_temp[1:12], (mri_temp_alt_lat[1:12] - 273.15), colors='grey')
mri_temp.collections[6].set_linewidth(2)
mri_temp.collections[6].set_color('white')
ax[1, 1].clabel(mri_temp, inline=1, fontsize=10)

ax[2, 0].contourf(cam5_lat, cam5_alt[:20], cam5_tciw_frac_alt_lat[:20], vmin=0, vmax=0.5)
ax[2, 0].set_title('f) CMIP5-CESM1-CAM5')
ax[2, 0].set_ylabel('Altitude (km)')
ax[2, 0].set_xlabel('Latitude')
cam5_temp = ax[2, 0].contour(cam5_lat, cam5_alt_temp[1:18], (cam5_temp_alt_lat[1:18] - 273.15), colors='grey')
cam5_temp.collections[6].set_linewidth(2)
cam5_temp.collections[6].set_color('white')
ax[2, 0].clabel(cam5_temp, inline=1, fontsize=10)

ax[2, 1].contourf(cam6_lat, cam6_alt[7:32], cam6_tciw_frac_alt_lat[7:32], vmin=0, vmax=0.5)
ax[2, 1].set_title('g) CMIP6-CESM2.1-CAM6')
ax[2, 1].set_xlabel('Latitude')
cam6_temp = ax[2, 1].contour(cam6_lat, cam6_alt[7:32], (cam6_temp_alt_lat[7:32] - 273.15), colors='grey')
cam6_temp.collections[6].set_linewidth(2)
cam6_temp.collections[6].set_color('white')
ax[2, 1].clabel(cam6_temp, inline=1, fontsize=10)


ax[1, 2].remove()  # don't display empty ax
ax[2, 2].remove()  # don't display empty ax

cbaxes = fig.add_axes([0.8, 0.2, 0.03, 0.3]) 
cbar = fig.colorbar(cont, cax=cbaxes)
cbar.set_clim(0, 0.5)
cbar.set_label('Cloud Ice Water Fraction')
plt.tight_layout(pad=0.4, w_pad=0.6, h_pad=3)
plt.savefig("2001_2005_contour_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""





#---Plot ECMWF Cloud Liquid Water Fraction---#

"""
plt.subplots()

plt.contourf(ecmwf_lat, ecmwf_alt[19:], ecmwf_tclw_frac_alt_lat[19:])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Liquid Cloud Fraction')

plt.title('2001 to 2005 ECMWF-ERA5 Cloud Liquid Water Fraction')

plt.savefig("2001_2005_contour_ecmwf_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot ECMWF Cloud Ice Water Fraction---#

"""
plt.subplots()

plt.contourf(ecmwf_lat, ecmwf_alt[9:], ecmwf_tciw_frac_alt_lat[9:])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Ice Cloud Fraction')

plt.title('2001 to 2005 ECMWF-ERA5 Cloud Ice Water Fraction')

plt.savefig("2001_2005_contour_ecmwf_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot GFDL-HIRAM-AMIP Cloud Liquid Water Fraction---#

"""
plt.subplots()

plt.contourf(gfdl_hiram_lat, gfdl_hiram_alt[0:17], gfdl_hiram_tclw_frac_alt_lat[0:17])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Liquid Cloud Fraction')

plt.title('2001 to 2005 GFDL-HIRAM-AMIP Cloud Liquid Water Fraction')

plt.savefig("2001_2005_contour_gfdl_hiram_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot GFDL-HIRAM-AMIP Cloud Ice Water Fraction---#

"""
plt.subplots()

plt.contourf(gfdl_hiram_lat, gfdl_hiram_alt[:26], gfdl_hiram_tciw_frac_alt_lat[:26])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Ice Cloud Fraction')

plt.title('2001 to 2005 GFDL-HIRAM-AMIP Cloud Ice Water Fraction')

plt.savefig("2001_2005_contour_gfdl_hiram.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot GFDL-AM4-AMIP Cloud Liquid Water Fraction---#

"""
plt.subplots()

plt.contourf(gfdl4_lat, gfdl4_alt[0:18], gfdl4_tclw_frac_alt_lat[0:18])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Liquid Cloud Fraction')

plt.title('2001 to 2005 GFDL-AM4-AMIP Cloud Liquid Water Fraction')

plt.savefig("2001_2005_contour_gfdl4_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot GFDL-AM4-AMIP Cloud Ice Water Fraction---#

"""
plt.subplots()

plt.contourf(gfdl4_lat, gfdl4_alt[:26], gfdl4_tciw_frac_alt_lat[:26])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Ice Cloud Fraction')

plt.title('2001 to 2005 GFDL-AM4-AMIP Cloud Ice Water Fraction')

plt.savefig("2001_2005_contour_gfdl4_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot MRI-CGCM3-AMIP Cloud Liquid Water Fraction---#

"""
plt.subplots()

plt.contourf(mri_cgcm_lat, mri_cgcm_alt[0:18], mri_cgcm_tclw_frac_alt_lat[0:18])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Liquid Cloud Fraction')

plt.title('2001 to 2005 MRI-CGCM3-AMIP Cloud Liquid Water Fraction')

plt.savefig("2001_2005_contour_mri_cgcm_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot MRI-CGCM3-AMIP Cloud Ice Water Fraction---#

"""
plt.subplots()
plt.contourf(mri_cgcm_lat, mri_cgcm_alt[:30], mri_cgcm_tciw_frac_alt_lat[:30])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Ice Cloud Fraction')

plt.title('2001 to 2005 MRI-CGCM3-AMIP Cloud Ice Water Fraction')

plt.savefig("2001_2005_contour_mri_cgcm_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot MRI_ESM2-AMIP Cloud Liquid Water Fraction---#

"""
plt.subplots()

plt.contourf(mri_lat, mri_alt[0:25], mri_tclw_frac_alt_lat[0:25])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Liquid Cloud Fraction')

plt.title('2001 to 2005 MRI_ESM2-AMIP Cloud Liquid Water Fraction')

plt.savefig("2001_2005_contour_mri_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot MRI_ESM2-AMIP Cloud Ice Water Fraction---#

"""
plt.subplots()

plt.contourf(mri_lat, mri_alt[:44], mri_tciw_frac_alt_lat[:44])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Ice Cloud Fraction')

plt.title('2001 to 2005 MRI_ESM2-AMIP Cloud Ice Water Fraction')

plt.savefig("2001_2005_contour_mri_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CESM2-CAM6-AMIP Cloud Liquid Water Fraction---#

"""
plt.subplots()
plt.contourf(cam6_lat, cam6_alt[18:32], cam6_tclw_frac_alt_lat[18:32])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Liquid Cloud Fraction')

plt.title('2001 to 2005 CESM2-CAM6-AMIP Cloud Liquid Water Fraction')

plt.savefig("2001_2005_contour_cam6_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CESM2-CAM6-AMIP Cloud Ice Water Fraction---#

"""
plt.subplots()

plt.contourf(cam6_lat, cam6_alt[7:32], cam6_tciw_frac_alt_lat[7:32])

plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Ice Cloud Fraction')

plt.title('2001 to 2005 CESM2-CAM6-AMIP Cloud Ice Water Fraction')

plt.savefig("2001_2005_contour_cam6_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CESM1-CAM5-AMIP Cloud Liquid Water Fraction---#

"""
plt.subplots()

plt.contourf(cam5_lat, cam5_alt[:13], cam5_tclw_frac_alt_lat[:13])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Liquid Cloud Fraction')

plt.title('2001 to 2005 CESM1-CAM5-AMIP Cloud Liquid Water Fraction')

plt.savefig("2001_2005_contour_cam5_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CESM1-CAM5-AMIP Cloud Ice Water Fraction---#

"""
plt.subplots()
plt.contourf(cam5_lat, cam5_alt[:20], cam5_tciw_frac_alt_lat[:20])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Ice Cloud Fraction')

plt.title('2001 to 2005 CESM1-CAM5-AMIP Cloud Ice Water Fraction')

plt.savefig("2001_2005_contour_cam5_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

