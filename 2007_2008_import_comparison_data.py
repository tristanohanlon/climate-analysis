# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This import the reduced datasets from 2007 to 2008 CCCM, ECMWF and gfdl4. 
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
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
#a = h5py.File('2007_2008_ecmwf.h5', 'r')

#GFDL-AM4-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('2007_2008_gfdl_am4.h5', 'r')

#GFDL-HIRAM-C360 Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL-HIRAM-C360/reduced_datasets')
h = h5py.File('2007_2008_gfdl_hiram.h5', 'r')

#CCCM Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('2007_2008_cccm.h5', 'r')

#MRI-ESM2-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
d = h5py.File('2007_2008_mri_esm2.h5', 'r')

#MRI-CGCM3-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/MRI-CGCM3-AMIP/reduced_datasets')
i = h5py.File('2007_2008_mri_cgcm3.h5', 'r')

#CESM2-CAM6-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets')
e = h5py.File('2007_2008_cesm2_cam6.h5', 'r')

#CAPLISO Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets')
f = h5py.File('2007_2008_calipso.h5', 'r')

#CERES Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets')
g = h5py.File('2007_2008_ceres.h5', 'r')
"""

# Home PC
#ECMWF Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
a = h5py.File('2007_2008_ecmwf_era5.h5', 'r')

#GFDL-AM4-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('2007_2008_gfdl_am4.h5', 'r')

#GFDL-HIRAM-C360 Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL-HIRAM-C360/reduced_datasets')
h = h5py.File('2007_2008_gfdl_hiram.h5', 'r')

#CCCM Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('2007_2008_cccm.h5', 'r')

#MRI-ESM2-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
d = h5py.File('2007_2008_mri_esm2.h5', 'r')

#MRI-CGCM3-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/MRI-CGCM3-AMIP/reduced_datasets')
i = h5py.File('2007_2008_mri_cgcm3.h5', 'r')

#CESM2-CAM6-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets')
e = h5py.File('2007_2008_cesm2_cam6.h5', 'r')

#CAPLISO Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets')
f = h5py.File('2007_2008_calipso.h5', 'r')

#CERES Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets')
g = h5py.File('2007_2008_ceres.h5', 'r')

"""
# Laptop
#ECMWF Data
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
#a = h5py.File('2007_2008_ecmwf.h5', 'r')

#GFDL-AM4-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('2007_2008_gfdl_am4.h5', 'r')

#GFDL-HIRAM-C360 Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/GFDL-HIRAM-C360/reduced_datasets')
h = h5py.File('2007_2008_gfdl_hiram.h5', 'r')

#CCCM Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('2007_2008_cccm.h5', 'r')

#MRI-ESM2-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
d = h5py.File('2007_2008_mri_esm2.h5', 'r')

#MRI-CGCM3-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/MRI-CGCM3-AMIP/reduced_datasets')
i = h5py.File('2007_2008_mri_cgcm3.h5', 'r')

#CESM2-CAM6-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets')
e = h5py.File('2007_2008_cesm2_cam6.h5', 'r')

#CAPLISO Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets')
f = h5py.File('2007_2008_calipso.h5', 'r')

#CERES Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets')
g = h5py.File('2007_2008_ceres.h5', 'r')
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

ecmwf_tclw_frac_alt_lat = (ecmwf_tclw_alt_lat / (ecmwf_tclw_alt_lat + ecmwf_tciw_alt_lat)) * (ecmwf_tcc_alt_lat)
ecmwf_tciw_frac_alt_lat = (ecmwf_tciw_alt_lat / (ecmwf_tclw_alt_lat + ecmwf_tciw_alt_lat)) * (ecmwf_tcc_alt_lat)


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




############################################################################### CCCM Data

#---CCCM Global Latitude Data---#

cccm_tcc_lat_g = c['tcc'][:] # 0-1
cccm_tcc_lat_g_enhanced  = c['tcc_enhanced'][:] # 0-1
cccm_tclw_lat_g = c['tclw'][:] #kgm^-2
cccm_tciw_lat_g = c['tciw'][:] #kgm^-2

#---CCCM Southern Ocean Latitude Data---#

cccm_tcc_lat_so = cccm_tcc_lat_g[cccm_tcc_lat_g[:,0]>=-70]
cccm_tcc_lat_so = cccm_tcc_lat_so[cccm_tcc_lat_so[:,0]<=-50] # 0-1

cccm_tclw_lat_so = cccm_tclw_lat_g[cccm_tclw_lat_g[:,0]>=-70]
cccm_tclw_lat_so = cccm_tclw_lat_so[cccm_tclw_lat_so[:,0]<=-50] #kgm^-2

cccm_tciw_lat_so = cccm_tciw_lat_g[cccm_tciw_lat_g[:,0]>=-70]
cccm_tciw_lat_so = cccm_tciw_lat_so[cccm_tciw_lat_so[:,0]<=-50] #kgm^-2

#---CCCM Phase Fractions---#

cccm_tclw_frac_lat_g = c['tclw_frac'][:]
cccm_tciw_frac_lat_g = c['tclw_frac'][:]

cccm_tclw_frac_lat_g_enhanced = c['tclw_frac_enhanced'][:]
cccm_tciw_frac_lat_g_enhanced = c['tclw_frac_enhanced'][:]

#---CCCM Global Profile---#

cccm_tcc_alt_g = c['cf'][:] # 0-1 4:101
cccm_tclw_alt_g = c['lw'][:113] #kg/kg 4:63
cccm_tciw_alt_g = c['iw'][:113] #kg/kg 4:93
cccm_temp_alt_g = c['temp'][:] #K
cccm_plevel_alt_g = c['pressure'][:] #hPa

cccm_tcc_temp_g = c['cf_t'][:] # 0-1
cccm_tclw_temp_g = c['lw_t'][:] #kg/kg
cccm_tciw_temp_g = c['iw_t'][:] #kg/kg

cccm_tclw_frac_alt_g = c['lw_frac'][:]
cccm_tciw_frac_alt_g = c['iw_frac'][:]

#---CCCM Temperature Phase Fractions---#

cccm_tclw_frac_temp_g = c['lw_frac_temp'][:]
cccm_tciw_frac_temp_g = c['iw_frac_temp'][:]

cccm_tclw_frac_temp_so = c['lw_frac_temp_so'][:]
cccm_tciw_frac_temp_so = c['iw_frac_temp_so'][:]


#---CCCM Southern Ocean Profile---#

cccm_tcc_alt_so = c['cf_so'][:] # 0-1
cccm_tclw_alt_so = c['lw_so'][:] #kg/kg
cccm_tciw_alt_so = c['iw_so'][:] #kg/kg
cccm_temp_alt_so = c['temp_so'][:] #K
cccm_plevel_alt_so = c['pressure_so'][:] #hPa

cccm_tcc_temp_so = c['cf_t_so'][:] # 0-1
cccm_tclw_temp_so = c['lw_t_so'][:] #kg/kg
cccm_tciw_temp_so = c['iw_t_so'][:] #kg/kg

cccm_tclw_frac_alt_so = c['lw_frac_so'][:]
cccm_tciw_frac_alt_so = c['iw_frac_so'][:]





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

gfdl4_tclw_frac_alt_lat = (gfdl4_tclw_alt_lat / (gfdl4_tclw_alt_lat + gfdl4_tciw_alt_lat)) * (gfdl4_tcc_alt_lat)
gfdl4_tciw_frac_alt_lat = (gfdl4_tciw_alt_lat / (gfdl4_tclw_alt_lat + gfdl4_tciw_alt_lat)) * (gfdl4_tcc_alt_lat)


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

gfdl_hiram_tclw_frac_alt_lat = (gfdl_hiram_tclw_alt_lat / (gfdl_hiram_tclw_alt_lat + gfdl_hiram_tciw_alt_lat)) * (gfdl_hiram_tcc_alt_lat)
gfdl_hiram_tciw_frac_alt_lat = (gfdl_hiram_tciw_alt_lat / (gfdl_hiram_tclw_alt_lat + gfdl_hiram_tciw_alt_lat)) * (gfdl_hiram_tcc_alt_lat)


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


mri_tclw_frac_alt_lat = (mri_tclw_alt_lat / (mri_tclw_alt_lat + mri_tciw_alt_lat)) * (mri_tcc_alt_lat)
mri_tciw_frac_alt_lat = (mri_tciw_alt_lat / (mri_tclw_alt_lat + mri_tciw_alt_lat)) * (mri_tcc_alt_lat)


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


mri_cgcm_tclw_frac_alt_lat = (mri_cgcm_tclw_alt_lat / (mri_cgcm_tclw_alt_lat + mri_cgcm_tciw_alt_lat)) * (mri_cgcm_tcc_alt_lat)
mri_cgcm_tciw_frac_alt_lat = (mri_cgcm_tciw_alt_lat / (mri_cgcm_tclw_alt_lat + mri_cgcm_tciw_alt_lat)) * (mri_cgcm_tcc_alt_lat)


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

cam6_tclw_frac_alt_lat = (cam6_tclw_alt_lat / (cam6_tclw_alt_lat + cam6_tciw_alt_lat)) * (cam6_tcc_alt_lat)
cam6_tciw_frac_alt_lat = (cam6_tciw_alt_lat / (cam6_tclw_alt_lat + cam6_tciw_alt_lat)) * (cam6_tcc_alt_lat)


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


############################################################################### CALIPSO Data

#---CALIPSO Global Latitude Data---#

calipso_tcc_lat_g = f['tcc'][:] # 0-1
calipso_tclw_frac_lat_g = f['tclw_frac'][:] # 0-1
calipso_tciw_frac_lat_g = f['tciw_frac'][:] # 0-1

#---CALIPSO Southern Ocean Latitude Data---#

calipso_tcc_lat_so = calipso_tcc_lat_g[calipso_tcc_lat_g[:,0]>=-70]
calipso_tcc_lat_so = calipso_tcc_lat_so[calipso_tcc_lat_so[:,0]<=-50] # 0-1

calipso_tclw_frac_lat_so = calipso_tclw_frac_lat_g[calipso_tclw_frac_lat_g[:,0]>=-70]
calipso_tclw_frac_lat_so = calipso_tclw_frac_lat_so[calipso_tclw_frac_lat_so[:,0]<=-50] # 0-1

calipso_tciw_frac_lat_so = calipso_tciw_frac_lat_g[calipso_tciw_frac_lat_g[:,0]>=-70]
calipso_tciw_frac_lat_so = calipso_tciw_frac_lat_so[calipso_tciw_frac_lat_so[:,0]<=-50] # 0-1


#---CALIPSO Global Profile---#

calipso_tcc_alt_g = f['cf'][:] # 0-1
calipso_tclw_frac_alt_g = f['lw_frac'][:] # 0-1
calipso_tciw_frac_alt_g = f['iw_frac'][:] # 0-1


#---CALIPSO lat-alt contour data---#

calipso_tcc_alt_lat = f['cf_alt_lat'][:]
calipso_tclw_frac_alt_lat = f['liq_frac_alt_lat'][:]
calipso_tciw_frac_alt_lat = f['ice_frac_alt_lat'][:]
calipso_lat = f['lat'][:]
calipso_lat = np.hstack(calipso_lat)
calipso_alt = f['alt'][:]
calipso_alt = np.hstack(calipso_alt)
calipso_temp = f['alt_t'][:]
calipso_temp = np.hstack(calipso_temp)


#---CALIPSO lat-temp contour data---#

calipso_cf_t_lat = f['cf_t_lat'][:]
calipso_lw_t_lat = f['lw_t_lat'][:]
calipso_iw_t_lat = f['iw_t_lat'][:]

#---CALIPSO temp profile data---#

calipso_tcc_temp_g = f['cf_t'][:] # 0-1
calipso_tclw_frac_temp_g = f['lw_t_frac'][:] # 0-1
calipso_tciw_frac_temp_g = f['iw_t_frac'][:] # 0-1


#---CALIPSO Southern Ocean Profile---#


calipso_tcc_alt_so = f['cf_so'][:] # 0-1
calipso_tclw_frac_alt_so = f['lw_frac_so'][:] # 0-1
calipso_tciw_frac_alt_so = f['iw_frac_so'][:] # 0-1

calipso_tcc_temp_so = f['cf_t_so'][:] # 0-1
calipso_tclw_frac_temp_so = f['lw_t_frac_so'][:] # 0-1
calipso_tciw_frac_temp_so = f['iw_t_frac_so'][:] # 0-1

############################################################################### CERES Data

#---CERES Global Latitude Data---#

ceres_tcc_lat_g = g['tcc'][:] # 0-1
ceres_tclw_lat_g = g['tclw'][:] #kgm^-2
ceres_tciw_lat_g = g['tciw'][:] #kgm^-2
ceres_tclw_frac_lat_g = g['tclw_frac'][:] # 0-1
ceres_tciw_frac_lat_g = g['tciw_frac'][:] # 0-1

#---CERES Southern Ocean Latitude Data---#

ceres_tcc_lat_so = ceres_tcc_lat_g[ceres_tcc_lat_g[:,0]>=-70]
ceres_tcc_lat_so = ceres_tcc_lat_so[ceres_tcc_lat_so[:,0]<=-50] # 0-1

ceres_tclw_lat_so = ceres_tclw_lat_g[ceres_tclw_lat_g[:,0]>=-70]
ceres_tclw_lat_so = ceres_tclw_lat_so[ceres_tclw_lat_so[:,0]<=-50] #kgm^-2

ceres_tciw_lat_so = ceres_tciw_lat_g[ceres_tciw_lat_g[:,0]>=-70]
ceres_tciw_lat_so = ceres_tciw_lat_so[ceres_tciw_lat_so[:,0]<=-50] #kgm^-2


############################################################################### End importing Data

end = time.time()
print('Importing data took:', end - start, 's')


"""
os.chdir('c:/Users/toha006/University/University/MSc/Models/Images/Meeting 1.7/2007 - 2008')
"""
os.chdir('e:/University/University/MSc/Models/Images')

############################################################################### Temperature Profiles

#---Plot Satellite Liquid Cloud Fraction with Temperature---#

"""

fig, ax = plt.subplots()


ax.plot(cccm_tclw_frac_temp_g[5:44,0],cccm_tclw_frac_temp_g[5:44,1], '-r', label='CCCM - Global')
ax.plot(cccm_tclw_frac_temp_so[2:40,0],cccm_tclw_frac_temp_so[2:40,1], '--r', label='CCCM - Southern Ocean')

ax.plot(calipso_tclw_frac_temp_g[20:34,0],calipso_tclw_frac_temp_g[20:34,1], '-b', label='CALIPSO - Global')
ax.plot(calipso_tclw_frac_temp_so[20:34,0],calipso_tclw_frac_temp_so[20:34,1], '--b', label='CALIPSO - Southern Ocean')
plt.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Liquid Cloud Fraction')

plt.title('2007 to 2008 Satellite Liquid Cloud Fraction vs Temperature')
plt.grid(True)
plt.savefig("2007_2008_satellite_liquid_T.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Liquid Fraction with Temperature---#



"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tclw_frac_temp_g[5:44,0],(cccm_tclw_frac_temp_g[5:44,1]/cccm_tcc_temp_g[5:44,1]), ':r', label='CCCM')
ax1.plot(calipso_tclw_frac_temp_g[17:34,0],(calipso_tclw_frac_temp_g[17:34,1]/calipso_tcc_temp_g[17:34,1]), ':b', label='CALIPSO')
ax1.plot(ecmwf_tclw_frac_temp_g[18:,0],(ecmwf_tclw_frac_temp_g[18:,1]/ecmwf_tcc_temp_g[18:,1]), '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tclw_frac_temp_g[:19,0],(gfdl_hiram_tclw_frac_temp_g[:19,1]/gfdl_hiram_tcc_temp_g[:19,1]), '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_g[:22,0],(mri_cgcm_tclw_frac_temp_g[:22,1]/mri_cgcm_tcc_temp_g[:22,1]), '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tclw_frac_temp_g[5:44,0],(cccm_tclw_frac_temp_g[5:44,1]/cccm_tcc_temp_g[5:44,1]), ':r', label='CCCM')
ax2.plot(calipso_tclw_frac_temp_g[17:34,0],(calipso_tclw_frac_temp_g[17:34,1]/calipso_tcc_temp_g[17:34,1]), ':b', label='CALIPSO')
ax2.plot(ecmwf_tclw_frac_temp_g[18:,0],(ecmwf_tclw_frac_temp_g[18:,1]/ecmwf_tcc_temp_g[18:,1]), '-k', label='ECMWF-ERA5')
ax2.plot(gfdl4_tclw_frac_temp_g[:20,0],(gfdl4_tclw_frac_temp_g[:20,1]/gfdl4_tcc_temp_g[:20,1]), '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_g[:28,0],(mri_tclw_frac_temp_g[:28,1]/mri_tcc_temp_g[:28,1]), '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_g[16:,0],(cam6_tclw_frac_temp_g[16:,1]/cam6_tcc_temp_g[16:,1]), '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Fraction')

#ax1.set_title('2007 to 2008 Global Liquid Fraction vs Temperature')

ax1.grid(True)
ax2.grid(True)
ax1.text(205, 1, 'a)')
ax2.text(205, 1, 'b)')
plt.savefig("2007_2008_global_liquid_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot Global Liquid Cloud Fraction with Temperature---#


"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tclw_frac_temp_g[5:44,0],cccm_tclw_frac_temp_g[5:44,1], ':r', label='CCCM')
ax1.plot(calipso_tclw_frac_temp_g[20:34,0],calipso_tclw_frac_temp_g[20:34,1], ':b', label='CALIPSO')
ax1.plot(ecmwf_tclw_frac_temp_so[18:,0],ecmwf_tclw_frac_temp_so[18:,1], '-k', label='ECMWF-ERA5')

ax1.plot(gfdl_hiram_tclw_frac_temp_g[:17,0],gfdl_hiram_tclw_frac_temp_g[:17,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_g[:18,0],mri_cgcm_tclw_frac_temp_g[:18,1], '-m', label='CMIP5-MRI_CGCM3-AMIP')


ax2.plot(cccm_tclw_frac_temp_g[5:44,0],cccm_tclw_frac_temp_g[5:44,1], ':r', label='CCCM')
ax2.plot(calipso_tclw_frac_temp_g[20:34,0],calipso_tclw_frac_temp_g[20:34,1], ':b', label='CALIPSO')
ax1.plot(ecmwf_tclw_frac_temp_so[18:,0],ecmwf_tclw_frac_temp_so[18:,1], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tclw_frac_temp_g[:18,0],gfdl4_tclw_frac_temp_g[:18,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_g[:26,0],mri_tclw_frac_temp_g[:26,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_g[19:,0],cam6_tclw_frac_temp_g[19:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

ax1.legend(loc='upper center', bbox_to_anchor=(1.5, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.5, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Cloud Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Cloud Fraction')

#ax1.set_title('2007 to 2008 Global Liquid Cloud Fraction vs Temperature')

plt.savefig("2007_2008_global_liquid_frac_T.svg", format="svg", bbox_inches='tight')
ax1.text(208, 0.3, 'a)')
ax2.text(208, 0.3, 'b)')

ax1.grid(True)
ax2.grid(True)
"""



#---Plot Southern Ocean Liquid Fraction with Temperature---#


"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tclw_frac_temp_so[2:40,0],(cccm_tclw_frac_temp_so[2:40,1]/cccm_tcc_temp_so[2:40,1]), ':r', label='CCCM')
ax1.plot(calipso_tclw_frac_temp_so[17:34,0],(calipso_tclw_frac_temp_so[17:34,1]/calipso_tcc_temp_so[17:34,1]), ':b', label='CALIPSO')
ax1.plot(ecmwf_tclw_frac_temp_so[18:,0],(ecmwf_tclw_frac_temp_so[18:,1]/ecmwf_tcc_temp_so[18:,1]), '-k', label='ECMWF-ERA5')

ax1.plot(gfdl_hiram_tclw_frac_temp_so[:19,0],(gfdl_hiram_tclw_frac_temp_so[:19,1]/gfdl_hiram_tcc_temp_so[:19,1]), '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_so[:22,0],(mri_cgcm_tclw_frac_temp_so[:22,1]/mri_cgcm_tcc_temp_so[:22,1]), '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tclw_frac_temp_so[2:40,0],(cccm_tclw_frac_temp_so[2:40,1]/cccm_tcc_temp_so[2:40,1]), ':r', label='CCCM')
ax2.plot(calipso_tclw_frac_temp_so[17:34,0],(calipso_tclw_frac_temp_so[17:34,1]/calipso_tcc_temp_so[17:34,1]), ':b', label='CALIPSO')
ax2.plot(ecmwf_tclw_frac_temp_so[18:,0],(ecmwf_tclw_frac_temp_so[18:,1]/ecmwf_tcc_temp_so[18:,1]), '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tclw_frac_temp_so[:20,0],(gfdl4_tclw_frac_temp_so[:20,1]/gfdl4_tcc_temp_so[:20,1]), '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_so[:28,0],(mri_tclw_frac_temp_so[:28,1]/mri_tcc_temp_so[:28,1]), '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_so[16:,0],(cam6_tclw_frac_temp_so[16:,1]/cam6_tcc_temp_so[16:,1]), '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Cloud Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Cloud Fraction')

#ax1.set_title('2007 to 2008 Southern Ocean Liquid Cloud Fraction vs Temperature')

ax1.text(200, 1, 'a)')
ax2.text(200, 1, 'b)')
ax1.grid(True)
ax2.grid(True)
plt.savefig("2007_2008_so_liquid_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Southern Ocean Liquid Cloud Fraction with Temperature---#


"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tclw_frac_temp_so[2:40,0],cccm_tclw_frac_temp_so[2:40,1], ':r', label='CCCM')
ax1.plot(calipso_tclw_frac_temp_so[20:34,0],calipso_tclw_frac_temp_so[20:34,1], ':b', label='CALIPSO')

ax1.plot(gfdl_hiram_tclw_frac_temp_so[:16,0],gfdl_hiram_tclw_frac_temp_so[:16,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_so[:18,0],mri_cgcm_tclw_frac_temp_so[:18,1], '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tclw_frac_temp_so[2:40,0],cccm_tclw_frac_temp_so[2:40,1], ':r', label='CCCM')
ax2.plot(calipso_tclw_frac_temp_so[20:34,0],calipso_tclw_frac_temp_so[20:34,1], ':b', label='CALIPSO')

ax2.plot(gfdl4_tclw_frac_temp_so[:18,0],gfdl4_tclw_frac_temp_so[:18,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_so[:26,0],mri_tclw_frac_temp_so[:26,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_so[19:,0],cam6_tclw_frac_temp_so[19:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

ax1.legend(loc='upper center', bbox_to_anchor=(1.5, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.5, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Cloud Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Cloud Fraction')

#ax1.set_title('2001 to 2005 Southern Ocean Liquid Cloud Fraction vs Temperature')

plt.savefig("2007_2008_so_liquid_frac_T.svg", format="svg", bbox_inches='tight')
ax1.text(210, 0.3, 'a)')
ax2.text(210, 0.3, 'b)')

ax1.grid(True)
ax2.grid(True)
"""
############################################################################### Global Latitude Plots

#---Plot Global Cloud Fraction with Latitude---#

"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tcc_lat_g_enhanced[:,0],cccm_tcc_lat_g_enhanced[:,1], ':r', label='CCCM')
ax1.plot(calipso_tcc_lat_g[:,0],calipso_tcc_lat_g[:,1], ':b', label='CAPLISO')
ax1.plot(ceres_tcc_lat_g[:,0],ceres_tcc_lat_g[:,1], '--r', label='CERES')
ax1.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-k', label='ECMWF-ERA5')

ax1.plot(gfdl_hiram_tcc_lat_g[:,0],gfdl_hiram_tcc_lat_g[:,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tcc_lat_g[:,0],mri_cgcm_tcc_lat_g[:,1], '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tcc_lat_g_enhanced[:,0],cccm_tcc_lat_g_enhanced[:,1], ':r', label='CCCM')
ax2.plot(calipso_tcc_lat_g[:,0],calipso_tcc_lat_g[:,1], ':b', label='CAPLISO')
ax2.plot(ceres_tcc_lat_g[:,0],ceres_tcc_lat_g[:,1], '--r', label='CERES')
ax2.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tcc_lat_g[:,0],gfdl4_tcc_lat_g[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_lat_g[:,0],mri_tcc_lat_g[:,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tcc_lat_g[:,0],cam6_tcc_lat_g[:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax1.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Cloud Fraction')
ax2.set_xlabel('Latitude')

ax1.set_title ('2007 - 2008 Global Cloud Fraction vs Latitude')

ax1.grid(True)
ax2.grid(True)
ax1.text(-140, 1, 'a)')
ax2.text(-140, 1, 'b)')

plt.grid(True)
plt.savefig("2007_2008_tcc_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot Global Comparison Liquid Cloud Fractions with Latitude---#

"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tclw_frac_lat_g_enhanced[:,0],cccm_tclw_frac_lat_g_enhanced[:,1], ':r', label='CCCM')
ax1.plot(calipso_tclw_frac_lat_g[:,0],calipso_tclw_frac_lat_g[:,1], ':b', label='CALIPSO')

ax1.plot(gfdl_hiram_tclw_frac_lat_g[:,0],gfdl_hiram_tclw_frac_lat_g[:,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_lat_g[:,0],mri_cgcm_tclw_frac_lat_g[:,1], '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tclw_frac_lat_g_enhanced[:,0],cccm_tclw_frac_lat_g_enhanced[:,1], ':r', label='CCCM')
ax2.plot(calipso_tclw_frac_lat_g[:,0],calipso_tclw_frac_lat_g[:,1], ':b', label='CALIPSO')

ax2.plot(gfdl4_tclw_frac_lat_g[:,0],gfdl4_tclw_frac_lat_g[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_lat_g[:,0],mri_tclw_frac_lat_g[:,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_lat_g[:,0],cam6_tclw_frac_lat_g[:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax1.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax1.set_ylabel('Cloud liquid Water Fraction')
ax2.set_ylabel('Cloud liquid Water Fraction')
ax2.set_xlabel('Latitude')

ax1.set_title ('2007 - 2008 Global Cloud Fraction vs Latitude')

ax1.grid(True)
ax2.grid(True)

plt.title('2007 to 2008 Global Liquid Fractions of Clouds vs Latitude')

plt.grid(True)
plt.savefig("2007_2008_tclw_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Global Comparison Ice Cloud Fractions with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_frac_lat_g[:,0],cccm_tciw_frac_lat_g[:,1], '--r', label='CCCM')
ax.plot(calipso_tclw_frac_lat_g[:,0],calipso_tclw_frac_lat_g[:,1], '--b', label='CALIPSO')
ax.plot(gfdl_hiram_tciw_frac_lat_g[:,0],gfdl_hiram_tciw_frac_lat_g[:,1], '--k', label='CMIP5-GFDL-HIRAM-AMIP')
ax.plot(gfdl4_tciw_frac_lat_g[:,0],gfdl4_tciw_frac_lat_g[:,1], '--g', label='CMIP6-GFDL-AM4-AMIP')
ax.plot(mri_cgcm_tciw_frac_lat_g[:,0],mri_cgcm_tciw_frac_lat_g[:,1], '--y', label='CMIP5-MRI_CGCM3-AMIP')
ax.plot(mri_tciw_frac_lat_g[:,0],mri_tciw_frac_lat_g[:,1], '--m', label='CMIP6-MRI_ESM2-AMIP')
ax.plot(cam6_tciw_frac_lat_g[:,0],cam6_tciw_frac_lat_g[:,1], '--c', label='CMIP6-CESM2-CAM6-AMIP')


ax.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Ice Water Fraction')

plt.title('2007 to 2008 Global Ice Fractions of Clouds vs Latitude')

plt.grid(True)
plt.savefig("2007_2008_tciw_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""




############################################################################### Global Altitude Plots

#---Plot Global Cloud Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tcc_alt_g[4:92,1],cccm_tcc_alt_g[4:92,0], ':r', label='CCCM')
ax1.plot(calipso_tcc_alt_g[:,1],calipso_tcc_alt_g[:,0], ':b', label='CALIPSO')
ax1.plot(ecmwf_tcc_alt_g[9:,1],ecmwf_tcc_alt_g[9:,0], '-k', label='ECMWF-ERA5')

ax1.plot(gfdl_hiram_tcc_alt_g[:23,1],gfdl_hiram_tcc_alt_g[:23,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tcc_alt_g[:25,1],mri_cgcm_tcc_alt_g[:25,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tcc_alt_g[4:92,1],cccm_tcc_alt_g[4:92,0], ':r', label='CCCM')
ax2.plot(calipso_tcc_alt_g[:,1],calipso_tcc_alt_g[:,0], ':b', label='CALIPSO')
ax2.plot(ecmwf_tcc_alt_g[9:,1],ecmwf_tcc_alt_g[9:,0], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tcc_alt_g[:23,1],gfdl4_tcc_alt_g[:23,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_alt_g[:42,1],mri_tcc_alt_g[:42,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tcc_alt_g[10:,1],cam6_tcc_alt_g[10:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.275));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.3));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Cloud Fraction')

#ax1.set_title('2007 to 2008 Global Cloud Fraction vs Altitude')

ax1.set_xlim(0, 0.35)
ax2.set_xlim(0, 0.35)

ax1.text(-0.1, 19, 'a)')
ax2.text(-0.05, 19, 'b)')

ax1.grid(True)
ax2.grid(True)

plt.savefig("2007_2008_tcc_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Cloud Liquid Water Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tclw_frac_alt_g[4:50,1],cccm_tclw_frac_alt_g[4:50,0], ':r', label='CCCM')
ax1.plot(calipso_tclw_frac_alt_g[:21,1],calipso_tclw_frac_alt_g[:21,0], ':b', label='CALIPSO')
ax1.plot(ecmwf_tclw_frac_alt_g[18:,1],ecmwf_tclw_frac_alt_g[18:,0], '-k', label='ECMWF-ERA5')

ax1.plot(gfdl_hiram_tclw_frac_alt_g[:18,1],gfdl_hiram_tclw_frac_alt_g[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_alt_g[:19,1],mri_cgcm_tclw_frac_alt_g[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tclw_frac_alt_g[4:50,1],cccm_tclw_frac_alt_g[4:50,0], ':r', label='CCCM')
ax2.plot(calipso_tclw_frac_alt_g[:21,1],calipso_tclw_frac_alt_g[:21,0], ':b', label='CALIPSO')
ax2.plot(ecmwf_tclw_frac_alt_g[18:,1],ecmwf_tclw_frac_alt_g[18:,0], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tclw_frac_alt_g[:18,1],gfdl4_tclw_frac_alt_g[:18,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_alt_g[:25,1],mri_tclw_frac_alt_g[:25,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_frac_alt_g[17:,1],cam6_tclw_frac_alt_g[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.275));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.3));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Liquid Water Fraction')
ax2.set_xlabel('Cloud Liquid Water Fraction')

#ax1.set_title('2007 to 2008 Cloud Liquid Water Fraction vs Altitude')

ax1.set_xlim(0, 0.35)
ax2.set_xlim(0, 0.35)

ax1.text(-0.1, 10, 'a)')
ax2.text(-0.05, 10, 'b)')

ax1.grid(True)

ax2.grid(True)
plt.savefig("2007_2008_tclw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""


#---Plot Global Cloud Liquid Water Content Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tclw_alt_g[4:50,1],cccm_tclw_alt_g[4:50,0], ':r', label='CCCM')
ax1.plot(calipso_tclw_alt_g[:21,1],calipso_tclw_alt_g[:21,0], ':b', label='CALIPSO')

ax1.plot(gfdl_hiram_tclw_alt_g[:18,1],gfdl_hiram_tclw_alt_g[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_alt_g[:19,1],mri_cgcm_tclw_alt_g[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tclw_alt_g[4:50,1],cccm_tclw_alt_g[4:50,0], ':r', label='CCCM')
ax2.plot(calipso_tclw_alt_g[:21,1],calipso_tclw_alt_g[:21,0], ':b', label='CALIPSO')

ax2.plot(gfdl4_tclw_alt_g[:18,1],gfdl4_tclw_alt_g[:18,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_alt_g[:25,1],mri_tclw_alt_g[:25,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_alt_g[17:,1],cam6_tclw_alt_g[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.225));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Liquid Water Content (g/kg)')
ax2.set_xlabel('Liquid Water Content (g/kg)')

#ax1.set_title('2007 to 2008 Cloud Liquid WaterContent vs Altitude')

ax1.set_xlim(0, 0.35)
ax2.set_xlim(0, 0.35)

ax1.text(-0.005, 9, 'a)')
ax2.text(-0.004, 9, 'b)')

ax1.grid(True)

ax2.grid(True)
plt.savefig("2007_2008_tclw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""


#---Plot Global Cloud Ice Water Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_frac_alt_g[4:92,1],cccm_tciw_frac_alt_g[4:92,0], '--r', label='CCCM')
ax.plot(calipso_tciw_frac_alt_g[:,1],calipso_tciw_frac_alt_g[:,0], '--b', label='CALIPSO')
ax.plot(ecmwf_tciw_frac_alt_g[10:,1],ecmwf_tciw_frac_alt_g[10:,0], '--k', label='ECMWF-ERA5')
ax.plot(gfdl4_tciw_frac_alt_g[:26,1],gfdl4_tciw_frac_alt_g[:26,0], '--g', label='CMIP6-GFDL-AM4-AMIP')
ax.plot(mri_tciw_frac_alt_g[:50,1],mri_tciw_frac_alt_g[:50,0], '--m', label='CMIP6-MRI-ESM2-AMIP')
ax.plot(cam6_tciw_frac_alt_g[7:,1],cam6_tciw_frac_alt_g[7:,0], '--c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Ice Water Fraction')

plt.title('2007 to 2008 Cloud Ice Water Fraction vs Altitude')

plt.grid(True)
plt.savefig("2007_2008_tciw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""


############################################################################### Southern Ocean Altitude Plots

#---Plot Southern Ocean Cloud Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tcc_alt_so[4:92,1],cccm_tcc_alt_so[4:92,0], ':r', label='CCCM')
ax1.plot(calipso_tcc_alt_so[:25,1],calipso_tcc_alt_so[:25,0], ':b', label='CALIPSO')
ax1.plot(ecmwf_tcc_alt_so[13:,1],ecmwf_tcc_alt_so[13:,0], '-k', label='ECMWF-ERA5')

ax1.plot(gfdl_hiram_tcc_alt_so[:23,1],gfdl_hiram_tcc_alt_so[:23,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tcc_alt_so[:25,1],mri_cgcm_tcc_alt_so[:25,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tcc_alt_so[4:92,1],cccm_tcc_alt_so[4:92,0], ':r', label='CCCM')
ax2.plot(calipso_tcc_alt_so[:25,1],calipso_tcc_alt_so[:25,0], ':b', label='CALIPSO')
ax2.plot(ecmwf_tcc_alt_so[13:,1],ecmwf_tcc_alt_so[13:,0], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tcc_alt_so[:23,1],gfdl4_tcc_alt_so[:23,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_alt_so[:42,1],mri_tcc_alt_so[:42,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tcc_alt_so[10:,1],cam6_tcc_alt_so[10:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.275));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.3));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Cloud Fraction')

#plt.title('2007 to 2008 Southern Ocean Cloud Fraction vs Altitude')

ax1.set_xlim(0, 0.65)
ax2.set_xlim(0, 0.65)

ax1.text(-0.11, 18.5, 'a)')
ax2.text(-0.09, 18.5, 'b)')

ax1.grid(True)
ax2.grid(True)
plt.savefig("2007_2008_tcc_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Cloud Liquid Water Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tclw_frac_alt_so[4:50,1],cccm_tclw_frac_alt_so[4:50,0], ':r', label='CCCM')
ax1.plot(calipso_tclw_frac_alt_so[:21,1],calipso_tclw_frac_alt_so[:21,0], ':b', label='CALIPSO')
ax1.plot(ecmwf_tclw_frac_alt_so[18:,1],ecmwf_tclw_frac_alt_so[18:,0], '-k', label='ECMWF-ERA5')

ax1.plot(gfdl_hiram_tclw_frac_alt_so[:18,1],gfdl_hiram_tclw_frac_alt_so[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_alt_so[:19,1],mri_cgcm_tclw_frac_alt_so[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tclw_frac_alt_so[4:50,1],cccm_tclw_frac_alt_so[4:50,0], ':r', label='CCCM')
ax2.plot(calipso_tclw_frac_alt_so[:21,1],calipso_tclw_frac_alt_so[:21,0], ':b', label='CALIPSO')
ax2.plot(ecmwf_tclw_frac_alt_so[18:,1],ecmwf_tclw_frac_alt_so[18:,0], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tclw_frac_alt_so[:18,1],gfdl4_tclw_frac_alt_so[:18,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_alt_so[:25,1],mri_tclw_frac_alt_so[:25,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_frac_alt_so[17:,1],cam6_tclw_frac_alt_so[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.275));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.3));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Liquid Water Fraction')
ax2.set_xlabel('Cloud Liquid Water Fraction')

#plt.title('2007 to 2008 Southern Ocean Cloud Fraction vs Altitude')

ax1.set_xlim(0, 0.6)
ax2.set_xlim(0, 0.6)

ax1.text(-0.10, 9, 'a)')
ax2.text(-0.05, 9, 'b)')

ax1.grid(True)
ax2.grid(True)
plt.savefig("2007_2008_tclw_frac_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot SO Cloud Liquid Water Content Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(cccm_tclw_alt_so[4:50,1],cccm_tclw_alt_so[4:50,0], ':r', label='CCCM')
ax1.plot(calipso_tclw_alt_so[:21,1],calipso_tclw_alt_so[:21,0], ':b', label='CALIPSO')

ax1.plot(gfdl_hiram_tclw_alt_so[:18,1],gfdl_hiram_tclw_alt_so[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_alt_so[:19,1],mri_cgcm_tclw_alt_so[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')

ax2.plot(cccm_tclw_alt_so[4:50,1],cccm_tclw_alt_so[4:50,0], ':r', label='CCCM')
ax2.plot(calipso_tclw_alt_so[:21,1],calipso_tclw_alt_so[:21,0], ':b', label='CALIPSO')

ax2.plot(gfdl4_tclw_alt_so[:18,1],gfdl4_tclw_alt_so[:18,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_alt_so[:25,1],mri_tclw_alt_so[:25,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_alt_so[17:,1],cam6_tclw_alt_so[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.225));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Liquid Water Content (g/kg)')
ax2.set_xlabel('Liquid Water Content (g/kg)')

plt.title('2007 to 2008 Southern Ocean Cloud Content vs Altitude')

ax1.set_xlim(0, 0.4)
ax2.set_xlim(0, 0.4)

ax1.text(-0.006, 9, 'a)')
ax2.text(-0.006, 9, 'b)')

ax1.grid(True)
ax2.grid(True)
plt.savefig("2007_2008_tclw_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""


#---Plot Southern Ocean Cloud Ice Water Content Fraction Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(calipso_tciw_frac_alt_so[:,1],calipso_tciw_frac_alt_so[:,0], '-r', label='CALIPSO')
ax.plot(ecmwf_tciw_frac_alt_so[:,1],ecmwf_tciw_frac_alt_so[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tciw_frac_alt_so[:,1],gfdl3_tciw_frac_alt_so[:,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tciw_frac_alt_so[:,1],gfdl4_tciw_frac_alt_so[:,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tciw_frac_alt_so[:,1],cam_tciw_frac_alt_so[:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tciw_frac_alt_so[:,1],cam6_tciw_frac_alt_so[:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Ice Water Fraction')

plt.title('2007 to 2008 Cloud Ice Water Fraction over the Southern Ocean vs Altitude')

plt.grid(True)
plt.savefig("2007_2008_tciw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""
############################################################################### Contour plots


#---Combined Grid tclw_frac---#
"""
fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(10, 10))

ax[0, 1].contourf(ecmwf_lat, ecmwf_alt[19:], ecmwf_tclw_frac_alt_lat[19:], vmin=0, vmax=0.5)
ax[0, 0].set_ylabel('Altitude (km)')
ax[0, 1].set_title('c) ECMWF-ERA5')
ecmwf_temp = ax[0, 1].contour(ecmwf_lat, ecmwf_alt[19:], (ecmwf_temp_alt_lat[19:] - 273.15), colors='grey')
ecmwf_temp.collections[5].set_linewidth(2)
ecmwf_temp.collections[5].set_color('white')
ax[0, 1].clabel(ecmwf_temp, inline=1, fontsize=10)

#ax[0, 0].contourf(cccm_lat, cccm_alt[0:17], cccm_tclw_frac_alt_lat[0:17], vmin=0, vmax=0.6)
#ax[0, 0].set_title('a) CCCM')
#cccm_temp = ax[0, 0].contour(cccm_lat, cccm_alt_temp[1:7], (cccm_temp_alt_lat[1:7] - 273.15), colors='grey')
#ax[0, 0].clabel(cccm_temp, inline=1, fontsize=10)

ax[0, 0].contourf(calipso_lat, calipso_alt[0:16], calipso_tclw_frac_alt_lat[0:16], vmin=0, vmax=0.5)
ecmwf_temp = ax[0, 0].contour(ecmwf_lat, ecmwf_alt[19:], (ecmwf_temp_alt_lat[19:] - 273.15), colors='grey')
ax[0, 0].set_title('a) CALIPSO-GOCCP')
ecmwf_temp.collections[5].set_linewidth(2)
ecmwf_temp.collections[5].set_color('white')
ax[0, 0].clabel(ecmwf_temp, inline=1, fontsize=10)

#calipso_temp = ax[0, 0].contour(calipso_lat, calipso_alt_temp[1:7], (calipso_temp_alt_lat[1:7] - 273.15), colors='grey')
#ax[0, 0].clabel(calipso_temp, inline=1, fontsize=10)

ax[1, 0].contourf(gfdl_hiram_lat, gfdl_hiram_alt[0:17], gfdl_hiram_tclw_frac_alt_lat[0:17], vmin=0, vmax=0.5)
ax[1, 0].set_title('d) CMIP5-GFDL-HIRAM')
ax[1, 0].set_ylabel('Altitude (km)')
gfdl_hiram_temp = ax[1, 0].contour(gfdl_hiram_lat, gfdl_hiram_alt_temp[1:7], (gfdl_hiram_temp_alt_lat[1:7] - 273.15), colors='grey')
gfdl_hiram_temp.collections[7].set_linewidth(2)
gfdl_hiram_temp.collections[7].set_color('white')
ax[1, 0].clabel(gfdl_hiram_temp, inline=1, fontsize=10)

ax[1, 1].contourf(gfdl4_lat, gfdl4_alt[0:18], gfdl4_tclw_frac_alt_lat[0:18], vmin=0, vmax=0.5)
ax[1, 1].set_title('e) CMIP6-GFDL-AM4')
gfdl4_temp = ax[1, 1].contour(gfdl4_lat, gfdl4_alt_temp[1:7], (gfdl4_temp_alt_lat[1:7] - 273.15), colors='grey')
gfdl4_temp.collections[5].set_linewidth(2)
gfdl4_temp.collections[5].set_color('white')
ax[1, 1].clabel(gfdl4_temp, inline=1, fontsize=10)

ax[2, 0].contourf(mri_cgcm_lat, mri_cgcm_alt[0:18], mri_cgcm_tclw_frac_alt_lat[0:18], vmin=0, vmax=0.5)
ax[2, 0].set_title('f) CMIP5-MRI-CGCM3')
ax[2, 0].set_xlabel('Latitude')
ax[2, 0].set_ylabel('Altitude (km)')
mri_cgcm_temp = ax[2, 0].contour(mri_cgcm_lat, mri_cgcm_alt_temp[1:7], (mri_cgcm_temp_alt_lat[1:7] - 273.15), colors='grey')
mri_cgcm_temp.collections[5].set_linewidth(2)
mri_cgcm_temp.collections[5].set_color('white')
ax[2, 0].clabel(mri_cgcm_temp, inline=1, fontsize=10)

cont=ax[2, 1].contourf(mri_lat, mri_alt[0:25], mri_tclw_frac_alt_lat[0:25], vmin=0, vmax=0.5)
ax[2, 1].set_title('g) CMIP6-MRI_ESM2')
mri_temp = ax[2, 1].contour(mri_lat, mri_alt_temp[1:7], (mri_temp_alt_lat[1:7] - 273.15), colors='grey')
mri_temp.collections[6].set_linewidth(2)
mri_temp.collections[6].set_color('white')
ax[2, 1].clabel(mri_temp, inline=1, fontsize=10)


ax[3, 1].contourf(cam6_lat, cam6_alt[17:32], cam6_tclw_frac_alt_lat[17:32], vmin=0, vmax=0.5)
ax[3, 1].set_title('h) CMIP6-CESM2.1-CAM6')
ax[3, 1].set_ylabel('Altitude (km)')
ax[3, 1].set_xlabel('Latitude')
cam6_temp = ax[3, 1].contour(cam6_lat, cam6_alt[17:32], (cam6_temp_alt_lat[17:32] - 273.15), colors='grey')
cam6_temp.collections[5].set_linewidth(2)
cam6_temp.collections[5].set_color('white')
ax[3, 1].clabel(cam6_temp, inline=1, fontsize=10)


ax[3, 0].remove()  # don't display empty ax
ax[0, 2].remove()  # don't display empty ax
ax[1, 2].remove()  # don't display empty ax
ax[2, 2].remove()  # don't display empty ax
ax[3, 2].remove()  # don't display empty ax

cbaxes = fig.add_axes([0.7, 0.5, 0.03, 0.3]) #(x-position, y-position, thickness, length)
cbar = fig.colorbar(cont, cax=cbaxes)
cbar.set_clim(0, 0.5)
cbar.set_label('Cloud liquid Water Fraction')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)
plt.savefig("2007_2008_contour_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Combined Grid tciw_frac---#

fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(10, 10))

ax[0, 1].contourf(ecmwf_lat, ecmwf_alt[9:], ecmwf_tciw_frac_alt_lat[9:], vmin=0, vmax=0.4)
ax[0, 1].set_xlabel('Latitude')
ax[0, 0].set_ylabel('Altitude (km)')
ax[0, 1].set_title('c) ECMWF-ERA5')
ecmwf_temp = ax[0, 1].contour(ecmwf_lat, ecmwf_alt[9:], (ecmwf_temp_alt_lat[9:] - 273.15), colors='grey')
ecmwf_temp.collections[6].set_linewidth(2)
ecmwf_temp.collections[6].set_color('white')
ax[0, 1].clabel(ecmwf_temp, inline=1, fontsize=10)

#ax[0, 0].contourf(cccm_lat, cccm_alt[0:17], cccm_tciw_frac_alt_lat[0:17], vmin=0, vmax=0.6)
#ax[0, 0].set_title('a) CCCM')
#cccm_temp = ax[0, 0].contour(cccm_lat, cccm_alt_temp[1:7], (cccm_temp_alt_lat[1:7] - 273.15), colors='grey')
#ax[0, 0].clabel(cccm_temp, inline=1, fontsize=10)

ax[0, 0].contourf(calipso_lat, calipso_alt, calipso_tciw_frac_alt_lat, vmin=0, vmax=0.4)
ax[0, 0].set_title('b) CALIPSO-GOCCP')
ecmwf_temp = ax[0, 0].contour(ecmwf_lat, ecmwf_alt[19:], (ecmwf_temp_alt_lat[19:] - 273.15), colors='grey')
ecmwf_temp.collections[5].set_linewidth(2)
ecmwf_temp.collections[5].set_color('white')
ax[0, 0].clabel(ecmwf_temp, inline=1, fontsize=10)

#calipso_temp = ax[0, 0].contour(calipso_lat, calipso_alt_temp[1:7], (calipso_temp_alt_lat[1:7] - 273.15), colors='grey')
#ax[0, 0].clabel(calipso_temp, inline=1, fontsize=10)

ax[1, 0].contourf(gfdl_hiram_lat, gfdl_hiram_alt[:26], gfdl_hiram_tciw_frac_alt_lat[:26], vmin=0, vmax=0.4)
ax[1, 0].set_title('a) CMIP5-GFDL-HIRAM')
gfdl_hiram_temp = ax[1, 0].contour(gfdl_hiram_lat, gfdl_hiram_alt_temp[1:13], (gfdl_hiram_temp_alt_lat[1:13] - 273.15), colors='grey')
gfdl_hiram_temp.collections[6].set_linewidth(2)
gfdl_hiram_temp.collections[6].set_color('white')
ax[1, 0].clabel(gfdl_hiram_temp, inline=1, fontsize=10)

ax[1, 1].contourf(gfdl4_lat, gfdl4_alt[:26], gfdl4_tciw_frac_alt_lat[:26], vmin=0, vmax=0.4)
ax[1, 1].set_title('b) CMIP6-GFDL-AM4')
gfdl4_temp = ax[1, 1].contour(gfdl4_lat, gfdl4_alt_temp[1:13], (gfdl4_temp_alt_lat[1:13] - 273.15), colors='grey')
gfdl4_temp.collections[6].set_linewidth(2)
gfdl4_temp.collections[6].set_color('white')
ax[1, 1].clabel(gfdl4_temp, inline=1, fontsize=10)

ax[2, 0].contourf(mri_cgcm_lat, mri_cgcm_alt[:30], mri_cgcm_tciw_frac_alt_lat[:30], vmin=0, vmax=0.4)
ax[2, 0].set_title('d) CMIP5-MRI-CGCM3')
ax[2, 0].set_ylabel('Altitude (km)')
mri_cgcm_temp = ax[2, 0].contour(mri_cgcm_lat, mri_cgcm_alt_temp[1:13], (mri_cgcm_temp_alt_lat[1:13] - 273.15), colors='grey')
mri_cgcm_temp.collections[6].set_linewidth(2)
mri_cgcm_temp.collections[6].set_color('white')
ax[2, 0].clabel(mri_cgcm_temp, inline=1, fontsize=10)

cont=ax[2, 1].contourf(mri_lat, mri_alt[:44], mri_tciw_frac_alt_lat[:44], vmin=0, vmax=0.4)
ax[2, 1].set_title('e) CMIP6-MRI_ESM2')
mri_temp = ax[2, 1].contour(mri_lat, mri_alt_temp[1:12], (mri_temp_alt_lat[1:12] - 273.15), colors='grey')
mri_temp.collections[6].set_linewidth(2)
mri_temp.collections[6].set_color('white')
ax[2, 1].clabel(mri_temp, inline=1, fontsize=10)

ax[3, 1].contourf(cam6_lat, cam6_alt[7:32], cam6_tciw_frac_alt_lat[7:32], vmin=0, vmax=0.4)
ax[3, 1].set_title('g) CMIP6-CESM2.1-CAM6')
ax[3, 0].set_ylabel('Altitude (km)')
ax[3, 1].set_xlabel('Latitude')
cam6_temp = ax[3, 1].contour(cam6_lat, cam6_alt[7:32], (cam6_temp_alt_lat[7:32] - 273.15), colors='grey')
cam6_temp.collections[6].set_linewidth(2)
cam6_temp.collections[6].set_color('white')
ax[3, 1].clabel(cam6_temp, inline=1, fontsize=10)


ax[1, 2].remove()  # don't display empty ax
ax[2, 2].remove()  # don't display empty ax
ax[3, 2].remove()  # don't display empty ax
ax[3, 0].remove()  # don't display empty ax
ax[0, 2].remove()  # don't display empty ax

cbaxes = fig.add_axes([0.7, 0.5, 0.03, 0.3]) 
cbar = fig.colorbar(cont, cax=cbaxes)
cbar.set_clim(0, 0.4)
cbar.set_label('Cloud Ice Water Fraction')
plt.tight_layout(pad=0.4, w_pad=0.6, h_pad=3)
plt.savefig("2007_2008_contour_tciw.svg", format="svg", bbox_inches='tight')
plt.show()












#---Plot CALIPSO Cloud Liquid Water Fraction---#

"""
plt.subplots()

plt.contourf(calipso_lat, calipso_alt[:16], calipso_tclw_frac_alt_lat[:16])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.title('2007 to 2008 CALIPSO Cloud Liquid Water Fraction')

plt.colorbar().set_label('Liquid Cloud Fraction')
plt.show()
plt.savefig("2007_2008_contour_calipso_tclw.svg", format="svg", bbox_inches='tight')


"""

#---Plot CALIPSO Cloud Ice Water Fraction---#

"""
plt.subplots()

plt.contourf(calipso_lat, calipso_alt, calipso_tciw_frac_alt_lat)
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.title('2007 to 2008 CALIPSO Cloud Ice Water Fraction')

plt.colorbar().set_label('Ice Cloud Fraction')
plt.show()
plt.savefig("2007_2008_contour_calipso_tciw.svg", format="svg", bbox_inches='tight')
"""



#---Plot ECMWF Cloud Liquid Water Fraction---#

"""
plt.subplots()

plt.contourf(ecmwf_lat, ecmwf_alt[0:18], ecmwf_tclw_frac_alt_lat[0:18])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Liquid Cloud Fraction')

plt.title('2007 to 2008 ECMWF-ERA5 Cloud Liquid Water Fraction')

plt.savefig("2007_2008_contour_ecmwf_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot ECMWF Cloud Ice Water Fraction---#

"""
plt.subplots()

plt.contourf(ecmwf_lat, ecmwf_alt[:26], ecmwf_tciw_frac_alt_lat[:26])
plt.clim(0, 0.5)

plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.colorbar().set_label('Ice Cloud Fraction')

plt.title('2007 to 2008 ECMWF-ERA5 Cloud Ice Water Fraction')

plt.savefig("2007_2008_contour_ecmwf_tciw.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 GFDL-HIRAM-AMIP Cloud Liquid Water Fraction')

plt.savefig("2007_2008_contour_gfdl_hiram_tclw.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 GFDL-HIRAM-AMIP Cloud Ice Water Fraction')

plt.savefig("2007_2008_contour_gfdl_hiram.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 GFDL-AM4-AMIP Cloud Liquid Water Fraction')

plt.savefig("2007_2008_contour_gfdl4_tclw.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 GFDL-AM4-AMIP Cloud Ice Water Fraction')

plt.savefig("2007_2008_contour_gfdl4_tciw.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 MRI-CGCM3-AMIP Cloud Liquid Water Fraction')

plt.savefig("2007_2008_contour_mri_cgcm_tclw.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 MRI-CGCM3-AMIP Cloud Ice Water Fraction')

plt.savefig("2007_2008_contour_mri_cgcm_tciw.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 MRI_ESM2-AMIP Cloud Liquid Water Fraction')

plt.savefig("2007_2008_contour_mri_tclw.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 MRI_ESM2-AMIP Cloud Ice Water Fraction')

plt.savefig("2007_2008_contour_mri_tciw.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 CESM2-CAM6-AMIP Cloud Liquid Water Fraction')

plt.savefig("2007_2008_contour_cam6_tclw.svg", format="svg", bbox_inches='tight')
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

plt.title('2007 to 2008 CESM2-CAM6-AMIP Cloud Ice Water Fraction')

plt.savefig("2007_2008_contour_cam6_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

