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
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
#a = h5py.File('2001_2005_ecmwf.h5', 'r')

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
"""


# Home PC
#ECMWF Data
#os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
#a = h5py.File('2001_2005_ecmwf.h5', 'r')

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

"""


############################################################################### ECMWF Data
"""
#---ECMWF Global Latitude Data---#

ecmwf_tcc_lat_g = a['tcc'][:] # 0-1
ecmwf_tclw_lat_g = a['tclw'][:] #kgm^-2
ecmwf_tciw_lat_g = a['tciw'][:] #kgm^-2

#---ECMWF Southern Ocean Latitude Data---#

ecmwf_tcc_lat_so = ecmwf_tcc_lat_g[ecmwf_tcc_lat_g[:,0]>=-70]
ecmwf_tcc_lat_so = ecmwf_tcc_lat_so[ecmwf_tcc_lat_so[:,0]<=-50] # 0-1

ecmwf_tclw_lat_so = ecmwf_tclw_lat_g[ecmwf_tclw_lat_g[:,0]>=-70]
ecmwf_tclw_lat_so = ecmwf_tclw_lat_so[ecmwf_tclw_lat_so[:,0]<=-50] #kgm^-2

ecmwf_tciw_lat_so = ecmwf_tciw_lat_g[ecmwf_tciw_lat_g[:,0]>=-70]
ecmwf_tciw_lat_so = ecmwf_tciw_lat_so[ecmwf_tciw_lat_so[:,0]<=-50] #kgm^-2

#---ECMWF Phase Fractions---#

ecmwf_tclw_frac_lat_g = (ecmwf_tclw_lat_g[:,1] / (ecmwf_tclw_lat_g[:,1] + ecmwf_tciw_lat_g[:,1])) * ecmwf_tcc_lat_g[:,1]
ecmwf_tciw_frac_lat_g = (ecmwf_tciw_lat_g[:,1] / (ecmwf_tclw_lat_g[:,1] + ecmwf_tciw_lat_g[:,1])) * ecmwf_tcc_lat_g[:,1]

ecmwf_tclw_frac_lat_g = np.vstack((ecmwf_tclw_lat_g[:,0], ecmwf_tclw_frac_lat_g)).T
ecmwf_tciw_frac_lat_g = np.vstack((ecmwf_tciw_lat_g[:,0], ecmwf_tciw_frac_lat_g)).T

#---ECMWF lat-alt contour data---#

ecmwf_tcc_alt_lat = b['cf_alt_lat'][:] #kg/kg
ecmwf_tclw_alt_lat = b['lw_alt_lat'][:] #kg/kg
ecmwf_tciw_alt_lat = b['iw_alt_lat'][:] #kg/kg
ecmwf_temp_alt_lat = b['temp_alt_lat'][:] #kg/kg
ecmwf_lat = b['lat'][:]
ecmwf_alt = b['alt'][:]

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

ecmwf_tclw_frac_temp_g = (ecmwf_tclw_temp_g[:,1] / (ecmwf_tclw_temp_g[:,1] + ecmwf_tciw_temp_g[:,1])) * ecmwf_tcc_temp_g[:,1]
ecmwf_tciw_frac_temp_g = (ecmwf_tciw_temp_g[:,1] / (ecmwf_tclw_temp_g[:,1] + ecmwf_tciw_temp_g[:,1])) * ecmwf_tcc_temp_g[:,1]

ecmwf_tclw_frac_temp_g = np.vstack((ecmwf_tclw_temp_g[:,0], ecmwf_tclw_frac_temp_g)).T
ecmwf_tciw_frac_temp_g = np.vstack((ecmwf_tciw_temp_g[:,0], ecmwf_tciw_frac_temp_g)).T


#---ECMWF Phase Profile Fractions---#

ecmwf_tclw_frac_alt_g = (ecmwf_tclw_alt_g[:,1] / (ecmwf_tclw_alt_g[:,1] + ecmwf_tciw_alt_g[:,1])) * ecmwf_tcc_alt_g[:,1]
ecmwf_tciw_frac_alt_g = (ecmwf_tciw_alt_g[:,1] / (ecmwf_tclw_alt_g[:,1] + ecmwf_tciw_alt_g[:,1])) * ecmwf_tcc_alt_g[:,1]

ecmwf_tclw_frac_alt_g = np.vstack((ecmwf_tclw_alt_g[:,0], ecmwf_tclw_frac_alt_g)).T
ecmwf_tciw_frac_alt_g = np.vstack((ecmwf_tciw_alt_g[:,0], ecmwf_tciw_frac_alt_g)).T

ecmwf_tcc_alt_g = a['cf'][8:] # 0-1
ecmwf_tclw_alt_g = a['lw'][14:] #kg/kg
ecmwf_tciw_alt_g = a['iw'][9:] #kg/kg


#---ECMWF Southern Ocean Profile---#

ecmwf_tcc_alt_so = a['cf_so'][11:] # 0-1
ecmwf_tclw_alt_so = a['lw_so'][18:] #kg/kg
ecmwf_tciw_alt_so = a['iw_so'][11:] #kg/kg
ecmwf_temp_alt_so = a['temp_so'][:] #K
ecmwf_plevel_alt_so = a['pressure_so'][:] #hPa

ecmwf_tcc_temp_so = a['cf_t_so'][:] # 0-1
ecmwf_tclw_temp_so = a['lw_t_so'][:] #kg/kg
ecmwf_tciw_temp_so = a['iw_t_so'][:] #kg/kg


ecmwf_tclw_frac_temp_so = (ecmwf_tclw_temp_so[:,1] / (ecmwf_tclw_temp_so[:,1] + ecmwf_tciw_temp_so[:,1])) * ecmwf_tcc_temp_so[:,1]
ecmwf_tciw_frac_temp_so = (ecmwf_tciw_temp_so[:,1] / (ecmwf_tclw_temp_so[:,1] + ecmwf_tciw_temp_so[:,1])) * ecmwf_tcc_temp_so[:,1]

ecmwf_tclw_frac_temp_so = np.vstack((ecmwf_tclw_temp_so[:,0], ecmwf_tclw_frac_temp_so)).T
ecmwf_tciw_frac_temp_so = np.vstack((ecmwf_tciw_temp_so[:,0], ecmwf_tciw_frac_temp_so)).T

"""

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

cam5_tclw_frac_alt_lat = (cam5_tclw_alt_lat / (cam5_tclw_alt_lat + cam5_tciw_alt_lat)) * (cam5_tcc_alt_lat)
cam5_tciw_frac_alt_lat = (cam5_tciw_alt_lat / (cam5_tclw_alt_lat + cam5_tciw_alt_lat)) * (cam5_tcc_alt_lat)


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


############################################################################### End importing Data

end = time.time()
print('Importing data took:', end - start, 's')


"""

os.chdir('c:/Users/toha006/University/University/MSc/Models/Images/Meeting 1.7/2001 - 2005')

os.chdir('E:/University/University/MSc/Models/Images')
"""
############################################################################### Temperature Profiles


#---Plot Liquid Cloud Fraction with Temperature---#


"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))


ax1.plot(gfdl_hiram_tclw_frac_temp_g[:17,0],gfdl_hiram_tclw_frac_temp_g[:17,1], '-g', label='Global - CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_g[:18,0],mri_cgcm_tclw_frac_temp_g[:18,1], '-m', label='Global - CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_temp_g[:12,0],cam5_tclw_frac_temp_g[:12,1], '-c', label='Global - CMIP5-CESM1-CAM5-AMIP')

ax2.plot(gfdl4_tclw_frac_temp_g[:18,0],gfdl4_tclw_frac_temp_g[:18,1], '-g', label='Global - CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_g[:26,0],mri_tclw_frac_temp_g[:26,1], '-m', label='Global - CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_g[19:,0],cam6_tclw_frac_temp_g[19:,1], '-c', label='Global - CMIP6-CESM2-CAM6-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Cloud Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Cloud Fraction')

ax1.set_title('2001 to 2005 Global Liquid Cloud Fraction vs Temperature')

ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_global_liquid_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Southern Ocean Liquid Cloud Fraction with Temperature---#


"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))


ax1.plot(gfdl_hiram_tclw_frac_temp_so[:17,0],gfdl_hiram_tclw_frac_temp_so[:17,1], '-g', label='Global - CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_temp_so[:18,0],mri_cgcm_tclw_frac_temp_so[:18,1], '-m', label='Global - CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_temp_so[:12,0],cam5_tclw_frac_temp_so[:12,1], '-c', label='Global - CMIP5-CESM1-CAM5-AMIP')

ax2.plot(gfdl4_tclw_frac_temp_so[:18,0],gfdl4_tclw_frac_temp_so[:18,1], '-g', label='Global - CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_temp_so[:26,0],mri_tclw_frac_temp_so[:26,1], '-m', label='Global - CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_temp_so[19:,0],cam6_tclw_frac_temp_so[19:,1], '-c', label='Global - CMIP6-CESM2-CAM6-AMIP')

ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Liquid Cloud Fraction')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Liquid Cloud Fraction')

ax1.set_title('2001 to 2005 Southern Ocean Liquid Cloud Fraction vs Temperature')

ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_so_liquid_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""

############################################################################### Global Latitude Plots

#---Plot Global Cloud Fraction with Latitude---#

"""
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True,figsize=(5, 6))

ax1.plot(gfdl_hiram_tcc_lat_g[:,0],gfdl_hiram_tcc_lat_g[:,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tcc_lat_g[:,0],mri_cgcm_tcc_lat_g[:,1], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tcc_lat_g[:,0],cam5_tcc_lat_g[:,1], '-c', label='CMIP5-CESM1-CAM5-AMIP')

ax2.plot(gfdl4_tcc_lat_g[:,0],gfdl4_tcc_lat_g[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_lat_g[:,0],mri_tcc_lat_g[:,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tcc_lat_g[:,0],cam6_tcc_lat_g[:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax1.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));
ax2.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Cloud Fraction')
ax2.set_xlabel('Latitude')

ax1.set_title ('2001 to 2005 Global Cloud Fraction vs Latitude')

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

ax2.plot(gfdl4_tclw_frac_lat_g[:,0],gfdl4_tclw_frac_lat_g[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_lat_g[:,0],mri_tclw_frac_lat_g[:,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax2.plot(cam6_tclw_frac_lat_g[:,0],cam6_tclw_frac_lat_g[:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')

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

ax1.plot(gfdl_hiram_tcc_alt_g[:23,1],gfdl_hiram_tcc_alt_g[:23,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tcc_alt_g[:25,1],mri_cgcm_tcc_alt_g[:25,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tcc_alt_g[:16,1],cam5_tcc_alt_g[:16,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')

ax2.plot(gfdl4_tcc_alt_g[:23,1],gfdl4_tcc_alt_g[:23,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_alt_g[:42,1],mri_tcc_alt_g[:42,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tcc_alt_g[10:,1],cam6_tcc_alt_g[10:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Cloud Fraction')

ax1.set_title('2001 to 2005 Global Cloud Fraction vs Altitude')

ax1.set_xlim(0, 0.35)
ax2.set_xlim(0, 0.35)

ax1.grid(True)
ax2.grid(True)

plt.savefig("2001_2005_tcc_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Cloud Liquid Water Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(gfdl_hiram_tclw_frac_alt_g[:18,1],gfdl_hiram_tclw_frac_alt_g[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_alt_g[:19,1],mri_cgcm_tclw_frac_alt_g[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_alt_g[:13,1],cam5_tclw_frac_alt_g[:13,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')

ax2.plot(gfdl4_tclw_frac_alt_g[:19,1],gfdl4_tclw_frac_alt_g[:19,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_alt_g[:26,1],mri_tclw_frac_alt_g[:26,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_frac_alt_g[17:,1],cam6_tclw_frac_alt_g[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Liquid Water Fraction')
ax2.set_xlabel('Cloud Liquid Water Fraction')

ax1.set_title('2001 to 2005 Global Cloud Liquid Water Fraction vs Altitude')

ax1.set_xlim(0, 0.35)
ax2.set_xlim(0, 0.35)

ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_tclw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Cloud Ice Water Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(ecmwf_tciw_frac_alt_g[10:,1],ecmwf_tciw_frac_alt_g[10:,0], '--k', label='ECMWF-ERA5')
ax.plot(gfdl4_tciw_frac_alt_g[:26,1],gfdl4_tciw_frac_alt_g[:26,0], '--g', label='CMIP6-GFDL-AM4-AMIP')
ax.plot(mri_tciw_frac_alt_g[:50,1],mri_tciw_frac_alt_g[:50,0], '--m', label='CMIP6-MRI-ESM2-AMIP')
ax.plot(cam6_tciw_frac_alt_g[7:,1],cam6_tciw_frac_alt_g[7:,0], '--c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Ice Water Fraction')

plt.title('2001 to 2005 Cloud Ice Water Fraction vs Altitude')

plt.grid(True)
plt.savefig("2001_2005_tciw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""


############################################################################### Southern Ocean Altitude Plots

#---Plot Southern Ocean Cloud Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(gfdl_hiram_tcc_alt_so[:23,1],gfdl_hiram_tcc_alt_so[:23,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tcc_alt_so[:25,1],mri_cgcm_tcc_alt_so[:25,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tcc_alt_so[:16,1],cam5_tcc_alt_so[:16,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')

ax2.plot(gfdl4_tcc_alt_so[:23,1],gfdl4_tcc_alt_so[:23,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_alt_so[:42,1],mri_tcc_alt_so[:42,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tcc_alt_so[10:,1],cam6_tcc_alt_so[10:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Cloud Fraction')

ax1.set_title('2001 to 2005 Southern Ocean Cloud Fraction vs Altitude')

ax1.set_xlim(0, 0.65)
ax2.set_xlim(0, 0.65)

ax1.grid(True)
ax2.grid(True)

plt.savefig("2001_2005_tcc_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Cloud Liquid Water Fraction Altitude Profile---#

"""
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(5, 6))

ax1.plot(gfdl_hiram_tclw_frac_alt_so[:18,1],gfdl_hiram_tclw_frac_alt_so[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_alt_so[:19,1],mri_cgcm_tclw_frac_alt_so[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(cam5_tclw_frac_alt_so[:13,1],cam5_tclw_frac_alt_so[:13,0], '-c', label='CMIP5-CESM1-CAM5-AMIP')

ax2.plot(gfdl4_tclw_frac_alt_so[:19,1],gfdl4_tclw_frac_alt_so[:19,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_alt_so[:26,1],mri_tclw_frac_alt_so[:26,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_frac_alt_so[17:,1],cam6_tclw_frac_alt_so[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.2));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.2));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Liquid Water Fraction')
ax2.set_xlabel('Cloud Liquid Water Fraction')

ax1.set_title('2001 to 2005 Southern Ocean Cloud Liquid Water Fraction vs Altitude')

ax1.set_xlim(0, 0.35)
ax2.set_xlim(0, 0.35)

ax1.grid(True)
ax2.grid(True)
plt.savefig("2001_2005_tclw_frac_alt_so.svg", format="svg", bbox_inches='tight')

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

plt.title('2001 to 2005 Cloud Ice Water Fraction over the Southern Ocean vs Altitude')

plt.grid(True)
plt.savefig("2001_2005_tciw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""
############################################################################### Contour plots


#---Plot ECMWF Cloud Liquid Water Fraction---#

"""
plt.subplots()

plt.contourf(ecmwf_lat, ecmwf_alt[0:18], ecmwf_tclw_frac_alt_lat[0:18])
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

plt.contourf(ecmwf_lat, ecmwf_alt[:26], ecmwf_tciw_frac_alt_lat[:26])
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

