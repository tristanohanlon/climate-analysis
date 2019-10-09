# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 10:06:47 2019

@author: toha006
"""
import numpy as np

model_dict = {
    "CMIP5-CESM1-CAM5" : "_Amon_CESM1-CAM5_amip_r1i1p1_197901-200512.nc",
    "CMIP5-GFDL-HIRAM-C360" : "_Amon_GFDL-HIRAM-C360_amip_r1i1p1_199901-200812.nc",
    "CMIP5-GISS-E2R" : "_Amon_GISS-E2-R_amip_r1i1p1_195101-201012.nc",
    "CMIP5-IPSL-CM5A-LR" : "_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc",
    "CMIP5-MIROC5" : "_Amon_MIROC5_amip_r1i1p1_199901-200812.nc",
    "CMIP5-MRI-CGCM3" : "_Amon_MRI-CGCM3_amip_r1i1p1_199901-201002.nc",
    
    "CMIP6-CESM2-CAM6" : "_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc",
    "CMIP6-GFDL-AM4" : "_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc",
    "CMIP6-GISS-E21G" : "_Amon_GISS-E2-1-G_amip_r1i1p1f1_gn_200101-201412.nc",
    "CMIP6-IPSL-CM6A-LR" : "_Amon_IPSL-CM6A-LR_amip_r1i1p1f1_gr_197901-201412.nc",
    "CMIP6-MIROC6" : "_Amon_MIROC6_amip_r1i1p1f1_gn_199901-201412.nc",
    "CMIP6-MRI-ESM2" : "_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_199901-201412.nc",
    
    "ECMWF" : "pressure_levels/2008_ECMWF_amon_plevels_T_cc_clw_ciw.nc"
   
        }

satellite_dict = {
    "CCCM" : "CCCM/CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905905.20060731.hdf",
    "MISR" : "MISR/MIL3YCFA/MISR_AM1_CFbA_2000.hdf",
    "CERES" : "CERES/CER_SSF1deg-Month_Terra-MODIS_Edition4A_400405.200607"
        }


home = 'E:/University/University/MSc/Models/'
uni = 'C:/Users/toha006/University/University/MSc/Models/'
hdd = 'D:/MSc/Models/'
laptop = 'C:/Users/tristan/University/University/MSc/Models/'

min_lat = -75.0
max_lat = 75.0
lat_division = 0.5
lat = np.arange(min_lat, max_lat, lat_division)

min_alt = 0.5
max_alt = 20
alt_division = 0.5
alt = np.arange(min_alt, max_alt, alt_division)

min_liq_alt = 0.5
max_liq_alt = 7
liq_alt_division = 0.25
liq_alt = np.arange(min_liq_alt, max_liq_alt, liq_alt_division)
