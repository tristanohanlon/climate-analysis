# -*- coding: utf-8 -*-
"""

@author: Tristan O'Hanlon

"""
import matplotlib.pyplot as plt
import h5py
import os
import numpy as np

#---Importing Data from Reduced Datasets---#

class Model:
    def __init__( self, path, name ):
        self.name = name
        a = h5py.File(path, 'r')
       
        #---universal data---#

        alt = a['alt'][:] # 0 to 20km
        liq_alt = a['liq_alt'][:] # 0 to 11km
        lat = a['lat'][:] # -75 to 75 degrees latitude
        
        #---temperature data---#

        self.ta_liq_g = a['ta_liq_g'][:] # global layer temperature corresponding to liq_alt
        self.ta_liq_so = a['ta_liq_so'][:] # southern ocean layer temperature corresponding to liq_alt

        #---zonal data---#

        self.clt = a['clt'][:] # total cloud fraction corresponding to lat

        #---alt-lat contour data---#

        self.cl_alt_lat = a['cl_alt_lat'][:] # total cloud fraction corresponding to alt and lat
        self.clw_alt_lat = a['clw_alt_lat'][:] # cloud liquid water fraction corresponding to liq_alt and lat
        self.cli_alt_lat = a['cli_alt_lat'][:] # cloud ice water fraction corresponding to alt and lat
        self.ta_alt_lat = a['ta_alt_lat'][:] # temperature corresponding to alt and lat

        #---global layer---#

        self.cl_g = a['cl_g'][:] # global layer total cloud fraction corresponding to alt
        self.clw_g = a['clw_g'][:] # global layer cloud liquid water fraction corresponding to liq_alt
        self.cli_g = a['cli_g'][:] # global layer cloud ice water fraction corresponding to alt

        #---southern ocean layer---#

        self.tcc_alt_so = a['cl_so'][:] # southern ocean layer total cloud fraction corresponding to alt
        self.tclw_alt_so = a['clw_so'][:] # southern ocean layer cloud liquid water fraction corresponding to liq_alt
        self.tciw_alt_so = a['ciw_so'][:] # southern ocean layer cloud ice water fraction corresponding to alt


models = {
    'CMIP5-CESM1-CAM5' : Model('CMIP5-CESM1-CAM5/07.2006_04.2011_ecmwf_era5.h5', 'CMIP5-CESM1-CAM5' ),
    'CMIP5-GFDL-HIRAM-C360' : Model('CMIP5-GFDL-HIRAM-C360/07.2006_04.2011_gfdl_am4.h5', 'CMIP5-GFDL-HIRAM-C360' ),
    'CMIP5-GISS-E2R' : Model('CMIP5-GISS-E2R/07.2006_04.2011_ecmwf_era5.h5', 'CMIP5-GISS-E2R' ),
    'CMIP5-IPSL-CM5A-LR' : Model('CMIP5-IPSL-CM5A-LR/07.2006_04.2011_gfdl_am4.h5', 'GFDL-HIRAM-C360' ),
    'CMIP5-MIROC5' : Model('CMIP5-MIROC5/07.2006_04.2011_ecmwf_era5.h5', 'CMIP5-MIROC5' ),
    'CMIP5-MRI-CGCM3' : Model('CMIP5-MRI-CGCM3/07.2006_04.2011_gfdl_am4.h5', 'CMIP5-MRI-CGCM3' ),
    
    'CMIP6-CESM2-CAM6' : Model('CMIP6-CESM2-CAM6/07.2006_04.2011_ecmwf_era5.h5', 'CMIP6-CESM2-CAM6' ),
    'CMIP6-GFDL-AM4' : Model('CMIP6-GFDL-AM4/07.2006_04.2011_gfdl_am4.h5', 'CMIP6-GFDL-AM4' ),
    'CMIP6-GISS-E21G' : Model('CMIP6-GISS-E21G/07.2006_04.2011_ecmwf_era5.h5', 'CMIP6-GISS-E21G' ),
    'CMIP6-IPSL-CM6A-LR' : Model('CMIP6-IPSL-CM6A-LR/07.2006_04.2011_gfdl_am4.h5', 'CMIP6-IPSL-CM6A-LR' ),
    'CMIP6-MIROC6' : Model('CMIP6-MIROC6/07.2006_04.2011_ecmwf_era5.h5', 'CMIP6-MIROC6' ),
    'CMIP6-MRI-ESM2' : Model('CMIP6-MRI-ESM2/07.2006_04.2011_gfdl_am4.h5', 'CMIP6-MRI-ESM2' )

}

#GFDL-HIRAM-C360 Data
h = Model('GFDL-HIRAM-C360/reduced_datasets/2001_2005_gfdl_hiram.h5', 'r')

#CCCM Data
c = Model('CCCM/reduced_datasets/07.2006_04.2011_cccm.h5', 'r')

#MRI-ESM2-AMIP Data
d = Model('MRI-ESM2-AMIP/reduced_datasets/07.2006_04.2011_mri_esm2.h5', 'r')

#MRI-CGCM3-AMIP Data
i = Model('MRI-CGCM3-AMIP/reduced_datasets/2001_2005_mri_cgcm3.h5', 'r')

#CESM1-CAM5-AMIP Data
p = Model('CESM1-CAM5-AMIP/reduced_datasets/2001_2005_cesm1_cam5.h5', 'r')

#CESM2-CAM6-AMIP Data
e = Model('CESM2-CAM6-AMIP/reduced_datasets/07.2006_04.2011_cesm2_cam6.h5', 'r')

#CAPLISO Data
f = Model('CALIPSO_GOCCP/reduced_datasets/07.2006_04.2011_calipso.h5', 'r')

#CERES Data
g = Model('CERES/reduced_datasets/07.2006_04.2011_ceres.h5', 'r')

#GISS-E2R-AMIP Data
j = Model('GISS-E2R-AMIP/reduced_datasets/2001_2005_giss_e2r.h5', 'r')

#GISS-E21G-AMIP Data
k = Model('GISS-E21G-AMIP/reduced_datasets/07.2006_04.2011_giss_e21g.h5', 'r')

#IPSL-CM5A-LR-AMIP Data
l = Model('IPSL-CM5A-LR-AMIP/reduced_datasets/2001_2005_ipsl_cm5a_lr.h5', 'r')

#IPSL-CM6A-LR-AMIP Data
m = Model('IPSL-CM6A-LR-AMIP/reduced_datasets/07.2006_04.2011_ipsl_cm6a_lr.h5', 'r')

#MIROC5-AMIP Data
n = Model('MIROC5-AMIP/reduced_datasets/2001_2005_miroc5.h5', 'r')

#MIROC6-AMIP Data
o = Model('MIROC6-AMIP/reduced_datasets/07.2006_04.2011_miroc6.h5', 'r')



############################################################################### Global Latitude Plots

#This needs to be at least as long as the number of models you will ever graph at the same time
colours = [':r', ':b', '-r' ]

def latitude_plot_models( models ):
    fig, (ax1) = plt.subplots()
    for model, col in zip( model_names, colours ):
        ax1.plot(model.tcc_lat_g_enhanced[:,0],model.tcc_lat_g_enhanced[:,1], col, label=model.name)
    plt.grid(True)
    plt.savefig("tcc_lat_g.svg", format="svg", bbox_inches='tight')
    plt.show()


latitude_plot_models( [ models['CWMWF-ERA5'], models[ 'CFDL-HIRAM-C360'] ])

latitude_plot_models( models.values() )

#---Plot Global Cloud Fraction with Latitude---#

"""
fig, (ax1) = plt.subplots()

ax1.plot(cccm_tcc_lat_g_enhanced[:,0],cccm_tcc_lat_g_enhanced[:,1], ':r', label='CCCM')
ax1.plot(calipso_tcc_lat_g[:,0],calipso_tcc_lat_g[:,1], ':b', label='CAPLISO')
ax1.plot(ceres_tcc_lat_g[:,0],ceres_tcc_lat_g[:,1], '-r', label='CERES')
ax1.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-k', label='ECMWF-ERA5')
ax1.plot(gfdl_hiram_tcc_lat_g[:,0],gfdl_hiram_tcc_lat_g[:,1], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(cam5_tcc_lat_g[:,0],cam5_tcc_lat_g[:,1], '-m', label='CMIP5-CESM1-CAM5-AMIP')



plt.grid(True)
plt.savefig("tcc_lat_g.svg", format="svg", bbox_inches='tight')
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
ax1.plot(miroc5_tcc_alt_g[:27,1],miroc5_tcc_alt_g[:27,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tcc_alt_g[:20,1],ipsl5_tcc_alt_g[:20,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tcc_alt_g[:24,1],giss5_tcc_alt_g[:24,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(cccm_tcc_alt_g[4:92,1],cccm_tcc_alt_g[4:92,0], ':r', label='CCCM')
ax2.plot(calipso_tcc_alt_g[:,1],calipso_tcc_alt_g[:,0], ':b', label='CALIPSO')
ax2.plot(ecmwf_tcc_alt_g[9:,1],ecmwf_tcc_alt_g[9:,0], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tcc_alt_g[:23,1],gfdl4_tcc_alt_g[:23,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_alt_g[:42,1],mri_tcc_alt_g[:42,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tcc_alt_g[10:,1],cam6_tcc_alt_g[10:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tcc_alt_g[:31,1],miroc6_tcc_alt_g[:31,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tcc_alt_g[:47,1],ipsl6_tcc_alt_g[:47,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tcc_alt_g[:24,1],giss6_tcc_alt_g[:24,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

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
ax1.plot(calipso_tclw_frac_alt_g[:20,1],calipso_tclw_frac_alt_g[:20,0], ':b', label='CALIPSO')
ax1.plot(ecmwf_tclw_frac_alt_g[18:,1],ecmwf_tclw_frac_alt_g[18:,0], '-k', label='ECMWF-ERA5')

ax1.plot(gfdl_hiram_tclw_frac_alt_g[:18,1],gfdl_hiram_tclw_frac_alt_g[:18,0], '-g', label='CMIP5-GFDL-HIRAM-AMIP')
ax1.plot(mri_cgcm_tclw_frac_alt_g[:19,1],mri_cgcm_tclw_frac_alt_g[:19,0], '-m', label='CMIP5-MRI_CGCM3-AMIP')
ax1.plot(miroc5_tclw_frac_alt_g[:20,1],miroc5_tclw_frac_alt_g[:20,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tclw_frac_alt_g[:15,1],ipsl5_tclw_frac_alt_g[:15,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tclw_frac_alt_g[:15,1],giss5_tclw_frac_alt_g[:15,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(cccm_tclw_frac_alt_g[4:50,1],cccm_tclw_frac_alt_g[4:50,0], ':r', label='CCCM')
ax2.plot(calipso_tclw_frac_alt_g[:20,1],calipso_tclw_frac_alt_g[:20,0], ':b', label='CALIPSO')
ax2.plot(ecmwf_tclw_frac_alt_g[18:,1],ecmwf_tclw_frac_alt_g[18:,0], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tclw_frac_alt_g[:18,1],gfdl4_tclw_frac_alt_g[:18,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_alt_g[:26,1],mri_tclw_frac_alt_g[:26,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_frac_alt_g[17:,1],cam6_tclw_frac_alt_g[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tclw_frac_alt_g[:21,1],miroc6_tclw_frac_alt_g[:21,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tclw_frac_alt_g[:40,1],ipsl6_tclw_frac_alt_g[:40,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tclw_frac_alt_g[:17,1],giss6_tclw_frac_alt_g[:17,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.3));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.325));

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
ax1.plot(miroc5_tcc_alt_so[:27,1],miroc5_tcc_alt_so[:27,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tcc_alt_so[:20,1],ipsl5_tcc_alt_so[:20,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tcc_alt_so[:24,1],giss5_tcc_alt_so[:24,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(cccm_tcc_alt_so[4:92,1],cccm_tcc_alt_so[4:92,0], ':r', label='CCCM')
ax2.plot(calipso_tcc_alt_so[:25,1],calipso_tcc_alt_so[:25,0], ':b', label='CALIPSO')
ax2.plot(ecmwf_tcc_alt_so[13:,1],ecmwf_tcc_alt_so[13:,0], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tcc_alt_so[:23,1],gfdl4_tcc_alt_so[:23,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tcc_alt_so[:42,1],mri_tcc_alt_so[:42,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tcc_alt_so[10:,1],cam6_tcc_alt_so[10:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tcc_alt_so[:31,1],miroc6_tcc_alt_so[:31,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tcc_alt_so[:47,1],ipsl6_tcc_alt_so[:47,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tcc_alt_so[:24,1],giss6_tcc_alt_so[:24,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.3));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.325));

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
ax1.plot(miroc5_tclw_frac_alt_so[:20,1],miroc5_tclw_frac_alt_so[:20,0], '-y', label='CMIP5-MIROC5-AMIP')
ax1.plot(ipsl5_tclw_frac_alt_so[:15,1],ipsl5_tclw_frac_alt_so[:15,0], '--m', label='CMIP5-IPSL-CM5A-LR-AMIP')
ax1.plot(giss5_tclw_frac_alt_so[:15,1],giss5_tclw_frac_alt_so[:15,0], '--g', label='CMIP5-NASA-GISS-E2R-AMIP')

ax2.plot(cccm_tclw_frac_alt_so[4:50,1],cccm_tclw_frac_alt_so[4:50,0], ':r', label='CCCM')
ax2.plot(calipso_tclw_frac_alt_so[:21,1],calipso_tclw_frac_alt_so[:21,0], ':b', label='CALIPSO')
ax2.plot(ecmwf_tclw_frac_alt_so[18:,1],ecmwf_tclw_frac_alt_so[18:,0], '-k', label='ECMWF-ERA5')

ax2.plot(gfdl4_tclw_frac_alt_so[:18,1],gfdl4_tclw_frac_alt_so[:18,0], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax2.plot(mri_tclw_frac_alt_so[:25,1],mri_tclw_frac_alt_so[:25,0], '-m', label='CMIP6-MRI-ESM2-AMIP')
ax2.plot(cam6_tclw_frac_alt_so[17:,1],cam6_tclw_frac_alt_so[17:,0], '-c', label='CMIP6-CESM2.1-CAM6-AMIP')
ax2.plot(miroc6_tclw_frac_alt_so[:21,1],miroc6_tclw_frac_alt_so[:21,0], '-y', label='CMIP6-MIROC6-AMIP')
ax2.plot(ipsl6_tclw_frac_alt_so[:40,1],ipsl6_tclw_frac_alt_so[:40,0], '--m', label='CMIP6-IPSL-CM6A-LR-AMIP')
ax2.plot(giss6_tclw_frac_alt_so[:17,1],giss6_tclw_frac_alt_so[:17,0], '--g', label='CMIP6-NASA-GISS-E21G-AMIP')

ax1.legend(loc='center', bbox_to_anchor=(0.3, -0.3));
ax2.legend(loc='center', bbox_to_anchor=(0.7, -0.325));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Liquid Water Fraction')
ax2.set_xlabel('Cloud Liquid Water Fraction')

#plt.title('2007 to 2008 Southern Ocean Cloud Fraction vs Altitude')

ax1.set_xlim(0, 0.6)
ax2.set_xlim(0, 0.6)

ax1.text(-0.17, 10, 'a)')
ax2.text(-0.08, 10, 'b)')

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
fig, ax = plt.subplots(nrows=7, ncols=3, figsize=(10, 20))

ax[0, 1].contourf(ecmwf_lat, ecmwf_alt[19:], ecmwf_tclw_frac_alt_lat[19:], vmin=0, vmax=0.5)
ax[0, 0].set_ylabel('Altitude (km)')
ax[0, 1].set_title('b) ECMWF-ERA5')
ecmwf_temp = ax[0, 1].contour(ecmwf_lat, ecmwf_alt[19:], (ecmwf_temp_alt_lat[19:] - 273.15), colors='grey')
ecmwf_temp.collections[5].set_linewidth(2)
ecmwf_temp.collections[5].set_color('white')
ax[0, 1].clabel(ecmwf_temp, inline=1, fontsize=10)

ax[0, 0].contourf(calipso_lat, calipso_alt[0:16], calipso_tclw_frac_alt_lat[0:16], vmin=0, vmax=0.5)
ecmwf_temp = ax[0, 0].contour(ecmwf_lat, ecmwf_alt[19:], (ecmwf_temp_alt_lat[19:] - 273.15), colors='grey')
ax[0, 0].set_title('a) CALIPSO-GOCCP')
ecmwf_temp.collections[5].set_linewidth(2)
ecmwf_temp.collections[5].set_color('white')
ax[0, 0].clabel(ecmwf_temp, inline=1, fontsize=10)

#calipso_temp = ax[0, 0].contour(calipso_lat, calipso_alt_temp[1:7], (calipso_temp_alt_lat[1:7] - 273.15), colors='grey')
#ax[0, 0].clabel(calipso_temp, inline=1, fontsize=10)

ax[1, 0].contourf(gfdl_hiram_lat, gfdl_hiram_alt[0:17], gfdl_hiram_tclw_frac_alt_lat[0:17], vmin=0, vmax=0.5)
ax[1, 0].set_title('c) CMIP5-GFDL-HIRAM')
ax[1, 0].set_ylabel('Altitude (km)')
gfdl_hiram_temp = ax[1, 0].contour(gfdl_hiram_lat, gfdl_hiram_alt_temp[1:7], (gfdl_hiram_temp_alt_lat[1:7] - 273.15), colors='grey')
gfdl_hiram_temp.collections[7].set_linewidth(2)
gfdl_hiram_temp.collections[7].set_color('white')
ax[1, 0].clabel(gfdl_hiram_temp, inline=1, fontsize=10)

ax[1, 1].contourf(gfdl4_lat, gfdl4_alt[0:18], gfdl4_tclw_frac_alt_lat[0:18], vmin=0, vmax=0.5)
ax[1, 1].set_title('d) CMIP6-GFDL-AM4')
gfdl4_temp = ax[1, 1].contour(gfdl4_lat, gfdl4_alt_temp[1:7], (gfdl4_temp_alt_lat[1:7] - 273.15), colors='grey')
gfdl4_temp.collections[6].set_linewidth(2)
gfdl4_temp.collections[6].set_color('white')
ax[1, 1].clabel(gfdl4_temp, inline=1, fontsize=10)

ax[2, 0].contourf(mri_cgcm_lat, mri_cgcm_alt[0:18], mri_cgcm_tclw_frac_alt_lat[0:18], vmin=0, vmax=0.5)
ax[2, 0].set_title('e) CMIP5-MRI-CGCM3')
ax[2, 0].set_xlabel('Latitude')
ax[2, 0].set_ylabel('Altitude (km)')
mri_cgcm_temp = ax[2, 0].contour(mri_cgcm_lat, mri_cgcm_alt_temp[1:7], (mri_cgcm_temp_alt_lat[1:7] - 273.15), colors='grey')
mri_cgcm_temp.collections[5].set_linewidth(2)
mri_cgcm_temp.collections[5].set_color('white')
ax[2, 0].clabel(mri_cgcm_temp, inline=1, fontsize=10)

cont=ax[2, 1].contourf(mri_lat, mri_alt[0:25], mri_tclw_frac_alt_lat[0:25], vmin=0, vmax=0.5)
ax[2, 1].set_title('f) CMIP6-MRI_ESM2')
mri_temp = ax[2, 1].contour(mri_lat, mri_alt_temp[1:7], (mri_temp_alt_lat[1:7] - 273.15), colors='grey')
mri_temp.collections[6].set_linewidth(2)
mri_temp.collections[6].set_color('white')
ax[2, 1].clabel(mri_temp, inline=1, fontsize=10)

ax[3, 0].contourf(cam5_lat, cam5_alt[:12], cam5_tclw_frac_alt_lat[:12], vmin=0, vmax=0.4)
ax[3, 0].set_title('g) CMIP5-CESM1-CAM5')
ax[3, 0].set_ylabel('Altitude (km)')
cam5_temp = ax[3, 0].contour(cam5_lat, cam5_alt_temp[1:7], (cam5_temp_alt_lat[1:7] - 273.15), colors='grey')
cam5_temp.collections[6].set_linewidth(2)
cam5_temp.collections[6].set_color('white')
ax[3, 0].clabel(cam5_temp, inline=1, fontsize=10)

ax[3, 1].contourf(cam6_lat, cam6_alt[19:32], cam6_tclw_frac_alt_lat[19:32], vmin=0, vmax=0.5)
ax[3, 1].set_title('h) CMIP6-CESM2.1-CAM6')
ax[3, 1].set_ylabel('Altitude (km)')
ax[3, 1].set_xlabel('Latitude')
cam6_temp = ax[3, 1].contour(cam6_lat, cam6_alt[19:32], (cam6_temp_alt_lat[19:32] - 273.15), colors='grey')
cam6_temp.collections[4].set_linewidth(2)
cam6_temp.collections[4].set_color('white')
ax[3, 1].clabel(cam6_temp, inline=1, fontsize=10)

ax[4, 0].contourf(miroc5_lat, miroc5_alt[1:19], miroc5_tclw_frac_alt_lat[1:19], vmin=0, vmax=0.4)
ax[4, 0].set_title('i) CMIP5-MIROC5')
ax[4, 0].set_ylabel('Altitude (km)')
miroc5_temp = ax[4, 0].contour(miroc5_lat, miroc5_alt_temp[1:7], (miroc5_temp_alt_lat[1:7] - 273.15), colors='grey')
miroc5_temp.collections[6].set_linewidth(2)
miroc5_temp.collections[6].set_color('white')
ax[4, 0].clabel(miroc5_temp, inline=1, fontsize=10)

ax[4, 1].contourf(miroc6_lat, miroc6_alt[1:19], miroc6_tclw_frac_alt_lat[1:19], vmin=0, vmax=0.4)
ax[4, 1].set_title('j) CMIP6-MIROC6')
miroc6_temp = ax[4, 1].contour(miroc6_lat, miroc6_alt_temp[1:7], (miroc6_temp_alt_lat[1:7] - 273.15), colors='grey')
miroc6_temp.collections[5].set_linewidth(2)
miroc6_temp.collections[5].set_color('white')
ax[4, 1].clabel(miroc6_temp, inline=1, fontsize=10)

ax[5, 0].contourf(giss5_lat, giss5_alt[:16], giss5_tclw_frac_alt_lat[:16], vmin=0, vmax=0.4)
ax[5, 0].set_title('k) CMIP5-GISS-E2R')
ax[5, 0].set_ylabel('Altitude (km)')
giss5_temp = ax[5, 0].contour(giss5_lat, giss5_alt_temp[1:7], (giss5_temp_alt_lat[1:7] - 273.15), colors='grey')
giss5_temp.collections[6].set_linewidth(2)
giss5_temp.collections[6].set_color('white')
ax[5, 0].clabel(giss5_temp, inline=1, fontsize=10)

ax[5, 1].contourf(giss6_lat, giss6_alt[:16], giss6_tclw_frac_alt_lat[:16], vmin=0, vmax=0.4)
ax[5, 1].set_title('l) CMIP6-GISS-E21G')
giss6_temp = ax[5, 1].contour(giss6_lat, giss6_alt_temp[1:7], (giss6_temp_alt_lat[1:7] - 273.15), colors='grey')
giss6_temp.collections[6].set_linewidth(2)
giss6_temp.collections[6].set_color('white')
ax[5, 1].clabel(giss6_temp, inline=1, fontsize=10)

ax[6, 0].contourf(ipsl5_lat, ipsl5_alt[:14], ipsl5_tclw_frac_alt_lat[:14], vmin=0, vmax=0.4)
ax[6, 0].set_title('m) CMIP5-IPSL-CM5A-LR')
ax[6, 0].set_ylabel('Altitude (km)')
ax[6, 0].set_xlabel('Latitude')
ipsl5_temp = ax[6, 0].contour(ipsl5_lat, ipsl5_alt_temp[1:7], (ipsl5_temp_alt_lat[1:7] - 273.15), colors='grey')
ipsl5_temp.collections[6].set_linewidth(2)
ipsl5_temp.collections[6].set_color('white')
ax[6, 0].clabel(ipsl5_temp, inline=1, fontsize=10)

ax[6, 1].contourf(ipsl6_lat, ipsl6_alt[:36], ipsl6_tclw_frac_alt_lat[:36], vmin=0, vmax=0.4)
ax[6, 1].set_title('n) CMIP6-IPSL-CM6A-LR')
ax[6, 1].set_xlabel('Latitude')
ipsl6_temp = ax[6, 1].contour(ipsl6_lat, ipsl6_alt_temp[1:7], (ipsl6_temp_alt_lat[1:7] - 273.15), colors='grey')
ipsl6_temp.collections[6].set_linewidth(2)
ipsl6_temp.collections[6].set_color('white')
ax[6, 1].clabel(ipsl6_temp, inline=1, fontsize=10)



ax[1, 2].remove()  # don't display empty ax
ax[2, 2].remove()  # don't display empty ax
ax[3, 2].remove()  # don't display empty ax
ax[4, 2].remove()  # don't display empty ax
ax[5, 2].remove()  # don't display empty ax
ax[0, 2].remove()  # don't display empty ax
ax[6, 2].remove()  # don't display empty ax


cbaxes = fig.add_axes([0.7, 0.5, 0.03, 0.3]) #(x-position, y-position, thickness, length)
cbar = fig.colorbar(cont, cax=cbaxes)
cbar.set_clim(0, 0.5)
cbar.set_label('Cloud liquid Water Fraction')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)
plt.savefig("contour_tclw.pdf", format="pdf", bbox_inches='tight')
plt.show()
"""
"""

#---Combined Grid tciw_frac---#

fig, ax = plt.subplots(nrows=7, ncols=3, figsize=(10, 20))

ax[0, 1].contourf(ecmwf_lat, ecmwf_alt[9:], ecmwf_tciw_frac_alt_lat[9:], vmin=0, vmax=0.4)
ax[0, 1].set_xlabel('Latitude')
ax[0, 0].set_ylabel('Altitude (km)')
ax[0, 1].set_title('b) ECMWF-ERA5')
ecmwf_temp = ax[0, 1].contour(ecmwf_lat, ecmwf_alt[9:], (ecmwf_temp_alt_lat[9:] - 273.15), colors='grey')
ecmwf_temp.collections[6].set_linewidth(2)
ecmwf_temp.collections[6].set_color('white')
ax[0, 1].clabel(ecmwf_temp, inline=1, fontsize=10)

ax[0, 0].contourf(calipso_lat, calipso_alt, calipso_tciw_frac_alt_lat, vmin=0, vmax=0.4)
ax[0, 0].set_title('a) CALIPSO-GOCCP')
ecmwf_temp = ax[0, 0].contour(ecmwf_lat, ecmwf_alt[19:], (ecmwf_temp_alt_lat[19:] - 273.15), colors='grey')
ecmwf_temp.collections[5].set_linewidth(2)
ecmwf_temp.collections[5].set_color('white')
ax[0, 0].clabel(ecmwf_temp, inline=1, fontsize=10)


ax[1, 0].contourf(gfdl_hiram_lat, gfdl_hiram_alt[:26], gfdl_hiram_tciw_frac_alt_lat[:26], vmin=0, vmax=0.4)
ax[1, 0].set_title('c) CMIP5-GFDL-HIRAM')
gfdl_hiram_temp = ax[1, 0].contour(gfdl_hiram_lat, gfdl_hiram_alt_temp[1:13], (gfdl_hiram_temp_alt_lat[1:13] - 273.15), colors='grey')
gfdl_hiram_temp.collections[6].set_linewidth(2)
gfdl_hiram_temp.collections[6].set_color('white')
ax[1, 0].clabel(gfdl_hiram_temp, inline=1, fontsize=10)

ax[1, 1].contourf(gfdl4_lat, gfdl4_alt[:26], gfdl4_tciw_frac_alt_lat[:26], vmin=0, vmax=0.4)
ax[1, 1].set_title('d) CMIP6-GFDL-AM4')
gfdl4_temp = ax[1, 1].contour(gfdl4_lat, gfdl4_alt_temp[1:13], (gfdl4_temp_alt_lat[1:13] - 273.15), colors='grey')
gfdl4_temp.collections[6].set_linewidth(2)
gfdl4_temp.collections[6].set_color('white')
ax[1, 1].clabel(gfdl4_temp, inline=1, fontsize=10)

ax[2, 0].contourf(mri_cgcm_lat, mri_cgcm_alt[:30], mri_cgcm_tciw_frac_alt_lat[:30], vmin=0, vmax=0.4)
ax[2, 0].set_title('e) CMIP5-MRI-CGCM3')
ax[2, 0].set_ylabel('Altitude (km)')
mri_cgcm_temp = ax[2, 0].contour(mri_cgcm_lat, mri_cgcm_alt_temp[1:13], (mri_cgcm_temp_alt_lat[1:13] - 273.15), colors='grey')
mri_cgcm_temp.collections[6].set_linewidth(2)
mri_cgcm_temp.collections[6].set_color('white')
ax[2, 0].clabel(mri_cgcm_temp, inline=1, fontsize=10)

cont=ax[2, 1].contourf(mri_lat, mri_alt[:44], mri_tciw_frac_alt_lat[:44], vmin=0, vmax=0.4)
ax[2, 1].set_title('f) CMIP6-MRI_ESM2')
mri_temp = ax[2, 1].contour(mri_lat, mri_alt_temp[1:12], (mri_temp_alt_lat[1:12] - 273.15), colors='grey')
mri_temp.collections[6].set_linewidth(2)
mri_temp.collections[6].set_color('white')
ax[2, 1].clabel(mri_temp, inline=1, fontsize=10)

ax[3, 1].contourf(cam6_lat, cam6_alt[7:32], cam6_tciw_frac_alt_lat[7:32], vmin=0, vmax=0.4)
ax[3, 1].set_title('g) CMIP6-CESM2.1-CAM6')
ax[3, 1].set_ylabel('Altitude (km)')
ax[3, 1].set_xlabel('Latitude')
cam6_temp = ax[3, 1].contour(cam6_lat, cam6_alt[7:32], (cam6_temp_alt_lat[7:32] - 273.15), colors='grey')
cam6_temp.collections[6].set_linewidth(2)
cam6_temp.collections[6].set_color('white')
ax[3, 1].clabel(cam6_temp, inline=1, fontsize=10)

ax[4, 0].contourf(miroc5_lat, miroc5_alt[1:32], miroc5_tciw_frac_alt_lat[1:32], vmin=0, vmax=0.4)
ax[4, 0].set_title('h) CMIP5-MIROC5')
ax[4, 0].set_ylabel('Altitude (km)')
miroc5_temp = ax[4, 0].contour(miroc5_lat, miroc5_alt_temp[1:13], (miroc5_temp_alt_lat[1:13] - 273.15), colors='grey')
miroc5_temp.collections[6].set_linewidth(2)
miroc5_temp.collections[6].set_color('white')
ax[4, 0].clabel(miroc5_temp, inline=1, fontsize=10)

ax[4, 1].contourf(miroc6_lat, miroc6_alt[1:37], miroc6_tciw_frac_alt_lat[1:37], vmin=0, vmax=0.4)
ax[4, 1].set_title('i) CMIP6-MIROC6')
miroc6_temp = ax[4, 1].contour(miroc6_lat, miroc6_alt_temp[1:13], (miroc6_temp_alt_lat[1:13] - 273.15), colors='grey')
miroc6_temp.collections[5].set_linewidth(2)
miroc6_temp.collections[5].set_color('white')
ax[4, 1].clabel(miroc6_temp, inline=1, fontsize=10)

ax[5, 0].contourf(giss5_lat, giss5_alt[:28], giss5_tciw_frac_alt_lat[:28], vmin=0, vmax=0.4)
ax[5, 0].set_title('j) CMIP5-GISS-E2R')
ax[5, 0].set_ylabel('Altitude (km)')
giss5_temp = ax[5, 0].contour(giss5_lat, giss5_alt_temp[1:13], (giss5_temp_alt_lat[1:13] - 273.15), colors='grey')
giss5_temp.collections[6].set_linewidth(2)
giss5_temp.collections[6].set_color('white')
ax[5, 0].clabel(giss5_temp, inline=1, fontsize=10)

ax[5, 1].contourf(giss6_lat, giss6_alt[:28], giss6_tciw_frac_alt_lat[:28], vmin=0, vmax=0.4)
ax[5, 1].set_title('k) CMIP6-GISS-E21G')
giss6_temp = ax[5, 1].contour(giss6_lat, giss6_alt_temp[1:13], (giss6_temp_alt_lat[1:13] - 273.15), colors='grey')
giss6_temp.collections[6].set_linewidth(2)
giss6_temp.collections[6].set_color('white')
ax[5, 1].clabel(giss6_temp, inline=1, fontsize=10)

ax[6, 0].contourf(ipsl5_lat, ipsl5_alt[:24], ipsl5_tciw_frac_alt_lat[:24], vmin=0, vmax=0.4)
ax[6, 0].set_title('l) CMIP5-IPSL-CM5A-LR')
ax[6, 0].set_ylabel('Altitude (km)')
ax[6, 0].set_xlabel('Latitude')
ipsl5_temp = ax[6, 0].contour(ipsl5_lat, ipsl5_alt_temp[1:13], (ipsl5_temp_alt_lat[1:13] - 273.15), colors='grey')
ipsl5_temp.collections[6].set_linewidth(2)
ipsl5_temp.collections[6].set_color('white')
ax[6, 0].clabel(ipsl5_temp, inline=1, fontsize=10)

ax[6, 1].contourf(ipsl6_lat, ipsl6_alt[:51], ipsl6_tciw_frac_alt_lat[:51], vmin=0, vmax=0.4)
ax[6, 1].set_title('m) CMIP6-IPSL-CM6A-LR')
ax[6, 1].set_xlabel('Latitude')
ipsl6_temp = ax[6, 1].contour(ipsl6_lat, ipsl6_alt_temp[1:13], (ipsl6_temp_alt_lat[1:13] - 273.15), colors='grey')
ipsl6_temp.collections[6].set_linewidth(2)
ipsl6_temp.collections[6].set_color('white')
ax[6, 1].clabel(ipsl6_temp, inline=1, fontsize=10)



ax[1, 2].remove()  # don't display empty ax
ax[2, 2].remove()  # don't display empty ax
ax[3, 2].remove()  # don't display empty ax
ax[4, 2].remove()  # don't display empty ax
ax[5, 2].remove()  # don't display empty ax
ax[3, 0].remove()  # don't display empty ax

ax[6, 2].remove()  # don't display empty ax

ax[0, 2].remove()  # don't display empty ax

cbaxes = fig.add_axes([0.7, 0.5, 0.03, 0.3]) 
cbar = fig.colorbar(cont, cax=cbaxes)
cbar.set_clim(0, 0.4)
cbar.set_label('Cloud Ice Water Fraction')
plt.tight_layout(pad=0.4, w_pad=0.6, h_pad=3)
plt.savefig("2007_2008_contour_tciw.svg", format="svg", bbox_inches='tight')
plt.show()

"""










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

