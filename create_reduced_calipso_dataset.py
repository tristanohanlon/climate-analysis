# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon - University of Auckland & Jonathan Rogers

These functions will create a standard model dataset when parsed:
    
reduce_dataset( directory, filename, save_location, start, end, cl_is_fractional=False )

example:
reduce_dataset( 'gfdl_am4', '_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'E:/University/University/MSc/Models/climate-analysis/GFDL-AM4', datetime.datetime( 2000, 1, 1 ), datetime.datetime( 2006, 1, 1 ) )

use ECMWF ta_alt_lat over contour data
Specify time range at the bottom
"""
import numpy as np
import os
from netCDF4 import Dataset
from netCDF4 import date2index
import h5py
import math
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
import datetime
from pprint import pprint
from sklearn.impute import SimpleImputer
import constants
from scipy import ndimage as nd

def fill(data, invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') 
    by the value of the nearest valid data cell

    Input:
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. 
                 data value are replaced where invalid is True
                 If None (default), use: invalid  = np.isnan(data)

    Output: 
        Return a filled array. 
    """    
    if invalid is None: invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid, 
                                    return_distances=False, 
                                    return_indices=True)
    return data[tuple(ind)]

location = constants.home # home, uni, hdd or laptop
os.chdir( location + 'Data/CALIPSO-GOCCP' )

f = Dataset('Map_OPAQ330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
raw_lat = np.array(f.variables['latitude'][4:85])
raw_lon = np.array(f.variables['longitude'][:]) + 180

#opaque (thick) clouds
cc = np.array(f.variables['cltcalipso_opaque'][:60]) # 7.2006 to 12.2020 - 54 months
cc[cc < 0] = None #set fill values to nan
cc = np.nanmean(cc, axis = 0) #average over time

#clt_lat_lon = np.nanmean(cc, axis = 0) #average over time

#thin clouds
tf = np.array(f.variables['cltcalipso_thin'][:60]) # 7.2006 to 12.2020 - 54 months
tf[tf < 0] = None #set fill values to nan
tf = np.nanmean(tf, axis = 0) #average over time

#total clouds
clt_lat_lon = (cc + tf)[4:85] #(lat, lon)

clt = np.nanmean(clt_lat_lon, axis = -1) #average over longitude (lat)

###########################---get lat-lon - phase fraction---###########################

f = Dataset('MapLowMidHigh_Phase330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
clivi = np.array(f.variables['cltcalipso_RPIC'][:60]) # 7.2006 to 12.2020 - 54 months

clivi[clivi < 0] = None #set fill values to nan

clivi_lat_lon = (np.nanmean(clivi, axis = 0))[4:85] #average over time
clwvi_lat_lon = 1 - clivi_lat_lon

clivi_lat_lon = clivi_lat_lon * clt_lat_lon
clwvi_lat_lon = clwvi_lat_lon * clt_lat_lon

clivi = np.nanmean(clivi_lat_lon, axis = -1) #average over longitude
clwvi = np.nanmean(clwvi_lat_lon, axis = -1) #average over longitude



###########################---get alt - cloud fraction---###########################

f = Dataset('3D_CloudFraction330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
raw_alt = np.array(f.variables['alt_mid'][:]) # km

cl = np.array(f.variables['clcalipso'][:60]) # 7.2006 to 12.2020 - 54 months

cl[cl < 0] = None #set fill values to nan

cl = np.nanmean(cl, axis = 0) #average over time
cl_alt_lat = np.nanmean(cl, axis = -1)[:,4:85] #average over longitude
cl_g = np.nanmean(cl_alt_lat, axis = -1) #average over latitude


###########################---get alt - phase fractions---###########################

f = Dataset('3D_CloudFraction_Phase330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')

#ice water fraction
ice = np.array(f.variables['clcalipso_RPIC'][0:60]) # 2.2011 to 4.2011 - 3 months
ice[ice < 0] = None #set fill values to nan
ice = np.nanmean(ice, axis = 0) #average over time

#ice water fraction with alt and lat

cli_alt_lat = (np.nanmean(ice, axis = -1))[:,4:85] #average over lon
clw_alt_lat = (1 - cli_alt_lat) * cl_alt_lat
cli_alt_lat = cli_alt_lat * cl_alt_lat 
#ice water fraction
cli_g = np.nanmean(cli_alt_lat, axis = -1) #average over lat
clw_g = np.nanmean(clw_alt_lat, axis = -1) #average over lat


###############---create southern ocean data---###############


cl_so = np.hstack((np.vstack(raw_lat), np.transpose(cl_alt_lat) )) #creates a (180,34) array
cl_so = cl_so[cl_so[:,0]>=-70]
cl_so = cl_so[cl_so[:,0]<=-50]
cl_so = cl_so[:,1:] #Split the combined array into just the tccl data, eliminating the first coloumn of latitude
cl_so = np.nanmean(cl_so, axis = 0)


clw_so = np.hstack((np.vstack(raw_lat), np.transpose(clw_alt_lat))) #creates a (180,34) array
clw_so = clw_so[clw_so[:,0]>=-70]
clw_so = clw_so[clw_so[:,0]<=-50]
clw_so = clw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
clw_so = np.nanmean(clw_so, axis = 0)


cli_so = np.hstack((np.vstack(raw_lat), np.transpose(cli_alt_lat))) #creates a (180,34) array
cli_so = cli_so[cli_so[:,0]>=-70]
cli_so = cli_so[cli_so[:,0]<=-50]
cli_so = cli_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
cli_so = np.nanmean(cli_so, axis = 0)


###########################---get alt - temp - phase fractions---###########################

f = Dataset('3D_CloudFraction_Temp330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
raw_ta = np.array(f.variables['temp_mid'][:]) + 273 # km (38)

cl_t = np.array(f.variables['cltemp'][0:60]) # 7.2006 to 12.2020 - 54 months
cl_t [cl_t  < 0] = None #set fill values to nan
cl_t  = np.nanmean(cl_t, axis = 0) #average over time
cl_t_lat =  np.nanmean(cl_t, axis = -1)[:,4:85] #average over longitude

cli_t = np.array(f.variables['cltemp_phase'][0:60]) # 7.2006 to 12.2020 - 54 months
cli_t[cli_t < 0] = None #set fill values to nan

cli_t = np.nanmean(cli_t, axis = 0) #average over time
cli_t_lat =  np.nanmean(cli_t, axis = -1)[:,4:85] #average over longitude

clw_t_lat =  (1 - cli_t_lat) * cl_t_lat   
clw_t_g = np.nanmean(clw_t_lat, axis = -1)

clw_t_so = np.hstack((np.vstack(raw_lat), np.transpose(clw_t_lat))) #creates a (180,34) array
clw_t_so = clw_t_so[clw_t_so[:,0]>=-70]
clw_t_so = clw_t_so[clw_t_so[:,0]<=-50]
clw_t_so = clw_t_so[:,1:] #Split the combined array into just the tcclw data, eliminating the first coloumn of latitude
clw_t_so = np.nanmean(clw_t_so, axis = 0)

full_clw_alt_lat = clw_alt_lat 

######################################

interpolated = interpolate.interp1d(raw_lat, clt, kind = 'cubic')
clt = interpolated(constants.lat)

interpolated = interpolate.interp1d(raw_lat, clwvi, kind = 'cubic')
clwvi = interpolated(constants.lat)

interpolated = interpolate.interp1d(raw_lat, clivi, kind = 'cubic')
clivi = interpolated(constants.lat)

interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clt_lat_lon ), kind = 'cubic')
clt_lat_lon = interpolated(constants.lat, constants.lon)
clt_lat_lon = np.transpose(clt_lat_lon)

interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clwvi_lat_lon ), kind = 'cubic')
clwvi_lat_lon = interpolated(constants.lat, constants.lon)
clwvi_lat_lon = np.transpose(clwvi_lat_lon)

interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clivi_lat_lon ), kind = 'cubic')
clivi_lat_lon = interpolated(constants.lat, constants.lon)
clivi_lat_lon = np.transpose(clivi_lat_lon)


interpolated = interpolate.interp1d(raw_alt, cl_g, kind = 'cubic', fill_value="extrapolate")
cl_g = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, cl_so, kind = 'cubic', fill_value="extrapolate")
cl_so = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, clw_g, kind = 'cubic', fill_value="extrapolate")
clw_g = interpolated(constants.liq_alt)

interpolated = interpolate.interp1d(raw_alt, cli_g, kind = 'cubic', fill_value="extrapolate")
cli_g = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, clw_so, kind = 'cubic', fill_value="extrapolate")
clw_so = interpolated(constants.liq_alt)

interpolated = interpolate.interp1d(raw_alt, cli_so, kind = 'cubic', fill_value="extrapolate")
cli_so = interpolated(constants.alt)


interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( cl_alt_lat ), kind = 'cubic')
cl_alt_lat = interpolated(constants.alt, constants.lat)

interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( clw_alt_lat ), kind = 'cubic')
full_clw_alt_lat = interpolated(constants.alt, constants.lat)

interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( clw_alt_lat ), kind = 'cubic')
clw_alt_lat = interpolated(constants.liq_alt, constants.lat)

interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( cli_alt_lat ), kind = 'cubic')
cli_alt_lat = interpolated(constants.alt, constants.lat)
cli_alt_lat [cli_alt_lat  < 0] = None #set fill values to nan



interpolated = interpolate.interp1d(raw_ta, clw_t_g, kind = 'cubic', fill_value="extrapolate")
clw_t_g = interpolated(constants.ta_g)

interpolated = interpolate.interp1d(raw_ta, clw_t_so, kind = 'cubic', fill_value="extrapolate")
clw_t_so = interpolated(constants.ta_so)


    #----------------------------#

os.chdir(location + 'climate-analysis/reduced_data')

with h5py.File( 'Jun_2006_Jun_2011_CALIPSO.h5', 'w' ) as p:
    
    p.create_dataset('clt', data=clt) # total cloud fraction corresponding to lat
    p.create_dataset('clt_lat_lon', data=clt_lat_lon ) # total cloud fraction corresponding to lat, lon
  
    p.create_dataset('clwvi', data=clwvi) # total cloud liquid water fraction corresponding to lat
    p.create_dataset('clwvi_lat_lon', data=clwvi_lat_lon ) # total cloud fraction corresponding to lat, lon

    p.create_dataset('clivi', data=clivi) # total cloud ice fraction corresponding to lat
    p.create_dataset('clivi_lat_lon', data=clivi_lat_lon ) # total cloud fraction corresponding to lat, lon
    
    p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
    p.create_dataset('clw_g', data=clw_g) # global layer cloud liquid water fraction corresponding to liq_alt
    p.create_dataset('cli_g', data=cli_g) # global layer cloud ice water fraction corresponding to alt
    
    p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
    p.create_dataset('clw_so', data=clw_so) # southern ocean layer cloud liquid water fraction corresponding to liq_alt
    p.create_dataset('cli_so', data=cli_so) # southern ocean layer cloud ice water fraction corresponding to alt
 
    p.create_dataset('clw_t_g', data=clw_t_g) # global layer cloud liquid water fraction corresponding to ta_g
    p.create_dataset('clw_t_so', data=clw_t_so) # global layer cloud liquid water fraction corresponding to ta_so
   
    p.create_dataset('cl_alt_lat', data=np.transpose(cl_alt_lat) ) # total cloud fraction corresponding to alt and lat
    p.create_dataset('cli_alt_lat', data=np.transpose(cli_alt_lat) ) # cloud ice water fraction corresponding to alt and lat
    p.create_dataset('clw_alt_lat', data=np.transpose(clw_alt_lat) )  # cloud liquid water fraction corresponding to liq_alt and lat
    p.create_dataset('full_clw_alt_lat', data=np.transpose(full_clw_alt_lat) )  # cloud liquid water fraction corresponding to liq_alt and lat


    p.close()


    # take from ECMWF - import first before saving
#f = h5py.File('Jun_2006_Jun_2011_CALIPSO.h5','a')
#f.create_dataset('full_ta_alt_lat', data=full_ta_alt_lat)
#f.create_dataset('ta_alt_lat', data=ta_alt_lat)
#f.close()

