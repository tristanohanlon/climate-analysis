# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon - University of Auckland, Samual Crookes & Jonathan Rogers


"""
from pprint import pprint
import time
import numpy as np
import os
from pyhdf import SD
import h5py
import matplotlib.pyplot as plt
from scipy import interpolate
import constants
from scipy import stats
import cartopy.crs as ccrs

###############################################################################
start = time.time()

#set location

location = constants.home
#lat = constants.lat

altitude_types = [ 'Irradiance layer center height profile', 'Layer center height profile (clouds and aerosol)', 'Irradiance level height profile' ]

class DataSet:
    def __init__( self, type_name, altitude_type ):
        self.type_name = type_name
        self.altitude_type = altitude_type
        self.data = np.zeros(( constants.lat.size, constants.alt.size ))
        self.data_counts = np.zeros( (constants.lat.size, constants.alt.size ))
        

# class grid_DataSet:
#     def __init__( self, type_name ):
#         self.type_name = type_name
#         self.data = np.zeros(( constants.lat.size, constants.lon.size ))
#         self.data_counts = np.zeros( (constants.lat.size, constants.lon.size ))


data_sets = [
        DataSet('Cloud fraction profile', 1 ),
        DataSet('Liquid water content profile used', 0 ), # kg/m3
        DataSet('Ice water content profile used', 0 ), # kg/m3
        DataSet('Temperature profile', 2 ),
        DataSet('Pressure profile', 2 )      
]        

# grid_data_sets = [
#         grid_DataSet('Cloud free area percent coverage (CALIPSO-CloudSat)' ),
#         grid_DataSet('Liquid water content profile used' ),
#         grid_DataSet('Ice water content profile used' ),
# ]        

# The directory where your HDF files are stored
os.chdir('E:/CCCM/test')  # Home PC


# Load every file in the directory
for filename in os.listdir():
    if os.path.isdir( filename ):
        continue
    pprint( filename )
    # Load the file
    f = SD.SD(filename)
    raw_lon = f.select('Longitude of subsatellite point at surface at observation').get()
    raw_lat = 90 - f.select('Colatitude of subsatellite point at surface at observation').get()
    altitudes = []
    for altitude_type in altitude_types:
        altitudes.append( f.select(altitude_type).get() )
    altitudes[0] /= 1000
    altitudes[2] /= 1000
    

    for data_set in data_sets:
        sel = f.select( data_set.type_name )
        data = sel.get()

        fill_value = sel.attributes()['_FillValue']
        for l_index, l in enumerate( raw_lat ):
            if l <= constants.min_lat or l >= constants.max_lat:
                continue
            lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
            for a_index, a in enumerate( altitudes[data_set.altitude_type] ):
                if a <= constants.min_alt or a >= constants.max_alt:
                    continue
                alt_bin = int( ( a - constants.min_alt ) / constants.alt_division )
                #print( lat_bin, alt_bin )
                val = data[l_index, a_index]
                if val == fill_value:
                    continue
                data_set.data[ lat_bin, alt_bin ] += val
                data_set.data_counts[ lat_bin, alt_bin ] += 1

#     for data_set in grid_data_sets:
#         sel = f.select( data_set.type_name )
#         data = sel.get()
#         fill_value = sel.attributes()['_FillValue']
#         data[ data == fill_value ] = None #set fill values to nan
#         if data_set.type_name == 'Cloud free area percent coverage (CALIPSO-CloudSat)':
#             data = (100 - data) / 100
#         else:
#             data = np.nansum(data, axis = 1)
#         for l_index, ( la, lo ) in enumerate( zip( raw_lat, raw_lon ) ):
#             if la <= constants.min_lat or la >= constants.max_lat:
#                 continue
            
#             lat_bin = int( ( la - constants.min_lat ) / constants.lat_division)
#             lon_bin = int( lo - 1 )
            
#             #print( lat_bin, alt_bin )
#             val = data[l_index]
#             if val == None:
#                 continue
#             data_set.data[ lat_bin, lon_bin ] += val
#             data_set.data_counts[ lat_bin, lon_bin ] += 1
            
# for data_set in grid_data_sets:
#     data_set.data /= data_set.data_counts
# #    pprint( data_set.data.shape )
  
# clt_lat_lon = grid_data_sets[0].data
# clwvi_lat_lon = grid_data_sets[1].data
# clivi_lat_lon = grid_data_sets[2].data

for data_set in data_sets:
    data_set.data /= data_set.data_counts
#    pprint( data_set.data.shape )
  
cl_alt_lat = data_sets[0].data # 0-1
full_clwc_alt_lat = data_sets[1].data # g/m^3
clic_alt_lat = data_sets[2].data # g/m^3
full_ta_alt_lat = data_sets[3].data # K
full_p_alt_lat = data_sets[4].data # hPa


# convert to cloud water mass fractions
#https://www.translatorscafe.com/unit-converter/en-US/calculator/altitude/

density_alt_lat = ( ( full_p_alt_lat * 100 ) / ( 287.052 * full_ta_alt_lat ) ) * 1000 # g/m^3
full_clw_alt_lat = full_clwc_alt_lat / density_alt_lat
cli_alt_lat = clic_alt_lat / density_alt_lat

#---create reduced altitude liquid mass fraction and temperatures---#

ta_alt_lat = full_ta_alt_lat[:,:constants.liq_alt_confine]
clw_alt_lat = full_clw_alt_lat[:,:constants.liq_alt_confine]
clwc_alt_lat = full_clwc_alt_lat[:,:constants.liq_alt_confine]

#---create fractions---#

full_clw_frac_alt_lat = ( full_clw_alt_lat / ( full_clw_alt_lat + cli_alt_lat ) ) * cl_alt_lat
cli_frac_alt_lat = ( cli_alt_lat / ( full_clw_alt_lat + cli_alt_lat ) ) * cl_alt_lat
clw_frac_alt_lat = full_clw_frac_alt_lat[:,:constants.liq_alt_confine]


#---create liquid and ice fractions---#

# clt = np.nanmean( clt_lat_lon, axis = -1 )
# clwvi = np.nanmean( clwvi_lat_lon, axis = -1 )
# clivi = np.nanmean( clivi_lat_lon, axis = -1 )


#---create southern ocean and global datasets---#
    
cl_so = constants.globalalt_latMean(np.transpose(cl_alt_lat[constants.so_idx_1:constants.so_idx_2]), constants.lat[constants.so_idx_1:constants.so_idx_2])
full_clw_so = constants.globalalt_latMean(np.transpose(full_clw_alt_lat[constants.so_idx_1:constants.so_idx_2]), constants.lat[constants.so_idx_1:constants.so_idx_2])
cli_so = constants.globalalt_latMean(np.transpose(cli_alt_lat[constants.so_idx_1:constants.so_idx_2]), constants.lat[constants.so_idx_1:constants.so_idx_2])

clwc_so = constants.globalalt_latMean(np.transpose(clwc_alt_lat[constants.so_idx_1:constants.so_idx_2]), constants.lat[constants.so_idx_1:constants.so_idx_2])
clic_so = constants.globalalt_latMean(np.transpose(clic_alt_lat[constants.so_idx_1:constants.so_idx_2]), constants.lat[constants.so_idx_1:constants.so_idx_2])

full_clw_frac_so = ( full_clw_so / ( full_clw_so + cli_so ) ) * cl_so
cli_frac_so = ( cli_so / ( full_clw_so + cli_so ) ) * cl_so
clw_frac_so = full_clw_frac_so[:constants.liq_alt_confine]
clw_so = full_clw_so[:constants.liq_alt_confine]

cl_g = constants.globalalt_latMean(np.transpose(cl_alt_lat[8:172]), constants.lat[8:172])
full_clw_g = constants.globalalt_latMean(np.transpose(full_clw_alt_lat[8:172]), constants.lat[8:172])
cli_g = constants.globalalt_latMean(np.transpose(cli_alt_lat[8:172]), constants.lat[8:172])

clwc_g = constants.globalalt_latMean(np.transpose(clwc_alt_lat[8:172]), constants.lat[8:172])
clic_g = constants.globalalt_latMean(np.transpose(clic_alt_lat[8:172]), constants.lat[8:172])

full_clw_frac_g = ( full_clw_g / ( full_clw_g + cli_g ) ) * cl_g
cli_frac_g = ( cli_g / ( full_clw_g + cli_g ) ) * cl_g
clw_frac_g = full_clw_frac_g[:constants.liq_alt_confine]
clw_g = full_clw_g[:constants.liq_alt_confine]

#----------------------------#

    ######## Binned Temperature Data ########
# values to bin: full_clw_alt_lat and full_ta_alt_lat
# binned into constants.ta_g array
# values in each bin to be summed
# call the summed values clw_t_g

stat = 'mean'
clwc_t_g, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat.flatten(), full_clwc_alt_lat.flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))
clwc_t_so, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat[constants.so_idx_1:constants.so_idx_2].flatten(), full_clwc_alt_lat[constants.so_idx_1:constants.so_idx_2].flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))

clw_t_g, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat.flatten(), full_clw_alt_lat.flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))
clw_t_so, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat[constants.so_idx_1:constants.so_idx_2].flatten(), full_clw_alt_lat[constants.so_idx_1:constants.so_idx_2].flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))

cli_t_g, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat.flatten(), cli_alt_lat.flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))
cli_t_so, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat[constants.so_idx_1:constants.so_idx_2].flatten(), cli_alt_lat[constants.so_idx_1:constants.so_idx_2].flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))

cl_t_g, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat.flatten(), cl_alt_lat.flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))
cl_t_so, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat[constants.so_idx_1:constants.so_idx_2].flatten(), cl_alt_lat[constants.so_idx_1:constants.so_idx_2].flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))

clw_frac_t_g = ( clw_t_g / ( clw_t_g + cli_t_g ) ) * cl_t_g
clw_frac_t_so = ( cli_t_g / ( clw_t_g + cli_t_g ) ) * cl_t_g

# fig, ax = plt.subplots()
# ax.plot( constants.ta, clw_t_g )
# ax.plot( constants.ta, clw_t_so )
# ax.axvline(x=273, label = '273K', color = 'black', linestyle='--')
# plt.grid(True)
# plt.show()
######################

#----------------------------#

os.chdir( location + '/climate-analysis/reduced_data' )

save_filename = 'Jun_2006_Apr_2011_CCCM.h5'

with h5py.File(save_filename, 'w') as p:

    # p.create_dataset('clt', data=clt)
    # p.create_dataset('clwvi', data=clwvi)
    # p.create_dataset('clivi', data=clivi)
    # p.create_dataset('clt_lat_lon', data=clt_lat_lon)
    # p.create_dataset('clwvi_lat_lon', data=clwvi_lat_lon)
    # p.create_dataset('clivi_lat_lon', data=clivi_lat_lon)

    p.create_dataset('density_alt_lat', data=np.transpose( density_alt_lat ) ) # global layer total cloud fraction corresponding to alt
    p.create_dataset('full_p_alt_lat', data=np.transpose( full_p_alt_lat ) ) # global layer total cloud fraction corresponding to alt

    p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
    p.create_dataset('clw_g', data=clw_g) # global layer cloud liquid water mass fraction in air ( kg/kg ) corresponding to liq_alt
    p.create_dataset('cli_g', data=cli_g) # global layer cloud ice water mass fraction in air ( kg/kg ) corresponding to alt
    p.create_dataset('clwc_g', data=clwc_g) # global layer cloud liquid water content ( g/m^3 ) corresponding to liq_alt
    p.create_dataset('clic_g', data=clic_g) # global layer cloud ice water content ( g/m^3 ) corresponding to alt
    p.create_dataset('clw_frac_g', data=clw_frac_g) # global layer cloud liquid water fraction corresponding to liq_alt
    p.create_dataset('cli_frac_g', data=cli_frac_g) # global layer cloud ice water fraction corresponding to alt
    
    p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
    p.create_dataset('clw_so', data=clw_so) # southern ocean layer cloud liquid water mass fraction in air ( kg/kg ) corresponding to liq_alt
    p.create_dataset('cli_so', data=cli_so) # southern ocean layer cloud ice water mass fraction in air ( kg/kg ) corresponding to alt
    p.create_dataset('clwc_so', data=clwc_so) # southern ocean layer cloud liquid water content ( g/m^3 ) corresponding to liq_alt
    p.create_dataset('clic_so', data=clic_so) # southern ocean layer cloud ice water content ( g/m^3 ) corresponding to alt
    p.create_dataset('clw_frac_so', data=clw_frac_so) # southern ocean layer cloud liquid water fraction corresponding to liq_alt
    p.create_dataset('cli_frac_so', data=cli_frac_so) # southern ocean layer cloud ice water fraction corresponding to alt

    p.create_dataset('cl_t_g', data=cl_t_g) # global layer cloud liquid water fraction corresponding to ta
    p.create_dataset('cl_t_so', data=cl_t_so) # global layer cloud liquid water fraction corresponding to ta

    p.create_dataset('clwc_t_g', data=clwc_t_g) # global layer cloud liquid water content ( g/m^3 ) corresponding to ta
    p.create_dataset('clwc_t_so', data=clwc_t_so) # global layer cloud liquid water content ( g/m^3 ) corresponding to ta

    p.create_dataset('clw_t_g', data=clw_t_g) # global layer cloud liquid water mass fraction in air ( kg/kg ) corresponding to ta
    p.create_dataset('clw_t_so', data=clw_t_so) # global layer cloud liquid water mass fraction in air ( kg/kg ) corresponding to ta

    p.create_dataset('clw_frac_t_g', data=clw_frac_t_g) # global layer cloud liquid water fraction corresponding to ta
    p.create_dataset('clw_frac_t_so', data=clw_frac_t_so) # global layer cloud liquid water fraction corresponding to ta

    p.create_dataset('full_clwc_alt_lat', data= np.transpose( full_clwc_alt_lat ) ) # cloud liquid water content ( g/m^3 ) corresponding to alt and lat
    p.create_dataset('clic_alt_lat', data= np.transpose( clic_alt_lat ) ) # cloud ice water content ( g/m^3 ) corresponding to alt and lat

    p.create_dataset('full_clw_frac_alt_lat', data= np.transpose( full_clw_frac_alt_lat ) ) # cloud liquid water fraction corresponding to alt and lat
    p.create_dataset('clw_frac_alt_lat', data= np.transpose( clw_frac_alt_lat ) ) # cloud liquid water fraction corresponding to liq_alt and lat
    p.create_dataset('cli_frac_alt_lat', data= np.transpose( cli_frac_alt_lat ) ) # cloud ice water fraction corresponding to alt and lat

    p.create_dataset('ta_alt_lat', data= np.transpose( ta_alt_lat ) ) # temperature corresponding to liq_alt and lat
    p.create_dataset('cl_alt_lat', data= np.transpose( cl_alt_lat ) ) # total cloud fraction corresponding to alt and lat
    p.create_dataset('clw_alt_lat', data= np.transpose( clw_alt_lat ) ) # cloud liquid water mass fraction in air ( kg/kg ) corresponding to liq_alt and lat
    p.create_dataset('cli_alt_lat', data= np.transpose( cli_alt_lat ) ) # cloud ice water mass fraction in air ( kg/kg ) corresponding to alt and lat

    p.create_dataset('full_ta_alt_lat', data= np.transpose( full_ta_alt_lat ) ) # temperature corresponding to alt and lat
    p.create_dataset('full_clw_alt_lat', data= np.transpose( full_clw_alt_lat ) ) # tcloud liquid water mass fraction in air ( kg/kg ) corresponding to alt and lat

    p.close()

end = time.time()
print('Create dataset took:', end - start, 's')



































