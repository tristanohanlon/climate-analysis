# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 14:44:25 2019

@author: Tristan O'Hanlon - University of Auckland

"""

#read hdf5
import h5py
import os
import matplotlib.pyplot as plt
import constants
import numpy as np
from scipy import interpolate
from scipy import stats


#--- Set Location, date period and model ---#

# specify model from the list above
model = 'CCCM' 

# specify location: home, uni, hdd, laptop
location = constants.home + 'climate-analysis/reduced_data/backup'

os.chdir( location )

min_alt = 0.25
max_alt = 20
alt_division = 0.25
raw_alt = np.arange(min_alt, max_alt, alt_division)


h5f = h5py.File( constants.date_cccm + '_' + model + '.h5', 'r')
cl_alt_lat = h5f['cl_alt_lat'][:]
cl_g = h5f['cl_g'][:]
cl_so = h5f['cl_so'][:]
cl_t_g = h5f['cl_t_g'][:]
cl_t_so = h5f['cl_t_so'][:]

cli_alt_lat = h5f['cli_alt_lat'][:]
cli_frac_alt_lat = h5f['cli_frac_alt_lat'][:]
cli_frac_g = h5f['cli_frac_g'][:]
cli_frac_so = h5f['cli_frac_so'][:]
cli_frac_t_g = h5f['cli_frac_t_g'][:]
cli_frac_t_so = h5f['cli_frac_t_so'][:]
cli_g = h5f['cli_g'][:]
cli_so = h5f['cli_so'][:]
clic_alt_lat = h5f['clic_alt_lat'][:]
clic_g = h5f['clic_g'][:]
clic_so = h5f['clic_so'][:]

clw_alt_lat = h5f['clw_alt_lat'][:]
clw_frac_alt_lat = h5f['clw_frac_alt_lat'][:]
clw_frac_g = h5f['clw_frac_g'][:]
clw_frac_so = h5f['clw_frac_so'][:]
clw_frac_t_g = h5f['clw_frac_t_g'][:]
clw_frac_t_so = h5f['clw_frac_t_so'][:]
clw_g = h5f['clw_g'][:]
clw_so = h5f['clw_so'][:]
clw_t_g = h5f['clw_t_g'][:]
clw_t_so = h5f['clw_t_so'][:]

clwc_g = h5f['clwc_g'][:]
clwc_so = h5f['clwc_so'][:]
clwc_t_g = h5f['clwc_t_g'][:]
clwc_t_so = h5f['clwc_t_so'][:]

full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
full_clw_frac_alt_lat = h5f['full_clw_frac_alt_lat'][:]
full_clwc_alt_lat = h5f['full_clwc_alt_lat'][:]

full_p_alt_lat = h5f['full_p_alt_lat'][:]
full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
ta_alt_lat = h5f['ta_alt_lat'][:]

density_alt_lat = h5f['density_alt_lat'][:]
# clt = h5f['clt'][:]
# clt_lat_lon = h5f['clt_lat_lon'][:]
# clwvi = h5f['clwvi'][:]
# clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
# clivi = h5f['clivi'][:]
# clivi_lat_lon = h5f['clivi_lat_lon'][:]
   
# print(constants.globalalt_latMeanVal(cl_alt_lat[:,constants.so_idx_1:constants.so_idx_2], constants.lat[constants.so_idx_1:constants.so_idx_2]))


interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt, cl_alt_lat[:,7:171] , kind = 'cubic')
cl_alt_lat = interpolated( constants.lat, constants.alt )
cl_alt_lat[:,0:7] = np.nan
cl_alt_lat[:,172:] = np.nan
cl_alt_lat[ cl_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp1d(raw_alt, cl_g, kind = 'cubic', fill_value = 'extrapolate')
cl_g = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, cl_so, kind = 'cubic', fill_value = 'extrapolate')
cl_so = interpolated(constants.alt)

interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt, cli_alt_lat[:,7:171] , kind = 'cubic')
cli_alt_lat = interpolated( constants.lat, constants.alt )
cli_alt_lat[:,0:7] = np.nan
cli_alt_lat[:,172:] = np.nan
cli_alt_lat[ cli_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt, clic_alt_lat[:,7:171] , kind = 'cubic')
clic_alt_lat = interpolated( constants.lat, constants.alt )
clic_alt_lat[:,0:7] = np.nan
clic_alt_lat[:,172:] = np.nan
clic_alt_lat[ clic_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp1d(raw_alt, cli_g, kind = 'cubic', fill_value = 'extrapolate')
cli_g = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, cli_so, kind = 'cubic', fill_value = 'extrapolate')
cli_so = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, clic_g, kind = 'cubic', fill_value = 'extrapolate')
clic_g = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, clic_so, kind = 'cubic', fill_value = 'extrapolate')
clic_so = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, cli_frac_g, kind = 'cubic', fill_value = 'extrapolate')
cli_frac_g = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, cli_frac_so, kind = 'cubic', fill_value = 'extrapolate')
cli_frac_so = interpolated(constants.alt)

interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt[:27], clw_alt_lat[:,7:171] , kind = 'cubic')
clw_alt_lat = interpolated( constants.lat, constants.alt[:constants.liq_alt_confine] )
clw_alt_lat[:,0:7] = np.nan
clw_alt_lat[:,172:] = np.nan
clw_alt_lat[ clw_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp1d(raw_alt[:27], clw_frac_g, kind = 'cubic', fill_value = 'extrapolate')
clw_frac_g = interpolated(constants.alt[:constants.liq_alt_confine]) / cl_g[:constants.liq_alt_confine]

interpolated = interpolate.interp1d(raw_alt[:27], clw_frac_so, kind = 'cubic', fill_value = 'extrapolate')
clw_frac_so = interpolated(constants.alt[:constants.liq_alt_confine]) / cl_so[:constants.liq_alt_confine]

interpolated = interpolate.interp1d(raw_alt[:27], clw_g, kind = 'cubic', fill_value = 'extrapolate')
clw_g = interpolated(constants.alt[:constants.liq_alt_confine])

interpolated = interpolate.interp1d(raw_alt[:27], clw_so, kind = 'cubic', fill_value = 'extrapolate')
clw_so = interpolated(constants.alt[:constants.liq_alt_confine])

interpolated = interpolate.interp1d(raw_alt[:27], clwc_g, kind = 'cubic', fill_value = 'extrapolate')
clwc_g = interpolated(constants.alt[:constants.liq_alt_confine])

interpolated = interpolate.interp1d(raw_alt[:27], clwc_so, kind = 'cubic', fill_value = 'extrapolate')
clwc_so = interpolated(constants.alt[:constants.liq_alt_confine])

interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt, full_clw_alt_lat[:,7:171] , kind = 'cubic')
full_clw_alt_lat = interpolated( constants.lat, constants.alt )
full_clw_alt_lat[:,0:7] = np.nan
full_clw_alt_lat[:,172:] = np.nan
full_clw_alt_lat[ full_clw_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt, full_clwc_alt_lat[:,7:171] , kind = 'cubic')
full_clwc_alt_lat = interpolated( constants.lat, constants.alt )
full_clwc_alt_lat[:,0:7] = np.nan
full_clwc_alt_lat[:,172:] = np.nan
full_clwc_alt_lat[ full_clwc_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt, full_p_alt_lat[:,7:171] , kind = 'cubic')
full_p_alt_lat = interpolated( constants.lat, constants.alt )
full_p_alt_lat[:,0:7] = np.nan
full_p_alt_lat[:,172:] = np.nan
full_p_alt_lat[ full_p_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt, full_ta_alt_lat[:,7:171] , kind = 'cubic')
full_ta_alt_lat = interpolated( constants.lat, constants.alt )
full_ta_alt_lat[:,0:7] = np.nan
full_ta_alt_lat[:,172:] = np.nan
full_ta_alt_lat[ full_ta_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt[:27], ta_alt_lat[:,7:171] , kind = 'cubic')
ta_alt_lat = interpolated( constants.lat, constants.alt[:constants.liq_alt_confine] )
ta_alt_lat[:,0:7] = np.nan
ta_alt_lat[:,172:] = np.nan
ta_alt_lat[ ta_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d(constants.lat[7:171], raw_alt, density_alt_lat[:,7:171] , kind = 'cubic')
density_alt_lat = interpolated( constants.lat, constants.alt )
density_alt_lat[:,0:7] = np.nan
density_alt_lat[:,172:] = np.nan
density_alt_lat[ density_alt_lat < 0 ] = np.nan

cli_frac_alt_lat = cli_alt_lat / (cli_alt_lat + full_clw_alt_lat)
full_clw_frac_alt_lat = full_clw_alt_lat / (cli_alt_lat + full_clw_alt_lat)

stat = 'mean'
cli_t_g, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat.flatten(), cli_alt_lat.flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))
cli_t_so, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat[:,constants.so_idx_1:constants.so_idx_2].flatten(), cli_alt_lat[:,constants.so_idx_1:constants.so_idx_2].flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))

clw_frac_t_g = clw_t_g / ( clw_t_g + cli_t_g )
clw_frac_t_so = clw_t_so / ( clw_t_so + cli_t_so )

cli_frac_t_g = cli_t_g / ( clw_t_g + cli_t_g )
cli_frac_t_so = cli_t_so / ( clw_t_so + cli_t_so )



for index,key in enumerate(h5f.keys()):
    print (index, key)


