"""
Created on Sun Oct  6 14:44:25 2019

@author: Tristan O'Hanlon - University of Auckland

"""

import h5py
import os
import matplotlib.pyplot as plt
import constants
import numpy as np
# import cartopy.crs as ccrs
import datetime
from netCDF4 import Dataset
from netCDF4 import date2index

#--- Set Location, date period and model ---#

# specify model from the list above
# specify location: home, uni, hdd, laptop
model = 'CMIP6-AMIP-GFDL-CM4' # see comments above for available models
location = constants.hdd # home, uni, hdd or laptop
os.chdir( location + 'Data/' + model )
variable = 'cl'

# Get the list of all files and directories in current dir
dir_list = os.listdir() 
  
# print the list 
print(dir_list) 

# search the file names and select the one that has the variable name

start = datetime.datetime( 2006, 1, 1 )
end = datetime.datetime( 2010, 12, 1 )

dataset_filenames = []
for f in dir_list:
        if f.startswith(variable + '_'):
                dataset_filenames.append( f )

if len( dataset_filenames ) != 1:
        print( "Too many datasets match this variable name")
        exit(1)

def string_to_datetime( s ):
        year = int( s[:4] )
        month = int( s[4:] )
        return datetime.datetime( year , month, 1 )


def decode_date_range_from_filename( f ):
        date_start_pos = f.rfind( '_')
        if date_start_pos == -1:
                return None,None
        date_start_pos += 1
        date_start_end = f.rfind( '-')
        if date_start_end == -1:
                return None,None
        date_end_pos = date_start_end + 1
        date_end_end = f.rfind('.')
        if date_end_end == -1:
                return None,None
        date_start = f[date_start_pos:date_start_end]
        date_end = f[date_end_pos:date_end_end]
        if len( date_start ) != 6 or len( date_end ) != 6:
                return None,None        
        return string_to_datetime( date_start ),string_to_datetime( date_end )
        

start,end = decode_date_range_from_filename( dataset_filenames[0])

with Dataset( dataset_filenames[0], 'r') as f:
        raw_lat = constants.extract_data( 'lat', f)
        raw_lon = constants.extract_data( 'lon', f)

        data = np.nanmean( constants.extract_data_over_time(variable, f, start, end ), axis = 0 ) # average over time
        data = constants.global3DMean(data, raw_lat)
        print(data)


fig, ax = plt.subplots()
ax.plot( data, raw_lev )
ax.set_ylabel('Altitude (km)')
plt.grid(True)
plt.show()





