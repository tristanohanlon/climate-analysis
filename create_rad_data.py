
from cftime import DatetimeNoLeap
import constants
import model_rad_dt # change to model_rad for pre average
from netCDF4 import date2index
import datetime


start = DatetimeNoLeap( 2007, 1, 1 )
end = DatetimeNoLeap( 2010, 12, 1 )
start_dt = datetime.datetime( 2007, 1, 1 )
end_dt = datetime.datetime( 2010, 12, 1 )
location = constants.home

#  Access keys by models.keys() and values by models.values()
#  Get both out by: for name, path in zip(models.keys(), models.values()):

models = {
    "CESM2-CAM6" : "CMIP6-AMIP-CESM2-CAM6",
    "GFDL-CM4" : "CMIP6-AMIP-GFDL-CM4",
    "BCC-ESM1" : "CMIP6-AMIP-BCC-ESM1",
    "IPSL-CM6A-LR" : "CMIP6-AMIP-IPSL-CM6A-LR",
    "MRI-ESM2" : "CMIP6-AMIP-MRI-ESM2",
    "MIROC6" : "CMIP6-AMIP-MIROC6",
}

label = '2'

# confine latitudes between:
lat_bnd_1 = -80
lat_bnd_2 = 80

save_outputs = True # save output graphs to a pdf and global mean data to excel

# set ice and liquid droplet radius in microns
liquid_r = 30 # below 60 microns
ice_r = 60 # above 15, below 130

model_rad_dt.radiation( start, 
                    end,
                    start_dt, 
                    end_dt, 
                    location, 
                    models, 
                    label, 
                    lat_bnd_1,
                    lat_bnd_2,
                    save_outputs,
                    liquid_r,
                    ice_r )