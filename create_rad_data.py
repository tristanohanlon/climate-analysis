
import datetime
import constants
import model_rad


start = datetime.datetime( 2007, 1, 1 )
end = datetime.datetime( 2010, 1, 1 )
location = constants.home
model = 'CMIP6-AMIP-GFDL-CM4'
label = '15'

# confine latitudes between:
min_lat = -70
max_lat = 70

use_aerosol_files = False
use_surface_albedo = True
plot_diagnostic_data = False
save_outputs = True # save output graphs and global mean data to excel
show_plots = False

# set ice and liquid droplet radius in microns
liquid_r = 30 # below 60 microns
ice_r = 130 # above 15, below 130

# convert specific_humidity to g/kg and check

model_rad.radiation( start, 
                    end, 
                    location, 
                    model, 
                    label, 
                    min_lat,
                    max_lat,
                    use_aerosol_files, 
                    use_surface_albedo, 
                    plot_diagnostic_data,
                    save_outputs,
                    liquid_r,
                    ice_r,
                    show_plots )