import datetime
import constants
import model_rad


start = datetime.datetime( 2007, 1, 1 )
end = datetime.datetime( 2010, 1, 1 )
location = constants.home
model = 'CMIP6-AMIP-GFDL-CM4'
label = '1'

use_aerosol_files = True
use_surface_albedo = True
plot_diagnostic_data = True
plot_outputs = True

# set ice and liquid droplet radius in microns
liquid_r = 60
ice_r = 60

model_rad.radiation( start, 
                    end, 
                    location, 
                    model, 
                    label, 
                    use_aerosol_files, 
                    use_surface_albedo, 
                    plot_diagnostic_data,
                    plot_outputs,
                    liquid_r,
                    ice_r )