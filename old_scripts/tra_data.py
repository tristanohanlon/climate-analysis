import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import os
import constants
import numpy as np
import pprint
import climlab
from climlab.radiation import RRTMG
from climlab.radiation import RRTMG_SW

location = constants.laptop + 'Data/'

# Import dataset from file
with Dataset(location + 'lb_mean_19820101_month.nc', 'r') as f: 
    print(f.variables.keys())
    Tatm = f.variables['temp'][0, :, 0, 0]
    Ts = f.variables['t_surf'][0, :, 0, 0]
    SH = f.variables['sphum'][0, :, 0, 0]
    p = f.variables['pfull'][:]

    # Model radiative fluxes from file
    rsdt = 280.004044102902
    rsutcs = f.variables['swup_toa_clr'][0, :, 0, 0]
    rsdscs = f.variables['swdn_sfc_clr'][0, :, 0, 0]
    rlutcs = f.variables['olr_clr'][0, :, 0, 0]
    rldscs = f.variables['lwdn_sfc_clr'][0, :, 0, 0]
    rsut = f.variables['swup_toa'][0, :, 0, 0]
    rlut = f.variables['olr'][0, :, 0, 0]
    model_heating_rate = f.variables['tdt_sw_clr'][0, :, 0, 0]


# Set absorbers
o2 = 0.
o3 = 0.
ccl4 = 0.
ch4 = 0.
ccl4 = 0.
cfc11 = 0.
cfc12 = 0.
cfc113 = 0.
cfc22 = 0.
ch4 = 0.
n2o = 0.
co2 = 400.e-6


absorbermean = {'O3': o3, 'CO2': co2, 'CH4':ch4, 'N2O':n2o, 'O2': o2,'CCL4':ccl4, 
    'CFC11':cfc11, 'CFC12':cfc12, 'CFC113':cfc113, 'CFC22':cfc22}

state = climlab.column_state(lev = p)
state.Tatm[:] = Tatm
state.Ts[:] = Ts

control_rad = RRTMG_SW(state=state,
            insolation = 280.004044102902,
            albedo = 0,
            absorber_vmr = absorbermean,
            specific_humidity = SH,
            verbose=False,
            )

control_rad.compute_diagnostics()

control_rsdt = control_rad.SW_flux_down_clr[0]
control_rsut = control_rad.SW_flux_up[0]
control_rsutcs = control_rad.SW_flux_up_clr[0]
control_rsdscs = control_rad.SW_flux_down_clr[-1]
heating_rate = control_rad.heating_rate
# control_rlut = control_rad.LW_flux_up[0]
# control_rlutcs = control_rad.LW_flux_up_clr[0]
# control_rldscs = control_rad.LW_flux_down_clr[-1]

#  Print and check mean data - single model only

print('===RRTMG OUTPUT===')
print('rsdt:')
print(control_rsdt)
print('rsutcs:')
print(control_rsutcs)
print('rsdscs:')
print(control_rsdscs)
print('rlutcs:')
print(control_rlutcs)
print('rldscs:')
print(control_rldscs)
print('rsut:')
print(control_rsut)
print('rlut:')
print(control_rlut)

print('===MODEL DATA===')
print('rsdt:')
print(rsdt)
print('rsutcs:')
print(rsutcs)
print('rsdscs:')
print(rsdscs)
print('rlutcs:')
print(rlutcs)
print('rldscs:')
print(rldscs)
print('rsut:')
print(rsut)
print('rlut:')
print(rlut)


