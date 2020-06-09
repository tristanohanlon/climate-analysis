import numpy as np
import matplotlib.pyplot as plt
import climlab
from climlab.radiation import RRTMG
import xarray as xr
import constants
import os
from netCDF4 import Dataset

location = constants.home
model = 'CMIP6-AMIP-GFDL-CM4'
os.chdir( location + 'Data/' + model )


#  Take global, annual average

with xr.open_dataset(constants.variable_to_filename( 'cl' ), decode_times=False) as cl:
    weight = np.cos(np.deg2rad(cl.lat)) / np.cos(np.deg2rad(cl.lat)).mean(dim='lat')
    CLglobal = (cl.cl * weight).mean(dim=('lat','lon','time')) / 100
    a = cl.ap
    b = cl.b

with xr.open_dataset(constants.variable_to_filename( 'ps' ), decode_times=False) as ps:
    ps = ps.ps

p = a + b*ps
plev = (p * weight).mean(dim=('lat','lon','time'))

with xr.open_dataset(constants.variable_to_filename( 'clw' ), decode_times=False) as clw:
    CLWglobal = (clw.clw * weight).mean(dim=('lat','lon','time'))

with xr.open_dataset(constants.variable_to_filename( 'cli' ), decode_times=False) as cli:
    CLIglobal = (cli.cli * weight).mean(dim=('lat','lon','time'))

with xr.open_dataset(constants.variable_to_filename( 'ta' ), decode_times=False) as ta:
    Tglobal = (ta.ta * weight).mean(dim=('lat','lon','time'))

with xr.open_dataset(constants.variable_to_filename( 'hus' ), decode_times=False) as hus:
    SHglobal = (hus.hus * weight).mean(dim=('lat','lon','time'))  # kg/kg

with xr.open_dataset(constants.variable_to_filename( 'o3' ), decode_times=False) as o3:
    O3global = (o3.o3 * weight).mean(dim=('lat','lon','time'))  # kg/kg

state = climlab.column_state(num_lev=50)
lev = state.Tatm.domain.axes['lev'].points * 100

# interpolate to model pressure levels
Tinterp = np.interp(lev, np.flipud(Tglobal.plev[:18]), np.flipud(Tglobal[:18]))
SHinterp = np.interp(lev, np.flipud(SHglobal.plev[:18]), np.flipud(SHglobal[:18]))
O3interp = np.interp(lev, np.flipud(O3global.plev[:18]), np.flipud(O3global[:18]))
CLinterp = np.interp(lev, np.flipud(plev), np.flipud(CLglobal))
CLWinterp = np.interp(lev, np.flipud(plev), np.flipud(CLWglobal))
CLIinterp = np.interp(lev, np.flipud(plev), np.flipud(CLIglobal))
#  Need to 'flipud' because the interpolation routine 
#  needs the pressure data to be in increasing order
#  Create a state dictionary with corresponsing number of cl data levels

#  Set the temperature to the observed values
state.Tatm[:] = Tinterp

r_liq = 14.  # Cloud water drop effective radius (microns)
r_ice = 60.  # Cloud water drop effective radius (microns)

#  Convert mixing ratio (kg/kg) into cloud water content in (g/m3)

alt = constants.p_to_alt( lev )  * 1000 # in km
air_density = ((lev) / (286.9 * Tinterp))
clwc = (CLWinterp * air_density) * 1000
ciwc = (CLIinterp * air_density) * 1000

prev_alt = alt[0] + (alt[0] - alt[1])
CLWPglobal = np.zeros_like(lev)
for i, ( row, altitude ) in enumerate( zip( clwc, alt ) ):
    CLWPglobal[i] = row * (prev_alt - altitude)
    prev_alt = altitude

prev_alt = alt[-1] + (alt[-1] - alt[-2])
CIWPglobal = np.zeros_like(lev)
for i, ( row, altitude ) in enumerate( zip( ciwc, alt ) ):
    CIWPglobal[i] = row * (prev_alt - altitude)
    prev_alt = altitude
#  Convert cloud water content to water path (g/m2)
#  WP = WC * deltaz
# CLWPglobal = constants.wc_to_wp( clwc, alt )
# CIWPglobal = constants.wc_to_wp( ciwc, alt )


#  Loop through all pressure levels
#  Set up a radiation model with the cloud layer at the current pressure level
#  Compute CRE and store the results
CRE_LW = {}
CRE_SW = {}
OLR = np.zeros_like(lev)
ASR = np.zeros_like(lev)
OLRclr = np.zeros_like(lev)
ASRclr = np.zeros_like(lev)
for i in range(lev.size):
    # Whole-column cloud characteristics
    # The cloud fraction is a Gaussian bump centered at the current level        
    rad = RRTMG(state=state, 
                albedo=0.2,
                O3 = O3interp,
                specific_humidity=SHinterp,
                cldfrac=CLinterp*np.exp(-(lev-lev[i])**2/(2*25.)**2),
                verbose=False,
                clwp = CLWPglobal *1000,
                ciwp = CIWPglobal *1000,
                r_liq = np.zeros_like(state.Tatm) + r_liq,
                r_ice = np.zeros_like(state.Tatm) + r_ice,               
                )
    rad.compute_diagnostics()
    OLR[i] = rad.OLR
    OLRclr[i] = rad.OLRclr
    ASR[i] = rad.ASR
    ASRclr[i] = rad.ASRclr
CRE_LW = -(OLR - OLRclr)
CRE_SW = (ASR - ASRclr)

#  Make some plots of the CRE dependence on cloud height
fig, axes = plt.subplots(1,3, figsize=(16,6))
ax = axes[0]
ax.plot(CRE_LW, lev)
ax.set_ylabel('Pressure (hPa)')
ax.set_xlabel('LW cloud radiative effect (W/m2)')

ax = axes[1]
ax.plot(CRE_SW, lev)
ax.set_xlabel('SW cloud radiative effect (W/m2)')

ax = axes[2]
ax.plot(CRE_SW + CRE_LW, lev)
ax.set_xlabel('Net cloud radiative effect (W/m2)')

for ax in axes:
    ax.invert_yaxis()
    ax.legend()
    ax.grid()
fig.suptitle('Cloud Radiative Effect as a function of the vertical height of the cloud layer', fontsize=16)
plt.show()