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

ta_file = constants.variable_to_filename( 'ta' )
hus_file = constants.variable_to_filename( 'hus' )
cl_file = constants.variable_to_filename( 'cl' )


with xr.open_dataset(ta_file, decode_times=False) as ta:
    weight = np.cos(np.deg2rad(ta.lat)) / np.cos(np.deg2rad(ta.lat)).mean(dim='lat')
    Tglobal = (ta.ta * weight).mean(dim=('lat','lon','time'))

with xr.open_dataset(hus_file, decode_times=False) as hus:
    SHglobal = (hus.hus * weight).mean(dim=('lat','lon','time'))  # kg/kg

with xr.open_dataset(cl_file, decode_times=False) as cl:
    CLglobal = (cl.cl * weight).mean(dim=('lat','lon','time')) / 100

#  Take global, annual average and convert to correct units (Kelvin and kg/kg)
#  Create a state dictionary with 50 levels
state = climlab.column_state(num_lev=50)
lev = state.Tatm.domain.axes['lev'].points * 100

# interpolate to model pressure levels
Tinterp = np.interp(lev, np.flipud(Tglobal.plev), np.flipud(Tglobal))
SHinterp = np.interp(lev, np.flipud(SHglobal.plev), np.flipud(SHglobal))
CLinterp = np.interp(lev, np.flipud(CLglobal.lev * 100000), np.flipud(CLglobal))
#  Need to 'flipud' because the interpolation routine 
#  needs the pressure data to be in increasing order

#  Set the temperature to the observed values
state.Tatm[:] = Tinterp
# cldfrac = np.zeros_like(state.Tatm) + CLinterp
#  Define some local cloud characteristics
#  We are going to repeat the calculation 
#   for three different types of clouds:
#   thin, medium, and thick
# cldfrac = 0.5  # layer cloud fraction
r_liq = 14.  # Cloud water drop effective radius (microns)
# in-cloud liquid water path (g/m2)
clwp = {'thin': 20.,
        'med': 60.,
        'thick': 200.,}

#  Loop through three types of cloud
#  for each type, loop through all pressure levels
#  Set up a radiation model with the cloud layer at the current pressure level
#  Compute CRE and store the results
CRE_LW = {}
CRE_SW = {}
for thickness in clwp:
    OLR = np.zeros_like(lev)
    ASR = np.zeros_like(lev)
    OLRclr = np.zeros_like(lev)
    ASRclr = np.zeros_like(lev)
    for i in range(lev.size):
        # Whole-column cloud characteristics
        #  The cloud fraction is a Gaussian bump centered at the current level        
        mycloud = {'clwp': np.zeros_like(state.Tatm) + clwp[thickness],
                   'r_liq': np.zeros_like(state.Tatm) + r_liq,}
        rad = RRTMG(state=state, 
                    albedo=0.2,
                    specific_humidity=SHinterp,
                    cldfrac=CLinterp*np.exp(-(lev-lev[i])**2/(2*25.)**2),
                    verbose=False,
                    **mycloud)
        rad.compute_diagnostics()
        OLR[i] = rad.OLR
        OLRclr[i] = rad.OLRclr
        ASR[i] = rad.ASR
        ASRclr[i] = rad.ASRclr
    CRE_LW[thickness] = -(OLR - OLRclr)
    CRE_SW[thickness] = (ASR - ASRclr)

#  Make some plots of the CRE dependence on cloud height
fig, axes = plt.subplots(1,3, figsize=(16,6))
ax = axes[0]
for thickness in clwp:
    ax.plot(CRE_LW[thickness], lev, label=thickness)
ax.set_ylabel('Pressure (hPa)')
ax.set_xlabel('LW cloud radiative effect (W/m2)')

ax = axes[1]
for thickness in clwp:
    ax.plot(CRE_SW[thickness], lev, label=thickness)
ax.set_xlabel('SW cloud radiative effect (W/m2)')

ax = axes[2]
for thickness in clwp:
    ax.plot(CRE_SW[thickness] + CRE_LW[thickness], lev, label=thickness)
ax.set_xlabel('Net cloud radiative effect (W/m2)')

for ax in axes:
    ax.invert_yaxis()
    ax.legend()
    ax.grid()
fig.suptitle('Cloud Radiative Effect as a function of the vertical height of the cloud layer', fontsize=16)
plt.show()