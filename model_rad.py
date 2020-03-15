

import numpy as np
import os
from netCDF4 import Dataset
from netCDF4 import date2index
import h5py
import math
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
from scipy import stats
import datetime
from pprint import pprint
from sklearn.impute import SimpleImputer
import constants
from scipy import ndimage as nd
import matplotlib.pyplot as plt
import climlab

start = datetime.datetime( 2007, 1, 1 )
end = datetime.datetime( 2010, 1, 1 )
location = constants.home
model = 'CMIP6-AMIP-GFDL-CM4'
os.chdir( location + 'Data/' + model )

#----------------------------- Aerosols ---------------------------#

o2=0.
ccl4=0.
cfc22=4.8743326488706363-11

# Mole Fraction of Ozone (time, lev, lat, lon)
variable = 'o3'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    o3 = 0.
else:
    with Dataset( check_file, 'r') as f:
        o3 = constants.extract_data_over_time( variable, f, start, end )
        o3[o3 > 1] = np.nan
        plev_o3 = constants.extract_data( 'plev', f ) # in Pa

# Global Mean Mole Fraction of CH4 (time)
variable = 'ch4global'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    ch4 = 1.557889297535934e-6
else:
    with Dataset( check_file, 'r') as f:
        ch4 = constants.extract_data_over_time( variable, f, start, end )
        ch4 = np.nanmean( ch4, axis = 0 ) # Average over time

# Global Mean Mole Fraction of N2O (time)
variable = 'n2oglobal'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    n2o = 3.007494866529774e-7
else:
    with Dataset( check_file, 'r') as f:
        n2o = constants.extract_data_over_time( variable, f, start, end )
        n2o = np.nanmean( n2o, axis = 0 ) # Average over time

# Global Mean Mole Fraction of CFC11 (time)
variable = 'cfc11global'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    cfc11 = 1.680e-10
else:
    with Dataset( check_file, 'r') as f:
        cfc11 = constants.extract_data_over_time( variable, f, start, end )
        cfc11 = np.nanmean( cfc11, axis = 0 ) # Average over time

# Global Mean Mole Fraction of CFC12 (time)
variable = 'cfc12global'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    cfc12 = 2.850e-10
else:
    with Dataset( check_file, 'r') as f:
        cfc12 = constants.extract_data_over_time( variable, f, start, end )
        cfc12 = np.nanmean( cfc12, axis = 0 ) # Average over time

# Global Mean Mole Fraction of CFC113 (time)
variable = 'cfc113global'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    cfc113 = 1.737268993839836e-11
else:
    with Dataset( check_file, 'r') as f:
        cfc113 = constants.extract_data_over_time( variable, f, start, end )
        cfc113 = np.nanmean( cfc113, axis = 0 ) # Average over time

# Total Atmospheric Mass (kg) of CO2 
variable = 'co2mass'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    co2mass = 400.e-6
else:
    with Dataset( check_file, 'r') as f:
        co2mass = constants.extract_data_over_time( variable, f, start, end )
        co2mass = np.nanmean( co2mass, axis = 0 ) # Average over time


#----------------------------- Cloud variables -----------------------------#

# Get layer cloud fraction and pressure level variables - (time, level, lat, lon)
with Dataset( constants.variable_to_filename( 'cl' ), 'r') as f:
    lat = constants.extract_data( 'lat', f )
    lon = constants.extract_data( 'lon', f )
    cl = constants.extract_data_over_time( 'cl', f, start, end ) / 100
    if 'CM4' in model or 'IPSL' in model:
        a = constants.extract_data( 'ap', f )
        b = constants.extract_data( 'b', f )
    else:
        a = constants.extract_data( 'a', f )
        b = constants.extract_data( 'b', f )
        p0 = np.array(f.variables['p0'][:])
        a = a*p0

# Get surface pressure - (time, lat, lon)
with Dataset( constants.variable_to_filename( 'ps' ), 'r') as f:
    ps = constants.extract_data_over_time( 'ps', f, start, end )
    ps = np.nanmean( ps, axis = 0 ) # average over time

# Get surface temperature - (time, lat, lon)
with Dataset( constants.variable_to_filename( 'ts' ), 'r') as f:
    ts = constants.extract_data_over_time( 'ts', f, start, end )

# Get specific humidity - (time, lev, lat, lon)
with Dataset( constants.variable_to_filename( 'hus' ), 'r') as f:
    hus = constants.extract_data_over_time( 'hus', f, start, end )
    hus[hus > 1] = np.nan
    plev_h = constants.extract_data( 'plev', f ) # in Pa

# Convert pressure level variables to pressure
a = np.swapaxes(np.swapaxes(np.tile(a, (np.shape(ps)[0], np.shape(ps)[1], 1) ), 0, 2), 1, 2)
b = np.swapaxes(np.swapaxes(np.tile(b, (np.shape(ps)[0], np.shape(ps)[1], 1) ), 0, 2), 1, 2)
ps = np.tile(ps, (np.shape(b)[0], 1, 1) )
p = (a + b*ps) # in Pa

# Get Cloud liquid water mass fraction in air (kg/kg) - (time, level, lat, lon)
with Dataset( constants.variable_to_filename( 'clw' ), 'r') as f:
    clw = constants.extract_data_over_time( 'clw', f, start, end )

# Get Cloud ice water mass fraction in air (kg/kg) - (time, level, lat, lon)
with Dataset( constants.variable_to_filename( 'cli' ), 'r') as f:
    cli = constants.extract_data_over_time( 'cli', f, start, end )

# Get temperature in atmosphere levels (time, level, lat, lon)
with Dataset( constants.variable_to_filename( 'ta' ), 'r' ) as f:
    ta = constants.extract_data_over_time( 'ta', f, start, end )
    ta[ta > 500] = np.nan
    plev_t = constants.extract_data( 'plev', f ) # in Pa




#-------------------------- Get top of atmosphere radiative flux data ------------------------#

# Incoming SW at TOA
with Dataset( constants.variable_to_filename( 'rsdt' ), 'r') as f:
    rsdt = constants.extract_data('rsdt', f )

# Outgoing SW at TOA
with Dataset( constants.variable_to_filename( 'rsut' ), 'r') as f:
    rsut = constants.extract_data('rsut', f )

# Outgoing SW at TOA assuming clear sky
with Dataset( constants.variable_to_filename( 'rsutcs' ), 'r') as f:
    rsutcs = constants.extract_data('rsutcs', f )

# Outgoing LW at TOA
with Dataset( constants.variable_to_filename( 'rlut' ), 'r') as f:
    rlut = constants.extract_data('rlut', f )

# Outgoing LW at TOA assuming clear sky
with Dataset( constants.variable_to_filename( 'rlutcs' ), 'r') as f:
    rlutcs = constants.extract_data('rlutcs', f )

# Net TOA flux
with Dataset( constants.variable_to_filename( 'rtmt' ), 'r') as f:
    rtmt = constants.extract_data('rtmt', f )


#-------------------------- Surface radiative flux data --------------------------#

# Incoming SW at surface
with Dataset( constants.variable_to_filename( 'rsds' ), 'r') as f:
    rsds = constants.extract_data('rsds', f )

# Incoming SW at surface assuming clear sky
with Dataset( constants.variable_to_filename( 'rsdscs' ), 'r') as f:
    rsdscs = constants.extract_data('rsdscs', f )

# Outgoing SW at surface
with Dataset( constants.variable_to_filename( 'rsus' ), 'r') as f:
    rsus = constants.extract_data('rsus', f )

# Outgoing SW at surface assuming clear sky
with Dataset( constants.variable_to_filename( 'rsuscs' ), 'r') as f:
    rsuscs = constants.extract_data('rsuscs', f )

# Outgoing LW at surface
with Dataset( constants.variable_to_filename( 'rlus' ), 'r') as f:
    rlus = constants.extract_data('rlus', f )

# Incoming LW at surface
with Dataset( constants.variable_to_filename( 'rlds' ), 'r') as f:
    rlds = constants.extract_data('rlds', f )

# Incoming LW at surface assuming clear sky
with Dataset( constants.variable_to_filename( 'rldscs' ), 'r') as f:
    rldscs = constants.extract_data('rldscs', f )


albedo_toa = rsut / rsdt
albedo_toa_cs = rsutcs / rsdt
albedo_surface = rsus / rsds
albedo_surface_cs = rsuscs / rsdscs


#-------------------------- State droplet and ice crystal sizes --------------------------#

r_liq = np.zeros((p.shape[0], lat.shape[0]))
r_ice = np.zeros((p.shape[0], lat.shape[0]))

r_liq[:] = 60
r_ice[:] = 60


#-------------------------- Average variables over time --------------------------#

cl = np.nanmean( cl, axis = 0 )
clw = np.nanmean( clw, axis = 0 )
cli = np.nanmean( cli, axis = 0 )
ta = np.nanmean( ta, axis = 0 )
ts = np.nanmean( ts, axis = 0 )
hus = np.nanmean( hus, axis = 0 )

rsdt = np.nanmean( rsdt, axis = 0 )
rsut = np.nanmean( rsut, axis = 0 )
rsutcs = np.nanmean( rsutcs, axis = 0 )
rlut = np.nanmean( rlut, axis = 0 )
rlutcs = np.nanmean( rlutcs, axis = 0 )
rtmt = np.nanmean( rtmt, axis = 0 )

rsds = np.nanmean( rsds, axis = 0 )
rsdscs = np.nanmean( rsdscs, axis = 0 )
rsus = np.nanmean( rsus, axis = 0 )
rsuscs = np.nanmean( rsuscs, axis = 0 )
rlus = np.nanmean( rlus, axis = 0 )
rlds = np.nanmean( rlds, axis = 0 )
rldscs = np.nanmean( rldscs, axis = 0 )

albedo_toa = np.nanmean( albedo_toa, axis = 0 )
albedo_toa_cs = np.nanmean( albedo_toa_cs, axis = 0 )
albedo_surface = np.nanmean( albedo_surface, axis = 0 )
albedo_surface_cs = np.nanmean( albedo_surface_cs, axis = 0 )
o3 = np.nanmean( o3, axis = 0 )

#-------------------------- Average variables over longitude --------------------------#

cl = np.nanmean( cl, axis = -1 )
clw = np.nanmean( clw, axis = -1 )
cli = np.nanmean( cli, axis = -1 )
ta = np.nanmean( ta, axis = -1 )
ts = np.nanmean( ts, axis = -1 )
hus = np.nanmean( hus, axis = -1 )

rsdt = np.nanmean( rsdt, axis = -1 )
rsut = np.nanmean( rsut, axis = -1 )
rsutcs = np.nanmean( rsutcs, axis = -1 )
rlut = np.nanmean( rlut, axis = -1 )
rlutcs = np.nanmean( rlutcs, axis = -1 )
rtmt = np.nanmean( rtmt, axis = -1 )

rsds = np.nanmean( rsds, axis = -1 )
rsdscs = np.nanmean( rsdscs, axis = -1 )
rsus = np.nanmean( rsus, axis = -1 )
rsuscs = np.nanmean( rsuscs, axis = -1 )
rlus = np.nanmean( rlus, axis = -1 )
rlds = np.nanmean( rlds, axis = -1 )
rldscs = np.nanmean( rldscs, axis = -1 )

albedo_toa = np.nanmean( albedo_toa, axis = -1 )
albedo_toa_cs = np.nanmean( albedo_toa_cs, axis = -1 )
albedo_surface = np.nanmean( albedo_surface, axis = -1 )
albedo_surface_cs = np.nanmean( albedo_surface_cs, axis = -1 )

o3 = np.nanmean( o3, axis = -1 )

# Get layer pressure levels (Pa)
p = constants.global3DMean(p, lat)

#-------------------------- Interpolate variables --------------------------#

imp = SimpleImputer(missing_values=np.nan, strategy='mean')
imp.fit(np.transpose(ta))  
ta = imp.transform(np.transpose(ta))
ta = np.transpose(ta)   
if plev_t.shape[0] > ta.shape[0]:
    plev = plev_t[plev_t.shape[0] - ta.shape[0]:] # reshape alt_temp if not equal
interpolated = interpolate.interp2d( lat, plev, ta, kind = 'linear')
ta = np.flip( interpolated( lat, p ), axis = 0 )

imp.fit(np.transpose(hus))  
hus = imp.transform(np.transpose(hus))
hus = np.transpose(hus)   
if plev_h.shape[0] > hus.shape[0]:
    plev_h = plev_h[plev_h.shape[0] - hus.shape[0]:] # reshape alt_temp if not equal
interpolated = interpolate.interp2d( lat, plev_h, hus, kind = 'linear')
hus = np.flip( interpolated( lat,p ), axis = 0 )

imp.fit(np.transpose(o3))  
o3 = imp.transform(np.transpose(o3))
o3 = np.transpose(o3)   
if plev_o3.shape[0] > o3.shape[0]:
    plev_o3 = plev_o3[plev_o3.shape[0] - o3.shape[0]:] # reshape alt_temp if not equal
interpolated = interpolate.interp2d( lat, plev_o3, o3, kind = 'linear')
o3 = np.flip( interpolated( lat,p ), axis = 0 )



#-------------------------- Input variables into radiative transfer code --------------------------#

sfc, atm = climlab.domain.zonal_mean_column(lat=lat, lev=p)

insolation = climlab.domain.Field(rsdt, domain=sfc)
albedo_surface = climlab.domain.Field(albedo_surface, domain=sfc)  
albedo_surface_cs = climlab.domain.Field(albedo_surface_cs, domain=sfc)  
albedo_toa = climlab.domain.Field(albedo_toa, domain=sfc)  
albedo_toa_cs = climlab.domain.Field(albedo_toa_cs, domain=sfc)  
ts = climlab.domain.Field(ts, domain=sfc)        
ta = climlab.domain.Field(np.transpose(ta), domain=atm)
specific_humidity = climlab.domain.Field(np.transpose(hus), domain=atm) # kg/kg
cldfrac = climlab.domain.Field(np.transpose(cl), domain=atm)
clwp = climlab.domain.Field(np.transpose(clw), domain=atm) # needs to be g/m2
ciwp = climlab.domain.Field(np.transpose(cli), domain=atm) # needs to be g/m2
o3 = climlab.domain.Field(np.transpose(o3), domain=atm)

r_liq = climlab.domain.Field(np.transpose(r_liq), domain=atm) # Cloud water drop effective radius (microns)
r_ice = climlab.domain.Field(np.transpose(r_ice), domain=atm) # Cloud ice particle effective size (microns)        

# dictionary of volumetric mixing ratios. Default values supplied if None
# If absorber_vmr = None then ozone will be interpolated to the model grid from a climatology file, or set to zero if ozone_file = None.
absorber = {'O3': o3, 'CO2': 400.e-6, 'CH4':ch4, 'N2O':n2o, 'O2': o2,'CCL4':ccl4, 
    'CFC11':cfc11, 'CFC12':cfc12, 'CFC113':cfc113, 'CFC22':cfc22}

#  State variables (Air and surface temperature)
state = {'Tatm': ta, 'Ts': ts}
h2o = climlab.radiation.ManabeWaterVapor(state=state)

rad = climlab.radiation.RRTMG(name='Radiation', 
                                state=state, 
                                specific_humidity=h2o.q, 
                                absorber_vmr=absorber, 
                                albedo=albedo_surface, 
                                insolation = insolation, 
                                cldfrac=cldfrac, 
                                r_liq=r_liq, 
                                r_ice=r_ice, 
                                clwp=clwp, 
                                ciwp=ciwp
                                )

rad.compute()

sw_down = rad.SW_flux_down
sw_down_clr = rad.SW_flux_down_clr
sw_up = rad.SW_flux_up
sw_up_clr = rad.SW_flux_up_clr

lw_down = rad.LW_flux_down
lw_up_clr = rad.LW_flux_up_clr
lw_down_clr = rad.LW_flux_down_clr

net_ASR = rad.ASR
net_OLR = rad.OLR

print(rad.SW_flux_down)

# rad_lw_dn[k,:,:]=rad_lw_dn[k,:,:]+rad.LW_flux_down.T/nt
# rad_lw_up_clr[k,:,:]=rad_lw_up_clr[k,:,:]+rad.LW_flux_up_clr.T/nt
# rad_lw_dn_clr[k,:,:]=rad_lw_dn_clr[k,:,:]+rad.LW_flux_down_clr.T/nt
# rad_sw_up[k,:,:]=rad_sw_up[k,:,:]+rad.SW_flux_up.T/nt
# rad_sw_dn[k,:,:]=rad_sw_dn[k,:,:]+rad.SW_flux_down.T/nt
# rad_sw_up_clr[k,:,:]=rad_sw_up_clr[k,:,:]+rad.SW_flux_up_clr.T/nt
# rad_sw_dn_clr[k,:,:]=rad_sw_dn_clr[k,:,:]+rad.SW_flux_down_clr.T/nt
