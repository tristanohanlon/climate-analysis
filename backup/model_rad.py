

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
print('ch4 = ' + str(ch4))

# Global Mean Mole Fraction of N2O (time)
variable = 'n2oglobal'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    n2o = 3.007494866529774e-7
else:
    with Dataset( check_file, 'r') as f:
        n2o = constants.extract_data_over_time( variable, f, start, end )
        n2o = np.nanmean( n2o, axis = 0 ) # Average over time
print('n2o = ' + str(n2o))

# Global Mean Mole Fraction of CFC11 (time)
variable = 'n2oglobal'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    cfc11 = 1.680e-10
else:
    with Dataset( check_file, 'r') as f:
        cfc11 = constants.extract_data_over_time( variable, f, start, end )
        cfc11 = np.nanmean( cfc11, axis = 0 ) # Average over time
print('cfc11 = ' + str(cfc11))

# Global Mean Mole Fraction of CFC12 (time)
variable = 'cfc12global'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    cfc12 = 2.850e-10
else:
    with Dataset( check_file, 'r') as f:
        cfc12 = constants.extract_data_over_time( variable, f, start, end )
        cfc12 = np.nanmean( cfc12, axis = 0 ) # Average over time
print('cfc12 = ' + str(cfc12))

# Global Mean Mole Fraction of CFC113 (time)
variable = 'cfc113global'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    cfc113 = 1.737268993839836e-11
else:
    with Dataset( check_file, 'r') as f:
        cfc113 = constants.extract_data_over_time( variable, f, start, end )
        cfc113 = np.nanmean( cfc113, axis = 0 ) # Average over time
print('cfc113 = ' + str(cfc113))

# Total Atmospheric Mass (kg) of CO2 
variable = 'co2mass'
check_file = constants.variable_to_filename( variable )
if check_file == None:
    co2mass = 400.e-6
else:
    with Dataset( check_file, 'r') as f:
        co2mass = constants.extract_data_over_time( variable, f, start, end )
        co2mass = np.nanmean( co2mass, axis = 0 ) # Average over time
print('co2mass = ' + str(co2mass))


#----------------------------- Cloud variables -----------------------------#

# Get layer cloud fraction and pressure level variables - (time, level, lat, lon)
with Dataset( constants.variable_to_filename( 'cl' ), 'r') as f:
    lat = constants.extract_data( 'lat', f )
    # define latitude confinement
    lat_confine_1 = np.abs(lat - (-72)).argmin()
    lat_confine_2 = np.abs(lat - (75)).argmin()
    lat = lat[lat_confine_1:lat_confine_2]

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
    cl = cl[:,:,lat_confine_1:lat_confine_2,:]

# Get surface pressure - (time, lat, lon)
with Dataset( constants.variable_to_filename( 'ps' ), 'r') as f:
    ps = constants.extract_data_over_time( 'ps', f, start, end )
    ps = np.nanmean( ps, axis = 0 ) # average over time

# Get surface temperature - (time, lat, lon)
with Dataset( constants.variable_to_filename( 'ts' ), 'r') as f:
    ts = constants.extract_data_over_time( 'ts', f, start, end )
    ts = ts[:,lat_confine_1:lat_confine_2,:]

# Get specific humidity - (time, lev, lat, lon)
with Dataset( constants.variable_to_filename( 'hus' ), 'r') as f:
    hus = constants.extract_data_over_time( 'hus', f, start, end )
    hus[hus > 1] = np.nan
    plev_h = constants.extract_data( 'plev', f ) # in Pa
    hus = hus[:,:,lat_confine_1:lat_confine_2,:]

# Convert pressure level variables to pressure
a = np.swapaxes(np.swapaxes(np.tile(a, (np.shape(ps)[0], np.shape(ps)[1], 1) ), 0, 2), 1, 2)
b = np.swapaxes(np.swapaxes(np.tile(b, (np.shape(ps)[0], np.shape(ps)[1], 1) ), 0, 2), 1, 2)
ps = np.tile(ps, (np.shape(b)[0], 1, 1) )
p = (a + b*ps) # in Pa
p = p[:,lat_confine_1:lat_confine_2,:]


# Get Cloud liquid water mass fraction in air (kg/kg) - (time, level, lat, lon)
with Dataset( constants.variable_to_filename( 'clw' ), 'r') as f:
    clw = constants.extract_data_over_time( 'clw', f, start, end )
    clw = clw[:,:,lat_confine_1:lat_confine_2,:]

# Get Cloud ice water mass fraction in air (kg/kg) - (time, level, lat, lon)
with Dataset( constants.variable_to_filename( 'cli' ), 'r') as f:
    cli = constants.extract_data_over_time( 'cli', f, start, end )
    cli = cli[:,:,lat_confine_1:lat_confine_2,:]

# Get temperature in atmosphere levels (time, level, lat, lon)
with Dataset( constants.variable_to_filename( 'ta' ), 'r' ) as f:
    ta = constants.extract_data_over_time( 'ta', f, start, end )
    ta[ta > 500] = np.nan
    plev_t = constants.extract_data( 'plev', f ) # in Pa
    ta = ta[:,:,lat_confine_1:lat_confine_2,:]




#-------------------------- Get top of atmosphere radiative flux data ------------------------#

# Incoming SW at TOA
with Dataset( constants.variable_to_filename( 'rsdt' ), 'r') as f:
    rsdt = constants.extract_data('rsdt', f )
    rsdt = rsdt[:,lat_confine_1:lat_confine_2,:]

# Outgoing SW at TOA
with Dataset( constants.variable_to_filename( 'rsut' ), 'r') as f:
    rsut = constants.extract_data('rsut', f )
    rsut = rsut[:,lat_confine_1:lat_confine_2,:]

# Outgoing SW at TOA assuming clear sky
with Dataset( constants.variable_to_filename( 'rsutcs' ), 'r') as f:
    rsutcs = constants.extract_data('rsutcs', f )
    rsutcs = rsutcs[:,lat_confine_1:lat_confine_2,:]

# Outgoing LW at TOA
with Dataset( constants.variable_to_filename( 'rlut' ), 'r') as f:
    rlut = constants.extract_data('rlut', f )
    rlut = rlut[:,lat_confine_1:lat_confine_2,:]

# Outgoing LW at TOA assuming clear sky
with Dataset( constants.variable_to_filename( 'rlutcs' ), 'r') as f:
    rlutcs = constants.extract_data('rlutcs', f )
    rlutcs = rlutcs[:,lat_confine_1:lat_confine_2,:]

# Net TOA flux
with Dataset( constants.variable_to_filename( 'rtmt' ), 'r') as f:
    rtmt = constants.extract_data('rtmt', f )
    rtmt = rtmt[:,lat_confine_1:lat_confine_2,:]


#-------------------------- Surface radiative flux data --------------------------#

# Incoming SW at surface
with Dataset( constants.variable_to_filename( 'rsds' ), 'r') as f:
    rsds = constants.extract_data('rsds', f )
    rsds = rsds[:,lat_confine_1:lat_confine_2,:]

# Incoming SW at surface assuming clear sky
with Dataset( constants.variable_to_filename( 'rsdscs' ), 'r') as f:
    rsdscs = constants.extract_data('rsdscs', f )
    rsdscs = rsdscs[:,lat_confine_1:lat_confine_2,:]

# Outgoing SW at surface
with Dataset( constants.variable_to_filename( 'rsus' ), 'r') as f:
    rsus = constants.extract_data('rsus', f )
    rsus = rsus[:,lat_confine_1:lat_confine_2,:]

# Outgoing SW at surface assuming clear sky
with Dataset( constants.variable_to_filename( 'rsuscs' ), 'r') as f:
    rsuscs = constants.extract_data('rsuscs', f )
    rsuscs = rsuscs[:,lat_confine_1:lat_confine_2,:]

# Outgoing LW at surface
with Dataset( constants.variable_to_filename( 'rlus' ), 'r') as f:
    rlus = constants.extract_data('rlus', f )
    rlus = rlus[:,lat_confine_1:lat_confine_2,:]

# Incoming LW at surface
with Dataset( constants.variable_to_filename( 'rlds' ), 'r') as f:
    rlds = constants.extract_data('rlds', f )
    rlds = rlds[:,lat_confine_1:lat_confine_2,:]

# Incoming LW at surface assuming clear sky
with Dataset( constants.variable_to_filename( 'rldscs' ), 'r') as f:
    rldscs = constants.extract_data('rldscs', f )
    rldscs = rldscs[:,lat_confine_1:lat_confine_2,:]


albedo_toa = rsut / rsdt
albedo_toa_cs = rsutcs / rsdt
albedo_surface = rsus / rsds
albedo_surface_cs = rsuscs / rsdscs


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
o3 = o3[:,lat_confine_1:lat_confine_2,:]

#-------------------------- Average variables over longitude --------------------------#

cl = np.nanmean( cl, axis = -1 )
clw = np.nanmean( clw, axis = -1 )
cli = np.nanmean( cli, axis = -1 )
ta = np.nanmean( ta, axis = -1 )
ts = np.nanmean( ts, axis = -1 )
hus = np.nanmean( hus, axis = -1 )
plev = np.nanmean( p, axis = -1 )

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
    plev_t = plev_t[plev_t.shape[0] - ta.shape[0]:] # reshape alt_temp if not equal
interpolated = interpolate.interp2d( lat, plev_t, ta, kind = 'linear')
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


#-------------------------- State droplet and ice crystal sizes --------------------------#

# convert mixing ratio (kg/kg) into cloud water content in (g/m3)
clwc = constants.mix_ratio_to_water_content( clw, ta, plev )
ciwc = constants.mix_ratio_to_water_content( cli, ta, plev )
alt = constants.p_to_alt( p ) # in km

# convert cloud water content to water path (g/m2)
clwp = constants.wc_to_wp( clwc, alt )
ciwp = constants.wc_to_wp( ciwc, alt )


r_liq = np.zeros((p.shape[0], lat.shape[0]))
r_ice = np.zeros((p.shape[0], lat.shape[0]))

r_liq[:] = 60
r_ice[:] = 60

#-------------------------- test plots --------------------------#


# fig, ax = plt.subplots()
# cont = ax.contourf( lat, alt, o3, cmap='coolwarm' )
# ax.set_xlabel('Latitude')
# ax.set_ylabel('Altitude (km)')
# ax.set_title( 'o3 mole fraction' )
# cbar = fig.colorbar(cont, orientation='horizontal')
# cbar.set_label('o3 mole fraction')
# plt.savefig( location + '/Images/RRTGM/' + "o3.svg", format="svg", bbox_inches='tight')
# plt.show()

# fig, ax = plt.subplots()
# cont = ax.contourf( lat, alt, ta, cmap='coolwarm' )
# ax.set_xlabel('Latitude')
# ax.set_ylabel('Altitude (km)')
# ax.set_title( 'atmosphere temperature' )
# cbar = fig.colorbar(cont, orientation='horizontal')
# cbar.set_label('ta (K)')
# plt.savefig( location + '/Images/RRTGM/' + "ta.svg", format="svg", bbox_inches='tight')
# plt.show()

# fig, ax = plt.subplots()
# cont = ax.contourf( lat, alt, hus, cmap='coolwarm' )
# ax.set_xlabel('Latitude')
# ax.set_ylabel('Altitude (km)')
# ax.set_title( 'atmosphere specific humidity' )
# cbar = fig.colorbar(cont, orientation='horizontal')
# cbar.set_label('atmosphere specific humidity (kg/kg)')
# plt.savefig( location + '/Images/RRTGM/' + "hus.svg", format="svg", bbox_inches='tight')
# plt.show()

# fig, ax = plt.subplots()
# cont = ax.contourf( lat, alt, cl, cmap='coolwarm' )
# ax.set_xlabel('Latitude')
# ax.set_ylabel('Altitude (km)')
# ax.set_title( 'cloud fraction' )
# cbar = fig.colorbar(cont, orientation='horizontal')
# cbar.set_label('cloud fraction')
# plt.savefig( location + '/Images/RRTGM/' + "cl.svg", format="svg", bbox_inches='tight')
# plt.show()

#-------------------------- Input variables into radiative transfer code --------------------------#

sfc, atm = climlab.domain.zonal_mean_column(lat=lat, lev=p)

# surface variables
insolation_field = climlab.domain.Field(rsdt, domain=sfc)
albedo_surface_field = climlab.domain.Field(albedo_surface, domain=sfc)  
albedo_surface_cs_field = climlab.domain.Field(albedo_surface_cs, domain=sfc)  
albedo_toa_field = climlab.domain.Field(albedo_toa, domain=sfc)  
albedo_toa_cs_field = climlab.domain.Field(albedo_toa_cs, domain=sfc)  
ts_field = climlab.domain.Field(ts, domain=sfc)

# atmosphere layer variables
ta_field = climlab.domain.Field(np.transpose(ta), domain=atm)
specific_humidity_field = climlab.domain.Field(np.transpose(hus), domain=atm) # kg/kg
cldfrac_field = climlab.domain.Field(np.transpose(cl), domain=atm)
clwp_field = climlab.domain.Field(np.transpose(clwp), domain=atm) # needs to be g/m2
ciwp_field = climlab.domain.Field(np.transpose(ciwp), domain=atm) # needs to be g/m2
o3_field = climlab.domain.Field(np.transpose(o3), domain=atm)
r_liq_field = climlab.domain.Field(np.transpose(r_liq), domain=atm) # Cloud water drop effective radius (microns)
r_ice_field = climlab.domain.Field(np.transpose(r_ice), domain=atm) # Cloud ice particle effective size (microns)        

# dictionary of volumetric mixing ratios. Default values supplied if None
# if absorber_vmr = None then ozone will be interpolated to the model grid from a climatology file, or set to zero if ozone_file = None.
absorber = {'O3': o3_field, 'CO2': 400.e-6, 'CH4':ch4, 'N2O':n2o, 'O2': o2,'CCL4':ccl4, 
    'CFC11':cfc11, 'CFC12':cfc12, 'CFC113':cfc113, 'CFC22':cfc22}
# absorber = {'O3': o3_field, 'CO2': 400.e-6, 'CH4':1.557889297535934e-6, 'N2O':3.007494866529774e-7, 'O2': o2,'CCL4':ccl4, 
#     'CFC11':1.680e-10, 'CFC12':2.850e-10, 'CFC113':1.737268993839836e-11, 'CFC22':cfc22}

#  state variables (Air and surface temperature)
state = {'Tatm': ta_field, 'Ts': ts_field}

#-------------------------- Execute Control Data --------------------------#

control_rad = climlab.radiation.RRTMG(name='Radiation', 
                                state=state, 
                                specific_humidity=specific_humidity_field, 
                                absorber_vmr=absorber,
                                albedo=albedo_surface_field, 
                                insolation = insolation_field, 
                                cldfrac=cldfrac_field, 
                                r_liq=r_liq_field, 
                                r_ice=r_ice_field, 
                                clwp=clwp_field, 
                                ciwp=ciwp_field
                                )

control_rad.compute()

sw_down = control_rad.SW_flux_down
sw_down[sw_down>1000]=np.nan
sw_down[sw_down<-1000]=np.nan

sw_down_clr = control_rad.SW_flux_down_clr
sw_down_clr[sw_down_clr>500]=np.nan
sw_down_clr[sw_down_clr<-1000]=np.nan

sw_net = control_rad.SW_flux_net
sw_net[sw_net>500]=np.nan
sw_net[sw_net<-1000]=np.nan

sw_up = control_rad.SW_flux_up
sw_up[sw_up>500]=np.nan
sw_up[sw_up<-1000]=np.nan

sw_up_clr = control_rad.SW_flux_up_clr
sw_up_clr[sw_up_clr>500]=np.nan
sw_up_clr[sw_up_clr<-500]=np.nan


lw_down = control_rad.LW_flux_down
lw_down_clr = control_rad.LW_flux_down_clr
lw_up_clr = control_rad.LW_flux_up_clr
lw_up = control_rad.LW_flux_up
lw_net = control_rad.LW_flux_net


net_ASR = control_rad.ASR
net_ASRcld = control_rad.ASRcld
net_ASRclr = control_rad.ASRclr   

net_OLR = control_rad.OLR
net_OLRcld = control_rad.OLRcld
net_OLRclr = control_rad.OLRclr   


print(control_rad)

#-------------------------- Plots --------------------------#

#--- Absorbed Solar Radiation ---#
fig, ax = plt.subplots()
ax.plot( lat, net_ASR, label='Net ASR', color = 'black', linestyle='-')
ax.plot( lat, net_ASRcld, label='ASR With Clouds', color = 'blue', linestyle='--' )
ax.plot( lat, net_ASRclr, label='ASR Clear Sky', color = 'black', linestyle=':' )
ax.set_ylabel('Radiation W/m^2')
ax.set_xlabel('Latitude')
ax.set_title ('Absorbed Solar Radiation')
ax.legend(loc='upper right');
plt.savefig( location + '/Images/RRTGM/' + "ASR.svg", format="svg", bbox_inches='tight')
plt.show()

#--- Outgoing Longwave Radiation ---#
fig, ax = plt.subplots()
ax.plot( lat, net_OLR, label='Net OLR', color = 'black', linestyle='-')
ax.plot( lat, net_OLRcld, label='OLR With Clouds', color = 'blue', linestyle='--' )
ax.plot( lat, net_OLRclr, label='OLR Clear Sky', color = 'black', linestyle=':' )
ax.set_ylabel('Radiation W/m^2')
ax.set_xlabel('Latitude')
ax.set_title ('Outgoing Longwave Radiation')
ax.legend(loc='upper right');
plt.savefig( location + '/Images/RRTGM/' + "OLR.svg", format="svg", bbox_inches='tight')
plt.show()

#--- SW Down - Clouds ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt[:15], np.transpose(sw_down[:,1:16]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('SW Down - Clouds W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "sw_down_cld.svg", format="svg", bbox_inches='tight')
plt.show()

#--- SW Down - Clear Sky ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt[:15], np.transpose(sw_down_clr[:,1:16]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('SW Down - Clear Sky W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "sw_down_clr.svg", format="svg", bbox_inches='tight')
plt.show()

#--- SW Up - Clouds ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt[:15], np.transpose(sw_up[:,1:16]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('SW Up - Clouds W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "sw_up_cld.svg", format="svg", bbox_inches='tight')
plt.show()

#--- SW Up - Clear Sky ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt[:20], np.transpose(sw_up_clr[:,1:21]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('SW Up - Clear Sky W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "sw_up_clr.svg", format="svg", bbox_inches='tight')
plt.show()

#--- SW - Net ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt[:15], np.transpose(sw_net[:,1:16]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('SW - Net W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "sw_net.svg", format="svg", bbox_inches='tight')
plt.show()

#--- LW Down - Clouds ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt, np.transpose(lw_down[:,1:]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('LW Down - Clouds W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "lw_down_cld.svg", format="svg", bbox_inches='tight')
plt.show()

#--- LW Down - Clear Sky ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt, np.transpose(lw_down_clr[:,1:]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('LW Down - Clear Sky W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "lw_down_clr.svg", format="svg", bbox_inches='tight')
plt.show()

#--- LW Up - Clouds ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt, np.transpose(lw_up[:,1:]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('LW Up - Clouds W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "lw_up_cld.svg", format="svg", bbox_inches='tight')
plt.show()

#--- LW Up - Clear Sky ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt, np.transpose(lw_up_clr[:,1:]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('LW Up - Clear Sky W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "lw_up_clr.svg", format="svg", bbox_inches='tight')
plt.show()

#--- LW - Net ---#
fig, ax = plt.subplots()
cont = ax.contourf( lat, alt, np.transpose(lw_net[:,1:]), cmap='coolwarm' )
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('LW - Net W/m^2')
plt.savefig( location + '/Images/RRTGM/' + "lw_net.svg", format="svg", bbox_inches='tight')
plt.show()
