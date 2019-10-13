# -*- coding: utf-8 -*-
"""

@author: Tristan O'Hanlon

"""
import h5py
import os
import matplotlib.pyplot as plt
import constants
import numpy as np
import cartopy.crs as ccrs

#---Importing Data from Reduced Datasets---#

class Model:
    def __init__( self, path, name ):
        self.name = name
        a = h5py.File(path, 'r')
       
        #---alt-lat contour data---#

        self.cl_alt_lat = a['cl_alt_lat'][:] # total cloud fraction corresponding to alt and lat
        self.clw_alt_lat = a['clw_alt_lat'][:] # cloud liquid water fraction corresponding to liq_alt and lat
        self.cli_alt_lat = a['cli_alt_lat'][:] # cloud ice water fraction corresponding to alt and lat
        self.ta_alt_lat = a['ta_alt_lat'][:] # temperature corresponding to liq_alt and lat
        self.full_clw_alt_lat = a['full_clw_alt_lat'][:] # cloud liquid water fraction corresponding to alt and lat
        self.full_ta_alt_lat = a['full_ta_alt_lat'][:] # temperature corresponding to alt and lat

        #---global profile---#

        self.cl_g = a['cl_g'][:] # global layer total cloud fraction corresponding to alt
        self.clw_g = a['clw_g'][:] # global layer cloud liquid water fraction corresponding to liq_alt
        self.cli_g = a['cli_g'][:] # global layer cloud ice water fraction corresponding to alt

        #---southern ocean profile---#

        self.cl_so = a['cl_so'][:] # southern ocean layer total cloud fraction corresponding to alt
        self.clw_so = a['clw_so'][:] # southern ocean layer cloud liquid water fraction corresponding to liq_alt
        self.cli_so = a['cli_so'][:] # southern ocean layer cloud ice water fraction corresponding to alt
        
        #---clould liquid water fraction maped to temperature---#
        
        self.clw_t_g = a['clw_t_g'][:] # corresponding to ta_g
        self.clw_t_so = a['clw_t_so'][:] # corresponding to ta_so
        
        #---zonal data---#
            
        self.clt = a['clt'][:] # total cloud fraction corresponding to lat
        self.clwvi = a['clwvi'][:] # total cloud liquid water fraction corresponding to lat
        self.clivi = a['clivi'][:] # total cloud ice water fraction corresponding to lat

        #---regional data---#        
        
        self.clt_lat_lon = a['clt_lat_lon'][:] # total cloud fraction corresponding to lat and lon
        self.clwvi_lat_lon = a['clwvi_lat_lon'][:] # total cloud liquid water fraction corresponding to lat and lon
        self.clivi_lat_lon = a['clivi_lat_lon'][:] # total cloud ice water fraction corresponding to lat and lon
       
        
class MISR_SAT:
    def __init__( self, path, name ):
        self.name = name
        a = h5py.File(path, 'r')
       
        #---global profile---#
        self.cl_g = a['cl_g'][:] # global layer total cloud fraction corresponding to alt

        #---southern ocean profile---#
        self.cl_so = a['cl_so'][:] # southern ocean layer total cloud fraction corresponding to alt


class CERES_SAT:
    def __init__( self, path, name ):
        self.name = name
        a = h5py.File(path, 'r')
       
        #---zonal data---#
            
        self.clt = a['clt'][:] # total cloud fraction corresponding to lat
        self.clwvi = a['clw'][:] # total cloud liquid water fraction corresponding to lat
        self.clivi = a['cli'][:] # total cloud ice water fraction corresponding to lat

        #---regional data---#        
        
#        self.clt_lat_lon = a['clt_lat_lon'][:] # total cloud fraction corresponding to lat and lon
#        self.clwvi_lat_lon = a['clwvi_lat_lon'][:] # total cloud liquid water fraction corresponding to lat and lon
#        self.clivi_lat_lon = a['clivi_lat_lon'][:] # total cloud ice water fraction corresponding to lat and lon
       
        #---regional radiation data---#        
        
        self.clr_toa_sw_reg = a['clr_toa_sw_reg'][:] # corresponding to lat and lon
        self.clr_toa_sw_reg = a['clr_toa_sw_reg'][:] # corresponding to lat and lon
        self.clr_toa_lw_reg = a['clr_toa_lw_reg'][:] # corresponding to lat and lon
        self.clr_toa_lw_reg = a['clr_toa_lw_reg'][:] # corresponding to lat and lon
        self.clr_toa_net_reg = a['clr_toa_net_reg'][:] # corresponding to lat and lon
        self.clr_toa_net_reg = a['clr_toa_net_reg'][:] # corresponding to lat and lon
        self.all_toa_sw_reg = a['all_toa_sw_reg'][:] # corresponding to lat and lon
        self.all_toa_sw_reg = a['all_toa_sw_reg'][:] # corresponding to lat and lon
        self.all_toa_lw_reg = a['all_toa_lw_reg'][:] # corresponding to lat and lon
        self.all_toa_lw_reg_so = a['all_toa_lw_reg_so'][:] # corresponding to lat and lon
        self.all_toa_net_reg = a['all_toa_net_reg'][:] # corresponding to lat and lon
        self.all_toa_net_reg_so = a['all_toa_net_reg_so'][:] # corresponding to lat and lon

class CCCM_SAT:
    def __init__( self, path, name ):
        self.name = name
        a = h5py.File(path, 'r')
       
        #---alt-lat contour data---#

        self.cl_alt_lat = a['cl_alt_lat'][:] # total cloud fraction corresponding to alt and lat
        self.clw_alt_lat = a['clw_alt_lat'][:] # cloud liquid water fraction corresponding to liq_alt and lat
        self.cli_alt_lat = a['cli_alt_lat'][:] # cloud ice water fraction corresponding to alt and lat
        self.ta_alt_lat = a['ta_alt_lat'][:] # temperature corresponding to liq_alt and lat
        self.full_clw_alt_lat = a['full_clw_alt_lat'][:] # cloud liquid water fraction corresponding to alt and lat
        self.full_ta_alt_lat = a['full_ta_alt_lat'][:] # temperature corresponding to alt and lat

        #---global profile---#

        self.cl_g = a['cl_g'][:] # global layer total cloud fraction corresponding to alt
        self.clw_g = a['clw_g'][:] # global layer cloud liquid water fraction corresponding to liq_alt
        self.cli_g = a['cli_g'][:] # global layer cloud ice water fraction corresponding to alt

        #---southern ocean profile---#

        self.cl_so = a['cl_so'][:] # southern ocean layer total cloud fraction corresponding to alt
        self.clw_so = a['clw_so'][:] # southern ocean layer cloud liquid water fraction corresponding to liq_alt
        self.cli_so = a['cli_so'][:] # southern ocean layer cloud ice water fraction corresponding to alt
        
        #---clould liquid water fraction maped to temperature---#
        
        self.clw_t_g = a['clw_t_g'][:] # corresponding to ta_g
        self.clw_t_so = a['clw_t_so'][:] # corresponding to ta_so
        
#        #---zonal data---#
#            
#        self.clt = a['clt'][:] # total cloud fraction corresponding to lat
#        self.clwvi = a['clwvi'][:] # total cloud liquid water fraction corresponding to lat
#        self.clivi = a['clivi'][:] # total cloud ice water fraction corresponding to lat
#
#        #---regional data---#        
#        
#        self.clt_lat_lon = a['clt_lat_lon'][:] # total cloud fraction corresponding to lat and lon
#        self.clwvi_lat_lon = a['clwvi_lat_lon'][:] # total cloud liquid water fraction corresponding to lat and lon
#        self.clivi_lat_lon = a['clivi_lat_lon'][:] # total cloud ice water fraction corresponding to lat and lon


date_cmip5 = 'Jan_2001_Dec_2005'
date_cmip6 = 'Jan_2006_Dec_2010'
date_ceres = 'Jun_2006_Jun_2011'
date_cccm = 'Jun_2006_Apr_2011'


models = {
    'CMIP5-CESM1-CAM5' : Model( location + 'climate-analysis/reduced_data/' + date_cmip5 + '_CMIP5-CESM1-CAM5.h5', 'CMIP5-CESM1-CAM5' ),
    'CMIP5-GFDL-HIRAM-C360' : Model( location + 'climate-analysis/reduced_data/' + date_cmip5 + '_CMIP5-GFDL-HIRAM-C360.h5', 'CMIP5-GFDL-HIRAM-C360' ),
    'CMIP5-GISS-E2R' : Model( location + 'climate-analysis/reduced_data/' + date_cmip5 + '_CMIP5-GISS-E2R.h5', 'CMIP5-GISS-E2R' ),
    'CMIP5-IPSL-CM5A-LR' : Model( location + 'climate-analysis/reduced_data/' + date_cmip5 + '_CMIP5-IPSL-CM5A-LR.h5', 'CMIP5-IPSL-CM5A-LR' ),
    'CMIP5-MIROC5' : Model( location + 'climate-analysis/reduced_data/' + date_cmip5 + '_CMIP5-MIROC5.h5', 'CMIP5-MIROC5' ),
    'CMIP5-MRI-CGCM3' : Model( location + 'climate-analysis/reduced_data/' + date_cmip5 + '_CMIP5-MRI-CGCM3.h5', 'CMIP5-MRI-CGCM3' ),
    
    'CMIP6-CESM2-CAM6' : Model( location + 'climate-analysis/reduced_data/' + date_cmip6 + '_CMIP6-CESM2-CAM6.h5', 'CMIP6-CESM2-CAM6' ),
    'CMIP6-GFDL-AM4' : Model( location + 'climate-analysis/reduced_data/' + date_cmip6 + '_CMIP6-GFDL-AM4.h5', 'CMIP6-GFDL-AM4' ),
    'CMIP6-GISS-E21G' : Model( location + 'climate-analysis/reduced_data/' + date_cmip6 + '_CMIP6-GISS-E21G.h5', 'CMIP6-GISS-E21G' ),
    'CMIP6-IPSL-CM6A-LR' : Model( location + 'climate-analysis/reduced_data/' + date_cmip6 + '_CMIP6-IPSL-CM6A-LR.h5', 'CMIP6-IPSL-CM6A-LR' ),
    'CMIP6-MIROC6' : Model( location + 'climate-analysis/reduced_data/' + date_cmip6 + '_CMIP6-MIROC6.h5', 'CMIP6-MIROC6' ),
    'CMIP6-MRI-ESM2' : Model( location + 'climate-analysis/reduced_data/' + date_cmip6 + '_CMIP6-MRI-ESM2.h5', 'CMIP6-MRI-ESM2' ),

    'ECMWF' : Model( location + 'climate-analysis/reduced_data/' + date_cmip6 + '_ECMWF.h5', 'ECMWF-ERA5' ),
    'CALIPSO' : Model( location + 'climate-analysis/reduced_data/' + date_ceres + '_CALIPSO.h5', 'CALIPSO-GOCCP' ),

    'MISR' : MISR_SAT( location + 'climate-analysis/reduced_data/' + date_cmip6 + '_MISR.h5', 'MISR - Satellite' ),
    'CCCM' : CCCM_SAT( location + 'climate-analysis/reduced_data/' + date_cccm + '_CCCM.h5', 'CCCM - Satellite' ),
    'CERES' : CERES_SAT( location + 'climate-analysis/reduced_data/' + date_ceres + '_CERES.h5', 'CERES - Satellite' ),
}


quantity = {
    'clt' : 'Cloud Fraction',
    'clwvi' : 'Cloud Liquid Water Fraction',
    'clivi' : 'Cloud ice Water Fraction',
    'cl' : 'Cloud Fraction',
    'clw' : 'Cloud Liquid Water Fraction',
    'cli' : 'Cloud ice Water Fraction'
}



#This needs to be at least as long as the number of models you will ever graph at the same time  
colours = ['-k', '-b', '-r', '-g', '-m', '-c', '--k', '--b', '--r', '--g', '--m', '--c', ':k', ':b', ':r', ':g', ':m', ':c' ]

############################################################################### global latitude plots

# plot zonal quantities: clt, clwvi or clivi with latitude

def latitude_plot_models( models, quantity, cmip5_range = True, cmip6_range = True ):
    
    # max and min range for all CMIP5 models - be able to turn this on or off
    cmip5_max = np.maximum.reduce( [ ( all CMIP5 models).quantity ] )
    cmip5_min = np.minimum.reduce( [ ( all CMIP5 models).quantity ] )
    
    # max and min range for all CMIP6 models - be able to turn this on or off
    cmip6_max = np.maximum.reduce( [ ( all CMIP6 models).quantity ] )
    cmip6_min = np.minimum.reduce( [ ( all CMIP6 models).quantity ] )

    fig, (ax) = plt.subplots()
    for model, col in zip( models, colours ):
        ax.plot( constants.lat, model.quantity, col, label=model.name )

    # plot the range limits from all models and fill between them.
    if cmip5_range = True:
        ax.fill_between( lat, cmip5_min, cmip5_max, facecolor='black', alpha=0.2,
                    label='CMIP5 Model Range' )
        
    if cmip6_range = True:
        ax.fill_between( lat, cmip6_min, cmip6_max, facecolor='yellow', alpha=0.4,
                    label='CMIP6 Model Range' )

    plt.grid( True )
    ax.set_ylabel( quantity.name )
    ax.set_xlabel( 'Latitude' )
    ax.set_title ( quantity.name + ' vs Latitude' )

    plt.savefig( location + '/Images/' + quantity + '_' + models.name + ".pdf", format="pdf", bbox_inches='tight' )
    plt.show()

location = constants.home
latitude_plot_models( [ models['ECMWF'], models[ 'GFDL-HIRAM-C360'] ], clt, cmip5_range = True, cmip6_range = False )
latitude_plot_models( models.values() )



############################################################################### altitude plots

# plot zonal quantities cl, clw or cli with altitude - 2 graphs: top = global, bottom = southern ocean
# cl and cli go with alt
# clw goes with liq_alt
# global = g, southern ocean = so

def altitude_plot_models( models, quantity, cmip5_range = True, cmip6_range = True ):

    # max and min range for all CMIP5 models - be able to turn this on or off
    cmip5_max = np.maximum.reduce( [ ( all CMIP5 models).quantity ] )
    cmip5_min = np.minimum.reduce( [ ( all CMIP5 models).quantity ] )

    cmip5_so_max = np.maximum.reduce( [ ( all CMIP5 models).quantity_so ] )
    cmip5_so_min = np.minimum.reduce( [ ( all CMIP5 models).quantity_so ] )
    
    # max and min range for all CMIP6 models - be able to turn this on or off
    cmip6_max = np.maximum.reduce( [ ( all CMIP6 models).quantity ] )
    cmip6_min = np.minimum.reduce( [ ( all CMIP6 models).quantity ] )

    cmip6_so_max = np.maximum.reduce( [ ( all CMIP6 models).quantity_so ] )
    cmip6_so_min = np.minimum.reduce( [ ( all CMIP6 models).quantity_so ] )
    
    
    fig, ( ax1, ax2 ) = plt.subplots()
    for model, col in zip( models, colours ):
        if quantity == 'clw':
            ax1.plot( model.quantity + '_g', constants.liq_alt, col, label=model.name )
        else:
            ax1.plot( model.quantity + '_g', constants.alt, col, label=model.name )
 
    for model, col in zip( models, colours ):
        if quantity == 'clw':
            ax2.plot( model.quantity + '_so', constants.liq_alt, col, label=model.name )
        else:
            ax2.plot( model.quantity + '_so', constants.alt, col, label=model.name )
        

    # plot the range limits from all models on both axes and fill between them.
    if cmip5_range = True:
        if quantity == 'clw':
            ax1.fill_betweenx( constants.liq_alt, cmip5_min, cmip5_max, facecolor='black', alpha=0.2,
                        label='CMIP5 Model Range' )
            ax2.fill_betweenx( constants.liq_alt, cmip5_so_min, cmip5_so_max, facecolor='black', alpha=0.2,
                        label='CMIP5 Model Range' )
        else:
            ax1.fill_betweenx( constants.alt, cmip5_min, cmip5_max, facecolor='black', alpha=0.2,
                        label='CMIP5 Model Range' )
            ax2.fill_betweenx( constants.alt, cmip5_so_min, cmip5_so_max, facecolor='black', alpha=0.2,
                        label='CMIP5 Model Range' )
           
    if cmip6_range = True:   
        if quantity == 'clw':
            ax1.fill_betweenx( constants.liq_alt, cmip6_min, cmip6_max, facecolor='yellow', alpha=0.4,
                        label='CMIP6 Model Range' )    
            ax2.fill_betweenx( constants.liq_alt, cmip6_so_min, cmip6_so_max, facecolor='yellow', alpha=0.4,
                        label='CMIP6 Model Range' )
        else:
            ax1.fill_betweenx( constants.alt, cmip6_min, cmip6_max, facecolor='yellow', alpha=0.4,
                        label='CMIP6 Model Range' )    
            ax2.fill_betweenx( constants.alt, cmip6_so_min, cmip6_so_max, facecolor='yellow', alpha=0.4,
                        label='CMIP6 Model Range' )


    ax.set_ylabel( quantity.name )
    ax.set_xlabel( 'Altitude (km)' )
    ax.set_title ( quantity.name + ' vs Altitude' )

    ax1.set_xlim( 0, 0.35 )
    ax2.set_xlim( 0, 0.35 )
    
    ax1.grid( True )
    ax2.grid( True )

    plt.savefig( location + '/Images/' + quantity + models + ".pdf", format="pdf", bbox_inches='tight' )
    plt.show()

location = constants.home
altitude_plot_models( [ models['ECMWF'], models[ 'GFDL-HIRAM-C360'] ], clw, cmip5_range = True, cmip6_range = True )
altitude_plot_models( models.values() )



############################################################################### temperature plots

# plot clw_t with temperature - 2 graphs: top = global, bottom = southern ocean
# global = g, southern ocean = so

def liq_temp_plot_models( models, cmip5_range = True, cmip6_range = True ):

    # max and min range for all CMIP5 models - be able to turn this on or off
    cmip5_max = np.maximum.reduce( [ ( all CMIP5 models).clw_t_g ] )
    cmip5_min = np.minimum.reduce( [ ( all CMIP5 models).clw_t_g ] )

    cmip5_so_max = np.maximum.reduce( [ ( all CMIP5 models).clw_t_so ] )
    cmip5_so_min = np.minimum.reduce( [ ( all CMIP5 models).clw_t_so ] )
    
    # max and min range for all CMIP6 models - be able to turn this on or off
    cmip6_max = np.maximum.reduce( [ ( all CMIP6 models).clw_t_g ] )
    cmip6_min = np.minimum.reduce( [ ( all CMIP6 models).clw_t_g ] )

    cmip6_so_max = np.maximum.reduce( [ ( all CMIP6 models).clw_t_so ] )
    cmip6_so_min = np.minimum.reduce( [ ( all CMIP6 models).clw_t_so ] )
    
    
    fig, ( ax1, ax2 ) = plt.subplots()
    for model, col in zip( models, colours ):
        ax1.plot( constants.ta_g, 'clw_t_g',  col, label=model.name )
 
    for model, col in zip( models, colours ):
        ax2.plot( constants.ta_so, 'clw_t_so', col, label=model.name )
        

    # plot the range limits from all models on both axes and fill between them.
    if cmip5_range = True:
        ax1.fill_between( constants.ta_g, cmip5_min, cmip5_max, facecolor='black', alpha=0.2,
                    label='CMIP5 Model Range' )
        ax2.fill_between( constants.ta_so, cmip5_so_min, cmip5_so_max, facecolor='black', alpha=0.2,
                    label='CMIP5 Model Range' )
           
    if cmip6_range = True:   
        ax1.fill_between( constants.ta_g, cmip6_min, cmip6_max, facecolor='yellow', alpha=0.4,
                    label='CMIP6 Model Range' )    
        ax2.fill_between( constants.ta_so, cmip6_so_min, cmip6_so_max, facecolor='yellow', alpha=0.4,
                    label='CMIP6 Model Range' )


    ax.set_ylabel( 'Cloud Liquid Water Fraction' )
    ax.set_xlabel( 'Temperature (C)' )
    ax.set_title ( 'Cloud Liquid Water Fraction vs Temperature (C)' )
   
    ax1.grid( True )
    ax2.grid( True )

    plt.savefig( location + '/Images/' + 'liq_temp' + '_' +  models + ".pdf", format="pdf", bbox_inches='tight' )
    plt.show()

location = constants.home
liq_temp_plot_models( [ models['CALIPSO'], models[ 'CMIP6-GFDL-AM4'] ], cmip5_range = True, cmip6_range = True )
liq_temp_plot_models( models.values() )



############################################################################### contour plots

def satellite_contours( models, zero_line ):
    
    fig, ax = plt.subplots( nrows=1, ncols=models.size )
    count = 0
    for model in models:
        if count = 0:
            prime = ax[ 0, count ].contourf( constants.lat, constants.liq_alt, model.clw_alt_lat )
        else:
            ax[ 0, count ].contourf( constants.lat, constants.liq_alt, model.clw_alt_lat )

        temp = ax[ 0, count ].contour( constants.lat, constants.liq_alt, ( model.ta_alt_lat - 273.15 ), colors='white' )
        temp.collections[ zero_line ].set_linewidth( 3 )
        temp.collections[ zero_line ].set_color( 'white' )
        ax.clabel( temp, inline=1, fontsize=10 )
        ax[ 0, count ].set_xlabel( 'Latitude' )
        ax[ 0, count ].set_title( model.name )
        count += 1

    # the common color bar is based off the model labeled 'prime'
    ax[ 0, 0 ].set_ylabel('Altitude (km)')
    cbaxes = fig.add_axes([0.7, 0.5, 0.03, 0.3]) #(x-position, y-position, thickness, length)
    cbar = fig.colorbar(prime, cax=cbaxes)
    cbar.set_label('Cloud liquid Water Fraction')    
    cbar.set_clim(0, 0.5)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)
    plt.savefig( location + '/Images/' + 'contours' + '_' +  models + ".pdf", format="pdf", bbox_inches='tight')
    plt.show()
   
location = constants.home
satellite_contours( [ models[ 'CALIPSO' ], 5, models[ 'CCCM' ], 6, models[ 'ECMWF' ], 6 ] )


def model_pair_contours( models, zero_line ):
    
    fig, ax = plt.subplots( nrows=1, ncols=models.size )
    count = 0
    for model in models:
        if count = 0:
            prime = ax[ 0, count ].contourf( constants.lat, constants.liq_alt, model.clw_alt_lat )
        else:
            ax[ 0, count ].contourf( constants.lat, constants.liq_alt, model.clw_alt_lat )

        temp = ax[ 0, count ].contour( constants.lat, constants.liq_alt, ( model.ta_alt_lat - 273.15 ), colors='white' )
        temp.collections[ zero_line ].set_linewidth( 3 )
        temp.collections[ zero_line ].set_color( 'white' )
        ax.clabel( temp, inline=1, fontsize=10 )
        ax[ 0, count ].set_xlabel( 'Latitude' )
        ax[ 0, count ].set_title( model.name )
        count += 1

    # the common color bar is based off the model labeled 'prime'
    ax[ 0, 0 ].set_ylabel('Altitude (km)')
    cbaxes = fig.add_axes([0.7, 0.5, 0.03, 0.3]) #(x-position, y-position, thickness, length)
    cbar = fig.colorbar(prime, cax=cbaxes, location=bottom, orientation='horizontal' )
    cbar.set_label('Cloud liquid Water Fraction')    
    cbar.set_clim(0, 0.5)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)
    plt.savefig( location + '/Images/' + 'contours' + '_' +  models + ".pdf", format="pdf", bbox_inches='tight')
    plt.show()


location = constants.home
model_pair_contours( [ models[ 'CMIP5-CESM1-CAM5' ], 5, models[ 'CMIP6-CESM2-CAM6' ], 6 )


############################################################################### regional plots

def region_plot( models, quantity ):
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.coastlines()
    count = 0
    for model in models:
        if count = 0:
            prime = ax.contourf(constants.lon, constants.lat, model.quantity_lat_lon, transform=ccrs.PlateCarree())
        else:
            ax.contourf(constants.lon, constants.lat, model.quantity_lat_lon, transform=ccrs.PlateCarree())
        count += 1

    cbar = plt.colorbar(prime, orientation='horizontal')
    cbar.set_label('Cloud Fraction')
    plt.show()


