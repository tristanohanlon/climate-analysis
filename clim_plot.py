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
location = constants.home

class Model:
    def __init__( self, path, name ):
        self.name = name
        a = h5py.File(path, 'r')

        self.data = {}
        for e in a.keys():
            self.data[ e ] = a[ e ][:]
        


models = {
    'CMIP5-AMIP-CESM1-CAM5' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip5 + '_CMIP5-AMIP-CESM1-CAM5.h5', 'CMIP5-AMIP-CESM1-CAM5' ),
    'CMIP5-AMIP-GFDL-CM3' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip5 + '_CMIP5-AMIP-GFDL-CM3.h5', 'CMIP5-AMIP-GFDL-CM3' ),
    'CMIP5-AMIP-MIROC5' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip5 + '_CMIP5-AMIP-MIROC5.h5', 'CMIP5-AMIP-MIROC5' ),
    'CMIP5-AMIP-MRI-CGCM3' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip5 + '_CMIP5-AMIP-MRI-CGCM3.h5', 'CMIP5-AMIP-MRI-CGCM3' ),
    'CMIP5-AMIP-GISS-E2R' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip5 + '_CMIP5-AMIP-GISS-E2R.h5', 'CMIP5-AMIP-GISS-E2R' ),
    'CMIP5-AMIP-IPSL-CM5A-LR' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip5 + '_CMIP5-AMIP-IPSL-CM5A-LR.h5', 'CMIP5-AMIP-IPSL-CM5A-LR' ),
    
    'CMIP6-AMIP-CESM2-CAM6' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip6 + '_CMIP6-AMIP-CESM2-CAM6.h5', 'CMIP6-AMIP-CESM2-CAM6' ),
    'CMIP6-AMIP-GFDL-CM4' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip6 + '_CMIP6-AMIP-GFDL-CM4.h5', 'CMIP6-AMIP-GFDL-CM4' ),
    'CMIP6-AMIP-MIROC6' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip6 + '_CMIP6-AMIP-MIROC6.h5', 'CMIP6-AMIP-MIROC6' ),
    'CMIP6-AMIP-MRI-ESM2' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip6 + '_CMIP6-AMIP-MRI-ESM2.h5', 'CMIP6-AMIP-MRI-ESM2' ),
    'CMIP6-AMIP-GISS-E21G' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip6 + '_CMIP6-AMIP-GISS-E21G.h5', 'CMIP6-AMIP-GISS-E21G' ),
    'CMIP6-AMIP-IPSL-CM6A-LR' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip6 + '_CMIP6-AMIP-IPSL-CM6A-LR.h5', 'CMIP6-AMIP-IPSL-CM6A-LR' ),

    # 'ECMWF' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cmip6 + '_ECMWF.h5', 'ECMWF-ERA5' ),
    'CALIPSO' : Model( location + 'climate-analysis/reduced_data/' + constants.date_ceres + '_CALIPSO.h5', 'CALIPSO-GOCCP' ),

    'CCCM' : Model( location + 'climate-analysis/reduced_data/' + constants.date_cccm + '_CCCM.h5', 'CCCM' ),
    'CERES' : Model( location + 'climate-analysis/reduced_data/' + constants.date_ceres + '_CERES.h5', 'CERES' ),
}

cosp_models = {
    'CMIP5-AMIP-CESM1-CAM5' : Model( location + 'climate-analysis/reduced_data/' + 'COSP_' + constants.date_cmip5 + '_CMIP5-AMIP-CESM1-CAM5.h5', 'CMIP5-AMIP-CESM1-CAM5-COSP' ),
    'CMIP5-AMIP-GFDL-CM3' : Model( location + 'climate-analysis/reduced_data/' + 'COSP_' + constants.date_cmip5 + '_CMIP5-AMIP-GFDL-CM3.h5', 'CMIP5-AMIP-GFDL-CM3-COSP' ),
    'CMIP5-AMIP-MIROC5' : Model( location + 'climate-analysis/reduced_data/' + 'COSP_' + constants.date_cmip5 + '_CMIP5-AMIP-MIROC5.h5', 'CMIP5-AMIP-MIROC5-COSP' ),
    'CMIP5-AMIP-MRI-CGCM3' : Model( location + 'climate-analysis/reduced_data/' + 'COSP_' + constants.date_cmip5 + '_CMIP5-AMIP-MRI-CGCM3.h5', 'CMIP5-AMIP-MRI-CGCM3-COSP' ),
    
    'CMIP6-AMIP-CESM2-CAM6' : Model( location + 'climate-analysis/reduced_data/' + 'COSP_' + constants.date_cmip6 + '_CMIP6-AMIP-CESM2-CAM6.h5', 'CMIP6-AMIP-CESM2-CAM6-COSP' ),
    'CMIP6-AMIP-GFDL-CM4' : Model( location + 'climate-analysis/reduced_data/' + 'COSP_' + constants.date_cmip6 + '_CMIP6-AMIP-GFDL-CM4.h5', 'CMIP6-AMIP-GFDL-CM4-COSP' ),
    'CMIP6-AMIP-MIROC6' : Model( location + 'climate-analysis/reduced_data/' + 'COSP_' + constants.date_cmip6 + '_CMIP6-AMIP-MIROC6.h5', 'CMIP6-AMIP-MIROC6-COSP' ),
    'CMIP6-AMIP-MRI-ESM2' : Model( location + 'climate-analysis/reduced_data/' + 'COSP_' + constants.date_cmip6 + '_CMIP6-AMIP-MRI-ESM2.h5', 'CMIP6-AMIP-MRI-ESM2-COSP' ),
}




quantity = {
    'clt' : 'Cloud Fraction',
    'clt_l' : 'Low Cloud Fraction',
    'clt_m' : 'Medium Cloud Fraction',
    'clt_h' : 'High Cloud Fraction',
   
    'clwvi' : 'Mean Cloud Liquid Water Path (kg/m^2)',
    'clivi' : 'Mean Cloud Ice Water Path (kg/m^2)',
    'cl' : 'Cloud Fraction',
    'clw' : 'Mean Cloud Liquid Water Mass Fraction in Air (g/kg)',
    'clwc' : 'Mean Cloud Liquid Water Content (g/m^3)',
    'clw_frac' : 'Mean Cloud Liquid Water Fraction',
    'cli_frac' : 'Mean Cloud Ice Water Fraction',    
    'cli' : 'Mean Cloud Ice Water Mass Fraction in Air (kg/kg)',
    'clt_lc' : 'Low Cloud Fraction',
    'clwvi_lc' : 'Mean Low Cloud Liquid Water Path (kg/m^2)',
    'clivi_lc' : 'Mean Low Cloud Ice Water Path (kg/m^2)',
    'loaddust' : 'The total dry mass of dust aerosol particles per unit area g/m^2',
    'loadss' : 'The total dry mass of sea salt aerosol particles per unit area g/m^2',
    'loadso4' : 'The total dry mass of sulphate aerosol particles per unit area g/m^2',
    'rsdt' : 'TOA Incoming Shortwave Flux W/m^2',
    'rsut' : 'TOA Outgoing Shortwave Flux W/m^2',
    'rsutcs' : 'TOA Outgoing Shortwave Flux Assuming Clear Sky W/m^2',
    'rtmt' : 'Net Downwards Radiative Flux at TOA W/m^2',
    'albedo_reg' : 'Albedo'
}



#This needs to be at least as long as the number of models you will ever graph at the same time  
colours = ['--b', '-b', '--r', '-r', '-k', '-g', '--k', '--b', '--r', '--g', '--m', '--c', ':k', ':b', ':r', ':g', ':m', ':c' ]
colours_cosp = ['-b', '-r', '-k', '-g']
############################################################################### global latitude plots

# plot zonal quantities: clt, clwvi or clivi with latitude

def latitude_plot_models( models_to_graph, quantity_to_graph, cmip5_range = True, cmip6_range = True ):
    

    fig, (ax) = plt.subplots()
    for model, col in zip( models_to_graph, colours ):
        ax.plot( constants.lat[constants.lat_confine_1:constants.lat_confine_2], model.data[quantity_to_graph][constants.lat_confine_1:constants.lat_confine_2], col, label=model.name )

    # plot the range limits from all models and fill between them.
    if cmip5_range == True:
        # max and min range for all CMIP5 models - be able to turn this on or off
        cmip5models = [ model.data[quantity_to_graph] for name, model in models.items() if name.startswith('CMIP5') ]
        cmip5_max = np.maximum.reduce( cmip5models )
        cmip5_min = np.minimum.reduce( cmip5models )

        ax.fill_between( constants.lat[constants.lat_confine_1:constants.lat_confine_2], cmip5_min[constants.lat_confine_1:constants.lat_confine_2], cmip5_max[constants.lat_confine_1:constants.lat_confine_2], facecolor='black', alpha=0.3,
                    label='CMIP5 Model Range' )
        
    if cmip6_range == True:
        cmip6models = [ model.data[quantity_to_graph] for name, model in models.items() if name.startswith('CMIP6') ]
        cmip6_max = np.maximum.reduce( cmip6models )
        cmip6_min = np.minimum.reduce( cmip6models )

        ax.fill_between( constants.lat[constants.lat_confine_1:constants.lat_confine_2], cmip6_min[constants.lat_confine_1:constants.lat_confine_2], cmip6_max[constants.lat_confine_1:constants.lat_confine_2], facecolor='green', alpha=0.3,
                    label='CMIP6 Model Range' )

    plt.grid( True )
    ax.legend(loc='upper center', bbox_to_anchor=(1.25, 1.0));
    ax.set_ylabel( quantity[quantity_to_graph] )
    ax.set_xlabel( 'Latitude' )
    ax.set_title ( quantity[quantity_to_graph] + ' vs Latitude' )

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])

    plt.savefig( location + '/Images/' + quantity_to_graph + '_' + all_model_names + ".svg", format="svg", bbox_inches='tight' )
    plt.show()



############################################################################### altitude plots

# plot zonal quantities cl, clw or cli with altitude - 2 graphs: top = global, bottom = southern ocean
# cl and cli go with alt
# clw goes with liq_alt
# global = g, southern ocean = so

def altitude_plot_models( models_to_graph, quantity_to_graph, cmip5_range = True, cmip6_range = True ):
 
    
    fig, ( ax1, ax2 ) = plt.subplots(1, 2)
    for model, col in zip( models_to_graph, colours ):
        if quantity_to_graph == 'clw'  or quantity_to_graph == 'clwc':
            ax1.plot( model.data[quantity_to_graph + '_g']*1000, constants.liq_alt, col, label=model.name )
        elif quantity_to_graph == 'clw_frac':
            ax1.plot( model.data[quantity_to_graph + '_g'], constants.liq_alt, col, label=model.name )
        else:
            ax1.plot( model.data[quantity_to_graph + '_g'], constants.alt, col, label=model.name )
 
    for model, col in zip( models_to_graph, colours ):
        if quantity_to_graph == 'clw' or quantity_to_graph == 'clwc':
            ax2.plot( model.data[quantity_to_graph + '_so']*1000, constants.liq_alt, col, label=model.name )
        elif quantity_to_graph == 'clw_frac':
            ax2.plot( model.data[quantity_to_graph + '_so'], constants.liq_alt, col, label=model.name )
        else:
            ax2.plot( model.data[quantity_to_graph + '_so'], constants.alt, col, label=model.name )
        

    # plot the range limits from all models on both axes and fill between them.
    if cmip5_range == True:
        cmip5models_g = [ model.data[quantity_to_graph + '_g'] for name, model in models.items() if name.startswith('CMIP5') ]
        cmip5models_so = [ model.data[quantity_to_graph + '_so'] for name, model in models.items() if name.startswith('CMIP5') ]
        cmip5_g_max = np.maximum.reduce( cmip5models_g )
        cmip5_g_min = np.minimum.reduce( cmip5models_g )
        cmip5_so_max = np.maximum.reduce( cmip5models_so )
        cmip5_so_min = np.minimum.reduce( cmip5models_so )

        if quantity_to_graph == 'clw' or quantity_to_graph == 'clwc':
            ax1.fill_betweenx( constants.liq_alt, cmip5_g_min*1000, cmip5_g_max*1000, facecolor='black', alpha=0.3,
                        label='CMIP5 Model Range' )
            ax2.fill_betweenx( constants.liq_alt, cmip5_so_min*1000, cmip5_so_max*1000, facecolor='black', alpha=0.3,
                        label='CMIP5 Model Range' )
        else:
            ax1.fill_betweenx( constants.alt, cmip5_g_min, cmip5_g_max, facecolor='black', alpha=0.3,
                        label='CMIP5 Model Range' )
            ax2.fill_betweenx( constants.alt, cmip5_so_min, cmip5_so_max, facecolor='black', alpha=0.3,
                        label='CMIP5 Model Range' )
           
    if cmip6_range == True:
        cmip6models_g = [ model.data[quantity_to_graph + '_g'] for name, model in models.items() if name.startswith('CMIP6') ]
        cmip6models_so = [ model.data[quantity_to_graph + '_so'] for name, model in models.items() if name.startswith('CMIP6') ]
        cmip6_g_max = np.maximum.reduce( cmip6models_g )
        cmip6_g_min = np.minimum.reduce( cmip6models_g )
        cmip6_so_max = np.maximum.reduce( cmip6models_so )
        cmip6_so_min = np.minimum.reduce( cmip6models_so ) 
   
        if quantity_to_graph == 'clw' or quantity_to_graph == 'clwc':
            ax1.fill_betweenx( constants.liq_alt, cmip6_g_min*1000, cmip6_g_max*1000, facecolor='green', alpha=0.3,
                        label='CMIP6 Model Range' )    
            ax2.fill_betweenx( constants.liq_alt, cmip6_so_min*1000, cmip6_so_max*1000, facecolor='green', alpha=0.3,
                        label='CMIP6 Model Range' )
        else:
            ax1.fill_betweenx( constants.alt, cmip6_g_min, cmip6_g_max, facecolor='green', alpha=0.3,
                        label='CMIP6 Model Range' )    
            ax2.fill_betweenx( constants.alt, cmip6_so_min, cmip6_so_max, facecolor='green', alpha=0.3,
                        label='CMIP6 Model Range' )


    ax1.set_xlabel( quantity[quantity_to_graph] )
    ax2.set_xlabel( quantity[quantity_to_graph] )
    ax1.set_ylabel( 'Altitude (km)' )
    ax1.set_title ( 'Global' )
    ax2.set_title ( 'Southern Ocean' )
    ax2.legend(loc='upper center', bbox_to_anchor=(1.55, 1.0));

    ax1.axhline(y=3, label = 'Low Cloud Boundary', color = 'black', linestyle='--')
    ax2.axhline(y=3, label = 'Low Cloud Boundary', color = 'black', linestyle='--')

    # ax1.set_xlim( 0, 0.6 )
    # ax2.set_xlim( 0, 0.6 )
    
    ax1.grid( True )
    ax2.grid( True )

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])

    plt.savefig( location + '/Images/' + quantity_to_graph + '_' + all_model_names + ".svg", format="svg", bbox_inches='tight' )
    plt.show()


############################################################################### temperature plots

# plot cl_t with temperature - 2 graphs: top = global, bottom = southern ocean
# global = g, southern ocean = so

def cloud_temp_plot_models( models_to_graph, cmip5_range = True, cmip6_range = True ):

    cmip5models_g = [ model.data['cl_t_g'] for name, model in models.items() if name.startswith('CMIP5') ]
    cmip6models_g = [ model.data['cl_t_g'] for name, model in models.items() if name.startswith('CMIP6') ]
    cmip5models_so = [ model.data['cl_t_so'] for name, model in models.items() if name.startswith('CMIP5') ]
    cmip6models_so = [ model.data['cl_t_so'] for name, model in models.items() if name.startswith('CMIP6') ]
 
    cmip5_g_max = np.maximum.reduce( cmip5models_g )
    cmip5_g_min = np.minimum.reduce( cmip5models_g )
    cmip6_g_max = np.maximum.reduce( cmip6models_g )
    cmip6_g_min = np.minimum.reduce( cmip6models_g )

    cmip5_so_max = np.maximum.reduce( cmip5models_so )
    cmip5_so_min = np.minimum.reduce( cmip5models_so )
    cmip6_so_max = np.maximum.reduce( cmip6models_so )
    cmip6_so_min = np.minimum.reduce( cmip6models_so )
    
    fig, ( ax1, ax2 ) = plt.subplots(2,1, figsize=(15, 7))
    for model, col in zip( models_to_graph, colours ):
        ax1.plot( constants.ta, model.data['cl_t_g'],  col, label=model.name )
 
    for model, col in zip( models_to_graph, colours ):
        ax2.plot( constants.ta, model.data['cl_t_so'], col, label=model.name )
        

    # plot the range limits from all models on both axes and fill between them.
    if cmip5_range == True:
        ax1.fill_between( constants.ta, cmip5_g_min, cmip5_g_max, facecolor='red', alpha=0.3,
                    label='CMIP5 Model Range' )
        ax2.fill_between( constants.ta, cmip5_so_min, cmip5_so_max, facecolor='red', alpha=0.3,
                    label='CMIP5 Model Range' )
           
    if cmip6_range == True:   
        ax1.fill_between( constants.ta, cmip6_g_min, cmip6_g_max, facecolor='blue', alpha=0.2,
                    label='CMIP6 Model Range' )    
        ax2.fill_between( constants.ta, cmip6_so_min, cmip6_so_max, facecolor='blue', alpha=0.2,
                    label='CMIP6 Model Range' )


    ax1.set_ylabel( 'Cloud Fraction' )
    ax2.set_ylabel( 'Cloud Fraction' )
    ax2.set_xlabel( 'Temperature (K)' )
    ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
    ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

    ax1.set_title ( 'Global' )
    ax2.set_title ( 'Southern Ocean' )
    ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0));

    ax1.grid( True )
    ax2.grid( True )

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])
    fig.tight_layout()
    plt.savefig( location + '/Images/' + 'cl_t' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight' )
    plt.show()

#---------------------------------------------------------------------#

# plot clw_t with temperature - 2 graphs: top = global, bottom = southern ocean
# global = g, southern ocean = so

def liq_temp_plot_models( models_to_graph, cmip5_range = True, cmip6_range = True ):

    cmip5models_g = [ model.data['clw_t_g'] for name, model in models.items() if name.startswith('CMIP5') ]
    cmip6models_g = [ model.data['clw_t_g'] for name, model in models.items() if name.startswith('CMIP6') ]
    cmip5models_so = [ model.data['clw_t_so'] for name, model in models.items() if name.startswith('CMIP5') ]
    cmip6models_so = [ model.data['clw_t_so'] for name, model in models.items() if name.startswith('CMIP6') ]
 
    cmip5_g_max = np.maximum.reduce( cmip5models_g )*1000
    cmip5_g_min = np.minimum.reduce( cmip5models_g )*1000
    cmip6_g_max = np.maximum.reduce( cmip6models_g )*1000
    cmip6_g_min = np.minimum.reduce( cmip6models_g )*1000

    cmip5_so_max = np.maximum.reduce( cmip5models_so )*1000
    cmip5_so_min = np.minimum.reduce( cmip5models_so )*1000
    cmip6_so_max = np.maximum.reduce( cmip6models_so )*1000
    cmip6_so_min = np.minimum.reduce( cmip6models_so )*1000
    
    fig, ( ax1, ax2 ) = plt.subplots(2,1, figsize=(10, 8))
    for model, col in zip( models_to_graph, colours ):
        ax1.plot( constants.ta, model.data['clw_t_g']*1000,  col, label=model.name )
 
    for model, col in zip( models_to_graph, colours ):
        ax2.plot( constants.ta, model.data['clw_t_so']*1000, col, label=model.name )
        

    # plot the range limits from all models on both axes and fill between them.
    if cmip5_range == True:
        ax1.fill_between( constants.ta, cmip5_g_min, cmip5_g_max, facecolor='red', alpha=0.3,
                    label='CMIP5 Model Range' )
        ax2.fill_between( constants.ta, cmip5_so_min, cmip5_so_max, facecolor='red', alpha=0.3,
                    label='CMIP5 Model Range' )
           
    if cmip6_range == True:   
        ax1.fill_between( constants.ta, cmip6_g_min, cmip6_g_max, facecolor='blue', alpha=0.3,
                    label='CMIP6 Model Range' )    
        ax2.fill_between( constants.ta, cmip6_so_min, cmip6_so_max, facecolor='blue', alpha=0.3,
                    label='CMIP6 Model Range' )


    ax1.set_ylabel( quantity['clw'] )
    ax2.set_ylabel( quantity['clw'] )
    ax2.set_xlabel( 'Temperature (K)' )
    ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
    ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

    ax1.set_title ( 'Global' )
    ax2.set_title ( 'Southern Ocean' )
    ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0));

    ax1.grid( True )
    ax2.grid( True )

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])
    fig.tight_layout()
    plt.savefig( location + '/Images/' + 'liq_temp' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight' )
    plt.show()

#---------------------------------------------------------------------#


# plot clwc_t with temperature - 2 graphs: top = global, bottom = southern ocean
# global = g, southern ocean = so

def clw_frac_temp_plot_models( models_to_graph, quantity_to_graph, cmip5_range = True, cmip6_range = True ):

    fig, ( ax1, ax2 ) = plt.subplots(2,1, figsize=(15, 7))
    for model, col in zip( models_to_graph, colours ):
        ax1.plot( constants.ta, model.data[quantity_to_graph + '_t_g'],  col, label=model.name )
        ax2.plot( constants.ta, model.data[quantity_to_graph + '_t_so'], col, label=model.name )
        

    # plot the range limits from all models on both axes and fill between them.
    if cmip5_range == True:
        cmip5models_g = [ model.data[quantity_to_graph + '_t_g'] for name, model in models.items() if name.startswith('CMIP5') ]
        cmip5models_so = [ model.data[quantity_to_graph + '_t_so'] for name, model in models.items() if name.startswith('CMIP5') ]
        cmip5_g_max = np.maximum.reduce( cmip5models_g )
        cmip5_g_min = np.minimum.reduce( cmip5models_g )
        cmip5_so_max = np.maximum.reduce( cmip5models_so )
        cmip5_so_min = np.minimum.reduce( cmip5models_so )

        ax1.fill_between( constants.ta, cmip5_g_min, cmip5_g_max, facecolor='red', alpha=0.3,
                    label='CMIP5 Model Range' )
        ax2.fill_between( constants.ta, cmip5_so_min, cmip5_so_max, facecolor='red', alpha=0.3,
                    label='CMIP5 Model Range' )
           
    if cmip6_range == True:   
        cmip6models_g = [ model.data[quantity_to_graph + '_t_g'] for name, model in models.items() if name.startswith('CMIP6') ]
        cmip6models_so = [ model.data[quantity_to_graph + '_t_so'] for name, model in models.items() if name.startswith('CMIP6') ]
        cmip6_g_max = np.maximum.reduce( cmip6models_g )
        cmip6_g_min = np.minimum.reduce( cmip6models_g )
        cmip6_so_max = np.maximum.reduce( cmip6models_so )
        cmip6_so_min = np.minimum.reduce( cmip6models_so )

        ax1.fill_between( constants.ta, cmip6_g_min, cmip6_g_max, facecolor='blue', alpha=0.3,
                    label='CMIP6 Model Range' )    
        ax2.fill_between( constants.ta, cmip6_so_min, cmip6_so_max, facecolor='blue', alpha=0.3,
                    label='CMIP6 Model Range' )


    ax1.set_ylabel( 'Mean Cloud Liquid Water Fraction' )
    ax2.set_ylabel( 'Mean Cloud Liquid Water Fraction' )
    ax2.set_xlabel( 'Temperature (K)' )
    ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
    ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

    ax1.set_title ( 'Global' )
    ax2.set_title ( 'Southern Ocean' )
    ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0));

    ax1.grid( True )
    ax2.grid( True )

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])
    fig.tight_layout()
    plt.savefig( location + '/Images/' + 'clw_frac_t' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight' )
    plt.show()

#---------------------------------------------------------------------#


# plot clwc_t with temperature - 2 graphs: top = global, bottom = southern ocean
# global = g, southern ocean = so

def cli_frac_temp_plot_models( models_to_graph, quantity_to_graph, cmip5_range = True, cmip6_range = True ):

    fig, ( ax1, ax2 ) = plt.subplots(2,1, figsize=(15, 7))
    for model, col in zip( models_to_graph, colours ):
        ax1.plot( constants.ta, model.data[quantity_to_graph + '_t_g'],  col, label=model.name )
        ax2.plot( constants.ta, model.data[quantity_to_graph + '_t_so'], col, label=model.name )
        

    # plot the range limits from all models on both axes and fill between them.
    if cmip5_range == True:
        cmip5models_g = [ model.data[quantity_to_graph + '_t_g'] for name, model in models.items() if name.startswith('CMIP5') ]
        cmip5models_so = [ model.data[quantity_to_graph + '_t_so'] for name, model in models.items() if name.startswith('CMIP5') ]
        cmip5_g_max = np.maximum.reduce( cmip5models_g )
        cmip5_g_min = np.minimum.reduce( cmip5models_g )
        cmip5_so_max = np.maximum.reduce( cmip5models_so )
        cmip5_so_min = np.minimum.reduce( cmip5models_so )

        ax1.fill_between( constants.ta, cmip5_g_min, cmip5_g_max, facecolor='red', alpha=0.3,
                    label='CMIP5 Model Range' )
        ax2.fill_between( constants.ta, cmip5_so_min, cmip5_so_max, facecolor='red', alpha=0.3,
                    label='CMIP5 Model Range' )
           
    if cmip6_range == True:   
        cmip6models_g = [ model.data[quantity_to_graph + '_t_g'] for name, model in models.items() if name.startswith('CMIP6') ]
        cmip6models_so = [ model.data[quantity_to_graph + '_t_so'] for name, model in models.items() if name.startswith('CMIP6') ]
        cmip6_g_max = np.maximum.reduce( cmip6models_g )
        cmip6_g_min = np.minimum.reduce( cmip6models_g )
        cmip6_so_max = np.maximum.reduce( cmip6models_so )
        cmip6_so_min = np.minimum.reduce( cmip6models_so )

        ax1.fill_between( constants.ta, cmip6_g_min, cmip6_g_max, facecolor='blue', alpha=0.3,
                    label='CMIP6 Model Range' )    
        ax2.fill_between( constants.ta, cmip6_so_min, cmip6_so_max, facecolor='blue', alpha=0.3,
                    label='CMIP6 Model Range' )


    ax1.set_ylabel( 'Mean Cloud Ice Water Fraction' )
    ax2.set_ylabel( 'Mean Cloud Ice Water Fraction' )
    ax2.set_xlabel( 'Temperature (K)' )
    ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
    ax2.axvline(x=273, label = '273K', color = 'black', linestyle='--')

    ax1.set_title ( 'Global' )
    ax2.set_title ( 'Southern Ocean' )
    ax1.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0));

    ax1.grid( True )
    ax2.grid( True )

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])
    fig.tight_layout()
    plt.savefig( location + '/Images/' + 'ice_temp' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight' )
    plt.show()

############################################################################### contour plots

def satellite_contours( models_to_graph, zero_lines ):
    assert( len( models_to_graph ) == len( zero_lines ) )

    fig, ax = plt.subplots( )
    for count, (model, zero_line ) in enumerate( zip( models_to_graph, zero_lines ) ):
        prime = ax.contourf( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.alt[:constants.liq_alt_confine], model.data['full_clw_frac_alt_lat'][:constants.liq_alt_confine,constants.lat_confine_1:constants.lat_confine_2], vmin=0, vmax=0.25, cmap='coolwarm' )
        temp = ax.contour( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.alt[:constants.liq_alt_confine], ( model.data[ 'full_ta_alt_lat' ][:constants.liq_alt_confine,constants.lat_confine_1:constants.lat_confine_2] - 273.15 ), colors='white' )
        temp.collections[ zero_line ].set_linewidth( 3 )
        temp.collections[ zero_line ].set_color( 'white' )
        ax.clabel( temp, inline=1, fontsize=10 )
        ax.set_xlabel( 'Latitude' )
        ax.set_ylabel( 'Altitude (km)' )
        ax.set_title( model.name )

    # the common color bar is based off the model labeled 'prime'
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)

    # the common color bar is based off the model labeled 'prime'
    cbar = fig.colorbar(prime, orientation='horizontal', fraction = 0.05 )
    cbar.set_label('Mean Cloud Liquid Water Fraction')    
    cbar.set_clim(0, 0.25)

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])

    plt.savefig( location + '/Images/' + 'sat_contours' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight')
    plt.show()
   

#---------------------------------------------------------------------#


def model_clw_contours( models_to_graph, zero_lines ):
    assert( len( models_to_graph ) == len( zero_lines ) )

    fig, ax = plt.subplots( nrows=1, ncols=len( models_to_graph ), figsize=(15, 5) )
    for count, (model, zero_line ) in enumerate( zip( models_to_graph, zero_lines ) ):
        prime = ax[ count ].contourf( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.liq_alt, model.data['clw_alt_lat'][:,constants.lat_confine_1:constants.lat_confine_2]*1000, vmin=0, vmax=0.06, cmap='coolwarm' )
        temp = ax[ count ].contour( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.liq_alt, ( model.data[ 'ta_alt_lat' ][:,constants.lat_confine_1:constants.lat_confine_2] - 273.15 ), colors='white' )
        temp.collections[ zero_line ].set_linewidth( 3 )
        temp.collections[ zero_line ].set_color( 'white' )
        ax[ count ].clabel( temp, inline=1, fontsize=10 )
        ax[ count ].set_xlabel( 'Latitude' )
        ax[ count ].set_ylabel( 'Altitude (km)' )
        ax[ count ].set_title( model.name )
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)

    # the common color bar is based off the model labeled 'prime'
    cbar = fig.colorbar(prime, ax=ax.ravel().tolist(), orientation='horizontal', fraction = 0.05 )
    cbar.set_label('Mean Cloud Liquid Water Mass Fraction in Air (g/kg)')    
    cbar.set_clim(0, 0.06)

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])

    plt.savefig( location + '/Images/' + 'clw_contours' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight')
    plt.show()

#---------------------------------------------------------------------#


def model_clw_frac_contours( models_to_graph, zero_lines ):
    assert( len( models_to_graph ) == len( zero_lines ) )

    fig, ax = plt.subplots( nrows=1, ncols=len( models_to_graph ), figsize=(15, 5) )
    for count, (model, zero_line ) in enumerate( zip( models_to_graph, zero_lines ) ):
        prime = ax[ count ].contourf( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.liq_alt, model.data['clw_frac_alt_lat'][:,constants.lat_confine_1:constants.lat_confine_2], vmin=0, vmax=0.25, cmap='coolwarm' )
        temp = ax[ count ].contour( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.liq_alt, ( model.data[ 'ta_alt_lat' ][:,constants.lat_confine_1:constants.lat_confine_2] - 273.15 ), colors='white' )
        temp.collections[ zero_line ].set_linewidth( 3 )
        temp.collections[ zero_line ].set_color( 'white' )
        ax[ count ].clabel( temp, inline=1, fontsize=10 )
        ax[ count ].set_xlabel( 'Latitude' )
        ax[ count ].set_ylabel( 'Altitude (km)' )
        ax[ count ].set_title( model.name )
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)

    # the common color bar is based off the model labeled 'prime'
    cbar = fig.colorbar(prime, ax=ax.ravel().tolist(), orientation='horizontal', fraction = 0.05 )
    cbar.set_label('Mean Cloud Liquid Water Fraction')    
    cbar.set_clim(0, 0.25)

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])

    plt.savefig( location + '/Images/' + 'clw_frac_contours' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight')
    plt.show()
#---------------------------------------------------------------------#

def model_clwc_contours( models_to_graph, zero_lines ):
    assert( len( models_to_graph ) == len( zero_lines ) )

    fig, ax = plt.subplots( nrows=1, ncols=len( models_to_graph ), figsize=(15, 5) )
    for count, (model, zero_line ) in enumerate( zip( models_to_graph, zero_lines ) ):
        prime = ax[ count ].contourf( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.liq_alt, model.data['clwc_alt_lat'][:,constants.lat_confine_1:constants.lat_confine_2], vmin=0, vmax=0.06, cmap='coolwarm' )
        temp = ax[ count ].contour( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.liq_alt, ( model.data[ 'ta_alt_lat' ][:,constants.lat_confine_1:constants.lat_confine_2] - 273.15 ), colors='white' )
        temp.collections[ zero_line ].set_linewidth( 3 )
        temp.collections[ zero_line ].set_color( 'white' )
        ax[ count ].clabel( temp, inline=1, fontsize=10 )
        ax[ count ].set_xlabel( 'Latitude' )
        ax[ count ].set_ylabel( 'Altitude (km)' )
        ax[ count ].set_title( model.name )
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)

    # the common color bar is based off the model labeled 'prime'
    cbar = fig.colorbar(prime, ax=ax.ravel().tolist(), orientation='horizontal', fraction = 0.05 )
    cbar.set_label('Mean Cloud Liquid Water Mass Content (g/m^3)')    
    cbar.set_clim(0, 0.06)

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])

    plt.savefig( location + '/Images/' + 'clwc_contours' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight')
    plt.show()


#---------------------------------------------------------------------#


def model_cl_contours( models_to_graph, zero_lines ):
    assert( len( models_to_graph ) == len( zero_lines ) )

    fig, ax = plt.subplots( nrows=1, ncols=len( models_to_graph ), figsize=(15, 5) )
    for count, (model, zero_line ) in enumerate( zip( models_to_graph, zero_lines ) ):
        prime = ax[ count ].contourf( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.alt, model.data['cl_alt_lat'][:,constants.lat_confine_1:constants.lat_confine_2], vmin=0, vmax=0.5, cmap='coolwarm' )
        temp = ax[ count ].contour( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.alt, ( model.data[ 'full_ta_alt_lat' ][:,constants.lat_confine_1:constants.lat_confine_2] - 273.15 ), colors='white' )
        temp.collections[ zero_line ].set_linewidth( 3 )
        temp.collections[ zero_line ].set_color( 'white' )
        ax[ count ].clabel( temp, inline=1, fontsize=10 )
        ax[ count ].set_xlabel( 'Latitude' )
        ax[ count ].set_ylabel( 'Altitude (km)' )
        ax[ count ].set_title( model.name )
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)

    # the common color bar is based off the model labeled 'prime'
    cbar = fig.colorbar(prime, ax=ax.ravel().tolist(), orientation='horizontal', fraction = 0.05 )
    cbar.set_label('Cloud Fraction')    
    cbar.set_clim(0, 0.5)

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])

    plt.savefig( location + '/Images/' + 'cl_contours' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight')
    plt.show()


#---------------------------------------------------------------------#


def model_cli_contours( models_to_graph, zero_lines ):
    assert( len( models_to_graph ) == len( zero_lines ) )

    fig, ax = plt.subplots( nrows=1, ncols=len( models_to_graph ), figsize=(15, 5) )
    for count, (model, zero_line ) in enumerate( zip( models_to_graph, zero_lines ) ):
        prime = ax[ count ].contourf( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.alt, model.data['cli_frac_alt_lat'][:,constants.lat_confine_1:constants.lat_confine_2], vmin=0, vmax=0.25, cmap='coolwarm' )
        temp = ax[ count ].contour( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.alt, ( model.data[ 'full_ta_alt_lat' ][:,constants.lat_confine_1:constants.lat_confine_2] - 273.15 ), colors='white' )
        temp.collections[ zero_line ].set_linewidth( 3 )
        temp.collections[ zero_line ].set_color( 'white' )
        ax[ count ].clabel( temp, inline=1, fontsize=10 )
        ax[ count ].set_xlabel( 'Latitude' )
        ax[ count ].set_ylabel( 'Altitude (km)' )
        ax[ count ].set_title( model.name )
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=3)

    # the common color bar is based off the model labeled 'prime'
    cbar = fig.colorbar(prime, ax=ax.ravel().tolist(), orientation='horizontal', fraction = 0.05 )
    cbar.set_label('Mean Cloud Ice Water Fraction')    
    cbar.set_clim(0, 0.25)

    all_model_names = '_'.join( [ model.name for model in models_to_graph ])

    plt.savefig( location + '/Images/' + 'contours' + '_' +  all_model_names + ".svg", format="svg", bbox_inches='tight')
    plt.show()

############################################################################### regional plots


def region_plot( models_to_graph, quantity_to_graph ):
    ax = plt.axes( projection=ccrs.Mollweide(  central_longitude=180 ) )
    ax.coastlines()
    for model in  models_to_graph:
        # if model == 'CALIPSO':
        # change to np.roll(model.data[ quantity_to_graph + '_lat_lon' ], 179)

        if quantity_to_graph == 'albedo_reg':
            prime = ax.contourf( constants.lon, constants.lat, model.data[ 'albedo_reg' ], transform=ccrs.PlateCarree(), cmap='coolwarm' )
        elif quantity_to_graph == 'albedo_bias':
            prime = ax.contourf( constants.lon, constants.lat, (models['CERES'].data[ 'albedo_reg' ] - model.data[ 'albedo_reg' ]), transform=ccrs.PlateCarree(), cmap='coolwarm' )
        else:
            prime = ax.contourf( constants.lon, constants.lat, model.data[ quantity_to_graph + '_lat_lon' ], transform=ccrs.PlateCarree(), cmap='coolwarm' )
        ax.set_title( model.name )

        cbar = plt.colorbar(prime, orientation='horizontal')

        if quantity_to_graph == 'albedo_reg':
            cbar.set_label( 'Albedo' )
        elif quantity_to_graph == 'albedo_bias':
            cbar.set_label( 'Albedo Bias (CERES observation - model data)' )
        else:
            cbar.set_label( quantity[quantity_to_graph] )

        # all_model_names = '_'.join( [ model.name for model in models_to_graph ])

        plt.savefig( location + '/Images/' + 'regional' + '_' +  quantity_to_graph + '_' +  model.name + ".svg", format="svg", bbox_inches='tight')
        plt.show()



