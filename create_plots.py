"""
Created on Sun Oct  6 14:44:25 2019

@author: Tristan O'Hanlon - University of Auckland, Jonathan Rogers

    CMIP5-AMIP-CESM1-CAM5 : 1979.01-2005.12 --- cl must be dived by 100
    CMIP5-AMIP-GFDL-CM3 : 1999.01-2008.12
    CMIP5-AMIP-GISS-E2R : 1951.01-2010.12 
    CMIP5-AMIP-IPSL-CM5A-LR : 1979.01-2009.12
    CMIP5-AMIP-MIROC5 : 1999.01-2008.12
    CMIP5-AMIP-MRI-CGCM3 : 1999.01-2010.02
    
    CMIP6-AMIP-CESM2-CAM6 : 1950.01-2014.12   
    CMIP6-AMIP-GFDL-CM4 : 1980.01-2014.12
    CMIP6-AMIP-GISS-E21G : 2001.01-2014.12
    CMIP6-AMIP-IPSL-CM6A-LR : 1979.01-2014.12
    CMIP6-AMIP-MIROC6 : 1999.01-2014.12
    CMIP6-AMIP-MRI-ESM2 : 1999.01-2014.12


    
    CCCM
    ECMWF
    CERES
    CALIPSO
"""

import constants
import clim_plot

location = constants.home

# plot zonal quantities: clt, clt_l, clwvi or clivi with latitude
# clim_plot.latitude_plot_models( [   clim_plot.cosp_models[ 'CMIP5-AMIP-GFDL-CM3'],
#                                     clim_plot.cosp_models[ 'CMIP6-AMIP-GFDL-CM4'],
#                                     clim_plot.cosp_models[ 'CMIP5-AMIP-CESM1-CAM5'],
#                                     clim_plot.cosp_models[ 'CMIP6-AMIP-CESM2-CAM6'],
#                                     clim_plot.models[ 'CALIPSO' ]], 'clt_l', cmip5_range = True, cmip6_range = True )

# plot layer quantities cl, clw or cli with altitude - 2 graphs: left = global, right = southern ocean
# clim_plot.altitude_plot_models( [ clim_plot.cosp_models[ 'CMIP5-AMIP-GFDL-CM3'], 
#                                     clim_plot.cosp_models[ 'CMIP6-AMIP-GFDL-CM4'],
#                                     clim_plot.cosp_models[ 'CMIP5-AMIP-CESM1-CAM5'], 
#                                     clim_plot.cosp_models[ 'CMIP6-AMIP-CESM2-CAM6'], 
#                                     clim_plot.models[ 'CALIPSO'] 
#                                     ], 'cl', cmip5_range = True, cmip6_range = True )

# plot cl_t with temperature - 2 graphs: top = global, bottom = southern ocean
# clim_plot.cloud_temp_plot_models( [ clim_plot.models[ 'CMIP5-AMIP-GFDL-CM3'], 
#                                     clim_plot.models[ 'CMIP6-AMIP-GFDL-CM4'], 
#                                     clim_plot.models[ 'CMIP5-AMIP-CESM1-CAM5'], 
#                                     clim_plot.models[ 'CMIP6-AMIP-CESM2-CAM6'] ], cmip5_range = True, cmip6_range = True )


# plot clw_frac_t with temperature - 2 graphs: top = global, bottom = southern ocean
# clim_plot.clw_frac_temp_plot_models( [ clim_plot.models[ 'CMIP5-AMIP-GFDL-CM3'],
#                                     clim_plot.models[ 'CMIP6-AMIP-GFDL-CM4'], 
#                                     clim_plot.models[ 'CMIP5-AMIP-CESM1-CAM5'],
#                                     clim_plot.models[ 'CMIP6-AMIP-CESM2-CAM6'] ], 'clw_frac', cmip5_range = True, cmip6_range = True )


# plot cli_frac_t with temperature - 2 graphs: top = global, bottom = southern ocean
# clim_plot.cli_frac_temp_plot_models( [ clim_plot.models[ 'CMIP5-AMIP-GFDL-CM3'],
#                                     clim_plot.models[ 'CMIP6-AMIP-GFDL-CM4'], 
#                                     clim_plot.models[ 'CMIP5-AMIP-CESM1-CAM5'],
#                                     clim_plot.models[ 'CMIP6-AMIP-CESM2-CAM6'] ], 'cli_frac', cmip5_range = True, cmip6_range = True )


# single cl contour plot. Zero lines specify the temperature contour of 0 degrees
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP5-AMIP-CESM1-CAM5' ] ], [ 5 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP5-AMIP-GFDL-CM3' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP5-AMIP-GISS-E2R' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP5-AMIP-IPSL-CM5A-LR' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP5-AMIP-MIROC5' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP5-AMIP-MRI-CGCM3' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-CESM2-CAM6' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-GFDL-CM4' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-IPSL-CM6A-LR' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-MIROC6' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-MRI-ESM2' ] ], [ 6 ] )
# clim_plot.cl_contour( [ clim_plot.models[ 'CCCM' ] ], [ 6 ] )

# single clw_frac_contour plot. Zero lines specify the temperature contour of 0 degrees
# clim_plot.clw_frac_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-CESM2-CAM6' ] ], [ 5 ] )
# clim_plot.clw_frac_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-GFDL-CM4' ] ], [ 5 ] )
# clim_plot.clw_frac_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-IPSL-CM6A-LR' ] ], [ 5 ] )
# clim_plot.clw_frac_contour( [ clim_plot.models[ 'CALIPSO' ] ], [ 5 ] )

# single cli_frac_contour plot. Zero lines specify the temperature contour of 0 degrees
# clim_plot.cli_frac_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-CESM2-CAM6' ] ], [ 5 ] )
# clim_plot.cli_frac_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-GFDL-CM4' ] ], [ 5 ] )
# clim_plot.cli_frac_contour( [ clim_plot.cosp_models[ 'CMIP6-AMIP-IPSL-CM6A-LR' ] ], [ 5 ] )
# clim_plot.cli_frac_contour( [ clim_plot.models[ 'CALIPSO' ] ], [ 5 ] )

# clw model pair contour plots. Zero lines specify the temperature contour of 0 degrees
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-AMIP-CESM1-CAM5' ], clim_plot.models[ 'CMIP6-AMIP-CESM2-CAM6' ] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-AMIP-GFDL-CM3'], clim_plot.models[ 'CMIP6-AMIP-GFDL-CM4'] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-AMIP-GISS-E2R' ], clim_plot.models[ 'CMIP6-AMIP-GISS-E21G' ] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-AMIP-IPSL-CM5A-LR'], clim_plot.models[ 'CMIP6-AMIP-IPSL-CM6A-LR'] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-AMIP-MIROC5' ], clim_plot.models[ 'CMIP6-AMIP-MIROC6' ] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-AMIP-MRI-CGCM3'], clim_plot.models[ 'CMIP6-AMIP-MRI-ESM2'] ], [ 5, 5 ] )


# cl model pair contour plots. Zero lines specify the temperature contour of 0 degrees
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-AMIP-CESM1-CAM5' ], clim_plot.models[ 'CMIP6-AMIP-CESM2-CAM6' ] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-AMIP-GFDL-CM3'], clim_plot.models[ 'CMIP6-AMIP-GFDL-CM4'] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-AMIP-GISS-E2R' ], clim_plot.models[ 'CMIP6-AMIP-GISS-E21G' ] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-AMIP-IPSL-CM5A-LR'], clim_plot.models[ 'CMIP6-AMIP-IPSL-CM6A-LR'] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-AMIP-MIROC5' ], clim_plot.models[ 'CMIP6-AMIP-MIROC6' ] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-AMIP-MRI-CGCM3'], clim_plot.models[ 'CMIP6-AMIP-MRI-ESM2'] ], [ 6, 6 ] )

# cli model pair contour plots. Zero lines specify the temperature contour of 0 degrees
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-AMIP-CESM1-CAM5' ], clim_plot.models[ 'CMIP6-AMIP-CESM2-CAM6' ] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-AMIP-GFDL-CM3'], clim_plot.models[ 'CMIP6-AMIP-GFDL-CM4'] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-AMIP-GISS-E2R' ], clim_plot.models[ 'CMIP6-AMIP-GISS-E21G' ] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-AMIP-IPSL-CM5A-LR'], clim_plot.models[ 'CMIP6-AMIP-IPSL-CM6A-LR'] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-AMIP-MIROC5' ], clim_plot.models[ 'CMIP6-AMIP-MIROC6' ] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-AMIP-MRI-CGCM3'], clim_plot.models[ 'CMIP6-AMIP-MRI-ESM2'] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CALIPSO'], clim_plot.models[ 'CALIPSO'] ], [ 6, 6 ] )

# plot regional clt, clw_frac, cli_frac, clt_l, clwvi_l, albedo_reg, loaddust, loadss, rtmt, rsut, rsdt, rsutcs quantities overlaying the world map
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP5-AMIP-CESM1-CAM5']], 'clt' )
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP5-AMIP-GFDL-CM3']], 'clt' )
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP5-AMIP-GISS-E2R']], 'clt' ) 
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP5-AMIP-IPSL-CM5A-LR']], 'clt' ) 
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP5-AMIP-MIROC5']], 'clt' )
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP5-AMIP-MRI-CGCM3']], 'clt' )
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP6-AMIP-CESM2-CAM6']], 'clt' )                                  
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP6-AMIP-GFDL-CM4']], 'clt' ) 
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP6-AMIP-GISS-E21G']], 'clt' ) 
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP6-AMIP-IPSL-CM6A-LR']], 'clt' ) 
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP6-AMIP-MIROC6']], 'clt' )
# clim_plot.region_plot( [ clim_plot.models[ 'CMIP6-AMIP-MRI-ESM2']], 'clt' ) 
# clim_plot.region_plot( [ clim_plot.models[ 'CALIPSO']], 'clt' ) 
                            
# plot regional clt, clw_frac, cli_frac bias quantities overlaying the world map
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP5-AMIP-CESM1-CAM5']], 'clt_l' )
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP5-AMIP-GFDL-CM3']], 'clt_l' )
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP5-AMIP-GISS-E2R']], 'clt_l' ) 
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP5-AMIP-IPSL-CM5A-LR']], 'clt_l' ) 
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP5-AMIP-MIROC5']], 'clt_l' )
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP5-AMIP-MRI-CGCM3']], 'clt_l' )
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP6-AMIP-CESM2-CAM6']], 'clw_frac' )                                  
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP6-AMIP-GFDL-CM4']], 'clw_frac' ) 
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP6-AMIP-IPSL-CM6A-LR']], 'clw_frac' ) 
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP6-AMIP-MIROC6']], 'clt_l' )
# clim_plot.calipso_region_bias_plot( [ clim_plot.cosp_models[ 'CMIP6-AMIP-MRI-ESM2']], 'clt_l' ) 
