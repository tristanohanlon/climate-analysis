"""
Created on Sun Oct  6 14:44:25 2019

@author: Tristan O'Hanlon - University of Auckland, Jonathan Rogers

    CMIP5-CESM1-CAM5
    CMIP5-GFDL-HIRAM-C360
    CMIP5-GISS-E2R
    CMIP5-IPSL-CM5A-LR
    CMIP5-MIROC5
    CMIP5-MRI-CGCM3
    
    CMIP6-CESM2-CAM6
    CMIP6-GFDL-AM4
    CMIP6-GISS-E21G
    CMIP6-IPSL-CM6A-LR
    CMIP6-MIROC6
    CMIP6-MRI-ESM2
    
    CCCM
    ECMWF
    CERES
    CALIPSO
"""

import constants
import clim_plot

location = constants.home


# plot zonal quantities: clt, clwvi or clivi with latitude
#clim_plot.latitude_plot_models( [ clim_plot.models[ 'CMIP5-GFDL-HIRAM-C360'], 
                                    # clim_plot.models[ 'CMIP6-GFDL-AM4'], 
                                    # clim_plot.models[ 'CMIP5-CESM1-CAM5'], 
                                    # clim_plot.models[ 'CMIP6-CESM2-CAM6'], 
                                    # clim_plot.models[ 'CERES'], 
                                    # clim_plot.models[ 'CALIPSO'], 
                                    # clim_plot.models[ 'ECMWF'] ], 'clwvi_lc', cmip5_range = True, cmip6_range = True )

# plot zonal quantities cl, clw or cli with altitude - 2 graphs: top = global, bottom = southern ocean
#clim_plot.altitude_plot_models( [ clim_plot.models[ 'CMIP5-GFDL-HIRAM-C360'], 
#                                     clim_plot.models[ 'CMIP6-GFDL-AM4'], 
#                                     clim_plot.models[ 'CMIP5-CESM1-CAM5'], 
#                                     clim_plot.models[ 'CMIP6-CESM2-CAM6'] ], 'clw', cmip5_range = True, cmip6_range = True )

# plot cl_t with temperature - 2 graphs: top = global, bottom = southern ocean
# clim_plot.cloud_temp_plot_models( [ clim_plot.models[ 'CMIP5-GFDL-HIRAM-C360'], 
#                                     clim_plot.models[ 'CMIP6-GFDL-AM4'], 
#                                     clim_plot.models[ 'CMIP5-CESM1-CAM5'], 
#                                     clim_plot.models[ 'CMIP6-CESM2-CAM6'], 
#                                     clim_plot.models[ 'CCCM'] ], cmip5_range = True, cmip6_range = True )


# plot clw_t with temperature - 2 graphs: top = global, bottom = southern ocean
clim_plot.liqc_temp_plot_models( [ clim_plot.models[ 'CMIP5-GFDL-HIRAM-C360'], 
                                    clim_plot.models[ 'CMIP6-GFDL-AM4'], 
                                    clim_plot.models[ 'CMIP5-CESM1-CAM5'], 
                                    clim_plot.models[ 'CMIP6-CESM2-CAM6'],
                                    clim_plot.models[ 'CCCM'] ], cmip5_range = True, cmip6_range = True )

# clw satellite contour plots. Zero lines specify the temperature contour of 0 degrees
# clim_plot.satellite_contours( [ clim_plot.models[ 'CALIPSO' ], clim_plot.models[ 'CCCM' ] ], [ 5, 5 ] )

# clw model pair contour plots. Zero lines specify the temperature contour of 0 degrees
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-CESM1-CAM5' ], clim_plot.models[ 'CMIP6-CESM2-CAM6' ] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-GFDL-HIRAM-C360'], clim_plot.models[ 'CMIP6-GFDL-AM4'] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-GISS-E2R' ], clim_plot.models[ 'CMIP6-GISS-E21G' ] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-IPSL-CM5A-LR'], clim_plot.models[ 'CMIP6-IPSL-CM6A-LR'] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-MIROC5' ], clim_plot.models[ 'CMIP6-MIROC6' ] ], [ 5, 5 ] )
# clim_plot.model_clw_frac_contours( [ clim_plot.models[ 'CMIP5-MRI-CGCM3'], clim_plot.models[ 'CMIP6-MRI-ESM2'] ], [ 5, 5 ] )


# cl model pair contour plots. Zero lines specify the temperature contour of 0 degrees
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-CESM1-CAM5' ], clim_plot.models[ 'CMIP6-CESM2-CAM6' ] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-GFDL-HIRAM-C360'], clim_plot.models[ 'CMIP6-GFDL-AM4'] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-GISS-E2R' ], clim_plot.models[ 'CMIP6-GISS-E21G' ] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-IPSL-CM5A-LR'], clim_plot.models[ 'CMIP6-IPSL-CM6A-LR'] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-MIROC5' ], clim_plot.models[ 'CMIP6-MIROC6' ] ], [ 6, 6 ] )
# clim_plot.model_cl_contours( [ clim_plot.models[ 'CMIP5-MRI-CGCM3'], clim_plot.models[ 'CMIP6-MRI-ESM2'] ], [ 6, 6 ] )

# cli model pair contour plots. Zero lines specify the temperature contour of 0 degrees
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-CESM1-CAM5' ], clim_plot.models[ 'CMIP6-CESM2-CAM6' ] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-GFDL-HIRAM-C360'], clim_plot.models[ 'CMIP6-GFDL-AM4'] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-GISS-E2R' ], clim_plot.models[ 'CMIP6-GISS-E21G' ] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-IPSL-CM5A-LR'], clim_plot.models[ 'CMIP6-IPSL-CM6A-LR'] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-MIROC5' ], clim_plot.models[ 'CMIP6-MIROC6' ] ], [ 6, 6 ] )
# clim_plot.model_cli_contours( [ clim_plot.models[ 'CMIP5-MRI-CGCM3'], clim_plot.models[ 'CMIP6-MRI-ESM2'] ], [ 6, 6 ] )

# plot regional clt, clwvi, clivi, clt_lc, clwvi_lc, aerosol_norm, mmrdust, mmroa, mmrso4, rtmt, rsut, rsdt, rsutcs quantities overlaying the world map
# clim_plot.region_plot( [ clim_plot.models['CALIPSO'] ], 'clt_lc' )