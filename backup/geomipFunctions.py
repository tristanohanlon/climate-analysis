# Module containing useful functions and classes for analyzing GeoMIP (and CMIP) model
# output. Some of these are adapted from "cesmFunctions.py" which was created to
# analyze CESM runs done for the Ackerman et al. "strategy" study. 
#
# Rick Russotto
# Started May 2, 2016

import numpy as np
import scipy.interpolate


#"Enhanced time mean" function to deal with different month lengths
# array: numpy array to take time mean of
# calendar: which calendar to use (Gregorian, 365-day or 360-day)
# leapYears: array of leap years if Gregorian calendar, assuming the start year is year 0.
# mask: whether to use the masked array average function (slows it down but works with masked arrays) 
#
# default of "calendar" is 360 for backward compatibility (no weighting according to month length)                      
#                      
# Assumes time is first dimension.                  
def enhancedTimeMean(array, calendar='360', leapYears = 'none', mask=False):
    if calendar == '360': 
        return np.mean(array, axis=0) #All 30-day months so no weighting necessary
    elif calendar == '365':
        weights = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        weights = np.tile(weights,np.shape(array)[0]/12) #Repeat month lengths by number of years
        if mask:
            return np.ma.average(array,weights=weights,axis=0) #Note this creates a masked array which may slow things down. 
                                         #But if I use np.average instead of np.ma.average, it will not take mask 
                                         #into account when I have masked arrays. 
        else:
            return np.average(array,weights=weights,axis=0)
    elif calendar == 'Gregorian': 
        weights = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        weights = np.tile(weights,np.shape(array)[0]/12) #Repeat month lengths by number of years
        weights[leapYears*12+1] = 29 #Change February days to 29 in the leap years
        if mask:
            return np.ma.average(array,weights=weights,axis=0)
        else:
            return np.average(array,weights=weights,axis=0)            
    else:
        print('Invalid calendar; nothing returned')
        return




# Function to calculate the vertical energy flux convergence into 
# the atmospheric column. 
#
# Inputs: 
# DataDict: dictionary of NetCDF4 Dataset objects containing:
# DataDict['rsdt']: Shortwave down at TOA
# DataDict['rsut']: Shortwave up at TOA
# DataDict['rlut']: Longwave up at TOA
# DataDict['rsds']: Shortwave down at surface
# DataDict['rlds']: Longwave down at surface
# DataDict['rsus']: Shortwave up at surface
# DataDict['rlus']: Longwave up at surface
# DataDict['hfls']: Latent heat up at surface
# DataDict['hfss']: Sensible heat up at surface
# firstYear, lastYear: boundaries for time mean (first year in the data is year 0)
# flipSensible: change sign of sensible heat flux (GISS model did this)
# calendar, leapYears: see "enhancedTimeMean"
# use_prsn: whether to consider term from latent heat of fusion of snow falling into ocean
# oceanOnly: if false, consider snow falling on land as well.
def AtmosEnergyVerticalFluxConvergence(DataDict, firstYear, lastYear, flipSensible=False, 
                                          calendar='360', leapYears='none', use_prsn=True, 
                                          landFracPercent = True, oceanOnly = True):
    #Calculate the energy flux into the column
    if not(flipSensible):
        fluxConv = (  DataDict['rsdt'].variables['rsdt'][firstYear*12:lastYear*12+12,:,:] 
                    - DataDict['rsut'].variables['rsut'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rlut'].variables['rlut'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rsds'].variables['rsds'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rlds'].variables['rlds'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['rsus'].variables['rsus'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['rlus'].variables['rlus'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['hfls'].variables['hfls'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['hfss'].variables['hfss'][firstYear*12:lastYear*12+12,:,:])
    else:
        fluxConv = (  DataDict['rsdt'].variables['rsdt'][firstYear*12:lastYear*12+12,:,:] 
                    - DataDict['rsut'].variables['rsut'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rlut'].variables['rlut'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rsds'].variables['rsds'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rlds'].variables['rlds'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['rsus'].variables['rsus'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['rlus'].variables['rlus'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['hfls'].variables['hfls'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['hfss'].variables['hfss'][firstYear*12:lastYear*12+12,:,:])
    #Add latent heat of fusion term in if specified
    if use_prsn:
        Lf = 3.34e5
        if oceanOnly == True: #Only consider snowfall over ocean
            if landFracPercent == True:                             #Remember land fraction was in %, except for MIROC
                fusionTerm = (DataDict['prsn'].variables['prsn'][firstYear*12:lastYear*12+12,:,:] * 
                              (1.0-DataDict['sftlf'].variables['sftlf'][:]/100))*Lf #1 minus land fraction to get ocean
            else:
                fusionTerm = (DataDict['prsn'].variables['prsn'][firstYear*12:lastYear*12+12,:,:] * 
                              (1.0-DataDict['sftlf'].variables['sftlf'][:]))*Lf #1 minus land fraction to get ocean
        else: #Consider snowfall over ocean and land
            fusionTerm = DataDict['prsn'].variables['prsn'][firstYear*12:lastYear*12+12,:,:]*Lf
        fluxConv = fluxConv + fusionTerm            
    #fluxConvMultiYearMean = np.mean(fluxConv,axis=0)  #Old algorithm--doesn't take different month length into account
    fluxConvMultiYearMean = enhancedTimeMean(fluxConv, calendar=calendar, leapYears=leapYears)            
    return fluxConvMultiYearMean
    
    
# Equivalent function for the ocean, using the surface flux output. (Maybe needs snowfall term too?)
#
# Uses same DataDict as above, except doesn't use the 3 TOA variables. It will still work if you
# pass it a dictionary that does include those. 
def OceanEnergyVerticalFluxConvergence(DataDict, firstYear, lastYear, flipSensible=False):
    #Calculate the energy flux into the ocean (signs flipped from atmosphere)
    if not(flipSensible):
        fluxConv = (  DataDict['rsds'].variables['rsds'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['rlds'].variables['rlds'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rsus'].variables['rsus'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rlus'].variables['rlus'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['hfls'].variables['hfls'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['hfss'].variables['hfss'][firstYear*12:lastYear*12+12,:,:])
    else:
        fluxConv = (  DataDict['rsds'].variables['rsds'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['rlds'].variables['rlds'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rsus'].variables['rsus'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['rlus'].variables['rlus'][firstYear*12:lastYear*12+12,:,:]
                    - DataDict['hfls'].variables['hfls'][firstYear*12:lastYear*12+12,:,:]
                    + DataDict['hfss'].variables['hfss'][firstYear*12:lastYear*12+12,:,:])
    fluxConvMultiYearMean = np.mean(fluxConv,axis=0)
    return fluxConvMultiYearMean        
    
# Function to calculate the convergence into the atmospheric column of just latent heat flux, based on 
# the precipitation and evaporation data. Really, this is flux of moisture, times Lv.
#
# Inputs: 
# data_pr: NetCDF Dataset object for precipitation (kg/m^2/s)
# data_hfls: NetCDF Dataset object for surface latent heat flux (upward)  (W m^-2)
# Lv: latent heat of vaporization (J/kg)
# Lf: latent heat of fusion (J/kg)
# firstYear, lastYear: boundaries for time mean (first year in the data is year 0)
# data_prsn: NetCDF Dataset object for snowfall. In some models, the latent heat of fusion from snowfall
#            was added to the "hfls" variable, so we need to subtract it to get back the evaporation. 
#            If not specified, not taken into account. (Maybe would be better to use evspsbl instead?)
def AtmosLatentHeatFluxConvergence(data_pr, data_hfls, firstYear, lastYear, Lv=2.5e6, Lf=3.34e5,
                                      calendar='360', leapYears='none', data_prsn='none'):
    hfls = data_hfls.variables['hfls'][firstYear*12:lastYear*12+12,:,:]
    pr = data_pr.variables['pr'][firstYear*12:lastYear*12+12,:,:]
    prLv = pr*Lv #This results in W m^-2 if units consistent with variable descriptions above
    if data_prsn == 'none':
        fluxConv = hfls - prLv
    else:
        sn = data_prsn.variables['prsn'][firstYear*12:lastYear*12+12,:,:]
        snLf = sn*Lf        
        fluxConv = hfls - snLf - prLv
    #fluxConvMultiYearMean = np.mean(fluxConv,axis=0)
    fluxConvMultiYearMean = enhancedTimeMean(fluxConv, calendar=calendar, leapYears=leapYears)   
    return fluxConvMultiYearMean        

# Function to calculate zonally integrated northward energy transport by the atmosphere based on 
# the vertical flux of energy into the atmospheric column. Specifically this is transport of Moist Static Energy.
# This should also work for the ocean despite its name (although correction factor form may be different, which may
# necessitate a separate ocean function later). 
# This will also work for the latent heat flux calculated using "AtmosLatentHeatFluxConvergence". 
#
# Inputs:
# fluxConv: net vertical energy flux into atmospheric column,
#           calculated using AtmosEnergyVerticalFluxConvergence() function
#           (This is 2d in lat, lon dimensions, having already been time averaged)
# lat: model latitudes (degrees)
# a: radius of earth (m) (default: 6.37e6)
# correct (default True): whether to subtract out a bias associated with northward transport at North Pole
# Must specify one of the following for grid cell area calculations: 
#
# areacella: CMIP5 reported grid cell areas (as a 2D array)
# lon: model longitudes (degrees)
#
# Result: a 1D vector of northward atmospheric energy transport 
def AtmosEnergyTransportNorthward(fluxConv, lat, a= 6.37e6, correct=True, areacella='none', lon='none'):
    if areacella=='none': #Didn't specify separate grid cell area file
        #First calculate grid cell area
        latDiff = lat[1]-lat[0]
        lonDiff = lon[1]-lon[0]
        gridCellAreas = a*a*(latDiff*np.pi/180.)*(lonDiff*np.pi/180.)*np.cos(lat*np.pi/180.) #Vector based on latitude
        #Now calculate total energy flux in grid boxes in Watts (not W/m^2)
        fluxWatts = fluxConv*gridCellAreas[:,None] #Dimensions of fluxConv are lat,lon;
                                               #"None" ensures weighting by lat, not lon, even if square grid.
    else: #Use specified areas
        gridCellAreas = areacella
        fluxWatts = fluxConv*gridCellAreas
        #print('Used areacella') 
    #Integrate over longitudes
    fluxWattsZonalSum = np.sum(fluxWatts,axis=1) #Result: vector, dimension: latitude
    #Now cumuliatively integrate over latitudes
    energyTransportNorthward = np.cumsum(fluxWattsZonalSum)
    # Often the calculated northward transport at the North Pole isn't exactly zero, because model is not conserving 
    # energy. Need to subtract out this bias so that both poles have zero meridional transport. 
    if correct: 
        N = energyTransportNorthward[len(energyTransportNorthward)-1] #Northward energy transport at NP
        E = N/2.*(1+np.sin(lat*np.pi/180.)) #See paper notes, 5 July 2016 for derivation. This form assumes each latitude band's
        energyTransportNorthward = energyTransportNorthward - E #...contribution to the error is proportional to its area.
    return energyTransportNorthward
    
    
# Standalone correction function for northward transport at North Pole, for a given energy transport array
# (also need model latitudes, in degrees, as input argument)
def correctNorthPoleTransport(transportArray, lats):
    N = transportArray[len(transportArray)-1] 
    E = N/2.*(1+np.sin(lats*np.pi/180.)) 
    correctedArray = transportArray - E 
    return correctedArray
    
# Version of correction function for when transport at southernmost latitude is not 0 
# (this happens in moist EBM results--partly because we got rid of any grid boxes at the pole)
# Subtract the residual from 0 at southernmost latitude from the entire dataset first, then do 
# the ramp correction.     
def correctNorthPoleTransportSouthOffset(transportArray, lats):
    S = transportArray[0]
    intermediateArray = transportArray - S
    N = intermediateArray[len(intermediateArray)-1]
    E = N/2.*(1+np.sin(lats*np.pi/180.)) 
    correctedArray = intermediateArray - E 
    return correctedArray
    
# Function to calculate northward energy transport across the equator
# (Special case of "AtmosEnergyTransportNorthward")
#
# Inputs:
# fluxConv: net vertical energy flux into atmospheric column,
#           calculated using AtmosEnergyVerticalFluxConvergence() function
#           (This is 2d in lat, lon dimensions, having already been time averaged)
# lat: model latitudes (degrees)
# lon: model longitudes (degrees)
# a: radius of earth (m)
#
# Result: a single number for total energy transport across the equator (positive northward)
def AtmosEnergyTransportCrossEquator(fluxConv, lat, lon, a):
    #First find the full energy transport vector for the whole earth
    #(Need cumulative sum because only at the pole can energy only go one direction)
    energyTransportNorthward = AtmosEnergyTransportNorthward(fluxConv, lat, lon, a)
    #Now find the cross-equatorial transport. 
    #Two cases based on odd or even number of latitudes
    numLats = len(lat)
    if not(np.mod(numLats,2)): #Even number of latitude bands--equator is between 2 boxes
        #Since transport is out of the grid box at north end, use the last element in 1st half of array
        energyTransportCrossEquator = energyTransportNorthward[numLats/2-1]
    else: #Odd number of latitude bands--equator is within one box
        #Rely on integer floor division here--don't want to mess up in Python 3 so use //
        minIndex = numLats//2-1 #e.g. if 5 latitude bands, indexed from 0 to 4, this gives 1
        maxIndex = numLats//2   #...and this gives 2; we want the average of [1] and [2] index
        energyTransportCrossEquator = (energyTransportNorthward[minIndex] + 
                                       energyTransportNorthward[maxIndex])/2. 
    return energyTransportCrossEquator
    
    


# Function to calculate the mean of a netCDF4 Dataset variable between minimum
# and maximum latitudes--returns a 1D array with the monthly means. 
#
# Inputs:
# Data: NetCDF4 dataset for the variable of interest
# varname: variable name
# minLat, maxLat: minimum and maximum latitudes to consider (degrees)
# latString: what latitude is called in the netCDF file
def bandMean(Data, varname, minLat, maxLat, latString = 'lat'):
    #Subsetting-indexing stuff
    #globalLats = Data.variables['lat'][:]
    globalLats = Data.variables[latString][:]
    a = globalLats >= minLat
    b = globalLats <= maxLat
    e=a*b
    bandLats = globalLats[e]
    foundIndexMinLat = False
    foundIndexMaxLat = False
    i = 0  #Loop counter
    while foundIndexMinLat == False:
        if e[i] == True:
            indexMinLat = i
            foundIndexMinLat = True
        i += 1
    i = 1  #With negative indices this starts at the end
    while foundIndexMaxLat == False:
        if e[-i] == True:
            indexMaxLat = len(e)-i
            foundIndexMaxLat = True
        i += 1
    #Now retrieving the data in the band--at all times
    globalData = Data.variables[varname][:]
    bandData = globalData[:, indexMinLat:indexMaxLat+1, :] 
    #Now, calculating weights and replicating them to the 3D array
    areaWeights = np.cos(bandLats*np.pi/180)
    areaWeights3D = np.swapaxes(np.tile(areaWeights,
                                (np.shape(bandData)[0],np.shape(bandData)[2],1)),
                                1,2)
    #Now calculating the weighted mean
    weighted3DMatrix = bandData*areaWeights3D
    sumWeighted = np.sum(weighted3DMatrix, axis=(1,2))
    sumWeights3D = np.sum(areaWeights3D, axis=(1,2))
    weightedMean = sumWeighted/sumWeights3D
    return weightedMean
    

# The following three functions calculate meridional transport of latent, 
# geopotential and thermal energy from daily output for meridional wind, 
# geopotential height, specific humidity, and air tempertaure, and take
# a mean over specified range of days. 
#
# Inputs: 
# data_va:  NetCDF Dataset object for the daily meridional wind output (m/s)
# data_ta:  NetCDF Dataset object for the daily temperature output (K)
# data_hus: NetCDF Dataset object for the daily specific humidity output (kg/kg)
# data_zg:  NetCDF Dataset object for the daily geopotential height output (m)
# firstDay: first day of the multi-day mean (measured from 0)
# lastDay:  last day of the multi-day mean (measured from 0)
#
# Outputs: 3D NumPy arrays (plev, lat, lon) of energy fluxes in J/kg*m/s
# 
# Input dimensions: time, plev, lat, lon (in NetCDF Dataset object)

def meridionalTransportLatentEnergyMultiDayMean(data_va, data_hus, firstDay, lastDay):
    #Physical constants
    Lv = 2.5e6 #Latent heat of vaporization at 0 deg C (J/kg (H2O))
    #Extract data from the NetCDF files
    va = data_va.variables['va'][firstDay:lastDay+1,:,:,:]  #m/s
    q = data_hus.variables['hus'][firstDay:lastDay+1,:,:,:]  #kg (H2O)/kg (air)
    #Multiply specific humidity by meridional velocity to get the moisture flux
    vq = va*q #kg/kg*m/s
    #Take time mean
    vqbar = np.mean(vq,axis=0)
    #Multiply by Lv to get in energy units
    Lvqbar = Lv*vqbar #J/kg*m/s
    return Lvqbar

def meridionalTransportGeopotentialEnergyMultiDayMean(data_va, data_zg, firstDay, lastDay):
    #Physical constants
    g = 9.81 #gravitational acceleration (m s^-2)
    #Extract data from the NetCDF files
    va = data_va.variables['va'][firstDay:lastDay+1,:,:,:]  #m/s
    zg = data_zg.variables['zg'][firstDay:lastDay+1,:,:,:]  #m
    #Multiply together to get flux
    vz = va*zg #m^2/s
    #Take time mean
    vzbar = np.mean(vz,axis=0)
    #Multiply by g to get in energy units
    gvzbar = g*vzbar #m^3/s^3     (J = kg * m^2/s^2 so this is J/kg*m/s)
    return gvzbar 
    

def meridionalTransportThermalEnergyMultiDayMean(data_va, data_ta, firstDay, lastDay):
    #Physical constants
    Cp = 1004. #Heat capacity of air at constant pressure (J kg^-1 K^-1)
    #Extract data from the NetCDF files
    va = data_va.variables['va'][firstDay:lastDay+1,:,:,:]  #m/s
    ta = data_ta.variables['ta'][firstDay:lastDay+1,:,:,:]  #K
    #Multiply together to get flux
    vt = va*ta #K*m/s
    #Take time mean
    vtbar = np.mean(vt,axis=0)
    #Multiply by Cp to get in energy units
    Cpvtbar = Cp*vtbar #J/kg*m/s
    return Cpvtbar
    
#Generic function to do zonal and vertical integration for the quantities computed by the above functions. 
#
# Inputs: 
# flux: 3D array of energy fluxes (units: J/kg*m/s) computed by above 3 functions (dimensions: plev, lat, lon)
# plev: the pressure levels in Pa
# lat: the latitudes in degrees
# lon: the longitudes in degrees
def meridionalTransportZonalVerticalIntegration(flux,plev,lat,lon):
    g = 9.81 #gravitational acceleration (m s^-2)
    a = 6.37e6 #Radius of earth (meters)
    #First do vertical integration  (= integral of flux*dp/g (see Overland et al., 1996, eq. 1))
    pdiff = -1.*np.diff(plev)       #Differences between the pressure levels (Pa) (listed descending so flip sign)
    fluxInterp = (flux[0:len(plev)-1,:,:]+flux[1:len(plev),:,:])/2. #Average flux values midway between pressure levels
    verticallyIntegratedFlux = np.sum(fluxInterp*pdiff[:,None,None]/g, axis=0)  #Vertically integrated flux (W/m) (lat,lon)
    #Now do zonal integration (= integral of vertically integrated flux * dx)
    lonDiff = lon[1]-lon[0] 
    gridCellWidth = (lonDiff*np.pi/180.)*a*np.cos(lat*np.pi/180.) #Width of each grid cell (m) (1-D array, dimension: lat)
    zonallyIntegratedFlux = np.sum(verticallyIntegratedFlux*gridCellWidth[:,None],axis=1) #Northward energy transport (W) (lat)
    return zonallyIntegratedFlux
    
# Alternative version of above function that uses original flux values and weights according to distance between 
# midpoints between pressure levels (to avoid mask "bleeding" into level above bottom, if bottom masked due to terrain)    
def meridionalTransportZonalVerticalIntegrationV2(flux,plev,lat,lon):
    g = 9.81 #gravitational acceleration (m s^-2)
    a = 6.37e6 #Radius of earth (meters)
    #First do vertical integration  (= integral of flux*dp/g (see Overland et al., 1996, eq. 1))
    pMidpoints = (plev[0:len(plev)-1]+plev[1:len(plev)])/2.
    pWeights = -1.*np.diff(pMidpoints) #Differences between the pressure levels (Pa) (listed descending so flip sign)
    pWeights = np.append(plev[0]-pMidpoints[0], pWeights)
    pWeights = np.append(pWeights, pMidpoints[-1])
    #print(plev) #Debug
    #print(pWeights) #Debug
    verticallyIntegratedFlux = np.sum(flux*pWeights[:,None,None]/g, axis=0)  #Vertically integrated flux (W/m) (lat,lon)
    #Now do zonal integration (= integral of vertically integrated flux * dx)
    lonDiff = lon[1]-lon[0] 
    gridCellWidth = (lonDiff*np.pi/180.)*a*np.cos(lat*np.pi/180.) #Width of each grid cell (m) (1-D array, dimension: lat)
    zonallyIntegratedFlux = np.sum(verticallyIntegratedFlux*gridCellWidth[:,None],axis=1) #Northward energy transport (W) (lat)
    return zonallyIntegratedFlux
    
# Function to calculate the multi-model mean when you have arrays that are functions of latitudes 
# (e.g. zonal mean, northward energy transport) for multiple models. This function first interpolates
# each model's output to a common grid, then averages over the different models. 
#
# Inputs: 
# arrayDict: dictionary containing the 1-D arrays that vary as a function of latitude
# latDict: dictionary containing the 1-D arrays of latitudes
# interpLats: latitudes to interpolate to
# 
# The 2 dictionaries should have identical keys designating which model the data is from.
# Should be both the same size.
def multiModelMeanLatitudeVarying(arrayDict, latDict, interpLats):
    #interpLats = np.linspace(-90.,90.,181) #Interpolate to 1-degree grid #No-make this an argument
    interpMat = np.zeros((len(arrayDict.keys()),len(interpLats))) #Fill this matrix one model at a time
    i = 0
    for key in arrayDict.keys():
        interpMat[i,:] = np.interp(interpLats, latDict[key], arrayDict[key])
        i = i + 1
    multiModelMean = np.mean(interpMat, axis=0)
    return multiModelMean
    
# Function to regrid 2D arrays with latitude-pressure dimensions to a common 
# set of latitudes (assuming the pressure levels are all the same).
#
# This is analogous to the function above called "multiModelMeanLatitudeVarying" 
# but does loop of np.interp over the different pressure levels.
#   
# Inputs: 
# arrayDict: dictionary containing the 2-D arrays that vary as a function of latitude and pressure
#            (dimensions in that order)
# latDict: dictionary containing the 1-D arrays of latitudes for each model
# interpLats: latitudes to interpolate to
# numLevels: number of pressure levels (too hard to extract from the dict without knowing variable name)
#            (default 17: this is the number of standard CMIP5 pressure levels)
#
# Returns:
# A 3D array with dimensions (lat, height, model) 
# Order of models not guaranteed.
#    
# Also doesn't do mean over the models (do that in script). 
# The "arraySignAgreement3D" function should work for output from this, as well 
# as for the lat-lon version.
def collocateModelsLatitudePressure(arrayDict, latDict, interpLats, numLevels=17):
    i = 0 #Loop counter for models
    for key in arrayDict.keys():
        print('Interpolating to common latitudes for model: ' + modelNames[key])
        if i == 0: 
            array3D = np.zeros((len(interpLats), numLevels, len(arrayDict.keys()))) #initialize array to return
        for j in np.arange(numLevels): #Loop over pressure levels
            #print(j)
#            print(np.shape(interpLats))
#            print(np.shape(latDict[key]))
#            print(np.shape(arrayDict[key]))
#            print(interpLats)
#            print(latDict[key])
#            print(arrayDict[key])
            array3D[:,j,i] = np.interp(interpLats, latDict[key], arrayDict[key][j,:])
        i = i + 1
    return array3D
        
    
    
# Functions for taking multi-model statistics (means, agreement on sign, etc.) 
# on lat-lon grid. 
# General strategy is to loop through the models, 
# regrid the model output to a specific grid, and save to a 3D array (lat, lon, model). 
# Return this in one function; then have other functions that calculate things like 
# agreement.
    
# Function to create a 3D array of regridded data, given dicts of model data, lats and lons, (which must all have same keys), 
# and lats and lons to interpolate to (lats and lons are 1D arrays).
# The order of the models will be dependent on the input dicts, so I shouldn't write any code that assumes the order matters. 
# This uses the scipi.interpolate.griddata function, with the default linear interpolation.
#
# Required some reshaping to make data "unstructured". See "griddataTest.py" to understand the logic of this. 
# This seems to run reasonably fast, for interpolating to 2-degree grid. 
# Input data dimensions must be (lat, lon).
#
# Inside this function, x is lats, y is lons.
# For "method" (linear, nearest, etc.): see griddata documentation
def collocateModels(modelData, modelLats, modelLons, interpLats, interpLons, method='linear'):
    array3D = np.zeros((len(interpLats), len(interpLons), len(modelData.keys()))) #Create empty array to fill in
    interpPoints_x, interpPoints_y = np.meshgrid(interpLats, interpLons)
    interpPoints_x = np.transpose(interpPoints_x)
    interpPoints_y = np.transpose(interpPoints_y)
    i = 0 #loop counter
    for key, value in modelData.items():
        print('Interpolating for model: ')
        print(key)
        #Need to reshape the modelData to be unstructured, and also create lists of lats and lons. 
        #For linear and cubic methods, gives NaN for points outside "convex hull" of the data, meaning
        #NaNs at poles and the prime meridian. Poles are fine, but I don't want white strip at 
        #prime meridian. Get around this by duplicating end values (i.e. "pad" the data and longitudes). 
        if not(method == 'nearest'):
            #Need to choose a new variable name because otherwise Python changes the original modelData and modelLons called by the function
            #instead off making a copy, which causes a bug if the function is called twice in one script for the same model. 
            modelDataPadded = np.append(modelData[key][:,-1:], modelData[key], axis=1) #Append the end value to the beginning
            print('modelDataPadded: ' + str(np.shape(modelDataPadded)))
            modelDataPadded = np.append(modelDataPadded, modelDataPadded[:,1:2], axis=1) #Append the (original) beginning value to the end
            print('modelDataPadded: ' + str(np.shape(modelDataPadded)))
            lonSpacing = modelLons[key][1] - modelLons[key][0] #All models have evenly spaced longitude grid
            modelLonsPadded = np.insert(modelLons[key], 0, modelLons[key][0]-1*lonSpacing) #Add additional longitude to beginning
            print('modelLonsPadded: ' + str(np.shape(modelLonsPadded)))            
            modelLonsPadded = np.insert(modelLonsPadded, len(modelLonsPadded), modelLonsPadded[-1]+lonSpacing) #Add additional longitude to end
            print('modelLonsPadded: ' + str(np.shape(modelLonsPadded)))
            flatData = modelDataPadded.flat[:]
            print('flatData: ' + str(np.shape(flatData)))
            flatLons = np.tile(modelLonsPadded, (1, len(modelLats[key]))).flat[:]
            print('flatLons: ' + str(np.shape(flatLons)))
            tiledLats = np.tile(modelLats[key], (len(modelLonsPadded),1))
            print('tiledLats: ' + str(np.shape(tiledLats)))
        else:
            flatData = modelData[key].flat[:]
            flatLons = np.tile(modelLons[key], (1, len(modelLats[key]))).flat[:] 
            tiledLats = np.tile(modelLats[key], (len(modelLons[key]), 1))
        transposedTiledLats = np.transpose(tiledLats)
        print('transposedTiledLats: ' + str(np.shape(transposedTiledLats)))
        flatLats = transposedTiledLats.flat[:]
        print('flatLats: ' + str(np.shape(flatLats)))
        points = np.zeros((len(flatLats),2))
        points[:,0] = flatLats
        points[:,1] = flatLons
        print(np.shape(value))
        print(np.shape(points))
        print(np.shape(flatData))
        print(np.shape(interpPoints_x))
        print(np.shape(interpPoints_y))
        array3D[:,:,i] = scipy.interpolate.griddata(points, flatData, (interpPoints_x, interpPoints_y), method=method)
        i = i + 1
    return array3D        

        
# Function to return true when no more than a certain number ("maxDisagreeing") 
# of elements on a 3d array, looking along a specific dimension, disagree on the sign from the overall consensus. 
# Anticipate using this along "model" dimension (dim=2) for the output from the "collocateModels" function.
def arraySignAgreement3D(array, maxDisagreeing, dim = 2):
    isPositiveArray = array > 0 #(Won't count zero as agreeing on sign of change; this is more conservative)
    sumPositives = np.sum(isPositiveArray, axis=dim)
    n = np.shape(array)[dim]
    # OK, this algorithm is kind of confusing, so here's an explanation of how it works, using example of 9
    # models, when I want at least 7 to agree on sign. 
    # I have at each grid point the number models with a positive sign. 
    # If this number is 0, 1 or 2, then most of the models are negative with at most 2 positives. 
    # If this number is 7, 8, or 9, then most of the models are positive with at most 2 negatives. 
    # I create a Boolean array for each of these cases separately, then use a bitwise or to combine them
    # (get an error if I try to do this all in one line). 
    positivesWhenMostNegative = sumPositives <= maxDisagreeing
    negativesWhenMostPositive = sumPositives >= n - maxDisagreeing
    returnArray = positivesWhenMostNegative | negativesWhenMostPositive
    return returnArray

# Function to calculate mean of a variable from start of one year to end of a later year,
# for multi-year or multidecadal means. Not taking a global mean, despite the name--intended for mapping. 
# Years measured from 0. Optional monthOffset allows for starting from a month other than the first month
# in the dataset (useful for model runs that don't start in January, for example).
# 
def multiYearMean(Data, varname, firstYear, lastYear, monthOffset=0, calendar='360', leapYears='none'):
    allData = np.squeeze(Data.variables[varname][:]) #Squeeze in case it's one slice of a 4D array
    #yearMeanData = np.mean(allData[firstYear*12+monthOffset:lastYear*12+12+monthOffset,:,:],axis=0)
    yearMeanData = enhancedTimeMean(allData[firstYear*12+monthOffset:lastYear*12+12+monthOffset,:,:],
                                       calendar=calendar, leapYears=leapYears)
    return yearMeanData        
    
# Function to calculate the global mean of a variable
# from the GeoMIP output from start of one year to end of a later year. 
def globalMeanMultiYear(Data, varname, firstYear, lastYear, monthOffset=0): #"data" is netCDF4 data object
    yearMeanData = multiYearMean(Data, varname, firstYear, lastYear, monthOffset)  #Take time mean; gives 2D array (lat,lon)
    #print(np.shape(yearMeanData)) #debug
    latitudes = Data.variables['lat'][:]
    areaWeights = np.cos(latitudes*np.pi/180)
    areaWeights2D = np.tile(areaWeights,(np.shape(yearMeanData)[1],1))
    #print(np.shape(areaWeights2D)) #debug
    areaWeights2D = np.swapaxes(areaWeights2D,0,1)
    weighted2DMatrix = yearMeanData*areaWeights2D
    sumWeighted = np.sum(weighted2DMatrix)
    sumWeights = np.sum(areaWeights2D)
    weightedMean = sumWeighted/sumWeights
    return weightedMean       
    
# Function to calculate zonal mean of a vairable from start of one year to end of a later year. 
# Use the "globalMeanMultiYear" function above and then do zonal average. 
def zonalMeanMultiYear(Data, varname, firstYear, lastYear, monthOffset=0, calendar='360', leapYears='none'):
    yearMeanData = multiYearMean(Data, varname, firstYear, lastYear, monthOffset, calendar, leapYears)    
    zonalMeanData = np.mean(yearMeanData, axis=1)  
    return zonalMeanData
    
# Function to calculate the latitude of the median precipitation in a zonal mean 
# precipitation profile--that is, latitude at which half of precipitation is 
# to the north and half is to the south.
#
# Inputs:
# zonalMeanData: the zonal mean of the precipitation (calculated e.g. by zonalMeanMultiYear() function)
# lats: the latitudes corresponding to the zonalMeanData profile
# minLat: the minimum latitude to consider when taking the median
# maxLat: the maximum latitude to consider when taking the median
def precipMedian(zonalMeanData, lats, minLat, maxLat):
    a = lats >= minLat
    b = lats <= maxLat
    c = a*b
    zonalMeanSubset = zonalMeanData[c]
    latsSubset = lats[c]
    areaWeightedPrecip = zonalMeanSubset*np.cos(latsSubset*np.pi/180)
    totalPrecip = np.sum(areaWeightedPrecip)
    precipFraction = areaWeightedPrecip / totalPrecip
    cumulativePrecipFraction = np.cumsum(precipFraction)
    #Find highest latitude < 50% of rain and lowest latitude > 50% of rain
    indexLowerBound = np.argmin(1./(cumulativePrecipFraction-0.5)) #1/ means largest neg. number is closest to .5 on the - side
    indexUpperBound = np.argmax(1./(cumulativePrecipFraction-0.5)) #1/ means largest pos. number is closest to .5 on the + side
    lowerLat = latsSubset[indexLowerBound]
    upperLat = latsSubset[indexUpperBound]
    #Interpolate between the 2 latitudes based on how close 50% is to either side
    medianLat = lowerLat + (upperLat-lowerLat)*(.5-cumulativePrecipFraction[indexLowerBound])/(
                            cumulativePrecipFraction[indexUpperBound]-cumulativePrecipFraction[indexLowerBound])
    return medianLat
    
    
#Function to calculate zonal mean of multiple years of a height-varying variable already 
#in pressure coordinates (e.g. "ta" in all the models, or "cl" in the GISS model).
#Assuming dimensions are (time, plev, lat, lon). Since we're averaging over longitudes and time, 
#result will be 2D array: plev, lat.
#
# Inputs: 
# Data: NetCDF4 Dataset object containing the variable of interest. 
# varname: variable name for the variable in the Data object
#
# firstYear: first year of the multi-year mean (measured from start of run)
# lastYear:  last    "  "   "    "     "    "   
def zonalMeanMultiYearPressureCoordinates(Data, varname, firstYear, lastYear):
    datavar = Data.variables[varname][firstYear*12:lastYear*12+12,:,:,:]   #time, plev, lat, lon
    #print(datavar[0,:,:,0]) #Debug: see how missing data are handled
    #print(datavar[0,:,:,50]) #Debug: see how missing data are handled
    #What happens is, it reads them in as masked arrays. So mean works fine, but this poses problems for 
    #saving output as npy files. 
    dataZonalMean = np.nanmean(datavar, axis=3)  #Not quite sure how missing data are handled
    dataTimeMean = np.nanmean(dataZonalMean, axis=0)
    return dataTimeMean
    
    
# Function to interpolate height-varying model output variables to desired pressure levels and calculate
# zonal mean and multi-year mean. Takes some inspiration from code by Dan Vimont at
# http://www.aos.wisc.edu/~dvimont/matlab/
#
# This version is for Hybrid Sigma vertical coordinates.
#
# Inputs: 
# Data: NetCDF4 Dataset object containing the variable of interest. 
# varname: variable name for the variable in the Data object
# pLevels: pressure levels to interpolate to, in hPa. 
# firstYear: first year of the multi-year mean (measured from start of run)
# lastYear:  last    "  "   "    "     "    "
# modelTag: 2-letter tag for the model (see elsewhere in this module), 
#           for model-specific variations in coordinate variable names, etc.
# DataPS: optional--NetCDF4 Dataset object containing the surface pressure data. This is 
#         optional because some of the other files ("cl", etc.) already have this. 
# varPS: optional--PS specified as array instead of Dataset object (for CSIRO model)
#
# Output:
# A 2-D Numpy array with latitude on one axis and pressure on the other, and NaNs where no data exist 
# due to terrain (i.e. over Antarctica). Dimension order: pressure, latitude
def zonalMeanMultiYearPressureVaryingVariable(Data, varname, pLevels, firstYear, lastYear, modelTag, DataPS='none', varPS = 'none'):
            
    # Extract the data from the objects as well as the surface pressure
    datavar = Data.variables[varname][firstYear*12:lastYear*12+12,:,:,:]   #time, lev, lat, lon
    if varPS == 'none': 
        if DataPS == 'none':
            PS = Data.variables['ps'][firstYear*12:lastYear*12+12,:,:]  #time, lat, lon; Surface pressure (Pa)
        else:
            PS = DataPS.variables['ps'][firstYear*12:lastYear*12+12,:,:]  #time, lat, lon; Surface pressure (Pa)
    else: #Specify PS as Python array (variable) defined outside the function, rather than NetCDF Dataset object
        PS = varPS
    
    # Calculate pressure at the points where the Data are located.
    # Conversion from hybrid sigma levels to pressure levels
    # based on equation: #P[i,j,k] = a[k]*p0 + b[k]*PS[i,j] 
    # But P and PS vary in time whereas a and b don't, so need to use broadcasting.    
    B = Data.variables['b'][:] #dimension: lev
    if(modelTag == 'ca' or modelTag == 'ip' or modelTag == 'mp'): #CanESM or IPSL or MPI model
        AkP0 = Data.variables['ap'][:] #Pa; dimension: lev
    else: #BNU or CCSM4 or MIROC model (as of June 27th, 2016) #June 1, 2017: CESM and NorESM also this way
        #P0 = Data.variables['p0'][:] #Reference pressure (Pa)
        P0 = 100000 #Reference pressure (Pa). It's 100000 in all the models, but CCSM and CESM mistakenly used hPa.
        #print(np.shape(P0))
        A = Data.variables['a'][:]    #Dimension: lev
        #print(np.shape(A))
        AkP0 = A*P0
        #print(np.shape(AkP0))
    BkPSij = B[None,:,None,None]*PS[:,None,:,:] #Result: 4D array (time, lev, lat, lon)
    print('Shape of 4D array to interpolate:')
    print(np.shape(BkPSij))
    pressureMat = (AkP0[None,:,None,None] + BkPSij)/100. #Divide by 100 to go from Pa to hPa
    
    # Okay, now the hard part:
    # Interpolate the data from the model's native pressure levels to the desired pressure levels. 
    # The model's pressure levels vary with time, latitude, and longitude, whereas we need consistent
    # pressure levels for time and zonal averaging. So we have to do a linear interpolation in the 
    # vertical coordinate that varies with latitude, longitude and time. But doing this in nested loops
    # is way too slow, so we need to use matrix operations. 
    #
    # Further complicating the picture is the fact that sometimes the desired pressure level lies outside
    # the range of the model pressure levels, e.g. due to terrain. To account for this we need to put
    # nans in these places where there is no data, and use nanmean at the end for zonal and time mean.  
    #
    # General strategy: only one loop, over the new pressure levels. At each desired pressure level, calculate 
    # the difference between the native pressures and the desired pressure, and find the vertical 
    # indices of the native pressure closest to the desired pressure above and below. Find the native pressures and 
    # the data values corresponding to these indices in order to do the linear interpolation. In order to put nans 
    # where there is no data, keep track of indices columns where all native pressures are either 
    # above or below the desired pressure, and slot in nans in the appropriate place right before the final 
    # interpolation calculation.

    # First: reshape the native pressures and data into 2D arrays, with time/latitude/longitude all on one axis
    # and vertical coordinate on the other axis. This makes it simpler to extract data at 
    # the particular vertical index we want (which varies with time, latitude, and longitude) later on. 
    # But we will still use the 4D arrays in other parts of the calculation.
    s = np.shape(pressureMat)
    numTimes = s[0]
    numLevs = s[1]
    numLats = s[2]
    numLons = s[3]    
    #Make vertical level the last axis so that columns remain intact when matrices are reshaped    
    pressure2D = np.swapaxes(pressureMat,1,3) #Now it's time, lon, lat, level
    pressure2D = np.reshape(pressure2D, (numTimes*numLats*numLons, numLevs))
    data2D = np.swapaxes(datavar,1,3)
    data2D = np.reshape(data2D, (numTimes*numLats*numLons, numLevs))
        
    # Preallocate array to hold interpolated data
    interpMat = np.empty((pressureMat.shape[0], len(pLevels), pressureMat.shape[2], pressureMat.shape[3]))
    
    # Now: loop over the desired pressure levels and interpolate data to each one
    for k in range(0, len(pLevels)):
        print('Interpolating to ' + str(pLevels[k]) + ' hPa')
        
        # This code block: find the upper boundaries of the native pressure and data for interpolation
        print('Positive side')
        pressureMatDiffPos = pressureMat - pLevels[k] #Result: 4D array of differences between native and desired pressure.
                                                      #Positive values: higher native pressure than the level we're filling
                                                      #Negative values: lower  native pressure than the level we're filling
        # We're only interested in positive values (higher pressure end) for now, 
        # so set negative ones to a very high number and then find the index associated with the minimum along the 
        # vertical coordinate. 
        pressureMatDiffPos[pressureMatDiffPos < 0] = 1.e10
        upperIndex = np.argmin(pressureMatDiffPos, axis=1)  #upperIndex is 3D array in time, lat, lon
        # Next: Record indices where we're trying to interpolate to greater than the maximum native pressure in the column.
        # If this is the case, every value of pressureMatDiffPos in the column will have been set to 1.e10, so the 
        # difference between the max and min in the column will be zero.
        # We'll create an array of boolean values that are true if this is one such column. They'll be used
        # later to slot in nans.
        nanUpperIndex = np.max(pressureMatDiffPos, axis=1) - np.min(pressureMatDiffPos,axis=1) #3D array: time, lat, lon
        nanUpperIndexBool = np.ones(np.shape(nanUpperIndex), dtype = bool)
        nanUpperIndexBool[nanUpperIndex >= 1.e-99] = False
        # Now, convert the 3D arrays containing vertical indices of interest and boolean values for nans to 1D vectors
        upperIndex1D = np.swapaxes(upperIndex, 1,2) #switched lat and lon to match reshaped matrices above        
        upperIndex1D = np.reshape(upperIndex1D, numTimes*numLats*numLons) #Convert to a vector
        nanUpperIndex1D = np.swapaxes(nanUpperIndexBool, 1,2)
        nanUpperIndex1D = np.reshape(nanUpperIndex1D, numTimes*numLats*numLons)
        # Now, extract the native pressure and data values at the upper bound indices we found
        # (I tested this method for extracting data using a sample 2D array in the command line)
        upperPressureBound = pressure2D[range(0,len(upperIndex1D)),upperIndex1D] #Result: 1D vector        
        upperDataBound = data2D[range(0,len(upperIndex1D)),upperIndex1D]
        # Set the pressure and data boundary data to nans where we are trying to interpolate to outside the data range
        upperPressureBound[nanUpperIndex1D] = np.nan        
        upperDataBound[nanUpperIndex1D] = np.nan 
        #print(upperDataBound) #debug output
    
        # This code block: same as above but for the lower boundaries. (far less comments)
        print('Negative side')
        pressureMatDiffNeg = pressureMat - pLevels[k]        
        pressureMatDiffNeg[pressureMatDiffNeg > 0] = 1.e10 #This time we are only interested in negative values
        lowerIndex = np.argmin(np.abs(pressureMatDiffNeg), axis=1) 
        nanLowerIndex = np.max(pressureMatDiffNeg, axis=1) - np.min(pressureMatDiffNeg,axis=1) 
        nanLowerIndexBool = np.ones(np.shape(nanLowerIndex), dtype = bool)
        nanLowerIndexBool[nanLowerIndex >= 1.e-99] = False
        lowerIndex1D = np.swapaxes(lowerIndex, 1,2) 
        lowerIndex1D = np.reshape(lowerIndex1D, numTimes*numLats*numLons) 
        nanLowerIndex1D = np.swapaxes(nanLowerIndexBool, 1,2)
        nanLowerIndex1D = np.reshape(nanLowerIndex1D, numTimes*numLats*numLons)        
        lowerPressureBound = pressure2D[range(0,len(lowerIndex1D)),lowerIndex1D]         
        lowerDataBound = data2D[range(0,len(lowerIndex1D)),lowerIndex1D]
        lowerPressureBound[nanLowerIndex1D] = np.nan        
        lowerDataBound[nanLowerIndex1D] = np.nan
        #print(lowerDataBound) #debug output
        
        # Now: linearly interpolate the data in log pressure space 
        # (interpolating in log pressure means interpolating linearly w.r.t. height)
        interpVec = lowerDataBound+(upperDataBound-lowerDataBound)*(np.log(pLevels[k])-np.log(lowerPressureBound))/(
                                                                    np.log(upperPressureBound)-np.log(lowerPressureBound))        
        
        # Finally: Reshape to 3D matrix to put in the later 4D matrix
        interpSlice = np.reshape(interpVec, (numTimes, numLons, numLats)) #time, lon, lat for consistency with above
        interpSlice = np.swapaxes(interpSlice,1,2) # switch lat and lon again
        interpMat[:,k,:,:] = interpSlice #Populate the returned matrix
        
        print('Finished k = ' + str(k))
    
    # Now we have a 4D matrix of interpolated data in time, pressure, latitude and longitude.
    # Now take zonal and time means
    interpMatZonalMean = np.nanmean(interpMat, axis=3)
    interpMatTimeMean = np.nanmean(interpMatZonalMean, axis=0)
    
    return interpMatTimeMean
    
#Similar function for hybrid height coordinates, for the HadGEM2-ES model. 

    
    
    
# Function similar to above but without the zonal or time means, to return a 4D array that I can do other 
# things with in wrapper functions. Subset in time (to save on computation) but don't take time mean yet. 
# This time will not have model as input; instead check for which format of hybrid sigma coordinates is used.
#
# Inputs: 
# Data: NetCDF4 Dataset object containing the variable of interest. 
# varname: variable name for the variable in the Data object
# pLevels: pressure levels to interpolate to, in hPa. 
# firstYear: first year of the time subset (measured from start of run, where start year is 0)
# lastYear:  last    "  "   "    "     "    "
# DataPS: optional--NetCDF4 Dataset object containing the surface pressure output. This is 
#         optional because some of the other files ("cl", etc.) already have this. 
#
# Output:
# A 4-D Numpy array with dimensions: time, pressure, lat, lon
# and NaNs where no data exist due to terrain (i.e. over Antarctica).    
#
def convertHybridSigmaToPressureCoords(Data, varname, pLevels, firstYear, lastYear, DataPS='none'):
    # Extract the output from the objects, and do time subset
    datavar = Data.variables[varname][firstYear*12:lastYear*12+12,:,:,:]   #time, lev, lat, lon
    if DataPS == 'none':
        PS = Data.variables['ps'][firstYear*12:lastYear*12+12,:,:]  #time, lat, lon; Surface pressure (Pa)
    else:
        PS = DataPS.variables['ps'][firstYear*12:lastYear*12+12,:,:]  #time, lat, lon; Surface pressure (Pa)
    
    #See above function for more detailed documentation.
    
    # Calculate pressure at the points where the Data are located.
    # Conversion from hybrid sigma levels to pressure levels
    # based on equation: #P[i,j,k] = a[k]*p0 + b[k]*PS[i,j] 
    # But P and PS vary in time whereas a and b don't, so need to use broadcasting.    
    B = Data.variables['b'][:] #dimension: lev
    if('ap' in Data.variables.keys()): #CanESM or IPSL or MPI model
        AkP0 = Data.variables['ap'][:] #Pa; dimension: lev
    else: #BNU or CCSM4 or MIROC model (as of June 27th, 2016) #June 1, 2017: CESM and NorESM also this way
        #P0 = Data.variables['p0'][:] #Reference pressure (Pa)
        P0 = 100000 #Reference pressure (Pa). It's 100000 in all the models, but CCSM and CESM mistakenly used hPa.
        A = Data.variables['a'][:]    #Dimension: lev
        AkP0 = A*P0
    BkPSij = B[None,:,None,None]*PS[:,None,:,:] #Result: 4D array (time, lev, lat, lon)
    pressureMat = (AkP0[None,:,None,None] + BkPSij)/100. #Divide by 100 to go from Pa to hPa
    
    #Reshape output to 2D array    
    s = np.shape(pressureMat)
    numTimes = s[0]
    numLevs = s[1]
    numLats = s[2]
    numLons = s[3]    
    #Make vertical level the last axis so that columns remain intact when matrices are reshaped    
    pressure2D = np.swapaxes(pressureMat,1,3) #Now it's time, lon, lat, level
    pressure2D = np.reshape(pressure2D, (numTimes*numLats*numLons, numLevs))
    data2D = np.swapaxes(datavar,1,3)
    data2D = np.reshape(data2D, (numTimes*numLats*numLons, numLevs))
    
    # Preallocate array to hold interpolated data
    interpMat = np.empty((pressureMat.shape[0], len(pLevels), pressureMat.shape[2], pressureMat.shape[3]))
    
    # Now: loop over the desired pressure levels and interpolate data to each one
    for k in range(0, len(pLevels)):
        print('Interpolating to ' + str(pLevels[k]) + ' hPa')
        
        # This code block: find the upper boundaries of the native pressure and data for interpolation
        print('Positive side')
        pressureMatDiffPos = pressureMat - pLevels[k] #Result: 4D array of differences between native and desired pressure.
        pressureMatDiffPos[pressureMatDiffPos < 0] = 1.e10
        upperIndex = np.argmin(pressureMatDiffPos, axis=1)  #upperIndex is 3D array in time, lat, lon
        nanUpperIndex = np.max(pressureMatDiffPos, axis=1) - np.min(pressureMatDiffPos,axis=1) #3D array: time, lat, lon
        nanUpperIndexBool = np.ones(np.shape(nanUpperIndex), dtype = bool)
        nanUpperIndexBool[nanUpperIndex >= 1.e-99] = False
        upperIndex1D = np.swapaxes(upperIndex, 1,2) #switched lat and lon to match reshaped matrices above        
        upperIndex1D = np.reshape(upperIndex1D, numTimes*numLats*numLons) #Convert to a vector
        nanUpperIndex1D = np.swapaxes(nanUpperIndexBool, 1,2)
        nanUpperIndex1D = np.reshape(nanUpperIndex1D, numTimes*numLats*numLons)
        upperPressureBound = pressure2D[range(0,len(upperIndex1D)),upperIndex1D] #Result: 1D vector        
        upperDataBound = data2D[range(0,len(upperIndex1D)),upperIndex1D]
        # Set the pressure and data boundary data to nans where we are trying to interpolate to outside the data range
        upperPressureBound[nanUpperIndex1D] = np.nan        
        upperDataBound[nanUpperIndex1D] = np.nan 
    
        # This code block: same as above but for the lower boundaries. 
        print('Negative side')
        pressureMatDiffNeg = pressureMat - pLevels[k]        
        pressureMatDiffNeg[pressureMatDiffNeg > 0] = 1.e10 #This time we are only interested in negative values
        lowerIndex = np.argmin(np.abs(pressureMatDiffNeg), axis=1) 
        nanLowerIndex = np.max(pressureMatDiffNeg, axis=1) - np.min(pressureMatDiffNeg,axis=1) 
        nanLowerIndexBool = np.ones(np.shape(nanLowerIndex), dtype = bool)
        nanLowerIndexBool[nanLowerIndex >= 1.e-99] = False
        lowerIndex1D = np.swapaxes(lowerIndex, 1,2) 
        lowerIndex1D = np.reshape(lowerIndex1D, numTimes*numLats*numLons) 
        nanLowerIndex1D = np.swapaxes(nanLowerIndexBool, 1,2)
        nanLowerIndex1D = np.reshape(nanLowerIndex1D, numTimes*numLats*numLons)        
        lowerPressureBound = pressure2D[range(0,len(lowerIndex1D)),lowerIndex1D]         
        lowerDataBound = data2D[range(0,len(lowerIndex1D)),lowerIndex1D]
        lowerPressureBound[nanLowerIndex1D] = np.nan        
        lowerDataBound[nanLowerIndex1D] = np.nan
        
        # Now: linearly interpolate the data in log pressure space 
        # (interpolating in log pressure means interpolating linearly w.r.t. height)
        interpVec = lowerDataBound+(upperDataBound-lowerDataBound)*(np.log(pLevels[k])-np.log(lowerPressureBound))/(
                                                                    np.log(upperPressureBound)-np.log(lowerPressureBound))        
        
        # Finally: Reshape to 3D matrix to put in the later 4D matrix
        interpSlice = np.reshape(interpVec, (numTimes, numLons, numLats)) #time, lon, lat for consistency with above
        interpSlice = np.swapaxes(interpSlice,1,2) # switch lat and lon again
        interpMat[:,k,:,:] = interpSlice #Populate the returned matrix
        
        print('Finished k = ' + str(k))
    
    return interpMat
    
# Function to obtain the low, middle and high cloud fraction from the "cl" output, based on output from above,
# averaged over time. Assumes time subset has already been done.
# Assumes random overlap (i.e. clouds in adjacent layers independent).
# Do overlap calcualtion before time mean. 
# Returns 3 2D arrays for low, middle and high cloud fraction, in a dict with keys 'low', 'middle', 'high'.
# Input: "clPressure": 4D array (time, pressure, lat, lon) of cloud fractions, e.g. obtained from above. 
#        "pLevels": the pressure levels at which the cloud fractions are specified
#        "divide100": Whetehr to divide cloud fraction by 100 first. This is true by default, but 2 models
#                     (BNU and CCSM4) didn't follow CMIP5 conventions and reported fraction instead of percent.        
# monthWeights: 12-member array of weights intended for use with "np.ma.average", associated with each month. 
#               Intended for seasonal means, e.g. DJF: [1,1,0,0,0,0,0,0,0,0,0,1]. 
#               Counted from first month in the "clPressure" array, so if first month is not January, need
#               to adjust input weights array. By default, not used and annual mean taken.
def calcLowMiddleHighClouds(clPressure, pLevels, divide100=True, monthWeights='none'):
    if divide100:
        clPressure = clPressure / 100.
    #Subset to low, middle and high 
    clLow = clPressure[:, pLevels >= 680, :, :]
    clMiddle = clPressure[:, np.logical_and(pLevels < 680, pLevels >= 440), :, :]
    clHigh = clPressure[:, pLevels < 440, :, :]
    #Calculate cloud fractions in low, middle and high clouds using random overlap assumption. 
    #The probability of clouds in any layer is one minus the probability of clouds in no layers. 
    #The probability of clouds in no layer is the product of 1 - cloud fraction in each layer, if 
    #the layers are independent.
    oneMinus_clLow = 1.0 - clLow
    oneMinus_clMiddle = 1.0 - clMiddle 
    oneMinus_clHigh = 1.0 - clHigh
    #Set nans to 1.0 before taking product
    oneMinus_clLow[np.isnan(oneMinus_clLow)] = 1.0
    oneMinus_clMiddle[np.isnan(oneMinus_clMiddle)] = 1.0
    oneMinus_clHigh[np.isnan(oneMinus_clHigh)] = 1.0
    #Do the multiplication
    clRandomLow = 1 - np.prod(oneMinus_clLow, axis=1) #should now have dimensions time, lat, lon
    #print(np.shape(clRandomLow)) #debug
    clRandomMiddle = 1 - np.prod(oneMinus_clMiddle, axis=1)
    clRandomHigh = 1 - np.prod(oneMinus_clHigh, axis=1)
    #Create dict to return, and do time averaging
    returnDict = dict()
    if monthWeights == 'none': #annual mean
        returnDict['low'] = np.mean(clRandomLow, axis=0)
        returnDict['middle'] = np.mean(clRandomMiddle, axis=0)
        returnDict['high'] = np.mean(clRandomHigh, axis=0)
    else: #seasonal mean
        tiledWeights = np.tile(monthWeights,np.shape(clRandomLow)[0]/12)
        returnDict['low'] = np.average(clRandomLow, weights=tiledWeights, axis=0)
        returnDict['middle'] = np.average(clRandomMiddle, weights=tiledWeights, axis=0)
        returnDict['high'] = np.average(clRandomHigh, weights=tiledWeights, axis=0)
    return returnDict
    
#Version of the above for models in height coordinates (i.e. HadGEM2-ES)
#Boundaries: 3250, 6500 m, corresponding to 440 and 680 hPa using standard atmosphere.
def calcLowMiddleHighCloudsHeight(clHeight, zLevels, divide100=True, monthWeights='none'):
    if divide100: 
        clHeight = clHeight / 100. 
    #subset to low, middle and high
    clLow = clHeight[:, zLevels <= 3250, :, :]
    clMiddle = clHeight[:, np.logical_and(zLevels > 3250, zLevels <= 6500), :, :]
    clHigh = clHeight[:, zLevels > 6500, :, :]
    #Below this line, same as the pressure version
    oneMinus_clLow = 1.0 - clLow
    oneMinus_clMiddle = 1.0 - clMiddle 
    oneMinus_clHigh = 1.0 - clHigh
    oneMinus_clLow[np.isnan(oneMinus_clLow)] = 1.0
    oneMinus_clMiddle[np.isnan(oneMinus_clMiddle)] = 1.0
    oneMinus_clHigh[np.isnan(oneMinus_clHigh)] = 1.0
    clRandomLow = 1 - np.prod(oneMinus_clLow, axis=1) 
    clRandomMiddle = 1 - np.prod(oneMinus_clMiddle, axis=1)
    clRandomHigh = 1 - np.prod(oneMinus_clHigh, axis=1)
    returnDict = dict()
    if monthWeights == 'none': #annual mean
        returnDict['low'] = np.mean(clRandomLow, axis=0)
        returnDict['middle'] = np.mean(clRandomMiddle, axis=0)
        returnDict['high'] = np.mean(clRandomHigh, axis=0)
    else: #seasonal mean
        tiledWeights = np.tile(monthWeights,np.shape(clRandomLow)[0]/12)
        returnDict['low'] = np.average(clRandomLow, weights=tiledWeights, axis=0)
        returnDict['middle'] = np.average(clRandomMiddle, weights=tiledWeights, axis=0)
        returnDict['high'] = np.average(clRandomHigh, weights=tiledWeights, axis=0)
    return returnDict

    
    
# Function to calculate Estimated Inversion Strength (Wood & Bretherton, 2006), 
# averaged over multiple years (but not doing any spatial averaging).
# Results only valid over the ocean. 
#
# Inputs:
# Data:      A dictionary of NetCDF Dataset objects containing the following variables, with keys as follows:
#            'tas': Near-surface air temperature
#            'ta': Air temperature at various heights
#            'huss': Near-surface specific humidity 
#            'ps': Surface pressure
# firstYear: First year of the multi-year averaging, measured from the start of the file starting at 0
# lastYear:  Last  year of the multi-year averaging, measured from the start of the file starting at 0
# index700:  Index of the 700 hPa pressure level, starting from 0. For standard CMIP5 pressure levels in 
#            order of decreasing pressure, this is 3.      
# monthWeights: see "calcLowMiddleHighClouds" documentation
#    
# Returns: a 2D array (lat, lon dimensions)
def multiYearMeanEIS(Data, firstYear, lastYear, index700=3, monthWeights='none'):
    #Define physical constants    
    Rd = 287.  # J kg^-1 K^-1
    Rv = 461.  # J kg^-1 K^-1
    Cp = 1004. # J kg^-1 K^-1
    Lv = 2.5e6 # J kg^-1
    g = 9.81   # m s^-2
    
    #Define the temperature, pressure and specific humidity variables 
    #(Internal variable names holdovers from CESM analysis script)
    #Also do time subset now, so not doing calculations for entire dataset, only the years needed.
    T700 = Data['ta'].variables['ta'][firstYear*12:lastYear*12+12,index700,:,:] #Temperature at 700 hPa in K (time, lat, lon)
    T0 = Data['tas'].variables['tas'][firstYear*12:lastYear*12+12, :,:] #Surface air temperature in K
    PS = Data['ps'].variables['ps'][firstYear*12:lastYear*12+12, :,:]/100. #Divide by 100 to go from Pa (CMIP5 standard) to hPa
    QS = Data['huss'].variables['huss'][firstYear*12:lastYear*12+12, :,:]  #Surface specific humidity, as a fraction (kg/kg)
    
    #Calculate the LTS (part of the EIS)
    Theta0 = T0*np.power((1000./PS),(Rd/Cp))
    Theta700 = T700*np.power((1000./700.),(Rd/Cp))
    LTS_3D = Theta700 - Theta0 #3D in the sense of time, lat, lon
    
    #Calculate the height of the Lifting Condensation Level (LCL)
    partialPressureWaterVapor = QS*(29./18.) 
    e = partialPressureWaterVapor*PS #Vapor pressure in hPa
    Td = 243.5*np.log(e/6.112)/(17.67-np.log(e/6.112)) #Dew point, Bolton, 1980, eq. 11 (rearr. to more commonly used form)
    Td = Td + 273.15     #Td was in degrees C, needed to convert to K.
    LCL = (T0 - Td)/.008 #Lamb & Verlinde, 2011, eq. 6.16  (The denominator is in deg C / m)
    
    #Calculate the height of the 700 hPa surface
    z700 = (Rd*T0/g)*np.log(PS/700.) #Wood & Bretherton, 2006, eq. 6
    
    #Calculate the moist-adiabatic potential temperature gradient
    T850 = (T0 + T700)/2.                                 #850 hPa temperature
    T850degC = T850 - 273.15                              #Convert to degrees C for following equation
    es850 = 6.112*np.exp(17.67*T850degC/(T850degC+243.5)) #Saturation vapor pressure, Bolton, 1980, eq. 10
    Qsat850 = (es850/850.)*(18./29.)                      #Saturation specific humidity at 850 hPa
    Gamma_m850 = (g/Cp)*(1-(1+Lv*Qsat850/(Rd*T850))/(1+Lv*Lv*Qsat850/(Cp*Rv*T850*T850))) #Wood & Bretherton, 2006, eq. 5
    
    #Put it all together to get EIS
    EIS_3D = LTS_3D - Gamma_m850*(z700-LCL) #Wood & Bretehrton, 2006, eq. 4
    
    #Now take multi-year mean (have already subsetted in time, so do this over entire time axis)
    if monthWeights == 'none':
        multiYearMeanData = np.mean(EIS_3D, axis=0)
    else:
        tiledWeights = np.tile(monthWeights,np.shape(EIS_3D)[0]/12)
        multiYearMeanData = np.average(EIS_3D, weights=tiledWeights, axis=0)
    return multiYearMeanData    
    

# Similar function to "multiYearMeanEIS" but instead returns LTS, the lower tropospheric stability 
# (difference in potential temperature between 1000 and 700 hPa)
# Same inputs as multiYearMeanEIS, except don't need the humidity stuff in the "Data" input dict.
# monthWeights: see "calcLowMiddleHighClouds" documentation
def multiYearMeanLTS(Data, firstYear, lastYear, index700=3, monthWeights='none'):
    #Define physical constants    
    Rd = 287.  # J kg^-1 K^-1
    Cp = 1004. # J kg^-1 K^-1
    
    #Define the temperature, pressure and specific humidity variables 
    #(Internal variable names holdovers from CESM analysis script)
    #Also do time subset now, so not doing calculations for entire dataset, only the years needed.
    T700 = Data['ta'].variables['ta'][firstYear*12:lastYear*12+12,index700,:,:] #Temperature at 700 hPa in K (time, lat, lon)
    T0 = Data['tas'].variables['tas'][firstYear*12:lastYear*12+12, :,:] #Surface air temperature in K
    PS = Data['ps'].variables['ps'][firstYear*12:lastYear*12+12, :,:]/100. #Divide by 100 to go from Pa (CMIP5 standard) to hPa
    
    #Calculate the LTS (part of the EIS)
    Theta0 = T0*np.power((1000./PS),(Rd/Cp))
    Theta700 = T700*np.power((1000./700.),(Rd/Cp))
    LTS_3D = Theta700 - Theta0 #3D in the sense of time, lat, lon
    
    #Now take multi-year mean (have already subsetted in time, so do this over entire time axis)
    if monthWeights == 'none':
        multiYearMeanData = np.mean(LTS_3D, axis=0)
    else:
        tiledWeights = np.tile(monthWeights,np.shape(LTS_3D)[0]/12)
        multiYearMeanData = np.average(LTS_3D, weights=tiledWeights, axis=0)
    return multiYearMeanData        
    
    
#Dictionary of colors for each model for plotting purposes
#Note the 2-digit codes for the model which I will consistently use as dictionary keys. 
modelColors = dict()
modelColors['bn'] = 'magenta'          #BNU-ESM
modelColors['ca'] = 'indigo'           #CanESM-2
modelColors['cc'] = 'cornflowerblue'   #CCSM4
modelColors['ce'] = 'blue'             #CESM-CAM5.1-FV
modelColors['cs'] = 'teal'             #CSIRO
modelColors['ec'] = 'springgreen'      #EC-EARTH-DMI
modelColors['gi'] = 'green'            #GISS-E2-R
modelColors['hc'] = 'lime'             #HadCM3
modelColors['hg'] = 'gold'             #HadGEM2-ES
modelColors['ip'] = 'darkgoldenrod'    #IPSL-CM5A-LR
modelColors['mi'] = 'orange'           #MIROC-ESM
modelColors['mp'] = 'red'              #MPI-ESM-LR
modelColors['no'] = 'maroon'           #NorESM1
modelColors['mm'] = 'k'                #Multi-Model Mean 

#On that note, need dictionary of names for the models, for plot legends. 
modelNames = dict()
modelNames['bn'] = 'BNU-ESM'
modelNames['ca'] = 'CanESM-2'
modelNames['cc'] = 'CCSM4'
modelNames['ce'] = 'CESM-CAM5.1-FV'
modelNames['cs'] = 'CSIRO-Mk3L-1-2'
modelNames['ec'] = 'EC-EARTH-DMI'
modelNames['gi'] = 'GISS-E2-R'
modelNames['hc'] = 'HadCM3'
modelNames['hg'] = 'HadGEM2-ES'
modelNames['ip'] = 'IPSL-CM5A-LR'
modelNames['mi'] = 'MIROC-ESM'
modelNames['mp'] = 'MPI-ESM-LR'
modelNames['no'] = 'NorESM1-M'
modelNames['mm'] = 'Multi-Model Mean'



#### Stuff for time series analysis ####


class MonthlyTimeSeries:

    def __init__(self, data, startMonth, startYear):
        self.data1D = data #raw data for plots of the full time series
        self.startMonth = startMonth
        self.startYear = startYear
        #Set up a 2D matrix with each column representing a month
        #while padding incomplete years at beginning and end with nans.
        if (len(data)+startMonth-1)%12 == 0:
            self.numYearsSpanned = (len(data)+startMonth-1)//12
        else:
            self.numYearsSpanned = (len(data)+startMonth-1)//12+1
        padBeginning = startMonth-1
        padEnd = 12-(len(data)+padBeginning)%12
        if padEnd == 12:
            padEnd = 0
        paddedData = np.pad(data, (padBeginning, padEnd), 
                            'constant', constant_values = (np.nan,np.nan)) 
        self.data2D = np.reshape(paddedData, (self.numYearsSpanned,12))
        self.years = np.arange(startYear, startYear+self.numYearsSpanned)

    #annual mean (incomplete years will be counted)
    def annualMean(self):
        return np.nanmean(self.data2D,axis=1)

    #mam (spring) mean
    def mamMean(self):
        return np.nanmean(self.data2D[:,2:5],axis=1)

    #jja (summer) mean
    def jjaMean(self):
        return np.nanmean(self.data2D[:,5:8],axis=1)

    #son (autumn) mean
    def sonMean(self):
        return np.nanmean(self.data2D[:,8:11],axis=1)

    #djf (winter) mean
    #How do I want to do this? Want to include December of previous year
    #reshape into vector, add 1 nan at beginning and delete the last one,
    #reshape back to matrix (now elements are shifted) and nanmean the 1st
    #3 columns. 
    def djfMean(self):
        vect = np.reshape(self.data2D, (self.data2D.size,1)) 
        vect=np.squeeze(vect)
        #vect = np.pad(vect, (1,0), 'constant', constant_values=(np.nan,np.nan))
        vect = np.concatenate((np.array([np.nan]),vect))
        vect = vect[0:len(vect)-1] #Delete last month so still divisible by 12
        mat = np.reshape(vect, (self.numYearsSpanned,12))
        return np.nanmean(mat[:,0:3],axis=1)


# Function to calculate the global mean of a variable
# from the GCM output. Now using the 
# latitudes taken from the netCDF file itself.
# This version is for all times (returns an array with the 
# means for every month). 
def globalMean(Data, varname): #"data" is netCDF4 data object
    latitudes = Data.variables['lat'][:]
    areaWeights = np.cos(latitudes*np.pi/180)
    datavar = np.squeeze(Data.variables[varname][:]) #Squeeze in case this is 1 slice from a height-varying array
    areaWeights3D = np.swapaxes(np.tile(areaWeights,
                                (np.shape(datavar)[0],np.shape(datavar)[2],1)),
                                1,2)  #Replicating the area weights 
    weighted3DMatrix = datavar*areaWeights3D
    sumWeighted = np.sum(weighted3DMatrix,axis=(1,2))
    sumWeights3D = np.sum(areaWeights3D,axis=(1,2))
    weightedMean = sumWeighted/sumWeights3D
    return weightedMean




####  Alternative versions of various functions which were used for testing  ####

#Northward energy flux functions in which time mean was done after integration 
#(didn't matter, and this version is slower so don't use)
def AtmosEnergyVerticalFluxConvergenceV2(DataDict, firstYear, lastYear):
    #Calculate the energy flux into the column
    fluxConv = (  DataDict['rsdt'].variables['rsdt'][firstYear*12:lastYear*12+12,:,:] 
                - DataDict['rsut'].variables['rsut'][firstYear*12:lastYear*12+12,:,:]
                - DataDict['rlut'].variables['rlut'][firstYear*12:lastYear*12+12,:,:]
                - DataDict['rsds'].variables['rsds'][firstYear*12:lastYear*12+12,:,:]
                - DataDict['rlds'].variables['rlds'][firstYear*12:lastYear*12+12,:,:]
                + DataDict['rsus'].variables['rsus'][firstYear*12:lastYear*12+12,:,:]
                + DataDict['rlus'].variables['rlus'][firstYear*12:lastYear*12+12,:,:]
                + DataDict['hfls'].variables['hfls'][firstYear*12:lastYear*12+12,:,:]
                + DataDict['hfss'].variables['hfss'][firstYear*12:lastYear*12+12,:,:])
    return fluxConv

def AtmosEnergyTransportNorthwardV2(fluxConv, lat, lon, a):
    #First calculate grid cell area
    latDiff = lat[1]-lat[0]
    lonDiff = lon[1]-lon[0]
    gridCellAreas = a*a*(latDiff*np.pi/180.)*(lonDiff*np.pi/180.)*np.cos(lat*np.pi/180.) #Vector based on latitude
        
    #Now calculate total energy flux in grid boxes in Watts (not W/m^2)
    fluxWatts = fluxConv*gridCellAreas[None, :,None] #Dimensions of fluxConv are time, lat,lon;
                                               #"None" ensures weighting by lat, not lon, even if square grid.

    #Now integrate over longitudes
    fluxWattsZonalSum = np.sum(fluxWatts,axis=2) #Result: 2D array, dimensions: time, lat
    #Now cumuliatively integrate over latitudes
    energyTransportNorthward = np.cumsum(fluxWattsZonalSum, axis=1) #result: 2D array, dimensions: time, lat
    #Finally, take time mean
    energyTransportNorthwardMultiYearMean = np.mean(energyTransportNorthward,axis=0)
    return energyTransportNorthwardMultiYearMean
