# -*- coding: utf-8 -*-
"""
Created on Mon May  6 10:30:55 2019

@author: Tristan O'Hanlon

Convert water content phase data (g/m^3) to specific phase content (kg/kg)
"""

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = (pressure[:,1] * 100) / (286.9 * temp[:,1]) * 1000
air_density = air_density[1:138]

iwc = iw[:,1] / air_density