# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

#---Import Latitude Variables---#

import h5py
import os
import numpy as np
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Uni Laptop
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Laptop


a = h5py.File('2006_CCCM.h5', 'r')
a_tcc = a['tcc'][:]
a_tclw = a['tclw'][:]
a_tciw = a['tciw'][:]
a_temp = a['temp'][:]
a_pressure = a['pressure'][:]

a_temp_so = a['temp_so'][:]
a_pressure_so = a['pressure_so'][:]

a_cf = a['cf'][:]
a_lw = a['lw'][:]
a_iw = a['iw'][:]

a_cf_so = a['cf_so'][:]
a_lw_so = a['lw_so'][:]
a_iw_so = a['iw_so'][:]

a_cf_t = a['cf_t'][:]
a_lw_t = a['lw_t'][:]
a_iw_t = a['iw_t'][:]

a_cf_t_so = a['cf_t_so'][:]
a_lw_t_so = a['lw_t_so'][:]
a_iw_t_so = a['iw_t_so'][:]

###############################################################################

b = h5py.File('2007_CCCM.h5', 'r')
b_tcc = b['tcc'][:]
b_tclw = b['tclw'][:]
b_tciw = b['tciw'][:]
b_temp = b['temp'][:]
b_pressure = b['pressure'][:]

b_temp_so = b['temp_so'][:]
b_pressure_so = b['pressure_so'][:]

b_cf = b['cf'][:]
b_lw = b['lw'][:]
b_iw = b['iw'][:]

b_cf_so = b['cf_so'][:]
b_lw_so = b['lw_so'][:]
b_iw_so = b['iw_so'][:]

b_cf_t = b['cf_t'][:]
b_lw_t = b['lw_t'][:]
b_iw_t = b['iw_t'][:]

b_cf_t_so = b['cf_t_so'][:]
b_lw_t_so = b['lw_t_so'][:]
b_iw_t_so = b['iw_t_so'][:]

###############################################################################

c = h5py.File('2008_CCCM.h5', 'r')
c_tcc = c['tcc'][:]
c_tclw = c['tclw'][:]
c_tciw = c['tciw'][:]
c_temp = c['temp'][:]
c_pressure = c['pressure'][:]

c_temp_so = c['temp_so'][:]
c_pressure_so = c['pressure_so'][:]

c_cf = c['cf'][:]
c_lw = c['lw'][:]
c_iw = c['iw'][:]

c_cf_so = c['cf_so'][:]
c_lw_so = c['lw_so'][:]
c_iw_so = c['iw_so'][:]

c_cf_t = c['cf_t'][:]
c_lw_t = c['lw_t'][:]
c_iw_t = c['iw_t'][:]

c_cf_t_so = c['cf_t_so'][:]
c_lw_t_so = c['lw_t_so'][:]
c_iw_t_so = c['iw_t_so'][:]

###############################################################################

d = h5py.File('2009_CCCM.h5', 'r')
d_tcc = d['tcc'][:]
d_tclw = d['tclw'][:]
d_tciw = d['tciw'][:]
d_temp = d['temp'][:]
d_pressure = d['pressure'][:]

d_temp_so = d['temp_so'][:]
d_pressure_so = d['pressure_so'][:]

d_cf = d['cf'][:]
d_lw = d['lw'][:]
d_iw = d['iw'][:]

d_cf_so = d['cf_so'][:]
d_lw_so = d['lw_so'][:]
d_iw_so = d['iw_so'][:]

d_cf_t = d['cf_t'][:]
d_lw_t = d['lw_t'][:]
d_iw_t = d['iw_t'][:]

d_cf_t_so = d['cf_t_so'][:]
d_lw_t_so = d['lw_t_so'][:]
d_iw_t_so = d['iw_t_so'][:]

###############################################################################

e = h5py.File('2010_CCCM.h5', 'r')
e_tcc = e['tcc'][:]
e_tclw = e['tclw'][:]
e_tciw = e['tciw'][:]
e_temp = e['temp'][:]
e_pressure = e['pressure'][:]

e_temp_so = e['temp_so'][:]
e_pressure_so = e['pressure_so'][:]

e_cf = e['cf'][:]
e_lw = e['lw'][:]
e_iw = e['iw'][:]

e_cf_so = e['cf_so'][:]
e_lw_so = e['lw_so'][:]
e_iw_so = e['iw_so'][:]

e_cf_t = e['cf_t'][:]
e_lw_t = e['lw_t'][:]
e_iw_t = e['iw_t'][:]

e_cf_t_so = e['cf_t_so'][:]
e_lw_t_so = e['lw_t_so'][:]
e_iw_t_so = e['iw_t_so'][:]

###############################################################################

f = h5py.File('2011_CCCM.h5', 'r')
f_tcc = f['tcc'][:]
f_tclw = f['tclw'][:]
f_tciw = f['tciw'][:]
f_temp = f['temp'][:]
f_pressure = f['pressure'][:]

f_temp_so = f['temp_so'][:]
f_pressure_so = f['pressure_so'][:]

f_cf = f['cf'][:]
f_lw = f['lw'][:]
f_iw = f['iw'][:]

f_cf_so = f['cf_so'][:]
f_lw_so = f['lw_so'][:]
f_iw_so = f['iw_so'][:]

f_cf_t = f['cf_t'][:]
f_lw_t = f['lw_t'][:]
f_iw_t = f['iw_t'][:]

f_cf_t_so = f['cf_t_so'][:]
f_lw_t_so = f['lw_t_so'][:]
f_iw_t_so = f['iw_t_so'][:]

############################################################################### tcc

combined = np.vstack((a_tcc, b_tcc, c_tcc, d_tcc, e_tcc, f_tcc))

#print("get unique lats")
unique = np.unique(a_tcc[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (lat, cloud ice water content) empty array
averages_total = unique.size
tcc = np.empty((averages_total,2),dtype=float)

# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud ice water content entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat, cloud ice water content) elements and subtotal the same lat values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_lat is None:
        current_lat = item[0];
    
    # If the lat is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_lat:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        tcc[i] = [current_lat, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current latitude
        current_lat = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
tcc[i] = [current_lat, average]

############################################################################### tclw

combined = np.vstack((a_tclw, b_tclw, c_tclw, d_tclw, e_tclw, f_tclw))

#print("get unique lats")
unique = np.unique(a_tclw[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (lat, cloud ice water content) empty array
averages_total = unique.size
tclw = np.empty((averages_total,2),dtype=float)

# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud ice water content entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat, cloud ice water content) elements and subtotal the same lat values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_lat is None:
        current_lat = item[0];
    
    # If the lat is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_lat:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        tclw[i] = [current_lat, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current latitude
        current_lat = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
tclw[i] = [current_lat, average]

############################################################################### tciw

combined = np.vstack((a_tciw, b_tciw, c_tciw, d_tciw, e_tciw, f_tciw))

#print("get unique lats")
unique = np.unique(a_tciw[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (lat, cloud ice water content) empty array
averages_total = unique.size
tciw = np.empty((averages_total,2),dtype=float)

# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud ice water content entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat, cloud ice water content) elements and subtotal the same lat values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_lat is None:
        current_lat = item[0];
    
    # If the lat is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_lat:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        tciw[i] = [current_lat, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current latitude
        current_lat = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
tciw[i] = [current_lat, average]

############################################################################### cf

combined = np.vstack((a_cf, b_cf, c_cf, d_cf, e_cf, f_cf))

#print("get unique alt")
unique = np.unique(a_cf[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
cf = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        cf[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
cf[i] = [current_alt, average]

############################################################################### cf_so

combined = np.vstack((a_cf_so, b_cf_so, c_cf_so, d_cf_so, e_cf_so, f_cf_so))

#print("get unique alt")
unique = np.unique(a_cf_so[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
cf_so = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        cf_so[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
cf_so[i] = [current_alt, average]

############################################################################### lw

combined = np.vstack((a_lw, b_lw, c_lw, d_lw, e_lw, f_lw))

#print("get unique alt")
unique = np.unique(a_lw[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
lw = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        lw[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
lw[i] = [current_alt, average]

############################################################################### lw_so

combined = np.vstack((a_lw_so, b_lw_so, c_lw_so, d_lw_so, e_lw_so, f_lw_so))

#print("get unique alt")
unique = np.unique(a_lw_so[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
lw_so = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        lw_so[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
lw_so[i] = [current_alt, average]

############################################################################### iw

combined = np.vstack((a_iw, b_iw, c_iw, d_iw, e_iw, f_iw))

#print("get unique alt")
unique = np.unique(a_iw[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
iw = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        iw[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
iw[i] = [current_alt, average]

############################################################################### iw_so

combined = np.vstack((a_iw_so, b_iw_so, c_iw_so, d_iw_so, e_iw_so, f_iw_so))

#print("get unique alt")
unique = np.unique(a_iw_so[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
iw_so = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        iw_so[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
iw_so[i] = [current_alt, average]

############################################################################### temp

combined = np.vstack((a_temp, b_temp, c_temp, d_temp, e_temp, f_temp))

#print("get unique alt")
unique = np.unique(a_temp[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
temp = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        temp[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
temp[i] = [current_alt, average]

############################################################################### temp_so

combined = np.vstack((a_temp_so, b_temp_so, c_temp_so, d_temp_so, e_temp_so, f_temp_so))

#print("get unique alt")
unique = np.unique(a_temp_so[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
temp_so = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        temp_so[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
temp_so[i] = [current_alt, average]

############################################################################### pressure

combined = np.vstack((a_pressure, b_pressure, c_pressure, d_pressure, e_pressure, f_pressure))

#print("get unique alt")
unique = np.unique(a_pressure[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
pressure = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        pressure[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
pressure[i] = [current_alt, average]

############################################################################### pressure_so

combined = np.vstack((a_pressure_so, b_pressure_so, c_pressure_so, d_pressure_so, e_pressure_so, f_pressure_so))

#print("get unique alt")
unique = np.unique(a_pressure_so[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (alt, cloud content) empty array
averages_total = unique.size
pressure_so = np.empty((averages_total,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction content entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt, cloud content) elements and subtotal the same alt values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

    if current_alt is None:
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        pressure_so[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
pressure_so[i] = [current_alt, average]

###############################################################################

#---Create new temperature datasets---#

cf_scale = cf[12:109]
cf_scale_so = cf_so[12:109]

t_scale_c = temp[:,1]
t_scale_c = t_scale_c[36:133]

t_scale_liw = temp[:,1]
t_scale_liw = t_scale_liw[0:137]

cf_t = np.vstack((t_scale_c, cf_scale[:,1])).T
cf_t_so = np.vstack((t_scale_c, cf_scale_so[:,1])).T
lw_t = np.vstack((t_scale_liw, lw[:,1])).T
lw_t_so = np.vstack((t_scale_liw, lw_so[:,1])).T
iw_t = np.vstack((t_scale_liw, iw[:,1])).T
iw_t_so = np.vstack((t_scale_liw, iw_so[:,1])).T


###############################################################################
    

import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Uni Laptop
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Laptop

with h5py.File('07.2006_04.2011_CCCM.h5', 'w') as p:
    
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    
    p.create_dataset('cf', data=cf)
    p.create_dataset('cf_so', data=cf_so)
    p.create_dataset('lw', data=lw)
    p.create_dataset('lw_so', data=lw_so)
    p.create_dataset('iw', data=iw)
    p.create_dataset('iw_so', data=iw_so)

    p.create_dataset('temp', data=temp)
    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure)
    p.create_dataset('pressure_so', data=pressure_so)

    p.create_dataset('cf_t', data=cf_t)
    p.create_dataset('cf_t_so', data=cf_t_so)
    p.create_dataset('lw_t', data=lw_t)
    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('iw_t', data=iw_t)
    p.create_dataset('iw_t_so', data=iw_t_so)

    p.close()






















