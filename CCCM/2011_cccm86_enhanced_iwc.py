# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon


"""

import time
import numpy as np
import os
from pyhdf import SD
import h5py
from scipy import integrate
from sklearn.impute import SimpleImputer

############################################################################### Get lwc (lat,alt)

cccm86_enhanced_iwc = []

start = time.time()

# The directory where your HDF files are stored
os.chdir('d:/Downloads/CCCM/2011')  # Home PC

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)

    cccm86_enhanced_iwc = cccm86_enhanced_iwc + f.select('Ice water content profile used').get().tolist()

end = time.time()
print('Importing data from files to lists took:', end - start, 's')

############################################################################### Convert lwc to np.array

start = time.time()

cccm86_enhanced_iwc = np.array(cccm86_enhanced_iwc)
end = time.time()
print('Create arrays took:', end - start, 's')


###############################################################################

os.chdir('e:/University/University/MSc/Models/Data/CCCM/raw_data')  # Home PC

with h5py.File('2011_raw_cccm86_enhanced_iwc.h5', 'w') as p:

    p.create_dataset('cccm86_enhanced_iwc', data=cccm86_enhanced_iwc)

    p.close()
    
###############################################################################  Load raw lwc
"""
os.chdir('e:/University/University/MSc/Models/Data/CCCM/raw_data')  # Home PC
f = h5py.File('2011_raw_cccm86_enhanced_iwc.h5', 'r')

cccm86_enhanced_iwc = f['cccm86_enhanced_iwc'][:]
f.close()

os.chdir('e:/University/University/MSc/Models/Data/CCCM/raw_data')  # Home PC
f = h5py.File('2011_lat.h5', 'r')

lat = f['lat'][:]
f.close()
"""
############################################################################### Reduce altitude levels to match cloud fraction profile (113)

#reduce altitude levels to match cloud fraction profile
cccm86_enhanced_iwc = np.transpose(cccm86_enhanced_iwc)
cccm86_enhanced_iwc = cccm86_enhanced_iwc[24:]
cccm86_enhanced_iwc = np.transpose(cccm86_enhanced_iwc)

############################################################################### Load lat and alt data


os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')  # Home PC
f = h5py.File('2011_CCCM_profile_variables.h5', 'r')

#lat = f['lat'][:]
cccm123_alt = f['alt'][:]
cccm121_alt = f['alt_c'][:]
cccm124_alt = f['alt_t'][:]
f.close()


############################################################################### fill data to nan and fit all nan values to average

#---Don't do this for vaed loaded data, only raw data---#
cccm86_enhanced_iwc[cccm86_enhanced_iwc > 20] = np.nan

imp = SimpleImputer(missing_values=np.nan, strategy='mean')
imp.fit(np.transpose(cccm86_enhanced_iwc))  
a = imp.transform(np.transpose(cccm86_enhanced_iwc))
cccm86_enhanced_iwc = np.transpose(a)

############################################################################### Get Southern Ocean Data


# Join the two lists as if they were two columns side by side, into a list of two elements each
a = np.vstack(lat)
combined = np.hstack((a, cccm86_enhanced_iwc))
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
combined = combined[combined[:,0]>=-70]
combined = combined[combined[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
cccm86_enhanced_iwc_so = combined[:,1:114]
#alt = alt[0:137] #scale alt if necessary


############################################################################### Reduce to separate lat and alt data

cccm86_enhanced_iwc_lat = np.nanmean(cccm86_enhanced_iwc, axis = 1)

cccm86_enhanced_iwc_alt = np.vstack((cccm121_alt, np.nanmean(cccm86_enhanced_iwc, axis = 0))).T
cccm86_enhanced_iwc_alt_so = np.vstack((cccm121_alt, np.nanmean(cccm86_enhanced_iwc_so, axis = 0))).T


############################################################################### Create global average latitude data set

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.vstack((lat, cccm86_enhanced_iwc_lat)).T
#print ("combined")
#print (combined)

#print("get unique lats")
unique = np.unique(lat)
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

# Averages of (lat, cloud ice water content) empty array
averages_total = unique.size
cccm_tciw_lat = np.empty((averages_total,2),dtype=float)

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
        """
        print("setting current_lat to item[0]")
        print("(current_lat == item[0]) = ", end='')
        print(current_lat == item[0]) 
        """
        current_lat = item[0];
    
    # If the lat is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_lat:
        
        # Find the average value.
        average = subtotal / number
        """
        print("--------")
        print("lat: ", end='')
        print(current_lat, end='')
        print(", avg: ", end='')
        print(average, end='')
        print(", subtotal: ", end='')
        print(subtotal, end='')
        print(", number: ", end='')
        print(number)
        """
        # Append the average
        cccm_tciw_lat[i] = [current_lat, average]
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
cccm_tciw_lat[i] = [current_lat, average]

cccm86_enhanced_iwc_lat = cccm_tciw_lat

############################################################################### Save data sets

os.chdir('//synthesis/e/University/University/MSc/Models/Data/CCCM/raw_data')
with h5py.File('2011_cccm86_enhanced_iwc.h5', 'w') as p:

    p.create_dataset('cccm86_enhanced_iwc', data=cccm86_enhanced_iwc)
    p.create_dataset('cccm86_enhanced_iwc_alt', data=cccm86_enhanced_iwc_alt)
    p.create_dataset('cccm86_enhanced_iwc_lat', data=cccm86_enhanced_iwc_lat)
    p.create_dataset('cccm86_enhanced_iwc_alt_so', data=cccm86_enhanced_iwc_alt_so)

    p.close()    


    
    