# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This import the reduced datasets from 2006 - 2011 CCCM. 
The code can select either global or southern ocean data.

"""
import time
import matplotlib.pyplot as plt
import h5py
import os

start = time.time()


#---Importing Data from Reduced Datasets---#

# Uni Laptop

#CCCM Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
a = h5py.File('2006_CCCM.h5', 'r')
b = h5py.File('2007_CCCM.h5', 'r')
c = h5py.File('2008_CCCM.h5', 'r')
d = h5py.File('2009_CCCM.h5', 'r')
e = h5py.File('2010_CCCM.h5', 'r')
f = h5py.File('2011_CCCM.h5', 'r')


"""
# Home PC
#ECMWF Data


#CCCM Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('2010_CCCM.h5', 'r')


"""



############################################################################### 2006

#---CCCM Global Latitude Data---#

a_cccm_tcc_lat_g = a['tcc'][:] # 0-1
a_cccm_tclw_lat_g = a['tclw'][:] #kg/kg
a_cccm_tciw_lat_g = a['tciw'][:] #kg/kg

#---CCCM Southern Ocean Latitude Data---#

a_cccm_tcc_lat_so = a_cccm_tcc_lat_g[a_cccm_tcc_lat_g[:,0]>=-70]
a_cccm_tcc_lat_so = a_cccm_tcc_lat_so[a_cccm_tcc_lat_so[:,0]<=-50] # 0-1

a_cccm_tclw_lat_so = a_cccm_tclw_lat_g[a_cccm_tclw_lat_g[:,0]>=-70]
a_cccm_tclw_lat_so = a_cccm_tclw_lat_so[a_cccm_tclw_lat_so[:,0]<=-50] # kg/kg

a_cccm_tciw_lat_so = a_cccm_tciw_lat_g[a_cccm_tciw_lat_g[:,0]>=-70]
a_cccm_tciw_lat_so = a_cccm_tciw_lat_so[a_cccm_tciw_lat_so[:,0]<=-50] # kg/kg

#---CCCM Global Profile---#

a_cccm_tcc_alt_g = a['cf'][12:109] # 0-1
a_cccm_tclw_alt_g = a['lw'][70:133] #kg/kg
a_cccm_tciw_alt_g = a['iw'][36:133] #kg/kg
a_cccm_temp_alt_g = a['temp'][:] #K
a_cccm_plevel_alt_g = a['pressure'][:] #hPa

a_cccm_tcc_temp_g = a['cf_t'][:] # 0-1
a_cccm_tclw_temp_g = a['lw_t'][:] #kg/kg
a_cccm_tciw_temp_g = a['iw_t'][:] #kg/kg


#---CCCM Southern Ocean Profile---#

a_cccm_tcc_alt_so = a['cf_so'][:] # 0-1
a_cccm_tclw_alt_so = a['lw_so'][53:] #kg/kg
a_cccm_tciw_alt_so = a['iw_so'][22:] #kg/kg
a_cccm_temp_alt_so = a['temp_so'][:] #K
a_cccm_plevel_alt_so = a['pressure_so'][:] #hPa

a_cccm_tcc_temp_so = a['cf_t_so'][:] # 0-1
a_cccm_tclw_temp_so = a['lw_t_so'][:] #kg/kg
a_cccm_tciw_temp_so = a['iw_t_so'][:] #kg/kg


############################################################################### 2007

#---CCCM Global Latitude Data---#

b_cccm_tcc_lat_g = b['tcc'][:] # 0-1
b_cccm_tclw_lat_g = b['tclw'][:] #kg/kg
b_cccm_tciw_lat_g = b['tciw'][:] #kg/kg

#---CCCM Southern Ocean Latitude Data---#

b_cccm_tcc_lat_so = b_cccm_tcc_lat_g[b_cccm_tcc_lat_g[:,0]>=-70]
b_cccm_tcc_lat_so = b_cccm_tcc_lat_so[b_cccm_tcc_lat_so[:,0]<=-50] # 0-1

b_cccm_tclw_lat_so = b_cccm_tclw_lat_g[b_cccm_tclw_lat_g[:,0]>=-70]
b_cccm_tclw_lat_so = b_cccm_tclw_lat_so[b_cccm_tclw_lat_so[:,0]<=-50] # kg/kg

b_cccm_tciw_lat_so = b_cccm_tciw_lat_g[b_cccm_tciw_lat_g[:,0]>=-70]
b_cccm_tciw_lat_so = b_cccm_tciw_lat_so[b_cccm_tciw_lat_so[:,0]<=-50] # kg/kg

#---CCCM Global Profile---#

b_cccm_tcc_alt_g = b['cf'][12:109] # 0-1
b_cccm_tclw_alt_g = b['lw'][70:133] #kg/kg
b_cccm_tciw_alt_g = b['iw'][36:133] #kg/kg
b_cccm_temp_alt_g = b['temp'][:] #K
b_cccm_plevel_alt_g = b['pressure'][:] #hPa

b_cccm_tcc_temp_g = b['cf_t'][:] # 0-1
b_cccm_tclw_temp_g = b['lw_t'][:] #kg/kg
b_cccm_tciw_temp_g = b['iw_t'][:] #kg/kg


#---CCCM Southern Ocean Profile---#

b_cccm_tcc_alt_so = b['cf_so'][:] # 0-1
b_cccm_tclw_alt_so = b['lw_so'][53:] #kg/kg
b_cccm_tciw_alt_so = b['iw_so'][22:] #kg/kg
b_cccm_temp_alt_so = b['temp_so'][:] #K
b_cccm_plevel_alt_so = b['pressure_so'][:] #hPa

b_cccm_tcc_temp_so = b['cf_t_so'][:] # 0-1
b_cccm_tclw_temp_so = b['lw_t_so'][:] #kg/kg
b_cccm_tciw_temp_so = b['iw_t_so'][:] #kg/kg

############################################################################### 2008

#---CCCM Global Latitude Data---#

c_cccm_tcc_lat_g = c['tcc'][:] # 0-1
c_cccm_tclw_lat_g = c['tclw'][:] #kg/kg
c_cccm_tciw_lat_g = c['tciw'][:] #kg/kg

#---CCCM Southern Ocean Latitude Data---#

c_cccm_tcc_lat_so = c_cccm_tcc_lat_g[c_cccm_tcc_lat_g[:,0]>=-70]
c_cccm_tcc_lat_so = c_cccm_tcc_lat_so[c_cccm_tcc_lat_so[:,0]<=-50] # 0-1

c_cccm_tclw_lat_so = c_cccm_tclw_lat_g[c_cccm_tclw_lat_g[:,0]>=-70]
c_cccm_tclw_lat_so = c_cccm_tclw_lat_so[c_cccm_tclw_lat_so[:,0]<=-50] # kg/kg

c_cccm_tciw_lat_so = c_cccm_tciw_lat_g[c_cccm_tciw_lat_g[:,0]>=-70]
c_cccm_tciw_lat_so = c_cccm_tciw_lat_so[c_cccm_tciw_lat_so[:,0]<=-50] # kg/kg

#---CCCM Global Profile---#

c_cccm_tcc_alt_g = c['cf'][12:109] # 0-1
c_cccm_tclw_alt_g = c['lw'][70:133] #kg/kg
c_cccm_tciw_alt_g = c['iw'][36:133] #kg/kg
c_cccm_temp_alt_g = c['temp'][:] #K
c_cccm_plevel_alt_g = c['pressure'][:] #hPa

c_cccm_tcc_temp_g = c['cf_t'][:] # 0-1
c_cccm_tclw_temp_g = c['lw_t'][:] #kg/kg
c_cccm_tciw_temp_g = c['iw_t'][:] #kg/kg


#---CCCM Southern Ocean Profile---#

c_cccm_tcc_alt_so = c['cf_so'][:] # 0-1
c_cccm_tclw_alt_so = c['lw_so'][53:] #kg/kg
c_cccm_tciw_alt_so = c['iw_so'][22:] #kg/kg
c_cccm_temp_alt_so = c['temp_so'][:] #K
c_cccm_plevel_alt_so = c['pressure_so'][:] #hPa

c_cccm_tcc_temp_so = c['cf_t_so'][:] # 0-1
c_cccm_tclw_temp_so = c['lw_t_so'][:] #kg/kg
c_cccm_tciw_temp_so = c['iw_t_so'][:] #kg/kg

############################################################################### 2009

#---CCCM Global Latitude Data---#

d_cccm_tcc_lat_g = d['tcc'][:] # 0-1
d_cccm_tclw_lat_g = d['tclw'][:] #kg/kg
d_cccm_tciw_lat_g = d['tciw'][:] #kg/kg

#---CCCM Southern Ocean Latitude Data---#

d_cccm_tcc_lat_so = d_cccm_tcc_lat_g[d_cccm_tcc_lat_g[:,0]>=-70]
d_cccm_tcc_lat_so = d_cccm_tcc_lat_so[d_cccm_tcc_lat_so[:,0]<=-50] # 0-1

d_cccm_tclw_lat_so = d_cccm_tclw_lat_g[d_cccm_tclw_lat_g[:,0]>=-70]
d_cccm_tclw_lat_so = d_cccm_tclw_lat_so[d_cccm_tclw_lat_so[:,0]<=-50] # kg/kg

d_cccm_tciw_lat_so = d_cccm_tciw_lat_g[d_cccm_tciw_lat_g[:,0]>=-70]
d_cccm_tciw_lat_so = d_cccm_tciw_lat_so[d_cccm_tciw_lat_so[:,0]<=-50] # kg/kg

#---CCCM Global Profile---#

d_cccm_tcc_alt_g = d['cf'][12:109] # 0-1
d_cccm_tclw_alt_g = d['lw'][70:133] #kg/kg
d_cccm_tciw_alt_g = d['iw'][36:133] #kg/kg
d_cccm_temp_alt_g = d['temp'][:] #K
d_cccm_plevel_alt_g = d['pressure'][:] #hPa

d_cccm_tcc_temp_g = d['cf_t'][:] # 0-1
d_cccm_tclw_temp_g = d['lw_t'][:] #kg/kg
d_cccm_tciw_temp_g = d['iw_t'][:] #kg/kg


#---CCCM Southern Ocean Profile---#

d_cccm_tcc_alt_so = d['cf_so'][:] # 0-1
d_cccm_tclw_alt_so = d['lw_so'][53:] #kg/kg
d_cccm_tciw_alt_so = d['iw_so'][22:] #kg/kg
d_cccm_temp_alt_so = d['temp_so'][:] #K
d_cccm_plevel_alt_so = d['pressure_so'][:] #hPa

d_cccm_tcc_temp_so = d['cf_t_so'][:] # 0-1
d_cccm_tclw_temp_so = d['lw_t_so'][:] #kg/kg
d_cccm_tciw_temp_so = d['iw_t_so'][:] #kg/kg

############################################################################### 2010

#---CCCM Global Latitude Data---#

e_cccm_tcc_lat_g = e['tcc'][:] # 0-1
e_cccm_tclw_lat_g = e['tclw'][:] #kg/kg
e_cccm_tciw_lat_g = e['tciw'][:] #kg/kg

#---CCCM Southern Ocean Latitude Data---#

e_cccm_tcc_lat_so = e_cccm_tcc_lat_g[e_cccm_tcc_lat_g[:,0]>=-70]
e_cccm_tcc_lat_so = e_cccm_tcc_lat_so[e_cccm_tcc_lat_so[:,0]<=-50] # 0-1

e_cccm_tclw_lat_so = e_cccm_tclw_lat_g[e_cccm_tclw_lat_g[:,0]>=-70]
e_cccm_tclw_lat_so = e_cccm_tclw_lat_so[e_cccm_tclw_lat_so[:,0]<=-50] # kg/kg

e_cccm_tciw_lat_so = e_cccm_tciw_lat_g[e_cccm_tciw_lat_g[:,0]>=-70]
e_cccm_tciw_lat_so = e_cccm_tciw_lat_so[e_cccm_tciw_lat_so[:,0]<=-50] # kg/kg

#---CCCM Global Profile---#

e_cccm_tcc_alt_g = e['cf'][12:109] # 0-1
e_cccm_tclw_alt_g = e['lw'][70:133] #kg/kg
e_cccm_tciw_alt_g = e['iw'][36:133] #kg/kg
e_cccm_temp_alt_g = e['temp'][:] #K
e_cccm_plevel_alt_g = e['pressure'][:] #hPa

e_cccm_tcc_temp_g = e['cf_t'][:] # 0-1
e_cccm_tclw_temp_g = e['lw_t'][:] #kg/kg
e_cccm_tciw_temp_g = e['iw_t'][:] #kg/kg


#---CCCM Southern Ocean Profile---#

e_cccm_tcc_alt_so = e['cf_so'][:] # 0-1
e_cccm_tclw_alt_so = e['lw_so'][53:] #kg/kg
e_cccm_tciw_alt_so = e['iw_so'][22:] #kg/kg
e_cccm_temp_alt_so = e['temp_so'][:] #K
e_cccm_plevel_alt_so = e['pressure_so'][:] #hPa

e_cccm_tcc_temp_so = e['cf_t_so'][:] # 0-1
e_cccm_tclw_temp_so = e['lw_t_so'][:] #kg/kg
e_cccm_tciw_temp_so = e['iw_t_so'][:] #kg/kg

############################################################################### 2011

#---CCCM Global Latitude Data---#

f_cccm_tcc_lat_g = f['tcc'][:] # 0-1
f_cccm_tclw_lat_g = f['tclw'][:] #kg/kg
f_cccm_tciw_lat_g = f['tciw'][:] #kg/kg

#---CCCM Southern Ocean Latitude Data---#

f_cccm_tcc_lat_so = f_cccm_tcc_lat_g[f_cccm_tcc_lat_g[:,0]>=-70]
f_cccm_tcc_lat_so = f_cccm_tcc_lat_so[f_cccm_tcc_lat_so[:,0]<=-50] # 0-1

f_cccm_tclw_lat_so = f_cccm_tclw_lat_g[f_cccm_tclw_lat_g[:,0]>=-70]
f_cccm_tclw_lat_so = f_cccm_tclw_lat_so[f_cccm_tclw_lat_so[:,0]<=-50] # kg/kg

f_cccm_tciw_lat_so = f_cccm_tciw_lat_g[f_cccm_tciw_lat_g[:,0]>=-70]
f_cccm_tciw_lat_so = f_cccm_tciw_lat_so[f_cccm_tciw_lat_so[:,0]<=-50] # kg/kg

#---CCCM Global Profile---#

f_cccm_tcc_alt_g = f['cf'][12:109] # 0-1
f_cccm_tclw_alt_g = f['lw'][70:133] #kg/kg
f_cccm_tciw_alt_g = f['iw'][36:133] #kg/kg
f_cccm_temp_alt_g = f['temp'][:] #K
f_cccm_plevel_alt_g = f['pressure'][:] #hPa

f_cccm_tcc_temp_g = f['cf_t'][:] # 0-1
f_cccm_tclw_temp_g = f['lw_t'][:] #kg/kg
f_cccm_tciw_temp_g = f['iw_t'][:] #kg/kg


#---CCCM Southern Ocean Profile---#

f_cccm_tcc_alt_so = f['cf_so'][:] # 0-1
f_cccm_tclw_alt_so = f['lw_so'][53:] #kg/kg
f_cccm_tciw_alt_so = f['iw_so'][22:] #kg/kg
f_cccm_temp_alt_so = f['temp_so'][:] #K
f_cccm_plevel_alt_so = f['pressure_so'][:] #hPa

f_cccm_tcc_temp_so = f['cf_t_so'][:] # 0-1
f_cccm_tclw_temp_so = f['lw_t_so'][:] #kg/kg
f_cccm_tciw_temp_so = f['iw_t_so'][:] #kg/kg

############################################################################### Global Latitude Plots

#---Plot Global Cloud Fraction with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(a_cccm_tcc_lat_g[:,0],a_cccm_tcc_lat_g[:,1], '-r', label='2006')
ax.plot(b_cccm_tcc_lat_g[:,0],b_cccm_tcc_lat_g[:,1], '-b', label='2007')
ax.plot(c_cccm_tcc_lat_g[:,0],c_cccm_tcc_lat_g[:,1], '-g', label='2008')
ax.plot(d_cccm_tcc_lat_g[:,0],d_cccm_tcc_lat_g[:,1], '-c', label='2009')
ax.plot(e_cccm_tcc_lat_g[:,0],e_cccm_tcc_lat_g[:,1], '-m', label='2010')
ax.plot(f_cccm_tcc_lat_g[:,0],f_cccm_tcc_lat_g[:,1], '-y', label='2011')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Fraction')

plt.title('Global CCCM Cloud Fraction vs Latitude')

plt.grid(True)
plt.savefig("tcc_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Global Specific Liquid Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(a_cccm_tclw_lat_g[:,0],a_cccm_tclw_lat_g[:,1], '-r', label='2006')
ax.plot(b_cccm_tclw_lat_g[:,0],b_cccm_tclw_lat_g[:,1], '-b', label='2007')
ax.plot(c_cccm_tclw_lat_g[:,0],c_cccm_tclw_lat_g[:,1], '-g', label='2008')
ax.plot(d_cccm_tclw_lat_g[:,0],d_cccm_tclw_lat_g[:,1], '-c', label='2009')
ax.plot(e_cccm_tclw_lat_g[:,0],e_cccm_tclw_lat_g[:,1], '-m', label='2010')
ax.plot(f_cccm_tclw_lat_g[:,0],f_cccm_tclw_lat_g[:,1], '-y', label='2011')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('Global CCCM Liquid Water Content vs Latitude')

plt.grid(True)
plt.savefig("tclw_lat_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Specific Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()


ax.plot(a_cccm_tciw_lat_g[:,0],a_cccm_tciw_lat_g[:,1], '--r', label='2006')
ax.plot(b_cccm_tciw_lat_g[:,0],b_cccm_tciw_lat_g[:,1], '--b', label='2007')
ax.plot(c_cccm_tciw_lat_g[:,0],c_cccm_tciw_lat_g[:,1], '--g', label='2008')
ax.plot(d_cccm_tciw_lat_g[:,0],d_cccm_tciw_lat_g[:,1], '--c', label='2009')
ax.plot(e_cccm_tciw_lat_g[:,0],e_cccm_tciw_lat_g[:,1], '--m', label='2010')
ax.plot(f_cccm_tciw_lat_g[:,0],f_cccm_tciw_lat_g[:,1], '--y', label='2011')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('Global CCCM Ice Water Content vs Latitude')

plt.grid(True)
plt.savefig("tciw_lat_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Comparison Specific Liquid and Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()
ax.plot(a_cccm_tclw_lat_g[:,0],a_cccm_tclw_lat_g[:,1], '-r', label='2006')
ax.plot(b_cccm_tclw_lat_g[:,0],b_cccm_tclw_lat_g[:,1], '-b', label='2007')
ax.plot(c_cccm_tclw_lat_g[:,0],c_cccm_tclw_lat_g[:,1], '-g', label='2008')
ax.plot(d_cccm_tclw_lat_g[:,0],d_cccm_tclw_lat_g[:,1], '-c', label='2009')
ax.plot(e_cccm_tclw_lat_g[:,0],e_cccm_tclw_lat_g[:,1], '-m', label='2010')
ax.plot(f_cccm_tclw_lat_g[:,0],f_cccm_tclw_lat_g[:,1], '-y', label='2011')
ax.plot(a_cccm_tciw_lat_g[:,0],a_cccm_tciw_lat_g[:,1], '--r', label='2006')
ax.plot(b_cccm_tciw_lat_g[:,0],b_cccm_tciw_lat_g[:,1], '--b', label='2007')
ax.plot(c_cccm_tciw_lat_g[:,0],c_cccm_tciw_lat_g[:,1], '--g', label='2008')
ax.plot(d_cccm_tciw_lat_g[:,0],d_cccm_tciw_lat_g[:,1], '--c', label='2009')
ax.plot(e_cccm_tciw_lat_g[:,0],e_cccm_tciw_lat_g[:,1], '--m', label='2010')
ax.plot(f_cccm_tciw_lat_g[:,0],f_cccm_tciw_lat_g[:,1], '--y', label='2011')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('Global CCCM Liquid and Ice Water Content vs Latitude')

plt.grid(True)
plt.savefig("tclw_tciw_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Global Cloud Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(a_cccm_tcc_alt_g[:,1],a_cccm_tcc_alt_g[:,0], '-r', label='2006')
ax.plot(b_cccm_tcc_alt_g[:,1],b_cccm_tcc_alt_g[:,0], '-b', label='2007')
ax.plot(c_cccm_tcc_alt_g[:,1],c_cccm_tcc_alt_g[:,0], '-g', label='2008')
ax.plot(d_cccm_tcc_alt_g[:,1],d_cccm_tcc_alt_g[:,0], '-c', label='2009')
ax.plot(e_cccm_tcc_alt_g[:,1],e_cccm_tcc_alt_g[:,0], '-m', label='2010')
ax.plot(f_cccm_tcc_alt_g[:,1],f_cccm_tcc_alt_g[:,0], '-y', label='2011')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Fraction')

plt.title('Global CCCM Cloud Fraction vs Altitude')

plt.grid(True)

plt.savefig("tcc_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Specific Liquid Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(a_cccm_tclw_alt_g[:,1]*10000,a_cccm_tclw_alt_g[:,0], '-r', label='2006')
ax.plot(b_cccm_tclw_alt_g[:,1]*10000,b_cccm_tclw_alt_g[:,0], '-b', label='2007')
ax.plot(c_cccm_tclw_alt_g[:,1]*10000,c_cccm_tclw_alt_g[:,0], '-g', label='2008')
ax.plot(d_cccm_tclw_alt_g[:,1]*10000,d_cccm_tclw_alt_g[:,0], '-c', label='2009')
ax.plot(e_cccm_tclw_alt_g[:,1]*10000,e_cccm_tclw_alt_g[:,0], '-m', label='2010')
ax.plot(f_cccm_tclw_alt_g[:,1]*10000,f_cccm_tclw_alt_g[:,0], '-y', label='2011')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('Global CCCM Specific Liquid Water Content vs Altitude')

plt.grid(True)
plt.savefig("tclw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Specific Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(a_cccm_tciw_alt_g[:,1]*10000,a_cccm_tciw_alt_g[:,0], '--r', label='2006')
ax.plot(b_cccm_tciw_alt_g[:,1]*10000,b_cccm_tciw_alt_g[:,0], '--b', label='2007')
ax.plot(c_cccm_tciw_alt_g[:,1]*10000,c_cccm_tciw_alt_g[:,0], '--g', label='2008')
ax.plot(d_cccm_tciw_alt_g[:,1]*10000,d_cccm_tciw_alt_g[:,0], '--c', label='2009')
ax.plot(e_cccm_tciw_alt_g[:,1]*10000,e_cccm_tciw_alt_g[:,0], '--m', label='2010')
ax.plot(f_cccm_tciw_alt_g[:,1]*10000,f_cccm_tciw_alt_g[:,0], '--y', label='2011')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('Global CCCM Specific Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("tciw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Specific Liquid and Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(a_cccm_tclw_alt_g[:,1]*10000,a_cccm_tclw_alt_g[:,0], '-r', label='2006')
ax.plot(b_cccm_tclw_alt_g[:,1]*10000,b_cccm_tclw_alt_g[:,0], '-b', label='2007')
ax.plot(c_cccm_tclw_alt_g[:,1]*10000,c_cccm_tclw_alt_g[:,0], '-g', label='2008')
ax.plot(d_cccm_tclw_alt_g[:,1]*10000,d_cccm_tclw_alt_g[:,0], '-c', label='2009')
ax.plot(e_cccm_tclw_alt_g[:,1]*10000,e_cccm_tclw_alt_g[:,0], '-m', label='2010')
ax.plot(f_cccm_tclw_alt_g[:,1]*10000,f_cccm_tclw_alt_g[:,0], '-y', label='2011')
ax.plot(a_cccm_tciw_alt_g[:,1]*10000,a_cccm_tciw_alt_g[:,0], '--r', label='2006')
ax.plot(b_cccm_tciw_alt_g[:,1]*10000,b_cccm_tciw_alt_g[:,0], '--b', label='2007')
ax.plot(c_cccm_tciw_alt_g[:,1]*10000,c_cccm_tciw_alt_g[:,0], '--g', label='2008')
ax.plot(d_cccm_tciw_alt_g[:,1]*10000,d_cccm_tciw_alt_g[:,0], '--c', label='2009')
ax.plot(e_cccm_tciw_alt_g[:,1]*10000,e_cccm_tciw_alt_g[:,0], '--m', label='2010')
ax.plot(f_cccm_tciw_alt_g[:,1]*10000,f_cccm_tciw_alt_g[:,0], '--y', label='2011')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('Global CCCM Specific Liquid and Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("tclw_tciw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""