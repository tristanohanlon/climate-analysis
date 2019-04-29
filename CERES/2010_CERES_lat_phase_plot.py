# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: toha006
"""
#!/usr/bin/python
import matplotlib.pyplot as plt

"""
Run 2010_tcc, 2010_tciw and 2010_tclw first to generate the datasets.
"""
        
plt.figure()
fig, ax = plt.subplots()

ax.plot(ceres_tcc_lat[:,0],ceres_tcc_lat[:,1], '-r', label='Total Cloud Fraction')
ax.plot(tclw[:,0],tclw[:,1], '-b', label='Liquid Water')
ax.plot(tciw[:,0],tciw[:,1], '--b', label='Ice Water')
#ax.axis('equal')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
          
ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Fraction', color='r')
plt.title('Total Cloud Fraction and Phase vs Latitude - CERES 2010')

plt.grid(True)
plt.show()