# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 09:47:57 2019

@author: Tristan
"""


fig, ax = plt.subplots()
ax.plot( constants.lat, (clt - cal_clt)*100, '-k' )
ax.plot( constants.lat, (clwvi - cal_clwvi)*100, '-b' )
ax.plot( constants.lat, (clivi - cal_clivi)*100, '-c' )

ax.axhline(y=0, color = 'black', linestyle='-')

ax.set_ylabel('Simulated - Observed %')
ax.set_xlabel('Latitude')
ax.set_title ('CAM5 CALIPSO Cloud Cover Bias')
plt.grid(True)
plt.show()


ax = plt.axes(projection=ccrs.Mollweide(central_longitude=180))
ax.coastlines()
p = ax.contourf(constants.lon, constants.lat, (clwvi_lat_lon - cal_clwvi_lat_lon) * 100, transform=ccrs.PlateCarree(), cmap='coolwarm', vmin = -45, vmax = 45)
cbar = plt.colorbar(p, orientation='horizontal')
ax.set_title ('CAM5 Liquid %')
cbar.set_clim(-45, 45)

cbar.set_label('CALIPSO Cloud Cover Bias %')
plt.show()
