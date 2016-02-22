# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:05:40 2016

@author: Charlie

Glacier code for modeling class
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

####################INITIALIZE

#material properties
rho_i = 917.
g = 9.81
a = 2.1e-16
slide_ratio = 0.05

#spatial domain
x_min = 0 #m
dx = 100 #m
x_max = 30000 #m
x = np.arange(x_min, x_max + dx, dx, dtype=np.float128)
ice_thickness = np.zeros((len(x)), dtype=np.float128)
valley_width = 3000 #m

#initial bedrock topography: plane for now
zb_max = 4000 #m
zb_slope = -0.03 #m/m
zb = np.zeros((len(x)), dtype=np.float128)
zb[:] = zb_max + zb_slope * x
surface_elev = zb + ice_thickness

#mass balance/climate
z_ELA = 3700 #m
dbdz = 0.01 #m/yr/m
b_cap = 1

#time domain
t_min = 0
dt = 0.001 #years
t_max = 800 #years
times = np.arange(t_min, t_max + dt, dt)

#plotting
t_plot = 10 #years

glacier_fig = plt.figure(figsize=(14,6)) #instantiate figure
glacier = plt.subplot()
plt.gcf().subplots_adjust(bottom=0.20)
plt.xlabel('Distance [km]')
plt.ylabel('Elevation [m]')
plt.ion()
plt.show()

edge_ice_thickness = np.zeros((len(x) - 1), dtype=np.float128)
edge_ice_surface_slopes = np.zeros((len(x) - 1), dtype=np.float128)
########################RUN
it = 0 #iteration counter
for t in range(len(times) - 1):
    it += 1
    current_time = times[t]
    #mass balance
    b = dbdz * (surface_elev-z_ELA)
    b[:] = b.clip(max = b_cap) #cap mass balance
    #print ice_thickness
    #interpolate ice thickness at cell edges
    edge_ice_thickness[:] = (ice_thickness[1:] + ice_thickness[0:-1]) / 2 #arithmetic mean
    edge_ice_surface_slopes[:] = -np.diff(surface_elev) / dx  
    
    #calculate ice flux
    q = np.zeros((len(x)-1), dtype=np.float128)
    q[:] = (a / 5) * np.power(rho_i * g * edge_ice_surface_slopes, 3) * np.power(edge_ice_thickness, 5)
    q = np.insert(q, 0, 0)
    q = np.append(q, 0)
    
    d_ice_thickness_dt = b - np.diff(q) / dx #conservation of mass
    ice_thickness += d_ice_thickness_dt * dt #step forwards in time
    ice_thickness[:] = ice_thickness.clip(min = 0) #eliminate negative ice thickness
    
    surface_elev = zb + ice_thickness
    
    if current_time % t_plot == 0: #plot stuff
        glacier.clear()
        glacier.plot(x / 1000, zb, color='k', linewidth = 2, label='Bedrock')
        glacier.plot(x / 1000, surface_elev, color='b', label='Ice')
        #glacier.fill_between(x / 1000, zb, surface_elev, where = surface_elev > zb, facecolor = 'b', interpolate=True)
        #self.schematic.plot((min(x_for_plotting), max(self.subducting_plate_x)), (rel_sl, rel_sl), color='b', label='Water')
        #self.schematic.fill_between(x_for_plotting, -4000, self.bedrock_height[0:len(plate_for_plot)], facecolor='.5', interpolate=True)
        #self.schematic.fill_between(self.subducting_plate_x/1000, self.bedrock_height+self.coral_height, rel_sl, where=rel_sl >self.bedrock_height+self.coral_height, facecolor='b', interpolate=True)
        #self.schematic.fill_between(x_for_plotting, self.bedrock_height[0:len(plate_for_plot)], self.bedrock_height[0:len(plate_for_plot)]+self.coral_height[0:len(plate_for_plot)], where=self.bedrock_height[0:len(plate_for_plot)]+self.coral_height[0:len(plate_for_plot)] >self.bedrock_height[0:len(plate_for_plot)], facecolor='hotpink', interpolate=True)
        glacier.set_xlim(0, x_max/1000)
        glacier.set_ylim(3000, 5000)
        plt.title('A Glacier')
        plt.xlabel('Distance [km]')
        plt.ylabel('Elevation [m]')
        plt.text(5, 3500, 'Time [yrs]: %.1f' % current_time)
        #glacier_fig.tight_layout()
        plt.pause(0.2)
    else:
        pass
########################FINALIZE