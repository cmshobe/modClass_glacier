# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:05:40 2016

@author: Charlie Shobe

Glacier code for modeling class:
-Combines glacial erosion with fluvial erosion where there is no ice.
-ELA can fluctuate sinusoidally if so chosen.
-VERY SLOW, requires a very tiny dt for stability because so nonlinear.
-Could be made fast with a combination Euler Backward/Newton-Raphson iteration.
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
dx = 500 #m
x_max = 50000 #m
x = np.arange(x_min, x_max + dx, dx, dtype=np.float128)
ice_thickness = np.zeros((len(x)), dtype=np.float128)
e_param = .005 #no idea what this should be, or its units.
fluvial_k = .0005 #fluvial erosion efficiency parameter

#initial bedrock topography: plane for now
zb_max = 4000 #m
zb_slope = -0.03 #m/m
zb = np.zeros((len(x)), dtype=np.float128)
zb[:] = zb_max + zb_slope * x
surface_elev = zb + ice_thickness

#time domain
t_min = 0
dt = 0.005 #years
t_max = 10000 #years
times = np.arange(t_min, t_max + dt, dt)

#mass balance/climate
z_ELA = 3700 #m
sigma_ELA = 0 #set to 0 to turn off oscillations, or amplitude of ELA
period = 300
z_ELA_array = z_ELA + sigma_ELA * np.cos(2 * np.pi * times / period)
dbdz = 0.01 #m/yr/m
b_cap = 2

#plotting
t_plot = 100 #years

glacier_fig = plt.figure(figsize=(14,6)) #instantiate figure
glacier = plt.subplot(211)
plt.gcf().subplots_adjust(bottom=0.20)
plt.xlabel('Distance [km]')
plt.ylabel('Elevation [m]')

discharge_profile = plt.subplot(212)
plt.xlabel('Distance [km]')
plt.ylabel('Ice Discharge [m2/yr]')
plt.ion()
plt.show()

edge_ice_thickness = np.zeros((len(x) - 1), dtype=np.float128)
edge_ice_surface_slopes = np.zeros((len(x) - 1), dtype=np.float128)
########################RUN
it = 0 #iteration counter
q_def = np.zeros((len(x)-1), dtype=np.float128)
q_slide = np.zeros((len(x)-1), dtype=np.float128)
deformation_speed = np.zeros((len(x)-1), dtype=np.float128)
for t in range(len(times) - 1):
    it += 1
    current_time = times[t]
    #mass balance
    b = dbdz * (surface_elev-z_ELA_array[t])
    b[:] = b.clip(max = b_cap) #cap mass balance
    
    #interpolate ice thickness at cell edges
    edge_ice_thickness[:] = (ice_thickness[1:] + ice_thickness[0:-1]) / 2 #arithmetic mean
    edge_ice_surface_slopes[:] = -np.diff(surface_elev) / dx  
    
    #calculate ice flux
    dens_grav_power = np.power(rho_i * g * edge_ice_surface_slopes, 3)
    edge_thick_power = np.power(edge_ice_thickness, 5)
    edge_thick_power_2 = np.power(edge_ice_thickness, 4)
    q_def[:] = (a / 5) * dens_grav_power * edge_thick_power #ice q due to internal deformation
    deformation_speed[:] = (a / 5) * dens_grav_power * edge_thick_power_2   
    q_slide[:] = (slide_ratio * deformation_speed) * edge_ice_thickness
    q = q_def + q_slide
    q = np.insert(q, 0, 0)
    q = np.append(q, 0)
    
    #erosion of the bedrock
    e = e_param * np.power(slide_ratio * deformation_speed, 2)    
    zb[:-1] -= e * dt
    
    d_ice_thickness_dt = b - np.diff(q) / dx #conservation of mass
    ice_thickness += d_ice_thickness_dt * dt #step forwards in time
    ice_thickness[:] = ice_thickness.clip(min = 0) #eliminate negative ice thickness
    
    #fluvial bedrock erosion where no ice
    fluv_domain_start = np.nonzero(ice_thickness)[0][-1] + 1 #first ice-free node
    fluv_domain = x[fluv_domain_start:]
    if len(fluv_domain) > 1: #if all covered, no fluvial erosion
        fluv_length = fluv_domain - x[fluv_domain_start]
        fluv_slope = -(zb[fluv_domain_start + 1 :] - zb[fluv_domain_start : -1]) / dx
        fluv_slope = np.append(fluv_slope, fluv_slope[-1])
        drainage_area = np.power(fluv_length / 1, 1 / .571) #hack constant assumed == 1
        fluv_erosion = fluvial_k * np.power(drainage_area, 1 / 2) * fluv_slope
        zb[fluv_domain_start:] -= fluv_erosion * dt    
    else:
        pass
    surface_elev = zb + ice_thickness
    
    if current_time % t_plot == 0: #plot stuff
        glacier.clear()
        glacier.plot(x / 1000, zb, color='k', linewidth = 2, label='Bedrock')
        glacier.plot(x / 1000, surface_elev, color='b', label='Ice')
        glacier.set_xlim(0, x_max/1000)
        glacier.set_ylim(1000, 5000)
        plt.xlabel('Distance [km]')
        glacier.set_ylabel('Elevation [m]')
        glacier.text(5, 2000, 'Time [yrs]: %.1f' % current_time)
        discharge_profile.plot(x / 1000, q[:-1])
        plt.xlabel('Distance [km]')
        plt.ylabel('Ice Discharge [m2/yr]')
        discharge_profile.set_ylim(0, 50000)
        plt.pause(0.001)
        glacier_fig.savefig('glacier'+str(current_time)+'.png')
    else:
        pass
########################FINALIZE