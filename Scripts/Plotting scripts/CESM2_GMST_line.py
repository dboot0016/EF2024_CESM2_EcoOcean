#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script for line plots of CESM2 output:
    1. Global Mean Surface Temperature (Fig. 1b, e)
    
@author: Amber A. Boot (a.a.boot@uu.nl)
"""

#%% Load in modules
import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import xesmf as xe
import cmocean as cm
import gsw

#%% Name some variables
year_start=2015
year_end1=2100
year_end2=2100
RUN1='control'
RUN2='hosing_05' # List of runs: control, hosing_05, 
quality=300 # quality when figure is saved in dpi
a=5         # amount of years averaged over
save_fig='yes'
years=np.arange(year_start,year_end2+1,1)

VARS='TREFHT','-'
UNIT='[$^{\circ}$C]','-'
CONV=1,1
DESCR='SAT','-'

#%%
data1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/'   # Location of dataset(s) 

load_var1 = xr.open_dataset(f'{data1}/area_gr.nc') # area on a 1x1 degree grid
area_gr=load_var1['areacello'][:,:].compute().squeeze()  
    
#%% Location CESM2 data per scenario and run 
data_snel1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_126/'
data_snel2='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_126/'
data_snel3='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_585/'
data_snel4='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_585/'  

#%%
n=0
#for n in range (len(VARS)-1):
var=VARS[n]
unit=UNIT[n]
conv=CONV[n]
descr=DESCR[n]

print(var)
         
#%% Load in data
load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end1)+'_0126.nc' )
VAR2=load_var2[var][:,:,:].compute().squeeze()

plt.plot(VAR2['time'])

load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end1)+'_0126.nc' )
VAR1=load_var1[var][:,:,:].compute().squeeze()

load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end2)+'_0585.nc' )
VAR4=load_var2[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end2)+'_0585.nc' )
VAR3=load_var1[var][:,:,:].compute().squeeze()

#%% From K to degrees C
VAR1_gr = VAR1 - 273.16
VAR2_gr = VAR2 - 273.16
VAR3_gr = VAR3 - 273.16
VAR4_gr = VAR4 - 273.16
            
#%% Average data globally
lat=VAR1_gr.lat
lon=VAR1_gr.lon

RE = 6.371e6  # [m] Earth radius

dy = 2*np.pi*RE*(lat[1]-lat[0]).values/360                              # grid size in y-direction
dx = 2*np.pi*RE*((lon[1]-lon[0]).values*np.cos(np.deg2rad(lat)))/360    # grid size in x-direction

VAR_lon=(VAR1_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat=(VAR_lon[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var1=VAR_lat/(area_gr.sum())

VAR_lon=(VAR2_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat=(VAR_lon[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var2=VAR_lat/(area_gr.sum())

VAR_lon=(VAR3_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat=(VAR_lon[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var3=VAR_lat/(area_gr.sum())

VAR_lon=(VAR4_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat=(VAR_lon[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var4=VAR_lat/(area_gr.sum())

#%% Moving mean
tm=a # Moving mean period is a (5 years)
Lv11 = pd.DataFrame(plt_var1)
rolling_windows = Lv11.rolling(tm)
LV11=np.squeeze(rolling_windows.mean())

Lv11 = pd.DataFrame(plt_var2)
rolling_windows = Lv11.rolling(tm)
LV12=np.squeeze(rolling_windows.mean())

Lv11 = pd.DataFrame(plt_var3)
rolling_windows = Lv11.rolling(tm)
LV13=np.squeeze(rolling_windows.mean())

Lv11 = pd.DataFrame(plt_var4)
rolling_windows = Lv11.rolling(tm)
LV14=np.squeeze(rolling_windows.mean())

#%% Plot results
datadir = '/Users/daan/Desktop/EF2024_CESM2_EcoOcean/Figures/'
A=([2025,2040,2055,2070,2085,2100])

FS=18
LW=4
t1=np.arange(2016,2100.5,1)  
t2=t1

fig = plt.figure(figsize=(7, 5))  
ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t2,LV13,linewidth=LW,color='tab:blue',linestyle='solid',label='CTL-585')
ax.plot(t2,LV14,linewidth=LW,color='tab:orange',linestyle='solid',label='HOS-585')
ax.plot(t1[:],LV11,linewidth=LW,color='tab:blue',linestyle='dashed',label='CTL-126')
ax.plot(t1,LV12,linewidth=LW,color='tab:orange',linestyle='dashed',label='HOS-126')
plt.grid()
plt.xlabel('Time',fontsize=FS)
plt.ylabel('GMST [$^{\circ}$C]',fontsize=FS)
plt.xlim([year_start+5,year_end2])
plt.xticks(A,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-4)

if save_fig == 'yes':
    plt.savefig(datadir+'Figure_1b.png', format='png', dpi=quality, transparent=True,bbox_inches='tight')

fig = plt.figure(figsize=(7, 5))  
ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t2,LV14-LV13,linewidth=LW,color='tab:green',label='585')
ax.plot(t1[:],np.array(LV12)[:]-LV11,linewidth=LW,color='tab:green',linestyle='dashed',label='126')
plt.grid()
plt.xlabel('Time',fontsize=FS)
plt.ylabel('$\Delta$GMST [$^{\circ}$C]',fontsize=FS)
plt.xlim([year_start+5,year_end2])
plt.ylim([-1.4,0.05])
plt.xticks(A,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-4)

plt.xlim([2020,2100])
if save_fig == 'yes':
    plt.savefig(datadir+'Figure_1e.png', format='png', dpi=quality, transparent=True,bbox_inches='tight')

