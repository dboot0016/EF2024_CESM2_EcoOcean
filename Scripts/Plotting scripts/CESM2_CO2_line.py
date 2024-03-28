#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script for line plots of CESM2 output:
    1. Atmospheric CO2 (Fig. 1a, d)
    
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
import regionmask
from cmip6_preprocessing.regionmask import merged_mask

#%% Name some variables
year_start=2015
year_end1=2100
year_end2=2100
RUN1='control' # List of runs: control, hosing_05, 
RUN2='hosing_05'
quality=300 # quality when figure is saved in dpi
a=1         # amount of years averaged over
tm=a
save_fig='yes'

A = ([2025, 2040, 2055, 2070, 2085, 2100])
LW=4 # Line width
FS=18 # Base font size

VARS='TMCO2','-'
UNIT= '[ppm]','-'
CONV= 1,1
DESCR='CO$_2$','-'

#%% Load in native ocean and 1x1 grid for regridder
data1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/'  # Location of dataset(s) 
load_var1 = xr.open_dataset(f'{data1}/area_gn.nc')          # Area on native grid
area_gn=load_var1['areacello'][:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')          # Area on 1 x 1 degree grid
area_gr=load_var1['areacello'][:,:].compute().squeeze()

#%% Location CESM2 data per scenario and run 
data_snel1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_126/'
data_snel2='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_126/'
data_snel3='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_585/'
data_snel4='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_585/'  

#%% Select variables, unit, description
n=0
var=VARS[n]
unit=UNIT[n]
conv=CONV[n]
descr=DESCR[n]
  
#%% Call on datasets (might need to change name of dataset)
load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end1)+'_0126.nc' )
VAR1_gr=load_var1[var][:,:,:].compute().squeeze()
    
load_var1 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end1)+'_0126.nc' )
VAR2_gr=load_var1[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end2)+'_0585.nc' )
VAR3_gr=load_var1[var][:,:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end2)+'_0585.nc' )
VAR4_gr=load_var1[var].compute().squeeze()

#%% Integrate globally
lat=VAR1_gr.lat
lon=VAR1_gr.lon

RE = 6.371e6  # [m] Earth radius

dy = 2*np.pi*RE*(lat[1]-lat[0]).values/360                              # grid size in y-direction
dx = 2*np.pi*RE*((lon[1]-lon[0]).values*np.cos(np.deg2rad(lat)))/360    # grid size in x-direction

VAR_lon=(VAR1_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat=(VAR_lon[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var1=VAR_lat/5.14e18*1.0e6 * 28.966 / 44.0 # Transform units from kg to ppm

VAR_lon2=(VAR2_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat2=(VAR_lon2[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var2=VAR_lat2/5.14e18*1.0e6 * 28.966 / 44.0

VAR_lon=(VAR3_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat=(VAR_lon[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var3=VAR_lat/5.14e18*1.0e6 * 28.966 / 44.0

VAR_lon2=(VAR4_gr*dx).sum(['lon'])            # Averaged over lon
VAR_lat2=(VAR_lon2[:,:]*dy).sum(['lat'])                   # Averaged over lat
plt_var4=VAR_lat2/5.14e18*1.0e6 * 28.966 / 44.0
    
#%% Take moving mean with period tm
plt_VAR1 = pd.DataFrame(plt_var1)
rolling_windows = plt_VAR1.rolling(tm)
PLT_VAR1=np.squeeze(rolling_windows.mean())

plt_VAR2 = pd.DataFrame(plt_var2)
rolling_windows = plt_VAR2.rolling(tm)
PLT_VAR2=np.squeeze(rolling_windows.mean())
    
plt_VAR1 = pd.DataFrame(plt_var3)
rolling_windows = plt_VAR1.rolling(tm)
PLT_VAR3=np.squeeze(rolling_windows.mean())

plt_VAR2 = pd.DataFrame(plt_var4)
rolling_windows = plt_VAR2.rolling(tm)
PLT_VAR4=np.squeeze(rolling_windows.mean())

#%% Plot results
t1=np.arange(2016,2100.5,1)  
t2=t1

datadir = '/Users/daan/Desktop/EF2024_CESM2_EcoOcean/Figures/'
 
fig = plt.figure(figsize=(7, 5))  
ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t2,PLT_VAR3,linewidth=LW,color='tab:blue',linestyle='solid',label='CTL-585')
ax.plot(t2,PLT_VAR4,linewidth=LW,color='tab:orange',linestyle='solid',label='HOS-585')
ax.plot(t1,PLT_VAR1,linewidth=LW,color='tab:blue',linestyle='dashed',label='CTL-126')
ax.plot(t1,PLT_VAR2[:],linewidth=LW,color='tab:orange',linestyle='dashed',label='HOS-126')

plt.grid()
plt.xlabel('Time',fontsize=FS)
plt.ylabel(descr+' ' + unit,fontsize=FS)
plt.xlim([year_start+5,year_end2])
plt.xticks(A,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-4)
plt.xlim([2020,2100])

if save_fig == 'yes':
    plt.savefig(datadir+'Figure_1a.png', format='png', dpi=quality, transparent=True,bbox_inches='tight')

fig = plt.figure(figsize=(7, 5))  
ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t2,PLT_VAR4-PLT_VAR3,linewidth=LW,color='tab:green',label='585')
ax.plot(t1,np.array(PLT_VAR2)[:]-PLT_VAR1,linewidth=LW,color='tab:green',linestyle='dashed',label='126')
plt.grid()
plt.xlabel('Time',fontsize=FS)
plt.ylabel('$\Delta$'+descr+' ' + unit,fontsize=FS)
plt.xlim([year_start+5,year_end2])
plt.xticks(A,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.xlim([2020,2100])
plt.legend(fontsize=FS-4)
plt.ylim([-4,6])

if save_fig == 'yes':
    plt.savefig(datadir+'Figure_1d.png', format='png', dpi=quality, transparent=True,bbox_inches='tight')
