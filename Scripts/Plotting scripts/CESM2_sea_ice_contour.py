#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script for contour plot of sea ice fraction fields from the CESM2 output:
    1. Sea ice fraction (Fig. S6)
    
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

#%%
year_start1 = 2016 # Begin average period beginning
year_start2 = 2095 # Begin average period end
year_end1 = 2020 # End average period beginning
year_end2 = 2099 # End average period end

RUN1='control' 
RUN2='hosing_05' 

quality = 300 # quality for saving figures in dp
a = 5         # amount of years averaged over
FS=20         # Base fontsize
save_fig='yes' # Put on 'yes' when figures need to be saved

VARS= 'ICEFRAC' ,'-'                        # Variables 
UNIT = '[-]','-'    # Units
DESCR = 'Ice fraction','-'                 # Description used in plots
CMAP = ([cm.cm.ice,'-'])                 # Used color maps
FIG_NR = 'S6','-'                                   # Figure number

# Extent of the plots
lat1 = -90
lat2 = 90
lon1 = -179-60
lon2 = 179-60   
    
#%% Location CESM2 data per scenario and run 
data_snel1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_126/'
data_snel2='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_126/'
data_snel3='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_585/'
data_snel4='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_585/' 

#%% Plotting function for the contour plots
def subplot(data1,data2,data3,type_i,scen,exp):
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
    
    ax.coastlines(resolution='110m',zorder=5)
    ax.add_feature(cfeature.LAND,zorder=4)
    
    # type_i == 0: Average 2016 - 2020 in absolute values
    # type_i == 1: Difference CTL
    # type_i == 2: Difference HOS - CTL
    
    if type_i == 0:
        im=plt.pcolormesh(lon,lat,data1,vmin=vmn0,vmax=vmx0,transform=ccrs.PlateCarree(),cmap=cmap)
        ax.set_title('Average 2016-2020 ('+str(exp)+'-' + str(scen)+')',fontsize=FS)
    elif type_i == 1:
        im=plt.pcolormesh(lon,lat,data2-data1,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        ax.set_title('Difference ('+str(exp)+'-' + str(scen)+')',fontsize=FS)
    elif type_i == 2:
        im=plt.pcolormesh(lon,lat,data3-data2,vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        ax.set_title('Difference HOS - CTL ('+str(scen)+')',fontsize=FS)
  
    # Colorbar specifics
    cbar=plt.colorbar(im,orientation='horizontal',)
    if type_i == 0:
        cbar.ax.set_xlabel(descr+' ' + unit, fontsize=16)
    else:
        cbar.ax.set_xlabel('$\Delta$'+descr+' ' + unit, fontsize=16)
        
    cbar.ax.set_yticks(fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.xaxis.offsetText.set_fontsize(16)
    
    # Grid lines and longitude and latitude notations
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)

#%% Load in native and 1x1 grid for regridder
data1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/'  # Location of dataset(s) 
load_var1 = xr.open_dataset(f'{data1}/area_gn.nc')          # Area on native grid
area_gn=load_var1['areacello'][:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')          # Area on 1 x 1 degree grid
area_gr=load_var1['areacello'][:,:].compute().squeeze()

#%% Prepare regridder
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(area_gn, ds_out, 'bilinear',periodic=True)

#%% Select variable, units, description, colormap, and figure number (for saving)
for var_i in range(1):
    # Loop to loop over the three different variables
    
    var=VARS[var_i]     # Get variable
    reg = 'global'
    
    unit=UNIT[var_i]    # Get units for plot
    descr=DESCR[var_i]  # Get description of variable
    cmap=CMAP[var_i]    # Get color map
    fig_nr=FIG_NR[var_i]# Get figure number

    print(var)
       
    #%% Load in data
    # VAR11: CTL-126 2016 - 2020
    # VAR12: CTL-126 2095 - 2099
    # VAR21: HOS-126 2016 - 2020
    # VAR22: HOS-126 2095 - 2099
    # VAR31: CTL-585 2016 - 2020
    # VAR32: CTL-585 2095 - 2099
    # VAR41: HOS-585 2016 - 2020
    # VAR42: HOS-585 2095 - 2099
    
    load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
    VAR21_gr=load_var2[var][:,:].compute().squeeze()

    load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
    VAR11_gr=load_var1[var][:,:].compute().squeeze()

    load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
    VAR41_gr=load_var2[var][:,:].compute().squeeze()
    
    load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
    VAR31_gr=load_var1[var][:,:].compute().squeeze()
    
    load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
    VAR22_gr=load_var2[var][:,:].compute().squeeze()

    load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
    VAR12_gr=load_var1[var][:,:].compute().squeeze()

    load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
    VAR42_gr=load_var2[var][:,:].compute().squeeze()
    
    load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
    VAR32_gr=load_var1[var][:,:].compute().squeeze()
    
    lat = VAR32_gr.lat
    lon = VAR32_gr.lon
        
    #%%
    # Directory for figures
    datadir = '/Users/daan/Desktop/EF2024_CESM2_EcoOcean/Figures/'
    
    # Select bounds for color scheme
    # Average 2016 - 2020
    vmn0=0
    vmx0=0.4
    
    # Difference CTL (average 2095-2099 minus average 2016-2020)
    vmn1=-0.65
    vmx1=-vmn1
    
    # Difference HOS - CTL (average 2095-2099)
    vmn2=-0.25
    vmx2=-vmn2
    
    # Plot figures 
    scen = '126' # Scenario
    exp = 'CTL' # Experiment
    type_i = 0 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(5.5, 5))
    subplot(VAR11_gr,VAR12_gr,VAR22_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'a.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '126' # Scenario
    exp = 'CTL' # Experiment
    type_i = 1 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(5.5, 5))
    subplot(VAR11_gr,VAR12_gr,VAR22_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'b.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126' # Scenario
    exp = 'CTL' # Experiment
    type_i = 2 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(5.5, 5))
    subplot(VAR11_gr,VAR12_gr,VAR22_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'c.png', format='png', dpi=quality,bbox_inches='tight')
        
    # Plot figures 
    scen = '585' # Scenario
    exp = 'CTL' # Experiment
    type_i = 0 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(5.5, 5))
    subplot(VAR31_gr,VAR32_gr,VAR42_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'d.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '585' # Scenario
    exp = 'CTL' # Experiment
    type_i = 1 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(5.5, 5))
    subplot(VAR31_gr,VAR32_gr,VAR42_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'e.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '585' # Scenario
    exp = 'CTL' # Experiment
    type_i = 2 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(5.5, 5))
    subplot(VAR31_gr,VAR32_gr,VAR42_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'f.png', format='png', dpi=quality,bbox_inches='tight')
        
