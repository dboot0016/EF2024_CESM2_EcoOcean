#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script for contour plots of EcoOcean output:
    1. Total System Biomass (Fig. 7)
    2. Total Community Biomass (Fig. 8)
    3. Commercial biomass (Fig. 9)
    4. Microzooplankton biomass (Fig. S21)
    5. Mesozooplankton biomass (Fig. S22)
    6. Large zooplankton biomass (Fig. S23)
    
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

VARS=  'TSB','TCB','COM','nsmz','nmdz','nlgz'                         # Variables 
UNIT = '[%]','[%]','[%]','[%]','[%]','[%]'     # Units
DESCR = 'TSB','TCB','COM','nsmz','nmdz','nlgz'                 # Description used in plots                # Used color maps
FIG_NR = '7','8','9','S21','S22','S23'                                   # Figure number

# Extent of the plots
lat1 = -90
lat2 = 90
lon1 = -179-60
lon2 = 179-60   
    
#%% Load in 1x1 grid for lat and lon
data1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/'  # Location of dataset(s) 

load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')          # Area on 1 x 1 degree grid
area_gr=load_var1['areacello'][:,:].compute().squeeze()

lat = area_gr.lat
lat = lat[5:] # EcoOcean grid starts at 84.5S while CESM starts at 89.5S
lon = area_gr.lon

#%% Location CESM2 data per scenario and run 
data_snel1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/CTL_126/'
data_snel2='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/HOS_126/'
data_snel3='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/CTL_585/'
data_snel4='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/HOS_585/' 

#%% Plotting function for the contour plots
def subplot(data1,data2,data3,data4,type_i,exp,scen):
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(-60))
    ax.set_extent([lon1, lon2, lat1, lat2], ccrs.PlateCarree())
        
    ax.coastlines(resolution='110m',zorder=5)
    ax.add_feature(cfeature.LAND,zorder=4)
    
    # type_i == 0: Difference CTL
    # type_i == 1: Difference HOS
    # type_i == 2: Difference HOS - CTL
    
    if type_i == 0:
        im=plt.pcolormesh(lon,lat,np.roll(np.flipud((data2-data1)/data1*100),180),vmin=vmn0,vmax=vmx0,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        ax.set_title(var+' ('+str(scen)+'-' + str(exp)+')',fontsize=FS)
    elif type_i == 1:
        im=plt.pcolormesh(lon,lat,np.roll(np.flipud((data4-data3)/data3*100),180),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        ax.set_title(var+' ('+str(scen)+'-' + str(exp)+')',fontsize=FS)
    elif type_i == 2:
        im=plt.pcolormesh(lon,lat,np.roll(np.flipud((data4-data3)/data3*100-(data2-data1)/data1*100),180),vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
        ax.set_title(var+' ($\Delta$-'+str(exp)+')',fontsize=FS)
  
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

    ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
    ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
    ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
    ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
    ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)
    
    ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
    ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
    ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
    ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

#%% Select variable, units, description, colormap, and figure number (for saving)
for var_i in range(len(VARS)):
    # Loop to loop over the different variables
    
    var=VARS[var_i]     # Get variable
    unit=UNIT[var_i]    # Get units for plot
    descr=DESCR[var_i]  # Get description of variable
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
    VAR21=load_var2[var][:,:].compute().squeeze()

    load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
    VAR11=load_var1[var][:,:].compute().squeeze()

    load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
    VAR41=load_var2[var][:,:].compute().squeeze()
    
    load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
    VAR31=load_var1[var][:,:].compute().squeeze()
    
    load_var2 = xr.open_dataset(f'{data_snel2}/'+var+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
    VAR22=load_var2[var][:,:].compute().squeeze()

    load_var1 = xr.open_dataset(f'{data_snel1}/'+var+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
    VAR12=load_var1[var][:,:].compute().squeeze()

    load_var2 = xr.open_dataset(f'{data_snel4}/'+var+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
    VAR42=load_var2[var][:,:].compute().squeeze()
    
    load_var1 = xr.open_dataset(f'{data_snel3}/'+var+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
    VAR32=load_var1[var][:,:].compute().squeeze()
    
    #%%
    # Directory for figures
    datadir = '/Users/daan/Desktop/EF2024_CESM2_EcoOcean/Figures/'
    
    # Select bounds for color scheme
    # Difference CTL simulation (average 2095 - 2099 minus average 2016 - 2020)
    vmn0=-70
    vmx0=-vmn0
    
    # Difference HOS simulation (average 2095 - 2099 minus average 2016 - 2020)
    vmn1=-70
    vmx1=-vmn1
    
    # Difference HOS - CTL simulation (averaged 2095 - 2099)
    vmn2=-70
    vmx2=-vmn2
    
    # Plot figures 
    scen = '126' # Scenario
    exp = 'CTL' # Experiment
    type_i = 0 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11,VAR12,VAR21,VAR22,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'a.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '126' # Scenario
    exp = 'HOS' # Experiment
    type_i = 1 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11,VAR12,VAR21,VAR22,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'b.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126' # Scenario
    exp = 'CTL' # Experiment
    type_i = 2 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11,VAR12,VAR21,VAR22,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'c.png', format='png', dpi=quality,bbox_inches='tight')
        
    # Plot figures 
    scen = '585' # Scenario
    exp = 'CTL' # Experiment
    type_i = 0 # Type of plot: 0 (difference CTL), 1 (difference HOS), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31,VAR32,VAR41,VAR42,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'d.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '585' # Scenario
    exp = 'HOS' # Experiment
    type_i = 1 # Type of plot: 0 (difference CTL), 1 (difference HOS), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31,VAR32,VAR41,VAR42,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'e.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '585' # Scenario
    exp = 'CTL' # Experiment
    type_i = 2 # Type of plot: 0 (difference CTL), 1 (difference HOS), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31,VAR32,VAR41,VAR42,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'f.png', format='png', dpi=quality,bbox_inches='tight')
        