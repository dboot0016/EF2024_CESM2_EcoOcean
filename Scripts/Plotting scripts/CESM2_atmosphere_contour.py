#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script for contour plots atmospheric fields from the CESM2 output:
    1. Surface air temperature (Fig. S1)
    2. Total precipitation (Fig. S2)
    3. Zonal wind stress (Fig. S3)

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

VARS= 'TREFHT','PRECT','TAUX'                               # Variables for SAT, precipitation and zonal wind stress
UNIT = '[$^{\circ}$C]','[mm day$^{-1}$]','[N m$^{-2}$]'     # Units
DESCR = 'SAT','Precipitation',r'$\tau_{x}$'                 # Description used in plots
CMAP = ([cm.cm.thermal,'BrBG','seismic'])                 # Used color maps
FIG_NR = 'S1','S2','S3'                                     # Figure number

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
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(-60))
    ax.set_extent([lon1, lon2, lat1, lat2], ccrs.PlateCarree())
        
    ax.coastlines(resolution='110m',zorder=5)
    
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
    
    if var == 'PRECT':
        # For total precipitation, convective and large scale precipitation are summed
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'PRECC_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARa1=load_var2['PRECC'][:,:].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'PRECC_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARb1=load_var1['PRECC'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'PRECL_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARa2=load_var2['PRECL'][:,:].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel1}/PRECL_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARb2=load_var1['PRECL'][:,:].compute().squeeze()
        
        VAR21 = (VARa1+VARa2)*1e3*86400 # Transform to mm/day
        VAR11 = (VARb1+VARb2)*1e3*86400 # Transform to mm/day
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'PRECC_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARa1=load_var2['PRECC'][:,:].compute().squeeze()
        #
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'PRECC_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARb1=load_var1['PRECC'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'PRECL_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARa2=load_var2['PRECL'][:,:].compute().squeeze()
        #
        load_var1 = xr.open_dataset(f'{data_snel3}/PRECL_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARb2=load_var1['PRECL'][:,:].compute().squeeze()
        
        VAR41 = (VARa1+VARa2)*1e3*86400 # Transform to mm/day
        VAR31 = (VARb1+VARb2)*1e3*86400 # Transform to mm/day
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'PRECC_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARa12=load_var2['PRECC'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'PRECC_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARb12=load_var1['PRECC'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'PRECL_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARa22=load_var2['PRECL'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/PRECL_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARb22=load_var1['PRECL'][:,:].compute().squeeze()
        
        VAR22 = (VARa12+VARa22)*1e3*86400 # Transform to mm/day
        VAR12 = (VARb12+VARb22)*1e3*86400 # Transform to mm/day
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'PRECC_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARa12=load_var2['PRECC'][:,:].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'PRECC_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARb12=load_var1['PRECC'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'PRECL_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARa22=load_var2['PRECL'][:,:].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/PRECL_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARb22=load_var1['PRECL'][:,:].compute().squeeze()
        
        VAR42 = (VARa12+VARa22)*1e3*86400 # Transform to mm/day
        VAR32 = (VARb12+VARb22)*1e3*86400 # Transform to mm/day
        
    else:
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

    #%% Manipulate data if necessary
    if var == 'TREFHT':
        # Transform from K to degrees C
        VAR11_gr = VAR11 - 273.16
        VAR21_gr = VAR21 - 273.16
        VAR31_gr = VAR31 - 273.16
        VAR41_gr = VAR41 - 273.16
        
        VAR12_gr = VAR12 - 273.16
        VAR22_gr = VAR22 - 273.16
        VAR32_gr = VAR32 - 273.16
        VAR42_gr = VAR42 - 273.16
        
        lat = VAR11.lat
        lon = VAR11.lon
        
    else:
        VAR11_gr = VAR11
        VAR21_gr = VAR21
        VAR31_gr = VAR31
        VAR41_gr = VAR41
        
        VAR12_gr = VAR12
        VAR22_gr = VAR22
        VAR32_gr = VAR32
        VAR42_gr = VAR42
        
        lat = VAR11.lat
        lon = VAR11.lon
        
    #%%
    # Directory for figures
    datadir = '/Users/daan/Desktop/EF2024_CESM2_EcoOcean/Figures/'
    
    # Select bounds for color scheme
    # Average 2016 - 2020
    Vmn0 = ([-10,0,-0.2])
    Vmx0 = ([35,10,0.2])
    
    # Difference CTL simulation (average 2095 - 2099 minus average 2016 - 2020)
    Vmn1 = ([-7,-5,-0.06])
    Vmx1 = ([7,5,0.06])
    
    # Difference HOS - CTL simulation (averaged 2095 - 2099)
    Vmn2 = ([-3.1,-5,-0.06])
    Vmx2 = ([3.1,5,0.06])
    
    vmn0=Vmn0[var_i]
    vmx0=Vmx0[var_i]
    
    vmn1=Vmn1[var_i]
    vmx1=Vmx1[var_i]
    
    vmn2=Vmn2[var_i]
    vmx2=Vmx2[var_i]
    
    # Plot figures 
    scen = '126' # Scenario
    exp = 'CTL' # Experiment
    type_i = 0 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11_gr,VAR12_gr,VAR22_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'a.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '126' # Scenario
    exp = 'CTL' # Experiment
    type_i = 1 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11_gr,VAR12_gr,VAR22_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'b.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '126' # Scenario
    exp = 'CTL' # Experiment
    type_i = 2 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR11_gr,VAR12_gr,VAR22_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'c.png', format='png', dpi=quality,bbox_inches='tight')
        
    # Plot figures 
    scen = '585' # Scenario
    exp = 'CTL' # Experiment
    type_i = 0 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31_gr,VAR32_gr,VAR42_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'d.png', format='png', dpi=quality,bbox_inches='tight')
        
    scen = '585' # Scenario
    exp = 'CTL' # Experiment
    type_i = 1 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31_gr,VAR32_gr,VAR42_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'e.png', format='png', dpi=quality,bbox_inches='tight')
    
    scen = '585' # Scenario
    exp = 'CTL' # Experiment
    type_i = 2 # Type of plot: 0 (average 2016-2020 - CTL), 1 (difference CTL), 2 (difference HOS - CTL)
    
    fig = plt.figure(figsize=(7, 5))
    subplot(VAR31_gr,VAR32_gr,VAR42_gr,type_i,scen,exp)
    
    if save_fig == 'yes':
        plt.savefig(datadir+'Figure_'+str(fig_nr)+'f.png', format='png', dpi=quality,bbox_inches='tight')
        

