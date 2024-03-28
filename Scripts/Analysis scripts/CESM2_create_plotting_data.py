#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script to create datasets for plotting CESM2 data.
Data necessary to run this script is available upon request.
This data is not put in the repository because of its size.
    
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

NAMES= 'TEMP','NO3','photoC_TOT','TEMP2','SALT','SALT2','WVEL','PO4','SiO3_150m','Fe_150m','POC_flux_100m','photoC_diat_zint','photoC_sp_zint','diazC','diatC','spC','TEMP_150m','TEMP_bot','TEMP_col','HMXL','TREFHT','TAUX','PRECC','PRECL','ICEFRAC','BSF'                         # Variables 
VARS= 'TEMP','__xarray_dataarray_variable__','photoC_TOT','TEMP2','SALT','SALT2','__xarray_dataarray_variable__','__xarray_dataarray_variable__','__xarray_dataarray_variable__','__xarray_dataarray_variable__','POC_FLUX_100m','__xarray_dataarray_variable__','__xarray_dataarray_variable__','diazC','diatC','spC','__xarray_dataarray_variable__','__xarray_dataarray_variable__','__xarray_dataarray_variable__','HMXL','TREFHT','TAUX','PRECC','PRECL','ICEFRAC','BSF'                         # Variables 
    
#%% Location where datasets should be saved to
data_snel1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_126/'
data_snel2='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_126/'
data_snel3='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_585/'
data_snel4='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_585/' 

#%% Location where original datasets are placed
data_snel5='/Users/daan/Documents/EcoOcean/CESM2 data/CTL-126/'
data_snel6='/Users/daan/Documents/EcoOcean/CESM2 data/HOS-126/'
data_snel7='/Users/daan/Documents/EcoOcean/CESM2 data/CTL-585/'
data_snel8='/Users/daan/Documents/EcoOcean/CESM2 data/HOS-585/'

#%%
for var_i in range(len(VARS)):
    var = VARS[var_i]
    name = NAMES[var_i]
    
    print(name)

    # Filenames on how to be saved
    filename_c126_1 = name + '_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc'
    filename_h126_1 = name + '_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc'
    filename_c585_1 = name + '_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc'
    filename_h585_1 = name + '_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc'
    
    filename_c126_2 = name + '_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc'
    filename_h126_2 = name + '_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc'
    filename_c585_2 = name + '_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc'
    filename_h585_2 = name + '_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc'
    
    # Load in data
    load_var2 = xr.open_dataset(f'{data_snel6}/'+name+'_yr_'+RUN2+'_'+str(2015)+'_'+str(2100)+'_0126.nc' )
    VAR2_1=load_var2[var].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
    VAR2_2=load_var2[var].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
    
    load_var1 = xr.open_dataset(f'{data_snel5}/'+name+'_yr_'+RUN1+'_'+str(2015)+'_'+str(2100)+'_0126.nc' )
    VAR1_1=load_var1[var].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
    VAR1_2=load_var1[var].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
    
    load_var2 = xr.open_dataset(f'{data_snel8}/'+name+'_yr_'+RUN2+'_'+str(2015)+'_'+str(2100)+'_0585.nc' )
    VAR4_1=load_var2[var].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
    VAR4_2=load_var2[var].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
    
    load_var1 = xr.open_dataset(f'{data_snel7}/'+name+'_yr_'+RUN1+'_'+str(2015)+'_'+str(2100)+'_0585.nc' )
    VAR3_1=load_var1[var].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
    VAR3_2=load_var1[var].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
    
    if not var == name:
        VAR1_1=VAR1_1.to_dataset(name = name)
        VAR1_2=VAR1_2.to_dataset(name = name)
        VAR2_1=VAR2_1.to_dataset(name = name)
        VAR2_2=VAR2_2.to_dataset(name = name)
        VAR3_1=VAR3_1.to_dataset(name = name)
        VAR3_2=VAR3_2.to_dataset(name = name)
        VAR4_1=VAR4_1.to_dataset(name = name)
        VAR4_2=VAR4_2.to_dataset(name = name)
        
    VAR1_1.to_netcdf(data_snel1+filename_c126_1)
    VAR1_2.to_netcdf(data_snel1+filename_c126_2)
    
    VAR2_1.to_netcdf(data_snel2+filename_h126_1)
    VAR2_2.to_netcdf(data_snel2+filename_h126_2)
    
    VAR3_1.to_netcdf(data_snel3+filename_c585_1)
    VAR3_2.to_netcdf(data_snel3+filename_c585_2)
    
    VAR4_1.to_netcdf(data_snel4+filename_h585_1)
    VAR4_2.to_netcdf(data_snel4+filename_h585_2)
    