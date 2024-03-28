#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script to create datasets for plotting EcoOcean data.
Data necessary to run this script is available upon request.
This data is not put in the repository because of its size.
    
@author: Amber A. Boot (a.a.boot@uu.nl)
"""

#%% Import modules
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
scen = 'ctl_126'

datadir = '/Users/daan/Downloads/Ctl-126/asc/'
datadir_c1 = '/Users/daan/Downloads/Ctl-126/'
datadir_c5 = '/Users/daan/Downloads/Ctl-585/'
datadir_h1 = '/Users/daan/Downloads/Hos-126/'
datadir_h5 = '/Users/daan/Downloads/Hos-585/'

figdir_c1 = '/Users/daan/Downloads/Figures Ecoocean/CTL_126/' 
figdir_c5 = '/Users/daan/Downloads/Figures Ecoocean/CTL_585/' 
figdir_h1 = '/Users/daan/Downloads/Figures Ecoocean/HOS_126/' 
figdir_h5 = '/Users/daan/Downloads/Figures Ecoocean/HOS_585/' 
figdir_d1 = '/Users/daan/Downloads/Figures Ecoocean/DIFF_126/' 
figdir_d5 = '/Users/daan/Downloads/Figures Ecoocean/DIFF_585/' 

figdir_reg = '/Users/daan/Downloads/Figures Ecoocean/Regional line plots/'
maskdir = '/Users/daan/Downloads/Figures Ecoocean/'

figdir = '/Users/daan/Downloads/Figures/Line plots/' 

names = np.array(pd.read_csv(datadir+'var_names_ecoocean.csv',sep=';',header=None))

area_1 = np.loadtxt(f'{datadir_c1}/'+"cellarea.asc", skiprows=6)
area = xr.DataArray(area_1)  

lat = np.linspace(-84,90,175)
lon =(np.arange(0,359.5,1))

#%%
varN = np.append(names[:,1],['TSB','TCB','B10','B30','Commercial'])

#%% Load in data for all functional groups
n = 0
var_n = 'all_species'
filename = var_n + '_2015_2100.nc'

load_var1 = xr.open_dataset(f'{datadir_c1}/'+filename)
VAR_c1=load_var1[var_n].compute().squeeze() 
load_var1 = xr.open_dataset(f'{datadir_c5}/'+filename)
VAR_c5=load_var1[var_n].compute().squeeze() 
load_var1 = xr.open_dataset(f'{datadir_h1}/'+filename)
VAR_h1=load_var1[var_n].compute().squeeze() 
load_var1 = xr.open_dataset(f'{datadir_h5}/'+filename)
VAR_h5=load_var1[var_n].compute().squeeze() 
    
#%%
VAR_c1 = VAR_c1.transpose("time","dim_0","dim_1","species")
VAR_c5 = VAR_c5.transpose("time","dim_0","dim_1","species")
VAR_h1 = VAR_h1.transpose("time","dim_0","dim_1","species")
VAR_h5 = VAR_h5.transpose("time","dim_0","dim_1","species")

#%% Determine specific groups CTL-126
TSB_c1 = VAR_c1.sum(['species'])
TCB_c1 = TSB_c1 - VAR_c1[:,:,:,8] - VAR_c1[:,:,:,18] - VAR_c1[:,:,:,38] - VAR_c1[:,:,:,39]
B10_c1 = TCB_c1 - VAR_c1[:,:,:,34] - VAR_c1[:,:,:,31] - VAR_c1[:,:,:,1] - VAR_c1[:,:,:,25] - VAR_c1[:,:,:,27] - VAR_c1[:,:,:,28] - VAR_c1[:,:,:,29] - VAR_c1[:,:,:,30] - VAR_c1[:,:,:,23] - VAR_c1[:,:,:,22] - VAR_c1[:,:,:,14] - VAR_c1[:,:,:,48] - VAR_c1[:,:,:,49] - VAR_c1[:,:,:,32] - VAR_c1[:,:,:,33]
B30_c1 = B10_c1 - VAR_c1[:,:,:,24] - VAR_c1[:,:,:,24] - VAR_c1[:,:,:,4] - VAR_c1[:,:,:,7] - VAR_c1[:,:,:,11] - VAR_c1[:,:,:,13] - VAR_c1[:,:,:,17] - VAR_c1[:,:,:,37]
COM_c1 = VAR_c1[:,:,:,2] + VAR_c1[:,:,:,5] + VAR_c1[:,:,:,6] +VAR_c1[:,:,:,7] + VAR_c1[:,:,:,9] + VAR_c1[:,:,:,10] + VAR_c1[:,:,:,11] + VAR_c1[:,:,:,13] + VAR_c1[:,:,:,15] + VAR_c1[:,:,:,16] + VAR_c1[:,:,:,17] + VAR_c1[:,:,:,20] + VAR_c1[:,:,:,21] + VAR_c1[:,:,:,24] + VAR_c1[:,:,:,31] + VAR_c1[:,:,:,35] + VAR_c1[:,:,:,36] + VAR_c1[:,:,:,37]
COM_c1 = COM_c1 + VAR_c1[:,:,:,41] + VAR_c1[:,:,:,42] + VAR_c1[:,:,:,43] + VAR_c1[:,:,:,44] + VAR_c1[:,:,:,45] + VAR_c1[:,:,:,46] + VAR_c1[:,:,:,47] + VAR_c1[:,:,:,48] 

#%% Determine specific groups CTL-585
TSB_c5 = VAR_c5.sum(['species'])
TCB_c5 = TSB_c5 - VAR_c5[:,:,:,8] - VAR_c5[:,:,:,18] - VAR_c5[:,:,:,38] - VAR_c5[:,:,:,39]
B10_c5 = TCB_c5 - VAR_c5[:,:,:,34] - VAR_c5[:,:,:,31] - VAR_c5[:,:,:,1] - VAR_c5[:,:,:,25] - VAR_c5[:,:,:,27] - VAR_c5[:,:,:,28] - VAR_c5[:,:,:,29] - VAR_c5[:,:,:,30] - VAR_c5[:,:,:,23] - VAR_c5[:,:,:,22] - VAR_c5[:,:,:,14] - VAR_c5[:,:,:,48] - VAR_c5[:,:,:,49] - VAR_c5[:,:,:,32] - VAR_c5[:,:,:,33]
B30_c5 = B10_c5 - VAR_c5[:,:,:,24] - VAR_c5[:,:,:,24] - VAR_c5[:,:,:,4] - VAR_c5[:,:,:,7] - VAR_c5[:,:,:,11] - VAR_c5[:,:,:,13] - VAR_c5[:,:,:,17] - VAR_c5[:,:,:,37]
COM_c5 = VAR_c5[:,:,:,2] + VAR_c5[:,:,:,5] + VAR_c5[:,:,:,6] +VAR_c5[:,:,:,7] + VAR_c5[:,:,:,9] + VAR_c5[:,:,:,10] + VAR_c5[:,:,:,11] + VAR_c5[:,:,:,13] + VAR_c5[:,:,:,15] + VAR_c5[:,:,:,16] + VAR_c5[:,:,:,17] + VAR_c5[:,:,:,20] + VAR_c5[:,:,:,21] + VAR_c5[:,:,:,24] + VAR_c5[:,:,:,31] + VAR_c5[:,:,:,35] + VAR_c5[:,:,:,36] + VAR_c5[:,:,:,37]
COM_c5 = COM_c5 + VAR_c5[:,:,:,41] + VAR_c5[:,:,:,42] + VAR_c5[:,:,:,43] + VAR_c5[:,:,:,44] + VAR_c5[:,:,:,45] + VAR_c5[:,:,:,46] + VAR_c5[:,:,:,47] + VAR_c5[:,:,:,48] 

#%% Determine specific groups HOS-126
TSB_h1 = VAR_h1.sum(['species'])
TCB_h1 = TSB_h1 - VAR_h1[:,:,:,8] - VAR_h1[:,:,:,18] - VAR_h1[:,:,:,38] - VAR_h1[:,:,:,39]
B10_h1 = TCB_h1 - VAR_h1[:,:,:,34] - VAR_h1[:,:,:,31] - VAR_h1[:,:,:,1] - VAR_h1[:,:,:,25] - VAR_h1[:,:,:,27] - VAR_h1[:,:,:,28] - VAR_h1[:,:,:,29] - VAR_h1[:,:,:,30] - VAR_h1[:,:,:,23] - VAR_h1[:,:,:,22] - VAR_h1[:,:,:,14] - VAR_h1[:,:,:,48] - VAR_h1[:,:,:,49] - VAR_h1[:,:,:,32] - VAR_h1[:,:,:,33]
B30_h1 = B10_h1 - VAR_h1[:,:,:,24] - VAR_h1[:,:,:,24] - VAR_h1[:,:,:,4] - VAR_h1[:,:,:,7] - VAR_h1[:,:,:,11] - VAR_h1[:,:,:,13] - VAR_h1[:,:,:,17] - VAR_h1[:,:,:,37]
COM_h1 = VAR_h1[:,:,:,2] + VAR_h1[:,:,:,5] + VAR_h1[:,:,:,6] +VAR_h1[:,:,:,7] + VAR_h1[:,:,:,9] + VAR_h1[:,:,:,10] + VAR_h1[:,:,:,11] + VAR_h1[:,:,:,13] + VAR_h1[:,:,:,15] + VAR_h1[:,:,:,16] + VAR_h1[:,:,:,17] + VAR_h1[:,:,:,20] + VAR_h1[:,:,:,21] + VAR_h1[:,:,:,24] + VAR_h1[:,:,:,31] + VAR_h1[:,:,:,35] + VAR_h1[:,:,:,36] + VAR_h1[:,:,:,37]
COM_h1 = COM_h1 + VAR_h1[:,:,:,41] + VAR_h1[:,:,:,42] + VAR_h1[:,:,:,43] + VAR_h1[:,:,:,44] + VAR_h1[:,:,:,45] + VAR_h1[:,:,:,46] + VAR_h1[:,:,:,47] + VAR_h1[:,:,:,48] 

#%% Determine specific groups HOS-585
TSB_h5 = VAR_h5.sum(['species'])
TCB_h5 = TSB_h5 - VAR_h5[:,:,:,8] - VAR_h5[:,:,:,18] - VAR_h5[:,:,:,38] - VAR_h5[:,:,:,39]
B10_h5 = TCB_h5 - VAR_h5[:,:,:,34] - VAR_h5[:,:,:,31] - VAR_h5[:,:,:,1] - VAR_h5[:,:,:,25] - VAR_h5[:,:,:,27] - VAR_h5[:,:,:,28] - VAR_h5[:,:,:,29] - VAR_h5[:,:,:,30] - VAR_h5[:,:,:,23] - VAR_h5[:,:,:,22] - VAR_h5[:,:,:,14] - VAR_h5[:,:,:,48] - VAR_h5[:,:,:,49] - VAR_h5[:,:,:,32] - VAR_h5[:,:,:,33]
B30_h5 = B10_h5 - VAR_h5[:,:,:,24] - VAR_h5[:,:,:,24] - VAR_h5[:,:,:,4] - VAR_h5[:,:,:,7] - VAR_h5[:,:,:,11] - VAR_h5[:,:,:,13] - VAR_h5[:,:,:,17] - VAR_h5[:,:,:,37]
COM_h5 = VAR_h5[:,:,:,2] + VAR_h5[:,:,:,5] + VAR_h5[:,:,:,6] +VAR_h5[:,:,:,7] + VAR_h5[:,:,:,9] + VAR_h5[:,:,:,10] + VAR_h5[:,:,:,11] + VAR_h5[:,:,:,13] + VAR_h5[:,:,:,15] + VAR_h5[:,:,:,16] + VAR_h5[:,:,:,17] + VAR_h5[:,:,:,20] + VAR_h5[:,:,:,21] + VAR_h5[:,:,:,24] + VAR_h5[:,:,:,31] + VAR_h5[:,:,:,35] + VAR_h5[:,:,:,36] + VAR_h5[:,:,:,37]
COM_h5 = COM_h5 + VAR_h5[:,:,:,41] + VAR_h5[:,:,:,42] + VAR_h5[:,:,:,43] + VAR_h5[:,:,:,44] + VAR_h5[:,:,:,45] + VAR_h5[:,:,:,46] + VAR_h5[:,:,:,47] + VAR_h5[:,:,:,48] 

#%% Make datasets TSB, TCB, COM
TSB1 = TSB_c1.to_dataset(name='TSB')
TSB2 = TSB_h1.to_dataset(name='TSB')
TSB3 = TSB_c5.to_dataset(name='TSB')
TSB4 = TSB_h5.to_dataset(name='TSB')

TCB1 = TCB_c1.to_dataset(name='TCB')
TCB2 = TCB_h1.to_dataset(name='TCB')
TCB3 = TCB_c5.to_dataset(name='TCB')
TCB4 = TCB_h5.to_dataset(name='TCB')

COM1 = COM_c1.to_dataset(name='COM')
COM2 = COM_h1.to_dataset(name='COM')
COM3 = COM_c5.to_dataset(name='COM')
COM4 = COM_h5.to_dataset(name='COM')

#%% Directories where saved data goes to
data_snel1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/CTL_126/'
data_snel2='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/HOS_126/'
data_snel3='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/CTL_585/'
data_snel4='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/HOS_585/' 

#%% Save TSB, TCB, COM
TSB1.to_netcdf(data_snel1+'TSB_2015_2100.nc')
TSB2.to_netcdf(data_snel2+'TSB_2015_2100.nc')
TSB3.to_netcdf(data_snel3+'TSB_2015_2100.nc')
TSB4.to_netcdf(data_snel4+'TSB_2015_2100.nc')

TCB1.to_netcdf(data_snel1+'TCB_2015_2100.nc')
TCB2.to_netcdf(data_snel2+'TCB_2015_2100.nc')
TCB3.to_netcdf(data_snel3+'TCB_2015_2100.nc')
TCB4.to_netcdf(data_snel4+'TCB_2015_2100.nc')

COM1.to_netcdf(data_snel1+'COM_2015_2100.nc')
COM2.to_netcdf(data_snel2+'COM_2015_2100.nc')
COM3.to_netcdf(data_snel3+'COM_2015_2100.nc')
COM4.to_netcdf(data_snel4+'COM_2015_2100.nc')

#%% Select data for TSB, TCB and COM contour plots
TSB1_1 = TSB1.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
TSB1_2 = TSB1.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
TSB2_1 = TSB2.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
TSB2_2 = TSB2.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
TSB3_1 = TSB3.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
TSB3_2 = TSB3.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
TSB4_1 = TSB4.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
TSB4_2 = TSB4.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()

TCB1_1 = TCB1.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
TCB1_2 = TCB1.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
TCB2_1 = TCB2.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
TCB2_2 = TCB2.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
TCB3_1 = TCB3.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
TCB3_2 = TCB3.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
TCB4_1 = TCB4.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
TCB4_2 = TCB4.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()

COM1_1 = COM1.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
COM1_2 = COM1.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
COM2_1 = COM2.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
COM2_2 = COM2.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
COM3_1 = COM3.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
COM3_2 = COM3.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
COM4_1 = COM4.sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
COM4_2 = COM4.sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()

#%% Get data for microzooplankton (nsmz), mesozooplankton (nmdz) and large zooplankton (nlgz)
nsmz1_1 = VAR_c1[:,:,:,30].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nsmz1_2 = VAR_c1[:,:,:,30].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
nsmz2_1 = VAR_h1[:,:,:,30].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nsmz2_2 = VAR_h1[:,:,:,30].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
nsmz3_1 = VAR_c5[:,:,:,30].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nsmz3_2 = VAR_c5[:,:,:,30].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
nsmz4_1 = VAR_h5[:,:,:,30].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nsmz4_2 = VAR_h5[:,:,:,30].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()

nmdz1_1 = VAR_c1[:,:,:,29].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nmdz1_2 = VAR_c1[:,:,:,29].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
nmdz2_1 = VAR_h1[:,:,:,29].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nmdz2_2 = VAR_h1[:,:,:,29].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
nmdz3_1 = VAR_c5[:,:,:,29].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nmdz3_2 = VAR_c5[:,:,:,29].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
nmdz4_1 = VAR_h5[:,:,:,29].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nmdz4_2 = VAR_h5[:,:,:,29].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()

nlgz1_1 = VAR_c1[:,:,:,23].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nlgz1_2 = VAR_c1[:,:,:,23].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
nlgz2_1 = VAR_h1[:,:,:,23].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nlgz2_2 = VAR_h1[:,:,:,23].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
nlgz3_1 = VAR_c5[:,:,:,23].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nlgz3_2 = VAR_c5[:,:,:,23].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()
nlgz4_1 = VAR_h5[:,:,:,23].sel(time=slice('2016-01-01','2020-12-31')).mean('time').compute().squeeze()
nlgz4_2 = VAR_h5[:,:,:,23].sel(time=slice('2095-01-01','2099-12-31')).mean('time').compute().squeeze()

#%% Save data
TSB1_1.to_netcdf(data_snel1+'TSB_yr_control_2016_2020_126.nc')
TSB2_1.to_netcdf(data_snel2+'TSB_yr_hosing_05_2016_2020_126.nc')
TSB3_1.to_netcdf(data_snel3+'TSB_yr_control_2016_2020_585.nc')
TSB4_1.to_netcdf(data_snel4+'TSB_yr_hosing_05_2016_2020_585.nc')
TSB1_2.to_netcdf(data_snel1+'TSB_yr_control_2095_2099_126.nc')
TSB2_2.to_netcdf(data_snel2+'TSB_yr_hosing_05_2095_2099_126.nc')
TSB3_2.to_netcdf(data_snel3+'TSB_yr_control_2095_2099_585.nc')
TSB4_2.to_netcdf(data_snel4+'TSB_yr_hosing_05_2095_2099_585.nc')

TCB1_1.to_netcdf(data_snel1+'TCB_yr_control_2016_2020_126.nc')
TCB2_1.to_netcdf(data_snel2+'TCB_yr_hosing_05_2016_2020_126.nc')
TCB3_1.to_netcdf(data_snel3+'TCB_yr_control_2016_2020_585.nc')
TCB4_1.to_netcdf(data_snel4+'TCB_yr_hosing_05_2016_2020_585.nc')
TCB1_2.to_netcdf(data_snel1+'TCB_yr_control_2095_2099_126.nc')
TCB2_2.to_netcdf(data_snel2+'TCB_yr_hosing_05_2095_2099_126.nc')
TCB3_2.to_netcdf(data_snel3+'TCB_yr_control_2095_2099_585.nc')
TCB4_2.to_netcdf(data_snel4+'TCB_yr_hosing_05_2095_2099_585.nc')

COM1_1.to_netcdf(data_snel1+'COM_yr_control_2016_2020_126.nc')
COM2_1.to_netcdf(data_snel2+'COM_yr_hosing_05_2016_2020_126.nc')
COM3_1.to_netcdf(data_snel3+'COM_yr_control_2016_2020_585.nc')
COM4_1.to_netcdf(data_snel4+'COM_yr_hosing_05_2016_2020_585.nc')
COM1_2.to_netcdf(data_snel1+'COM_yr_control_2095_2099_126.nc')
COM2_2.to_netcdf(data_snel2+'COM_yr_hosing_05_2095_2099_126.nc')
COM3_2.to_netcdf(data_snel3+'COM_yr_control_2095_2099_585.nc')
COM4_2.to_netcdf(data_snel4+'COM_yr_hosing_05_2095_2099_585.nc')

#%% Rename variable in dataarray
nsmz1_1 = nsmz1_1.to_dataset(name='nsmz')
nsmz1_2 = nsmz1_2.to_dataset(name='nsmz')
nsmz2_1 = nsmz2_1.to_dataset(name='nsmz')
nsmz2_2 = nsmz2_2.to_dataset(name='nsmz')
nsmz3_1 = nsmz3_1.to_dataset(name='nsmz')
nsmz3_2 = nsmz3_2.to_dataset(name='nsmz')
nsmz4_1 = nsmz4_1.to_dataset(name='nsmz')
nsmz4_2 = nsmz4_2.to_dataset(name='nsmz')

nmdz1_1 = nmdz1_1.to_dataset(name='nmdz')
nmdz1_2 = nmdz1_2.to_dataset(name='nmdz')
nmdz2_1 = nmdz2_1.to_dataset(name='nmdz')
nmdz2_2 = nmdz2_2.to_dataset(name='nmdz')
nmdz3_1 = nmdz3_1.to_dataset(name='nmdz')
nmdz3_2 = nmdz3_2.to_dataset(name='nmdz')
nmdz4_1 = nmdz4_1.to_dataset(name='nmdz')
nmdz4_2 = nmdz4_2.to_dataset(name='nmdz')

nlgz1_1 = nlgz1_1.to_dataset(name='nlgz')
nlgz1_2 = nlgz1_2.to_dataset(name='nlgz')
nlgz2_1 = nlgz2_1.to_dataset(name='nlgz')
nlgz2_2 = nlgz2_2.to_dataset(name='nlgz')
nlgz3_1 = nlgz3_1.to_dataset(name='nlgz')
nlgz3_2 = nlgz3_2.to_dataset(name='nlgz')
nlgz4_1 = nlgz4_1.to_dataset(name='nlgz')
nlgz4_2 = nlgz4_2.to_dataset(name='nlgz')

#%% Save data
nsmz1_1.to_netcdf(data_snel1+'nsmz_yr_control_2016_2020_126.nc')
nsmz2_1.to_netcdf(data_snel2+'nsmz_yr_hosing_05_2016_2020_126.nc')
nsmz3_1.to_netcdf(data_snel3+'nsmz_yr_control_2016_2020_585.nc')
nsmz4_1.to_netcdf(data_snel4+'nsmz_yr_hosing_05_2016_2020_585.nc')
nsmz1_2.to_netcdf(data_snel1+'nsmz_yr_control_2095_2099_126.nc')
nsmz2_2.to_netcdf(data_snel2+'nsmz_yr_hosing_05_2095_2099_126.nc')
nsmz3_2.to_netcdf(data_snel3+'nsmz_yr_control_2095_2099_585.nc')
nsmz4_2.to_netcdf(data_snel4+'nsmz_yr_hosing_05_2095_2099_585.nc')

nmdz1_1.to_netcdf(data_snel1+'nmdz_yr_control_2016_2020_126.nc')
nmdz2_1.to_netcdf(data_snel2+'nmdz_yr_hosing_05_2016_2020_126.nc')
nmdz3_1.to_netcdf(data_snel3+'nmdz_yr_control_2016_2020_585.nc')
nmdz4_1.to_netcdf(data_snel4+'nmdz_yr_hosing_05_2016_2020_585.nc')
nmdz1_2.to_netcdf(data_snel1+'nmdz_yr_control_2095_2099_126.nc')
nmdz2_2.to_netcdf(data_snel2+'nmdz_yr_hosing_05_2095_2099_126.nc')
nmdz3_2.to_netcdf(data_snel3+'nmdz_yr_control_2095_2099_585.nc')
nmdz4_2.to_netcdf(data_snel4+'nmdz_yr_hosing_05_2095_2099_585.nc')

nlgz1_1.to_netcdf(data_snel1+'nlgz_yr_control_2016_2020_126.nc')
nlgz2_1.to_netcdf(data_snel2+'nlgz_yr_hosing_05_2016_2020_126.nc')
nlgz3_1.to_netcdf(data_snel3+'nlgz_yr_control_2016_2020_585.nc')
nlgz4_1.to_netcdf(data_snel4+'nlgz_yr_hosing_05_2016_2020_585.nc')
nlgz1_2.to_netcdf(data_snel1+'nlgz_yr_control_2095_2099_126.nc')
nlgz2_2.to_netcdf(data_snel2+'nlgz_yr_hosing_05_2095_2099_126.nc')
nlgz3_2.to_netcdf(data_snel3+'nlgz_yr_control_2095_2099_585.nc')
nlgz4_2.to_netcdf(data_snel4+'nlgz_yr_hosing_05_2095_2099_585.nc')


