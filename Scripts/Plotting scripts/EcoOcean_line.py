#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script for line plots of EcoOcean output:
    1. Total System Biomass (Fig. 6a, d)
    2. Total Consumer Biomass (Fig. 6b, e)
    3. Commercial biomass (Fig. 6c, f)
    
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

#%% Name data directories for loading data and saving figures
datadir = '/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/'
data_snel1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/CTL_126/'
data_snel2='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/HOS_126/'
data_snel3='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/CTL_585/'
data_snel4='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/EcoOcean data/HOS_585/' 

figdir = '/Users/daan/Desktop/EF2024_CESM2_EcoOcean/Figures/' 

#%% Area of EcoOcean grid
area_1 = np.loadtxt(f'{datadir}/'+"cellarea.asc", skiprows=6)
area = xr.DataArray(area_1)  

lat = np.linspace(-84,90,175)
lon =(np.arange(0,359.5,1))

#%% Plotting variables
LW = 4      # Linewidth
FS = 20     # Font size
quality = 300   # Quality figures in dpi

#%% Variables
varN = (['TSB','TCB','COM'])
FIGNR1 = 'a','b','c'
FIGNR2 = 'd','e','f'

#%% Function to integrate data
def prepare_data(VAR):
    var_n = VAR
    var = var_n.where(var_n != var_n[0,-1,-1])
    var_area = var*area
    var_tot = np.sum(np.sum(var_area,axis=1),axis=1)
    var_tot_rel = (var_tot-var_tot[0])/var_tot[0]*100
    
    return var_area, var_tot_rel

#%% Function to produce line plot
def line_plot(data1,data2,data3,data4,var,type_i):
    
    # type_i == 0: The four scenarios
    # type_i == 1: Difference HOS - CTL
    
    if type_i == 0:
        plt.plot(np.arange(2015,2100,1),data3,linewidth=LW,color='tab:blue',label='CTL-585')
        plt.plot(np.arange(2015,2100,1),data4,linewidth=LW,color='tab:orange',label='HOS-585')
        plt.plot(np.arange(2015,2100,1),data1,linewidth=LW,linestyle='dashed',color='tab:blue',label='CTL-126')
        plt.plot(np.arange(2015,2100,1),data2,linewidth=LW,linestyle='dashed',color='tab:orange',label='HOS-126')
        plt.ylim([-16,2])
    elif type_i == 1:
        plt.plot(np.arange(2015,2100,1),data4-data3,linewidth=LW,color='tab:green',label='SSP5-8.5')
        plt.plot(np.arange(2015,2100,1),data2-data1,linewidth=LW,linestyle='dashed',color='tab:green',label='SSP1-2.6')
        plt.ylim([-5,2])
        
    plt.xlabel('Time',fontsize=FS-2)
    plt.ylabel('Relative change [%]',fontsize=FS-2)
    plt.xticks([2025,2040,2055,2070,2085],fontsize=FS-4)
    plt.yticks(fontsize=FS-4)
    
    if var == 'COM':
        plt.title('Commercial',fontsize=FS)
    else:
        plt.title(str(var),fontsize=FS)
    plt.grid()
    
    plt.legend(fontsize=FS-6)
    plt.xlim([2020,2099])

#%%
for var_i in range(len(varN)):

    filename = varN[var_i] + '_2015_2100.nc'
    var_n = varN[var_i]
    fignr1 = FIGNR1[var_i]
    fignr2 = FIGNR2[var_i]
    
    load_var1 = xr.open_dataset(f'{data_snel1}/'+filename)
    VAR_c1=load_var1[var_n].compute().squeeze() 
    load_var1 = xr.open_dataset(f'{data_snel3}/'+filename)
    VAR_c5=load_var1[var_n].compute().squeeze() 
    load_var1 = xr.open_dataset(f'{data_snel2}/'+filename)
    VAR_h1=load_var1[var_n].compute().squeeze() 
    load_var1 = xr.open_dataset(f'{data_snel4}/'+filename)
    VAR_h5=load_var1[var_n].compute().squeeze() 
        
    #%%
    VAR_c1 = VAR_c1.transpose("time","dim_0","dim_1")
    VAR_c5 = VAR_c5.transpose("time","dim_0","dim_1")
    VAR_h1 = VAR_h1.transpose("time","dim_0","dim_1")
    VAR_h5 = VAR_h5.transpose("time","dim_0","dim_1")
    
#%%
    cont_c1, line_c1 = prepare_data(VAR_c1)
    cont_c5, line_c5 = prepare_data(VAR_c5)
    cont_h1, line_h1 = prepare_data(VAR_h1)
    cont_h5, line_h5 = prepare_data(VAR_h5)
    
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_axes([0.25,0.15,0.67,0.8])
    line_plot(line_c1,line_h1,line_c5,line_h5,var_n,0)
    plt.savefig(figdir+'Figure_6'+fignr1+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_axes([0.25,0.15,0.67,0.8])
    line_plot(line_c1,line_h1,line_c5,line_h5,var_n,1)
    plt.savefig(figdir+'Figure_6'+fignr2+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')