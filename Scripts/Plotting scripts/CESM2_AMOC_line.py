#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script for line plots of CESM2 output:
    1. AMOC strength (Fig. 1c, f)
    
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
import matplotlib.colors as colors
    
#%%
year_start=2015
year_end1=2100
year_end2=2100
RUN1='control'
RUN2='hosing_05' # List of runs: control, hosing_05, 
quality=300 # quality when figure is saved in dpi
a=5         # amount of years averaged over
tm=a
save_fig='yes'
years=np.arange(year_start,year_end2+1,1)
dyr=year_end2-year_start
A=([2025,2040,2055,2070,2085,2100])

FS=18
LW=6

VARS='MOC'    
unit='[Sv]'
CONV= 1
descr='AMOC'

#%% Location CESM2 data per scenario and run 
data_snel1='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_126/'
data_snel2='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_126/'
data_snel3='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/CTL_585/'
data_snel4='/Users/daan/Desktop/EF2024_CESM2_EcoOcean/CESM2 data/HOS_585/'  

#%% Call on datasets (might need to change name of dataset)
load_var1 = xr.open_dataset(f'{data_snel1}/'+VARS+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end1)+'_0126.nc' )
VAR1_gr=load_var1[VARS][:,1,0,:,:].squeeze()
l1=load_var1[VARS][:,1,0,:,274].squeeze()

depth=-load_var1['moc_z']/100
lat=load_var1['lat_aux_grid']

load_var1 = xr.open_dataset(f'{data_snel2}/'+VARS+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end1)+'_0126.nc' )
VAR2_gr=load_var1[VARS][:,1,0,:,:].squeeze()
l2=load_var1[VARS][:,1,0,:,274].squeeze()

plt.plot(VAR2_gr['time'])

#%%
load_var1 = xr.open_dataset(f'{data_snel3}/'+VARS+'_yr_'+RUN1+'_'+str(year_start)+'_'+str(year_end2)+'_0585.nc' )
VAR3_gr=load_var1[VARS][:,1,0,:,:].squeeze()
l3=load_var1[VARS][:,1,0,:,274].squeeze()

depth=-load_var1['moc_z']/100
lat=load_var1['lat_aux_grid']

load_var1 = xr.open_dataset(f'{data_snel4}/'+VARS+'_yr_'+RUN2+'_'+str(year_start)+'_'+str(year_end2)+'_0585.nc' )
VAR4_gr=load_var1[VARS][:,1,0,:,:].squeeze()
l4=load_var1[VARS][:,1,0,:,274].squeeze()

#%% Determine maximum
len_1 = len(VAR2_gr['time'])

VAR1=np.zeros((len_1,))
VAR2=np.zeros((len_1,))

for i in range(len_1):
    VAR1[i]=max(l1[i,:])
    VAR2[i]=max(l2[i,:])

#%% Determine maximum
len_3 = len(VAR4_gr['time'])

VAR3=np.zeros((len_3,))
VAR4=np.zeros((len_3,))

for i in range(len_3):
    VAR3[i]=max(l3[i,:])
    VAR4[i]=max(l4[i,:])
    
#%% Turn into xarray and take moving average with period tm
VAR1=VAR1*xr.ones_like(l1[:,0].squeeze())
VAR_gr1=VAR1.rolling(time=tm, center=False).mean()

VAR2=VAR2*xr.ones_like(l1[:,0].squeeze())
VAR_gr2=VAR2.rolling(time=tm, center=False).mean()

VAR3=VAR3*xr.ones_like(l3[:,0].squeeze())
VAR_gr3=VAR3.rolling(time=tm, center=False).mean()

VAR4=VAR4*xr.ones_like(l3[:,0].squeeze())
VAR_gr4=VAR4.rolling(time=tm, center=False).mean()

#%% Figure  
datadir = '/Users/daan/Desktop/EF2024_CESM2_EcoOcean/Figures/'
t2=np.arange(2016,2100.5,1) 
t1=t2 
LW=4
FS=18

fig = plt.figure(figsize=(7, 5))
ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t2,VAR_gr3,linewidth=LW,color='tab:blue',linestyle='solid',label='CTL-585')
ax.plot(t2,VAR_gr4,linewidth=LW,color='tab:orange',linestyle='solid',label='HOS-585')
ax.plot(t1,VAR_gr1,linewidth=LW,color='tab:blue',linestyle='dashed',label='CTL-126')
ax.plot(t1,VAR_gr2,linewidth=LW,color='tab:orange',linestyle='dashed',label='HOS-126')
plt.xlabel('Time',fontsize=FS)
plt.ylabel('AMOC (26.5$^{\circ}$N) [Sv]',fontsize=FS)
plt.grid()
plt.xlim([2020,year_end2])
plt.xticks(A,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-4)

if save_fig == 'yes':
    plt.savefig(datadir+'Figure_1c.png', format='png', dpi=quality, transparent=True, bbox_inches = 'tight')

#%% Figure
fig = plt.figure(figsize=(7, 5))
ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t2,VAR_gr4-VAR_gr3,linewidth=LW,color='tab:green',label='585')
ax.plot(t1,VAR_gr2-VAR_gr1,linewidth=LW,color='tab:green',linestyle='dashed',label='126')
plt.xlabel('Time',fontsize=FS)
plt.ylabel('AMOC (26.5$^{\circ}$N) [Sv]',fontsize=FS)
plt.grid()
plt.xlim([2020,year_end2])
plt.ylim([-9,0.5])
plt.xticks(A,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-4)

if save_fig == 'yes':
    plt.savefig(datadir+'Figure_1f.png', format='png', dpi=quality, transparent=True, bbox_inches = 'tight')
