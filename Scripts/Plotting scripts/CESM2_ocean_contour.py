#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    
'Global marine ecosystem response to a strong AMOC weakening under low and high future emission scenarios'
Amber A. Boot, J. Steenbeek, M. Coll, A.S. von der Heydt and H.A. Dijkstra
Submitted to Earth's Future March 2024

Script for contour plots of ocean fields from the CESM2 output:
    1. Sea surface temperature (Fig. 2)
    2. Nitrate concentrations (Fig. 3)
    3. Net Primary Production (Fig. 4)
    4. Total phytoplankton biomass (Fig. 5)
    5. Stratification (Fig. S4)
    6. Upwelling velocities (Fig. S5)
    7. Phosphate concentrations (Fig. S7)
    8. Silicate concentrations (Fig. S8)
    9. Dissolved iron concentrations (Fig. S9)
    10. Export production at 100m (Fig. S10)
    11. Diatom net primary production (Fig. S11)
    12. Small phytoplankton net primary production (Fig. S12)
    13. Diazotroph biomass (Fig. S13)
    14. Diatom biomass (Fig. S14)
    15. Small phytoplankton biomass (Fig. S15)
    16. Ocean temperature averaged top 150m (Fig. S16)
    17. Ocean temperatuee averaged column (Fig. S17)
    18. Ocean temperature bottom (Fig. S18)
    19. Barotropic stream function (Fig. S19)
    20. Maximum mixed layer depth (Fig. S20)
    
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

VARS= 'TEMP','NO3','photoC_TOT','phyto','Strat','WVEL','PO4','SiO3_150m','Fe_150m','POC_flux_100m','photoC_diat_zint','photoC_sp_zint','diazC','diatC','spC','TEMP_150m','TEMP_col','TEMP_bot','BSF','HMXL'                         # Name of dataset
UNIT = '[$^{\circ}$C]','[mol m$^{-2}$]','[mol C m$^{-2}$ s$^{-1}$]','[mol C m$^{-2}$]','[kg m$^{-3}$]','[m day$^{-1}$]','[mol m$^{-2}$]','[mol m$^{-2}$]','[mmol m$^{-2}$]','[mol C m$^{-2}$ s$^{-1}$]','[mol C m$^{-2}$ s$^{-1}$]','[mol C m$^{-2}$ s$^{-1}$]','[mol C m$^{-2}$]','[mol C m$^{-2}$]','[mol C m$^{-2}$]','[$^{\circ}$C]','[$^{\circ}$C]','[$^{\circ}$C]','[Sv]','[m]'     # Units
DESCR = 'SST','NO$_3^-$ (0-150m)','NPP','Total phytoplankton','Stratification','Vertical velocity (z = 150m)','PO$_4^{3-}$ (0-150m)','SiO$_3$ (0-150m)','Fe (0-150m)','EP (z = 100m)','NPP (diat; 0-150m)','NPP (sp; 0-150m)','Diazotrophs (0-150m)','Diatom (0-150m)','Small phytoplankton (0-150m)','Temperature (0-150m)','Temperature (column)','Temperature (bottom)','Barotropic $\psi$','MLD (max)'                 # Description used in plots
CMAP = ([cm.cm.thermal,cm.cm.algae,cm.cm.algae_r,cm.cm.tempo_r,cm.cm.ice,'PuOr_r',cm.cm.algae,cm.cm.algae,cm.cm.algae,cm.cm.matter_r,cm.cm.algae_r,cm.cm.algae_r,cm.cm.tempo_r,cm.cm.tempo_r,cm.cm.tempo_r,cm.cm.thermal,cm.cm.thermal,cm.cm.thermal,'PiYG',cm.cm.deep_r])                 # Used color maps
FIG_NR = '2','3','4','5','S4','S5','S7','S8','S9','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20'                                   # Figure number

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

    ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
    ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
    ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
    ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
    ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)
    
    ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
    ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
    ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
    ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

#%% Function to calculate density
def dens(salt,temp,lev):
    CT1=gsw.CT_from_pt(salt,temp)
    p1=gsw.p_from_z(lev,salt['lat'])
    
    rho=gsw.density.rho(salt,CT1,p1)

    return rho

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
n = 0 
for var_i in range(len(VARS)-n):
    # Loop to loop over the different variables
    
    var=VARS[var_i+n]     # Get variable
    name=VARS[var_i+n]    # Get name of dataset # Redundant now because var == name by default
    unit=UNIT[var_i+n]    # Get units for plot
    descr=DESCR[var_i+n]  # Get description of variable
    cmap=CMAP[var_i+n]    # Get color map
    fig_nr=FIG_NR[var_i+n]# Get figure number

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
    
    if var == 'Strat':
        # Stratificatin is defined as the density difference between z = 200m and the surface
        # To determine the density we need temperature and salinity 
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'TEMP_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARa1=load_var2['TEMP'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'TEMP_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARb1=load_var1['TEMP'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'SALT_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARa2=load_var2['SALT'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/SALT_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARb2=load_var1['SALT'][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'TEMP2_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARc1=load_var2['TEMP2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'TEMP2_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARd1=load_var1['TEMP2'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'SALT2_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARc2=load_var2['SALT2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/SALT2_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VARd2=load_var1['SALT2'][:,:].compute().squeeze()  
        
        VAR21 = dens(VARc2,VARc1,-200)-dens(VARa2,VARa1,0)
        VAR11 = dens(VARd2,VARd1,-200)-dens(VARb2,VARb1,0)
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'TEMP_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARa1=load_var2['TEMP'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'TEMP_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARb1=load_var1['TEMP'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'SALT_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARa2=load_var2['SALT'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/SALT_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARb2=load_var1['SALT'][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'TEMP2_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARc1=load_var2['TEMP2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'TEMP2_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARd1=load_var1['TEMP2'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'SALT2_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARc2=load_var2['SALT2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/SALT2_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VARd2=load_var1['SALT2'][:,:].compute().squeeze()  
        
        VAR41 = dens(VARc2,VARc1,-200)-dens(VARa2,VARa1,0)
        VAR31 = dens(VARd2,VARd1,-200)-dens(VARb2,VARb1,0)
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'TEMP_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARa1=load_var2['TEMP'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'TEMP_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARb1=load_var1['TEMP'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'SALT_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARa2=load_var2['SALT'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/SALT_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARb2=load_var1['SALT'][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'TEMP2_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARc1=load_var2['TEMP2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/'+'TEMP2_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARd1=load_var1['TEMP2'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+'SALT2_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARc2=load_var2['SALT2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel1}/SALT2_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VARd2=load_var1['SALT2'][:,:].compute().squeeze()  
        
        VAR22 = dens(VARc2,VARc1,-200)-dens(VARa2,VARa1,0)
        VAR12 = dens(VARd2,VARd1,-200)-dens(VARb2,VARb1,0)
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'TEMP_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARa1=load_var2['TEMP'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'TEMP_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARb1=load_var1['TEMP'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'SALT_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARa2=load_var2['SALT'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/SALT_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARb2=load_var1['SALT'][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'TEMP2_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARc1=load_var2['TEMP2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/'+'TEMP2_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARd1=load_var1['TEMP2'][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+'SALT2_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARc2=load_var2['SALT2'][:,:].compute().squeeze()
        #%%
        load_var1 = xr.open_dataset(f'{data_snel3}/SALT2_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VARd2=load_var1['SALT2'][:,:].compute().squeeze()  
        
        VAR42 = dens(VARc2,VARc1,-200)-dens(VARa2,VARa1,0)
        VAR32 = dens(VARd2,VARd1,-200)-dens(VARb2,VARb1,0)
    
    elif var == 'phyto':
        # Total phytoplankton biomass is sum of diazotrophs, diatoms and small phytoplankton
        # Load in diazotrophs
        var = 'diazC'
        name = 'diazC'
            
        load_var2 = xr.open_dataset(f'{data_snel2}/'+name+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VAR211=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel1}/'+name+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VAR111=load_var1[var].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+name+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VAR411=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+name+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VAR311=load_var1[var].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+name+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VAR221=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel1}/'+name+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VAR121=load_var1[var].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+name+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VAR421=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+name+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VAR321=load_var1[var].compute().squeeze()
        
        # Load in diatoms
        var = 'diatC'
        name = 'diatC'
            
        load_var2 = xr.open_dataset(f'{data_snel2}/'+name+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VAR212=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel1}/'+name+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VAR112=load_var1[var].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+name+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VAR412=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+name+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VAR312=load_var1[var].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+name+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VAR222=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel1}/'+name+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VAR122=load_var1[var].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+name+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VAR422=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+name+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VAR322=load_var1[var].compute().squeeze()
        
        # Load in small phytoplankton
        var = 'spC'
        name = 'spC'
            
        load_var2 = xr.open_dataset(f'{data_snel2}/'+name+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VAR213=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel1}/'+name+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VAR113=load_var1[var].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+name+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VAR413=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+name+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VAR313=load_var1[var].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+name+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VAR223=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel1}/'+name+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VAR123=load_var1[var].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel4}/'+name+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VAR423=load_var2[var].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+name+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VAR323=load_var1[var].compute().squeeze()
            
        # Sum up the three phytoplankton types
        VAR11 = VAR111 + VAR112 + VAR113
        VAR21 = VAR211 + VAR212 + VAR213
        VAR31 = VAR311 + VAR312 + VAR313
        VAR41 = VAR411 + VAR412 + VAR413
        
        VAR12 = VAR121 + VAR122 + VAR123
        VAR22 = VAR221 + VAR222 + VAR223
        VAR32 = VAR321 + VAR322 + VAR323
        VAR42 = VAR421 + VAR422 + VAR423

    else:
        load_var2 = xr.open_dataset(f'{data_snel2}/'+name+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VAR21=load_var2[var][:,:].compute().squeeze()
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+name+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_126.nc' )
        VAR11=load_var1[var][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+name+'_yr_'+RUN2+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VAR41=load_var2[var][:,:].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+name+'_yr_'+RUN1+'_'+str(year_start1)+'_'+str(year_end1)+'_585.nc' )
        VAR31=load_var1[var][:,:].compute().squeeze()
        
        load_var2 = xr.open_dataset(f'{data_snel2}/'+name+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VAR22=load_var2[var][:,:].compute().squeeze()
    
        load_var1 = xr.open_dataset(f'{data_snel1}/'+name+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_126.nc' )
        VAR12=load_var1[var][:,:].compute().squeeze()
    
        load_var2 = xr.open_dataset(f'{data_snel4}/'+name+'_yr_'+RUN2+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VAR42=load_var2[var][:,:].compute().squeeze()
        
        load_var1 = xr.open_dataset(f'{data_snel3}/'+name+'_yr_'+RUN1+'_'+str(year_start2)+'_'+str(year_end2)+'_585.nc' )
        VAR32=load_var1[var][:,:].compute().squeeze()

    #%% Manipulate data if necessary
    if name == 'phyto' or name == 'spC_150m' or name == 'diatC_150m' or name == 'diazC_150m' or name == 'TEMP_150m' or name == 'TEMP_bot' or name == 'TEMP_col':
        # These datasets are already regridded
        VAR11_gr = VAR11
        VAR21_gr = VAR21 
        VAR31_gr = VAR31
        VAR41_gr = VAR41
        
        VAR12_gr = VAR12
        VAR22_gr = VAR22
        VAR32_gr = VAR32
        VAR42_gr = VAR42
        
    else:
        VAR11_gr = regridder(VAR11)
        VAR11_gr=np.roll(VAR11_gr,-180)   
        VAR21_gr = regridder(VAR21)
        VAR21_gr=np.roll(VAR21_gr,-180) 
        VAR31_gr = regridder(VAR31)
        VAR31_gr=np.roll(VAR31_gr,-180)   
        VAR41_gr = regridder(VAR41)
        VAR41_gr=np.roll(VAR41_gr,-180) 
         
        VAR12_gr = regridder(VAR12)
        VAR12_gr=np.roll(VAR12_gr,-180)   
        VAR22_gr = regridder(VAR22)
        VAR22_gr=np.roll(VAR22_gr,-180) 
        VAR32_gr = regridder(VAR32)
        VAR32_gr=np.roll(VAR32_gr,-180)   
        VAR42_gr = regridder(VAR42)
        VAR42_gr=np.roll(VAR42_gr,-180)
    
    if var == 'WVEL':
        # Transform to m/day
        VAR11_gr = VAR11_gr * 1e-2 * 86400
        VAR21_gr = VAR21_gr * 1e-2 * 86400
        VAR31_gr = VAR31_gr * 1e-2 * 86400
        VAR41_gr = VAR41_gr * 1e-2 * 86400
        
        VAR12_gr = VAR12_gr * 1e-2 * 86400
        VAR22_gr = VAR22_gr * 1e-2 * 86400
        VAR32_gr = VAR32_gr * 1e-2 * 86400
        VAR42_gr = VAR42_gr * 1e-2 * 86400
        
    lat = area_gr.lat
    lon = area_gr.lon
        
    #%%
    # Directory for figures
    datadir = '/Users/daan/Desktop/EF2024_CESM2_EcoOcean/Figures/'
    
    # Select bounds for color scheme
    # Average 2016 - 2020
    Vmn0 = ([0,0,0,0,0,-0.45,0,0,0,0,0,0,0,0,0,0,0,0,-30,0])
    Vmx0 = ([35,2.5,6e-7,0.15,6,0.45,0.25,2.5,0.21,1e-7,4.5e-7,4.5e-7,0.007,0.15,0.07,25,15,12,30,200])
    
    # Difference CTL simulation (average 2095 - 2099 minus average 2016 - 2020)
    Vmn1 = ([-5,-0.6,-2.1e-7,-0.07,-5,-0.21,-0.078,-1.1,-0.05,-5e-8,-3.1e-7,-2.1e-7,-0.005,-0.11,-0.07,-5,-5,-2,-16,-80])
    
    # Difference HOS - CTL simulation (averaged 2095 - 2099)
    Vmn2 = ([-3.1,-0.6,-1.6e-7,-0.041,-5,-0.21,-0.078,-1.1,-0.05,-5e-8,-3.1e-7,-2.1e-7,-0.005,-0.11,-0.07,-5,-5,-2,-16,-80])
    
    vmn0=Vmn0[var_i+n]
    vmx0=Vmx0[var_i+n]
    
    vmn1=Vmn1[var_i+n]
    vmx1=-vmn1
    
    vmn2=Vmn2[var_i+n]
    vmx2=-vmn2
    
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
        
