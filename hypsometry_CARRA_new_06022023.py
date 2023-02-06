#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 12:04:02 2021
@author: AD and JEB
"""
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
from glob import glob
import pickle
# from netCDF4 import Dataset
import pandas as pd
from datetime import datetime 
from scipy.spatial import cKDTree
import matplotlib as mpl 
import metpy.calc
from metpy.units import units
import xarray as xr
import gzip
AD=1
# if AD:
    # os.environ['PROJ_LIB'] = r'C:/Users/Armin/Anaconda3/pkgs/proj4-5.2.0-ha925a31_1/Library/share' #Armin needed to not get an error with import basemap: see https://stackoverflow.com/questions/52295117/basemap-import-error-in-pycharm-keyerror-proj-lib
# from mpl_toolkits.basemap import Basemap

# CARRA grid info
# Lambert_Conformal()
#     grid_mapping_name: lambert_conformal_conic
#     standard_parallel: 72.0
#     longitude_of_central_meridian: -36.0
#     latitude_of_projection_origin: 72.0
#     earth_radius: 6367470.0
#     false_easting: 1334211.3405653758
#     false_northing: 1584010.8994621644
#     longitudeOfFirstGridPointInDegrees: 302.903
#     latitudeOfFirstGridPointInDegrees: 55.81


path='/Users/jason/Dropbox/CARRA/CARRA_rain/'
if AD:
    path='/home/rmean/Dokumente/Work/GEUS/CARRA_rain/'
    raw_path='/home/rmean/Dokumente/Work/GEUS/CARRA/'
    tool_path='/home/rmean/Dokumente/Work/GEUS/CARRA_tools/'
os.chdir(path)

#global plot settings
th=1
fs=10 # font size
ms=15# marker size
mult=1
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams["font.size"] = fs
plt.rcParams["mathtext.default"]='regular' #no italics on titles etc.

#define Box
LLlat=58
LLlon=-75+360      
URlat=85
URlon=-6+360

#define upscaling RACMO 1km
# def upscaleRACMOtoCARRA(data):
#     threshold=np.sqrt(2) * 2500
#     fn=tool_path+'resampling_key_RACMO1km_to_CARRA_dist.npy'
#     dist1=np.load(fn)
#     fn=tool_path+'resampling_key_RACMO1km_to_CARRA_idx.npy'
#     idx1=np.load(fn)
#     raw_values = data.flatten()[idx1]
#     raw_values[dist1 > threshold] = np.nan
#     upsampled_values = np.nanmean(raw_values, axis=1)
    
#     fn=tool_path+'resampling_key_RACMO1km_to_CARRA_idx1.npy'
#     idx2=np.load(fn)

#     tp_up=upsampled_values[idx2]
#     tpR=tp_up.reshape(niC,njC)
#     return tpR


fig = plt.figure(figsize=(10,5))
# ax1 = fig.add_subplot(111)
# n_cols=2 ; n_rows=1 ; 
# fig, (ax1, ax3) = plt.subplots(n_rows,n_cols,figsize=(10,4))
gs = fig.add_gridspec(1, 2, hspace=0, wspace=0)
(ax1, ax3) = gs.subplots(sharex='col', sharey='row')
ax2 = ax1.twiny()  #second x-axis
ax4 = ax3.twiny()  #second x-axis
mark=['-', '--', ':']
color=['g', 'brown', 'b', 'm', 'orange','purple']
label=['CARRA', 'ERA5', 'MAR 15km', 'NHM-SMAP', 'RACMO 5.5km', 'JRA-55']
xx0=0.1 ; yy0=1.5 ; dy2=-0.028 ; cc=0

#define resampling method
resampling='upscale'
# resampling='interpolation'


for jj in range(0,len(label)-1):  #len(label)#CARRA and ERA5

    if jj==0:  #CARRA
    
        ni=1269 ; nj=1069
        niC=1269 ; njC=1069
        fn='./ancil/2.5km_CARRA_west_lat_1269x1069.npy'
        lat=np.fromfile(fn, dtype=np.float32)
        lat_mat=lat.reshape(ni, nj)
        
        fn='./ancil/2.5km_CARRA_west_lon_1269x1069.npy'
        lon=np.fromfile(fn, dtype=np.float32)
        lon_mat=lon.reshape(ni, nj)
           
       
        # read ice mask
        fn='./ancil/CARRA_W_domain_ice_mask.nc'
        nc2 = xr.open_dataset(fn)
        
        mask = np.array(nc2.z)
        
        fn='./ancil/mask_peri_glaciers_ice_caps_1269x1069.npy'
        mask_peri = np.load(fn)
        mask_peri=np.rot90(mask_peri.T)
        # mask_fuzzy_PG=np.array(mask)
        # mask_fuzzy_PG[mask_peri==0]=np.nan
       
        LLlat=lat_mat[0,0]
        LLlon=lon_mat[0,0]-360
        lon0=lon_mat[int(round(ni/2)),int(round(nj/2))]-360
        lat0=lat_mat[int(round(ni/2)),int(round(nj/2))]         
        URlat=lat_mat[ni-1,nj-1]
        URlon=lon_mat[ni-1,nj-1]
        
        # m = Basemap(llcrnrlon=LLlon, llcrnrlat=LLlat, urcrnrlon=URlon, urcrnrlat=URlat, lat_0=lat0, lon_0=lon0, resolution='h', projection='lcc')

        mask_iceland=1
        mask_svalbard=1
        mask_jan_mayen=1
               
        if mask_jan_mayen:
            mask[((lon_mat-360>-15)&(lat_mat>66.6)&(lat_mat<75))]=0
        if mask_iceland:
            mask[((lon_mat-360>-30)&(lat_mat<66.6))]=0
        if mask_svalbard:
            mask[0:300,800:]=0
        
        mask[mask==0]=np.nan
        maskC=mask.copy()
        
        
        #---------------------------------------------------------- read in data
        # fn='./ancil/CARRA_W_domain_elev.nc'
        # ds=xr.open_dataset(fn)
        # elevC=ds.z
        
        fn='./ancil/CARRA_W_elev_lat_lon.nc'
        dsC=xr.open_dataset(fn)
        elevC=dsC.z.values  
        # elevC[mask_peri==0]=0
        
               
        years=np.arange(1998,2021).astype('str')
        n_years=len(years)
        wo=1
        sumx=np.zeros((ni, nj))
        varnams=['rf','tp','sf']
        cc=0
        for yy,year in enumerate(years):
            if yy>=0:
                fn='./output_annual/tp_'+year+'_'+str(ni)+'x'+str(nj)+'_float32.npy'
                print(year)
                tp=np.fromfile(fn, dtype=np.float32)
                tp=tp.reshape(ni, nj)
                sumx+=tp
                cc+=1
        
        tpC=sumx/cc 
        # tpC[mask_peri==0]=0

    
    if jj==1:  #ERA5
        #load Greenland land mask created with a watershed algorithm (to get rid of other land pixels)
        ni=1440 ; nj=721

        fn='./ancil/ERA5_regional_masks_raster_1440x721.npy'
        mask = np.fromfile(fn, dtype='int16')
        mask=mask.reshape(nj, ni)
        mask[((mask>1))]=0 #set 0 to non-greenland values, Greenland has a mask value of 1
        
        #load ERA5 to CARRA grid resampling key
        fn=tool_path+'resampling_key_ERA5_to_CARRA.pkl'
        infile = open(fn,'rb')
        df_res=pickle.load(infile)

        #ERA 5 elevation data
        fn='./ancil/ERA5_mask_and_terrain.nc'
        ds=xr.open_dataset(fn)
        dsz=ds.z[0,:,:].values *mask #read and mask data
        #slice dataset to bbox
        # z_sl=dsz.sel(longitude=slice(LLlon, URlon), latitude=slice(URlat, LLlat))
 
        if resampling=='upscale':   #ERA elev data with CARRA resolution
            print(resampling)
            elev_res= dsz[df_res.col_e, df_res.row_e]
            elev_res=elev_res.reshape(niC, njC)
        
        if resampling=='interpolation':   #interpolate into CARRA grid
            print(resampling)
            ds.load()   
            elev_res1 = ds.interp(latitude=dsC.latitude, longitude=dsC.longitude)
            elev_res=elev_res1.z[0,:,:].values  #read and mask data
          
        #from geopot to elevation
        z_elev2=elev_res* units['meter ** 2 / second ** 2']
        elevE=np.array(metpy.calc.geopotential_to_height(z_elev2)) #convert geopot into elevation
        # elevE[mask_peri==0]=0
        
        
        # ERA5 tp data
        fn=raw_path+'./ERA5_tp_1998_2020_monthly.grib'
        ds=xr.open_dataset(fn,engine='cfgrib')
        # ds = xr.load_dataset(fn, engine="cfgrib")
        dstp=ds.tp[:,:,:]*mask  
        
        if resampling=='interpolation':  #interpolate into CARRA grid
            print(resampling)
            ds.load() 
            ds1= ds.interp(latitude=dsC.latitude, longitude=dsC.longitude)
            dstp=ds1.tp
             
        #units conversion
        average_days_per_month = (7*31+4*30+28+(5/23))/12
        tp_E=np.sum(dstp.values, axis=0)*1000 *average_days_per_month/23  #tp mm per year #from m/month to mm/month of Greenland *31 since monthly data is in m/day
        # tpE[mask_peri==0]=0
        
        if resampling=='upscale':    #ERA tp data with CARRA resolution
            print(resampling)
            tpE= tp_E[df_res.col_e, df_res.row_e]
            tpE=tpE.reshape(niC, njC)
        

    if jj==2:  #MAR
        MAR_BS15=1 #choose with BS (1) or without (0) -> only with calculated for newest version 3.13
        # ni=446; nj=240 #6km
        # ni=179; nj=96 #15km old version
        ni=181; nj=101  #new version 
        
        # load MAR to CARRA grid resampling key
        # fn=tool_path+'resampling_key_MAR6km_to_CARRA.pkl' #6km
        # fn=tool_path+'resampling_key_MAR15km_to_CARRA.pkl' #15km old version
        fn=tool_path+'resampling_key_MAR15km_BS_to_CARRA.pkl' #15km new version
        infile = open(fn,'rb')
        df_res=pickle.load(infile)
        
        #MAR data
        # fn=raw_path+'MARv3.11.5-6km.nc' #6km
        # fn=raw_path+'MARv3.11.5-15km.nc' #15km  old version
        # fn=raw_path+'withBS/MARv3.13.0-withBS-15km-yearly-ERA5-1998-2020.nc' #15km new version (version 3.13.0 on 6.2.23)
        if MAR_BS15: fn=raw_path+'withBS/MARv3.13.0-withBS-15km-yearly-ERA5-1998-2020.nc'  #15km new version (version 3.13.0 on 6.2.23)
        ds=xr.open_dataset(fn)
        
        #MAR elevation
        elev=np.array(ds.z[:, :]) #surface height in m
        # elev=np.array(ds.SH[:, :]) #surface height in m -> old version
        
        #MAR tp 1998-2020
        tp=np.array(ds.tp[:, :]) #tp in mmWE/yr
        tp_M=np.mean(tp, axis=0)

        #old version  (original dataset is 1950-2020)
        # sf=np.array(ds.SF2[48:, :, :]) #snowfall in mmWE/yr
        # rf=np.array(ds.RF2[48:, :, :]) #rainfall in mmWE/yr
        # tp=sf+rf
        
        if resampling=='upscale':    #MAR tp data with CARRA resolution
            print(resampling)
            tpM= tp_M[df_res.row_m, df_res.col_m]
            tpM=tpM.reshape(niC, njC)
        
            elevM=elev[df_res.row_m, df_res.col_m]
            elevM=elevM.reshape(niC, njC)
      
        
    if jj==3:  #NHM-SMAP
        ni=550; nj=450
        
        #load NHM to CARRA grid resampling key
        fn=tool_path+'resampling_key_NHM_to_CARRA.pkl'
        infile = open(fn,'rb')
        df_res=pickle.load(infile)
        
        #NHM data
        fn=raw_path+'NHM-SMAP_v2_9820.nc'
        ds=xr.open_dataset(fn)
        
        #NHM elevation
        elev=np.array(ds.HEIGHT[0,:,:]) #surface height in m
        
        #NHM average tp 1998-2020 
        tp_N=np.array(ds.PRECIP[0,:,:]) #rainfall in mmWE

        
        if resampling=='upscale':    #NHM tp data with CARRA resolution
            print(resampling)
            tpN= tp_N[df_res.row_n, df_res.col_n]
            tpN=tpN.reshape(niC, njC)
        
            elevN=elev[df_res.row_n, df_res.col_n]
            elevN=elevN.reshape(niC, njC)
            
    if jj==4:  #RACMO
        ni=2700; nj=1496
        
        #load RACMO 5km to CARRA grid resampling key
        fn=tool_path+'resampling_key_RACMO5km_to_CARRA.pkl'
        infile = open(fn,'rb')
        df_res=pickle.load(infile)
        
        #RACMO data
        fn=raw_path+'precip.1958-2020.FGRN055_BN_RACMO2.3p2_ERA5_3h_FGRN055.YY.nc.gz'
        with gzip.open(fn) as gz:
            ds=xr.open_dataset(gz,decode_times=False)
            #RACMO average tp 1998-2020 
            tp_data=np.array(ds.precip[:,0,:,:]) #total precip in mmWE
            tp_R=np.mean(tp_data[40:,:,:], axis=0)
        
        #RACMO elevation
        fn=raw_path+'FGRN055_Masks_5.5km.nc.gz'
        with gzip.open(fn) as gz:
            ds=xr.open_dataset(gz,decode_times=False)
            elev=np.array(ds.Topography[:,:]) #surface height in m
        
        
        if resampling=='upscale':    #RACMO tp data with CARRA resolution
            print(resampling)
            
            #direct match resampling for RACMO 5km
            tpR= tp_R[df_res.row_r, df_res.col_r]
            tpR=tpR.reshape(niC, njC)
        
            elevR=elev[df_res.row_r, df_res.col_r]
            elevR=elevR.reshape(niC, njC)
            
            #upscale method for RACMO 1km (needed since RACMO resolution (1km) < CARRA resolution (2.5km))
            # tpR=upscaleRACMOtoCARRA(tp_R)
            # elevR=upscaleRACMOtoCARRA(elev)
      
    if jj==5:  #JRA
        ni=73; nj=145
        
        #load JRA to CARRA grid resampling key
        fn=tool_path+'resampling_key_JRA_to_CARRA.pkl'
        infile = open(fn,'rb')
        df_res=pickle.load(infile)
        
        #JRA data
        fn=raw_path+'JRA-55/jra55prec9820.nc'
        ds=xr.open_dataset(fn)
        
        #JRA average tp 1998-2020 
        tp_J=np.array(ds.prec[:,:]) #rainfall in mmWE
        
        #JRA elevation
        fn=raw_path+'JRA-55/jra55elev9820.nc'
        ds=xr.open_dataset(fn)
        elev=np.array(ds.elev[:,:]) #surface height in m
        elev[elev<0]=0
        # plt.imshow(elev)
               
        if resampling=='upscale':    #JRA tp data with CARRA resolution
            print(resampling)
            tpJ= tp_J[df_res.col_j, df_res.row_j]
            tpJ=tpJ.reshape(niC, njC)
        
            elevJ=elev[df_res.col_j, df_res.row_j]
            elevJ=elevJ.reshape(niC, njC)

    #-----------------------------------------------calculate area and hypso
    
    #### tp was masked for consistency!! (meaning that fuzzy values decrease in tp)
    if jj==0: meanx=tpC.copy()*maskC; elev=elevC.copy() #don't mask elevation with *maskC since that is decreasing elevation values (due to range between 0 and 1)
    if jj==1: meanx=tpE.copy()*maskC; elev=elevE.copy()
    if jj==2: meanx=tpM.copy()*maskC; elev=elevM.copy()
    if jj==3: meanx=tpN.copy()*maskC; elev=elevN.copy()
    if jj==4: meanx=tpR.copy()*maskC; elev=elevR.copy()
    if jj==5: meanx=tpJ.copy()*maskC; elev=elevJ.copy()
    binx=100
    elevs=np.arange(0,3250,binx)
    n=len(elevs)
    
    areax=2.5e3**2
    
    hypso=np.zeros(n)*np.nan
    areas=np.zeros(n)*np.nan
    tp_cum=np.zeros(n)*np.nan
    
    for i in range(n):
        v=[((elev>elevs[i])&(elev<=elevs[i]+binx)&(maskC>0))] #only values on IS (mask) and within elevation range
        hypso[i]=np.nanmean(meanx[tuple(v)])#*areax/1e12    --> mean of tp
        areas[i]=np.nansum(v)
        # if i== 8: 
        #     plt.imshow(np.array(v)[0,:,:])
        #     # plt.xlim([1000,1400])
        #     # plt.ylim([0,200])
        #     plt.colorbar()
        # print(i,elevs[i],np.sum(v),hypso[i],areas[i])
        # if i== 55: 
        #     plt.imshow(np.array(v)[0,:,:])
    areas*=areax
    
    #statistics
    if jj==0: areas_all=np.zeros((n,len(label)-1))
    areas_all[:,jj]=areas
    print(hypso[11])
    
    
    #plot parameters
    ax1.plot(hypso,elevs, mark[0], color=color[jj]) # hypsometric curve
    ax2.plot(areas/1e6,elevs, mark[1], color=color[jj]) # area-altitude distribution
    ax1.set_ylabel('elevation, m')
    ax1.set_xlabel('average total precipitation \n 1998-2020, mm y$^{-1}$', color='k') #tab:blue
    ax2.set_xlabel('area, km${^2}$', color='dimgrey') #tab:orange
    ax2.spines['bottom'].set_color('k') #tab:blue
    ax2.spines['top'].set_color('dimgrey') #tab:orange
    ax1.tick_params(axis='x', colors='k') #C0
    ax2.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    ax2.tick_params(axis='x', colors='dimgrey') #tab:orange
    ax1.xaxis.grid(color='k', alpha=0.8) # vertical lines
    ax2.xaxis.grid(linestyle='--') # vertical lines
    # Hide x labels and tick labels 
    ax1.label_outer()

# print('Ratio ERA/CARRA', (np.nansum(tpE*maskC)/np.nansum(tpC*maskC)))
#add legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color=color[0], linestyle=mark[0], lw=1),
                Line2D([0], [0], color=color[1], linestyle=mark[0], lw=1), 
                Line2D([0], [0], color=color[2], linestyle=mark[0], lw=1),
                Line2D([0], [0], color=color[3], linestyle=mark[0], lw=1),
                Line2D([0], [0], color=color[4], linestyle=mark[0], lw=1)]
custom_lines2 = [Line2D([0], [0], color='k', linestyle=mark[0], lw=1),
                Line2D([0], [0], color='k', linestyle=mark[1], lw=1)]
# ax1.legend(custom_lines, label,  loc='center left', frameon=False)
# ax2.legend(custom_lines2, ['tp','area'],  loc='upper right', frameon=False)

#save file
plt.savefig('./Figs/tp1998to2020_hypso_elev', bbox_inches='tight', dpi=300)


#---------------------- cumulative map -------------


for jj in range(0,len(label)-1):  #len(label)#CARRA and ERA5
    if jj==0: meanx=tpC.copy()*maskC; elev=elevC.copy() #don't mask elevation with *maskC since that is decreasing elevation values (due to range between 0 and 1)
    if jj==1: meanx=tpE.copy()*maskC; elev=elevE.copy()
    if jj==2: meanx=tpM.copy()*maskC; elev=elevM.copy()
    if jj==3: meanx=tpN.copy()*maskC; elev=elevN.copy()
    if jj==4: meanx=tpR.copy()*maskC; elev=elevR.copy()
    if jj==5: meanx=tpJ.copy()*maskC; elev=elevJ.copy()
    
    for i in range(n):
        v=[((elev>elevs[i])&(elev<=elevs[i]+binx)&(maskC>0))] #only values on IS (mask) and within elevation range
        tp_cum[i]=np.nansum(meanx[tuple(v)])*areax/1e12    #tuple(v)
        areas[i]=np.nansum(v)
    areas*=areax
    
    ax3.plot(np.cumsum(tp_cum), elevs, mark[0], color=color[jj])
    ax4.plot(np.cumsum(areas)/1e6, elevs, mark[1], color=color[jj])
    # ax3.plot(tp_cum, elevs, mark[0], color=color[jj])  #tp_cum (sum of tp per bin)
    # ax4.plot(areas/1e6, elevs, mark[1], color=color[jj])
    # ax2.spines['bottom'].set_color('r') #tab:blue
    # ax1.tick_params(axis='x', colors='r') #C0
    # plt.legend(handles=handles, frameon=False, fontsize=fs)

    idx1 = (np.abs(np.cumsum(tp_cum) - (np.median(tp_cum)*len(tp_cum)/2))).argmin() #find elev nearest to 50% mark > median
    # idx1 = (np.abs(np.cumsum(tp_cum) - (np.cumsum(tp_cum)[-1]/2))).argmin() #find elev nearest to 50% mark ->mean
    idx2 = (np.abs(np.cumsum(areas)/1e6 - (np.median(areas)*len(areas)/(2*1e6)))).argmin() #find elev nearest to 50% mark -> median
    # idx2 = (np.abs(np.cumsum(areas)/1e6 - (np.cumsum(areas)[-1]/(2*1e6)))).argmin() #find elev nearest to 50% mark -> mean
    ax3.axhline(elevs[idx1],  linestyle=mark[0], color=color[jj], alpha=0.3) #50%
    ax4.axhline(elevs[idx2],  linestyle=mark[1], color=color[jj], alpha=0.3) #50%
    # print(jj, idx1,idx2, tp_cum[idx1], tp_cum[idx2] )

    

ax4.tick_params(axis='x', colors='dimgrey') 
ax4.tick_params(axis='x', colors='dimgrey')  
ax4.spines['top'].set_color('dimgrey') #tab:orange 
ax3.set_xlabel("cumulative average total precipitation \n 1998-2020, Gt y$^{-1}$", fontsize=fs)
ax4.set_xlabel("cumulative area, $km^2$", fontsize=fs, color='dimgrey')
ax3.set_ylabel("elevation, m", fontsize=fs)
ax3.xaxis.grid(color='k', alpha=0.8) # vertical lines
ax4.xaxis.grid(linestyle='--') # vertical lines
ax4.set_xlim([-110000,2490000])
ax3.label_outer()


ax3.legend(custom_lines, label,  loc='lower right', frameon=True)
ax4.legend(custom_lines2, ['tp','area'],  loc='upper left', frameon=True)


mult=1.4
#text (a) and (b)
ax1.text(-0.15,1.13,'(a)',transform=ax1.transAxes, fontsize=fs*mult,
    verticalalignment='top',rotation=0,color='k', rotation_mode="anchor")
ax3.text(1.01,1.13,'(b)',transform=ax1.transAxes, fontsize=fs*mult,
    verticalalignment='top',rotation=0,color='k', rotation_mode="anchor")
mult=1
#text median
ax3.text(1.01,.66,'50th percentiles',transform=ax1.transAxes, fontsize=fs*mult,
    verticalalignment='top',rotation=0,color='k', rotation_mode="anchor")


# # #add in-situ elevation data
fn='./ancil/Greenland_snow_pit_SWE_v20210626.xlsx'
df_insitu = pd.read_excel(fn, engine='openpyxl')
df_insitu = df_insitu.drop(df_insitu[df_insitu["End Year"] <1997].index)
df_insitu = df_insitu.drop(df_insitu[df_insitu.Name=='Basin5'].index)
df_insitu = df_insitu.drop(df_insitu[df_insitu.Name=='JAR2'].index)
df_insitu = df_insitu.drop_duplicates(subset="Name")
thresh_elev=1900
frac1=df_insitu[df_insitu.Elevation>thresh_elev]
frac2=df_insitu[df_insitu.Elevation<thresh_elev]
print(len(frac1)/len(df_insitu), len(frac2)/len(df_insitu))

# #get in-situ data per elevation bin
# pits=np.zeros(n)*np.nan
# for i in range(n):
#     v=[((df_insitu.Elevation>elevs[i])&(df_insitu.Elevation<=elevs[i]+binx))] 
#     pits[i]=np.nansum(v)
# #plot 
# ax3 = fig.add_subplot(121)
# ax3.plot(pits, elevs)
# for i, v in enumerate(pits):
#     ax3.text(v, elevs[i], str(v), color='k', fontsize=fs)
# ax3.axes.xaxis.set_ticklabels([])
# ax3.axes.xaxis.set_ticks([])


#%% tp vs area relative impact
r_area=areas/np.nanmax(areas)
# r_tp=hypso/np.nanmax(hypso) #mean values of tp per bin
r_tp=tp_cum/np.nanmax(tp_cum) #absolute Gt of tp mass per bin

plt.plot(r_area, elevs, label='ratio of area')
plt.plot(r_tp, elevs, label='ratio of tp')
plt.plot(r_area*r_tp, elevs, label='ratio of area*tp')

plt.xlabel("relative impact of precipitation ", fontsize=fs, color='k')
plt.ylabel("elevation, m", fontsize=fs, color='k')

plt.legend()


#%% statistics, area model agreement above certain elevation
elev_thresh=750
areas_all2=areas_all.copy()
areas_all2[elevs<elev_thresh]=np.nan
areas_all2[-1]=np.nan #remove last line values
areas_mean=areas_all2.mean(axis=1)
areas_perc=np.zeros((n,len(label)-1)) 
for i in range(len(areas_all2)):  #calculate deviation (%) from mean
    # areas_diff[i]=areas_all[i,:]-areas_mean[i]
    areas_perc[i]=areas_all2[i,:]/areas_mean[i]
max_perc=np.nanmax(areas_perc)
min_perc=np.nanmin(areas_perc)
median_perc=np.nanmedian(areas_perc)
print('max: +',(max_perc-1)*100, 'min: -', (1-min_perc)*100, 'median: ', median_perc )
    
#%% Elevation/TP model difference with CARRA

#what to plot
varnam=['elevation', 'total precipitation']
varnam_short=['elev', 'tp']
units=['m', 'mm']
models=['ERA 5', 'MAR', 'NHM-SMAP', 'RACMO']

# in-situ locations
fn='./ancil/Table_1.csv'
gdf = pd.read_csv(fn)
gdf = gdf.drop(gdf[gdf.Program == "Kjær et al, unpublished 2021"].index)
gdf = gdf.drop(gdf[gdf.Program == "Kjær"].index)
gdf = gdf.reset_index(drop=True)

means=np.zeros((len(models)))
medians=np.zeros((len(models)))
rmss=np.zeros((len(models)))

for jj in range(len(varnam_short)): 
    fig, ax = plt.subplots(1, 4)
    # fig = plt.figure(figsize=(8,15.25)) # Notice aspect ratio (not equal) #8,9.75
    for ii in range(len(models)):
        # graph=221+ii #graph box number, 2 by 2 graph
        graph=141+ii #graph box number, 4 in a row graph
        plt.subplot(graph) 
        if jj==0: threshold=200
        if jj==1: threshold=500
        #read data
        varCARRA=eval(varnam_short[jj]+'C') #CARRA data
        varMODEL=eval(varnam_short[jj]+list(models[ii])[0]) #model data
        diff=varMODEL-varCARRA  #model-CARRA
        #plot data
        cmap = mpl.cm.get_cmap("bwr").copy()
        p3=m.imshow(diff*maskC.copy(), vmax=threshold, vmin=-threshold, cmap=cmap)
        p3.cmap.set_under('purple')
        p3.cmap.set_over('yellow')
        m.drawcoastlines(color='k',linewidth=.08)
        
        #add in-situ stations
        lon1, lat1 = m(*(gdf.lon, gdf.lat)) #transform into m map projection
        m.scatter(lon1, lat1,c='k', s=.5)
        
        #statistics
        stat=diff.copy()
        stat[np.isnan(maskC.copy())]=np.nan #mask
        # stat[elevC>2000]=np.nan #ice edge only filter
        rms = np.sqrt(np.nanmean(np.square(stat)))
        print(models[ii], ': mean=', "%0.0f" % np.nanmean(stat), ', median=', "%0.0f" % np.nanmedian(stat), ', rms=', "%0.0f" % rms)
        means[ii]='%.0f'% np.nanmean(stat)
        medians[ii]='%.0f'% np.nanmedian(stat)
        rmss[ii]='%.0f'% np.sqrt(np.nanmean(np.square(stat)))
        
        #add text model
        xx0=.7; yy0=.07
        # xx0=1500000; yy0=350000
        mult=0.7
        plt.text(xx0, yy0, models[ii]+'\n minus \n CARRA', fontsize=fs*mult,ha='center', transform=ax[ii].transAxes, color='k') 

        # zooming in
        xlim0=400000 ; xlim1= 1950000 ; ylim0=120000; ylim1=3000000 ; 
        plt.xlim(xlim0,xlim1)        
        plt.ylim(ylim0,ylim1)
        


    # plt.suptitle(varnam[jj])
    
    #add colorbar
    cb_ax = fig.add_axes([.126, 0.2, 0.775, .02])
    # cb_ax = fig.add_axes([.92, 0.233, 0.01, .54])
    # cb_ax = fig.add_axes([.92, 0.125, 0.02, .755]) #2by2 graph
    bounds=[-threshold, -int(threshold/2), 0, int(threshold/2), threshold]
    cbar2 = plt.colorbar(p3, orientation='horizontal', cax=cb_ax, extend='both', ticks=bounds)
    # cbar2.ax.set_xticklabels(bounds, fontsize=fs*mult)
    cbar2.outline.set_linewidth(.2)
    cbar2.ax.tick_params(width=.2)
    cbar2.ax.tick_params(labelsize=fs*mult)
    cbar2.ax.set_xlabel(varnam[jj]+' difference, '+units[jj], fontsize=fs*mult)
    plt.subplots_adjust(wspace=0, hspace=0) #no white sace between subplots
    
    #set outline box linewidth
    for axis in ['top','bottom','left','right']:
        ax[0].spines[axis].set_linewidth(0.2)
        ax[1].spines[axis].set_linewidth(0.2)
        ax[2].spines[axis].set_linewidth(0.2)
        ax[3].spines[axis].set_linewidth(0.2)
        
# statistics to csv
df_stat = pd.DataFrame(models)
df_stat['elev_mean']= means
df_stat['elev_median']=medians
df_stat['elev_rms']=rmss
df_stat.to_csv(path+'output_csv/hypsometry_statistics.csv')

#%% statistics
# test=diff.copy()
# test[np.isnan(maskC)]=np.nan
# test[elevN>2000]=np.nan

# # plt.imshow(ttest)
# # plt.colorbar()
# print(models[ii], ': mean=', np.nanmean(ttest), ', median=', np.nanmedian(ttest))

m1=tpN.copy()
mC=tpC.copy()
#mask
m1*=maskC.copy()
mC*=maskC.copy()
# m1[np.isnan(maskC)]=np.nan
# mC[np.isnan(maskC)]=np.nan
print("%0.0f" % np.nanmean(m1), "%0.0f" % np.nanmean(mC), 1-(np.nanmean(m1)/np.nanmean(mC)) )

#in Gt
areax=2.5e3**2
v=[~np.isnan(maskC)]
areas_number=np.sum(v)

m1_Gt=np.nansum(m1)*areax /1e12
mC_Gt=np.nansum(mC)*areax /1e12
print("%0.0f" % np.nanmean(m1_Gt), "%0.0f" % np.nanmean(mC_Gt), (np.nanmean(m1_Gt)-np.nanmean(mC_Gt)) )

#%% Compare two random models

#what to plot
model1='CARRA'
model2='MAR'
varnam=['elevation', 'total precipitation']
# var=[elevC, elevE, tpC, tpE]
var=[eval('elev'+list(model1)[0]), eval('elev'+list(model2)[0]), eval('tp'+list(model1)[0]), eval('tp'+list(model2)[0])]

#start plotting
label=[model1, model2, model2+' - '+model1]
dd=0
for jj in range(len(varnam)): 
    fig, ax = plt.subplots(1, 3)
        
    if jj==0: threshold=3250
    if jj==0: threshold=2000
    plt.subplot(131) #model1
    p1=m.imshow(var[jj+dd]*maskC, vmax=threshold, vmin=0)
    m.drawcoastlines(color='k',linewidth=.08)
    
    plt.subplot(132) #model2
    p2=m.imshow(var[jj+dd+1]*maskC, vmax=threshold, vmin=.08)
    m.drawcoastlines(color='k',linewidth=.1)
    
    plt.subplot(133) #model2-model1
    threshold=200
    diff=var[jj+dd+1]-var[jj+dd]
    p3=m.imshow(diff*maskC, vmax=threshold, vmin=-threshold, cmap='bwr')
    m.drawcoastlines(color='k',linewidth=.08)
    # plt.title('ERA-CARRA total precipitation')
    
    
    #add colorbars
    mult=0.7
    ax_cbar1 = fig.add_axes([0.126, 0.24, .51, 0.02])
    if jj==0: bounds=[1,1000, 2000, 3000]
    if jj==1: bounds=[1,500,1000,1500, 2000]
    ex=['neither','max']
    cbar1 = plt.colorbar(p2, orientation='horizontal', cax=ax_cbar1, extend=ex[jj], ticks=bounds)
    if jj==0: cbar1.ax.set_xlabel('elevation, m', fontsize = fs*mult)
    if jj==1: cbar1.ax.set_xlabel('average total precipitation 1998-2020, mm', fontsize = fs*mult)
    bounds[0]=0
    cbar1.ax.set_xticklabels(bounds, fontsize=fs*mult)
    cbar1.ax.tick_params(width=.4)
    # cb_ax2 = fig.add_axes([.92, 0.275, 0.015, 0.455])
    cb_ax2 = fig.add_axes([0.6475, 0.24, 0.25, 0.02])
    bounds=[-200, -100, 0,100, 200]
    cbar2 = plt.colorbar(p3, orientation='horizontal', cax=cb_ax2, extend='both', ticks=bounds)
    cbar2.ax.set_xticklabels(bounds, fontsize=fs*mult)
    if jj==0: cbar2.ax.set_xlabel('difference, m', fontsize = fs*mult)
    if jj==1: cbar2.ax.set_xlabel('difference, mm', fontsize = fs*mult)
    cbar2.ax.tick_params(width=.4)
    p2.cmap.set_over('magenta')
    p3.cmap.set_under('purple')
    p3.cmap.set_over('yellow')
    cbar1.outline.set_linewidth(.4)
    cbar2.outline.set_linewidth(.4)
    
    #set outline box linewidth
    for axis in ['top','bottom','left','right']:
        ax[0].spines[axis].set_linewidth(0.5)
        ax[1].spines[axis].set_linewidth(0.5)
        ax[2].spines[axis].set_linewidth(0.5)
    
    #add text model
    for i in range(3): 
        xx0=.7; yy0=.10
        plt.text(xx0, yy0, label[i], fontsize=fs*mult,ha='center', transform=ax[i].transAxes, color='k') 
    
    plt.subplots_adjust(wspace=0, hspace=0) #no white sace between subplots
    
    dd+=1
    
    #save file
    plt.savefig('./Figs/'+varnam[jj]+'CARRAvsERA5', bbox_inches='tight', dpi=300)

#%% mask peri 4 areas

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()  #second x-axis

#define model dataset and mask accordingly
tp=tpC*mask_peri
elev=elevC.copy()
elev[mask_peri==0]=0


#split in 4 regions
thresh_lat=70
thresh_lon=-44+360

tp_sw=tp.copy()
tp_nw=tp.copy()
tp_se=tp.copy()
tp_ne=tp.copy()
elev_sw=elev.copy()
elev_nw=elev.copy()
elev_se=elev.copy()
elev_ne=elev.copy()
#southwest
v=np.logical_or(dsC.latitude>thresh_lat, dsC.longitude<thresh_lon) #what to cut
tp_sw[v]=np.nan
elev_sw[v]=np.nan
#northwest
v=np.logical_or(dsC.latitude<thresh_lat, dsC.longitude<thresh_lon) #what to cut
tp_nw[v]=np.nan
elev_nw[v]=np.nan
#southeast
v=np.logical_or(dsC.latitude>thresh_lat, dsC.longitude>thresh_lon) #what to cut
tp_se[v]=np.nan
elev_se[v]=np.nan
#northeast
v=np.logical_or(dsC.latitude<thresh_lat, dsC.longitude>thresh_lon) #what to cut
tp_ne[v]=np.nan
elev_ne[v]=np.nan

tp_all=[tp_sw,tp_se, tp_nw, tp_ne]
elev_all=[elev_sw, elev_se, elev_nw, elev_ne]
# plt.imshow(tp_ne*maskC)

mark=['-', '--', ':', '-.']

for k in range(4):
    meanx=tp_all[k].copy()*mask_peri; 
    elev=elev_all[k].copy() #don't mask elevation with *maskC since that is decreasing elevation values (due to range between 0 and 1)
    binx=50
    elevs=np.arange(0,3250,binx)
    n=len(elevs)
    
    
    areax=2.5e3**2
    # if jj==1: areax=9e3**2
    
    hypso=np.zeros(n)*np.nan
    areas=np.zeros(n)*np.nan
    
    for i in range(n):
        v=[((elev>elevs[i])&(elev<=elevs[i]+binx)&(maskC>0))] #only values on IS (mask) and within elevation range
        hypso[i]=np.mean(meanx[v])#*areax/1e6)    #tuple(v)
        areas[i]=np.sum(v)

    areas*=areax
    
    ax1.plot(hypso,elevs, mark[k], color='C0') # hypsometric curve
    ax2.plot(areas/1e6,elevs, mark[k], color='tab:orange') # area-altitude distribution
    ax1.set_ylabel('elevation, m')
    ax1.set_xlabel('average total precipitation 1998-2020, mm y$^{-1}$', color='tab:blue')
    ax2.set_xlabel('area, km${^2}$', color='tab:orange')
    ax2.spines['bottom'].set_color('tab:blue')
    ax2.spines['top'].set_color('tab:orange')
    ax1.tick_params(axis='x', colors='C0')
    ax2.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    ax2.tick_params(axis='x', colors='tab:orange')

#add legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='k', linestyle=mark[0], lw=1),
                Line2D([0], [0], color='k', linestyle=mark[1], lw=1),
                Line2D([0], [0], color='k', linestyle=mark[2], lw=1),
                Line2D([0], [0], color='k', linestyle=mark[3], lw=1)]
label=['southwest', 'southeast','northwest',  'northeast']
ax1.legend(custom_lines, label, loc='upper right', frameon=False)

#%% produce new CARRA elev nc file
# from netCDF4 import Dataset,num2date
# ni=1269 ; nj=1069
# outpath='C:/Users/Armin/Documents/Work/GEUS/Github/CARRA/'
# ofile=outpath+'elev_lat_lon1.nc'
# test=CARRAz.values
# n_days=len(CARRAz.values)
# print("start making .nc file")
# # os.system("/bin/rm "+ofile)
# ncfile = Dataset(ofile,mode='w',format='NETCDF4_CLASSIC')
# lat_dim = ncfile.createDimension('lat', nj)     # latitude axis
# lon_dim = ncfile.createDimension('lon', ni)    # longitude axis
# # time_dim = ncfile.createDimension('time', n_days) # unlimited axis (can be appended to)

# ncfile.subtitle="subtitle"


# latitude = ncfile.createVariable('latitude', np.float32, ('lon','lat'))
# latitude.units = 'degrees_north'
# latitude.long_name = 'latitude'
# longitude = ncfile.createVariable('longitude', np.float32, ('lon','lat'))
# longitude.units = 'degrees_east'
# longitude.long_name = 'longitude'
# # time = ncfile.createVariable('time', np.float64, ('time',))
# # time.units = 'days since '+year+'-01-01'
# # time.units = 'days'
# # time.long_name = 'time'
# # Define a 3D variable to hold the data
# print("compressing")
# temp = ncfile.createVariable('z',np.float32,('lon','lat'),zlib=True,least_significant_digit=3) # note: unlimited dimension is leftmost
# temp.units = 'm' # degrees Kelvin
# temp.standard_name = 'z' # this is a CF standard name

# nlats = len(lat_dim); nlons = len(lon_dim); ntimes = 3

# latitude[:,:]=lat_mat
# longitude[:,:]=lon_mat
# temp[:,:] = test  # Appends data along unlimited dimension


# print("-- Wrote data, temp.shape is now ", temp.shape)
# print(ofile)
   
# ncfile.close(); print('Dataset is closed!')