#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 08:52:21 2021

@author: AD
"""
import numpy as np
import os
from glob import glob
from netCDF4 import Dataset,num2date
import pandas as pd
from datetime import datetime 
import calendar
import netCDF4
import xarray as xr
from cftime import num2date, date2num
import matplotlib.pyplot as plt

years=np.arange(1998,2021).astype('str')

ni=1269
nj=1069

niM=181
njM=101

   
AD=1
path='/Users/jason/Dropbox/CARRA/prog/map_CARRA_west/'
if AD:
    path='C:/Users/Armin/Documents/Work/GEUS/Github/CARRA_rain/'
    raw_path='C:/Users/Armin/Documents/Work/GEUS/Github/CARRA/'
os.chdir(path)


tp_M_all=np.zeros((len(years), niM, njM))

for yy,year in enumerate(years):
    fn=raw_path+'MARv3.12.0.4-noBS-15km-daily-ERA5-'+year+'.nc'
    ds=xr.open_dataset(fn)
    
    # sf=np.array(ds.SF[:, :, :]) #snowfall in mmWE/day
    # rf=np.array(ds.RF[:, :, :]) #snowfall in mmWE/day
    smb=np.array(ds.SMB[:, 0,:, :]) #surface mass balance in mmWE/day
    su=np.array(ds.SU[:, 0,:, :]) #sublimation form snow/soil in mmWE/day
    ru=np.array(ds.RU[:, 0,:, :]) #runoff of meltwater and rainwater in mmWE/day
    tp=smb+ru+su  #definition of tp through smb
    tp_M=np.sum(tp, axis=0) #mmWE/year
    tp_M_all[yy,:,:]=tp_M


#other parameters
elev=np.array(ds.SH[:, :]) #surface height in m
latM=np.array(ds.LAT[:, :]) #latitude
lonM=np.array(ds.LON[:, :]) #longitude


#%% visualistaion check
plt.imshow(tp_M[:,:])
plt.colorbar()


#%% produce new MAR tp1998-2020 and elev nc file

ni=niM ; nj=njM
varnam='tp'

outpath='C:/Users/Armin/Documents/Work/GEUS/Github/CARRA/'
ofile=outpath+'MARv3.12.0.4-noBS-15km-yearly-ERA5-1998-2020.nc'
tp_out=tp_M_all
n_years=len(tp_M_all)
print("start making .nc file")
# os.system("/bin/rm "+ofile)
ncfile = Dataset(ofile,mode='w',format='NETCDF4_CLASSIC')
lat_dim = ncfile.createDimension('lat', nj)     # latitude axis
lon_dim = ncfile.createDimension('lon', ni)    # longitude axis
time_dim = ncfile.createDimension('time', n_years) # unlimited axis (can be appended to)

ncfile.subtitle="subtitle"


latitude = ncfile.createVariable('latitude', np.float32, ('lon','lat'))
latitude.units = 'degrees_north'
latitude.long_name = 'latitude'
longitude = ncfile.createVariable('longitude', np.float32, ('lon','lat'))
longitude.units = 'degrees_east'
longitude.long_name = 'longitude'
time = ncfile.createVariable('time', np.float64, ('time',))
time.units = 'years'
time.long_name = 'time'
# Define a 3D variable to hold the data
print("compressing")
temp = ncfile.createVariable(varnam,np.float32,('time','lon','lat'),zlib=True,least_significant_digit=3) # note: unlimited dimension is leftmost
temp.units = 'mm/year'
temp.standard_name = varnam # this is a CF standard name

temp2 = ncfile.createVariable('z',np.float32,('lon','lat'),zlib=True,least_significant_digit=3) # note: unlimited dimension is leftmost
temp2.units = 'm'
temp2.standard_name = 'z' # this is a CF standard name


nlats = len(lat_dim); nlons = len(lon_dim); ntimes = 3

latitude[:,:]=latM
longitude[:,:]=lonM
time=years
temp[:,:,:] = tp_out  # Appends data along unlimited dimension
temp2[:,:] = elev  # Appends data along unlimited dimension


print("-- Wrote data, temp.shape is now ", temp.shape)
print(ofile)
   
ncfile.close(); print('Dataset is closed!')
 
