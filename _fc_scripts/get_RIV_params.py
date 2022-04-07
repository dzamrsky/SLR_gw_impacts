# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 14:19:34 2019

@author: daniel


Script that loads various datasets and calculates the RIV stage, bot elev and COND.

"""

import os
import dask.array as da
from IPython.display import Image
import xarray as xa
from dask.diagnostics import ProgressBar
from netCDF4 import Dataset
import numpy as np


#   define the different raster files and databases that will be used to derive the various parameter values
flo1k_dir = r'g:\_ORIGINAL_DATA\FLO1K\FLO1K_qav_mean_1km.nc'













qav_dir_5km = r'd:\_ORIGINAL_DATA\FLO1K\FLO1K.5min.ts.1960.2015.qav.nc'
data = xa.open_dataset(qav_dir_5km, chunks={'lat' : 1080, 'lon' : 2160}, engine = 'netcdf4', mask_and_scale = -0.0)
data = xa.where(data['qav']>0, data['qav'], np.nan)

with ProgressBar():
    mean_qav_5km = data.mean(dim = ('time'), skipna = True)
    new_arr_5km = mean_qav_5km.values


qav_dir = r'd:\_ORIGINAL_DATA\FLO1K\FLO1K.5min.ts.1960.2015.qav.nc'
data = xa.open_dataset(qav_dir, chunks={'lat' : 2160, 'lon' : 4320})

mean_qav_2 = data.mean(dim = ('time'))
new_arr_2 = mean_qav_2.qav.values

xa_name = r'd:\_ORIGINAL_DATA\FLO1K\FLO1K_5km_qav_mean.nc'

ncout = Dataset(xa_name,'w','NETCDF4'); # using netCDF3 for output format 
ncout.createDimension('lon', 4320)
ncout.createDimension('lat', 2160)
lonvar = ncout.createVariable('lon','float32',('lon'))
lonvar[:] = data['lon'].values
latvar = ncout.createVariable('lat','float32',('lat'))
latvar[:] = data['lat'].values

qav_mean = ncout.createVariable('qav_mean','float32',('lat','lon'))
qav_mean.setncattr('units','m3/s')
qav_mean[:] = new_arr_2

ncout.close()

new_arr_2.to_netcdf(xa_name)      








