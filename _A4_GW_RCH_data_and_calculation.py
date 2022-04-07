# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 10:55:12 2019

@author: daniel
"""

import pandas as pd
from sklearn import linear_model
import statsmodels.api as sm
import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import gdal
from osgeo import osr
import xarray as xr
import os
import math
from sklearn.metrics import mean_squared_error
from math import sqrt

from mpl_toolkits.basemap import Basemap
import matplotlib
from matplotlib.colorbar import ColorbarBase

##  First the GW_Rch data provided by Chinchu.
gw_rch_dir = r'g:\_ORIGINAL_DATA\_GW_Recharge_Chinchu_Mohan\Recharge_bestModel.mat'
csv_in_dir = r'g:\_ORIGINAL_DATA\_GW_Recharge_Chinchu_Mohan\_gw_rch_input_data.csv'

#   First transform the .mat grid into tiff 
gw_rch_best_model = scipy.io.loadmat(gw_rch_dir)['Recharge_bestModel']
gw_rch_avg = np.mean(gw_rch_best_model, axis = 2)
gw_rch_dir_out = r'g:\_ORIGINAL_DATA\_GW_Recharge_Chinchu_Mohan\Recharge_bestModel.tiff'

xmin, ymin, xmax, ymax = -180., -90., 180., 90.
nrows, ncols = np.shape(gw_rch_avg)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(gw_rch_dir_out, gw_rch_avg.shape[1], gw_rch_avg.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(gw_rch_avg)
output_raster.FlushCache()
output_raster = None


##  Use the historical land use data (.NC) files to extract land use at different time intervals in history.
##  Those will be at 1500 AD, 1600AD, 1700AD, 1800AD, 1900AD and 2000AD. The land use datasets contain fractions
##  of each land use in the given raster cell (50km * 50km). We use the u2 dataset that includes urban areas as well.
##  For the paleo gw recharge we will use the u1 dataset without urban areas - (the u1 at 1500 will simulate the 
##  land use conditions at 1400AD and before that). 

out_dir = r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\_A4_data'

u2_crop = xr.open_dataset(r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\data\LUHa_u2.v1_gcrop.nc4', decode_times = False)
u2_primary_land = xr.open_dataset(r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\data\LUHa_u2.v1_gothr.nc4', decode_times = False)
u2_pasture = xr.open_dataset(r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\data\LUHa_u2.v1_gpast.nc4', decode_times = False)
u2_secondary_land = xr.open_dataset(r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\data\LUHa_u2.v1_gsecd.nc4', decode_times = False)
u2_urban = xr.open_dataset(r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\data\LUHa_u2.v1_gurbn.nc4', decode_times = False)

##  the land use coefficient will be calculated as follows (based on Chinchu's method and classification)
##  lu_coef = 5 * crop_pct + 4 * (past_pct + primary_pct) + 3 * secondary_pct + 2 * urb_pct + 1 * (100 - crop_pct - past_pct - second_pct - primary_pct - urb_pct) 
##  during these time steps:
time_out = [0, 100, 200, 300, 400, 500]

for t in time_out:
    #   select the right arrays from all the u2 datasets
    crop = u2_crop.sel(time = t)['prop_crop'].values
    primary_land = u2_primary_land.sel(time = t)['prop_primary'].values
    pasture = u2_pasture.sel(time = t)['prop_past'].values
    secondary_land = u2_secondary_land.sel(time = t)['prop_secd'].values
    urban = u2_urban.sel(time = t)['prop_urbn'].values
    #   calculate the LU coeffecient array and save it as tiff file
    lu_coef_arr = 5 * crop + 4 * pasture + 3 * (secondary_land + primary_land) + 2 * urban + 1 * (1 - crop - pasture - secondary_land - primary_land - urban) 
    lu_out_name = os.path.join(out_dir, 'LU_coeff_' + str(1500 + t) + '.tiff')

    #print(t, crop[120, 360], primary_land[120, 360], pasture[120, 360], secondary_land[120, 360], urban[120, 360])
    xmin, ymin, xmax, ymax = -180., -90., 180., 90.
    nrows, ncols = np.shape(lu_coef_arr)
    xres = (xmax - xmin)/float(ncols)
    yres = (ymax - ymin)/float(nrows)
    geotransform = (xmin, xres, 0, ymax, 0, -yres)  
    output_raster = gdal.GetDriverByName('GTiff').Create(lu_out_name, lu_coef_arr.shape[1], lu_coef_arr.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(4326)      
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(lu_coef_arr)
    output_raster.FlushCache()
    output_raster = None

##  create the 1400 dataset - which is the 1500 year without urban pct
u1_crop = xr.open_dataset(r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\data\LUHa_t1.v1_gcrop.nc4', decode_times = False)
u1_primary_land = xr.open_dataset(r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\data\LUHa_t1.v1_gothr.nc4', decode_times = False)
u1_pasture = xr.open_dataset(r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\data\LUHa_t1.v1_gpast.nc4', decode_times = False)
u1_secondary_land = xr.open_dataset(r'g:\_ORIGINAL_DATA\_Harmonized_Global_Land_Use_1500_2100\data\LUHa_t1.v1_gsecd.nc4', decode_times = False)

crop = u1_crop.sel(time = 0)['prop_crop'].values
primary_land = u1_primary_land.sel(time = 0)['prop_primary'].values
pasture = u1_pasture.sel(time = 0)['prop_pasture'].values
secondary_land = u1_secondary_land.sel(time = 0)['prop_secd'].values
#   calculate the LU coeffecient array and save it as tiff file
lu_coef_arr = 5 * crop + 4 * pasture + 3 * (secondary_land + primary_land) + 2 * urban + 1 * (1 - crop - pasture - secondary_land - primary_land - urban) 
lu_out_name = os.path.join(out_dir, 'LU_coeff_' + str(1400) + '.tiff')

#print(t, crop[120, 360], primary_land[120, 360], pasture[120, 360], secondary_land[120, 360], urban[120, 360])
xmin, ymin, xmax, ymax = -180., -90., 180., 90.
nrows, ncols = np.shape(lu_coef_arr)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(lu_out_name, lu_coef_arr.shape[1], lu_coef_arr.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(lu_coef_arr)
output_raster.FlushCache()
output_raster = None




##  Open the CRU climatological datasets - the PET and P data that were used by Chinchu to calculate her GW recharge estimates.
##  We will create annual average values for each of these 3 variables and store them in one netcdf file. This will be used to 
##  recreate the GW recharge values from Chinchu (to double check that our approach is correct) and also to calculate the T/PET 
##  ratio that will be later used as input into Rens's equation to calculate the paleo PET data. 
out_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\CRU_climate_data.nc'

pet_80s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.1981.1990.pet.dat.nc') 
pet_90s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.1991.2000.pet.dat.nc') 
pet_00s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.2001.2010.pet.dat.nc') 
pet_10s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.2011.2018.pet.dat.nc') 

prec_80s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.1981.1990.pre.dat.nc') 
prec_90s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.1991.2000.pre.dat.nc') 
prec_00s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.2001.2010.pre.dat.nc') 
prec_10s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.2011.2018.pre.dat.nc') 

tmp_80s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.1981.1990.tmp.dat.nc') 
tmp_90s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.1991.2000.tmp.dat.nc') 
tmp_00s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.2001.2010.tmp.dat.nc') 
tmp_10s = xr.open_dataset(r'g:\_ORIGINAL_DATA\_CRU_climate_data\cru_ts4.03.2011.2018.tmp.dat.nc') 

#   create empty lists that will be filled with arrays of each variable, to be stored in the final nc file
pet_lst, tmp_lst, tmp_K_lst, pre_lst = [], [], [], []

#   start with the 80s data
years = [1984, 1985, 1986, 1987, 1988, 1989, 1990]
for year in years:
    pet = pet_80s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['pet'].values
    tmp = tmp_80s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['tmp'].values
    pre = prec_80s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['pre'].values
    #   calculate the averages and transform from C to K 
    pet_avg = np.mean(pet, axis = 0)
    tmp_avg = np.mean(tmp, axis = 0)
    pre_avg = np.mean(pre, axis = 0)
    tmp_K = tmp_avg + 273.15
    pet_lst.append(pet_avg)
    tmp_lst.append(tmp_avg)
    pre_lst.append(pre_avg)
    tmp_K_lst.append(tmp_K)

#   the 90s data
years = [1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000]
for year in years:
    pet = pet_90s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['pet'].values
    tmp = tmp_90s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['tmp'].values
    pre = prec_90s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['pre'].values
    #   calculate the averages and transform from C to K 
    pet_avg = np.mean(pet, axis = 0)
    tmp_avg = np.mean(tmp, axis = 0)
    pre_avg = np.mean(pre, axis = 0)
    tmp_K = tmp_avg + 273.15
    pet_lst.append(pet_avg)
    tmp_lst.append(tmp_avg)
    pre_lst.append(pre_avg)
    tmp_K_lst.append(tmp_K)

#   the 00s data
years = [2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010]
for year in years:
    pet = pet_00s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['pet'].values
    tmp = tmp_00s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['tmp'].values
    pre = prec_00s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['pre'].values
    #   calculate the averages and transform from C to K 
    pet_avg = np.mean(pet, axis = 0)
    tmp_avg = np.mean(tmp, axis = 0)
    pre_avg = np.mean(pre, axis = 0)
    tmp_K = tmp_avg + 273.15
    pet_lst.append(pet_avg)
    tmp_lst.append(tmp_avg)
    pre_lst.append(pre_avg)
    tmp_K_lst.append(tmp_K)

#   the 10s data
years = [2011, 2012, 2013, 2014]
for year in years:
    pet = pet_10s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['pet'].values
    tmp = tmp_10s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['tmp'].values
    pre = prec_10s.sel(time = slice(str(year) + '-01-01', str(year) + '-12-31'))['pre'].values
    #   calculate the averages and transform from C to K 
    pet_avg = np.mean(pet, axis = 0)
    tmp_avg = np.mean(tmp, axis = 0)
    pre_avg = np.mean(pre, axis = 0)
    tmp_K = tmp_avg + 273.15
    pet_lst.append(pet_avg)
    tmp_lst.append(tmp_avg)
    pre_lst.append(pre_avg)
    tmp_K_lst.append(tmp_K)

all_years = [1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000,\
             2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014]
nc_out = xr.Dataset(data_vars = {'PET_annual_avg' : (('year', 'lat', 'lon'), pet_lst),
                                 'PREC_annual_avg' : (('year', 'lat', 'lon'), pre_lst),
                                 'TMP_annual_avg' : (('year', 'lat', 'lon'), tmp_lst),
                                 'TMP_K_annual_avg' : (('year', 'lat', 'lon'), tmp_K_lst)},
                    coords = {'year' : all_years, 'lat' : pet_80s.coords['lat'].values, 'lon' : pet_80s.coords['lon'].values})
nc_out.to_netcdf(out_dir)    

##  based on Rens's suggestion and paper by Gardner et. al (2009), the relation between T and PET is :
##
##      (1) PET (mm/yr) = 1.2e10 * exp(-4620 / T (annual average T in Kelvins)) 
##
##  from this we can introduce the ratio between now and the past
##
##      (2) PET_past / PET_now = exp(-4620 / T_past) / exp(-4620 / T_now) 
##
##  so to calculate the past PET value we simply use the following formula
##
##      (3) PET_past = PET_now *[exp(-4620 / T_past) / exp(-4620 / T_now)] 

#   for the T_now in Kelvins take the average over 30 years (1984 - 2014)
t_now_K_mean = np.mean(np.array(tmp_K_lst), axis = 0)
pet_mean = np.mean(np.array(pet_lst) * 365.25, axis = 0)

pet_mean_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\pet_mean.nc'
xmin, ymin, xmax, ymax = -180., -90., 180., 90.
nrows, ncols = np.shape(pet_mean)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(pet_mean_dir, pet_mean.shape[1], pet_mean.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(pet_mean[::-1])
output_raster.FlushCache()
output_raster = None

pet_calc_mean = np.zeros((360, 720))
for i in range(360):
    for j in range(720):
        pet_calc_mean[i, j] = 1.2e10 * math.exp(-4620 / t_now_K_mean[i, j])

pet_calc_mean_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\pet_calc_mean.nc'
xmin, ymin, xmax, ymax = -180., -90., 180., 90.
nrows, ncols = np.shape(pet_calc_mean)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(pet_calc_mean_dir, pet_calc_mean.shape[1], pet_calc_mean.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(pet_calc_mean[::-1])
output_raster.FlushCache()
output_raster = None

pet_diff = pet_calc_mean - pet_mean
pet_diff_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\pet_DIFF.nc'
xmin, ymin, xmax, ymax = -180., -90., 180., 90.
nrows, ncols = np.shape(pet_diff)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(pet_diff_dir, pet_diff.shape[1], pet_diff.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(pet_diff[::-1])
output_raster.FlushCache()
output_raster = None





##  read in the comparison of different PET sources - the MODIS annual average (2000 - 2013), CRU average (1984 - 2014)
##  and the PET calculated based on average annual temperature (1984 - 1984).

csv_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\PET_comparison_GEBCO_elev_2.csv'
df = pd.read_csv(csv_dir)

#   get the different PET values at each point
pet_modis = list(df['MOD16A3_PE'].values / 10.)
pet_chinchu = list(df['PET'].values)
pet_calc_T = list(df['pet_calc_m'].values)
pet_cru = list(df['pet_mean'].values)
lat_chinchu = list(df['Latitude'].values)
temp_chinchu = list(df['T'].values)
gebco_elev = list(df['GEBCO_2014'].values / 1000.)

rms_modis = sqrt(mean_squared_error(pet_modis, pet_chinchu))
rms_cru = sqrt(mean_squared_error(pet_cru, pet_chinchu))
rms_cacl_T = sqrt(mean_squared_error(pet_calc_T, pet_chinchu))

#   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
fig = plt.figure(figsize=(11.69,8.27))
#   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
ax1 = plt.subplot2grid((3, 4), (0, 0))  #  Modis vs observed   
ax2 = plt.subplot2grid((3, 4), (0, 1))  #  Modis error vs latitude
ax3 = plt.subplot2grid((3, 4), (1, 0))  #  CRU annual mean vs observed 
ax4 = plt.subplot2grid((3, 4), (1, 1))  #  CRU error vs latitude
ax5 = plt.subplot2grid((3, 4), (2, 0))  #  Calculated by Temp vs observed
ax6 = plt.subplot2grid((3, 4), (2, 1))  #  Calculated by Temp error vs latitude
ax7 = plt.subplot2grid((3, 4), (0, 2))  #  CRU error vs temp
ax8 = plt.subplot2grid((3, 4), (1, 2))  #  Calculated by Temp vs temp
ax9 = plt.subplot2grid((3, 4), (2, 2))  #  Calculated by Temp error vs temp
ax10 = plt.subplot2grid((3, 4), (0, 3))  #  CRU error vs elevation (GEBCO)
ax11 = plt.subplot2grid((3, 4), (1, 3))  #  Calculated by Temp vs elevation (GEBCO)
ax12 = plt.subplot2grid((3, 4), (2, 3))  #  Calculated by Temp error vs elevation (GEBCO)
#   set the exact positions for each axis
ax1.set_position([0.05, 0.675, 0.21, 0.25]) # [left, bottom, width, height]
ax2.set_position([0.315, 0.675, 0.18, 0.25])
ax3.set_position([0.05, 0.375, 0.21, 0.25])        
ax4.set_position([0.315, 0.375, 0.18, 0.25]) # [left, bottom, width, height]
ax5.set_position([0.05, 0.075, 0.21, 0.25])
ax6.set_position([0.315, 0.075, 0.18, 0.25])   
ax7.set_position([0.55, 0.675, 0.18, 0.25]) # [left, bottom, width, height]
ax8.set_position([0.55, 0.375, 0.18, 0.25])      
ax9.set_position([0.55, 0.075, 0.18, 0.25])
ax10.set_position([0.785, 0.675, 0.18, 0.25])
ax11.set_position([0.785, 0.375, 0.18, 0.25]) # [left, bottom, width, height]
ax12.set_position([0.785, 0.075, 0.18, 0.25])     

ax1.scatter(pet_chinchu, pet_modis, color='green', marker='o', s = 5)
ax1.set_xlim([0, math.ceil(max(pet_chinchu) / 1000.) * 1000.])
ax1.set_ylim([0, math.ceil(max(pet_modis) / 1000.) * 1000.])
ax1.plot([0, 3000], [0, 3000], 'r--')
#ax1.text(1 , 1, 'RMSE = ' + str(round(rms_modis, 1)))
ax1.text(0.75, 0.1,'RMSE = ' + str(round(rms_modis, 1)), ha='center', va='center', transform = ax1.transAxes)
ax1.set_xlabel('Observed PET (mm/year)', fontsize = 8)
ax1.set_ylabel('MODIS PET (mm/year)', fontsize = 8)
ax1.tick_params(labelsize = 8)        

ax2.scatter([x1 - x2 for (x1, x2) in zip(pet_modis, pet_chinchu)], lat_chinchu, color = 'red', marker='o', s = 5)
ax2.set_xlim([0, math.ceil(max([x1 - x2 for (x1, x2) in zip(pet_modis, pet_chinchu)]) / 1000.) * 1000.])
ax2.set_ylim([-90, 90])
ax2.set_xlabel('MODIS - Observed PET (mm/year)', fontsize = 8)
ax2.set_ylabel('Latitude (degrees)', fontsize = 8)
ax2.tick_params(labelsize = 8)        

ax3.scatter(pet_chinchu, pet_cru, color='green', marker='o', s = 5)
ax3.set_xlim([0, math.ceil(max(pet_chinchu) / 1000.) * 1000.])
ax3.set_ylim([0, math.ceil(max(pet_cru) / 1000.) * 1000.])
ax3.plot([0, 3000], [0, 3000], 'r--')
ax3.text(0.75, 0.1,'RMSE = ' + str(round(rms_cru, 1)), ha='center', va='center', transform = ax3.transAxes)
ax3.set_xlabel('Observed PET (mm/year)', fontsize = 8)
ax3.set_ylabel('CRU PET (mm/year)', fontsize = 8)
ax3.tick_params(labelsize = 8)        

ax4.scatter([x1 - x2 for (x1, x2) in zip(pet_cru, pet_chinchu)], lat_chinchu, color = 'red', marker='o', s = 5)
ax4.set_xlim([0, math.ceil(max([x1 - x2 for (x1, x2) in zip(pet_cru, pet_chinchu)]) / 1000.) * 1000.])
ax4.set_ylim([-90, 90])
ax4.set_xlabel('CRU - Observed PET (mm/year)', fontsize = 8)
ax4.set_ylabel('Latitude (degrees)', fontsize = 8)
ax4.tick_params(labelsize = 8)        

ax5.scatter(pet_chinchu, pet_calc_T, color='green', marker='o', s = 5)
ax5.set_xlim([0, math.ceil(max(pet_chinchu) / 1000.) * 1000.])
ax5.set_ylim([0, math.ceil(max(pet_calc_T) / 1000.) * 1000.])
ax5.plot([0, 3000], [0, 3000], 'r--')
ax5.text(0.75, 0.1,'RMSE = ' + str(round(rms_cacl_T, 1)), ha='center', va='center', transform = ax5.transAxes)
ax5.set_xlabel('Observed PET (mm/year)', fontsize = 8)
ax5.set_ylabel('PET (mm/year) calculated via T(K)', fontsize = 8)
ax5.tick_params(labelsize = 8)        

ax6.scatter([x1 - x2 for (x1, x2) in zip(pet_calc_T, pet_chinchu)], lat_chinchu, color = 'red', marker='o', s = 5)
ax6.set_xlim([0, math.ceil(max([x1 - x2 for (x1, x2) in zip(pet_calc_T, pet_chinchu)]) / 1000.) * 1000.])
ax6.set_ylim([-90, 90])
ax6.set_xlabel('PET via T(K) - Observed PET (mm/year)', fontsize = 8)
ax6.set_ylabel('Latitude (degrees)', fontsize = 8)
ax6.tick_params(labelsize = 8)        

ax7.scatter([x1 - x2 for (x1, x2) in zip(pet_modis, pet_chinchu)], temp_chinchu, color = 'red', marker='o', s = 5)
ax7.set_xlim([0, math.ceil(max([x1 - x2 for (x1, x2) in zip(pet_modis, pet_chinchu)]) / 1000.) * 1000.])
ax7.set_ylim([-5, 50])
ax7.set_xlabel('MODIS - Observed PET (mm/year)', fontsize = 8)
ax7.set_ylabel('Mean annual Temperature (C)', fontsize = 8)
ax7.tick_params(labelsize = 8)        

ax8.scatter([x1 - x2 for (x1, x2) in zip(pet_cru, pet_chinchu)], temp_chinchu, color = 'red', marker='o', s = 5)
ax8.set_xlim([0, math.ceil(max([x1 - x2 for (x1, x2) in zip(pet_cru, pet_chinchu)]) / 1000.) * 1000.])
ax8.set_ylim([-5, 50])
ax8.set_xlabel('CRU - Observed PET (mm/year)', fontsize = 8)
ax8.set_ylabel('Mean annual Temperature (C)', fontsize = 8)
ax8.tick_params(labelsize = 8)        

ax9.scatter([x1 - x2 for (x1, x2) in zip(pet_calc_T, pet_chinchu)], temp_chinchu, color = 'red', marker='o', s = 5)
ax9.set_xlim([0, math.ceil(max([x1 - x2 for (x1, x2) in zip(pet_calc_T, pet_chinchu)]) / 1000.) * 1000.])
ax9.set_ylim([-5, 50])
ax9.set_xlabel('PET via T(K) - Observed PET (mm/year)', fontsize = 8)
ax9.set_ylabel('Mean annual Temperature (C)', fontsize = 8)
ax9.tick_params(labelsize = 8)        

ax10.scatter([x1 - x2 for (x1, x2) in zip(pet_modis, pet_chinchu)], gebco_elev, color = 'red', marker='o', s = 5)
ax10.set_xlim([0, math.ceil(max([x1 - x2 for (x1, x2) in zip(pet_modis, pet_chinchu)]) / 1000.) * 1000.])
ax10.set_ylim([0, 3.5])
ax10.set_xlabel('MODIS - Observed PET (mm/year)', fontsize = 8)
ax10.set_ylabel('GEBCO elevation (km asl.)', fontsize = 8)
ax10.tick_params(labelsize = 8)        

ax11.scatter([x1 - x2 for (x1, x2) in zip(pet_cru, pet_chinchu)], gebco_elev, color = 'red', marker='o', s = 5)
ax11.set_xlim([0, math.ceil(max([x1 - x2 for (x1, x2) in zip(pet_cru, pet_chinchu)]) / 1000.) * 1000.])
ax11.set_ylim([0, 3.5])
ax11.set_xlabel('CRU - Observed PET (mm/year)', fontsize = 8)
ax11.set_ylabel('GEBCO elevation (km asl.)', fontsize = 8)
ax11.tick_params(labelsize = 8)        

ax12.scatter([x1 - x2 for (x1, x2) in zip(pet_calc_T, pet_chinchu)], gebco_elev, color = 'red', marker='o', s = 5)
ax12.set_xlim([0, math.ceil(max([x1 - x2 for (x1, x2) in zip(pet_calc_T, pet_chinchu)]) / 1000.) * 1000.])
ax12.set_ylim([0, 3.5])
ax12.set_xlabel('PET via T(K) - Observed PET (mm/year)', fontsize = 8)
ax12.set_ylabel('GEBCO elevation (km asl.)', fontsize = 8)
ax12.tick_params(labelsize = 8)        

plt.savefig(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\PET_comparison.png', dpi = 300)
plt.close(fig)

"""
df['MODIS_min_obs'] = pd.Series([x1 - x2 for (x1, x2) in zip(pet_modis, pet_chinchu)], index=df.index)
df['CRU_min_obs'] = pd.Series([x1 - x2 for (x1, x2) in zip(pet_cru, pet_chinchu)], index=df.index)
df['Calc_PET_min_obs'] = pd.Series([x1 - x2 for (x1, x2) in zip(pet_calc_T, pet_chinchu)], index=df.index)
df.to_csv(csv_dir)
"""


##  Next step is to create paleo PET datasets based on temperature. Time scale is from LGM till 2000. 
##  There are multiple datasets available to estimate paleo recharge:
##      - we use the Worldclim paleo data (P and T) that give average annual values at LGM and mid-holocene.
##      - additionally, to estimate GW recharge at the intermediate time steps, based on temperature changes in the last ~30ka

lu_paleo_dir = lu_out_name
lu_2000_dir = os.path.join(out_dir, 'LU_coeff_' + str(2000) + '.tiff')
prec_lgm_dir = r'g:\_ORIGINAL_DATA\_WorldClim\ccsm4_LGM_mean_an_prec.tif'   # units are in mm
temp_lgm_dir = r'g:\_ORIGINAL_DATA\_WorldClim\ccsm4_LGM_mean_an_temp.tif'   # units are in K
#prec_mid_hol_dir = r'g:\_ORIGINAL_DATA\_WorldClim\ccmidbi12.tif'
#temp_mid_hol_dir = r'g:\_ORIGINAL_DATA\_WorldClim\ccmidbi1.tif'
temp_mid_hol_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\temp_midhol_resamp.tif'
prec_mid_hol_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\prec_midhol_resamp.tif'
#prec_present_dir = r'g:\_ORIGINAL_DATA\_WorldClim\wc2.0_bio_2.5m_12.tif'
prec_present_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\prec_present_resamp.tif'
temp_present_dir = r'g:\_ORIGINAL_DATA\_WorldClim\wc2.0_bio_2.5m_01.tif '   # units are in deg. C
ne_coastline = r'g:\_ORIGINAL_DATA\natural_earth\ne_10m_land'

#   calculate the current PET based on temperature, and then make the plot with MODIS, CRU and calculated PET
t_now_K_mean = np.mean(np.array(tmp_K_lst), axis = 0)
pet_cru_mean = np.mean(np.array(pet_lst) * 365.25, axis = 0)

pet_mean_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\pet_mean.nc'
xmin, ymin, xmax, ymax = -180., -90., 180., 90.
nrows, ncols = np.shape(pet_cru_mean)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(pet_mean_dir, pet_cru_mean.shape[1], pet_cru_mean.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(pet_mean[::-1])
output_raster.FlushCache()
output_raster = None

pet_calc_mean = np.zeros((360, 720))
for i in range(360):
    for j in range(720):
        pet_calc_mean[i, j] = 1.2e10 * math.exp(-4620 / t_now_K_mean[i, j])

pet_calc_mean_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\pet_calc_mean.nc'
xmin, ymin, xmax, ymax = -180., -90., 180., 90.
nrows, ncols = np.shape(pet_calc_mean)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(pet_calc_mean_dir, pet_calc_mean.shape[1], pet_calc_mean.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(pet_calc_mean[::-1])
output_raster.FlushCache()
output_raster = None

"""
#   read in the MODIS PET dataset for plotting
pet_modis_tif_dir = r'g:\_ORIGINAL_DATA\MODIS_evapotranspiration\MOD16A3_PET_2000_to_2013_mean.tif'
ds = gdal.Open(pet_modis_tif_dir)
ds_band = ds.GetRasterBand(1)
data = ds_band.ReadAsArray()
ds_gtrans = ds.GetGeoTransform()
"""

nc_out_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_SEAWAT_GW_RCH_input_data.nc'

temp_lgm = gdal.Open(temp_lgm_dir)
temp_lgm_band = temp_lgm.GetRasterBand(1)
temp_lgm_data = temp_lgm_band.ReadAsArray()
temp_midhol = gdal.Open(temp_mid_hol_dir)
temp_midhol_band = temp_midhol.GetRasterBand(1)
temp_midhol_data = temp_midhol_band.ReadAsArray()
temp_present = gdal.Open(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\temp_avg_annual_present.tif')
temp_present_band = temp_present.GetRasterBand(1)
temp_present_data = temp_present.ReadAsArray()

prec_lgm = gdal.Open(prec_lgm_dir)
prec_lgm_band = prec_lgm.GetRasterBand(1)
prec_lgm_data = prec_lgm_band.ReadAsArray()
prec_midhol = gdal.Open(prec_mid_hol_dir)
prec_midhol_band = prec_midhol.GetRasterBand(1)
prec_midhol_data = prec_midhol_band.ReadAsArray()
prec_present = gdal.Open(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\prec_present_resamp.tif')
prec_present_band = prec_present.GetRasterBand(1)
prec_present_data = prec_present.ReadAsArray()

#   read the annomaly (compared to long term averaged present values) for the last millenia
prec_LMR_dir = r'g:\_ORIGINAL_DATA\_Last_Millennium_Reanalysis_LMR\prate_MCruns_ensemble_mean_LMRv2.1.nc'
prec_LMR = xr.open_dataset(prec_LMR_dir, decode_times = False)
temp_LMR_dir = r'g:\_ORIGINAL_DATA\_Last_Millennium_Reanalysis_LMR\air_MCruns_ensemble_mean_LMRv2.1.nc' 
temp_LMR = xr.open_dataset(temp_LMR_dir, decode_times = False)

##  Create the directory name and the time steps that will be filled with different input data - P, PET, T and the final calculated GW_RCH
temp_nc_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_GW_rch_input_30kaBP_to_05kaAP.nc' 
time_steps = np.arange(-30000, -1000, 1000).tolist()
time_steps_2 = np.arange(-1900, 500, 100).tolist()
time_steps = time_steps + time_steps_2   #  the time step is the start of the SP (years BP ~ 2000CE)

x_coord_lst = np.linspace(-180., 180., 4320)
y_coord_lst = np.linspace(-60, 90, 1800)

nc_out_temp_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_TEMP_input_data_rounded_float32.nc'
nc_out_prec_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_PREC_input_data_rounded_float32.nc'
nc_out_pet_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_PET_input_data_rounded_float32.nc'
nc_out_lu_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_LU_input_data_rounded_float32.nc'


"""     --------------------------------------------------------------------    
                                TEMPERATURE
        --------------------------------------------------------------------    """

if not os.path.exists(nc_out_temp_dir):
    #   calculate the temperature for each time step until present - using the LGM, midhol and present T (in C)
    temp_arr_lst = []
    
    temp_lgm_data_nan = temp_lgm_data.astype('float64') 
    temp_lgm_data_nan[temp_lgm_data_nan == -32768] = np.nan
    
    temp_midhol_data_nan = temp_midhol_data
    temp_midhol_data_nan[temp_midhol_data_nan <= -32768] = np.nan
    
    #   the LGM occurred at approx -21ka BP so write constant temperature record until that time step 
    for a in range(10):
        temp_arr_lst.append(temp_lgm_data_nan / 10.)
    
    #   LGM and midhol coastal extent is different (LGM = lower sea level) and therefore the coverage of 
    #   the temperature datasets at these time steps is different. In cells where the the temperature is NaN 
    #   during the mid holocene change the temperatures graduall by an average temperature change in the surrounding area.
    temp_lgm_min_midhol = (temp_lgm_data_nan - temp_midhol_data_nan) 
    temp_lgm_min_midhol_per_ann = temp_lgm_min_midhol / 14.
    
    #   create boolean arrays to map areas where lgm is non nan and holocene is - these will be approximated in the next calculaitons
    lgm_nonnan = ~np.isnan(temp_lgm_data_nan)
    midhol_nonnan = np.isnan(temp_midhol_data_nan)
    cells_nonnan = np.logical_and(lgm_nonnan, midhol_nonnan)
    
    add_lgm_min_midhol = np.zeros((1800, 4320))
    for i in range(1800):
        for j in range(4320):
            if cells_nonnan[i, j] == True:
                avg_temp = temp_lgm_min_midhol_per_ann[max(0, i - 180) : min(1800, i + 180), max(0, j - 432) : min(4320, j + 432)]
                avg_temp_val = round(np.nanmean(avg_temp), 2)
                add_lgm_min_midhol[i, j] = avg_temp_val / 10.
    
    temp_lgm_min_midhol_per_ann_no_nan = np.nan_to_num(temp_lgm_min_midhol_per_ann, nan = 0)   
    lgm_min_midhol = temp_lgm_min_midhol_per_ann_no_nan / 10. + add_lgm_min_midhol
    
    for b in range(14):
        new_temp_arr = temp_arr_lst[-1] - lgm_min_midhol
        temp_arr_lst.append(new_temp_arr)
    
    #   do the same for the temperature between the midhol and 1ky BP = in total 5 time steps
    #   for the years 6000 BP and 2000 BP use the midhol data and the worldclim - LMR annomaly data respectively  
    temp_yr_MC_avg_0 = np.mean(temp_LMR['air']['time' == 0].values, axis = 0)
    temp_yr_MC_avg_0_resized = np.resize(temp_yr_MC_avg_0, temp_present_data.shape)
    
    temp_2000_BP = temp_present_data / 10. + temp_yr_MC_avg_0_resized
    temp_2000_BP_nan = temp_2000_BP
    temp_2000_BP_nan[temp_2000_BP_nan <= -990] = np.nan
    
    temp_midhol_min_2000_BP = (temp_midhol_data_nan / 10. - temp_2000_BP_nan) 
    temp_midhol_min_2000_BP_per_ann = temp_midhol_min_2000_BP / 4.
    
    for b in range(4):
        new_temp_arr = temp_arr_lst[-1] - temp_midhol_min_2000_BP_per_ann
        temp_arr_lst.append(new_temp_arr)
    
    #   final step is to create the temperature records for the centuries between 2ky BP and 0.5ky AP
    #   first we create an average centennial annomaly compared to the current temperature, do that for 20 time steps (until present)
    temp_present_data_nan = temp_present_data / 10
    temp_present_data_nan[temp_present_data_nan < -990] = np.nan
    
    for c in range(20):
        #   select the arrays from the original netcdf file
        temp_arr = temp_LMR['air'].sel(time = slice((365. * 100 * c), (365. * 100 * (c + 1)))).values
        temp_annomalies = np.mean(temp_arr, axis = 1)
        temp_annomalies_avg = np.mean(temp_annomalies, axis = 0)
        temp_annomalies_avg_resized = np.resize(temp_annomalies_avg, temp_present_data.shape)
        temp_new_arr = temp_present_data_nan + temp_annomalies_avg_resized
        temp_arr_lst.append(temp_new_arr)
    
    #   assume the same present day conditions for the final 5 centuries
    for d in range(5):    
        temp_arr_lst.append(temp_present_data_nan)
    
    #min_t, max_t = np.nanmin(temp_arr_lst), np.nanmax(temp_arr_lst)
    temp_zones = [-40., -10., -2.5, 0.0, 5.0, 10.0, 20., 32.]
    
    #   make plots for each time step
    temp_plots_out_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_temp_timesteps'
    
    for i in range(len(time_steps)):
        temp_arr = temp_arr_lst[i]
        time = time_steps[i]
        if time <= 0:
            out_name = os.path.join(temp_plots_out_dir, '_temp_deg_C_' + str(abs(time)) + '_BP.png')
        else:
            out_name = os.path.join(temp_plots_out_dir, '_temp_deg_C_' + str(abs(time)) + '_AP.png')
            
        print(out_name)
        ## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
        fig = plt.figure(figsize = (11.69, 8.27))
        #   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
        ax1 = plt.subplot2grid((4, 1), (0, 0))  #  LGM
        ax2 = plt.subplot2grid((4, 1), (3, 0))  #  Colorbar
        #   set the exact positions for each axis
        ax1.set_position([0.05, 0.05, 0.9, 0.75]) # [left, bottom, width, height]     
        ax2.set_position([0.25, 0.025, 0.5, 0.05]) # [left, bottom, width, height]
    
        #   the LGM
        m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
                    llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
        # draw parallels and meridians.
        m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
        m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
        data = temp_arr
        x = np.linspace(-180., 180., data.shape[1])
        y = np.linspace(-60, 90, data.shape[0])
        xx, yy = np.meshgrid(x, y)
        unique_n = temp_zones
        unique_n_label = [str(i) for i in list(unique_n)]
        #   define the colormap
        cmap = plt.cm.viridis_r
        norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
        m.contourf(xx, yy, data[::-1], levels = unique_n, cmap = cmap, norm = norm)
    
        cb = ColorbarBase(ax2, cmap = cmap, norm = norm, orientation = 'horizontal')
        cb.ax.set_xlabel('Average annual temperature (deg Celsius)', fontsize = 9, y = 1.025)
        cb.ax.tick_params(labelsize = 9)
        
        plt.savefig(out_name, dpi = 400, facecolor = 'w', edgecolor = 'w',
            orientation = 'landscape', papertype = None, format = None,
            transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
            frameon = None)
        plt.close()

    temp_arr_round_lst = []
    for i in range(len(temp_arr_lst)):
        new_arr = np.asarray(np.around(temp_arr_lst[i], 1), dtype = np.float32)
        temp_arr_round_lst.append(new_arr)
    
    xa_sum = xr.Dataset(data_vars = {'avg_annual_temp_deg_c' : (('time', 'y', 'x'), np.array(temp_arr_round_lst))},
                        coords = {'time' : time_steps,
                                  'x' : x_coord_lst,
                                  'y' : y_coord_lst})
    xa_sum.to_netcdf(nc_out_temp_dir)    
    del xa_sum
else:
    open_temp_nc = xr.open_dataset(nc_out_temp_dir) 
    temp_arr_lst = open_temp_nc['avg_annual_temp_deg_c'].values


"""     --------------------------------------------------------------------    
                                PRECIPITATION 
        --------------------------------------------------------------------    """

if not os.path.exists(nc_out_prec_dir):
    #   calculate the precipitation for each time step until present - using the LGM, midhol and present P (mm per year)
    prec_arr_lst = []
    prec_lgm_data_nan = prec_lgm_data.astype('float64') 
    prec_lgm_data_nan[prec_lgm_data_nan == -32768] = np.nan
    prec_midhol_data_nan = prec_midhol_data
    prec_midhol_data_nan[prec_midhol_data_nan <= -32768] = np.nan
    
    #   the LGM occurred at approx -21ka BP so write constant precipitation record until that time step 
    for a in range(10):
        prec_arr_lst.append(prec_lgm_data_nan)    
        
    prec_lgm_min_midhol = (prec_lgm_data_nan - prec_midhol_data_nan) 
    prec_lgm_min_midhol_per_ann = prec_lgm_min_midhol / 14.    
        
    #   create boolean arrays to map areas where lgm is non nan and holocene is - these will be approximated in the next calculaitons
    lgm_nonnan = ~np.isnan(prec_lgm_data_nan)
    midhol_nonnan = np.isnan(prec_midhol_data_nan)
    cells_nonnan = np.logical_and(lgm_nonnan, midhol_nonnan)    
        
    add_lgm_min_midhol = np.zeros((1800, 4320))
    for i in range(1800):
        for j in range(4320):
            if cells_nonnan[i, j] == True:
                avg_prec = prec_lgm_min_midhol_per_ann[max(0, i - 180) : min(1800, i + 180), max(0, j - 432) : min(4320, j + 432)]
                avg_prec_val = round(np.nanmean(avg_prec), 2)
                add_lgm_min_midhol[i, j] = avg_prec_val 
    
    lgm_min_midhol = prec_lgm_min_midhol_per_ann + add_lgm_min_midhol

    prec_lgm_min_midhol_per_ann_no_nan = np.nan_to_num(prec_lgm_min_midhol_per_ann, nan = 0)   
    lgm_min_midhol = prec_lgm_min_midhol_per_ann_no_nan / 10. + add_lgm_min_midhol
    
    for b in range(14):
        new_prec_arr = prec_arr_lst[-1] - lgm_min_midhol
        prec_arr_lst.append(new_prec_arr)
        
    #   do the same for the temperature between the midhol and 1ky BP = in total 5 time steps
    #   for the years 6000 BP and 2000 BP use the midhol data and the worldclim - LMR annomaly data respectively  
    prec_yr_MC_avg_0 = np.mean(prec_LMR['prate']['time' == 0].values, axis = 0)
    prec_yr_MC_avg_0_mm_yr = prec_yr_MC_avg_0 * 86400 * 365.
    prec_yr_MC_avg_0_resized = np.resize(prec_yr_MC_avg_0_mm_yr, prec_present_data.shape)
    
    prec_2000_BP = prec_present_data + prec_yr_MC_avg_0_resized
    prec_2000_BP_nan = prec_2000_BP
    prec_2000_BP_nan[prec_2000_BP_nan <= -990] = np.nan
    
    prec_midhol_min_2000_BP = (prec_midhol_data_nan - prec_2000_BP_nan) 
    prec_midhol_min_2000_BP_per_ann = prec_midhol_min_2000_BP / 4.
    
    for b in range(4):
        new_prec_arr = prec_arr_lst[-1] - prec_midhol_min_2000_BP_per_ann
        prec_arr_lst.append(new_prec_arr)
    
    #   final step is to create the temperature records for the centuries between 2ky BP and 0.5ky AP
    #   first we create an average centennial annomaly compared to the current temperature, do that for 20 time steps (until present)
    prec_present_data_nan = prec_present_data
    prec_present_data_nan[prec_present_data_nan < -990] = np.nan
    
    for c in range(20):
        #   select the arrays from the original netcdf file
        prec_arr = temp_LMR['air'].sel(time = slice((365. * 100 * c), (365. * 100 * (c + 1)))).values
        prec_annomalies = np.mean(prec_arr, axis = 1)
        prec_annomalies_avg = np.mean(prec_annomalies, axis = 0)
        prec_annomalies_avg_resized = np.resize(prec_annomalies_avg, prec_present_data.shape)
        prec_new_arr = prec_present_data_nan + prec_annomalies_avg_resized
        prec_arr_lst.append(prec_new_arr)
    
    #   assume the same present day conditions for the final 5 centuries
    for d in range(5):    
        prec_arr_lst.append(prec_present_data_nan)
    
    min_prec, max_prec = np.nanmin(prec_arr_lst), np.nanmax(prec_arr_lst)
    prec_zones = [0, 100, 250, 500, 1000, 2500, 5000, 15000]
    
    #   make plots for each time step
    prec_plots_out_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_prec_timesteps'
    
    for i in range(len(time_steps)):
        prec_arr = prec_arr_lst[i]
        time = time_steps[i]
        if time <= 0:
            out_name = os.path.join(prec_plots_out_dir, '_prec_mm_yr_' + str(abs(time)) + '_BP.png')
        else:
            out_name = os.path.join(prec_plots_out_dir, '_prec_mm_yr_' + str(abs(time)) + '_AP.png')
            
        print(out_name)
        ## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
        fig = plt.figure(figsize = (11.69, 8.27))
        #   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
        ax1 = plt.subplot2grid((4, 1), (0, 0))  #  LGM
        ax2 = plt.subplot2grid((4, 1), (3, 0))  #  Colorbar
        #   set the exact positions for each axis
        ax1.set_position([0.05, 0.05, 0.9, 0.75]) # [left, bottom, width, height]     
        ax2.set_position([0.25, 0.025, 0.5, 0.05]) # [left, bottom, width, height]
    
        #   the LGM
        m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
                    llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
        # draw parallels and meridians.
        m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
        m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
        data = prec_arr
        x = np.linspace(-180., 180., data.shape[1])
        y = np.linspace(-60, 90, data.shape[0])
        xx, yy = np.meshgrid(x, y)
        unique_n = [0, 100, 250, 500, 1000, 2500, 5000, 15000]#prec_zones
        unique_n_label = [str(i) for i in list(unique_n)]
        #   define the colormap
        cmap = plt.cm.plasma_r
        norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
        m.contourf(xx, yy, data[::-1], levels = unique_n, interpolation = None, cmap = cmap, norm = norm)
    
        cb = ColorbarBase(ax2, cmap = cmap, norm = norm, orientation = 'horizontal')
        cb.ax.set_xlabel('Average annual precipitation (mm/yr)', fontsize = 9, y = 1.025)
        cb.ax.tick_params(labelsize = 9)
        
        plt.savefig(out_name, dpi = 400, facecolor = 'w', edgecolor = 'w',
            orientation = 'landscape', papertype = None, format = None,
            transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
            frameon = None)
        plt.close()

    prec_arr_round_lst = []
    for i in range(len(prec_arr_lst)):
        new_arr = np.asarray(np.around(prec_arr_lst[i], 1), dtype = np.float32)
        prec_arr_round_lst.append(new_arr)
    
    xa_sum = xr.Dataset(data_vars = {'avg_annual_prec_mm_yr' : (('time', 'y', 'x'), np.array(prec_arr_round_lst))},
                        coords = {'time' : time_steps,
                                  'x' : x_coord_lst,
                                  'y' : y_coord_lst})
    xa_sum.to_netcdf(nc_out_prec_dir)    
else:
    open_prec_nc = xr.open_dataset(nc_out_prec_dir) 
    prec_arr_lst = open_prec_nc['avg_annual_prec_mm_yr'].values

temp_pres_interp_dir_out = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\temp_pres_interpolated.tiff'
if not os.path.exists(temp_pres_interp_dir_out):
    temp_present_arr = temp_present_data / 10
    #temp_present_arr[temp_present_arr == -999.9] = np.nan
    buff_row, buff_col = 180, 432 
    for i in range(1800):
        print(i)
        for j in range(4320):
                if temp_present_arr[i, j] != -999.9:
                    pass
                else: 
                    if temp_lgm_data[i, j] != -32768:
                        try:
                            avg_temp_non_nan = temp_present_arr[max(0, i - buff_row) : min(1800, i + buff_row), max(0, j - buff_col) : min(4320, j + buff_col)]
                            bol_arr = avg_temp_non_nan > -990.
                            avg_temp = avg_temp_non_nan[bol_arr]
                            avg_temp_val = round(np.mean(avg_temp), 1)
                            temp_present_arr[i, j] = avg_temp_val
                        except ValueError:
                            temp_present_arr[i, j] = np.nan
                    else:
                        temp_present_arr[i, j] = np.nan
    xmin, ymin, xmax, ymax = -180., -90., 180., 90.
    nrows, ncols = np.shape(gw_rch_avg)
    xres = (xmax - xmin)/float(ncols)
    yres = (ymax - ymin)/float(nrows)
    geotransform = (xmin, xres, 0, ymax, 0, -yres)  
    output_raster = gdal.GetDriverByName('GTiff').Create(temp_pres_interp_dir_out, temp_present_arr.shape[1], temp_present_arr.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(4326)      
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(temp_present_arr)
    output_raster.FlushCache()
    output_raster = None
else:
    temp_to_open = gdal.Open(temp_pres_interp_dir_out)
    temp_to_open_band = temp_to_open.GetRasterBand(1)
    temp_present_arr = temp_to_open_band.ReadAsArray()


"""     --------------------------------------------------------------------    
                                PET 
        --------------------------------------------------------------------    """

#   do the same about the PET calculations
pet_calc_mean_resamp = gdal.Open(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\pet_calc_mean_resamp.tif')
pet_calc_mean_resamp_band = pet_calc_mean_resamp.GetRasterBand(1)
pet_calc_mean_resamp_data = pet_calc_mean_resamp_band.ReadAsArray()
pet_calc_mean_resamp_data = pet_calc_mean_resamp_data[: 1800, :]

pet_calc_mean_resamp_data_interp = pet_calc_mean_resamp_data
pet_dir_out = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\PET_pres_interpolated.tiff'
if not os.path.exists(pet_dir_out):
    
    buff_row, buff_col = 180, 432 
    for i in range(1800):
        print(i)
        for j in range(4320):
                if not math.isnan(pet_calc_mean_resamp_data_interp[i, j]):
                    pass
                else: 
                    if temp_lgm_data[i, j] != -32768:
                        try:
                            avg_pet_non_nan = pet_calc_mean_resamp_data[max(0, i - buff_row) : min(1800, i + buff_row), max(0, j - buff_col) : min(4320, j + buff_col)]
                            avg_pet_val = round(np.nanmean(avg_pet_non_nan), 1)
                            pet_calc_mean_resamp_data_interp[i, j] = avg_pet_val
                        except ValueError:
                            pet_calc_mean_resamp_data_interp[i, j] = np.nan
                    else:
                        pet_calc_mean_resamp_data_interp[i, j] = np.nan
                        
    xmin, ymin, xmax, ymax = -180., -90., 180., 90.
    nrows, ncols = np.shape(gw_rch_avg)
    xres = (xmax - xmin)/float(ncols)
    yres = (ymax - ymin)/float(nrows)
    geotransform = (xmin, xres, 0, ymax, 0, -yres)  
    output_raster = gdal.GetDriverByName('GTiff').Create(pet_dir_out, pet_calc_mean_resamp_data_interp.shape[1], pet_calc_mean_resamp_data_interp.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(4326)      
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(pet_calc_mean_resamp_data_interp)
    output_raster.FlushCache()
    output_raster = None
else:
    pet_to_open = gdal.Open(pet_dir_out)
    pet_to_open_band = pet_to_open.GetRasterBand(1)
    pet_present_arr = pet_to_open_band.ReadAsArray()


    
if not os.path.exists(nc_out_pet_dir):    
    #   make a list of arrays, with calculated PET for each time step
    pet_arr_lst = []
    
    for k in range(len(temp_arr_lst)):
        
        temp_arr_ts = temp_arr_lst[k]
        pet_arr = np.zeros((1800, 4320))
        
        for i in range(1800):
            print(i)
            for j in range(4320):
                
                if math.isnan(temp_arr_ts[i, j]):
                    pet_arr[i, j] = np.nan
                else:
                    pet_arr[i, j] = pet_present_arr[i, j] * math.exp(-4620 / (273.15 + temp_arr_ts[i, j])) / math.exp(-4620 / (273.15 + temp_present_arr[i, j])) 
    
        pet_arr_lst.append(pet_arr)
        
    #   make plots for each time step
    prec_plots_out_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_pet_timesteps'
        
    for i in range(len(time_steps)):
        pet_arr = pet_arr_lst[i]
        time = time_steps[i]
        if time <= 0:
            out_name = os.path.join(prec_plots_out_dir, '_PET_mm_yr_' + str(abs(time)) + '_BP.png')
        else:
            out_name = os.path.join(prec_plots_out_dir, '_PET_mm_yr_' + str(abs(time)) + '_AP.png')
            
        print(out_name)
        ## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
        fig = plt.figure(figsize = (11.69, 8.27))
        #   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
        ax1 = plt.subplot2grid((4, 1), (0, 0))  #  LGM
        ax2 = plt.subplot2grid((4, 1), (3, 0))  #  Colorbar
        #   set the exact positions for each axis
        ax1.set_position([0.05, 0.05, 0.9, 0.75]) # [left, bottom, width, height]     
        ax2.set_position([0.25, 0.025, 0.5, 0.05]) # [left, bottom, width, height]
    
        #   the LGM
        m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
                    llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
        # draw parallels and meridians.
        m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
        m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
        data = pet_arr
        x = np.linspace(-180., 180., data.shape[1])
        y = np.linspace(-60, 90, data.shape[0])
        xx, yy = np.meshgrid(x, y)
        unique_n = [0, 100, 250, 500, 1000, 2500, 5000]#prec_zones
        unique_n_label = [str(i) for i in list(unique_n)]
        #   define the colormap
        cmap = plt.cm.plasma_r
        norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
        m.contourf(xx, yy, data[::-1], levels = unique_n, interpolation = None, cmap = cmap, norm = norm)
    
        cb = ColorbarBase(ax2, cmap = cmap, norm = norm, orientation = 'horizontal')
        cb.ax.set_xlabel('Average annual PET (mm/yr)', fontsize = 9, y = 1.025)
        cb.ax.tick_params(labelsize = 9)
        
        plt.savefig(out_name, dpi = 400, facecolor = 'w', edgecolor = 'w',
            orientation = 'landscape', papertype = None, format = None,
            transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
            frameon = None)
        plt.close()       

    pet_arr_round_lst = []
    for i in range(len(pet_arr_lst)):
        new_arr = np.asarray(np.around(pet_arr_lst[i], 1), dtype = np.float32)
        pet_arr_round_lst.append(new_arr)
    
    xa_sum = xr.Dataset(data_vars = {'avg_pet_mm_yr' : (('time', 'y', 'x'), np.array(pet_arr_round_lst))},
                        coords = {'time' : time_steps,
                                  'x' : x_coord_lst,
                                  'y' : y_coord_lst})
    xa_sum.to_netcdf(nc_out_pet_dir)    
    
else:
    open_pet_nc = xr.open_dataset(nc_out_pet_dir) 
    pet_arr_lst = open_pet_nc['avg_pet_mm_yr'].values    
    
    
"""     --------------------------------------------------------------------    
                                LAND USE 
        --------------------------------------------------------------------    """
        
#   make a list of arrays with the LU for each time step. Assume the same LU for time steps before 12000 BP
#   using the HYDE 3.2.1. database and its classification - lumping multiple classes into one so it matches
#   the classes used by Chinchu et. al. 

##  11, 12 - Urban 
##  21, 22, 23, 24 - originally classified as Village, here put into the Urban category
##  31, 32, 33, 34 - Cropland
##  41, 42, 43 - Pasture
##  51, 52, 53, 54, 61 - Forest
##  62, 63 - barren     

hyde_LU_dir = r'g:\_ORIGINAL_DATA\_HYDE_3_2_1\anthromes\_unzipped'    

#   first deal with the LGM LU estimation - which again has to be interpolated to cover the uncovered areas during the LGM low sea level.
lu_ts_file = os.path.join(hyde_LU_dir, 'anthromes10000BC.asc')
lu_to_open = gdal.Open(lu_ts_file)
lu_to_open_band = lu_to_open.GetRasterBand(1)
lu_arr = lu_to_open_band.ReadAsArray()[:1800, :]

lu_arr_nan = lu_arr
lu_arr_nan[lu_arr_nan == -9999.] = np.nan

prec_lgm_data_nan = prec_lgm_data.astype('float64') 
prec_lgm_data_nan[prec_lgm_data_nan == -32768] = np.nan
prec_midhol_data_nan = prec_midhol_data
prec_midhol_data_nan[prec_midhol_data_nan <= -32768] = np.nan

prec_lgm_min_midhol = (prec_lgm_data_nan - prec_midhol_data_nan) 
prec_lgm_min_midhol_per_ann = prec_lgm_min_midhol / 14.    
    
#   create boolean arrays to map areas where lgm is non nan and holocene is - these will be approximated in the next calculaitons
lgm_nonnan = ~np.isnan(prec_lgm_data_nan)
midhol_nonnan = np.isnan(prec_midhol_data_nan)
cells_nonnan = np.logical_and(lgm_nonnan, midhol_nonnan)  

lu_lgm_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\LU_hyde_LGM_interpolated.tiff'
if not os.path.exists(lu_lgm_dir):
    buff_row, buff_col = 180, 432 
    lu_lgm = lu_arr
    for i in range(1800):
        print(i)
        for j in range(4320):
            if cells_nonnan[i, j] == True:
                
                most_occ_arr = lu_arr_nan[max(0, i - 180) : min(1800, i + 180), max(0, j - 432) : min(4320, j + 432)]
                int_arr = np.asarray(most_occ_arr, np.int64)
                
                #   the count there will also have zero
                unique_elements, counts_elements = np.unique(int_arr, return_counts=True)
                unique_elements = unique_elements[1:].tolist()
                counts_elements = counts_elements[1:].tolist()
                lu_val = unique_elements[counts_elements.index(max(counts_elements))]
                
                lu_lgm[i, j] = lu_val 
                        
    xmin, ymin, xmax, ymax = -180., -90., 180., 90.
    nrows, ncols = np.shape(lu_lgm)
    xres = (xmax - xmin)/float(ncols)
    yres = (ymax - ymin)/float(nrows)
    geotransform = (xmin, xres, 0, ymax, 0, -yres)  
    output_raster = gdal.GetDriverByName('GTiff').Create(lu_lgm_dir, lu_lgm.shape[1], lu_lgm.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(4326)      
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(lu_lgm)
    output_raster.FlushCache()
    output_raster = None
else:
    lu_lgm_to_open = gdal.Open(lu_lgm_dir)
    lu_lgm_to_open_band = lu_lgm_to_open.GetRasterBand(1)
    lu_lgm = lu_lgm_to_open_band.ReadAsArray()

##  11, 12 - Urban 
##  21, 22, 23, 24 - originally classified as Village, here put into the Urban category
##  31, 32, 33, 34 - Cropland
##  41, 42, 43, 54, 62 - Pasture
##  51, 52, 53, 61 - Forest
##  63 - barren     

if not os.path.exists(nc_out_lu_dir):
    #   make a list of arrays, with calculated PET for each time step
    lu_arr_lst = []
    
    lu_lgm_reclass = lu_lgm
    lu_lgm_reclass[(lu_lgm_reclass <= 24)] = 2
    lu_lgm_reclass[(lu_lgm_reclass >= 31) & (lu_lgm_reclass <= 34)] = 5
    lu_lgm_reclass[(lu_lgm_reclass >= 41) & (lu_lgm_reclass <= 43)] = 4
    lu_lgm_reclass[(lu_lgm_reclass == 54)] = 4
    lu_lgm_reclass[(lu_lgm_reclass == 62)] = 4
    lu_lgm_reclass[(lu_lgm_reclass >= 51) & (lu_lgm_reclass <= 53)] = 3
    lu_lgm_reclass[(lu_lgm_reclass == 61)] = 3
    lu_lgm_reclass[(lu_lgm_reclass == 63)] = 1
    lu_lgm_reclass[(lu_lgm_reclass == 70)] = 1
    
    #   assign the LGM (HYDE at 10000BC to the stress periods)
    for i in range(19):
        lu_arr_lst.append(lu_lgm_reclass)
    
    #   next deal with time steps until 0AD, that are by millenia. 
    for i in range(9):
        #   first deal with the LGM LU estimation - which again has to be interpolated to cover the uncovered areas during the LGM low sea level.
        file_name = 'anthromes0' + str(9000 - i * 1000) + 'BC.asc'
        print(file_name)
        lu_ts_file = os.path.join(hyde_LU_dir, 'anthromes10000BC.asc')
        lu_to_open = gdal.Open(lu_ts_file)
        lu_to_open_band = lu_to_open.GetRasterBand(1)
        lu_arr = lu_to_open_band.ReadAsArray()[:1800, :]
        lu_arr_nan = lu_arr
        lu_arr_nan[lu_arr_nan == -9999.] = np.nan
    
        lu_arr_2 = lu_arr
        for i in range(1800):
            print(i)
            for j in range(4320):
                if cells_nonnan[i, j] == True:
                    
                    most_occ_arr = lu_arr_nan[max(0, i - 180) : min(1800, i + 180), max(0, j - 432) : min(4320, j + 432)]
                    int_arr = np.asarray(most_occ_arr, np.int64)
                    
                    #   the count there will also have zero
                    unique_elements, counts_elements = np.unique(int_arr, return_counts=True)
                    unique_elements = unique_elements[1:].tolist()
                    counts_elements = counts_elements[1:].tolist()
                    lu_val = unique_elements[counts_elements.index(max(counts_elements))]
                    
                    lu_arr_2[i, j] = lu_val 
    
        lu_arr_reclass = lu_arr_2
        lu_arr_reclass[(lu_arr_reclass <= 24)] = 2
        lu_arr_reclass[(lu_arr_reclass >= 31) & (lu_arr_reclass <= 34)] = 5
        lu_arr_reclass[(lu_arr_reclass >= 41) & (lu_arr_reclass <= 43)] = 4
        lu_arr_reclass[(lu_arr_reclass == 54)] = 4
        lu_arr_reclass[(lu_arr_reclass == 62)] = 4
        lu_arr_reclass[(lu_arr_reclass >= 51) & (lu_arr_reclass <= 53)] = 3
        lu_arr_reclass[(lu_arr_reclass == 61)] = 3
        lu_arr_reclass[(lu_arr_reclass == 63)] = 1
        lu_arr_reclass[(lu_arr_reclass == 70)] = 1
    
        lu_arr_lst.append(lu_arr_reclass)
    
    #   deal with 0AD
    lu_ts_file = os.path.join(hyde_LU_dir, 'anthromes0AD.asc')
    lu_to_open = gdal.Open(lu_ts_file)
    lu_to_open_band = lu_to_open.GetRasterBand(1)
    lu_arr = lu_to_open_band.ReadAsArray()[:1800, :]  
        
    lu_arr_reclass = lu_arr
    lu_arr_reclass[(lu_arr_reclass <= 24)] = 2
    lu_arr_reclass[(lu_arr_reclass >= 31) & (lu_arr_reclass <= 34)] = 5
    lu_arr_reclass[(lu_arr_reclass >= 41) & (lu_arr_reclass <= 43)] = 4
    lu_arr_reclass[(lu_arr_reclass == 54)] = 4
    lu_arr_reclass[(lu_arr_reclass == 62)] = 4
    lu_arr_reclass[(lu_arr_reclass >= 51) & (lu_arr_reclass <= 53)] = 3
    lu_arr_reclass[(lu_arr_reclass == 61)] = 3
    lu_arr_reclass[(lu_arr_reclass == 63)] = 1
    lu_arr_reclass[(lu_arr_reclass == 70)] = 1
    
    lu_arr_lst.append(lu_arr_reclass)
    
    for j in range(1, 20):
        #   first deal with the LGM LU estimation - which again has to be interpolated to cover the uncovered areas during the LGM low sea level.
        file_name = 'anthromes' + str(j * 100) + 'AD.asc'
        print(file_name)
    
        lu_ts_file = os.path.join(hyde_LU_dir, file_name)
        lu_to_open = gdal.Open(lu_ts_file)
        lu_to_open_band = lu_to_open.GetRasterBand(1)
        lu_arr = lu_to_open_band.ReadAsArray()[:1800, :]  
            
        lu_arr_reclass = lu_arr
        lu_arr_reclass[(lu_arr_reclass <= 24)] = 2
        lu_arr_reclass[(lu_arr_reclass >= 31) & (lu_arr_reclass <= 34)] = 5
        lu_arr_reclass[(lu_arr_reclass >= 41) & (lu_arr_reclass <= 43)] = 4
        lu_arr_reclass[(lu_arr_reclass == 54)] = 4
        lu_arr_reclass[(lu_arr_reclass == 62)] = 4
        lu_arr_reclass[(lu_arr_reclass >= 51) & (lu_arr_reclass <= 53)] = 3
        lu_arr_reclass[(lu_arr_reclass == 61)] = 3
        lu_arr_reclass[(lu_arr_reclass == 63)] = 1
        lu_arr_reclass[(lu_arr_reclass == 70)] = 1
        
        lu_arr_lst.append(lu_arr_reclass)
    
    for j in range(5):
        #   first deal with the LGM LU estimation - which again has to be interpolated to cover the uncovered areas during the LGM low sea level.
        file_name = 'anthromes2017AD.asc'
        print(file_name)
    
        lu_ts_file = os.path.join(hyde_LU_dir, file_name)
        lu_to_open = gdal.Open(lu_ts_file)
        lu_to_open_band = lu_to_open.GetRasterBand(1)
        lu_arr = lu_to_open_band.ReadAsArray()[:1800, :]  
            
        lu_arr_reclass = lu_arr
        lu_arr_reclass[(lu_arr_reclass <= 24)] = 2
        lu_arr_reclass[(lu_arr_reclass >= 31) & (lu_arr_reclass <= 34)] = 5
        lu_arr_reclass[(lu_arr_reclass >= 41) & (lu_arr_reclass <= 43)] = 4
        lu_arr_reclass[(lu_arr_reclass == 54)] = 4
        lu_arr_reclass[(lu_arr_reclass == 62)] = 4
        lu_arr_reclass[(lu_arr_reclass >= 51) & (lu_arr_reclass <= 53)] = 3
        lu_arr_reclass[(lu_arr_reclass == 61)] = 3
        lu_arr_reclass[(lu_arr_reclass == 63)] = 1
        lu_arr_reclass[(lu_arr_reclass == 70)] = 1
        
        lu_arr_lst.append(lu_arr_reclass)
    
    xa_sum = xr.Dataset(data_vars = {'LU_class' : (('time', 'y', 'x'), np.array(lu_arr_lst))},
                        coords = {'time' : time_steps,
                                  'x' : x_coord_lst,
                                  'y' : y_coord_lst})
    xa_sum.to_netcdf(nc_out_lu_dir)    
    
else:
    open_lu_nc = xr.open_dataset(nc_out_lu_dir) 
    lu_arr_lst = open_lu_nc['LU_class'].values  



#   calculate the groundwater recharge for each time step using the equation from Chinchu
cell_size = xr.open_dataset(r'g:\_ORIGINAL_DATA\_PCR_GLOB_waterdemand\cellsize05min_correct.nc')['cellsize05min_correct_map'].values[:-360, :]
cell_size_km2 = cell_size / 1000000
        

#   read and transform the Ksoil raster
ksoil_dir = r'g:\_ORIGINAL_DATA\Hydraul_Param_SoilGrids_Schaap_0\Ks_cm_d.tif'
ksoil_to_open = gdal.Open(ksoil_dir)
ksoil_to_open_band = ksoil_to_open.GetRasterBand(1)
ksoil_arr = ksoil_to_open_band.ReadAsArray()[:600, :]  
ksoil_arr_resize = np.zeros((1800, 4320))
for a in range(ksoil_arr.shape[0]):
    for b in range(ksoil_arr.shape[1]):
        if ksoil_arr[a, b] < -1000:
            ksoil_arr_resize[a * 3 : a * 3 + 3, b * 3 : b * 3 + 3] = 0
        else:
            ksoil_arr_resize[a * 3 : a * 3 + 3, b * 3 : b * 3 + 3] = ksoil_arr[a, b]

xmin, ymin, xmax, ymax = -180., -60., 180., 90.
nrows, ncols = np.shape(ksoil_arr_resize)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(r'g:\_ORIGINAL_DATA\Hydraul_Param_SoilGrids_Schaap_0\Ks_cm_d_1km.tif', ksoil_arr_resize.shape[1], ksoil_arr_resize.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(ksoil_arr_resize)
output_raster.FlushCache()
output_raster = None   


#   read in the clay content dataset
clay_cont_dir = r'G:\_ORIGINAL_DATA\soilgrids\CLYPPT_M_sl6_250m_ll_resamp_10km.tif'
clay_cont = gdal.Open(clay_cont_dir)
clay_cont_band = clay_cont.GetRasterBand(1)
clay_cont_arr = clay_cont_band.ReadAsArray()
#   calculate the number of empty rows in the Northern latitudes
# 72 rows missing in the north
for i in range(72):
    clay_cont_arr = np.insert(clay_cont_arr, [0], [4320 * [np.nan]], axis = 0)
# 48 rows missing in the south 
for i in range(48):
    clay_cont_arr = np.append(clay_cont_arr, [4320 * [np.nan]], axis = 0)

where_are_NaNs = np.isnan(clay_cont_arr)
clay_cont_arr[where_are_NaNs] = 0

gw_rch_plot_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH_Ksoil\_GW_RCH_timesteps'
gw_rch_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH_Ksoil\_tif_GW_RCH_timesteps'

ts_lst = open_lu_nc['time'].values.tolist()
gw_rch_arr_lst = []

gw_ts_max, gw_ts_min, gw_ts_mean = [], [], []
pet_ts_max, pet_ts_min, pet_ts_mean = [], [], []
p_ts_max, p_ts_min, p_ts_mean = [], [], []

rch_total_km3_yr = []
rch_avg_mm_yr = []

for ts in ts_lst:
    if ts < 0:
        suff = '_BP'
    elif ts > 0:
        suff = '_AP'
        
    out_name = 'gw_rch_ts_' + str(abs(ts)) + suff
    print(out_name)
    
    p_val = open_prec_nc.sel(time = ts)['avg_annual_prec_mm_yr'].values
    pet_val = open_pet_nc.sel(time = ts)['avg_pet_mm_yr'].values
    lu_val = open_lu_nc.sel(time = ts)['LU_class'].values
    
    #rch_val = (5.3543 + 0.0081 * p_val - 0.0043 * pet_val + 0.9567 * lu_val ) ** 2 # 
    #rch_val_old = (5.3543 + 0.0081 * p_val - 0.0043 * pet_val + 0.9567 * lu_val ) ** 2 # 
    rch_val = (6.3781 + 0.0086 * p_val - 0.0044 * pet_val + 1.0335 * lu_val - clay_cont_arr * 0.0606) ** 2 
    
    #rch_val_old_km = rch_val_old / 1000000
    rch_val_km = rch_val / 1000000
    
    rch_km3_yr = np.multiply(rch_val_km, cell_size_km2)    
    
    rch_avg_mm_yr.append(np.nanmean(rch_val))
    rch_total_km3_yr.append(np.nansum(rch_km3_yr))
    
    """
    gw_rch_arr_lst.append(rch_val)
    
    gw_ts_max.append(np.nanmax(rch_val))
    gw_ts_min.append(np.nanmin(rch_val))
    gw_ts_mean.append(np.nanmean(rch_val))

    pet_ts_max.append(np.nanmax(pet_val))
    pet_ts_min.append(np.nanmin(pet_val))
    pet_ts_mean.append(np.nanmean(pet_val))

    p_ts_max.append(np.nanmax(p_val))
    p_ts_min.append(np.nanmin(p_val))
    p_ts_mean.append(np.nanmean(p_val))
    
    rch_km3_yr = np.multiply(rch_val_km, cell_size_km2)    
    """
    plt_dir = os.path.join(gw_rch_plot_dir, out_name + '.png')
    """
    tif_dir = os.path.join(gw_rch_dir, out_name + '.tif')
    
    xmin, ymin, xmax, ymax = -180., -60., 180., 90.
    nrows, ncols = np.shape(rch_val)
    xres = (xmax - xmin)/float(ncols)
    yres = (ymax - ymin)/float(nrows)
    geotransform = (xmin, xres, 0, ymax, 0, -yres)  
    output_raster = gdal.GetDriverByName('GTiff').Create(tif_dir, rch_val.shape[1], rch_val.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(4326)      
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(rch_val)
    output_raster.FlushCache()
    output_raster = None   
    """
    ## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
    fig = plt.figure(figsize = (11.69, 7.27))
    #   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
    ax1 = plt.subplot2grid((2, 1), (0, 0))  #  LGM
    ax2 = plt.subplot2grid((2, 1), (1, 0))  #  Colorbar
    #   set the exact positions for each axis
    ax1.set_position([0.05, 0.05, 0.9, 0.75]) # [left, bottom, width, height]     
    ax2.set_position([0.4, 0.005, 0.2, 0.025]) # [left, bottom, width, height]
    
    #   the LGM
    m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
                llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
    # draw parallels and meridians.
    m.drawparallels(np.arange(-45.,76.,15.), labels = [True, False, False, False], dashes = [2,2])
    m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, True, False, False], dashes = [2,2])
    data = rch_val
    x = np.linspace(-180., 180., data.shape[1])
    y = np.linspace(-60, 90, data.shape[0])
    xx, yy = np.meshgrid(x, y)
    unique_n = [0, 100, 250, 500, 1000, 2500, 5000]#prec_zones
    unique_n_label = [str(i) for i in list(unique_n)]
    #   define the colormap
    cmap = plt.cm.plasma_r
    norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
    m.contourf(xx, yy, data[::-1], levels = unique_n, interpolation = None, cmap = cmap, norm = norm)

    cb = ColorbarBase(ax2, cmap = cmap, norm = norm, orientation = 'horizontal')
    cb.ax.set_xlabel('Average annual GW recharge (mm/yr)', fontsize = 9, y = 1.025)
    cb.ax.tick_params(labelsize = 8)
    
    plt.savefig(plt_dir, dpi = 300, facecolor = 'w', edgecolor = 'w',
        orientation = 'landscape', format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()           

nc_out_gw_rch_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH_Ksoil\_GW_RCH_rounded_float32.nc'
xa_sum = xr.Dataset(data_vars = {'gw_rch_mm_yr' : (('time', 'y', 'x'), np.array(gw_rch_arr_lst))},
                    coords = {'time' : time_steps,
                              'x' : x_coord_lst,
                              'y' : y_coord_lst})
xa_sum.to_netcdf(nc_out_gw_rch_dir)    



"""
rch_total_km3_yr = [i / 1000 for i in rch_total_km3_yr]
ts_lst = [i / 1000 for i in ts_lst]

fig, ax1 = plt.subplots()

ax1.set_xlabel('time (Ka BP)', fontsize = 8)
ax1.set_ylabel('Total GW recharge (1000 * km3/yr)', labelpad = 10, fontsize = 8)
ax1.plot(ts_lst, rch_total_km3_yr, color = 'green', label = 'total GW recharge', linewidth = 1)
ax1.tick_params(axis='y', labelsize = 8)
ax1.set_xlim(-20, 0.4)
ax1.tick_params(axis='x', labelsize = 8)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Average GW recharge (mm/yr)', rotation = 270, labelpad = 20, fontsize = 8)  # we already handled the x-label with ax1
ax2.plot(ts_lst, rch_avg_mm_yr, color=color, label = 'average GW recharge', linewidth = 1)
ax2.tick_params(axis='y', labelsize = 8)


# ask matplotlib for the plotted objects and their labels
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=3, prop={"size":8})

x_major_locator = matplotlib.ticker.FixedLocator([-17.5, -15., -12.5, -10., -7.5, -5, -2.5, 0], nbins = None)   
ax1.xaxis.set_major_locator(x_major_locator)
ax1.grid(axis='x', color = 'grey', alpha = 0.5, linewidth = 0.25)

#fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.savefig(r'g:\Water_Nexus\_A4_paper\_figures\gw_rch.png', dpi = 400, facecolor = 'w', edgecolor = 'w',
    orientation = 'landscape', papertype = None, format = None,
    transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
    frameon = None)
plt.close()    

import pylab
pylab.figure()
pylab.hist(rch_val[~np.isnan(rch_val)])
pylab.show()


"""






gw_rch_plot_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH_diff\_GW_RCH_timesteps'
gw_rch_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH_diff\_tif_GW_RCH_timesteps'

for ts in ts_lst:
    if ts < 0:
        suff = '_BP'
    elif ts > 0:
        suff = '_AP'
        
    out_name = 'gw_rch_ts_DIFF_' + str(abs(ts)) + suff
    print(out_name)
    
    p_val = open_prec_nc.sel(time = ts)['avg_annual_prec_mm_yr'].values
    pet_val = open_pet_nc.sel(time = ts)['avg_pet_mm_yr'].values
    lu_val = open_lu_nc.sel(time = ts)['LU_class'].values
    
    #rch_val = (5.3543 + 0.0081 * p_val - 0.0043 * pet_val + 0.9567 * lu_val ) ** 2 # 
    rch_val_old = (5.3543 + 0.0081 * p_val - 0.0043 * pet_val + 0.9567 * lu_val ) ** 2 # 
    rch_val = (6.3781 + 0.0086 * p_val - 0.0044 * pet_val + 1.0335 * lu_val - clay_cont_arr * 0.0606) ** 2 
    
    plt_dir = os.path.join(gw_rch_plot_dir, out_name + '.png')
    tif_dir = os.path.join(gw_rch_dir, out_name + '.tif')
    
    xmin, ymin, xmax, ymax = -180., -60., 180., 90.
    nrows, ncols = np.shape(rch_val)
    xres = (xmax - xmin)/float(ncols)
    yres = (ymax - ymin)/float(nrows)
    geotransform = (xmin, xres, 0, ymax, 0, -yres)  
    output_raster = gdal.GetDriverByName('GTiff').Create(tif_dir, rch_val.shape[1], rch_val.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(4326)      
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(rch_val_old - rch_val)
    output_raster.FlushCache()
    output_raster = None   
    
    ## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
    fig = plt.figure(figsize = (11.69, 8.27))
    #   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
    ax1 = plt.subplot2grid((2, 1), (0, 0))  #  LGM
    ax2 = plt.subplot2grid((2, 1), (1, 0))  #  Colorbar
    #   set the exact positions for each axis
    ax1.set_position([0.05, 0.05, 0.9, 0.75]) # [left, bottom, width, height]     
    ax2.set_position([0.25, 0.025, 0.5, 0.025]) # [left, bottom, width, height]
    
    #   the LGM
    m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
                llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
    # draw parallels and meridians.
    m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
    m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
    data = rch_val_old - rch_val
    x = np.linspace(-180., 180., data.shape[1])
    y = np.linspace(-60, 90, data.shape[0])
    xx, yy = np.meshgrid(x, y)
    unique_n = [-1000, -100, -50, -25, -10, 0, 10, 25, 50, 100, 1000]#prec_zones
    unique_n_label = [str(i) for i in list(unique_n)]
    #   define the colormap
    cmap = plt.cm.plasma_r
    norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
    m.contourf(xx, yy, data[::-1], levels = unique_n, interpolation = None, cmap = cmap, norm = norm)

    cb = ColorbarBase(ax2, cmap = cmap, norm = norm, orientation = 'horizontal')
    cb.ax.set_xlabel('Average annual GW recharge (mm/yr)', fontsize = 9, y = 1.025)
    cb.ax.tick_params(labelsize = 9)
    
    plt.savefig(plt_dir, dpi = 400, facecolor = 'w', edgecolor = 'w',
        orientation = 'landscape', papertype = None, format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
        frameon = None)
    plt.close()    







#else:
#    open_pet_nc = xr.open_dataset(nc_out_pet_dir) 
#    pet_arr_lst = open_pet_nc['avg_pet_mm_yr'].values   



crop = u1_crop.sel(time = 0)['prop_crop'].values
primary_land = u1_primary_land.sel(time = 0)['prop_primary'].values
pasture = u1_pasture.sel(time = 0)['prop_pasture'].values
secondary_land = u1_secondary_land.sel(time = 0)['prop_secd'].values

#   calculate the LU coeffecient array and save it as tiff file
lu_coef_arr = 5 * crop + 4 * pasture + 3 * (secondary_land + primary_land) + 2 * urban + 1 * (1 - crop - pasture - secondary_land - primary_land - urban) 
lu_out_name = os.path.join(out_dir, 'LU_coeff_' + str(1400) + '.tiff')


















    
    
    
    
    
    
pcr_tot_abs = xr.open_dataset(r'g:\_ORIGINAL_DATA\_PCR_GLOB_waterdemand\totalAbstraction_annuaTot_output_2015-12-31_to_2015-12-31.nc')    
pcr_tot_abs_arr = pcr_tot_abs['total_abstraction'].values     
    
    
plt.imshow(pcr_tot_abs_arr[0, :, :])
    









pet_lgm_arr = np.zeros((1800, 4320))
buff_row, buff_col = 180, 432 

for i in range(1800):
    for j in range(4320):
        #if math.isnan(pet_calc_mean_resamp_data[i, j]) or temp_lgm_data[i, j] == -32768:
        if temp_lgm_data[i, j] == -32768:
            pet_lgm_arr[i, j] = np.nan
        else:
            if temp_present_data[i, j] != -9999:
                pet_lgm_arr[i, j] = pet_calc_mean_resamp_data[i, j] * (math.exp(-4620 / (273.15 + temp_lgm_data[i, j] / 10)) / math.exp(-4620 / (273.15 + temp_present_data[i, j] / 10)))  
            else: 
                try:
                    avg_temp_non_nan = temp_present_data[max(0, i - buff_row) : min(1800, i + buff_row), max(0, j - buff_col) : min(4320, j + buff_col)]
                    bol_arr = avg_temp_non_nan > -1000
                    avg_temp = avg_temp_non_nan[bol_arr]
                    avg_temp_val = int(round(np.mean(avg_temp), 0))
                    
                    pet_non_nan = pet_calc_mean_resamp_data[max(0, i - buff_row) : min(1800, i + buff_row), max(0, j - buff_col) : min(4320, j + buff_col)]
                    bol_arr = pet_non_nan > -1000
                    avg_pet = pet_non_nan[bol_arr]
                    avg_pet_val = int(round(np.mean(avg_pet), 0))
                              
                    pet_lgm_arr[i, j] = avg_pet_val * (math.exp(round(-4620 / (273.15 + (temp_lgm_data[i, j] / 10)), 2)) / math.exp(round(-4620 / (273.15 + (avg_temp_val / 10)), 2)))

                    #print(i, j, pet_lgm_arr[i, j], pet_calc_mean_resamp_data[i, j],273.15 + (temp_lgm_data[i, j] / 10), 273.15 + (avg_temp_val / 10))

                except ValueError:
                    #print(i, j, '5000', pet_calc_mean_resamp_data[i, j], temp_lgm_data[i, j] / 10, temp_present_data[i, j] / 10)
                    pet_lgm_arr[i, j] = np.nan


pet_midhol_arr = np.zeros((1800, 4320))
for i in range(1800):
    for j in range(4320):
        #if math.isnan(pet_calc_mean_resamp_data[i, j]) or temp_lgm_data[i, j] == -32768:
        if temp_midhol_data[i, j] < -32768:
            pet_midhol_arr[i, j] = np.nan
        else:
            if temp_present_data[i, j] != -9999:
                pet_midhol_arr[i, j] = pet_calc_mean_resamp_data[i, j] * (math.exp(-4620 / (273.15 + temp_midhol_data[i, j] / 10)) / math.exp(-4620 / (273.15 + temp_present_data[i, j] / 10)))  
            else: 
                try:
                    avg_temp_non_nan = temp_present_data[max(0, i - buff_row) : min(1800, i + buff_row), max(0, j - buff_col) : min(4320, j + buff_col)]
                    bol_arr = avg_temp_non_nan > -1000
                    avg_temp = avg_temp_non_nan[bol_arr]
                    avg_temp_val = int(round(np.mean(avg_temp), 0))
                    
                    pet_non_nan = pet_calc_mean_resamp_data[max(0, i - buff_row) : min(1800, i + buff_row), max(0, j - buff_col) : min(4320, j + buff_col)]
                    bol_arr = pet_non_nan > -1000
                    avg_pet = pet_non_nan[bol_arr]
                    avg_pet_val = int(round(np.mean(avg_pet), 0))
                              
                    pet_midhol_arr[i, j] = avg_pet_val * (math.exp(round(-4620 / (273.15 + (temp_midhol_data[i, j] / 10)), 2)) / math.exp(round(-4620 / (273.15 + (avg_temp_val / 10)), 2)))

                    #print(i, j, pet_lgm_arr[i, j], pet_calc_mean_resamp_data[i, j],273.15 + (temp_lgm_data[i, j] / 10), 273.15 + (avg_temp_val / 10))

                except ValueError:
                    #print(i, j, '5000', pet_calc_mean_resamp_data[i, j], temp_lgm_data[i, j] / 10, temp_present_data[i, j] / 10)
                    pet_midhol_arr[i, j] = np.nan















plt.imshow(prec_lgm_data_nan)






    
pet_arr_lst = 53 * [np.zeros((1800, 4320))]
lu_arr_lst = 53 * [np.zeros((1800, 4320))]
gw_rch_arr_lst = 53 * [np.zeros((1800, 4320))]
    





xa_sum = xr.Dataset(data_vars = {'avg_annual_temp_deg_c' : (('time', 'y', 'x'), np.array(temp_arr_lst)),
                                 'avg_annual_prec_mm_yr' : (('time', 'y', 'x'), np.array(prec_arr_lst))},
                    coords = {'time' : time_steps,
                              'x' : x_coord_lst,
                              'y' : y_coord_lst})
xa_sum.to_netcdf(nc_out_dir)    



xa_sum = xr.Dataset(data_vars = {'avg_annual_temp_deg_c' : (('time', 'y', 'x'), np.array(temp_arr_lst)),
                                 'avg_annual_prec_mm_yr' : (('time', 'y', 'x'), np.array(prec_arr_lst)),
                                 'PET_calculated' : (('time', 'y', 'x'), np.array(pet_arr_lst)),
                                 'land_use' : (('time', 'y', 'x'), np.array(lu_arr_lst)),
                                 'GW_RCH_mm_day' : (('time', 'y', 'x'), np.array(gw_rch_arr_lst))},
                    coords = {'time' : time_steps,
                              'x' : x_coord_lst,
                              'y' : y_coord_lst})
xa_sum.to_netcdf(nc_out_dir)           

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



    
    
phyda_nc = xr.open_dataset(r'g:\_ORIGINAL_DATA\_PHYDA\da_hydro_AprMar_r.1-2000_d.05-Jan-2018.nc', decode_times = False)


phyda_nc['tas_sg'].values.shape



plt.imshow(add_lgm_min_midhol)










            """
            try:
                pet_lgm_arr[i, j] = pet_calc_mean_resamp_data[i, j] * (math.exp(-4620 / (273.15 + temp_lgm_data[i, j] / 10)) / math.exp(-4620 / (273.15 + temp_present_data[i, j] / 10)))             
            except TypeError:
                avg_temp_non_nan = temp_present_data[max(0, i - 50) : min(1800, i + 50), max(0, j - 50) : min(4320, j + 50)]
                bol_arr = avg_temp_non_nan > -1000
                avg_temp = avg_temp_non_nan[bol_arr]
                avg_temp_val = int(round(np.mean(avg_temp), 0))
                
                pet_non_nan = pet_calc_mean_resamp_data[max(0, i - 50) : min(1800, i + 50), max(0, j - 50) : min(4320, j + 50)]
                bol_arr = pet_non_nan > -1000
                avg_pet = pet_non_nan[bol_arr]
                avg_pet_val = int(round(np.mean(avg_pet), 0))
                                
                pet_lgm_arr[i, j] = 5000
                #PET_LGM[i, j] = avg_pet_val * (math.exp(-4620 / (273.15 + (temp_lgm_data[i, j] / 10))) / (math.exp(-4620 / (273.15 + (avg_temp_val / 10)))))
            """





## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
fig = plt.figure(figsize = (8.27, 11.69))
#   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
ax1 = plt.subplot2grid((4, 1), (0, 0))  #  LGM
ax2 = plt.subplot2grid((4, 1), (1, 0))  #  Mid_holocene
ax3 = plt.subplot2grid((4, 1), (2, 0))  #  Present
ax4 = plt.subplot2grid((4, 1), (3, 0))  #  Colorbar
#   set the exact positions for each axis
ax1.set_position([0.05, 0.65, 0.9, 0.23]) # [left, bottom, width, height]
ax2.set_position([0.05, 0.4, 0.9, 0.23])
ax3.set_position([0.05, 0.15, 0.9, 0.23])        
ax4.set_position([0.25, 0.08, 0.5, 0.02]) # [left, bottom, width, height]

#   the LGM
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
data = pet_lgm_arr
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-60, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
unique_n = np.arange(0, 3500, 500)
unique_n_label = [str(i) for i in list(unique_n)]
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
m.contourf(xx, yy, data[::-1], levels = unique_n, cmap = cmap)

#   the mid_holocene
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax2)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
data = pet_midhol_arr
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-60, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
unique_n = np.arange(0, 3500, 500)
unique_n_label = [str(i) for i in list(unique_n)]
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
m.contourf(xx, yy, data[::-1], levels = unique_n, cmap = cmap)

#   the current condition
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax3)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
data = pet_calc_mean_resamp_data#[:1800, :]
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-60, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
unique_n = np.arange(0, 3500, 500)
unique_n_label = [str(i) for i in list(unique_n)]
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 

data = data.astype('float64') 
#data[data == -1.70000000e+308] = np.nan
m.contourf(xx, yy, data[::-1], levels = unique_n, cmap = cmap)#, norm = norm)

cb = ColorbarBase(ax4, cmap = cmap, norm = norm, orientation = 'horizontal')
cb.ax.set_xlabel('Average PET (mm/year)', fontsize = 11, y  = 1.025)
cb.ax.tick_params(labelsize = 9)

plt.savefig(os.path.join(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH', 'PET_calculated_LGM_midhol_current.png'), dpi = 400, facecolor = 'w', edgecolor = 'w',
    orientation = 'portrait', papertype = None, format = None,
    transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
    frameon = None)
plt.close()
















## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
fig = plt.figure(figsize = (8.27, 11.69))
#   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
ax1 = plt.subplot2grid((4, 1), (0, 0))  #  LGM
ax2 = plt.subplot2grid((4, 1), (1, 0))  #  Mid_holocene
ax3 = plt.subplot2grid((4, 1), (2, 0))  #  Present
ax4 = plt.subplot2grid((4, 1), (3, 0))  #  Colorbar
#   set the exact positions for each axis
ax1.set_position([0.05, 0.65, 0.9, 0.23]) # [left, bottom, width, height]
ax2.set_position([0.05, 0.4, 0.9, 0.23])
ax3.set_position([0.05, 0.15, 0.9, 0.23])        
ax4.set_position([0.25, 0.08, 0.5, 0.02]) # [left, bottom, width, height]

#   the LGM
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds = gdal.Open(prec_lgm_dir)
data = ds.ReadAsArray()
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-60, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
unique_n = [0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2., 2.5, 5., 10., 20.]
unique_n_label = ['0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.5', '2', '2.5', '5', '10', '20']
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
data = data.astype('float64') 
data[data == -32768] = np.nan
data_mm_d = data/ 365.25
m.contourf(xx, yy, data_mm_d[::-1], levels = unique_n, cmap = cmap)

#   the mid_holocene
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax2)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds = gdal.Open(prec_mid_hol_dir)
data = ds.ReadAsArray()
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-60, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
#unique_n = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2., 2.5, 3., 3.5, 4., 4.5, 5., 10., 20.]
#unique_n_label = ['0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1.25', '1.5', '1.75', '2', '2.5', '3', '3.5', '4', '4.5', '5', '10', '20']
unique_n = [0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2., 2.5, 5., 10., 20.]
unique_n_label = ['0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.5', '2', '2.5', '5', '10', '20']
#   define the colormap
#cmap = plt.cm.viridis_r
#norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
data = data.astype('float64') 
data[data == -32768] = np.nan
data_mm_d = data/ 365.25
m.contourf(xx, yy, data_mm_d[::-1], levels = unique_n, cmap = cmap)#, norm = norm)

#   the current condition
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax3)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds = gdal.Open(prec_present_dir)
data = ds.ReadAsArray()
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-90, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
unique_n = [0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2., 2.5, 5., 10., 20.]
unique_n_label = ['0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.5', '2', '2.5', '5', '10', '20']
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 

data = data.astype('float64') 
data[data == -1.70000000e+308] = np.nan
data_mm_d = data/ 365.25
m.contourf(xx, yy, data_mm_d[::-1], levels = unique_n, cmap = cmap)#, norm = norm)

cb = ColorbarBase(ax4, cmap = cmap, norm = norm, orientation = 'horizontal')
cb.ax.set_xlabel('Average daily precipitation (mm)', fontsize = 11, y  = 1.025)
cb.ax.tick_params(labelsize = 9)

plt.savefig(os.path.join(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH', 'worldclim_prec_comparison.png'), dpi = 400, facecolor = 'w', edgecolor = 'w',
    orientation = 'portrait', papertype = None, format = None,
    transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
    frameon = None)
plt.close()





## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
fig = plt.figure(figsize = (8.27, 11.69))
#   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
ax1 = plt.subplot2grid((4, 1), (0, 0))  #  Mid_holocene - LGM
ax2 = plt.subplot2grid((4, 1), (1, 0))  #  Present - Mid_holocene
ax3 = plt.subplot2grid((4, 1), (3, 0))  #  Colorbar
#   set the exact positions for each axis
ax1.set_position([0.05, 0.575, 0.9, 0.375]) # [left, bottom, width, height]
ax2.set_position([0.05, 0.15, 0.9, 0.375]) 
ax3.set_position([0.25, 0.08, 0.5, 0.02]) # [left, bottom, width, height]

#   the LGM
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds_lgm = gdal.Open(prec_lgm_dir)
ds_midhol = gdal.Open(prec_mid_hol_dir)
data_lgm = ds_lgm.ReadAsArray()
data_midhol = ds_midhol.ReadAsArray()
x = np.linspace(-180., 180., data_lgm.shape[1])
y = np.linspace(-60, 90, data_lgm.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
#unique_n = np.arange(-20, 40, 5).tolist()#[-10., -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10]
unique_n = [-20, -10, 0, 2.5, 5, 10, 35]
unique_n_label = [str(i) for i in unique_n]
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
data_lgm = data_lgm.astype('float64') 
data_lgm[data_lgm == -32768] = np.nan
data_midhol = data_midhol.astype('float64') 
data_midhol[data_midhol < -10000] = np.nan

data_mm_d_diff = (data_midhol / 365.25) - (data_lgm / 365.25)
m.contourf(xx, yy, data_mm_d_diff[::-1], levels = unique_n, cmap = cmap)


#   the mid_holocene
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax2)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds_pres = gdal.Open(prec_present_dir)
data_pres = ds_pres.ReadAsArray()
x = np.linspace(-180., 180., data_pres.shape[1])
y = np.linspace(-60, 90, data_pres.shape[0])
xx, yy = np.meshgrid(x, y)
unique_n = [-20, -10, 0, 2.5, 5, 10, 35]
unique_n_label = [str(i) for i in unique_n]
#unique_n_label = ['0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.5', '2', '2.5', '5', '10', '20']
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
data_pres = data_pres.astype('float64') 
data_pres[data_pres < -10000] = np.nan

data_mm_d_diff = (data_pres / 365.25) - (data_midhol / 365.25)
m.contourf(xx, yy, data_mm_d_diff[::-1], levels = unique_n, cmap = cmap)#, norm = norm)

cb = ColorbarBase(ax3, cmap = cmap, norm = norm, orientation = 'horizontal')
cb.ax.set_xlabel('Average daily precipitation (mm)', fontsize = 11, y  = 1.025)
cb.ax.tick_params(labelsize = 9)

plt.savefig(os.path.join(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH', 'worldclim_prec_differences_timesteps.png'), dpi = 400, facecolor = 'w', edgecolor = 'w',
    orientation = 'portrait', papertype = None, format = None,
    transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
    frameon = None)
plt.close()



np.nanmin(data_mm_d_diff), np.nanmax(data_mm_d_diff)




























## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
fig = plt.figure(figsize = (8.27, 11.69))
#   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
ax1 = plt.subplot2grid((4, 1), (0, 0))  #  LGM
ax2 = plt.subplot2grid((4, 1), (1, 0))  #  Mid_holocene
ax3 = plt.subplot2grid((4, 1), (2, 0))  #  Present
ax4 = plt.subplot2grid((4, 1), (3, 0))  #  Colorbar
#   set the exact positions for each axis
ax1.set_position([0.05, 0.65, 0.9, 0.23]) # [left, bottom, width, height]
ax2.set_position([0.05, 0.4, 0.9, 0.23])
ax3.set_position([0.05, 0.15, 0.9, 0.23])        
ax4.set_position([0.25, 0.08, 0.5, 0.02]) # [left, bottom, width, height]

#   the LGM
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds = gdal.Open(temp_lgm_dir)
data = ds.ReadAsArray()
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-60, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
unique_n = np.arange(-50, 55, 5)
unique_n_label = [str(i) for i in list(unique_n)]
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
data = data.astype('float64') 
data[data == -32768] = np.nan
m.contourf(xx, yy, data[::-1] / 10., levels = unique_n, cmap = cmap)

#   the mid_holocene
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax2)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds = gdal.Open(temp_mid_hol_dir)
data = ds.ReadAsArray()
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-60, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
#unique_n = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2., 2.5, 3., 3.5, 4., 4.5, 5., 10., 20.]
#unique_n_label = ['0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1.25', '1.5', '1.75', '2', '2.5', '3', '3.5', '4', '4.5', '5', '10', '20']
unique_n = np.arange(-50, 55, 5)
unique_n_label = [str(i) for i in list(unique_n)]
#   define the colormap
#cmap = plt.cm.viridis_r
#norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
data = data.astype('float64') 
data[data == -32768] = np.nan
m.contourf(xx, yy, data[::-1] / 10., levels = unique_n, cmap = cmap)#, norm = norm)

#   the current condition
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax3)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds = gdal.Open(temp_present_dir)
data = ds.ReadAsArray()
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-90, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
unique_n = np.arange(-50, 55, 5)
unique_n_label = [str(i) for i in list(unique_n)]
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 

data = data.astype('float64') 
data[data == -1.70000000e+308] = np.nan
m.contourf(xx, yy, data[::-1], levels = unique_n, cmap = cmap)#, norm = norm)

cb = ColorbarBase(ax4, cmap = cmap, norm = norm, orientation = 'horizontal')
cb.ax.set_xlabel('Average annual temprature (C)', fontsize = 11, y  = 1.025)
cb.ax.tick_params(labelsize = 9)

plt.savefig(os.path.join(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH', 'worldclim_temperature_comparison.png'), dpi = 400, facecolor = 'w', edgecolor = 'w',
    orientation = 'portrait', papertype = None, format = None,
    transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
    frameon = None)
plt.close()































## first plot the P, T and then PET at the 3 timesteps based on the Worldclim datasets - LGM (-21ka), Mid-Holocene (-6ka), Present (0ka)
fig = plt.figure(figsize = (8.27, 11.69))
#   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
ax1 = plt.subplot2grid((4, 1), (0, 0))  #  LGM
ax2 = plt.subplot2grid((4, 1), (1, 0))  #  Mid_holocene
ax3 = plt.subplot2grid((4, 1), (2, 0))  #  Present
ax4 = plt.subplot2grid((4, 1), (3, 0))  #  Colorbar
#   set the exact positions for each axis
ax1.set_position([0.05, 0.65, 0.9, 0.23]) # [left, bottom, width, height]
ax2.set_position([0.05, 0.4, 0.9, 0.23])
ax3.set_position([0.05, 0.15, 0.9, 0.23])        
ax4.set_position([0.25, 0.08, 0.5, 0.02]) # [left, bottom, width, height]

#   the LGM
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax1)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds = gdal.Open(temp_lgm_dir)
data = ds.ReadAsArray()
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-60, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
unique_n = np.arange(-50, 55, 5)
unique_n_label = [str(i) for i in list(unique_n)]
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
data = data.astype('float64') 
data[data == -32768] = np.nan
m.contourf(xx, yy, data[::-1] / 10., levels = unique_n, cmap = cmap)

#   the mid_holocene
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax2)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds = gdal.Open(temp_mid_hol_dir)
data = ds.ReadAsArray()
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-60, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
#unique_n = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2., 2.5, 3., 3.5, 4., 4.5, 5., 10., 20.]
#unique_n_label = ['0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1.25', '1.5', '1.75', '2', '2.5', '3', '3.5', '4', '4.5', '5', '10', '20']
unique_n = np.arange(-50, 55, 5)
unique_n_label = [str(i) for i in list(unique_n)]
#   define the colormap
#cmap = plt.cm.viridis_r
#norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
data = data.astype('float64') 
data[data == -32768] = np.nan
m.contourf(xx, yy, data[::-1] / 10., levels = unique_n, cmap = cmap)#, norm = norm)

#   the current condition
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 90,\
            llcrnrlon = -180, urcrnrlon =180, lat_ts = 20, resolution = 'i', ax = ax3)
# draw parallels and meridians.
m.drawparallels(np.arange(-45.,76.,15.), labels = [True, True, False, False], dashes = [2,2])
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], dashes = [2,2])
#m.readshapefile(ne_coastline, 'areas', linewidth = 0.5)
ds = gdal.Open(temp_present_dir)
data = ds.ReadAsArray()
x = np.linspace(-180., 180., data.shape[1])
y = np.linspace(-90, 90, data.shape[0])
xx, yy = np.meshgrid(x, y)
#unique_n = [0., 100., 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 10000]
#unique_n_label = ['0.', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '10000']
unique_n = np.arange(-50, 55, 5)
unique_n_label = [str(i) for i in list(unique_n)]
#   define the colormap
cmap = plt.cm.viridis_r
norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 

data = data.astype('float64') 
data[data == -1.70000000e+308] = np.nan
m.contourf(xx, yy, data[::-1], levels = unique_n, cmap = cmap)#, norm = norm)

cb = ColorbarBase(ax4, cmap = cmap, norm = norm, orientation = 'horizontal')
cb.ax.set_xlabel('Average annual temprature (C)', fontsize = 11, y  = 1.025)
cb.ax.tick_params(labelsize = 9)

plt.savefig(os.path.join(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH', 'current_PET_comparison.png'), dpi = 400, facecolor = 'w', edgecolor = 'w',
    orientation = 'portrait', papertype = None, format = None,
    transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
    frameon = None)
plt.close()











































#   read in the csv file with P, T, and PET values
df = pd.read_csv(csv_in_dir)
sample = df.sample(572)   #   572 is 80% of the data points
test_sample = df.loc[~df.index.isin(sample.index)]

#X = sample[['T', 'abs_lat']]
X = sample[['T']]
y = sample['PET']

# sklearn
regr = linear_model.LinearRegression()
regr.fit(X, y)

print('Intercept: \n', regr.intercept_)
print('Coefficients: \n', regr.coef_)

predicted_pet = []
predicted_t = []

for index, row in test_sample.iterrows():
    #predicted_pet.append(regr.predict([row[['T']]])[0])
    predicted_pet.append(regr.intercept_ + row['T'] * regr.coef_[0])
    predicted_t.append(row['T'])
    #print(row['T'], row['PET'], regr.predict([row[['T', 'abs_lat']]])[0] - row['PET'])


# Note the difference in argument order
X = sm.add_constant(X) # adding a constant
model = sm.OLS(y, X).fit()
predictions = model.predict(X) # make the predictions by the model

# Print out the statistics
model.summary()



# plot the data
plt.plot(predicted_t, predicted_pet)








