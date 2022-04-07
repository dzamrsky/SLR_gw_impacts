# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 11:25:59 2021

@author: daniel
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:14:55 2019

@author: daniel
"""

"""
This script creates ARP for each sub-region. It takes the sub-region shapefile and loops through it. 
The only additional dataset that is read in is the new MERIT DEM, whose values are read in every 25, 50 and 100 meters 
to allow for various model resolutions. Both the value in the cross-section point and average cross-section point value
are measured . This is done by taking a point value each 25, 50 or 100m for a distance of 1km in line perpendicular to 
the cross-section 

"""

import numpy as np
import geopandas as gpd
import os
import gdal
import sys
import utm
import xarray as xr
import pandas as pd
import math
import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.ioff()
import osr
from shapely.geometry import Point, MultiPoint
from shapely.ops import nearest_points
##  directory paths on my laptop


sys_dir = r'g:\Water_Nexus\_A4_GUM\_scripts\_fc_scripts'
input_files_dir = r'g:\Water_Nexus\_A4_GUM\_input_data_COSCAT_REG'
main_dir = r'g:\Water_Nexus\_A4_GUM\_MODEL_input_files'

#   define all necessary directories and database table names
offshore_sed_dir = r'g:\_ORIGINAL_DATA\seabed_lithology\seabed_lithology_v1.nc'
wtd_dir = r'g:\_ORIGINAL_DATA\water_table_depth\wtd_world_3.tif'       
soil_type_dir = r'g:\_CREATED_DATA\soilgrids\merge_all_LZW.tif'
soil_thk_dir = r'g:\_ORIGINAL_DATA\soilgrids\absolute_depth_to_bedrock_BDTICM_M_1km_ll.tif'
glhymps_top_dir = r'g:\_CREATED_DATA\_A2_data\GLHYMPS\GLHYMPS_2_upper_layer_NODATA.tif'
glhymps_bot_dir = r'g:\_CREATED_DATA\_A2_data\GLHYMPS\GLHYMPS_1_bottom_layer_compressed.tif'
k_soil_dir = r'g:\_ORIGINAL_DATA\Hydraul_Param_SoilGrids_Schaap_0\ks_100cm_mm_d.tiff'
drn_rate_dir = r'g:\_ORIGINAL_DATA\LCS_D_global_drainage\LCS_Dd_global.tif'
gw_rch_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH\_tif_GW_RCH_timesteps'
gw_rch_clay_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_GW_RCH_Ksoil\_tif_GW_RCH_timesteps'
gebco_2014_dir = r'g:\_ORIGINAL_DATA\GEBCO 2014\GEBCO_2014_2D.nc'
merit_dem_folder_dir = r'g:\_ORIGINAL_DATA\MERIT_DEM'
coastal_dem_folder_dir = r'g:\_ORIGINAL_DATA\CoastalDEMv1.1_MSL\CoastalDEM90v1.1'
"""

sys_dir = r'/home/dzamrsky/_A4/_fc_scripts'
input_files_dir = r'/projects/0/qt16165/_dzamrsky/_A4_input_data'
main_dir = r'/projects/0/qt16165/_dzamrsky/_A4_MODEL_input_files'

offshore_sed_dir = os.path.join(input_files_dir, 'seabed_lithology', 'seabed_lithology_v1.nc') 
wtd_dir = os.path.join(input_files_dir, 'water_table_depth', 'wtd_world_3.tif')     
soil_type_dir = os.path.join(input_files_dir, 'soilgrids', 'merge_all_LZW.tif')   
soil_thk_dir = os.path.join(input_files_dir, 'soilgrids', 'absolute_depth_to_bedrock_BDTICM_M_1km_ll.tif')   
glhymps_top_dir = os.path.join(input_files_dir, 'GLHYMPS', 'GLHYMPS_2_upper_layer_NODATA.tif')  
glhymps_bot_dir = os.path.join(input_files_dir, 'GLHYMPS', 'GLHYMPS_1_bottom_layer_compressed.tif')  
k_soil_dir = os.path.join(input_files_dir, 'Hydraul_Param_SoilGrids_Schaap_0', 'ks_100cm_mm_d.tiff')  
drn_rate_dir = os.path.join(input_files_dir, 'LCS_D_global_drainage', 'LCS_Dd_global.tif')
gw_rch_dir = os.path.join(input_files_dir, '_GW_RCH', '_tif_GW_RCH_timesteps')  
gw_rch_clay_dir = os.path.join(input_files_dir, '_GW_RCH_Ksoil', '_tif_GW_RCH_timesteps')  
gebco_2014_dir = os.path.join(input_files_dir, 'GEBCO_2014', 'GEBCO_2014_2D.nc')  
merit_dem_folder_dir = os.path.join(input_files_dir, 'MERIT_DEM')  
coastal_dem_folder_dir = os.path.join(input_files_dir, 'CoastalDEMv1.1_MSL', 'CoastalDEM90v1.1') 
coastal_dem_folder_corr_res_dir = os.path.join(input_files_dir, 'CoastalDEMv1.1_MSL', 'corr_res')
"""
sys.path.append(sys_dir)    
import DTBSE_functions_py3 as db_fc
import _create_ARP_RIVBAS_tools as rivbas
import ws_masterscript_py3 as ws_ms
#import ws_datatools_py3 as ws               # connect to dbase, run sqls etc.
#from modeltools_py3 import cs_model 

#   the id loop is the only input into the script, the rest will be found in the lookup table
#id_loop = int(sys.argv[1])
id_loop = 16
print(id_loop)

#   open the csv file and read the row with the id_loop
#   change the status of the model to -1, the values are as follows
#       0  - not finished nor started
#       -1 - not finished, but started
#       1  - finished
#df_in = pd.read_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb.csv'))
df_in = pd.read_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb_ADDED_SRMs.csv'))
id_loop_row = df_in.loc[df_in['id_loop'] == id_loop]
status = df_in.loc[df_in['id_loop'] == id_loop, 'finished'].values[0]
if status != 1:#== 0:
    df_in.loc[df_in['id_loop'] == id_loop, 'finished'] = -1
    df_in.update(df_in)
    #df_in.to_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb.csv'), sep = ',', encoding = 'utf-8', index = False)
    df_in.to_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb_ADDED_SRMs.csv'), sep = ',', encoding = 'utf-8', index = False)
del df_in    

#   read in the original raw shapefile with coastal points and coastal types
#cs_regs_shp_dir = os.path.join(input_files_dir, 'cs_reg_id.shp')
cs_regs_shp_dir = os.path.join(input_files_dir, 'cs_reg_id_v2.shp')
cs_regs_raw = gpd.read_file(cs_regs_shp_dir)

coscat_id = id_loop_row['coscat_id'].values[0]
reg_id = id_loop_row['cst_reg_id'].values[0]

#   define the right string from the coscat id number
if coscat_id < 10:
    id_cs_str = '000' + str(int(coscat_id))
elif coscat_id >= 10 and coscat_id < 100:
    id_cs_str = '00' + str(int(coscat_id))
elif coscat_id >= 100 and coscat_id < 1000:
    id_cs_str = '0' + str(int(coscat_id))
else:
    id_cs_str = str(int(coscat_id))  

if reg_id < 10:
    str_id_subreg = '00' + str(int(reg_id))
elif reg_id >= 10 and reg_id < 100:
    str_id_subreg = '0' + str(int(reg_id))
else:
    str_id_subreg = str(int(reg_id))  

if reg_id != -1:
    #   define all the necessary sub directories
    coscat_dir = subreg_dir = os.path.join(main_dir, id_cs_str)
    subreg_id_lst = list(set(cs_regs_raw.loc[(cs_regs_raw['coscat'] == coscat_id)]['cst_reg_id'].values.tolist()))
    
    #   define the model directory where the model files and figures will be stored
    id_riv_bas_dir = os.path.join(subreg_dir, id_cs_str + '_REG_' + str_id_subreg)
    rch_pic_dir = os.path.join(id_riv_bas_dir, '_RCH_hist')
    geometry_pic_dir = os.path.join(id_riv_bas_dir, '_GEOMETRY_hist')
    geology_pic_dir = os.path.join(id_riv_bas_dir, '_GEOLOGY_hist')     
    
    if not os.path.exists(id_riv_bas_dir):
        os.makedirs(id_riv_bas_dir)
        os.makedirs(rch_pic_dir)
        os.makedirs(geometry_pic_dir)
        os.makedirs(geology_pic_dir)
    else:
        print('Directories for RIV BASIN id: ' + str_id_subreg + ' already exists.')

    #  determine how many points to create (25 means point per 1km) and the distance of the farthest point from point_5km
    n_points = 400
    max_dist = 200000
    #   define the number of points and distance from the cross-section point
    n_points_avg = 5
    avg_dist = 500
    
    #   define the list that will serve as output to csv (if necessary)
    to_csv = [[],[],[],[],[]]
    
    """ 1)  Open raster GEBCO 2014, connect to database and create new tables """
    ##  first open the raster file via GDAL
    raster_in = gdal.Open(gebco_2014_dir)
    
    ##  assign the geotransform, band and noval
    gt = raster_in.GetGeoTransform()
    rb = raster_in.GetRasterBand(1)
    noval = rb.GetNoDataValue()
    rx = raster_in.RasterXSize
    ry = raster_in.RasterYSize
    
    #   open the CSV file with the coastal points and other files
    df_points_5km = pd.read_csv(os.path.join(input_files_dir, '_csv_input_files', 'ne_10m_coastline_points_5km.csv'))
    df_sed_thk_v2 = pd.read_csv(os.path.join(input_files_dir, '_csv_input_files', 'cs_sed_thick_est_v2.csv'))
    df_cst_type = pd.read_csv(os.path.join(input_files_dir, '_csv_input_files', 'cs_coastal_types.csv'))
    
    #   open the nodes CSV file and transform it into geopandas dataframe
    df_nodes = pd.read_csv(os.path.join(input_files_dir, '_csv_input_files', 'ne_coastline_nodes_test.csv')) 
    geometry_nodes = [Point(xy) for xy in zip(df_nodes.x, df_nodes.y)]
    df_nodes = df_nodes.drop(['x', 'y'], axis = 1)
    crs = {'init': 'epsg:4326'}
    gdf_nodes = gpd.GeoDataFrame(df_nodes, crs = crs, geometry = geometry_nodes)
    nodes_geom = MultiPoint(geometry_nodes)
    
    csv_dir_ts = os.path.join(subreg_dir, '_SUBREG_representative_models_summary.csv')      
    csv_dir_ts_input_data_geom = os.path.join(subreg_dir, '_SUBREG_representative_models_GEOMETRY_stats.csv')   
    csv_dir_ts_input_data_rch = os.path.join(subreg_dir, '_SUBREG_representative_models_RCH_stats.csv')    
    csv_dir_ts_input_data_geol = os.path.join(subreg_dir, '_SUBREG_representative_models_GEOLOGY_stats.csv')    
    
    if os.path.exists(csv_dir_ts):
        pass
    else:
        f = open(csv_dir_ts,'w')
        f.write('id_COSCAT, id_SUBREG, n, fos, ls_elev, ls_elev_50') #write the header
        #   n = number of profiles in the HYBAS regions
        #   fos = FOS point found in the averaged Henry profile, if no = 0, if yes = distance from coast in km
        #               this is not based on the find_FOS function as the coastline can be shifted depending on the top_elev list,
        #               instead it is defined based on the last column in the IBOUND_array and the number of layers there 
        #   ls_elev = either 0 or 1 (false or true), if all the elevation points are > -120m offshore then = 0
        #   lst_elev_50 = either 0 or 1 (false or true), if all more than half of the the elevation points are > -120m offshore then = 0
        f.close() 
    
        f2 = open(csv_dir_ts_input_data_geom,'w')
        f2.write('id_COSCAT, id_SUBREG, cs_thk_mu, cs_thk_std, cs_width_mu, cs_width_std')
        f2.close()
        
        f3 = open(csv_dir_ts_input_data_rch,'w')
        f3.write('id_COSCAT, id_SUBREG, cs_pcr_rch_mu, cs_pcr_rch_std, cs_watergap_rch_mu, cs_watergap_rch_std, cs_p_min_et_rch_mu, cs_p_min_et_rch_std,\
        all_rch_mu, all_rch_std') #write the header
        f3.close()
        
        f4 = open(csv_dir_ts_input_data_geol,'w')
        f4.write('id_COSCAT, id_SUBREG, soil_type_mu, soil_type_std, soil_thk_mu, soil_thk_std,\
        glhymps_top_mu, glhymps_top_std, glhymps_bot_mu, glhymps_bot_std, soil_k_mu')         
        f4.close()


"""
cs_pts_dist = 100.
dir_out = id_riv_bas_dir
cont_dir = cont_direction
id_cs = 43818
"""

def merit_dem_topo(id_cs, cs_pts_dist, dir_out, cont_dir):
    #   get the x and y coordinates 
    #pts_info = [i for i in points_5km if i[0] == id_cs]
    pts_info = df_points_5km.loc[df_points_5km['id_cs'] == id_cs]
    
    #if pts_info != []:
    
    point_lat, point_lon = pts_info.x.values.tolist()[0], pts_info.y.values.tolist()[0]

    #   get the closest node to the coastal point
    nearest_geoms = nearest_points(Point(point_lat, point_lon), nodes_geom)
    near_idx1 = nearest_geoms[1]
    closest_node_lat, closest_node_lon = near_idx1.x, near_idx1.y        
    
    ##  calculate the distance between the point_5m and closest node
    ##  to do this we convert the lat/lon coordinates to UTM coordinates (using the utm package)
    point_5km_utm = utm.from_latlon(point_lon, point_lat)
    closest_node_utm = utm.from_latlon(closest_node_lon, closest_node_lat)
    zone_num, zone_let = point_5km_utm[2], point_5km_utm[3]

    ##  check that zone_num is between 1 and 60, case for areas that are split between
    ##  the zones 1 and 60 (e.g. certain islands in Pacific..)
    if zone_num > 60:
        zone_num = zone_num - 60

    ##  extract value on the coastline point - for GEBCO and GLIM
    gebco_val_cs_point = db_fc.read_raster_val(rb, gt, rx, ry, point_lat, point_lon)
    
    #   define the csv headers and list
    csv_headers = ['dist_from_cst', 'x_wgs84', 'y_wgs_84', 'gebco', 'gebco_avg', 'merit', 'merit_avg', 'merit_gebco', 'merit_gebco_avg']
    csv_lst = []
    
    """
    point_lat = c_coord_wgs_x_land
    point_lon = c_coord_wgs_y_land
    """
    
    def get_merit_val(point_lat, point_lon):

        #   first decide on the numbers 
        """
        merit_dem_lat = round(point_lat / 5.) * 5
        if point_lat < 0:
            merit_dem_lat = abs(merit_dem_lat) + 5
    
        merit_dem_lon = round(point_lon / 5.) * 5
        if point_lon < 0:
            merit_dem_lon = abs(merit_dem_lon) + 5    
        """
        
        if point_lat < 0:
            merit_dem_lat = math.ceil(point_lat / 5.) * 5
            merit_dem_lat = abs(merit_dem_lat) + 5
        else:
            merit_dem_lat = math.floor(point_lat / 5.) * 5
            
        if point_lon < 0:
            merit_dem_lon = math.floor(point_lon / 5.) * 5
            merit_dem_lon = abs(merit_dem_lon)# + 5    
        else:
            merit_dem_lon = math.floor(point_lon / 5.) * 5   
    
        #   transform to string
        merit_dem_lat = str(merit_dem_lat)
        merit_dem_lon = str(merit_dem_lon)
        if len(merit_dem_lon) == 1:
            merit_dem_lon = '0' + merit_dem_lon
        if len(merit_dem_lat) == 1:
            merit_dem_lat = '00' + merit_dem_lat
        elif len(merit_dem_lat) == 2:
            merit_dem_lat = '0' + merit_dem_lat
    
        #   add north/south and east/west based on the coordinates
        if point_lon < 0:
            merit_dem_lon = 's' + merit_dem_lon
        else:
            merit_dem_lon = 'n' + merit_dem_lon
        if point_lat < 0:
            merit_dem_lat = 'w' + merit_dem_lat
        else:
            merit_dem_lat = 'e' + merit_dem_lat
    
        #   open the MERIT DEM raster
        def find(name, path):
            for root, dirs, files in os.walk(path):
                if name in files:
                    return os.path.join(root, name)    
        merit_dir = find(merit_dem_lon + merit_dem_lat + '_dem.tif', merit_dem_folder_dir)
        
        if merit_dir:
            raster_in_merit = gdal.Open(merit_dir)
            ##  assign the geotransform, band and noval
            gt_merit = raster_in_merit.GetGeoTransform()
            rb_merit = raster_in_merit.GetRasterBand(1)
            rx_merit = raster_in_merit.RasterXSize
            ry_merit = raster_in_merit.RasterYSize
            merit_val_cs_point = db_fc.read_raster_val(rb_merit, gt_merit, rx_merit, ry_merit, point_lat, point_lon)

        else:
            merit_val_cs_point = np.nan
        
        return merit_val_cs_point

    #  determine how many points to create (25 means point per 1km) and the distance of the farthest point from point_5km
    max_dist = 200000
    n_points = int(max_dist / cs_pts_dist)
    x_coords = np.arange(-max_dist / 1000., max_dist / 1000. + cs_pts_dist / 1000., cs_pts_dist / 1000.)
    x_coords = [round(i, 2) for i in x_coords]
    #   define the number of points and distance from the cross-section point
    avg_dist = 500
    n_points_avg = int(avg_dist / cs_pts_dist)
    
    #gebco_avg_land_lst, gebco_avg_sea_lst = np.ones(len(x_coords)) * np.nan, np.ones(len(x_coords)) * np.nan
    #merit_avg_land_lst, merit_avg_sea_lst = np.ones(len(x_coords)) * np.nan, np.ones(len(x_coords)) * np.nan
    
    merit_val_cs_point = get_merit_val(point_lat, point_lon)

    if merit_val_cs_point < -100. or math.isnan(merit_val_cs_point):
        merit_val_cs_point = gebco_val_cs_point
    
    gebco_avg_lst = np.ones(len(x_coords)) * np.nan
    merit_avg_lst = np.ones(len(x_coords)) * np.nan
    gebco_lst = np.ones(len(x_coords)) * np.nan
    merit_lst = np.ones(len(x_coords)) * np.nan
    merit_gebco_lst = np.ones(len(x_coords)) * np.nan
    merit_gebco_avg_lst = np.ones(len(x_coords)) * np.nan
    
    gebco_lst[int(len(x_coords) / 2)] = gebco_val_cs_point
    gebco_avg_lst[int(len(x_coords) / 2)] = gebco_val_cs_point
    merit_lst[int(len(x_coords) / 2)] = merit_val_cs_point
    merit_avg_lst[int(len(x_coords) / 2)] = merit_val_cs_point
    merit_gebco_lst[int(len(x_coords) / 2)] = merit_val_cs_point
    merit_gebco_avg_lst[int(len(x_coords) / 2)] = merit_val_cs_point
    
    csv_lst.append([0, point_lat, point_lon, gebco_val_cs_point, gebco_val_cs_point, merit_val_cs_point, merit_val_cs_point])
    
    ##  run the loop for all the cross-section points
    for a in range(n_points):
        try:
            ##  calculate the b_length
            b_length = max_dist - a * cs_pts_dist
            c_coords_utm = db_fc.coord_from_triangle(point_5km_utm[0], point_5km_utm[1], closest_node_utm[0], closest_node_utm[1], b_length)
            #print(b_length)

            ##  transform the coordinates back to wgs84, position of the equator is taken into account
            c_coords_wgs = db_fc.equator_position_angles(c_coords_utm[1], c_coords_utm[0], c_coords_utm[3], c_coords_utm[2], zone_num, zone_let)
            c_coord_wgs_x_land, c_coord_wgs_y_land, c_coord_wgs_x_sea, c_coord_wgs_y_sea = \
            c_coords_wgs[0][1], c_coords_wgs[0][0], c_coords_wgs[1][1], c_coords_wgs[1][0]

            ##  read values from the GEBCO raster
            gebco_val_land = db_fc.read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_land, c_coord_wgs_y_land)
            gebco_val_sea = db_fc.read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_sea, c_coord_wgs_y_sea)
            
            ##  read values from the MERIT DEM raster
            merit_val_land = round(get_merit_val(c_coord_wgs_x_land, c_coord_wgs_y_land), 2)
            merit_val_sea = round(get_merit_val(c_coord_wgs_x_sea, c_coord_wgs_y_sea) ,2)
            
            ##  define sum and count of avg_values for final average value
            sum_avg_land_gebco, sum_avg_sea_gebco, sum_avg_land_merit, sum_avg_sea_merit = [], [], [], []
            ##  calculate the average value for each cross-section point
            for b in range(n_points_avg):
                avg_point_dist = avg_dist - b * cs_pts_dist

                avg_point_coord_land = db_fc.avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[1], c_coords_utm[0], avg_point_dist)
                avg_point_coord_sea = db_fc.avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[3], c_coords_utm[2], avg_point_dist)

                ##  transform to wgs coordinates
                avg_point_coord_wgs_land = db_fc.equator_position_angles(avg_point_coord_land[1], avg_point_coord_land[0], avg_point_coord_land[3], avg_point_coord_land[2], zone_num, zone_let)
                avg_point_coord_wgs_sea = db_fc.equator_position_angles(avg_point_coord_sea[1], avg_point_coord_sea[0], avg_point_coord_sea[3], avg_point_coord_sea[2], zone_num, zone_let)

                ##  GEBCO get the pixel value for the cross-section avg point for the land pixel
                avg_pixel_val_land_1 = db_fc.read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_land[0][1], avg_point_coord_wgs_land[0][0])
                sum_avg_land_gebco.append(avg_pixel_val_land_1)
                avg_pixel_val_land_2 = db_fc.read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_land[1][1], avg_point_coord_wgs_land[1][0])
                sum_avg_land_gebco.append(avg_pixel_val_land_2)
                ##  get the pixel value for the cross-section avg point for the sea pixel
                avg_pixel_val_sea_1 = db_fc.read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_sea[0][1], avg_point_coord_wgs_sea[0][0])
                sum_avg_sea_gebco.append(avg_pixel_val_sea_1)
                avg_pixel_val_sea_2 = db_fc.read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_sea[1][1], avg_point_coord_wgs_sea[1][0])
                sum_avg_sea_gebco.append(avg_pixel_val_sea_2)
            
                ##  MERIT get the pixel value for the cross-section avg point for the land pixel
                avg_pixel_val_land_merit_1 = get_merit_val(avg_point_coord_wgs_land[0][1], avg_point_coord_wgs_land[0][0])
                sum_avg_land_merit.append(avg_pixel_val_land_merit_1)
                avg_pixel_val_land_merit_2 = get_merit_val(avg_point_coord_wgs_land[1][1], avg_point_coord_wgs_land[1][0])
                sum_avg_land_merit.append(avg_pixel_val_land_merit_2)
                ##  get the pixel value for the cross-section avg point for the sea pixel
                avg_pixel_val_sea_merit_1 = get_merit_val(avg_point_coord_wgs_sea[0][1], avg_point_coord_wgs_sea[0][0])
                sum_avg_sea_merit.append(avg_pixel_val_sea_merit_1)
                avg_pixel_val_sea_merit_2 = get_merit_val(avg_point_coord_wgs_sea[1][1], avg_point_coord_wgs_sea[1][0])
                sum_avg_sea_merit.append(avg_pixel_val_sea_merit_2)            
            
            #   get rid of potential non-value data in the MERIT list
            sum_avg_land_merit = [i for i in sum_avg_land_merit if i > -100.]
            sum_avg_sea_merit = [i for i in sum_avg_sea_merit if i > -100.]
            
            gebco_avg_land, gebco_avg_sea = int(sum(sum_avg_land_gebco) / len(sum_avg_land_gebco)), int(sum(sum_avg_sea_gebco) / len(sum_avg_sea_gebco))
            try:
                merit_avg_land = round(sum(sum_avg_land_merit) / len(sum_avg_land_merit), 2)
            except ZeroDivisionError:
                merit_avg_land = np.nan
            try: 
                merit_avg_sea = round(sum(sum_avg_sea_merit) / len(sum_avg_sea_merit), 2)
            except ZeroDivisionError:
                merit_avg_sea = np.nan
                
            gebco_avg_lst[a]= gebco_avg_sea
            gebco_avg_lst[n_points + n_points - a] = gebco_avg_land            
            merit_avg_lst[a]= merit_avg_sea
            merit_avg_lst[n_points + n_points - a] = merit_avg_land

            gebco_lst[a]= gebco_val_sea
            gebco_lst[n_points + n_points - a] = gebco_val_land            
            merit_lst[a]= merit_val_sea
            merit_lst[n_points + n_points - a] = merit_val_land

            #   merge the MERIT and GEBCO data into the same list - so it can be potentially tested in the model. Take all the values 
            #   of MERIT wherever possible and fill in the non-data values by the GEBCO data.
            
            if merit_val_land < -100. or math.isnan(merit_val_land):
                merit_gebco_land = gebco_val_land
            else:
                merit_gebco_land = merit_val_land
            if merit_val_sea < -100. or math.isnan(merit_val_sea):
                merit_gebco_sea = gebco_val_sea
            else:
                merit_gebco_sea = merit_val_sea

            if merit_avg_land < -100. or math.isnan(merit_avg_land):
                merit_gebco_avg_land = gebco_avg_land
            else:
                merit_gebco_avg_land = merit_avg_land
            if merit_avg_sea < -100. or math.isnan(merit_avg_sea):
                merit_gebco_avg_sea = gebco_avg_sea
            else:
                merit_gebco_avg_sea = merit_avg_sea            
            
            merit_gebco_lst[a] = merit_gebco_sea
            merit_gebco_avg_lst[a] = merit_gebco_avg_sea
            merit_gebco_lst[n_points + n_points - a] = merit_gebco_land
            merit_gebco_avg_lst[n_points + n_points - a] = merit_gebco_avg_land

            csv_lst.append([-b_length, c_coord_wgs_x_land, c_coord_wgs_y_land, gebco_val_land, gebco_avg_land,\
                            merit_val_land, merit_avg_land, merit_gebco_land, merit_gebco_avg_land])
            csv_lst.append([b_length, c_coord_wgs_x_sea, c_coord_wgs_y_sea, gebco_val_sea, gebco_avg_sea,\
                            merit_val_sea, merit_avg_sea, merit_gebco_sea, merit_gebco_avg_sea])            
    
        ##  ZeroDivisionError happens when point is located on top of the node (closest node)
        ##  this is the case in small islands where equidistant points are created..
        except ZeroDivisionError:
            continue
            #print('Division by zero')    
    
    #   still make sure that there are no -9999m points in the merit lists
    merit_lst = [np.nan if i < -100. else round(i, 2) for i in merit_lst]   
    merit_avg_lst = [np.nan if i < -100. else round(i, 2) for i in merit_avg_lst]       
    merit_gebco_lst = [np.nan if i == -9999. else round(i, 2) for i in merit_gebco_lst]   
    merit_gebco_avg_lst = [np.nan if i == -9999. else round(i, 2) for i in merit_gebco_avg_lst]    
    
    #   create a netcdf file that will be returned as output
    xa_sum = xr.Dataset(data_vars = {'GEBCO_elev' : (('x'), gebco_lst),
                                     'GEBCO_elev_AVG' : (('x'), gebco_avg_lst),
                                     'MERIT_elev' : (('x'), merit_lst),
                                     'MERIT_elev_AVG' : (('x'), merit_avg_lst),
                                     'MERIT_GEBCO_elev' : (('x'), merit_gebco_lst),
                                     'MERIT_GEBCO_elev_AVG' : (('x'), merit_gebco_avg_lst)},
                        coords = {'x' : x_coords})

    xa_name = 'cs_' + str(id_cs) + '_pt_dist_' + str(int(cs_pts_dist)) + '_topo_merit.nc'
    xa_sum.to_netcdf(os.path.join(dir_out, xa_name))           

    #   save the csv file
    my_df = pd.DataFrame(csv_lst)
    csv_name = 'cs_' + str(id_cs) + '_pt_dist_' + str(int(cs_pts_dist)) + '_topo_merit.csv'
    my_df.to_csv(os.path.join(dir_out, csv_name), index = False, header = csv_headers)

    #   also plot the topography in one plot
    fig, ax = plt.subplots(figsize = (15, 8))
    
    #   calculate the indexes so the zoom area goes from 10km offshore to 50km inland
    x_st_idx = int(200000. / (cs_pts_dist) - 50000. / cs_pts_dist)
    x_end_idx = int(200000. / (cs_pts_dist) + 10000. / cs_pts_dist)
    
    if cont_dir == 'right':
        ax.plot(x_coords[x_st_idx : x_end_idx], list(reversed(gebco_lst[x_st_idx : x_end_idx])), c = 'deepskyblue', label = 'GEBCO', lw = 1)
        ax.plot(x_coords[x_st_idx : x_end_idx], list(reversed(gebco_avg_lst[x_st_idx : x_end_idx])), c = 'blue', label = 'GEBCO_avg', lw = 2)
        ax.plot(x_coords[x_st_idx : x_end_idx], list(reversed(merit_lst[x_st_idx : x_end_idx])), c = 'gold', label = 'MERIT', lw = 1)
        ax.plot(x_coords[x_st_idx : x_end_idx], list(reversed(merit_avg_lst[x_st_idx : x_end_idx])), c = 'darkorange', label = 'MERIT_avg', lw = 2)
    
    else:   
        ax.plot(x_coords[x_st_idx : x_end_idx], gebco_lst[x_st_idx : x_end_idx], c = 'deepskyblue', label = 'GEBCO', lw = 1)
        ax.plot(x_coords[x_st_idx : x_end_idx], gebco_avg_lst[x_st_idx : x_end_idx], c = 'blue', label = 'GEBCO_avg', lw = 2)
        ax.plot(x_coords[x_st_idx : x_end_idx], merit_lst[x_st_idx : x_end_idx], c = 'gold', label = 'MERIT', lw = 1)
        ax.plot(x_coords[x_st_idx : x_end_idx], merit_avg_lst[x_st_idx : x_end_idx], c = 'darkorange', label = 'MERIT_avg', lw = 2)
        
    ax.legend()
    ax.set(xlabel = 'Distance from coast (km)', ylabel = 'Elevation (m asl.)', title = 'Comparison of MERIT and GEBCO elevation profiles for CS = ' + str(id_cs))
    ax.grid()
    plt_name_cst_zoom = 'cs_' + str(id_cs) + '_pt_dist_' + str(int(cs_pts_dist)) + '_topo_CST_ZOOM.png'
    fig.savefig(os.path.join(dir_out, plt_name_cst_zoom), dpi = 300, facecolor = 'w', edgecolor = 'w', orientation = 'portrait', 
                papertype = None, format = None, transparent = False, bbox_inches = 'tight', pad_inches = 0.3, frameon = None)

    return gebco_lst, gebco_avg_lst, merit_lst, merit_avg_lst, merit_gebco_lst, merit_gebco_avg_lst

"""
cs_pts_dist = 100.
dir_out = id_riv_bas_dir

cont_dir = cont_direction

"""
def coastal_dem_topo(id_cs, cs_pts_dist, dir_out, cont_dir):
    #   get the x and y coordinates 
    pts_info = df_points_5km.loc[df_points_5km['id_cs'] == id_cs]
    
    #if pts_info != []:
    point_lat, point_lon = pts_info.x.values.tolist()[0], pts_info.y.values.tolist()[0]

    #   get the closest node to the coastal point
    nearest_geoms = nearest_points(Point(point_lat, point_lon), nodes_geom)
    near_idx1 = nearest_geoms[1]
    closest_node_lat, closest_node_lon = near_idx1.x, near_idx1.y        
    
    ##  calculate the distance between the point_5m and closest node
    ##  to do this we convert the lat/lon coordinates to UTM coordinates (using the utm package)
    point_5km_utm = utm.from_latlon(point_lon, point_lat)
    closest_node_utm = utm.from_latlon(closest_node_lon, closest_node_lat)
    zone_num, zone_let = point_5km_utm[2], point_5km_utm[3]

    ##  check that zone_num is between 1 and 60, case for areas that are split between
    ##  the zones 1 and 60 (e.g. certain islands in Pacific..)
    if zone_num > 60:
        zone_num = zone_num - 60

    ##  extract value on the coastline point - for GEBCO and GLIM
    gebco_val_cs_point = db_fc.read_raster_val(rb, gt, rx, ry, point_lat, point_lon)
    
    #   define the csv headers and list
    csv_headers = ['dist_from_cst', 'x_wgs84', 'y_wgs_84', 'gebco', 'gebco_avg', 'merit', 'merit_avg', 'merit_gebco', 'merit_gebco_avg']
    csv_lst = []
    
    """
    point_lat = point_lat
    point_lon = point_lon
    point_lat = avg_point_coord_wgs_land[0][1]
    point_lon = avg_point_coord_wgs_land[0][0]
    """

    def get_coastaldem_val(point_lat, point_lon):

        #   first decide on the numbers 
        """
        coastal_dem_lat = math.floor(point_lat / 1.) * 1
        if point_lat < 0:
            coastal_dem_lat = abs(coastal_dem_lat) + 1
    
        coastal_dem_lon = math.floor(point_lon / 1.) * 1
        if point_lon < 0:
            coastal_dem_lon = abs(coastal_dem_lon) + 1    
        """

        if point_lat < 0:
            coastal_dem_lat = math.ceil(point_lat / 1.) * 1
            coastal_dem_lat = abs(coastal_dem_lat) + 1
        else:
            coastal_dem_lat = math.floor(point_lat / 1.) * 1
            
        if point_lon < 0:
            coastal_dem_lon = math.floor(point_lon / 1.) * 1
            coastal_dem_lon = abs(coastal_dem_lon) #+ 1    
        else:
            coastal_dem_lon = math.floor(point_lon / 1.) * 1       
    
        #   transform to string
        coastal_dem_lat = str(coastal_dem_lat)
        coastal_dem_lon = str(coastal_dem_lon)
        if len(coastal_dem_lon) == 1:
            coastal_dem_lon = '0' + coastal_dem_lon
        if len(coastal_dem_lat) == 1:
            coastal_dem_lat = '00' + coastal_dem_lat
        elif len(coastal_dem_lat) == 2:
            coastal_dem_lat = '0' + coastal_dem_lat
    
        #   add north/south and east/west based on the coordinates
        if point_lon < 0:
            coastal_dem_lon = 'S' + coastal_dem_lon
        else:
            coastal_dem_lon = 'N' + coastal_dem_lon
        if point_lat < 0:
            coastal_dem_lat = 'W' + coastal_dem_lat
        else:
            coastal_dem_lat = 'E' + coastal_dem_lat
    
        #   open the MERIT DEM raster
        def find(name, path):
            for root, dirs, files in os.walk(path):
                if name in files:
                    return os.path.join(root, name)    
        coastal_dir = find(coastal_dem_lon + coastal_dem_lat + '.tif', coastal_dem_folder_dir)
        #print(coastal_dir)
        if coastal_dir:
            #print('dir exists')
            raster_in_coastal = gdal.Open(coastal_dir)
            ##  assign the geotransform, band and noval
            gt_coastal_orig = raster_in_coastal.GetGeoTransform()

            #   check that the pixel size is not 8.333333333333452e-05 because thats the wrong resoltion, has to be one order of magnitude lower
            if gt_coastal_orig[1] < 1e-04:
                #   reproject the raster after moving the previous dir into a new folder
                new_dir = os.path.join(coastal_dem_folder_corr_res_dir, coastal_dem_lon + coastal_dem_lat + '.tif')
                
                if not os.path.exists(new_dir):
                    rb_coastal_orig = raster_in_coastal.GetRasterBand(1)
                    xmin, ymax = gt_coastal_orig[0], gt_coastal_orig[3] + 1200 * gt_coastal_orig[-1] + 1.
                    xres = 8.333333333333452e-04
                    yres = 8.333333333333452e-04
                    geotransform = (xmin, xres, 0, ymax, 0, -yres)  
                    output_raster = gdal.GetDriverByName('GTiff').Create(new_dir, 1200, 1200, 1 ,gdal.GDT_Float32)  # Open the file
                    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
                    srs = osr.SpatialReference()                 # Establish its coordinate encoding
                    srs.ImportFromEPSG(4326)      
                    output_raster.SetProjection(srs.ExportToWkt())
                    output_raster.GetRasterBand(1).WriteArray(rb_coastal_orig.ReadAsArray())
                    output_raster.FlushCache()
                    output_raster = None

                raster_in_coastal = gdal.Open(new_dir)
                gt_coastal = raster_in_coastal.GetGeoTransform()
                rb_coastal = raster_in_coastal.GetRasterBand(1)
                rx_coastal = raster_in_coastal.RasterXSize
                ry_coastal = raster_in_coastal.RasterYSize

            else:
                gt_coastal = gt_coastal_orig
                rb_coastal = raster_in_coastal.GetRasterBand(1)
                rx_coastal = raster_in_coastal.RasterXSize
                ry_coastal = raster_in_coastal.RasterYSize

            coastal_val_cs_point = db_fc.read_raster_val(rb_coastal, gt_coastal, rx_coastal, ry_coastal, point_lat, point_lon)

        else:
            coastal_val_cs_point = np.nan
        
        return coastal_val_cs_point
    
    
    #  determine how many points to create (25 means point per 1km) and the distance of the farthest point from point_5km
    max_dist = 200000
    n_points = int(max_dist / cs_pts_dist)
    x_coords = np.arange(-max_dist / 1000., max_dist / 1000. + cs_pts_dist / 1000., cs_pts_dist / 1000.)
    x_coords = [round(i, 2) for i in x_coords]
    #   define the number of points and distance from the cross-section point
    avg_dist = 500
    n_points_avg = int(avg_dist / cs_pts_dist)
    
    #gebco_avg_land_lst, gebco_avg_sea_lst = np.ones(len(x_coords)) * np.nan, np.ones(len(x_coords)) * np.nan
    #merit_avg_land_lst, merit_avg_sea_lst = np.ones(len(x_coords)) * np.nan, np.ones(len(x_coords)) * np.nan
    
    coastaldem_val_cs_point = get_coastaldem_val(point_lat, point_lon)

    if coastaldem_val_cs_point == 0. or math.isnan(coastaldem_val_cs_point):
        coastaldem_val_cs_point = gebco_val_cs_point

    gebco_avg_lst = np.ones(len(x_coords)) * np.nan
    coastaldem_avg_lst = np.ones(len(x_coords)) * np.nan
    gebco_lst = np.ones(len(x_coords)) * np.nan
    coastaldem_lst = np.ones(len(x_coords)) * np.nan
    coastaldem_gebco_lst = np.ones(len(x_coords)) * np.nan
    coastaldem_gebco_avg_lst = np.ones(len(x_coords)) * np.nan
    
    gebco_lst[int(len(x_coords) / 2)] = gebco_val_cs_point
    gebco_avg_lst[int(len(x_coords) / 2)] = gebco_val_cs_point
    coastaldem_lst[int(len(x_coords) / 2)] = coastaldem_val_cs_point
    coastaldem_avg_lst[int(len(x_coords) / 2)] = coastaldem_val_cs_point
    coastaldem_gebco_lst[int(len(x_coords) / 2)] = coastaldem_val_cs_point
    coastaldem_gebco_avg_lst[int(len(x_coords) / 2)] = coastaldem_val_cs_point
    
    csv_lst.append([0, point_lat, point_lon, gebco_val_cs_point, gebco_val_cs_point, coastaldem_val_cs_point, coastaldem_val_cs_point])
    
    ##  run the loop for all the cross-section points
    for a in range(n_points):
        try:
            ##  calculate the b_length
            b_length = max_dist - a * cs_pts_dist
            c_coords_utm = db_fc.coord_from_triangle(point_5km_utm[0], point_5km_utm[1], closest_node_utm[0], closest_node_utm[1], b_length)
            #print(b_length)

            ##  transform the coordinates back to wgs84, position of the equator is taken into account
            c_coords_wgs = db_fc.equator_position_angles(c_coords_utm[1], c_coords_utm[0], c_coords_utm[3], c_coords_utm[2], zone_num, zone_let)
            c_coord_wgs_x_land, c_coord_wgs_y_land, c_coord_wgs_x_sea, c_coord_wgs_y_sea = \
            c_coords_wgs[0][1], c_coords_wgs[0][0], c_coords_wgs[1][1], c_coords_wgs[1][0]

            ##  read values from the GEBCO raster
            gebco_val_land = db_fc.read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_land, c_coord_wgs_y_land)
            gebco_val_sea = db_fc.read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_sea, c_coord_wgs_y_sea)
            
            ##  read values from the MERIT DEM raster
            coastaldem_val_land = round(get_coastaldem_val(c_coord_wgs_x_land, c_coord_wgs_y_land), 2)
            coastaldem_val_sea = round(get_coastaldem_val(c_coord_wgs_x_sea, c_coord_wgs_y_sea) ,2)
            
            ##  define sum and count of avg_values for final average value
            sum_avg_land_gebco, sum_avg_sea_gebco, sum_avg_land_coastaldem, sum_avg_sea_coastaldem = [], [], [], []
            ##  calculate the average value for each cross-section point
            for b in range(n_points_avg):
                avg_point_dist = avg_dist - b * cs_pts_dist

                avg_point_coord_land = db_fc.avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[1], c_coords_utm[0], avg_point_dist)
                avg_point_coord_sea = db_fc.avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[3], c_coords_utm[2], avg_point_dist)

                ##  transform to wgs coordinates
                avg_point_coord_wgs_land = db_fc.equator_position_angles(avg_point_coord_land[1], avg_point_coord_land[0], avg_point_coord_land[3], avg_point_coord_land[2], zone_num, zone_let)
                avg_point_coord_wgs_sea = db_fc.equator_position_angles(avg_point_coord_sea[1], avg_point_coord_sea[0], avg_point_coord_sea[3], avg_point_coord_sea[2], zone_num, zone_let)

                ##  GEBCO get the pixel value for the cross-section avg point for the land pixel
                avg_pixel_val_land_1 = db_fc.read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_land[0][1], avg_point_coord_wgs_land[0][0])
                sum_avg_land_gebco.append(avg_pixel_val_land_1)
                avg_pixel_val_land_2 = db_fc.read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_land[1][1], avg_point_coord_wgs_land[1][0])
                sum_avg_land_gebco.append(avg_pixel_val_land_2)
                ##  get the pixel value for the cross-section avg point for the sea pixel
                avg_pixel_val_sea_1 = db_fc.read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_sea[0][1], avg_point_coord_wgs_sea[0][0])
                sum_avg_sea_gebco.append(avg_pixel_val_sea_1)
                avg_pixel_val_sea_2 = db_fc.read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_sea[1][1], avg_point_coord_wgs_sea[1][0])
                sum_avg_sea_gebco.append(avg_pixel_val_sea_2)
            
                ##  MERIT get the pixel value for the cross-section avg point for the land pixel
                avg_pixel_val_land_coastaldem_1 = get_coastaldem_val(avg_point_coord_wgs_land[0][1], avg_point_coord_wgs_land[0][0])
                sum_avg_land_coastaldem.append(avg_pixel_val_land_coastaldem_1)
                avg_pixel_val_land_coastaldem_2 = get_coastaldem_val(avg_point_coord_wgs_land[1][1], avg_point_coord_wgs_land[1][0])
                sum_avg_land_coastaldem.append(avg_pixel_val_land_coastaldem_2)
                ##  get the pixel value for the cross-section avg point for the sea pixel
                avg_pixel_val_sea_coastaldem_1 = get_coastaldem_val(avg_point_coord_wgs_sea[0][1], avg_point_coord_wgs_sea[0][0])
                sum_avg_sea_coastaldem.append(avg_pixel_val_sea_coastaldem_1)
                avg_pixel_val_sea_coastaldem_2 = get_coastaldem_val(avg_point_coord_wgs_sea[1][1], avg_point_coord_wgs_sea[1][0])
                sum_avg_sea_coastaldem.append(avg_pixel_val_sea_coastaldem_2)            
            
            #   get rid of potential non-value data in the MERIT list
            sum_avg_land_coastaldem = [i for i in sum_avg_land_coastaldem if i > -100.]
            sum_avg_sea_coastaldem = [i for i in sum_avg_sea_coastaldem if i > -100.]
            
            gebco_avg_land, gebco_avg_sea = int(sum(sum_avg_land_gebco) / len(sum_avg_land_gebco)), int(sum(sum_avg_sea_gebco) / len(sum_avg_sea_gebco))
            try:
                coastaldem_avg_land = round(sum(sum_avg_land_coastaldem) / len(sum_avg_land_coastaldem), 2)
            except ZeroDivisionError:
                coastaldem_avg_land = np.nan
            try: 
                coastaldem_avg_sea = round(sum(sum_avg_sea_coastaldem) / len(sum_avg_sea_coastaldem), 2)
            except ZeroDivisionError:
                coastaldem_avg_sea = np.nan
                
            gebco_avg_lst[a]= gebco_avg_sea
            gebco_avg_lst[n_points + n_points - a] = gebco_avg_land            
            coastaldem_avg_lst[a]= coastaldem_avg_sea
            coastaldem_avg_lst[n_points + n_points - a] = coastaldem_avg_land

            gebco_lst[a]= gebco_val_sea
            gebco_lst[n_points + n_points - a] = gebco_val_land            
            coastaldem_lst[a]= coastaldem_val_sea
            coastaldem_lst[n_points + n_points - a] = coastaldem_val_land

            #   merge the MERIT and GEBCO data into the same list - so it can be potentially tested in the model. Take all the values 
            #   of MERIT wherever possible and fill in the non-data values by the GEBCO data.
            
            if coastaldem_val_land == 0. or math.isnan(coastaldem_val_land):
                coastaldem_gebco_land = gebco_val_land
            else:
                coastaldem_gebco_land = coastaldem_val_land
            if coastaldem_val_sea == 0. or math.isnan(coastaldem_val_sea):
                coastaldem_gebco_sea = gebco_val_sea
            else:
                coastaldem_gebco_sea = coastaldem_val_sea

            if coastaldem_avg_land == 0. or math.isnan(coastaldem_avg_land):
                coastaldem_gebco_avg_land = gebco_avg_land
            else:
                coastaldem_gebco_avg_land = coastaldem_avg_land
            if coastaldem_avg_sea == 0. or math.isnan(coastaldem_avg_sea):
                coastaldem_gebco_avg_sea = gebco_avg_sea
            else:
                coastaldem_gebco_avg_sea = coastaldem_avg_sea            
            
            coastaldem_gebco_lst[a] = coastaldem_gebco_sea
            coastaldem_gebco_avg_lst[a] = coastaldem_gebco_avg_sea
            coastaldem_gebco_lst[n_points + n_points - a] = coastaldem_gebco_land
            coastaldem_gebco_avg_lst[n_points + n_points - a] = coastaldem_gebco_avg_land

            csv_lst.append([-b_length, c_coord_wgs_x_land, c_coord_wgs_y_land, gebco_val_land, gebco_avg_land,\
                            coastaldem_val_land, coastaldem_avg_land, coastaldem_gebco_land, coastaldem_gebco_avg_land])
            csv_lst.append([b_length, c_coord_wgs_x_sea, c_coord_wgs_y_sea, gebco_val_sea, gebco_avg_sea,\
                            coastaldem_val_sea, coastaldem_avg_sea, coastaldem_gebco_sea, coastaldem_gebco_avg_sea])            
    
        ##  ZeroDivisionError happens when point is located on top of the node (closest node)
        ##  this is the case in small islands where equidistant points are created..
        except ZeroDivisionError:
            continue
            #print('Division by zero')    
    
    #   still make sure that there are no -9999m points in the merit lists
    coastaldem_lst = [np.nan if i == 0. else round(i, 2) for i in coastaldem_lst]   
    coastaldem_avg_lst = [np.nan if i == 0. else round(i, 2) for i in coastaldem_avg_lst]       
    coastaldem_gebco_lst = [np.nan if i == -9999. else round(i, 2) for i in coastaldem_gebco_lst]   
    coastaldem_gebco_avg_lst = [np.nan if i == -9999. else round(i, 2) for i in coastaldem_gebco_avg_lst]    
    
    #   create a netcdf file that will be returned as output
    xa_sum = xr.Dataset(data_vars = {'GEBCO_elev' : (('x'), gebco_lst),
                                     'GEBCO_elev_AVG' : (('x'), gebco_avg_lst),
                                     'CSDEM_elev' : (('x'), coastaldem_lst),
                                     'CSDEM_elev_AVG' : (('x'), coastaldem_avg_lst),
                                     'CSDEM_GEBCO_elev' : (('x'), coastaldem_gebco_lst),
                                     'CSDEM_GEBCO_elev_AVG' : (('x'), coastaldem_gebco_avg_lst)},
                        coords = {'x' : x_coords})

    xa_name = 'cs_' + str(id_cs) + '_pt_dist_' + str(int(cs_pts_dist)) + '_topo_csdem.nc'
    xa_sum.to_netcdf(os.path.join(dir_out, xa_name))           

    #   save the csv file
    my_df = pd.DataFrame(csv_lst)
    csv_name = 'cs_' + str(id_cs) + '_pt_dist_' + str(int(cs_pts_dist)) + '_topo_csdem.csv'
    my_df.to_csv(os.path.join(dir_out, csv_name), index = False, header = csv_headers)

    #   also plot the topography in one plot
    fig, ax = plt.subplots(figsize = (15, 8))
    
    #   calculate the indexes so the zoom area goes from 10km offshore to 50km inland
    x_st_idx = int(200000. / (cs_pts_dist) - 50000. / cs_pts_dist)
    x_end_idx = int(200000. / (cs_pts_dist) + 10000. / cs_pts_dist)
    
    
    if cont_dir == 'right':
        ax.plot(x_coords[x_st_idx : x_end_idx], list(reversed(gebco_lst[x_st_idx : x_end_idx])), c = 'deepskyblue', label = 'GEBCO', lw = 1)
        ax.plot(x_coords[x_st_idx : x_end_idx], list(reversed(gebco_avg_lst[x_st_idx : x_end_idx])), c = 'blue', label = 'GEBCO_avg', lw = 2)
        ax.plot(x_coords[x_st_idx : x_end_idx], list(reversed(coastaldem_lst[x_st_idx : x_end_idx])), c = 'limegreen', label = 'CoastalDEM', lw = 1)
        ax.plot(x_coords[x_st_idx : x_end_idx], list(reversed(coastaldem_avg_lst[x_st_idx : x_end_idx])), c = 'green', label = 'CoastalDEM_avg', lw = 2)
    
    else:   
        ax.plot(x_coords[x_st_idx : x_end_idx], gebco_lst[x_st_idx : x_end_idx], c = 'deepskyblue', label = 'GEBCO', lw = 1)
        ax.plot(x_coords[x_st_idx : x_end_idx], gebco_avg_lst[x_st_idx : x_end_idx], c = 'blue', label = 'GEBCO_avg', lw = 2)
        ax.plot(x_coords[x_st_idx : x_end_idx], coastaldem_lst[x_st_idx : x_end_idx], c = 'limegreen', label = 'CoastalDEM', lw = 1)
        ax.plot(x_coords[x_st_idx : x_end_idx], coastaldem_avg_lst[x_st_idx : x_end_idx], c = 'green', label = 'CoastalDEM_avg', lw = 2)
    
    ax.legend()
    ax.set(xlabel = 'Distance from coast (km)', ylabel = 'Elevation (m asl.)', title = 'Comparison of MERIT and GEBCO elevation profiles for CS = ' + str(id_cs))
    ax.grid()
    plt_name_cst_zoom = 'cs_' + str(id_cs) + '_pt_dist_' + str(int(cs_pts_dist)) + '_topo_CoastalDEM_CST_ZOOM.png'
    fig.savefig(os.path.join(dir_out, plt_name_cst_zoom), dpi = 300, facecolor = 'w', edgecolor = 'w', orientation = 'portrait', 
                papertype = None, format = None, transparent = False, bbox_inches = 'tight', pad_inches = 0.3, frameon = None)

    return gebco_lst, gebco_avg_lst, coastaldem_lst, coastaldem_avg_lst, coastaldem_gebco_lst, coastaldem_gebco_avg_lst




"""
function that reads in gw_rch for all time steps for each cross-section
"""
def get_gw_rch(id_cs, cs_pts_dist, dir_out, gw_rch_dir):
    
    pts_info = df_points_5km.loc[df_points_5km['id_cs'] == id_cs]

    time_steps = np.arange(-30000, -1000, 1000).tolist()
    time_steps_2 = np.arange(-1900, 500, 100).tolist()
    time_steps = time_steps + time_steps_2 
    
    #if pts_info != []:
    point_lat, point_lon = pts_info.x.values.tolist()[0], pts_info.y.values.tolist()[0]

    #   get the closest node to the coastal point
    nearest_geoms = nearest_points(Point(point_lat, point_lon), nodes_geom)
    near_idx1 = nearest_geoms[1]
    closest_node_lat, closest_node_lon = near_idx1.x, near_idx1.y        
    
    ##  calculate the distance between the point_5m and closest node
    ##  to do this we convert the lat/lon coordinates to UTM coordinates (using the utm package)
    point_5km_utm = utm.from_latlon(point_lon, point_lat)
    closest_node_utm = utm.from_latlon(closest_node_lon, closest_node_lat)
    zone_num, zone_let = point_5km_utm[2], point_5km_utm[3]

    ##  check that zone_num is between 1 and 60, case for areas that are split between
    ##  the zones 1 and 60 (e.g. certain islands in Pacific..)
    if zone_num > 60:
        zone_num = zone_num - 60

    gw_rch_all_ts_lst, gw_rch_avg_all_ts_lst = [], []
    
    """
    point_lat = point_lat
    point_lon = point_lon
    """

    def get_gw_rch_val(rb_coastal, gt_coastal, rx_coastal, ry_coastal, point_lat, point_lon, raster_in):

        #   first decide on the numbers 
        coastal_dem_lat = math.floor(point_lat / 1.) * 1
        if point_lat < 0:
            coastal_dem_lat = abs(coastal_dem_lat) + 1
    
        coastal_dem_lon = math.floor(point_lon / 1.) * 1
        if point_lon < 0:
            coastal_dem_lon = abs(coastal_dem_lon) + 1    
    
        #   transform to string
        coastal_dem_lat = str(coastal_dem_lat)
        coastal_dem_lon = str(coastal_dem_lon)
        if len(coastal_dem_lon) == 1:
            coastal_dem_lon = '0' + coastal_dem_lon
        if len(coastal_dem_lat) == 1:
            coastal_dem_lat = '00' + coastal_dem_lat
        elif len(coastal_dem_lat) == 2:
            coastal_dem_lat = '0' + coastal_dem_lat
    
        #   add north/south and east/west based on the coordinates
        if point_lon < 0:
            coastal_dem_lon = 'S' + coastal_dem_lon
        else:
            coastal_dem_lon = 'N' + coastal_dem_lon
        if point_lat < 0:
            coastal_dem_lat = 'W' + coastal_dem_lat
        else:
            coastal_dem_lat = 'E' + coastal_dem_lat

        if raster_in:
            coastal_val_cs_point = db_fc.read_raster_val(rb_coastal, gt_coastal, rx_coastal, ry_coastal, point_lat, point_lon)

        else:
            coastal_val_cs_point = np.nan
        
        return coastal_val_cs_point
    

    #   loop through different times
    for ts in time_steps:
        #print(ts)
        if ts < 0:
            suff = '_BP'
        elif ts > 0:
            suff = '_AP'

        gw_rch_path = os.path.join(gw_rch_dir, 'gw_rch_ts_' + str(abs(ts)) + suff + '.tif')

        raster_in_gw_rch = gdal.Open(gw_rch_path)
        ##  assign the geotransform, band and noval
        gt_gw_rch = raster_in_gw_rch.GetGeoTransform()
        rb_gw_rch = raster_in_gw_rch.GetRasterBand(1)
        rx_gw_rch = raster_in_gw_rch.RasterXSize
        ry_gw_rch = raster_in_gw_rch.RasterYSize

        #  determine how many points to create (25 means point per 1km) and the distance of the farthest point from point_5km
        max_dist = 200000
        n_points = int(max_dist / cs_pts_dist)
        x_coords = np.arange(-max_dist / 1000., max_dist / 1000. + cs_pts_dist / 1000., cs_pts_dist / 1000.)
        x_coords = [round(i, 2) for i in x_coords]
        #   define the number of points and distance from the cross-section point
        avg_dist = 500
        n_points_avg = int(avg_dist / cs_pts_dist)
        
        #gebco_avg_land_lst, gebco_avg_sea_lst = np.ones(len(x_coords)) * np.nan, np.ones(len(x_coords)) * np.nan
        #merit_avg_land_lst, merit_avg_sea_lst = np.ones(len(x_coords)) * np.nan, np.ones(len(x_coords)) * np.nan
        
        gw_rch_val_cs_point = get_gw_rch_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, point_lat, point_lon, raster_in_gw_rch)
        
        gw_rch_val_lst = np.ones(len(x_coords)) * np.nan
        gw_rch_val_avg_lst = np.ones(len(x_coords)) * np.nan
        
        gw_rch_val_lst[int(len(x_coords) / 2)] = gw_rch_val_cs_point
        gw_rch_val_avg_lst[int(len(x_coords) / 2)] = gw_rch_val_cs_point
        
        ##  run the loop for all the cross-section points
        for a in range(n_points):
            try:
                ##  calculate the b_length
                b_length = max_dist - a * cs_pts_dist
                c_coords_utm = db_fc.coord_from_triangle(point_5km_utm[0], point_5km_utm[1], closest_node_utm[0], closest_node_utm[1], b_length)
    
                ##  transform the coordinates back to wgs84, position of the equator is taken into account
                c_coords_wgs = db_fc.equator_position_angles(c_coords_utm[1], c_coords_utm[0], c_coords_utm[3], c_coords_utm[2], zone_num, zone_let)
                c_coord_wgs_x_land, c_coord_wgs_y_land, c_coord_wgs_x_sea, c_coord_wgs_y_sea = \
                c_coords_wgs[0][1], c_coords_wgs[0][0], c_coords_wgs[1][1], c_coords_wgs[1][0]

                ##  read values from the gw rch raster
                val_land = round(get_gw_rch_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, c_coord_wgs_x_land, c_coord_wgs_y_land, raster_in_gw_rch), 2)
                val_sea = round(get_gw_rch_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, c_coord_wgs_x_sea, c_coord_wgs_y_sea, raster_in_gw_rch), 2)
                
                ##  define sum and count of avg_values for final average value
                sum_avg_land, sum_avg_sea = [], []
                ##  calculate the average value for each cross-section point
                for b in range(n_points_avg):
                    avg_point_dist = avg_dist - b * cs_pts_dist
    
                    avg_point_coord_land = db_fc.avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[1], c_coords_utm[0], avg_point_dist)
                    avg_point_coord_sea = db_fc.avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[3], c_coords_utm[2], avg_point_dist)
    
                    ##  transform to wgs coordinates
                    avg_point_coord_wgs_land = db_fc.equator_position_angles(avg_point_coord_land[1], avg_point_coord_land[0], avg_point_coord_land[3], avg_point_coord_land[2], zone_num, zone_let)
                    avg_point_coord_wgs_sea = db_fc.equator_position_angles(avg_point_coord_sea[1], avg_point_coord_sea[0], avg_point_coord_sea[3], avg_point_coord_sea[2], zone_num, zone_let)
    
                    ##  GEBCO get the pixel value for the cross-section avg point for the land pixel
                    avg_pixel_val_land_1 = db_fc.read_raster_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, avg_point_coord_wgs_land[0][1], avg_point_coord_wgs_land[0][0])
                    sum_avg_land.append(avg_pixel_val_land_1)
                    avg_pixel_val_land_2 = db_fc.read_raster_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, avg_point_coord_wgs_land[1][1], avg_point_coord_wgs_land[1][0])
                    sum_avg_land.append(avg_pixel_val_land_2)
                    ##  get the pixel value for the cross-section avg point for the sea pixel
                    avg_pixel_val_sea_1 = db_fc.read_raster_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, avg_point_coord_wgs_sea[0][1], avg_point_coord_wgs_sea[0][0])
                    sum_avg_sea.append(avg_pixel_val_sea_1)
                    avg_pixel_val_sea_2 = db_fc.read_raster_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, avg_point_coord_wgs_sea[1][1], avg_point_coord_wgs_sea[1][0])
                    sum_avg_sea.append(avg_pixel_val_sea_2)
                
                #   get rid of potential non-value data in the MERIT list
                sum_avg_land = [i for i in sum_avg_land if i > -100.]
                sum_avg_sea = [i for i in sum_avg_sea if i > -100.]
                
                try:
                    avg_land = round(sum(sum_avg_land) / len(sum_avg_land), 2)
                except ZeroDivisionError:
                    avg_land = np.nan
                try: 
                    avg_sea = round(sum(sum_avg_sea) / len(sum_avg_sea), 2)
                except ZeroDivisionError:
                    avg_sea = np.nan
                        
                gw_rch_val_avg_lst[a]= avg_sea
                gw_rch_val_avg_lst[n_points + n_points - a] = avg_land
   
                gw_rch_val_lst[a]= val_sea
                gw_rch_val_lst[n_points + n_points - a] = val_land
    
            ##  ZeroDivisionError happens when point is located on top of the node (closest node)
            ##  this is the case in small islands where equidistant points are created..
            except ZeroDivisionError:
                continue
                print('Division by zero')    
    
        gw_rch_all_ts_lst.append(gw_rch_val_lst)
        gw_rch_avg_all_ts_lst.append(gw_rch_val_avg_lst)

    #   create a netcdf file that will be returned as output
    xa_sum = xr.Dataset(data_vars = {'gw_rch' : (('time', 'x'), gw_rch_all_ts_lst),
                                     'gw_rch_AVG' : (('time', 'x'), gw_rch_avg_all_ts_lst)},
                        coords = {'time' : time_steps,
                                  'x' : x_coords})

    xa_name = 'cs_' + str(id_cs) + '_pt_dist_' + str(int(cs_pts_dist)) + '_gw_rch.nc'
    xa_sum.to_netcdf(os.path.join(dir_out, xa_name))           


"""
function that reads in gw_rch for all time steps for each cross-section
"""
def get_gw_rch_clay(id_cs, cs_pts_dist, dir_out, gw_rch_dir):
    
    pts_info = df_points_5km.loc[df_points_5km['id_cs'] == id_cs]

    time_steps = np.arange(-30000, -1000, 1000).tolist()
    time_steps_2 = np.arange(-1900, 500, 100).tolist()
    time_steps = time_steps + time_steps_2 
    
    #if pts_info != []:
    point_lat, point_lon = pts_info.x.values.tolist()[0], pts_info.y.values.tolist()[0]

    #   get the closest node to the coastal point
    nearest_geoms = nearest_points(Point(point_lat, point_lon), nodes_geom)
    near_idx1 = nearest_geoms[1]
    closest_node_lat, closest_node_lon = near_idx1.x, near_idx1.y        
    
    ##  calculate the distance between the point_5m and closest node
    ##  to do this we convert the lat/lon coordinates to UTM coordinates (using the utm package)
    point_5km_utm = utm.from_latlon(point_lon, point_lat)
    closest_node_utm = utm.from_latlon(closest_node_lon, closest_node_lat)
    zone_num, zone_let = point_5km_utm[2], point_5km_utm[3]

    ##  check that zone_num is between 1 and 60, case for areas that are split between
    ##  the zones 1 and 60 (e.g. certain islands in Pacific..)
    if zone_num > 60:
        zone_num = zone_num - 60

    gw_rch_all_ts_lst, gw_rch_avg_all_ts_lst = [], []
    
    """
    point_lat = point_lat
    point_lon = point_lon
    """

    def get_gw_rch_val(rb_coastal, gt_coastal, rx_coastal, ry_coastal, point_lat, point_lon, raster_in):

        #   first decide on the numbers 
        coastal_dem_lat = math.floor(point_lat / 1.) * 1
        if point_lat < 0:
            coastal_dem_lat = abs(coastal_dem_lat) + 1
    
        coastal_dem_lon = math.floor(point_lon / 1.) * 1
        if point_lon < 0:
            coastal_dem_lon = abs(coastal_dem_lon) + 1    
    
        #   transform to string
        coastal_dem_lat = str(coastal_dem_lat)
        coastal_dem_lon = str(coastal_dem_lon)
        if len(coastal_dem_lon) == 1:
            coastal_dem_lon = '0' + coastal_dem_lon
        if len(coastal_dem_lat) == 1:
            coastal_dem_lat = '00' + coastal_dem_lat
        elif len(coastal_dem_lat) == 2:
            coastal_dem_lat = '0' + coastal_dem_lat
    
        #   add north/south and east/west based on the coordinates
        if point_lon < 0:
            coastal_dem_lon = 'S' + coastal_dem_lon
        else:
            coastal_dem_lon = 'N' + coastal_dem_lon
        if point_lat < 0:
            coastal_dem_lat = 'W' + coastal_dem_lat
        else:
            coastal_dem_lat = 'E' + coastal_dem_lat

        if raster_in:
            coastal_val_cs_point = db_fc.read_raster_val(rb_coastal, gt_coastal, rx_coastal, ry_coastal, point_lat, point_lon)

        else:
            coastal_val_cs_point = np.nan
        
        return coastal_val_cs_point
    

    #   loop through different times
    for ts in time_steps:
        #print(ts)
        if ts < 0:
            suff = '_BP'
        elif ts > 0:
            suff = '_AP'

        gw_rch_path = os.path.join(gw_rch_dir, 'gw_rch_ts_' + str(abs(ts)) + suff + '.tif')

        raster_in_gw_rch = gdal.Open(gw_rch_path)
        ##  assign the geotransform, band and noval
        gt_gw_rch = raster_in_gw_rch.GetGeoTransform()
        rb_gw_rch = raster_in_gw_rch.GetRasterBand(1)
        rx_gw_rch = raster_in_gw_rch.RasterXSize
        ry_gw_rch = raster_in_gw_rch.RasterYSize

        #  determine how many points to create (25 means point per 1km) and the distance of the farthest point from point_5km
        max_dist = 200000
        n_points = int(max_dist / cs_pts_dist)
        x_coords = np.arange(-max_dist / 1000., max_dist / 1000. + cs_pts_dist / 1000., cs_pts_dist / 1000.)
        x_coords = [round(i, 2) for i in x_coords]
        #   define the number of points and distance from the cross-section point
        avg_dist = 500
        n_points_avg = int(avg_dist / cs_pts_dist)
        
        #gebco_avg_land_lst, gebco_avg_sea_lst = np.ones(len(x_coords)) * np.nan, np.ones(len(x_coords)) * np.nan
        #merit_avg_land_lst, merit_avg_sea_lst = np.ones(len(x_coords)) * np.nan, np.ones(len(x_coords)) * np.nan
        
        gw_rch_val_cs_point = get_gw_rch_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, point_lat, point_lon, raster_in_gw_rch)
        
        gw_rch_val_lst = np.ones(len(x_coords)) * np.nan
        gw_rch_val_avg_lst = np.ones(len(x_coords)) * np.nan
        
        gw_rch_val_lst[int(len(x_coords) / 2)] = gw_rch_val_cs_point
        gw_rch_val_avg_lst[int(len(x_coords) / 2)] = gw_rch_val_cs_point
        
        ##  run the loop for all the cross-section points
        for a in range(n_points):
            try:
                ##  calculate the b_length
                b_length = max_dist - a * cs_pts_dist
                c_coords_utm = db_fc.coord_from_triangle(point_5km_utm[0], point_5km_utm[1], closest_node_utm[0], closest_node_utm[1], b_length)
    
                ##  transform the coordinates back to wgs84, position of the equator is taken into account
                c_coords_wgs = db_fc.equator_position_angles(c_coords_utm[1], c_coords_utm[0], c_coords_utm[3], c_coords_utm[2], zone_num, zone_let)
                c_coord_wgs_x_land, c_coord_wgs_y_land, c_coord_wgs_x_sea, c_coord_wgs_y_sea = \
                c_coords_wgs[0][1], c_coords_wgs[0][0], c_coords_wgs[1][1], c_coords_wgs[1][0]

                ##  read values from the gw rch raster
                val_land = round(get_gw_rch_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, c_coord_wgs_x_land, c_coord_wgs_y_land, raster_in_gw_rch), 2)
                val_sea = round(get_gw_rch_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, c_coord_wgs_x_sea, c_coord_wgs_y_sea, raster_in_gw_rch), 2)
                
                ##  define sum and count of avg_values for final average value
                sum_avg_land, sum_avg_sea = [], []
                ##  calculate the average value for each cross-section point
                for b in range(n_points_avg):
                    avg_point_dist = avg_dist - b * cs_pts_dist
    
                    avg_point_coord_land = db_fc.avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[1], c_coords_utm[0], avg_point_dist)
                    avg_point_coord_sea = db_fc.avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[3], c_coords_utm[2], avg_point_dist)
    
                    ##  transform to wgs coordinates
                    avg_point_coord_wgs_land = db_fc.equator_position_angles(avg_point_coord_land[1], avg_point_coord_land[0], avg_point_coord_land[3], avg_point_coord_land[2], zone_num, zone_let)
                    avg_point_coord_wgs_sea = db_fc.equator_position_angles(avg_point_coord_sea[1], avg_point_coord_sea[0], avg_point_coord_sea[3], avg_point_coord_sea[2], zone_num, zone_let)
    
                    ##  GEBCO get the pixel value for the cross-section avg point for the land pixel
                    avg_pixel_val_land_1 = db_fc.read_raster_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, avg_point_coord_wgs_land[0][1], avg_point_coord_wgs_land[0][0])
                    sum_avg_land.append(avg_pixel_val_land_1)
                    avg_pixel_val_land_2 = db_fc.read_raster_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, avg_point_coord_wgs_land[1][1], avg_point_coord_wgs_land[1][0])
                    sum_avg_land.append(avg_pixel_val_land_2)
                    ##  get the pixel value for the cross-section avg point for the sea pixel
                    avg_pixel_val_sea_1 = db_fc.read_raster_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, avg_point_coord_wgs_sea[0][1], avg_point_coord_wgs_sea[0][0])
                    sum_avg_sea.append(avg_pixel_val_sea_1)
                    avg_pixel_val_sea_2 = db_fc.read_raster_val(rb_gw_rch, gt_gw_rch, rx_gw_rch, ry_gw_rch, avg_point_coord_wgs_sea[1][1], avg_point_coord_wgs_sea[1][0])
                    sum_avg_sea.append(avg_pixel_val_sea_2)
                
                #   get rid of potential non-value data in the MERIT list
                sum_avg_land = [i for i in sum_avg_land if i > -100.]
                sum_avg_sea = [i for i in sum_avg_sea if i > -100.]
                
                try:
                    avg_land = round(sum(sum_avg_land) / len(sum_avg_land), 2)
                except ZeroDivisionError:
                    avg_land = np.nan
                try: 
                    avg_sea = round(sum(sum_avg_sea) / len(sum_avg_sea), 2)
                except ZeroDivisionError:
                    avg_sea = np.nan
                        
                gw_rch_val_avg_lst[a]= avg_sea
                gw_rch_val_avg_lst[n_points + n_points - a] = avg_land
   
                gw_rch_val_lst[a]= val_sea
                gw_rch_val_lst[n_points + n_points - a] = val_land
    
            ##  ZeroDivisionError happens when point is located on top of the node (closest node)
            ##  this is the case in small islands where equidistant points are created..
            except ZeroDivisionError:
                continue
                print('Division by zero')    
    
        gw_rch_all_ts_lst.append(gw_rch_val_lst)
        gw_rch_avg_all_ts_lst.append(gw_rch_val_avg_lst)

    #   create a netcdf file that will be returned as output
    xa_sum = xr.Dataset(data_vars = {'gw_rch' : (('time', 'x'), gw_rch_all_ts_lst),
                                     'gw_rch_AVG' : (('time', 'x'), gw_rch_avg_all_ts_lst)},
                        coords = {'time' : time_steps,
                                  'x' : x_coords})

    xa_name = 'cs_' + str(id_cs) + '_pt_dist_' + str(int(cs_pts_dist)) + '_gw_rch_clay.nc'
    xa_sum.to_netcdf(os.path.join(dir_out, xa_name))     



def plot_topo_comparison(lst_gebco, lst_gebco_avg, lst_merit, lst_merit_avg, lst_coastal, lst_coastal_avg, id_cs, cs_pts_dist, out_dir, cont_dir):

    max_dist = 200000
    x_coords = np.arange(-max_dist / 1000., max_dist / 1000. + cs_pts_dist / 1000., cs_pts_dist / 1000.)
    x_coords = [round(i, 2) for i in x_coords]
    
    #if cont_dir == 'left':
    #    x_coords = x_coords[::-1]
    
        #   also plot the topography in one plot
    fig, ax = plt.subplots(figsize = (15, 8))
    
    #   calculate the indexes so the zoom area goes from 10km offshore to 50km inland
    x_st_idx = int(200000. / (cs_pts_dist) - 50000. / cs_pts_dist)
    x_end_idx = int(200000. / (cs_pts_dist) + 10000. / cs_pts_dist)
    
    ax.plot(x_coords[x_st_idx : x_end_idx], lst_gebco[x_st_idx : x_end_idx], c = 'deepskyblue', label = 'GEBCO', lw = 1)
    ax.plot(x_coords[x_st_idx : x_end_idx], lst_gebco_avg[x_st_idx : x_end_idx], c = 'blue', label = 'GEBCO_avg', lw = 2)
    ax.plot(x_coords[x_st_idx : x_end_idx], lst_merit[x_st_idx : x_end_idx], c = 'gold', label = 'MERIT', lw = 1)
    ax.plot(x_coords[x_st_idx : x_end_idx], lst_merit_avg[x_st_idx : x_end_idx], c = 'darkorange', label = 'MERIT_avg', lw = 2)
    ax.plot(x_coords[x_st_idx : x_end_idx], lst_coastal[x_st_idx : x_end_idx], c = 'limegreen', label = 'CoastalDEM', lw = 1)
    ax.plot(x_coords[x_st_idx : x_end_idx], lst_coastal_avg[x_st_idx : x_end_idx], c = 'green', label = 'CoastalDEM_avg', lw = 2)
    
    ax.legend()
    ax.set(xlabel = 'Distance from coast (km)', ylabel = 'Elevation (m asl.)', title = 'DEM datasets comparison in coastal profile ID = ' + str(id_cs))
    ax.grid()
    plt_name_cst_zoom = 'cs_' + str(id_cs) + '_pt_dist_' + str(int(cs_pts_dist)) + '_ALL_DEM_COMPARISON.png'
    fig.savefig(os.path.join(out_dir, plt_name_cst_zoom), dpi = 300, facecolor = 'w', edgecolor = 'w', orientation = 'portrait', 
                papertype = None, format = None, transparent = False, bbox_inches = 'tight', pad_inches = 0.3, frameon = None)



if os.path.isfile(os.path.join(id_riv_bas_dir, '_input_lists.npy')):
    #   change status to 1, means finished succesfully
    #df_in = pd.read_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb.csv'))
    df_in = pd.read_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb_ADDED_SRMs.csv'))
    df_in.loc[df_in['id_loop'] == id_loop, 'finished'] = 1
    df_in.update(df_in)
    df_in.to_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb_ADDED_SRMs.csv'), sep = ',', encoding = 'utf-8', index = False)
    #df_in.to_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb.csv'), sep = ',', encoding = 'utf-8', index = False)
    del df_in        


print(status, reg_id, os.path.join(id_riv_bas_dir, '_input_lists.npy'))

if status != 1 and reg_id != -1 and not os.path.isfile(os.path.join(id_riv_bas_dir, '_input_lists.npy')):
    #   select the profiles and creat a list of id_cs
    id_cs_lst = cs_regs_raw.loc[(cs_regs_raw['coscat'] == coscat_id) & (cs_regs_raw['cst_reg_id'] == reg_id)]['id_cs'].values.tolist()
    
    print(id_cs_lst)
    
    csv_out_lst = []    
    #   define string and integer 
    #id_subreg_str = "ID_sub_reg = " + str(reg_id)
    #   add 0s in front of the name to have an ordered list of COSCAT ids in the folder 

    #   define the model dimensions and empty lists to be filled 
    del_col = 100.
    del_lay = 10.
    cs_points_dist = 500.
    
    lst_cs_width, lst_cs_thk = [], []
    lst_anchor_dist, lst_anchor_depth, lst_wtd, lst_topo, lst_soil_type, lst_soil_thk   = [], [], [], [], [], [] 
    lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch, lst_cs_p_min_et, lst_cs_k_soil = [], [], [], [], [], [] 
    lst_cs_drn_rate, lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay = [], [], []
      
    lst_topo_gebco_10m, lst_topo_gebco_avg_10m, lst_topo_merit_10m, lst_topo_merit_avg_10m, lst_gebco_merit_10m, lst_gebco_merit_avg_10m = [], [], [], [], [], []     
    lst_topo_gebco_25m, lst_topo_gebco_avg_25m, lst_topo_merit_25m, lst_topo_merit_avg_25m, lst_gebco_merit_25m, lst_gebco_merit_avg_25m = [], [], [], [], [], []    
    lst_topo_gebco_50m, lst_topo_gebco_avg_50m, lst_topo_merit_50m, lst_topo_merit_avg_50m, lst_gebco_merit_50m, lst_gebco_merit_avg_50m = [], [], [], [], [], []
    lst_topo_gebco_100m, lst_topo_gebco_avg_100m, lst_topo_merit_100m, lst_topo_merit_avg_100m, lst_gebco_merit_100m, lst_gebco_merit_avg_100m = [], [], [], [], [], [] 
    lst_topo_gebco_250m, lst_topo_gebco_avg_250m, lst_topo_merit_250m, lst_topo_merit_avg_250m, lst_gebco_merit_250m, lst_gebco_merit_avg_250m = [], [], [], [], [], [] 
    lst_topo_gebco_500m, lst_topo_gebco_avg_500m, lst_topo_merit_500m, lst_topo_merit_avg_500m, lst_gebco_merit_500m, lst_gebco_merit_avg_500m = [], [], [], [], [], [] 
    lst_topo_cstdem_10m, lst_topo_cstdem_25m, lst_topo_cstdem_50m, lst_topo_cstdem_100m, lst_topo_cstdem_250m, lst_topo_cstdem_500m = [], [], [], [], [], []
    lst_topo_cstdem_avg_10m, lst_topo_cstdem_avg_25m, lst_topo_cstdem_avg_50m, lst_topo_cstdem_avg_100m, lst_topo_cstdem_avg_250m, lst_topo_cstdem_avg_500m = [], [], [], [], [], []
    lst_topo_cstdem_gebco_10m, lst_topo_cstdem_gebco_25m, lst_topo_cstdem_gebco_50m, lst_topo_cstdem_gebco_100m, lst_topo_cstdem_gebco_250m, lst_topo_cstdem_gebco_500m = [], [], [], [], [], []
    lst_topo_cstdem_gebco_avg_10m, lst_topo_cstdem_gebco_avg_25m, lst_topo_cstdem_gebco_avg_50m, lst_topo_cstdem_gebco_avg_100m, lst_topo_cstdem_gebco_avg_250m, lst_topo_cstdem_gebco_avg_500m = [], [], [], [], [], []
       
    #   set the counters of different coastal types
    pts_count = 0
    fos = -999
    ls_elev = -999
    ls_elev_50 = -999   
    
    #   loop through the cross-sections and select required info from all the necessary tables
    for id_num in id_cs_lst:
        id_cs = int(id_num)
         
        #   info about coastal plain width and the anchor point
        cs_info = df_sed_thk_v2.loc[df_sed_thk_v2['id_cs'] == id_cs]
        cs_class = df_cst_type.loc[df_cst_type['id_cs'] == id_cs]
        coastal_class = cs_class.class_fin.values[0]
        
        #   create the groundwater recharge netcdf file
        if not os.path.isfile(os.path.join(id_riv_bas_dir, 'cs_' + str(id_cs) + '_pt_dist_100_gw_rch.nc')):
            get_gw_rch(id_cs, 100., id_riv_bas_dir, gw_rch_dir)  
        #   create the groundwater recharge netcdf file
        if not os.path.isfile(os.path.join(id_riv_bas_dir, 'cs_' + str(id_cs) + '_pt_dist_100_gw_rch_clay.nc')):
            get_gw_rch_clay(id_cs, 100., id_riv_bas_dir, gw_rch_clay_dir)  
            
        pts_info = df_points_5km.loc[df_points_5km['id_cs'] == id_cs]
        point_lat, point_lon = pts_info.x.values.tolist()[0], pts_info.y.values.tolist()[0]
        #   get the closest node to the coastal point
        nearest_geoms = nearest_points(Point(point_lat, point_lon), nodes_geom)
        near_idx1 = nearest_geoms[1]
        closest_node_lat, closest_node_lon = near_idx1.x, near_idx1.y   
        
        try:
            cs_plain_width = cs_info.cst_plain_width.values[0]
            anchor_dist = cs_info.nasa_point_dist.values[0]
            anchor_depth = cs_info.nasa_point_depth.values[0]
            avg_depth = cs_info.overall_avg.values[0]        
        except IndexError:
            continue
        
        #   if it is delta then and the anchor depth is more than 200km then assign 200km
        #if anchor_dist < 0.:
        #    anchor_dist = 200.0 + anchor_dist
        if 200. - anchor_dist < 0:
            anchor_dist = 200.
        else:
            anchor_dist = 200. - anchor_dist
            
        cs_model_point = ws_ms.point(None, None, 400, 200000, id_cs = id_cs)
        
        #   happens when there are no data extracted for the id_cs
        try:
            cs_model_point.get_input_data(adjust_to_ocean = True)    
            cs_wtd = cs_model_point.dta_wtd_depth 
            cs_topo = cs_model_point.dta_topo
            cs_soil_type = cs_model_point.dta_soil_type
            cs_soil_thk = cs_model_point.dta_soil_thk
            cs_offshore = cs_model_point.dta_offshore_sed
            cs_offshore_flt = [float(i) for i in cs_offshore]
            cs_nasa = cs_model_point.dta_thk_pel
            cs_nasa_flt = [float(i) for i in cs_nasa]    
            #cs_pcr_rch = cs_model_point.dta_pcr_rch
            #cs_watergap_rch = cs_model_point.dta_watergap_rch
            #cs_p_min_et = cs_model_point.dta_p_min_et
            cs_k_soil = cs_model_point.dta_k_soil
            cs_drn_rate = cs_model_point.dta_drn_rate
            cs_glhymps_top_lay = cs_model_point.dta_glhymps_top_lay                   
            cs_glhymps_bot_lay = cs_model_point.dta_glhymps_bot_lay  
    
            lst_cs_thk.append(avg_depth)
            lst_cs_width.append(cs_plain_width)   
            lst_anchor_dist.append(anchor_dist)
            lst_anchor_depth.append(anchor_depth)
            lst_wtd.append(cs_wtd)
            lst_topo.append(cs_topo)
            lst_soil_type.append(cs_soil_type)
            lst_soil_thk.append(cs_soil_thk)
            lst_cs_nasa.append(cs_nasa_flt)
            lst_cs_offshore.append(cs_offshore_flt)                
            #lst_cs_pcr_rch.append(cs_pcr_rch)
            #lst_cs_watergap_rch.append(cs_watergap_rch)
            #lst_cs_p_min_et.append(cs_p_min_et)
            lst_cs_k_soil.append(cs_k_soil)
            lst_cs_drn_rate.append(cs_drn_rate)  
            lst_cs_glhymps_top_lay.append(cs_glhymps_top_lay)
            lst_cs_glhymps_bot_lay.append(cs_glhymps_bot_lay)                
    
            cont_direction = cs_model_point.land_pos
            
            pts_count += 1             
            
        except IndexError:
            
            try:
                pts_count += 1    
                #continue
                df_cs_topo = pd.read_csv(os.path.join(input_files_dir, '_csv_input_files', 'cs_gebco_2014_avg_DELTAS.csv'))
                cs_topo = df_cs_topo.loc[df_cs_topo['id_cs'] == id_cs].values.tolist()[0][400:1201]
                del df_cs_topo
                df_cs_nasa = pd.read_csv(os.path.join(input_files_dir,'_csv_input_files', 'cs_nasa_thick_avg_DELTAS.csv'))
                cs_nasa = df_cs_nasa.loc[df_cs_nasa['id_cs'] == id_cs].values.tolist()[0][400:1201]
                del df_cs_nasa
                df_cs_pcr = pd.read_csv(os.path.join(input_files_dir,'_csv_input_files', 'cs_aq_thick_DELTAS.csv'))
                cs_pcr = df_cs_pcr.loc[df_cs_pcr['id_cs'] == id_cs].values.tolist()[0][400:1201]
                del df_cs_pcr
                df_cs_glim = pd.read_csv(os.path.join(input_files_dir,'_csv_input_files', 'cs_glim_litho_DELTAS.csv'))
                cs_glim = df_cs_glim.loc[df_cs_glim['id_cs'] == id_cs].values.tolist()[0][400:1201]
                cs_nasa_flt = cs_nasa
                del df_cs_glim            
                
                cs_offshore = rivbas.read_cs_raster_data_no_sql(offshore_sed_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)
                cs_offshore_flt  = [float(i) if i is not None else -32768.0 for i in cs_offshore]  #[float(i) for i in cs_offshore]  
                cs_wtd = rivbas.read_cs_raster_data_no_sql(wtd_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_soil_thk = rivbas.read_cs_raster_data_no_sql(soil_thk_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_soil_type = rivbas.read_cs_raster_data_no_sql(soil_type_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'int', avg = False)[1:]
                cs_glhymps_top_lay = rivbas.read_cs_raster_data_no_sql(glhymps_top_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_glhymps_bot_lay = rivbas.read_cs_raster_data_no_sql(glhymps_bot_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_k_soil = rivbas.read_cs_raster_data_no_sql(k_soil_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_drn_rate = rivbas.read_cs_raster_data_no_sql(drn_rate_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
    
                """
                cs_offshore = rivbas.read_cs_raster_data(offshore_sed_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_offshore_sed', num_type = 'int', avg = False)[1:]
                cs_offshore_flt  = [float(i) if i is not None else -32768.0 for i in cs_offshore]  #[float(i) for i in cs_offshore]  
                cs_wtd = rivbas.read_cs_raster_data(wtd_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_wtd_depth', num_type = 'real', avg = False)[1:]
                #cs_pcr_rch = rivbas.read_cs_raster_data(pcr_rch_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_rch_pcr', num_type = 'real', avg = False)[1:]
                #cs_p_min_et = rivbas.read_cs_raster_data(p_min_et_rch_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_p_min_et', num_type = 'real', avg = False)[1:]
                cs_soil_thk = rivbas.read_cs_raster_data(soil_thk_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_soil_thick', num_type = 'real', avg = False)[1:]
                cs_soil_type = rivbas.read_cs_raster_data(soil_type_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_soil_type', num_type = 'int', avg = False)[1:]
                #cs_watergap_rch = rivbas.read_cs_raster_data(watergap_rch_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_rch_watergap', num_type = 'real', avg = False)[1:]
                cs_glhymps_top_lay = rivbas.read_cs_raster_data(glhymps_top_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_glhymps_layer_1', num_type = 'real', avg = False)[1:]
                cs_glhymps_bot_lay = rivbas.read_cs_raster_data(glhymps_bot_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_glhymps_layer_2', num_type = 'real', avg = False)[1:]
                cs_k_soil = rivbas.read_cs_raster_data(k_soil_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_k_soil', num_type = 'real', avg = False)[1:]
                cs_drn_rate = rivbas.read_cs_raster_data(drn_rate_dir, id_cs, 400, 200000, 2, 2500, db_name, db_name_world, 'cs_drn_rate', num_type = 'real', avg = False)[1:]
                """
                
                n_pts = 400
                cs_point_gebco = cs_topo[n_pts]
                if cs_point_gebco < 0:
                    ##  loop through the neighbouring points and check if they are equal to or larger than 0
                    ##  also check in which direction from the mid-list value the coastal point actually is,
                    ##  assign that to cont_direction (continental) - either left or right (east or west)
                    for c in range(len(cs_topo) // 2):
                        next_val_left = cs_topo[n_pts - c]
                        next_val_right = cs_topo[n_pts + c]
                        if next_val_left >= 0:
                            coastal_point_indx = n_pts - c
                            cont_direction = 'left'
                            break
                        elif next_val_right >= 0:
                            coastal_point_indx = n_pts + c
                            cont_direction = 'right'
                            break
                        else:
                            pass
                ##  if the coastal point itself is = 0 assign it
                elif cs_point_gebco == 0:
                    coastal_point_indx = n_pts
                    ##  check which direction is the continental part - 5 pixels from the coastal point
                    next_val_left = cs_topo[coastal_point_indx - 25 : coastal_point_indx]
                    next_val_right = cs_topo[coastal_point_indx : coastal_point_indx + 25]
                    if float(sum(next_val_left) / 25.) < 0:
                        cont_direction = 'right'
                    elif float(sum(next_val_right) / 25.) < 0:
                        cont_direction = 'left'
                    else:
                        if float(sum(next_val_left) / 25.) > float(sum(next_val_right) / 25.):
                            cont_direction = 'left'
                        else:
                            cont_direction = 'right'
            
                else:
                    ##  loop through the values to see which one is closest to 0 but still positive
                    for z in range(len(cs_topo) // 2):
                        next_val_left = cs_topo[n_pts - z]
                        next_val_right = cs_topo[n_pts + z]
                        if next_val_left < 0:
                            coastal_point_indx = n_pts - z
                            cont_direction = 'right'
                            break
                        elif next_val_right < 0:
                            coastal_point_indx = n_pts + z
                            cont_direction = 'left'
                            break
                        else:
                            pass
        
                #   check the position of the land, use the function from the DTBSE_sed_thick_functions
                #   transform to position of the ocean
                if cont_direction == 'right':
                    #   reverse the list of values so it seems like ocean is on the right hand side
                    cs_topo = list(reversed(cs_topo))
                    cs_nasa = list(reversed(cs_nasa))
                    cs_pcr = list(reversed(cs_pcr))
                    cs_glim = list(reversed(cs_glim))
                    cs_nasa_flt = list(reversed(cs_nasa_flt))
                    cs_offshore = list(reversed(cs_offshore))
                    cs_offshore_flt = list(reversed(cs_offshore_flt))
                    cs_wtd = list(reversed(cs_wtd))
                    #cs_pcr_rch = list(reversed(cs_pcr_rch))
                    #cs_p_min_et = list(reversed(cs_p_min_et))
                    cs_soil_thk = list(reversed(cs_soil_thk))
                    cs_soil_type = list(reversed(cs_soil_type))
                    #cs_watergap_rch = list(reversed(cs_watergap_rch))
                    cs_glhymps_top_lay = list(reversed(cs_glhymps_top_lay)) 
                    cs_glhymps_bot_lay = list(reversed(cs_glhymps_bot_lay))
                    cs_k_soil = list(reversed(cs_k_soil))
                    cs_drn_rate = list(reversed(cs_drn_rate))
            
                else:
                    pass     
        
                lst_cs_thk.append(avg_depth)
                lst_cs_width.append(cs_plain_width)
                lst_anchor_dist.append(anchor_dist)
                lst_anchor_depth.append(anchor_depth)     
                lst_cs_nasa.append(cs_nasa_flt)
                lst_wtd.append(cs_wtd)
                lst_topo.append(cs_topo)
                lst_soil_type.append(cs_soil_type)
                lst_soil_thk.append(cs_soil_thk)
                lst_cs_nasa.append(cs_nasa_flt)
                lst_cs_offshore.append(cs_offshore_flt)                
                #lst_cs_pcr_rch.append(cs_pcr_rch)
                #lst_cs_watergap_rch.append(cs_watergap_rch)
                #lst_cs_p_min_et.append(cs_p_min_et)
                lst_cs_k_soil.append(cs_k_soil)
                lst_cs_drn_rate.append(cs_drn_rate)  
                lst_cs_glhymps_top_lay.append(cs_glhymps_top_lay)
                lst_cs_glhymps_bot_lay.append(cs_glhymps_bot_lay)
                
            except IndexError:
                continue

        x_cnt, y_cnt = 43200, 21600
        x_st, y_st = -180, 90 #  lower left corner
        px_x, px_y = abs((2 * x_st) / x_cnt), abs((2 * y_st) / y_cnt) * -1
    
        #   create the topography lists and files
        #   but first check if the csv files already exist, in that case just load the lists
        if os.path.isfile(os.path.join(id_riv_bas_dir, 'cs_' + str(id_cs) + '_pt_dist_100_topo_merit.csv')):
            topo_100m_csv = pd.read_csv(os.path.join(id_riv_bas_dir, 'cs_' + str(id_cs) + '_pt_dist_100_topo_merit.csv'))
            if cont_direction == 'left':
                topo_100m_csv = topo_100m_csv.sort_values(by=['dist_from_cst'], axis = 0, ascending=False, inplace=False, kind='quicksort', na_position='last')
            else:
                topo_100m_csv = topo_100m_csv.sort_values(by=['dist_from_cst'], axis = 0, ascending=True, inplace=False, kind='quicksort', na_position='last')
            topo_100m = []
            topo_100m.append(topo_100m_csv['gebco'].values.tolist())
            topo_100m.append(topo_100m_csv['gebco_avg'].values.tolist())
            topo_100m.append(topo_100m_csv['merit'].values.tolist())
            topo_100m.append(topo_100m_csv['merit_avg'].values.tolist())            
            topo_100m.append(topo_100m_csv['merit_gebco'].values.tolist())
            topo_100m.append(topo_100m_csv['merit_gebco_avg'].values.tolist())                        
        else:
            topo_100m = merit_dem_topo(id_cs, 100., id_riv_bas_dir, cont_direction)
            
        if os.path.isfile(os.path.join(id_riv_bas_dir, 'cs_' + str(id_cs) + '_pt_dist_100_topo_csdem.csv')):
            topo_coastaldem_100m_csv = pd.read_csv(os.path.join(id_riv_bas_dir, 'cs_' + str(id_cs) + '_pt_dist_100_topo_csdem.csv'))
            if cont_direction == 'left':
                topo_coastaldem_100m_csv = topo_coastaldem_100m_csv.sort_values(by=['dist_from_cst'], axis = 0, ascending=False, inplace=False, kind='quicksort', na_position='last')
            else:
                topo_coastaldem_100m_csv = topo_coastaldem_100m_csv.sort_values(by=['dist_from_cst'], axis = 0, ascending=True, inplace=False, kind='quicksort', na_position='last')
            topo_coastaldem_100m = []
            topo_coastaldem_100m.append(topo_coastaldem_100m_csv['gebco'].values.tolist())
            topo_coastaldem_100m.append(topo_coastaldem_100m_csv['gebco_avg'].values.tolist())
            topo_coastaldem_100m.append(topo_coastaldem_100m_csv['merit'].values.tolist())
            topo_coastaldem_100m.append(topo_coastaldem_100m_csv['merit_avg'].values.tolist())            
            topo_coastaldem_100m.append(topo_coastaldem_100m_csv['merit_gebco'].values.tolist())
            topo_coastaldem_100m.append(topo_coastaldem_100m_csv['merit_gebco_avg'].values.tolist())                        
        else:
            topo_coastaldem_100m = coastal_dem_topo(id_cs, 100., id_riv_bas_dir, cont_direction)            
            
        if topo_100m is not None and topo_coastaldem_100m is not None:
    
            gebco_lst_100m = topo_100m[0]
            gebco_avg_lst_100m = topo_100m[1]
            merit_lst_100m = topo_100m[2]
            merit_avg_lst_100m = topo_100m[3]
            merit_gebco_lst_100m = topo_100m[4]
            merit_gebco_avg_lst_100m = topo_100m[5]
            coastaldem_lst_100m = topo_coastaldem_100m[2]
            coastaldem_avg_lst_100m = topo_coastaldem_100m[3]
            coastaldem_gebco_lst_100m = topo_coastaldem_100m[4]
            coastaldem_gebco_avg_lst_100m = topo_coastaldem_100m[5]            
        
            gebco_lst_500m = gebco_lst_100m[::5]
            gebco_avg_lst_500m = gebco_avg_lst_100m[::5]
            merit_lst_500m = merit_lst_100m[::5]
            merit_avg_lst_500m = merit_avg_lst_100m[::5]
            merit_gebco_lst_500m = merit_gebco_lst_100m[::5]
            merit_gebco_avg_lst_500m = merit_gebco_avg_lst_100m[::5]
            coastaldem_lst_500m = coastaldem_lst_100m[::5]
            coastaldem_avg_lst_500m = coastaldem_avg_lst_100m[::5]
            coastaldem_gebco_lst_500m = coastaldem_gebco_lst_100m[::5]
            coastaldem_gebco_avg_lst_500m = coastaldem_gebco_avg_lst_100m[::5]                        

            x_lst_250 = np.arange(-200., 200 + 500. / 1000., 500. / 1000.)
            x_refined_lst_250 = np.arange(-200., 200 + 250. / 1000., 250. / 1000.)

            gebco_lst_250m = np.interp(x_refined_lst_250, x_lst_250, gebco_lst_500m).tolist()
            gebco_avg_lst_250m = np.interp(x_refined_lst_250, x_lst_250, gebco_avg_lst_500m).tolist()
            merit_lst_250m = np.interp(x_refined_lst_250, x_lst_250, merit_lst_500m).tolist()
            merit_avg_lst_250m = np.interp(x_refined_lst_250, x_lst_250, merit_avg_lst_500m).tolist()
            merit_gebco_lst_250m = np.interp(x_refined_lst_250, x_lst_250, merit_gebco_lst_500m).tolist()
            merit_gebco_avg_lst_250m = np.interp(x_refined_lst_250, x_lst_250, merit_gebco_avg_lst_500m).tolist()
            coastaldem_lst_250m = np.interp(x_refined_lst_250, x_lst_250, coastaldem_lst_500m).tolist()
            coastaldem_avg_lst_250m = np.interp(x_refined_lst_250, x_lst_250, coastaldem_avg_lst_500m).tolist()
            coastaldem_gebco_lst_250m = np.interp(x_refined_lst_250, x_lst_250, coastaldem_gebco_lst_500m).tolist()
            coastaldem_gebco_avg_lst_250m = np.interp(x_refined_lst_250, x_lst_250, coastaldem_gebco_avg_lst_500m).tolist()      
        
            x_lst_50 = np.arange(-200., 200 + 100. / 1000., 100. / 1000.)
            x_refined_lst_50 = np.arange(-200., 200 + 50. / 1000., 50. / 1000.)
        
            gebco_lst_50m = np.interp(x_refined_lst_50, x_lst_50, gebco_lst_100m).tolist()
            gebco_avg_lst_50m = np.interp(x_refined_lst_50, x_lst_50, gebco_avg_lst_100m).tolist()
            merit_lst_50m = np.interp(x_refined_lst_50, x_lst_50, merit_lst_100m).tolist()
            merit_avg_lst_50m = np.interp(x_refined_lst_50, x_lst_50, merit_avg_lst_100m).tolist()
            merit_gebco_lst_50m = np.interp(x_refined_lst_50, x_lst_50, merit_gebco_lst_100m).tolist()
            merit_gebco_avg_lst_50m = np.interp(x_refined_lst_50, x_lst_50, merit_gebco_avg_lst_100m).tolist()
            coastaldem_lst_50m = np.interp(x_refined_lst_50, x_lst_50, coastaldem_lst_100m).tolist()
            coastaldem_avg_lst_50m = np.interp(x_refined_lst_50, x_lst_50, coastaldem_avg_lst_100m).tolist()
            coastaldem_gebco_lst_50m = np.interp(x_refined_lst_50, x_lst_50, coastaldem_gebco_lst_100m).tolist()
            coastaldem_gebco_avg_lst_50m = np.interp(x_refined_lst_50, x_lst_50, coastaldem_gebco_avg_lst_100m).tolist()     
        
            x_lst_25 = np.arange(-200., 200 + 100. / 1000., 100. / 1000.)
            x_refined_lst_25 = np.arange(-200., 200 + 25. / 1000., 25. / 1000.)
            
            gebco_lst_25m = np.interp(x_refined_lst_25, x_lst_25, gebco_lst_100m).tolist()
            gebco_avg_lst_25m = np.interp(x_refined_lst_25, x_lst_25, gebco_avg_lst_100m).tolist()
            merit_lst_25m = np.interp(x_refined_lst_25, x_lst_25, merit_lst_100m).tolist()
            merit_avg_lst_25m = np.interp(x_refined_lst_25, x_lst_25, merit_avg_lst_100m).tolist()
            merit_gebco_lst_25m = np.interp(x_refined_lst_25, x_lst_25, merit_gebco_lst_100m).tolist()
            merit_gebco_avg_lst_25m = np.interp(x_refined_lst_25, x_lst_25, merit_gebco_avg_lst_100m).tolist()
            coastaldem_lst_25m = np.interp(x_refined_lst_25, x_lst_25, coastaldem_lst_100m).tolist()
            coastaldem_avg_lst_25m = np.interp(x_refined_lst_25, x_lst_25, coastaldem_avg_lst_100m).tolist()
            coastaldem_gebco_lst_25m = np.interp(x_refined_lst_25, x_lst_25, coastaldem_gebco_lst_100m).tolist()
            coastaldem_gebco_avg_lst_25m = np.interp(x_refined_lst_25, x_lst_25, coastaldem_gebco_avg_lst_100m).tolist()     
            
            x_lst_10 = np.arange(-200., 200 + 100. / 1000., 100. / 1000.)
            x_refined_lst_10 = np.arange(-200., 200 + 10. / 1000., 10. / 1000.)
            
            gebco_lst_10m = np.interp(x_refined_lst_10, x_lst_10, gebco_lst_100m).tolist()
            gebco_avg_lst_10m = np.interp(x_refined_lst_10, x_lst_10, gebco_avg_lst_100m).tolist()
            merit_lst_10m = np.interp(x_refined_lst_10, x_lst_10, merit_lst_100m).tolist()
            merit_avg_lst_10m = np.interp(x_refined_lst_10, x_lst_10, merit_avg_lst_100m).tolist()
            merit_gebco_lst_10m = np.interp(x_refined_lst_10, x_lst_10, merit_gebco_lst_100m).tolist()
            merit_gebco_avg_lst_10m = np.interp(x_refined_lst_10, x_lst_10, merit_gebco_avg_lst_100m).tolist()    
            coastaldem_lst_10m = np.interp(x_refined_lst_10, x_lst_10, coastaldem_lst_100m).tolist()
            coastaldem_avg_lst_10m = np.interp(x_refined_lst_10, x_lst_10, coastaldem_avg_lst_100m).tolist()
            coastaldem_gebco_lst_10m = np.interp(x_refined_lst_10, x_lst_10, coastaldem_gebco_lst_100m).tolist()
            coastaldem_gebco_avg_lst_10m = np.interp(x_refined_lst_10, x_lst_10, coastaldem_gebco_avg_lst_100m).tolist()     
            
            if cont_direction == 'right':
                gebco_lst_10m = list(reversed(gebco_lst_10m))
                gebco_avg_lst_10m = list(reversed(gebco_avg_lst_10m))
                merit_lst_10m = list(reversed(merit_lst_10m))
                merit_avg_lst_10m = list(reversed(merit_avg_lst_10m))
                merit_gebco_lst_10m = list(reversed(merit_gebco_lst_10m))
                merit_gebco_avg_lst_10m = list(reversed(merit_gebco_avg_lst_10m))
                coastaldem_lst_10m = list(reversed(coastaldem_lst_10m))
                coastaldem_avg_lst_10m = list(reversed(coastaldem_avg_lst_10m))
                coastaldem_gebco_lst_10m = list(reversed(coastaldem_gebco_lst_10m))
                coastaldem_gebco_avg_lst_10m = list(reversed(coastaldem_gebco_avg_lst_10m))   
                
                gebco_lst_25m = list(reversed(gebco_lst_25m))
                gebco_avg_lst_25m = list(reversed(gebco_avg_lst_25m))
                merit_lst_25m = list(reversed(merit_lst_25m))
                merit_avg_lst_25m = list(reversed(merit_avg_lst_25m))
                merit_gebco_lst_25m = list(reversed(merit_gebco_lst_25m))
                merit_gebco_avg_lst_25m = list(reversed(merit_gebco_avg_lst_25m))
                coastaldem_lst_25m = list(reversed(coastaldem_lst_25m))
                coastaldem_avg_lst_25m = list(reversed(coastaldem_avg_lst_25m))
                coastaldem_gebco_lst_25m = list(reversed(coastaldem_gebco_lst_25m))
                coastaldem_gebco_avg_lst_25m = list(reversed(coastaldem_gebco_avg_lst_25m))   
                
                gebco_lst_50m = list(reversed(gebco_lst_50m))
                gebco_avg_lst_50m = list(reversed(gebco_avg_lst_50m))
                merit_lst_50m = list(reversed(merit_lst_50m))
                merit_avg_lst_50m = list(reversed(merit_avg_lst_50m))
                merit_gebco_lst_50m = list(reversed(merit_gebco_lst_50m))
                merit_gebco_avg_lst_50m = list(reversed(merit_gebco_avg_lst_50m))
                coastaldem_lst_50m = list(reversed(coastaldem_lst_50m))
                coastaldem_avg_lst_50m = list(reversed(coastaldem_avg_lst_50m))
                coastaldem_gebco_lst_50m = list(reversed(coastaldem_gebco_lst_50m))
                coastaldem_gebco_avg_lst_50m = list(reversed(coastaldem_gebco_avg_lst_50m))   
                
                gebco_lst_100m = list(reversed(gebco_lst_100m))
                gebco_avg_lst_100m = list(reversed(gebco_avg_lst_100m))
                merit_lst_100m = list(reversed(merit_lst_100m))
                merit_avg_lst_100m = list(reversed(merit_avg_lst_100m))
                merit_gebco_lst_100m = list(reversed(merit_gebco_lst_100m))
                merit_gebco_avg_lst_100m = list(reversed(merit_gebco_avg_lst_100m))
                coastaldem_lst_100m = list(reversed(coastaldem_lst_100m))
                coastaldem_avg_lst_100m = list(reversed(coastaldem_avg_lst_100m))
                coastaldem_gebco_lst_100m = list(reversed(coastaldem_gebco_lst_100m))
                coastaldem_gebco_avg_lst_100m = list(reversed(coastaldem_gebco_avg_lst_100m))   

                gebco_lst_250m = list(reversed(gebco_lst_250m))
                gebco_avg_lst_250m = list(reversed(gebco_avg_lst_250m))
                merit_lst_250m = list(reversed(merit_lst_250m))
                merit_avg_lst_250m = list(reversed(merit_avg_lst_250m))
                merit_gebco_lst_250m = list(reversed(merit_gebco_lst_250m))
                merit_gebco_avg_lst_250m = list(reversed(merit_gebco_avg_lst_250m))
                coastaldem_lst_250m = list(reversed(coastaldem_lst_250m))
                coastaldem_avg_lst_250m = list(reversed(coastaldem_avg_lst_250m))
                coastaldem_gebco_lst_250m = list(reversed(coastaldem_gebco_lst_250m))
                coastaldem_gebco_avg_lst_250m = list(reversed(coastaldem_gebco_avg_lst_250m))  
                
                gebco_lst_500m = list(reversed(gebco_lst_500m))
                gebco_avg_lst_500m = list(reversed(gebco_avg_lst_500m))
                merit_lst_500m = list(reversed(merit_lst_500m))
                merit_avg_lst_500m = list(reversed(merit_avg_lst_500m))
                merit_gebco_lst_500m = list(reversed(merit_gebco_lst_500m))
                merit_gebco_avg_lst_500m = list(reversed(merit_gebco_avg_lst_500m))
                coastaldem_lst_500m = list(reversed(coastaldem_lst_500m))
                coastaldem_avg_lst_500m = list(reversed(coastaldem_avg_lst_500m))
                coastaldem_gebco_lst_500m = list(reversed(coastaldem_gebco_lst_500m))
                coastaldem_gebco_avg_lst_500m = list(reversed(coastaldem_gebco_avg_lst_500m))              

            #   plot the DEM comparison
            plot_topo_comparison(gebco_lst_100m, gebco_avg_lst_100m, merit_lst_100m, merit_avg_lst_100m, coastaldem_lst_100m,\
                                 coastaldem_avg_lst_100m, id_cs, 100., id_riv_bas_dir, cont_direction)
            
            #   last step, make sure that there are no non-data values in the merit_gebco lists ( = -9999.), also round values to 2 decimals
            merit_gebco_lst_500m = [round(i, 2) for i in merit_gebco_lst_500m]
            merit_gebco_lst_250m = [round(i, 2) for i in merit_gebco_lst_250m]
            merit_gebco_lst_100m = [round(i, 2) for i in merit_gebco_lst_100m]
            merit_gebco_lst_50m = [round(i, 2) for i in merit_gebco_lst_50m]
            merit_gebco_lst_25m = [round(i, 2) for i in merit_gebco_lst_25m]
            merit_gebco_lst_10m = [round(i, 2) for i in merit_gebco_lst_10m]
            
            merit_gebco_avg_lst_500m = [round(i, 2) for i in merit_gebco_avg_lst_500m]
            merit_gebco_avg_lst_250m = [round(i, 2) for i in merit_gebco_avg_lst_250m]
            merit_gebco_avg_lst_100m = [round(i, 2) for i in merit_gebco_avg_lst_100m]
            merit_gebco_avg_lst_50m = [round(i, 2) for i in merit_gebco_avg_lst_50m]
            merit_gebco_avg_lst_25m = [round(i, 2) for i in merit_gebco_avg_lst_25m]
            merit_gebco_avg_lst_10m = [round(i, 2) for i in merit_gebco_avg_lst_10m]
            
            coastaldem_gebco_lst_500m = [round(i, 2) for i in coastaldem_gebco_lst_500m]
            coastaldem_gebco_lst_250m = [round(i, 2) for i in coastaldem_gebco_lst_250m]
            coastaldem_gebco_lst_100m = [round(i, 2) for i in coastaldem_gebco_lst_100m]
            coastaldem_gebco_lst_50m = [round(i, 2) for i in coastaldem_gebco_lst_50m]
            coastaldem_gebco_lst_25m = [round(i, 2) for i in coastaldem_gebco_lst_25m]
            coastaldem_gebco_lst_10m = [round(i, 2) for i in coastaldem_gebco_lst_10m]
            
            coastaldem_gebco_avg_lst_500m = [round(i, 2) for i in coastaldem_gebco_avg_lst_500m]
            coastaldem_gebco_avg_lst_250m = [round(i, 2) for i in coastaldem_gebco_avg_lst_250m]
            coastaldem_gebco_avg_lst_100m = [round(i, 2) for i in coastaldem_gebco_avg_lst_100m]
            coastaldem_gebco_avg_lst_50m = [round(i, 2) for i in coastaldem_gebco_avg_lst_50m]
            coastaldem_gebco_avg_lst_25m = [round(i, 2) for i in coastaldem_gebco_avg_lst_25m]
            coastaldem_gebco_avg_lst_10m = [round(i, 2) for i in coastaldem_gebco_avg_lst_10m]            
            
            lst_topo_gebco_10m.append(gebco_lst_10m)
            lst_topo_gebco_avg_10m.append(gebco_avg_lst_10m)
            lst_topo_merit_10m.append(merit_lst_10m)
            lst_topo_merit_avg_10m.append(merit_avg_lst_10m) 
            lst_gebco_merit_10m.append(merit_gebco_lst_10m)
            lst_gebco_merit_avg_10m.append(merit_gebco_avg_lst_10m)   
            lst_topo_cstdem_10m.append(coastaldem_lst_10m)
            lst_topo_cstdem_avg_10m.append(coastaldem_avg_lst_10m)
            lst_topo_cstdem_gebco_10m.append(coastaldem_gebco_lst_10m)
            lst_topo_cstdem_gebco_avg_10m.append(coastaldem_gebco_avg_lst_10m)
            
            lst_topo_gebco_25m.append(gebco_lst_25m)
            lst_topo_gebco_avg_25m.append(gebco_avg_lst_25m)
            lst_topo_merit_25m.append(merit_lst_25m)
            lst_topo_merit_avg_25m.append(merit_avg_lst_25m)
            lst_gebco_merit_25m.append(merit_gebco_lst_25m)
            lst_gebco_merit_avg_25m.append(merit_gebco_avg_lst_25m)
            lst_topo_cstdem_25m.append(coastaldem_lst_25m)
            lst_topo_cstdem_avg_25m.append(coastaldem_avg_lst_25m)
            lst_topo_cstdem_gebco_25m.append(coastaldem_gebco_lst_25m)
            lst_topo_cstdem_gebco_avg_25m.append(coastaldem_gebco_avg_lst_25m)
            
            lst_topo_gebco_50m.append(gebco_lst_50m)
            lst_topo_gebco_avg_50m.append(gebco_avg_lst_50m)
            lst_topo_merit_50m.append(merit_lst_50m)
            lst_topo_merit_avg_50m.append(merit_avg_lst_50m)
            lst_gebco_merit_50m.append(merit_gebco_lst_50m)
            lst_gebco_merit_avg_50m.append(merit_gebco_avg_lst_50m) 
            lst_topo_cstdem_50m.append(coastaldem_lst_50m)
            lst_topo_cstdem_avg_50m.append(coastaldem_avg_lst_50m)
            lst_topo_cstdem_gebco_50m.append(coastaldem_gebco_lst_50m)
            lst_topo_cstdem_gebco_avg_50m.append(coastaldem_gebco_avg_lst_50m)
            
            lst_topo_gebco_100m.append(gebco_lst_100m)
            lst_topo_gebco_avg_100m.append(gebco_avg_lst_100m)
            lst_topo_merit_100m.append(merit_lst_100m)
            lst_topo_merit_avg_100m.append(merit_avg_lst_100m)
            lst_gebco_merit_100m.append(merit_gebco_lst_100m)
            lst_gebco_merit_avg_100m.append(merit_gebco_avg_lst_100m)
            lst_topo_cstdem_100m.append(coastaldem_lst_100m)
            lst_topo_cstdem_avg_100m.append(coastaldem_avg_lst_100m)
            lst_topo_cstdem_gebco_100m.append(coastaldem_gebco_lst_100m)
            lst_topo_cstdem_gebco_avg_100m.append(coastaldem_gebco_avg_lst_100m)

            lst_topo_gebco_250m.append(gebco_lst_250m)
            lst_topo_gebco_avg_250m.append(gebco_avg_lst_250m)
            lst_topo_merit_250m.append(merit_lst_250m)
            lst_topo_merit_avg_250m.append(merit_avg_lst_250m)
            lst_gebco_merit_250m.append(merit_gebco_lst_250m)
            lst_gebco_merit_avg_250m .append(merit_gebco_avg_lst_250m)
            lst_topo_cstdem_250m.append(coastaldem_lst_250m)
            lst_topo_cstdem_avg_250m.append(coastaldem_avg_lst_250m)
            lst_topo_cstdem_gebco_250m.append(coastaldem_gebco_lst_250m)
            lst_topo_cstdem_gebco_avg_250m.append(coastaldem_gebco_avg_lst_250m)
            
            lst_topo_gebco_500m.append(gebco_lst_500m)
            lst_topo_gebco_avg_500m.append(gebco_avg_lst_500m)
            lst_topo_merit_500m.append(merit_lst_500m)
            lst_topo_merit_avg_500m.append(merit_avg_lst_500m)
            lst_gebco_merit_500m.append(merit_gebco_lst_500m)
            lst_gebco_merit_avg_500m .append(merit_gebco_avg_lst_500m)
            lst_topo_cstdem_500m.append(coastaldem_lst_500m)
            lst_topo_cstdem_avg_500m.append(coastaldem_avg_lst_500m)
            lst_topo_cstdem_gebco_500m.append(coastaldem_gebco_lst_500m)
            lst_topo_cstdem_gebco_avg_500m.append(coastaldem_gebco_avg_lst_500m)
    
    cs_thk_mu, cs_thk_std = -1, -1                
    cs_width_mu, cs_width_std = -1, -1  
    pcr_rch_mu, pcr_rch_std = -1, -1  
    watergap_rch_mu, watergap_rch_std = -1, -1  
    p_min_et_rch_mu, p_min_et_rch_std = -1, -1  
    all_rch_mu, all_rch_std = -1, -1  
    soil_type_mu, soil_type_std = -1, -1  
    soil_thk_mu, soil_thk_std = -1, -1  
    glhymps_top_mu, glhymps_top_std = -1, -1  
    glhymps_bot_mu, glhymps_bot_std = -1, -1  
    k_soil_mu, k_soil_std = -1, -1         

    #   create plots for each coastal type 
    if lst_cs_thk == []:
        ls_elev = -1
        ls_elev_50 = -1
        pass
    else:
        try:
            
            #   save all the topo lists into netcdf files
            rivbas.create_topo_nc(lst_gebco_merit_avg_10m, lst_gebco_merit_avg_25m, lst_gebco_merit_avg_50m, lst_gebco_merit_avg_100m,\
                                  lst_gebco_merit_avg_250m, lst_gebco_merit_avg_500m, lst_topo_gebco_avg_10m, lst_topo_gebco_avg_25m,\
                                  lst_topo_gebco_avg_50m, lst_topo_gebco_avg_100m, lst_topo_gebco_avg_250m, lst_topo_gebco_avg_500m,\
                                  lst_topo_cstdem_gebco_avg_10m, lst_topo_cstdem_gebco_avg_25m, lst_topo_cstdem_gebco_avg_50m,\
                                  lst_topo_cstdem_gebco_avg_100m, lst_topo_cstdem_gebco_avg_250m, lst_topo_cstdem_gebco_avg_500m, id_riv_bas_dir)

            #   save all the input data lists into a numpy dictionary file so the script doesnt need to be rerun all the time if you want to change the IBOUND etc.
            input_dict = dict()
            input_dict['lst_cs_thk'] = lst_cs_thk
            input_dict['lst_cs_width'] = lst_cs_width
            input_dict['lst_anchor_dist'] = lst_anchor_dist
            input_dict['lst_anchor_depth'] = lst_anchor_depth
            input_dict['lst_wtd'] = lst_wtd
            input_dict['lst_gebco_merit_avg_100m'] = lst_gebco_merit_avg_100m
            input_dict['lst_topo_gebco_avg_500m'] = lst_topo_gebco_avg_500m
            input_dict['lst_soil_type'] = lst_soil_type
            input_dict['lst_soil_thk'] = lst_soil_thk
            input_dict['lst_cs_nasa'] = lst_cs_nasa
            input_dict['lst_cs_offshore'] = lst_cs_offshore
            input_dict['lst_cs_pcr_rch'] = lst_cs_pcr_rch
            input_dict['lst_cs_watergap_rch'] = lst_cs_watergap_rch
            input_dict['lst_cs_p_min_et'] = lst_cs_p_min_et
            input_dict['lst_cs_k_soil'] = lst_cs_k_soil
            input_dict['lst_cs_drn_rate'] = lst_cs_drn_rate
            input_dict['lst_cs_glhymps_top_lay'] = lst_cs_glhymps_top_lay
            input_dict['lst_cs_glhymps_bot_lay'] = lst_cs_glhymps_bot_lay
        
            dict_save_dir = os.path.join(id_riv_bas_dir, '_input_lists.npy')
            np.save(dict_save_dir, input_dict)
            """
            #   the topo data to use in models is - GEBCO_MERIT_100m as that combines updated inland DEM and 
            repr_model = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_gebco_merit_avg_100m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 100., 100., 10., id_cs, id_riv_bas_dir, '_merit_avg_100m')

            repr_model_gm_500 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_gebco_merit_avg_500m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 500., 100., 50., id_cs, id_riv_bas_dir, '_merit_avg_500m')

            repr_model_gm_250 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_gebco_merit_avg_250m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 250., 100., 25., id_cs, id_riv_bas_dir, '_merit_avg_250m')

            repr_model_gm_50 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_gebco_merit_avg_50m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 50., 50., 5., id_cs, id_riv_bas_dir, '_merit_avg_050m')
            
            repr_model_gm_25 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_gebco_merit_avg_25m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 25., 25., 2.5, id_cs, id_riv_bas_dir, '_merit_avg_025m')            

            repr_model_gm_10 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_gebco_merit_avg_10m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 10., 10., 1., id_cs, id_riv_bas_dir, '_merit_avg_010m')            
            
            repr_model_g_500 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_gebco_avg_500m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 500., 100., 50., id_cs, id_riv_bas_dir, '_gebco_avg_500m')

            repr_model_g_250 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_gebco_avg_250m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 250., 100., 25., id_cs, id_riv_bas_dir, '_gebco_avg_250m')

            repr_model_g_100 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_gebco_avg_100m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 100., 100., 10., id_cs, id_riv_bas_dir, '_gebco_avg_100m')

            repr_model_g_50 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_gebco_avg_50m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 50., 50., 5., id_cs, id_riv_bas_dir, '_gebco_avg_050m')

            repr_model_g_25 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_gebco_avg_25m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 25., 25., 2.5, id_cs, id_riv_bas_dir, '_gebco_avg_025m')

            repr_model_g_10 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_gebco_avg_10m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 10., 10., 1., id_cs, id_riv_bas_dir, '_gebco_avg_010m')

            repr_model_gcst_500 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_cstdem_gebco_avg_500m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 500., 100., 50., id_cs, id_riv_bas_dir, '_coastal_avg_500m')

            repr_model_gcst_250 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_cstdem_gebco_avg_250m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 250., 100., 25., id_cs, id_riv_bas_dir, '_coastal_avg_250m')

            repr_model_gcst_100 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_cstdem_gebco_avg_100m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 100., 100., 10., id_cs, id_riv_bas_dir, '_coastal_avg_100m')

            repr_model_gcst_50 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_cstdem_gebco_avg_50m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 50., 50., 5., id_cs, id_riv_bas_dir, '_coastal_avg_050m')

            repr_model_gcst_25 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_cstdem_gebco_avg_25m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 25., 25., 2.5, id_cs, id_riv_bas_dir, '_coastal_avg_025m')

            repr_model_gcst_10 = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo_cstdem_gebco_avg_10m, lst_topo_gebco_avg_500m, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate,\
                                                 lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, 10., 10., 1., id_cs, id_riv_bas_dir, '_coastal_avg_010m')

            fos = repr_model[-1]
            ls_elev = repr_model[-3]
            ls_elev_50 = repr_model[-2]  
            
            #   plot all the topographical profiles in one plot
            rivbas.plot_all_topo_profiles(lst_gebco_merit_avg_500m, lst_cs_thk, fos, 'MERIT_GEBCO_avg_500m', id_cs_str + '_REG_' + str_id_subreg, id_riv_bas_dir)
            rivbas.plot_all_topo_profiles(lst_topo_gebco_avg_500m, lst_cs_thk, fos, 'GEBCO_avg_500m', id_cs_str + '_REG_' + str_id_subreg, id_riv_bas_dir)
            rivbas.plot_all_topo_profiles(lst_topo_cstdem_gebco_avg_500m, lst_cs_thk, fos, 'CoastalDEM_GEBCO_avg_500m', id_cs_str + '_REG_' + str_id_subreg, id_riv_bas_dir)
            """
            #cs_thk = rivbas.get_normal_dist(lst_cs_thk, os.path.join(geometry_pic_dir, '_cs_thk.png'),\
            #x_axis_lbl = 'Log normal distribution of coastal thickness', nan_val = None, no_negative_vals = True, plot = True)
            #cs_thk_mu, cs_thk_std, cs_thk_data = cs_thk[0], cs_thk[1], cs_thk[2]
            
            #cs_width = rivbas.get_normal_dist(lst_cs_width, os.path.join(geometry_pic_dir, '_cs_width.png'),\
            #x_axis_lbl = 'Log normal distribution of coastal width', nan_val = None, no_negative_vals = True, plot = True)
            #cs_width_mu, cs_width_std, cs_width_data = cs_width[0], cs_width[1], cs_width[2]                

            #pcr_rch = rivbas.get_normal_dist(lst_cs_pcr_rch, os.path.join(rch_pic_dir, '_pcr_rch.png'),\
            #x_axis_lbl = 'Log normal distribution of PCR-GLOBWB recharge', nan_val = 0.0, no_negative_vals = True, plot = True)
            #pcr_rch_mu, pcr_rch_std, pcr_rch_data = pcr_rch[0], pcr_rch[1], pcr_rch[2]

            #watergap_rch = rivbas.get_normal_dist(lst_cs_watergap_rch, os.path.join(rch_pic_dir, '_watergap_rch.png'),\
            #x_axis_lbl = 'Log normal distribution of WATERGAP recharge', nan_val = 0.0, no_negative_vals = True, plot = True)
            #watergap_rch_mu, watergap_rch_std, watergap_rch_data = watergap_rch[0], watergap_rch[1], watergap_rch[2]

            #p_min_et_rch = rivbas.get_normal_dist(lst_cs_p_min_et, os.path.join(rch_pic_dir, '_p_min_et_rch.png'),\
            #x_axis_lbl = 'Log normal distribution of P - ET recharge', nan_val = 0.0, no_negative_vals = True, plot = True)
            #p_min_et_rch_mu, p_min_et_rch_std, p_min_et_rch_data = p_min_et_rch[0], p_min_et_rch[1], p_min_et_rch[2]
            
            #all_rch_lst = np.concatenate((pcr_rch, watergap_rch_data, p_min_et_rch_data), axis=0).tolist()
            #all_rch = rivbas.get_normal_dist(all_rch_lst, os.path.join(rch_pic_dir, '_all_rch.png'),\
            #x_axis_lbl = 'Log normal distribution of all inputs (summed) recharge', nan_val = 0.0, no_negative_vals = True, plot = True)
            #all_rch_mu, all_rch_std, all_data = all_rch[0], all_rch[1], all_rch[2]

            #soil_type = rivbas.get_normal_dist(lst_soil_type, os.path.join(geology_pic_dir, '_soil_type.png'),\
            #x_axis_lbl = 'Log normal distribution of soil type', nan_val = -1, no_negative_vals = True, plot = True)
            #soil_type_mu, soil_type_std, soil_type_data = soil_type[0], soil_type[1], soil_type[2]

            #soil_thk = rivbas.get_normal_dist(lst_soil_thk, os.path.join(geology_pic_dir, '_soil_thk.png'),\
            #x_axis_lbl = 'Log normal distribution of soil thickness', nan_val = -32768.0, no_negative_vals = True, plot = True)
            #soil_thk_mu, soil_thk_std, soil_thk_data = soil_thk[0], soil_thk[1], soil_thk[2]

            #glhymps_top = rivbas.get_normal_dist(lst_cs_glhymps_top_lay, os.path.join(geology_pic_dir, '_glhymps_top.png'),\
            #x_axis_lbl = 'Log normal distribution of GLHYMPS 2 Hk values', nan_val = -9999.0, no_negative_vals = True, plot = True)
            #glhymps_top_mu, glhymps_top_std, glhymps_top_data = glhymps_top[0], glhymps_top[1], glhymps_top[2]

            #glhymps_bot = rivbas.get_normal_dist(lst_cs_glhymps_bot_lay, os.path.join(geology_pic_dir, '_glhymps_bot.png'),\
            #x_axis_lbl = 'Log normal distribution of GLHYMPS 1 Hk values', nan_val = -9999.0, no_negative_vals = True, plot = True)
            #glhymps_bot_mu, glhymps_bot_std, glhymps_bot_data = glhymps_bot[0], glhymps_bot[1], glhymps_bot[2]

            #k_soil = rivbas.get_normal_dist(lst_cs_k_soil, os.path.join(geology_pic_dir, '_k_soil.png'),\
            #x_axis_lbl = 'Log normal distribution of Hk soil values', nan_val = -9999.0, no_negative_vals = True, plot = True)
            #k_soil_mu, k_soil_std, k_soil_data = k_soil[0], k_soil[1], k_soil[2]

            #   change status to 1, means finished succesfully
            #df_in = pd.read_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb.csv'))
            df_in = pd.read_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb_ADDED_SRMs.csv'))
            df_in.loc[df_in['id_loop'] == id_loop, 'finished'] = 1
            df_in.update(df_in)
            #df_in.to_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb.csv'), sep = ',', encoding = 'utf-8', index = False)
            df_in.to_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb_ADDED_SRMs.csv'), sep = ',', encoding = 'utf-8', index = False)
            del df_in        
        
        except TypeError:
            print('Something is fucked')
            #   change status to 1, means finished succesfully
            #df_in = pd.read_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb.csv'))
            df_in = pd.read_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb_ADDED_SRMs.csv'))
            df_in.loc[df_in['id_loop'] == id_loop, 'finished'] = -9
            df_in.update(df_in)
            #df_in.to_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb.csv'), sep = ',', encoding = 'utf-8', index = False)
            df_in.to_csv(os.path.join(input_files_dir, '_coscat_reg_main_tb_ADDED_SRMs.csv'), sep = ',', encoding = 'utf-8', index = False)
            del df_in
            pass

    """
    f = open(csv_dir_ts, 'a')
    f.write("\n")
    f.write(str(coscat_id) + ',' + str(reg_id) + ',' + str(len(id_cs_lst)) + ',' + str(fos) + ',' + str(ls_elev) + ',' + str(ls_elev_50)) # write the header
    f.close()    

    #f2 = open(csv_dir_ts_input_data_geom, 'a')
    #f2.write("\n")
    #f2.write(str(coscat_id) + ',' + str(reg_id) + ',' + str(cs_thk_mu) + ',' + str(cs_thk_std) + ',' + str(cs_width_mu) + ',' + str(cs_width_std))
    #f2.close()
    
    #f3 = open(csv_dir_ts_input_data_rch, 'a')
    #f3.write("\n")
    #f3.write(str(coscat_id) + ',' + str(reg_id) + ',' + str(pcr_rch_mu) + ',' + str(pcr_rch_std) + ',' + str(watergap_rch_mu) + ',' + str(watergap_rch_std)\
    #         + ',' + str(p_min_et_rch_mu) + ',' + str(p_min_et_rch_std) + ',' + str(all_rch_mu) + ',' + str(all_rch_std))              
    #f3.close()
            
    f4 = open(csv_dir_ts_input_data_geol, 'a')
    f4.write("\n")    
    f4.write(str(coscat_id) + ',' + str(reg_id) + ',' + str(soil_type_mu) + ',' + str(soil_type_std) + ',' + str(soil_thk_mu) + ',' + str(soil_thk_std)\
             + ',' + str(glhymps_top_mu) + ',' + str(glhymps_top_std) + ',' + str(glhymps_bot_mu) + ',' + str(glhymps_bot_std)\
             + ',' + str(k_soil_mu) + ',' + str(k_soil_std))    
    f4.close()
    """    
    
    
    
    
    
    
    
    
    
    
"""
Function to build an IBOUND array and other MODFLOW topographical lists necessary to build a model. Different topographical
lists (GEBCO or MERIT) can be filled in to test different DEM inputs. Another difference compared to models in A1, A2, A3 is the
possibility to change the resolution of the model. To implement also possibility to refine the grid at certain areas. 



pt_topo_dist = 100.
del_col = 100.
del_lay = 10.


thk_lst = lst_cs_thk
width_lst = lst_cs_width
anchor_dist_lst = lst_anchor_dist
anchor_depth_lst = lst_anchor_depth
wtd_lst = lst_wtd

topo_res = 100.
topo_lst = lst_gebco_merit_avg_100m
topo_gebco_500m_lst = lst_topo_gebco_avg_500m

soil_type_lst = lst_soil_type
soil_thk_lst = lst_soil_thk
nasa_lst = lst_cs_nasa
offshore_lst = lst_cs_offshore 
pcr_lst = lst_cs_pcr_rch
watergap_lst = lst_cs_watergap_rch
p_min_et_lst = lst_cs_p_min_et
k_soil_lst = lst_cs_k_soil
drn_rate_lst = lst_cs_drn_rate
riv_cond_01_lst = lst_cs_riv_cond_01
riv_cond_05_lst = lst_cs_riv_cond_05
riv_cond_10_lst = lst_cs_riv_cond_10
riv_width_lst = lst_cs_riv_width
riv_bot_elev_lst = lst_cs_riv_bot_elev
riv_head_elev_lst = lst_cs_riv_head_elev
glhymps_top_lay_lst = lst_cs_glhymps_top_lay
glhymps_bot_lay_lst = lst_cs_glhymps_bot_lay
id_coscat = id_cs
out_dir = id_riv_bas_dir      
topo_diff = 10.
bot_limit = -10000.
ibound_act_cells_limit = 2
offshore_dist_limit = 200000.
"""



    
    
    
    
    
    
    
