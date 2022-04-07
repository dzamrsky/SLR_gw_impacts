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
import fiona
import geopandas as gpd
import pandas as pd


import sys
sys.path.append(r'g:\Water_Nexus\_A4\_scripts\_fc_scripts')
import ws_masterscript_py3 as ws_ms
import ws_datatools_py3 as ws               # connect to dbase, run sqls etc.
from modeltools_py3 import cs_model 
import _create_ARP_RIVBAS_tools as rivbas

import gdal
import utm
import math
import rasterio

#   define the different raster files and databases that will be used to derive the various parameter values
#   flo1k = mean, maximum and minimum annual streamflow at 1 km resolution from 1960 through 2015
flo1k_dir = r'g:\_ORIGINAL_DATA\FLO1K\qav_mean_1km.tif'
cond_10_dir = r'e:\Water_Nexus\_A2\_data\_MODFLOW_input\_RIV_cond_M_val_1.0_noval_compressed.tif' 
cond_05_dir = r'e:\Water_Nexus\_A2\_data\_MODFLOW_input\_RIV_cond_M_val_0.5_noval_compressed.tif'
cond_01_dir = r'e:\Water_Nexus\_A2\_data\_MODFLOW_input\_RIV_cond_M_val_0.1_noval_compressed.tif'
pcr_cond_dir = r'e:\Water_Nexus\_A2\_data\_MODFLOW_input\use_this\riv_bed_conductance.map'
pcr_stage_dir = r'e:\Water_Nexus\_A2\_data\_MODFLOW_input\use_this\riv_head_surface_water_elevation.map'
pcr_rivbot_dir = r'e:\Water_Nexus\_A2\_data\_MODFLOW_input\use_this\riv_bottom_surface_water_bed_elevation.map'
gaia_width_dir = r'e:\Water_Nexus\_A2\_data\_GAIA_river_depth_width\_GAIA_river_width_1km.tif'
gaia_depth_dir = r'e:\Water_Nexus\_A2\_data\_GAIA_river_depth_width\_GAIA_river_depth_1km_compressed.tif'

ne_coastline = r'g:\_ORIGINAL_DATA\natural_earth\ne_10m_land'
coscat_polys = r'g:\_CREATED_DATA\_A2_data\coscat_analysis\coscat_dissolved'
# Read the original GLIM_su Shapefile
input = fiona.open('g:/_CREATED_DATA/_A2_data/GLIM_su_only.shp', 'r')

hybas_dir = r'g:\Water_Nexus\_A4\_HYBAS_input_files'            #   output folder for individual HYBAS models (Hydrosheds Basins)
plot_dir = r'g:\Water_Nexus\_A4\_HYBAS_input_files\_plots'

#   connect to the databases
#db_name = "'_coastal_dbase_sed_thick_valid'"
db_name = "'cs_db_deltas'"
db_name_world = "'cs_db_v1'"
db_host = "'localhost'"
db_user = "'postgres'"
db_pass = "'postgres'"

#table_db = "coastline_coscat"
table_db = "cs_coastal_types"

dbase_connect = ws.dbase_tools(db_name, db_user, db_host, db_pass)
dbase_connect_world = ws.dbase_tools(db_name_world, db_user, db_host, db_pass)

#   select all unique COSCAT id numbers 
riv_pts_dir = r'g:\Water_Nexus\_A4\_GISDATA\_AF_test_region\_TEST_AF_RIVBAS_PTS.shp' 
riv_pts_raw = gpd.read_file(riv_pts_dir)

#coscat_ids = dbase_connect_world.get_ids_coscat("coscat_id", "cs_coastal_types")
#coscat_ids = dbase_connect.get_ids_coscat("coscat_id", "cs_coastal_types")

csv_dir_ts = os.path.join(hybas_dir, '_COSCAT_representative_models_summary_ATE_MEAN_STDEV.csv')  

#   define the table in the database that contains the coastline points shapefile
#   this is not the cs table where the info is written to!
coastline_points_raw_table = 'ne_10m_coastline_points_5km'#'ne_points_5km'#'ne_points_5km_test_raster_edge'##'ne_coatline_points_5km_test_avg'#'ne_coatline_points_5km_test_avg'
#  determine how many points to create (25 means point per 1km)
#  and the distance of the farthest point from point_5km
n_points = 799
max_dist = 399500
#   define the number of points and distance from the cross-section point
n_points_avg = 5
avg_dist = 2500

csv_out_hdrs = ['id_cs', 'direction', 'dist_to_cst', 'X_wgs84', 'Y_wgs84', 'flo1k_val', 'flo1k_val_AVG', 'COND_m_0.1', 'COND_m_0.5', 'COND_m_1.0',\
                'R_cond_PCR', 'R_head_elev_PCR', 'R_riv_bot_elev_PCR', 'R_width_GAIA', 'R_depth_GAIA']





csv_dir_ts = os.path.join(hybas_dir, '_HYBAS_representative_models_summary.csv')           
f = open(csv_dir_ts,'w')
f.write('id_HYBAS, n, fos, ls_elev, ls_elev_50') #write the header
#   n = number of profiles in the HYBAS regions
#   fos = FOS point found in the averaged Henry profile, if no = 0, if yes = distance from coast in km
#               this is not based on the find_FOS function as the coastline can be shifted depending on the top_elev list,
#               instead it is defined based on the last column in the IBOUND_array and the number of layers there 
#   ls_elev = either 0 or 1 (false or true), if all the elevation points are > -120m offshore then = 0
#   lst_elev_50 = either 0 or 1 (false or true), if all more than half of the the elevation points are > -120m offshore then = 0
f.close() 

csv_dir_ts_input_data_geom = os.path.join(hybas_dir, '_HYBAS_representative_models_GEOMETRY_stats.csv')           
f2 = open(csv_dir_ts_input_data_geom,'w')
f2.write('id_HYBAS, cs_thk_mu, cs_thk_std, cs_width_mu, cs_width_std')
f2.close()

csv_dir_ts_input_data_rch = os.path.join(hybas_dir, '_HYBAS_representative_models_RCH_stats.csv')           
f3 = open(csv_dir_ts_input_data_rch,'w')
f3.write('id_HYBAS, cs_pcr_rch_mu, cs_pcr_rch_std, cs_watergap_rch_mu, cs_watergap_rch_std, cs_p_min_et_rch_mu, cs_p_min_et_rch_std,\
all_rch_mu, all_rch_std') #write the header
f3.close()

csv_dir_ts_input_data_geol = os.path.join(hybas_dir, '_HYBAS_representative_models_GEOLOGY_stats.csv')           
f4 = open(csv_dir_ts_input_data_geol,'w')
f4.write('id_HYBAS, soil_type_mu, soil_type_std, soil_thk_mu, soil_thk_std,\
glhymps_top_mu, glhymps_top_std, glhymps_bot_mu, glhymps_bot_std, soil_k_mu')         
f4.close()











#   make list of all BASIN id numbers from the input shapefile
id_bas_lst = list(set(list(riv_pts_raw['id_basin'])))
#   loop through the BASIN id numbers and select all the individual cross-section points
for id_bas in id_bas_lst:
    #id_bas = 79190
    
    csv_out_lst = []
    id_cs_lst = list(riv_pts_raw.loc[riv_pts_raw['id_basin'] == id_bas]['id_cs'])
    
    #   define string and integer 
    id_bas_str = "ID_riv_bas = " + str(id_bas)
    id_bas_int = int(id_bas)

    #   add 0s in front of the name to have an ordered list of COSCAT ids in the folder 
    if id_bas_int < 10:
        str_id_riv_bas = '00000' + str(id_bas_int)
    elif id_bas_int >= 10 and id_bas_int < 100:
        str_id_riv_bas = '0000' + str(id_bas_int)
    elif id_bas_int >= 100 and id_bas_int < 1000:
        str_id_riv_bas = '000' + str(id_bas_int)    
    elif id_bas_int >= 1000 and id_bas_int < 10000:
        str_id_riv_bas = '00' + str(id_bas_int)    
    elif id_bas_int >= 10000 and id_bas_int < 100000:
        str_id_riv_bas = '0' + str(id_bas_int)    
    else:
        str_id_riv_bas = str(id_bas_int)  

    #   define the model directory where the model files and figures will be stored
    id_riv_bas_dir = os.path.join(hybas_dir, str_id_riv_bas)
    rch_pic_dir = os.path.join(id_riv_bas_dir, '_RCH_hist')
    geometry_pic_dir = os.path.join(id_riv_bas_dir, '_GEOMETRY_hist')
    geology_pic_dir = os.path.join(id_riv_bas_dir, '_GEOLOGY_hist')    

    if not os.path.exists(id_riv_bas_dir):
        os.makedirs(id_riv_bas_dir)
        os.makedirs(rch_pic_dir)
        os.makedirs(geometry_pic_dir)
        os.makedirs(geology_pic_dir)
    else:
        print('Directories for RIV BASIN id: ' + str_id_riv_bas + ' already exists.')

    #   define the model dimensions and empty lists to be filled 
    del_col = 100.
    del_lay = 10.
    cs_points_dist = 500.
    
    lst_cs_width, lst_cs_thk = [], []
    lst_anchor_dist, lst_anchor_depth, lst_wtd, lst_topo, lst_soil_type, lst_soil_thk   = [], [], [], [], [], [] 
    lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch, lst_cs_p_min_et, lst_cs_k_soil = [], [], [], [], [], [] 
    lst_cs_drn_rate, lst_cs_riv_cond_01, lst_cs_riv_cond_05, lst_cs_riv_cond_10, lst_cs_riv_width, lst_cs_riv_bot_elev = [], [], [], [], [], [] 
    lst_cs_riv_head_elev, lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay = [], [], []
    
    #   set the counters of different coastal types
    pts_count = 0
    fos = -999
    ls_elev = -999
    ls_elev_50 = -999   

    #   loop through the cross-sections and select required info from all the necessary tables
    for id_num in id_cs_lst:
        id_cs = int(id_num)
        #   info about coastal plain width and the anchor point
        cs_info = dbase_connect.get_data_condition("cst_plain_width, nasa_point_dist, nasa_point_depth, overall_avg",\
                                                   "id_cs = " + str(id_cs), "cs_sed_thick_est_v2")
        #   also get the classification/delta of the coastal point
        cond_str_id_cs = "id_cs = '" + str(id_cs) + "'"
        
        cs_class = dbase_connect_world.get_data_condition("class_fin", cond_str_id_cs, "cs_coastal_types")
        coastal_class = cs_class[0][0]

        try:
            cs_plain_width, anchor_dist, anchor_depth, avg_depth = float(cs_info[0][0]), 200.0 - float(cs_info[0][1]), float(cs_info[0][2]), float(cs_info[0][3]) 
        except IndexError:
            continue
        
        #   if it is delta then and the anchor depth is more than 200km then assign 200km
        if anchor_dist < 0.:
            anchor_dist = 200.0 + anchor_dist

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
            cs_pcr_rch = cs_model_point.dta_pcr_rch
            cs_watergap_rch = cs_model_point.dta_watergap_rch
            cs_p_min_et = cs_model_point.dta_p_min_et
            cs_k_soil = cs_model_point.dta_k_soil
            cs_drn_rate = cs_model_point.dta_drn_rate
            cs_riv_cond_01 = cs_model_point.dta_riv_cond_01       
            cs_riv_cond_05 = cs_model_point.dta_riv_cond_05                      
            cs_riv_cond_10 = cs_model_point.dta_riv_cond_10                   
            cs_riv_width = cs_model_point.dta_riv_width       
            cs_riv_bot_elev = cs_model_point.dta_riv_bot_elev                     
            cs_riv_head_elev = cs_model_point.dta_riv_head_elev   
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
            lst_cs_pcr_rch.append(cs_pcr_rch)
            lst_cs_watergap_rch.append(cs_watergap_rch)
            lst_cs_p_min_et.append(cs_p_min_et)
            lst_cs_k_soil.append(cs_k_soil)
            lst_cs_drn_rate.append(cs_drn_rate)  
            lst_cs_riv_cond_01.append(cs_riv_cond_01)
            lst_cs_riv_cond_05.append(cs_riv_cond_05)
            lst_cs_riv_cond_10.append(cs_riv_cond_10)
            lst_cs_riv_width.append(cs_riv_width)
            lst_cs_riv_bot_elev.append(cs_riv_bot_elev)
            lst_cs_riv_head_elev.append(cs_riv_head_elev)
            lst_cs_glhymps_top_lay.append(cs_glhymps_top_lay)
            lst_cs_glhymps_bot_lay.append(cs_glhymps_bot_lay)                

            cont_direction = cs_model_point.land_pos
    
            pts_count += 1             
            
        except IndexError:

            pts_count += 1    
            #continue
            
            cs_topo = dbase_connect.get_data_condition("*", "id_cs = " + str(id_cs), "cs_gebco_2014_avg")
            cs_topo = cs_topo[0][399:1200]
            cs_nasa = dbase_connect.get_data_condition("*", "id_cs = " + str(id_cs), "nasa_thick_avg")
            cs_nasa = cs_nasa[0][399:1200]
            cs_pcr = dbase_connect.get_data_condition("*", "id_cs = " + str(id_cs), "cs_aq_thick")
            cs_pcr = cs_pcr[0][399:1200]
            cs_glim = dbase_connect.get_data_condition("*", "id_cs = " + str(id_cs), "cs_glim_litho")
            cs_glim = cs_glim[0][399:1200]
            cs_nasa_flt = cs_nasa#[float(i) for i in cs_nasa]    

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

            else:
                pass     

            lst_cs_thk.append(avg_depth)
            lst_cs_width.append(cs_plain_width)
            lst_anchor_dist.append(anchor_dist)
            lst_anchor_depth.append(anchor_depth)     
            lst_cs_nasa.append(cs_nasa_flt)
            lst_topo.append(cs_topo)

        x_cnt, y_cnt = 43200, 21600
        x_st, y_st = -180, 90 #  lower left corner
        px_x, px_y = abs((2 * x_st) / x_cnt), abs((2 * y_st) / y_cnt) * -1
        cs_id_out_lst = rivbas.get_riv_input(id_cs, flo1k_dir, cond_10_dir, cond_05_dir, cond_01_dir, pcr_cond_dir, pcr_stage_dir, pcr_rivbot_dir,\
                                             gaia_width_dir, gaia_depth_dir, x_st, y_st, px_x, px_y, cs_plain_width, cont_direction)
        csv_out_lst.append(cs_id_out_lst[:])
        
    fin_riv_bas_csv = pd.DataFrame.from_records([item for sublist in csv_out_lst for item in sublist], columns = csv_out_hdrs)
    fin_riv_bas_csv.to_csv(os.path.join(id_riv_bas_dir, '_flo1k_vals.csv'), sep = ',')

    riv_dict_input = rivbas.get_riv_info(os.path.join(id_riv_bas_dir, '_flo1k_vals.csv'), id_cs_lst, max(lst_cs_width) * 1000.)

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
            repr_model = rivbas.save_model_input_files_cst_type_topo('_', lst_cs_thk, lst_cs_width, lst_anchor_dist, lst_anchor_depth,\
                                                 lst_wtd, lst_topo, lst_soil_type, lst_soil_thk,\
                                                 lst_cs_nasa, lst_cs_offshore, lst_cs_pcr_rch, lst_cs_watergap_rch,\
                                                 lst_cs_p_min_et, lst_cs_k_soil, lst_cs_drn_rate, lst_cs_riv_cond_01,\
                                                 lst_cs_riv_cond_05, lst_cs_riv_cond_10, lst_cs_riv_width, lst_cs_riv_bot_elev,\
                                                 lst_cs_riv_head_elev, lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay, id_cs, id_riv_bas_dir)
            
            fos = repr_model[-1]
            ls_elev = repr_model[-3]
            ls_elev_50 = repr_model[-2]  
            fos_pt = rivbas.find_FOS(lst_topo, id_riv_bas_dir, 'avg_avg', avg = True)
            
            rivbas.plot_all_topo_profiles(id_cs, lst_topo, lst_cs_thk, fos_pt, '_', id_riv_bas_dir)
            
            cs_thk = rivbas.get_normal_dist(lst_cs_thk, os.path.join(geometry_pic_dir, '_cs_thk.png'),\
            x_axis_lbl = 'Log normal distribution of coastal thickness', nan_val = None, no_negative_vals = True, plot = True)
            cs_thk_mu, cs_thk_std, cs_thk_data = cs_thk[0], cs_thk[1], cs_thk[2]
            
            cs_width = rivbas.get_normal_dist(lst_cs_width, os.path.join(geometry_pic_dir, '_cs_width.png'),\
            x_axis_lbl = 'Log normal distribution of coastal width', nan_val = None, no_negative_vals = True, plot = True)
            cs_width_mu, cs_width_std, cs_width_data = cs_width[0], cs_width[1], cs_width[2]                

            pcr_rch = rivbas.get_normal_dist(lst_cs_pcr_rch, os.path.join(rch_pic_dir, '_pcr_rch.png'),\
            x_axis_lbl = 'Log normal distribution of PCR-GLOBWB recharge', nan_val = 0.0, no_negative_vals = True, plot = True)
            pcr_rch_mu, pcr_rch_std, pcr_rch_data = pcr_rch[0], pcr_rch[1], pcr_rch[2]

            watergap_rch = rivbas.get_normal_dist(lst_cs_watergap_rch, os.path.join(rch_pic_dir, '_watergap_rch.png'),\
            x_axis_lbl = 'Log normal distribution of WATERGAP recharge', nan_val = 0.0, no_negative_vals = True, plot = True)
            watergap_rch_mu, watergap_rch_std, watergap_rch_data = watergap_rch[0], watergap_rch[1], watergap_rch[2]

            p_min_et_rch = rivbas.get_normal_dist(lst_cs_p_min_et, os.path.join(rch_pic_dir, '_p_min_et_rch.png'),\
            x_axis_lbl = 'Log normal distribution of P - ET recharge', nan_val = 0.0, no_negative_vals = True, plot = True)
            p_min_et_rch_mu, p_min_et_rch_std, p_min_et_rch_data = p_min_et_rch[0], p_min_et_rch[1], p_min_et_rch[2]
            
            #all_rch_lst = np.concatenate((pcr_rch, watergap_rch_data, p_min_et_rch_data), axis=0).tolist()
            #all_rch = rivbas.get_normal_dist(all_rch_lst, os.path.join(rch_pic_dir, '_all_rch.png'),\
            #x_axis_lbl = 'Log normal distribution of all inputs (summed) recharge', nan_val = 0.0, no_negative_vals = True, plot = True)
            #all_rch_mu, all_rch_std, all_data = all_rch[0], all_rch[1], all_rch[2]

            soil_type = rivbas.get_normal_dist(lst_soil_type, os.path.join(geology_pic_dir, '_soil_type.png'),\
            x_axis_lbl = 'Log normal distribution of soil type', nan_val = -1, no_negative_vals = True, plot = True)
            soil_type_mu, soil_type_std, soil_type_data = soil_type[0], soil_type[1], soil_type[2]

            soil_thk = rivbas.get_normal_dist(lst_soil_thk, os.path.join(geology_pic_dir, '_soil_thk.png'),\
            x_axis_lbl = 'Log normal distribution of soil thickness', nan_val = -32768.0, no_negative_vals = True, plot = True)
            soil_thk_mu, soil_thk_std, soil_thk_data = soil_thk[0], soil_thk[1], soil_thk[2]

            glhymps_top = rivbas.get_normal_dist(lst_cs_glhymps_top_lay, os.path.join(geology_pic_dir, '_glhymps_top.png'),\
            x_axis_lbl = 'Log normal distribution of GLHYMPS 2 Hk values', nan_val = -9999.0, no_negative_vals = True, plot = True)
            glhymps_top_mu, glhymps_top_std, glhymps_top_data = glhymps_top[0], glhymps_top[1], glhymps_top[2]

            glhymps_bot = rivbas.get_normal_dist(lst_cs_glhymps_bot_lay, os.path.join(geology_pic_dir, '_glhymps_bot.png'),\
            x_axis_lbl = 'Log normal distribution of GLHYMPS 1 Hk values', nan_val = -9999.0, no_negative_vals = True, plot = True)
            glhymps_bot_mu, glhymps_bot_std, glhymps_bot_data = glhymps_bot[0], glhymps_bot[1], glhymps_bot[2]

            k_soil = rivbas.get_normal_dist(lst_cs_k_soil, os.path.join(geology_pic_dir, '_k_soil.png'),\
            x_axis_lbl = 'Log normal distribution of Hk soil values', nan_val = -9999.0, no_negative_vals = True, plot = True)
            k_soil_mu, k_soil_std, k_soil_data = k_soil[0], k_soil[1], k_soil[2]
            
        except TypeError:
            pass

    f = open(csv_dir_ts, 'a')
    f.write("\n")
    f.write(str(id_bas) + ',' + str(len(id_cs_lst)) + ',' + str(fos) + ',' + str(ls_elev) + ',' + str(ls_elev_50)) # write the header
    f.close()    

    f2 = open(csv_dir_ts_input_data_geom, 'a')
    f2.write("\n")
    f2.write(str(id_bas) + ',' + str(cs_thk_mu) + ',' + str(cs_thk_std) + ',' + str(cs_width_mu) + ',' + str(cs_width_std))
    f2.close()
    
    f3 = open(csv_dir_ts_input_data_rch, 'a')
    f3.write("\n")
    f3.write(str(id_bas) + ',' + str(pcr_rch_mu) + ',' + str(pcr_rch_std) + ',' + str(watergap_rch_mu) + ',' + str(watergap_rch_std)\
             + ',' + str(p_min_et_rch_mu) + ',' + str(p_min_et_rch_std) + ',' + str(all_rch_mu) + ',' + str(all_rch_std))              
    f3.close()
            
    f4 = open(csv_dir_ts_input_data_geol, 'a')
    f4.write("\n")    
    f4.write(str(id_bas) + ',' + str(soil_type_mu) + ',' + str(soil_type_std) + ',' + str(soil_thk_mu) + ',' + str(soil_thk_std)\
             + ',' + str(glhymps_top_mu) + ',' + str(glhymps_top_std) + ',' + str(glhymps_bot_mu) + ',' + str(glhymps_bot_std)\
             + ',' + str(k_soil_mu) + ',' + str(k_soil_std))    
    f4.close()





































"""

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

"""














