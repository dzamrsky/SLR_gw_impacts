# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 10:24:58 2020

@author: daniel
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

#main_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models'
out_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models_OUT_files'
input_dir = r'/projects/0/qt16165/_dzamrsky/_A4_MODEL_input_files'
master_csv_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models/_coscat_reg_main_tb_MODELRUNS.csv' 
#id_loop = int(sys.argv[1])


#main_dir = r'g:\Water_Nexus\_A4_models\_SLR_models'
out_dir = r'g:\Water_Nexus\_A4_models\_SLR_models_OUT_files'
input_dir = r'g:\Water_Nexus\_A4_models\_SLR_models_OUT_files\_input_files'
master_csv_dir = r'g:\Water_Nexus\_A4_models\_SLR_models\_coscat_reg_main_tb_MODELRUNS.csv' 
#id_loop = 3

os.makedirs(out_dir, exist_ok = True)
os.makedirs(os.path.join(out_dir, '_inl_fresh_csv'), exist_ok = True)
os.makedirs(os.path.join(out_dir, '_shlf_fresh_csv'), exist_ok = True)
os.makedirs(os.path.join(out_dir, '_avg_nc_files'), exist_ok = True)
os.makedirs(os.path.join(out_dir, '_input_files'), exist_ok = True)
os.makedirs(os.path.join(out_dir, '_conc_plots_files'), exist_ok = True)

#   the final output csv with global results 
#       IFGW - Inland Fresh GroundWater % - in the first 10km of the coastline!!
#       SWW - Salt Water Wedge (km inland reach)
#       DEMs -  me (MERIT), co (CoastalDEM), ge (GEBCO)
out_cols = ['fid', 'coscat', 'srm', 'cs_srm_id', 'inl_ext_km',\
            'IFGW_me_2000', 'IFGW_co_2000', 'IFGW_ge_2000',\
            'SWW_me_2000', 'SWW_co_2000', 'SWW_ge_2000',\
            'IFGW_me_2050', 'IFGW_co_2050', 'IFGW_ge_2050',\
            'SWW_me_2050', 'SWW_co_2050', 'SWW_ge_2050',\
            'IFGW_me_2100', 'IFGW_co_2100', 'IFGW_ge_2100',\
            'SWW_me_2100', 'SWW_co_2100', 'SWW_ge_2100',\
            'IFGW_me_2200', 'IFGW_co_2200', 'IFGW_ge_2200',\
            'SWW_me_2200', 'SWW_co_2200', 'SWW_ge_2200',\
            'IFGW_me_2300', 'IFGW_co_2300', 'IFGW_ge_2300',\
            'SWW_me_2300', 'SWW_co_2300', 'SWW_ge_2300',\
            'IFGW_me_2400', 'IFGW_co_2400', 'IFGW_ge_2400',\
            'SWW_me_2400', 'SWW_co_2400', 'SWW_ge_2400',\
            'IFGW_me_2500', 'IFGW_co_2500', 'IFGW_ge_2500',\
            'SWW_me_2500', 'SWW_co_2500', 'SWW_ge_2500']

#   create the csv files with only the column names, if they already exist dont create them
if not os.path.exists(os.path.join(out_dir, 'RCP_26_results.csv')):
    df = pd.DataFrame(columns = out_cols)
    df.set_index('fid')
    df.to_csv(os.path.join(out_dir, 'RCP_26_results.csv') , encoding = 'utf-8', index=False)
if not os.path.exists(os.path.join(out_dir, 'RCP_45_results.csv')):
    df = pd.DataFrame(columns = out_cols)
    df.set_index('fid')
    df.to_csv(os.path.join(out_dir, 'RCP_45_results.csv') , encoding = 'utf-8', index=False)
if not os.path.exists(os.path.join(out_dir, 'RCP_85_results.csv')):
    df = pd.DataFrame(columns = out_cols)
    df.set_index('fid')
    df.to_csv(os.path.join(out_dir, 'RCP_85_results.csv') , encoding = 'utf-8', index=False)

for id_loop in range(1, 225):

    try:
        #   get the region and coscat ID numbers
        df_in = pd.read_csv(master_csv_dir)
        id_loop_row = df_in.loc[df_in['id_loop'] == id_loop]
        coscat_id = int(id_loop_row['coscat_id'].values[0])
        subreg_id = int(id_loop_row['cst_reg_id'].values[0])
        
        #   add 0s in front of the name to have an ordered list of COSCAT ids in the folder 
        if subreg_id < 10:
            str_id_subreg = '00' + str(subreg_id)
        elif subreg_id >= 10 and subreg_id < 100:
            str_id_subreg = '0' + str(subreg_id)
        else:
            str_id_subreg = str(subreg_id)  
        
        #   define the right string from the coscat id number
        if coscat_id < 10:
            id_cs_str = '000' + str(coscat_id)
        elif coscat_id >= 10 and coscat_id < 100:
            id_cs_str = '00' + str(coscat_id)
        elif coscat_id >= 100 and coscat_id < 1000:
            id_cs_str = '0' + str(coscat_id)
        else:
            id_cs_str = str(coscat_id)  
        
        """
        rcp_26_lst, rcp_45_lst, rcp_85_lst = [], [], []
        for (dirpath, dirnames, filenames) in os.walk(main_dir):
            for file in filenames:
                if 'RCP_26' in file and '_MERGED.nc' in file and 'sc' not in file:
                    rcp_26_lst.append(os.path.join(dirpath, file))
                elif 'RCP_45' in file and '_MERGED.nc' in file and 'sc' not in file:
                    rcp_45_lst.append(os.path.join(dirpath, file))
                elif 'RCP_85' in file and '_MERGED.nc' in file and 'sc' not in file:
                    rcp_85_lst.append(os.path.join(dirpath, file))
        """
        
        #   function to calculate the % fresh in the 10km inland domain
        #conc_arr = conc_merit_2050
        #x_coords = x_lst.tolist()
        def fresh_inland(conc_arr, x_coords):
            #   find indexes of coastline (x = 0.05) and the inland extendt up to 10 km 
            x_coords = x_coords.tolist()
            x_0 = x_coords.index([i for i in x_coords if i > -10.][0])
            x_1 = x_coords.index([i for i in x_coords if i < 0.][-1])
            #   select the part of the conc array
            frsh_arr = conc_arr[:, x_0 : x_1 + 1]
            frsh_cnt, cells_tot = 0, 0
            for i in range(frsh_arr.shape[0]):
                for j in range(frsh_arr.shape[1]):
                    if frsh_arr[i, j] < 0.1:
                        frsh_cnt += 1
                    if not math.isnan(frsh_arr[i, j]):
                        cells_tot += 1
            #   return the final percentage of freshwater inland
            pct_out = round((frsh_cnt / cells_tot) * 100., 2)
            return pct_out

        #   function that find the maximum inland extent of saline/brackish
        #conc_arr = conc_coastal_2000
        #x_coords = x_lst.tolist()
        def sww_inland(conc_arr, idx_cst):
            #   find indexes of coastline (x = 0.05) and the inland extendt up to 10 km 
            frsh_arr = conc_arr[:, : idx_cst]
            sww = 0.0
            for i in range(frsh_arr.shape[1] - 1, 0, -1):
                #   get all layers in that column
                col_lst = frsh_arr[:, i].tolist()
                non_frsh = [j for j in col_lst if j > 17.5]
                if len(non_frsh) != 0:
                    #print(i)
                    sww += 0.1
                else:
                    pass
            return round(sww, 1)            
            
            """
            x_coords = x_coords.tolist()
            x_0 = x_coords.index([i for i in x_coords if i > -10.][0])
            x_1 = x_coords.index([i for i in x_coords if i < 0.][-1])
            #   select the part of the conc array
            frsh_arr = conc_arr[:, x_0 : x_1 + 1]
            for i in range(frsh_arr.shape[1]):
                #   get all layers in that column
                col_lst = frsh_arr[:, i].tolist()
                non_frsh = [j for j in col_lst if j > 34.5]
                if len(non_frsh) != 0:
                    sww = (x_1 - i) / 10.
                    break
                else:
                    sww = (x_1 - i) / 10.
            return sww
            """
            
        #   for each rcp in the list below create the final plots and add to the final CSV file 
        rcp_names = ['RCP_26', 'RCP_45', 'RCP_85']
        
        #   go through the list of nc files and get all values for each column + create finel output pictures
        for i in range(len(rcp_names)):
            rcp_name = rcp_names[i]
            #   get the coscat and srm number 
            cs_srm_str = id_cs_str + '_SRM_' + str_id_subreg
            coscat_str = int(cs_srm_str.split('_')[0].replace('0', ''))
            srm_str = int(cs_srm_str.split('SRM_')[1].replace('0', ''))
            cs_srm_id = cs_srm_str.split('_SRM_')[0] + '_' + str(srm_str)
            
            #   load the netcdf files, if the average files are not created then try to load the individual files 
            #sc_folder = os.path.join(main_dir, '_summary_' + cs_srm_str)
            avg_nc_folder = os.path.join(out_dir, '_avg_nc_files')
            
            #   load in the different DEM topographies
            try:
                input_dict = np.load(os.path.join(input_dir, cs_srm_str.split('_')[0], cs_srm_str.split('_')[0] + '_REG_' + cs_srm_str.split('SRM_')[1], '_gebco_100_IBOUND_IBOUND.npy'), allow_pickle = True).item()
                cst_idx = input_dict[('_%s_' + '%s_IBOUND') % ('gebco', int(100))]['cst_idx']
                cst_offset = input_dict[('_%s_' + '%s_IBOUND') % ('gebco', int(100))]['cst_offset']
                topo_dict = xr.open_dataset(os.path.join(input_dir, cs_srm_str.split('_')[0], cs_srm_str.split('_')[0] + '_REG_' + cs_srm_str.split('SRM_')[1], '_topo_100m_avg.nc'))
            except OSError:
                continue
        
            try:
                nc_in_coastal = xr.open_dataset(os.path.join(avg_nc_folder, cs_srm_str + '_' + rcp_name + '_coastal_MERGED.nc'))
                conc_coastal_2000 = nc_in_coastal.sel(time = 30000)['solute concentration'].values[:, :]
                conc_coastal_2050 = nc_in_coastal.sel(time = 30050)['solute concentration'].values[:, :]
                conc_coastal_2100 = nc_in_coastal.sel(time = 30100)['solute concentration'].values[:, :]
                conc_coastal_2200 = nc_in_coastal.sel(time = 30200)['solute concentration'].values[:, :]
                conc_coastal_2300 = nc_in_coastal.sel(time = 30300)['solute concentration'].values[:, :]
                conc_coastal_2400 = nc_in_coastal.sel(time = 30400)['solute concentration'].values[:, :]
                conc_coastal_2500 = nc_in_coastal.sel(time = 30500)['solute concentration'].values[:, :]
        
                #   read in the extent inland 
                inl_ext_km = nc_in_coastal['x'].values[0]
                #   get the x and y coordinates
                x_lst =  nc_in_coastal.coords['x'].values
                y_lst = np.arange(nc_in_coastal.coords['y'].values[0], nc_in_coastal.coords['y'].values[0] - conc_coastal_2000.shape[0] * 10., -10.)
                #y_lst =  nc_in.coords['y'].values
                if x_lst[0] > 0.:
                    x_st = max(cst_idx * -0.1, - 10.05)
                    x_lst = np.arange(cst_idx * -0.1 + 0.05, cst_idx * -0.1 + 0.05 + len(x_lst) * 0.1, 0.1)
                else:
                    x_st = max(x_lst[0], - 10.05)
                x_end = 5.05
                
                y_top_lst = conc_coastal_2000[:, 0].tolist()
                st = 1
                while len([x for x in y_top_lst if not math.isnan(x)]) == 0:
                    y_top_lst = conc_coastal_2000[:, st].tolist()
                    st += 1
                top = y_lst[y_top_lst.index(next(x for x in y_top_lst if not math.isnan(x)))] + 5.0
                    
                y_bot_lst = conc_coastal_2000[:, cst_idx].tolist()
                
                #   check how many nan values are in the beginning of the list
                bot_idx_1 = y_bot_lst.index([x for x in y_bot_lst if not math.isnan(x)][0])
                bot_idx_2 = len([x for x in y_bot_lst if not math.isnan(x)])
                try:
                    bot = y_lst[bot_idx_1 + bot_idx_2] + 5.0
                except IndexError:
                    bot = y_lst[-7] + 5.0
                    
                #top = y_lst[y_top_lst.index(next(x for x in y_top_lst if not math.isnan(x)))] + 5.0
                #y_bot_lst = nc_in_coastal.sel(x = x_end)['solute concentration'].values[-1].tolist()
                #bot = y_lst[y_bot_lst.index([x for x in y_bot_lst if not math.isnan(x)][-1])] + 5.0
                
                del nc_in_coastal
        
            except (FileNotFoundError, OSError):
                """
                try:
                    conc_coastal_2000 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30000.nc'))['solute concentration'].values[:, :]
                    conc_coastal_2050 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30050.nc'))['solute concentration'].values[:, :]
                    conc_coastal_2100 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30100.nc'))['solute concentration'].values[:, :]
                    conc_coastal_2200 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30200.nc'))['solute concentration'].values[:, :]
                    conc_coastal_2300 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30300.nc'))['solute concentration'].values[:, :]
                    conc_coastal_2400 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30400.nc'))['solute concentration'].values[:, :]
                    conc_coastal_2500 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30500.nc'))['solute concentration'].values[:, :]           
        
                    #   read in the extent inland 
                    inl_ext_km = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30000.nc'))['x'].values[0]
                    #   get the x and y coordinates
                    x_lst =  xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30000.nc')).coords['x'].values
                    y_lst = np.arange(xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30000.nc')).coords['y'].values[0],\
                                      xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30000.nc')).coords['y'].values[0] - conc_coastal_2000.shape[0] * 10., -10.)
                    #y_lst =  nc_in.coords['y'].values
                    x_st = max(x_lst[0], - 10.05)
                    x_end = 5.05
                    y_top_lst = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30000.nc')).sel(x = x_st)['solute concentration'].values[-1].tolist()
                    top = y_lst[y_top_lst.index(next(x for x in y_top_lst if not math.isnan(x)))]
                    
                    y_bot_lst = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'coastal', '_average_30000.nc')).sel(x = x_end)['solute concentration'].values[-1].tolist()
                    bot = y_lst[y_bot_lst.index(next(x for x in y_bot_lst if not math.isnan(x)))]
        
                except (FileNotFoundError, OSError):     
                """
                conc_coastal_2000 = None
                conc_coastal_2050 = None
                conc_coastal_2100 = None
                conc_coastal_2200 = None
                conc_coastal_2300 = None
                conc_coastal_2400 = None
                conc_coastal_2500 = None
                
            try:
                nc_in_merit = xr.open_dataset(os.path.join(avg_nc_folder, cs_srm_str + '_' + rcp_name + '_merit_MERGED.nc'))
                conc_merit_2000 = nc_in_merit.sel(time = 30000)['solute concentration'].values[:, :]
                conc_merit_2050 = nc_in_merit.sel(time = 30050)['solute concentration'].values[:, :]
                conc_merit_2100 = nc_in_merit.sel(time = 30100)['solute concentration'].values[:, :]
                conc_merit_2200 = nc_in_merit.sel(time = 30200)['solute concentration'].values[:, :]
                conc_merit_2300 = nc_in_merit.sel(time = 30300)['solute concentration'].values[:, :]
                conc_merit_2400 = nc_in_merit.sel(time = 30400)['solute concentration'].values[:, :]
                conc_merit_2500 = nc_in_merit.sel(time = 30500)['solute concentration'].values[:, :]
        
                #   read in the extent inland 
                if top is None:
                    inl_ext_km = nc_in_merit['x'].values[0]
                    #   get the x and y coordinates
                    x_lst =  nc_in_merit.coords['x'].values
                    y_lst = np.arange(nc_in_merit.coords['y'].values[0], nc_in_merit.coords['y'].values[0] - conc_merit_2000.shape[0] * 10., -10.)
                    #y_lst =  nc_in.coords['y'].values
                    if x_lst[0] > 0.:
                        x_st = max(cst_idx * -0.1, - 10.05)
                        x_lst = np.arange(cst_idx * -0.1 + 0.05, cst_idx * -0.1 + 0.05 + len(x_lst) * 0.1, 0.1)
                    else:
                        x_st = max(x_lst[0], - 10.05)
                    x_end = 5.05
                    
                    y_top_lst = conc_merit_2000[:, 0].tolist()
                    st = 1
                    while len([x for x in y_top_lst if not math.isnan(x)]) == 0:
                        y_top_lst = conc_merit_2000[:, st].tolist()
                        st += 1
                    top = y_lst[y_top_lst.index(next(x for x in y_top_lst if not math.isnan(x)))] + 5.0
                        
                    #y_bot_lst = conc_merit_2000[:, int(abs(round(x_st * 10, 1)))].tolist()   
                    y_bot_lst = conc_merit_2000[:, cst_idx].tolist()
                    bot_idx_1 = y_bot_lst.index([x for x in y_bot_lst if not math.isnan(x)][0])
                    bot_idx_2 = len([x for x in y_bot_lst if not math.isnan(x)])
                    bot = y_lst[bot_idx_1 + bot_idx_2] + 5.0
                        
                    #top = y_lst[y_top_lst.index(next(x for x in y_top_lst if not math.isnan(x)))] + 5.0
                    #y_bot_lst = nc_in_coastal.sel(x = x_end)['solute concentration'].values[-1].tolist()
                    #bot = y_lst[y_bot_lst.index([x for x in y_bot_lst if not math.isnan(x)][-1])] + 5.0
                
                del nc_in_merit
        
            except (FileNotFoundError, OSError):
                """
                try:
                    conc_merit_2000 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30000.nc'))['solute concentration'].values[:, :]
                    conc_merit_2050 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30050.nc'))['solute concentration'].values[:, :]
                    conc_merit_2100 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30100.nc'))['solute concentration'].values[:, :]
                    conc_merit_2200 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30200.nc'))['solute concentration'].values[:, :]
                    conc_merit_2300 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30300.nc'))['solute concentration'].values[:, :]
                    conc_merit_2400 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30400.nc'))['solute concentration'].values[:, :]
                    conc_merit_2500 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30500.nc'))['solute concentration'].values[:, :]           
        
                    #   read in the extent inland 
                    inl_ext_km = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30000.nc'))['x'].values[0]
                    #   get the x and y coordinates
                    x_lst =  xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30000.nc')).coords['x'].values
                    y_lst = np.arange(xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30000.nc')).coords['y'].values[0],\
                                      xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30000.nc')).coords['y'].values[0] - conc_merit_2000.shape[0] * 10., -10.)
                    #y_lst =  nc_in.coords['y'].values
                    x_st = max(x_lst[0], - 10.05)
                    x_end = 5.05
                    y_top_lst = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30000.nc')).sel(x = x_st)['solute concentration'].values[-1].tolist()
                    top = y_lst[y_top_lst.index(next(x for x in y_top_lst if not math.isnan(x)))]
                    
                    y_bot_lst = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'merit', '_average_30000.nc')).sel(x = x_end)['solute concentration'].values[-1].tolist()
                    bot = y_lst[y_bot_lst.index(next(x for x in y_bot_lst if not math.isnan(x)))]
        
                except (FileNotFoundError, OSError):
                """
                conc_merit_2000 = None
                conc_merit_2050 = None
                conc_merit_2100 = None
                conc_merit_2200 = None
                conc_merit_2300 = None
                conc_merit_2400 = None
                conc_merit_2500 = None
        
            try:
                nc_in_gebco = xr.open_dataset(os.path.join(avg_nc_folder, cs_srm_str + '_' + rcp_name + '_gebco_MERGED.nc'))
                conc_gebco_2000 = nc_in_gebco.sel(time = 30000)['solute concentration'].values[:, :]
                conc_gebco_2050 = nc_in_gebco.sel(time = 30050)['solute concentration'].values[:, :]
                conc_gebco_2100 = nc_in_gebco.sel(time = 30100)['solute concentration'].values[:, :]
                conc_gebco_2200 = nc_in_gebco.sel(time = 30200)['solute concentration'].values[:, :]
                conc_gebco_2300 = nc_in_gebco.sel(time = 30300)['solute concentration'].values[:, :]
                conc_gebco_2400 = nc_in_gebco.sel(time = 30400)['solute concentration'].values[:, :]
                conc_gebco_2500 = nc_in_gebco.sel(time = 30500)['solute concentration'].values[:, :]
        
                #   read in the extent inland 
                if top is None:
                    inl_ext_km = nc_in_gebco['x'].values[0]
                    #   get the x and y coordinates
                    x_lst =  nc_in_gebco.coords['x'].values
                    y_lst = np.arange(nc_in_gebco.coords['y'].values[0], nc_in_gebco.coords['y'].values[0] - conc_gebco_2000.shape[0] * 10., -10.)
                    #y_lst =  nc_in.coords['y'].values
                    if x_lst[0] > 0.:
                        x_st = max(cst_idx * -0.1, - 10.05)
                        x_lst = np.arange(cst_idx * -0.1 + 0.05, cst_idx * -0.1 + 0.05 + len(x_lst) * 0.1, 0.1)
                    else:
                        x_st = max(x_lst[0], - 10.05)
                    x_end = 5.05
                    
                    y_top_lst = conc_gebco_2000[:, 0].tolist()
                    st = 1
                    while len([x for x in y_top_lst if not math.isnan(x)]) == 0:
                        y_top_lst = conc_gebco_2000[:, st].tolist()
                        st += 1
                    top = y_lst[y_top_lst.index(next(x for x in y_top_lst if not math.isnan(x)))] + 5.0
                        
                    y_bot_lst = conc_gebco_2000[:, cst_idx].tolist()
                    bot_idx_1 = y_bot_lst.index([x for x in y_bot_lst if not math.isnan(x)][0])
                    bot_idx_2 = len([x for x in y_bot_lst if not math.isnan(x)])
                    bot = y_lst[bot_idx_1 + bot_idx_2] + 5.0                    
                    #bot = y_lst[y_bot_lst.index([x for x in y_bot_lst if not math.isnan(x)][-1])] + 5.0
                    #top = y_lst[y_top_lst.index(next(x for x in y_top_lst if not math.isnan(x)))] + 5.0
                    #y_bot_lst = nc_in_coastal.sel(x = x_end)['solute concentration'].values[-1].tolist()
                    #bot = y_lst[y_bot_lst.index([x for x in y_bot_lst if not math.isnan(x)][-1])] + 5.0                del nc_in_gebco
                
            except (FileNotFoundError, OSError):
                """
                try:
                    conc_gebco_2000 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30000.nc'))['solute concentration'].values[:, :]
                    conc_gebco_2050 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30050.nc'))['solute concentration'].values[:, :]
                    conc_gebco_2100 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30100.nc'))['solute concentration'].values[:, :]
                    conc_gebco_2200 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30200.nc'))['solute concentration'].values[:, :]
                    conc_gebco_2300 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30300.nc'))['solute concentration'].values[:, :]
                    conc_gebco_2400 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30400.nc'))['solute concentration'].values[:, :]
                    conc_gebco_2500 = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30500.nc'))['solute concentration'].values[:, :]           
        
                    #   read in the extent inland 
                    inl_ext_km = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30000.nc'))['x'].values[0]
                    #   get the x and y coordinates
                    x_lst =  xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30000.nc')).coords['x'].values
                    y_lst = np.arange(xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30000.nc')).coords['y'].values[0],\
                                      xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30000.nc')).coords['y'].values[0] - conc_gebco_2000.shape[0] * 10., -10.)
                    #y_lst =  nc_in.coords['y'].values
                    x_st = max(x_lst[0], - 10.05)
                    x_end = 5.05
                    y_top_lst = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30000.nc')).sel(x = x_st)['solute concentration'].values[-1].tolist()
                    top = y_lst[y_top_lst.index(next(x for x in y_top_lst if not math.isnan(x)))]
                    
                    y_bot_lst = xr.open_dataset(os.path.join(sc_folder, 'avg_sim', '_nc', 'gebco', '_average_30000.nc')).sel(x = x_end)['solute concentration'].values[-1].tolist()
                    bot = y_lst[y_bot_lst.index(next(x for x in y_bot_lst if not math.isnan(x)))]
        
                except (FileNotFoundError, OSError):    
                """
                conc_gebco_2000 = None
                conc_gebco_2050 = None
                conc_gebco_2100 = None
                conc_gebco_2200 = None
                conc_gebco_2300 = None
                conc_gebco_2400 = None
                conc_gebco_2500 = None
         
            topo_lst_m = [round(i, 2) for i in topo_dict['gebco_merit_avg'].values.tolist()]
            topo_lst_g = [round(i, 2) for i in topo_dict['gebco_avg'].values.tolist()]
            topo_lst_c = [round(i, 2) for i in topo_dict['gebco_coastal_avg'].values.tolist()]
            mid_idx = int(2000 + cst_offset * 10)
            m_topo = topo_lst_m[mid_idx - min(int(abs(round(x_st, 1) * 10)), cst_idx) : mid_idx] + topo_lst_g[mid_idx : mid_idx + int(10 * round(x_end, 1))]
            c_topo = topo_lst_c[mid_idx - min(int(abs(round(x_st, 1) * 10)), cst_idx) : mid_idx] + topo_lst_g[mid_idx : mid_idx + int(10 * round(x_end, 1))]
            g_topo = topo_lst_g[mid_idx - min(int(abs(round(x_st, 1) * 10)), cst_idx) : mid_idx + int(10 * round(x_end, 1))] 
            
            """
            #   read and transofrm the CONC array 
            curr_conc = nc_in_coastal.sel(time = 30000)['solute concentration'].values
            curr_conc = (curr_conc * 1000).astype(int)
            curr_conc[curr_conc < -1] = -1
            curr_conc[curr_conc > 35000] = 35000
            """
        
            #   define the constants for the colobrar
            cmap = plt.cm.jet
            cmaplist = [cmap(i) for i in range(cmap.N)]
            cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
            # define the bins and normalize
            bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]            
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
            cbar_title = 'Salinity (ppt)'
            
            fig = plt.Figure(figsize = (11, 16))  
                                                                       
            #   first row - DEM names                                                                                                                                      
            ax1_1 = plt.subplot2grid((10, 5), (0, 1), fig = fig)             
            ax1_2 = plt.subplot2grid((10, 5), (0, 2), fig = fig)             
            ax1_3 = plt.subplot2grid((10, 5), (0, 3), fig = fig)             
        
            #   second till eitghth row - conc profiles
            ax2_1 = plt.subplot2grid((10, 5), (1, 1), fig = fig)             
            ax2_2 = plt.subplot2grid((10, 5), (1, 2), fig = fig)             
            ax2_3 = plt.subplot2grid((10, 5), (1, 3), fig = fig)          
        
            ax3_1 = plt.subplot2grid((10, 5), (2, 1), fig = fig)             
            ax3_2 = plt.subplot2grid((10, 5), (2, 2), fig = fig)             
            ax3_3 = plt.subplot2grid((10, 5), (2, 3), fig = fig)          
        
            ax4_1 = plt.subplot2grid((10, 5), (3, 1), fig = fig)             
            ax4_2 = plt.subplot2grid((10, 5), (3, 2), fig = fig)             
            ax4_3 = plt.subplot2grid((10, 5), (3, 3), fig = fig)          
        
            ax5_1 = plt.subplot2grid((10, 5), (4, 1), fig = fig)             
            ax5_2 = plt.subplot2grid((10, 5), (4, 2), fig = fig)             
            ax5_3 = plt.subplot2grid((10, 5), (4, 3), fig = fig)          
        
            ax6_1 = plt.subplot2grid((10, 5), (5, 1), fig = fig)             
            ax6_2 = plt.subplot2grid((10, 5), (5, 2), fig = fig)             
            ax6_3 = plt.subplot2grid((10, 5), (5, 3), fig = fig)          
        
            ax7_1 = plt.subplot2grid((10, 5), (6, 1), fig = fig)             
            ax7_2 = plt.subplot2grid((10, 5), (6, 2), fig = fig)             
            ax7_3 = plt.subplot2grid((10, 5), (6, 3), fig = fig)          
        
            ax8_1 = plt.subplot2grid((10, 5), (7, 1), fig = fig)             
            ax8_2 = plt.subplot2grid((10, 5), (7, 2), fig = fig)             
            ax8_3 = plt.subplot2grid((10, 5), (7, 3), fig = fig)   
            
            #   9th row is the x_axis label
            ax9_1 = plt.subplot2grid((10, 5), (8, 1), colspan = 3, fig = fig)         
        
            #   10th row is the legend and the DEM comparison graph
            ax10_1 = plt.subplot2grid((10, 5), (9, 1), fig = fig)         
            ax10_2 = plt.subplot2grid((10, 5), (9, 2), colspan = 2, fig = fig)         
            
            #   left column is the time label and the right column is the y axis label
            ax11_1 = plt.subplot2grid((10, 5), (1, 0), rowspan = 7, fig = fig)      
            ax11_2 = plt.subplot2grid((10, 5), (1, 4), rowspan = 7, fig = fig) 
            
            #   specify the location of each of these figure parts
            ax1_1.set_position([0.1, 0.95, 0.266, 0.05])                        # [left, bottom, width, height]
            ax1_2.set_position([0.366, 0.95, 0.266, 0.05])
            ax1_3.set_position([0.632, 0.95, 0.266, 0.05])
            
            ax2_1.set_position([0.1, 0.85, 0.266, 0.1])                      
            ax2_2.set_position([0.366, 0.85, 0.266, 0.1])
            ax2_3.set_position([0.632, 0.85, 0.266, 0.1])    
            
            ax3_1.set_position([0.1, 0.75, 0.266, 0.1])                      
            ax3_2.set_position([0.366, 0.75, 0.266, 0.1])
            ax3_3.set_position([0.632, 0.75, 0.266, 0.1])        
            
            ax4_1.set_position([0.1, 0.65, 0.266, 0.1])                      
            ax4_2.set_position([0.366, 0.65, 0.266, 0.1])
            ax4_3.set_position([0.632, 0.65, 0.266, 0.1])    
            
            ax5_1.set_position([0.1, 0.55, 0.266, 0.1])                      
            ax5_2.set_position([0.366, 0.55, 0.266, 0.1])
            ax5_3.set_position([0.632, 0.55, 0.266, 0.1])            
            
            ax6_1.set_position([0.1, 0.45, 0.266, 0.1])                      
            ax6_2.set_position([0.366, 0.45, 0.266, 0.1])
            ax6_3.set_position([0.632, 0.45, 0.266, 0.1])        
            
            ax7_1.set_position([0.1, 0.35, 0.266, 0.1])                      
            ax7_2.set_position([0.366, 0.35, 0.266, 0.1])
            ax7_3.set_position([0.632, 0.35, 0.266, 0.1])    
            
            ax8_1.set_position([0.1, 0.25, 0.266, 0.1])                      
            ax8_2.set_position([0.366, 0.25, 0.266, 0.1])
            ax8_3.set_position([0.632, 0.25, 0.266, 0.1])         
            
            ax9_1.set_position([0.1, 0.2, 0.8, 0.05])       
                       
            ax10_1.set_position([0.1, 0.05, 0.025, 0.15])
            ax10_2.set_position([0.25, 0.05, 0.65, 0.15])        
            
            ax11_1.set_position([0.025, 0.25, 0.07, 0.7])
            ax11_2.set_position([0.9, 0.25, 0.07, 0.7])      
            
            #   add the DEM strings 
            txt = ax1_1.text(0.5, 0.1, 'CoastalDEM', fontsize = 11, fontweight = 'bold', ha = 'center')
            txt.set_clip_on(False)
            ax1_1.axis('off')
            txt = ax1_2.text(0.5, 0.1, 'GEBCO ', fontsize = 11, fontweight = 'bold', ha = 'center')
            txt.set_clip_on(False)
            ax1_2.axis('off')    
            txt = ax1_3.text(0.5, 0.1, 'Merit', fontsize = 11, fontweight = 'bold', ha = 'center')
            txt.set_clip_on(False)
            ax1_3.axis('off')    
            txt = ax9_1.text(0.5, 0.3, 'Distance from current coastline (km)', fontsize = 11, fontweight = 'bold', ha = 'center')
            txt.set_clip_on(False)
            ax9_1.axis('off')  
            txt = ax11_1.text(0.25, 0.5, '2500 AD                 2400 AD                 2300 AD                 2200 AD                 2100 AD                 2050 AD                 2000 AD', fontsize = 11, fontweight = 'bold', rotation = 90, va = 'center')
            txt.set_clip_on(False)
            ax11_1.axis('off')    
            txt = ax11_2.text(0.75, 0.5, 'Elavation (m bsl)', fontsize = 11, fontweight = 'bold', rotation = 270, va = 'center')
            txt.set_clip_on(False)
            ax11_2.axis('off')   
        
            frsh_c_2000, frsh_g_2000, frsh_m_2000, frsh_c_2050, frsh_g_2050, frsh_m_2050, frsh_c_2100, frsh_g_2100, frsh_m_2100 = -1, -1, -1, -1, -1, -1, -1, -1, -1
            frsh_c_2200, frsh_g_2200, frsh_m_2200, frsh_c_2300, frsh_g_2300, frsh_m_2300, frsh_c_2400, frsh_g_2400, frsh_m_2400 = -1, -1, -1, -1, -1, -1, -1, -1, -1
            frsh_c_2500, frsh_g_2500, frsh_m_2500 = -1, -1, -1

            sww_c_2000, sww_g_2000, sww_m_2000, sww_c_2050, sww_g_2050, sww_m_2050, sww_c_2100, sww_g_2100, sww_m_2100 = -1, -1, -1, -1, -1, -1, -1, -1, -1
            sww_c_2200, sww_g_2200, sww_m_2200, sww_c_2300, sww_g_2300, sww_m_2300, sww_c_2400, sww_g_2400, sww_m_2400 = -1, -1, -1, -1, -1, -1, -1, -1, -1
            sww_c_2500, sww_g_2500, sww_m_2500 = -1, -1, -1
            
            #   draw the concentration profile for given time step
            if conc_coastal_2000 is not None and np.count_nonzero(conc_coastal_2000) != 0:
                im2_1 = ax2_1.imshow(conc_coastal_2000, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst)  + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax2_1.set_xlim([x_st, x_end])
                ax2_1.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_c_2000 = fresh_inland(conc_coastal_2000, x_lst)
                sww_c_2000 = sww_inland(conc_coastal_2000, cst_idx)
            
            if conc_gebco_2000 is not None and np.count_nonzero(conc_gebco_2000) != 0:
                im2_2 = ax2_2.imshow(conc_gebco_2000, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax2_2.set_xlim([x_st, x_end])
                ax2_2.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_g_2000 = fresh_inland(conc_gebco_2000, x_lst)
                sww_g_2000 = sww_inland(conc_gebco_2000, cst_idx)
        
            if conc_merit_2000 is not None and np.count_nonzero(conc_merit_2000) != 0:    
                im2_3 = ax2_3.imshow(conc_merit_2000, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax2_3.set_xlim([x_st, x_end])
                ax2_3.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_m_2000 = fresh_inland(conc_merit_2000, x_lst)
                sww_m_2000 = sww_inland(conc_merit_2000, cst_idx)
        
            if conc_coastal_2050 is not None and np.count_nonzero(conc_coastal_2050) != 0:   
                im3_1 = ax3_1.imshow(conc_coastal_2050, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax3_1.set_xlim([x_st, x_end])
                ax3_1.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_c_2050 = fresh_inland(conc_coastal_2050, x_lst)
                sww_c_2050 = sww_inland(conc_coastal_2050, cst_idx)
                
            if conc_gebco_2050 is not None and np.count_nonzero(conc_gebco_2050) != 0:       
                im3_2 = ax3_2.imshow(conc_gebco_2050, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax3_2.set_xlim([x_st, x_end])
                ax3_2.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_g_2050 = fresh_inland(conc_gebco_2050, x_lst)
                sww_g_2050 = sww_inland(conc_gebco_2050, cst_idx)
                
            if conc_merit_2050 is not None and np.count_nonzero(conc_merit_2050) != 0:  
                im3_3 = ax3_3.imshow(conc_merit_2050, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax3_3.set_xlim([x_st, x_end])
                ax3_3.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_m_2050 = fresh_inland(conc_merit_2050, x_lst)
                sww_m_2050 = sww_inland(conc_merit_2050, cst_idx)
                
            if conc_coastal_2100 is not None and np.count_nonzero(conc_coastal_2100) != 0:  
                im4_1 = ax4_1.imshow(conc_coastal_2100, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax4_1.set_xlim([x_st, x_end])
                ax4_1.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_c_2100 = fresh_inland(conc_coastal_2100, x_lst)
                sww_c_2100 = sww_inland(conc_coastal_2100, cst_idx)
                
            if conc_gebco_2100 is not None and np.count_nonzero(conc_gebco_2100) != 0:      
                im4_2 = ax4_2.imshow(conc_gebco_2100, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax4_2.set_xlim([x_st, x_end])
                ax4_2.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])  
                frsh_g_2100 = fresh_inland(conc_gebco_2100, x_lst)
                sww_g_2100 = sww_inland(conc_gebco_2100, cst_idx)
                
            if conc_merit_2100 is not None and np.count_nonzero(conc_merit_2100) != 0:  
                im4_3 = ax4_3.imshow(conc_merit_2100, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax4_3.set_xlim([x_st, x_end])
                ax4_3.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_m_2100 = fresh_inland(conc_merit_2100, x_lst) 
                sww_m_2100 = sww_inland(conc_merit_2100, cst_idx)
               
            if conc_coastal_2200 is not None and np.count_nonzero(conc_coastal_2200) != 0:  
                im5_1 = ax5_1.imshow(conc_coastal_2200, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax5_1.set_xlim([x_st, x_end])
                ax5_1.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_c_2200 = fresh_inland(conc_coastal_2200, x_lst) 
                sww_c_2200 = sww_inland(conc_coastal_2200, cst_idx)
                
            if conc_gebco_2200 is not None and np.count_nonzero(conc_gebco_2200) != 0:      
                im5_2 = ax5_2.imshow(conc_gebco_2200, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax5_2.set_xlim([x_st, x_end])
                ax5_2.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])    
                frsh_g_2200 = fresh_inland(conc_gebco_2200, x_lst) 
                sww_g_2200 = sww_inland(conc_gebco_2200, cst_idx)
                    
            if conc_merit_2200 is not None and np.count_nonzero(conc_merit_2200) != 0:  
                im5_3 = ax5_3.imshow(conc_merit_2200, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax5_3.set_xlim([x_st, x_end])
                ax5_3.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_m_2200 = fresh_inland(conc_merit_2200, x_lst) 
                sww_m2200 = sww_inland(conc_merit_2200, cst_idx)
                
            if conc_coastal_2300 is not None and np.count_nonzero(conc_coastal_2300) != 0:  
                im6_1 = ax6_1.imshow(conc_coastal_2300, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax6_1.set_xlim([x_st, x_end])
                ax6_1.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_c_2300 = fresh_inland(conc_coastal_2300, x_lst) 
                sww_c_2300 = sww_inland(conc_coastal_2300, cst_idx)
                
            if conc_gebco_2300 is not None and np.count_nonzero(conc_gebco_2300) != 0:      
                im6_2 = ax6_2.imshow(conc_gebco_2300, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax6_2.set_xlim([x_st, x_end])
                ax6_2.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_g_2300 = fresh_inland(conc_gebco_2300, x_lst) 
                sww_g_2300 = sww_inland(conc_gebco_2300, cst_idx)
                
            if conc_merit_2300 is not None and np.count_nonzero(conc_merit_2300) != 0:  
                im6_3 = ax6_3.imshow(conc_merit_2300, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax6_3.set_xlim([x_st, x_end])
                ax6_3.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])   
                frsh_m_2300 = fresh_inland(conc_merit_2300, x_lst) 
                sww_m_2300 = sww_inland(conc_merit_2300, cst_idx)
                
            if conc_coastal_2400 is not None and np.count_nonzero(conc_coastal_2400) != 0:  
                im7_1 = ax7_1.imshow(conc_coastal_2400, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax7_1.set_xlim([x_st, x_end])
                ax7_1.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_c_2400 = fresh_inland(conc_coastal_2400, x_lst)
                sww_c_2400 = sww_inland(conc_coastal_2400, cst_idx)
                
            if conc_gebco_2400 is not None and np.count_nonzero(conc_gebco_2400) != 0:      
                im7_2 = ax7_2.imshow(conc_gebco_2400, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax7_2.set_xlim([x_st, x_end])
                ax7_2.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_g_2400 = fresh_inland(conc_gebco_2400, x_lst)
                sww_g_2400 = sww_inland(conc_gebco_2400, cst_idx)
                
            if conc_merit_2400 is not None and np.count_nonzero(conc_merit_2400) != 0:          
                im7_3 = ax7_3.imshow(conc_merit_2400, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax7_3.set_xlim([x_st, x_end])
                ax7_3.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])  
                frsh_m_2400 = fresh_inland(conc_merit_2400, x_lst) 
                sww_m_2400 = sww_inland(conc_merit_2400, cst_idx)
                
            if conc_coastal_2500 is not None and np.count_nonzero(conc_coastal_2500) != 0:  
                im8_1 = ax8_1.imshow(conc_coastal_2500, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax8_1.set_xlim([x_st, x_end])
                ax8_1.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_c_2500 = fresh_inland(conc_coastal_2500, x_lst) 
                sww_c_2500 = sww_inland(conc_coastal_2500, cst_idx)
                
            if conc_gebco_2500 is not None and np.count_nonzero(conc_gebco_2500) != 0:      
                im8_2 = ax8_2.imshow(conc_gebco_2500, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax8_2.set_xlim([x_st, x_end])
                ax8_2.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])
                frsh_g_2500 = fresh_inland(conc_gebco_2500, x_lst) 
                sww_g_2500 = sww_inland(conc_gebco_2500, cst_idx)
                
            if conc_merit_2500 is not None and np.count_nonzero(conc_merit_2500) != 0:          
                im8_3 = ax8_3.imshow(conc_merit_2500, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), vmin = 0, vmax = 35.0, animated = True)
                #   set the limits to the concentration plots
                ax8_3.set_xlim([x_st, x_end])
                ax8_3.set_ylim([bot - 25.0 + 5.0, max(y_lst) + 5.0])        
                frsh_m_2500 = fresh_inland(conc_merit_2500, x_lst) 
                sww_m_2500 = sww_inland(conc_merit_2500, cst_idx)
                
            #x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 5.) * 5., math.ceil(x_end), 5.0)[1:], nbins = None)    
            #x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 1.) * 1., math.ceil(x_end), 1.0), nbins = None)    
            #y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 50.) * 50., math.ceil(top), 50.0)[1:], nbins = None)    
            #y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 10.) * 10. , math.ceil(top), 10.0), nbins = None)    
            
            y_top = 50 * round(float(math.ceil(top)) / 50) + 50 + 5.0
            
            x_major_ticks = np.arange(math.floor(x_st / 5.) * 5., math.ceil(x_end), 5.0)[1:]
            x_minor_ticks = np.arange(math.floor(x_st / 1.) * 1., math.ceil(x_end), 1.0)
            y_major_ticks = np.arange(math.floor(bot / 50.) * 50., y_top + 1.0, 50.0)[1:]
            y_minor_ticks = np.arange(math.floor(bot / 10.) * 10., y_top + 1.0, 10.0) 
        
            ax2_1.set_xticks(x_major_ticks)
            ax2_1.set_xticks(x_minor_ticks, minor=True)
            ax2_1.set_yticks(y_major_ticks)
            ax2_1.set_yticks(y_minor_ticks, minor=True)      
            ax2_1.grid(which='minor', alpha=0.2)
            ax2_1.grid(which='major', alpha=0.5)
        
            ax2_2.set_xticks(x_major_ticks)
            ax2_2.set_xticks(x_minor_ticks, minor=True)
            ax2_2.set_yticks(y_major_ticks)
            ax2_2.set_yticks(y_minor_ticks, minor=True)      
            ax2_2.grid(which='minor', alpha=0.2)
            ax2_2.grid(which='major', alpha=0.5)
            
            ax2_3.set_xticks(x_major_ticks)
            ax2_3.set_xticks(x_minor_ticks, minor=True)
            ax2_3.set_yticks(y_major_ticks)
            ax2_3.set_yticks(y_minor_ticks, minor=True)      
            ax2_3.grid(which='minor', alpha=0.2)
            ax2_3.grid(which='major', alpha=0.5)
        
            ax3_1.set_xticks(x_major_ticks)
            ax3_1.set_xticks(x_minor_ticks, minor=True)
            ax3_1.set_yticks(y_major_ticks)
            ax3_1.set_yticks(y_minor_ticks, minor=True)      
            ax3_1.grid(which='minor', alpha=0.2)
            ax3_1.grid(which='major', alpha=0.5)
        
            ax3_2.set_xticks(x_major_ticks)
            ax3_2.set_xticks(x_minor_ticks, minor=True)
            ax3_2.set_yticks(y_major_ticks)
            ax3_2.set_yticks(y_minor_ticks, minor=True)      
            ax3_2.grid(which='minor', alpha=0.2)
            ax3_2.grid(which='major', alpha=0.5)
            
            ax3_3.set_xticks(x_major_ticks)
            ax3_3.set_xticks(x_minor_ticks, minor=True)
            ax3_3.set_yticks(y_major_ticks)
            ax3_3.set_yticks(y_minor_ticks, minor=True)      
            ax3_3.grid(which='minor', alpha=0.2)
            ax3_3.grid(which='major', alpha=0.5)
        
            ax4_1.set_xticks(x_major_ticks)
            ax4_1.set_xticks(x_minor_ticks, minor=True)
            ax4_1.set_yticks(y_major_ticks)
            ax4_1.set_yticks(y_minor_ticks, minor=True)      
            ax4_1.grid(which='minor', alpha=0.2)
            ax4_1.grid(which='major', alpha=0.5)
        
            ax4_2.set_xticks(x_major_ticks)
            ax4_2.set_xticks(x_minor_ticks, minor=True)
            ax4_2.set_yticks(y_major_ticks)
            ax4_2.set_yticks(y_minor_ticks, minor=True)      
            ax4_2.grid(which='minor', alpha=0.2)
            ax4_2.grid(which='major', alpha=0.5)
            
            ax4_3.set_xticks(x_major_ticks)
            ax4_3.set_xticks(x_minor_ticks, minor=True)
            ax4_3.set_yticks(y_major_ticks)
            ax4_3.set_yticks(y_minor_ticks, minor=True)      
            ax4_3.grid(which='minor', alpha=0.2)
            ax4_3.grid(which='major', alpha=0.5)
        
            ax5_1.set_xticks(x_major_ticks)
            ax5_1.set_xticks(x_minor_ticks, minor=True)
            ax5_1.set_yticks(y_major_ticks)
            ax5_1.set_yticks(y_minor_ticks, minor=True)      
            ax5_1.grid(which='minor', alpha=0.2)
            ax5_1.grid(which='major', alpha=0.5)
        
            ax5_2.set_xticks(x_major_ticks)
            ax5_2.set_xticks(x_minor_ticks, minor=True)
            ax5_2.set_yticks(y_major_ticks)
            ax5_2.set_yticks(y_minor_ticks, minor=True)      
            ax5_2.grid(which='minor', alpha=0.2)
            ax5_2.grid(which='major', alpha=0.5)
            
            ax5_3.set_xticks(x_major_ticks)
            ax5_3.set_xticks(x_minor_ticks, minor=True)
            ax5_3.set_yticks(y_major_ticks)
            ax5_3.set_yticks(y_minor_ticks, minor=True)      
            ax5_3.grid(which='minor', alpha=0.2)
            ax5_3.grid(which='major', alpha=0.5)
        
            ax6_1.set_xticks(x_major_ticks)
            ax6_1.set_xticks(x_minor_ticks, minor=True)
            ax6_1.set_yticks(y_major_ticks)
            ax6_1.set_yticks(y_minor_ticks, minor=True)      
            ax6_1.grid(which='minor', alpha=0.2)
            ax6_1.grid(which='major', alpha=0.5)
        
            ax6_2.set_xticks(x_major_ticks)
            ax6_2.set_xticks(x_minor_ticks, minor=True)
            ax6_2.set_yticks(y_major_ticks)
            ax6_2.set_yticks(y_minor_ticks, minor=True)      
            ax6_2.grid(which='minor', alpha=0.2)
            ax6_2.grid(which='major', alpha=0.5)
            
            ax6_3.set_xticks(x_major_ticks)
            ax6_3.set_xticks(x_minor_ticks, minor=True)
            ax6_3.set_yticks(y_major_ticks)
            ax6_3.set_yticks(y_minor_ticks, minor=True)      
            ax6_3.grid(which='minor', alpha=0.2)
            ax6_3.grid(which='major', alpha=0.5)
        
            ax7_1.set_xticks(x_major_ticks)
            ax7_1.set_xticks(x_minor_ticks, minor=True)
            ax7_1.set_yticks(y_major_ticks)
            ax7_1.set_yticks(y_minor_ticks, minor=True)      
            ax7_1.grid(which='minor', alpha=0.2)
            ax7_1.grid(which='major', alpha=0.5)
        
            ax7_2.set_xticks(x_major_ticks)
            ax7_2.set_xticks(x_minor_ticks, minor=True)
            ax7_2.set_yticks(y_major_ticks)
            ax7_2.set_yticks(y_minor_ticks, minor=True)      
            ax7_2.grid(which='minor', alpha=0.2)
            ax7_2.grid(which='major', alpha=0.5)
            
            ax7_3.set_xticks(x_major_ticks)
            ax7_3.set_xticks(x_minor_ticks, minor=True)
            ax7_3.set_yticks(y_major_ticks)
            ax7_3.set_yticks(y_minor_ticks, minor=True)      
            ax7_3.grid(which='minor', alpha=0.2)
            ax7_3.grid(which='major', alpha=0.5)
        
            ax8_1.set_xticks(x_major_ticks)
            ax8_1.set_xticks(x_minor_ticks, minor=True)
            ax8_1.set_yticks(y_major_ticks)
            ax8_1.set_yticks(y_minor_ticks, minor=True)      
            ax8_1.grid(which='minor', alpha=0.2)
            ax8_1.grid(which='major', alpha=0.5)
        
            ax8_2.set_xticks(x_major_ticks)
            ax8_2.set_xticks(x_minor_ticks, minor=True)
            ax8_2.set_yticks(y_major_ticks)
            ax8_2.set_yticks(y_minor_ticks, minor=True)      
            ax8_2.grid(which='minor', alpha=0.2)
            ax8_2.grid(which='major', alpha=0.5)
            
            ax8_3.set_xticks(x_major_ticks)
            ax8_3.set_xticks(x_minor_ticks, minor=True)
            ax8_3.set_yticks(y_major_ticks)
            ax8_3.set_yticks(y_minor_ticks, minor=True)      
            ax8_3.grid(which='minor', alpha=0.2)
            ax8_3.grid(which='major', alpha=0.5)
              
            ax2_3.yaxis.tick_right()
            ax2_3.yaxis.set_label_position("right")
            ax3_3.yaxis.tick_right()
            ax3_3.yaxis.set_label_position("right")
            ax4_3.yaxis.tick_right()
            ax4_3.yaxis.set_label_position("right")
            ax5_3.yaxis.tick_right()
            ax5_3.yaxis.set_label_position("right")
            ax6_3.yaxis.tick_right()
            ax6_3.yaxis.set_label_position("right")
            ax7_3.yaxis.tick_right()
            ax7_3.yaxis.set_label_position("right")
            ax8_3.yaxis.tick_right()
            ax8_3.yaxis.set_label_position("right")
        
            ax2_1.tick_params(labelleft = False, labelbottom = False)   
            ax2_2.tick_params(labelleft = False, labelbottom = False) 
            ax2_3.tick_params(labelbottom = False)   
            ax3_1.tick_params(labelleft = False, labelbottom = False)   
            ax3_2.tick_params(labelleft = False, labelbottom = False)   
            ax3_3.tick_params(labelbottom = False)   
            ax4_1.tick_params(labelleft = False, labelbottom = False)   
            ax4_2.tick_params(labelleft = False, labelbottom = False)   
            ax4_3.tick_params(labelbottom = False)   
            ax5_1.tick_params(labelleft = False, labelbottom = False)   
            ax5_2.tick_params(labelleft = False, labelbottom = False)   
            ax5_3.tick_params(labelbottom = False)   
            ax6_1.tick_params(labelleft = False, labelbottom = False)   
            ax6_2.tick_params(labelleft = False, labelbottom = False)   
            ax6_3.tick_params(labelbottom = False)   
            ax7_1.tick_params(labelleft = False, labelbottom = False)   
            ax7_2.tick_params(labelleft = False, labelbottom = False)   
            ax7_3.tick_params(labelbottom = False)   
            ax8_1.tick_params(labelleft = False)   
            ax8_2.tick_params(labelleft = False)    
            
            #   plot the colorbar
            cbar = plt.colorbar(im6_2, cax = ax10_1, cmap = cmap, pad=-0.2,  shrink=0.9, aspect=30, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
            #cbar.ax.set_title(cbar_title, fontsize = 9, y  = 1.025)
            cbar.ax.set_ylabel(cbar_title, rotation = 90)
            cbar.ax.yaxis.set_label_position('left')
            cbar.ax.tick_params(labelsize = 9)    
        
            #   plot the different DEM elevations
            if conc_merit_2000 is not None:
                ibound_arr = conc_merit_2000
                ibound_arr[ibound_arr < 40] = 0
            elif conc_gebco_2000 is not None:
                ibound_arr = conc_gebco_2000
                ibound_arr[ibound_arr < 40] = 0
            elif conc_coastal_2000 is not None:
                ibound_arr = conc_coastal_2000
                ibound_arr[ibound_arr < 40] = 0
        
            cmap = matplotlib.colors.ListedColormap(['silver'])
            bounds=[0,5]
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        
            im5_3 = ax10_2.imshow(ibound_arr, aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                       extent = (x_lst[0], x_lst[-1], min(y_lst) + 5.0, max(y_lst) + 5.0), animated = True, alpha = 0.5)
            
            #   check that the dimenstions are the same
            x_lst_topo = [round(x, 2) for x in np.arange(x_st, x_st + len(m_topo) * 0.1, 0.1).tolist()]
            if len(x_lst_topo) > len(m_topo):
                x_lst_topo = x_lst_topo[:len(m_topo)]

            ax10_2.plot(x_lst_topo, m_topo, lw = 1.5, color = 'forestgreen', label = 'MERIT')
            ax10_2.plot(x_lst_topo, c_topo, lw = 1.5, color = 'royalblue', label = 'CoastalDEM')
            ax10_2.plot(x_lst_topo, g_topo, lw = 1.5, color = 'red', label = 'GEBCO')

            #   set the limits to the concentration plots
            ax10_2.set_xlim([x_st, 1.])
            ax10_2.set_ylim([- 25.0 + 5.0, max(y_lst) + 5.0])
            #ax10_2.set_xticks(x_major_ticks)
            #ax10_2.set_xticks(x_minor_ticks, minor=True)
            #ax10_2.set_yticks(y_major_ticks)
            #ax10_2.set_yticks(y_minor_ticks, minor=True)  
            ax10_2.set_xlabel('Distance from current coastline (km)')  
            ax10_2.set_ylabel('Elavation (m bsl)')  
            ax10_2.legend(loc="lower left")
            ax10_2.grid(which='minor', alpha=0.2)
            ax10_2.grid(which='major', alpha=0.5)
            
            canvas = FigureCanvas(fig)
            canvas.print_figure(os.path.join(out_dir, '_conc_plots_files', cs_srm_str + '_' + rcp_name + '_conc_results.png'), dpi=500)
        
            #   get all the pct fresh for the RCP scenario
            df_in = pd.read_csv(os.path.join(out_dir, rcp_name + '_results.csv'))
            try:
                id_loop_row = df_in.loc[df_in[' cs_srm_id'] == cs_srm_id]    
            except KeyError:
                del df_in
                #   create the output CSV file, the headers will match the headers of the usual individual csv file
                f = open(os.path.join(out_dir, rcp_name + '_results.csv'),'a')
                #f.write("\n")
                f.write(' ,' + str(coscat_str) + ',' + str(srm_str) + ',' + cs_srm_id + ',' + str(x_st) + ',' + str(frsh_m_2000) + ',' + str(frsh_c_2000)\
                         + ',' + str(frsh_g_2000) + ',' + str(sww_m_2000) + ',' + str(sww_c_2000) + ',' + str(sww_g_2000)\
                         + ',' + str(frsh_m_2050) + ',' + str(frsh_c_2050) + ',' + str(frsh_g_2050) + ',' + str(sww_m_2050) + ',' + str(sww_c_2050) + ',' + str(sww_g_2050)\
                         + ',' + str(frsh_m_2100) + ',' + str(frsh_c_2100) + ',' + str(frsh_g_2100) + ',' + str(sww_m_2100) + ',' + str(sww_c_2100) + ',' + str(sww_g_2100)\
                         + ',' + str(frsh_m_2200) + ',' + str(frsh_c_2200) + ',' + str(frsh_g_2200) + ',' + str(sww_m_2200) + ',' + str(sww_c_2200) + ',' + str(sww_g_2200)\
                         + ',' + str(frsh_m_2300) + ',' + str(frsh_c_2300) + ',' + str(frsh_g_2300) + ',' + str(sww_m_2300) + ',' + str(sww_c_2300) + ',' + str(sww_g_2300)\
                         + ',' + str(frsh_m_2400) + ',' + str(frsh_c_2400) + ',' + str(frsh_g_2400) + ',' + str(sww_m_2400) + ',' + str(sww_c_2400) + ',' + str(sww_g_2400)\
                         + ',' + str(frsh_m_2500) + ',' + str(frsh_c_2500) + ',' + str(frsh_g_2500) + ',' + str(sww_m_2500) + ',' + str(sww_c_2500) + ',' + str(sww_g_2500))
                f.write("\n")
                f.close()                         

    except FileNotFoundError:
        print(id_loop)
        pass
    
    
    
    
    
    
    
    
    
    
    
