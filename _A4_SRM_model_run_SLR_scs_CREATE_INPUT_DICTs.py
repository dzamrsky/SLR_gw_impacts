import numpy as np
import os
import collections
import time
import gc
import flopy
import random
import sys

import math
#import geopandas as gpd
import xarray as xr
import pandas as pd
import statistics
from statistics import StatisticsError
import geopandas as gpd

#import scipy
import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.ioff()
#from scipy.stats import norm
#   the id loop is the only input into the script, the rest will be found in the lookup table
#id_loop = int(sys.argv[1])
#param_combo = int(sys.argv[2])
id_loop = 1986
#param_combo = 1

"""   -----------------------------------------------------------------------------------   """
"""               define all paths and directories, load additional packages                """ 
"""   -----------------------------------------------------------------------------------   """

#main_dir = r'g:\Water_Nexus\_A4_GUM\_SRM_models'
master_csv_dir = r'g:\Water_Nexus\_A4_GUM\_coscat_reg_main_tb_MODELRUNS_old.csv' 
model_sc_csv_dir = r'g:\Water_Nexus\_A4_models\_SLR_models\_model_scenarios_GEO.csv'
input_dir = r'g:\Water_Nexus\_A4_MODEL_input_files'
swat_exe_dir = r'g:\Water_Nexus\Modelling\swt_v4_00_05\exe\swt_v4x64.exe' 
scripts_dir = r'g:\Water_Nexus\_A4_GUM\_scripts\_fc_scripts'
#gis_data_dir =  r'g:\Water_Nexus\_A4_GUM\_GIS_data'
cs_regs_shp_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\cs_reg_id_v2.shp'
model_out_dir = r'g:\Water_Nexus\_A4_models\_SLR_models'
init_dir = r'g:\Water_Nexus\_A4\_COSCAT_input_files'

"""
master_csv_dir = r'/gpfs/work4/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models/_coscat_reg_main_tb_MODELRUNS.csv'
model_sc_csv_dir = r'/gpfs/work4/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models/_model_scenarios_GEO.csv'
input_dir = r'/gpfs/work4/0/qt16165/_dzamrsky/_A4_MODEL_input_files'
cs_regs_shp_dir = r'/gpfs/work4/0/qt16165/_dzamrsky/_A4_input_data/cs_reg_id_v2.shp'
model_out_dir = r'/gpfs/work4/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models'
init_dir = r'/gpfs/work4/0/qt16165/_dzamrsky/A3_final_results/coscat_files'
swat_exe_dir = r'/home/dzamrsky/swtv4'
scripts_dir = r'/home/dzamrsky/_A4/_fc_scripts'
"""
sys.path.append(scripts_dir)  
#   add the path with all additional tailor made functions to import the scripts
import model_geotools_py3 as mgt
import modelfunctions_py3 as mfc
import model_geotools_fresh_py3 as mgt_frsh
from modeltools_py3 import cs_model 
import _create_ARP_RIVBAS_tools as rivbas
#import ws_masterscript_py3 as ws_ms
#import model_plotting_tools as mpt

#   define the location of the SEAWAT executable
mf_exe_dir = swat_exe_dir
mt3d_exe_dir = swat_exe_dir


#   function that fills missing elements in lists.. who the fuck knows what went wrong when creating these lists in the first place
def fix_lst(in_lst):
    #   make sure all the lists are of same length
    lst = []
    for j in range(len(in_lst)):
        lst.append(len(in_lst[j]))
    max_lst_len = max(lst)
    for j in range(len(in_lst)):
        while len(in_lst[j]) < max_lst_len:
            in_lst[j].append(in_lst[j][-1])
    lst_new = []
    for j in range(len(in_lst)):
        lst_new.append(in_lst[j])
    out_lst = lst_new
    return out_lst


for param_combo in [1, 9, 17]:

    #   get the region and coscat ID numbers
    df_in = pd.read_csv(master_csv_dir)
    id_loop_row = df_in.loc[df_in['id_loop'] == id_loop]
    coscat_id = int(id_loop_row['coscat_id'].values[0])
    subreg_id = int(id_loop_row['cst_reg_id'].values[0])
    print('Running model for : ' + '\n' + '    id_loop   ' + str(id_loop) + '\n' + '    scenario  ' + str(param_combo)\
          + '\n' + '    COSCAT    ' + str(coscat_id) + '\n' + '    SRM       ' + str(subreg_id))
    
    df_model_sc_in = pd.read_csv(model_sc_csv_dir)
    id_model_sc_row = df_model_sc_in.loc[df_model_sc_in['sc'] == 'sc_' + str(param_combo)]
    
    #   change the status of the model to -1, the values are as follows
    #       0  - not finished nor started
    #       -1 - not finished, but started
    #       1  - finished
    status = int(df_in.loc[df_in['id_loop'] == id_loop, 'sc_' + str(param_combo)].values[0])
    
    if status != 1:
        df_in.loc[df_in['id_loop'] == id_loop, 'sc_' + str(param_combo)] = -1
        #df_in.at[id_loop_row.index[0], 'sc_' + str(param_combo)] = -1
        df_in.update(df_in)
        df_in.to_csv(master_csv_dir, sep = ',', encoding = 'utf-8', index = False)
    
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
    
    #   read in the original raw shapefile with coastal points and coastal types, and check the coastal type of the SRM
    cs_regs_shp_dir = os.path.join(cs_regs_shp_dir)
    cs_regs_raw = gpd.read_file(cs_regs_shp_dir)
    cs_type = cs_regs_raw.loc[cs_regs_raw['cst_reg_id'] == subreg_id]
    cs_type = cs_type.loc[cs_type['coscat_id'] == coscat_id]['cst_type'].values.tolist()
    cs_type = max(set(cs_type), key=cs_type.count)
    
    #   1 means delta
    #   2 means simple
    #   3 means island
    #   first check if it is delta type, if the delta folder doesnt exist, then it must be called _all_
    if cs_type == 1:
        if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_delta', 'avg_sim', 'avg_sim_MERGED.nc')):
            init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_delta', 'avg_sim', 'avg_sim_MERGED.nc')
        else:
            if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_all', 'avg_sim', 'avg_sim_MERGED.nc')):
                init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_all', 'avg_sim', 'avg_sim_MERGED.nc')
            else:
                if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim', 'avg_sim_MERGED.nc')):
                    init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim', 'avg_sim_MERGED.nc')
                else:
                    if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_other', 'avg_sim', 'avg_sim_MERGED.nc')):
                        init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_other', 'avg_sim', 'avg_sim_MERGED.nc')
                    else:
                        init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry', 'avg_sim', 'avg_sim_MERGED.nc')
    
    elif cs_type == 2:
        if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_henry', 'avg_sim', 'avg_sim_MERGED.nc')):
            init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry', 'avg_sim', 'avg_sim_MERGED.nc')
        else:
            if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim', 'avg_sim_MERGED.nc')):
                init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim', 'avg_sim_MERGED.nc')
            else:
                init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_all', 'avg_sim', 'avg_sim_MERGED.nc')
    
    elif cs_type == 3:
        if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_other', 'avg_sim', 'avg_sim_MERGED.nc')):
            init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_other', 'avg_sim', 'avg_sim_MERGED.nc')
        else:
            if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim', 'avg_sim_MERGED.nc')):
                init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim', 'avg_sim_MERGED.nc')
            else:
                init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_all', 'avg_sim', 'avg_sim_MERGED.nc')
    
    print(init_dict_dir)
    
    """   -----------------------------------------------------------------------------------   """
    """               Define the model constants and also the SEWAT directory                   """ 
    """   -----------------------------------------------------------------------------------   """
    
    cs_points_dist = 500.
    max_depth = -2000. 
    
    #   Specify the COSCAT ID and also the model types to be modelled
    model_sc_type = 'salt'
    del_col = 100.  #   the initial column width resolution to be applied during the initial model build
    del_lay = 10.   #   the initial layer thickness resolution to be applied during the initial model build
    
    topo_type = id_model_sc_row['topo'].values[0]   
    topo_100_type = topo_type # id_loop_row['topo_100'].values[0]   
    #del_col_res = id_model_sc_row['del_col'].values[0]   
    #del_lay_res = id_model_sc_row['del_lay'].values[0]   
    #rch_type = id_model_sc_row['gw_rch'].values[0]   
    #ld_rate = id_model_sc_row['LD'].values[0]   
    
    sand_lay_n = id_model_sc_row['sand_lay_n'].values[0]   
    p_fact = id_model_sc_row['p_fact'].values[0]   
    clay_cap_thk = id_model_sc_row['clay_cap_thk'].values[0]   
    
    #slr_sc = id_model_sc_row['slr'].values[0]   
    del_col_res = 100
    del_lay_res = 10
    rch_type = 'paleo'
    ld_rate = 1
    
    #del_col = del_col_res
    #del_lay = del_lay_res
    
    print('COSCAT_id = ' + str(coscat_id) + '  |  ' + 'SUBREG_id = ' + str(subreg_id), ' |  st_cond = ' + model_sc_type)
    print('Model resolution :  column width = ' + str(del_col_res) + '  |   layer thickness = ' + str(del_lay_res))
    print('Recharge typae = ' + rch_type)
    
    #  define all the inputs for the DIS package
    nrow = 1                                    # number of rows in model domain (2D so = 1)
    delc = 1.0                                  # an array of spacings along a column, using the list created earlier (2D so = 1)
    laycbd = 0                                  # 0 indicates no confining bed... -> this needs to be checked!!
    #   define input for the LPF package
    laytyp_val = 0
    #   PCG
    hclose = 1e-4 
    rclose = 1e+1
    #   OC
    th_oc_time_step_val = 1
    #   BTN
    porosity = 0.3
    dt0 = 365.25 #50.
    nprs = 1
    ifmtcn = 0
    chkmas = False
    nprmas = 10
    nprobs = 10
    #   ADV
    mixelm = 0
    #   DSP
    dmcoef = 0.0000864    # effective molecular diffusion coefficient [M2/D]
    al = 1.
    trpt = 0.1
    trpv = 0.1
    #   GCG
    iter1 = 500
    mxiter = 1
    isolve = 1
    cclose = 1e-7
    #   VDF
    iwtable = 0
    densemin = 1000.
    densemax = 1025.
    denseref = 1000.
    denseslp = 0.7143
    firstdt = 0.001
    
    """   -----------------------------------------------------------------------------------   """
    """              Read in the scanrio parameter values and the topography                    """ 
    """   -----------------------------------------------------------------------------------   """
    if del_col_res >= 100:
        del_col_res_str = str(del_col_res)
    else:
        del_col_res_str = '0' + str(del_col_res)
    
    subreg_input_dir = os.path.join(input_dir, id_cs_str, '%s') % (id_cs_str + '_REG_' + str_id_subreg)
    topo_nc_dir = os.path.join(subreg_input_dir, '_topo_%sm_avg.nc') % (del_col_res)
    topo_500m_nc_dir = os.path.join(subreg_input_dir, '_topo_500m_avg.nc')
    #ibound_dict_dir = os.path.join(subreg_input_dir, '_%s_avg' + '_%sm_IBOUND.npy') % (topo_100_type, 100)
    #ibound_dict_dir_res = os.path.join(subreg_input_dir, '_%s_avg' + '_%sm_IBOUND.npy') % (topo_100_type, del_col_res_str)
    #ibound_dict_dir = os.path.join(subreg_input_dir, '_%s_avg' + '_025m_IBOUND.npy') % (topo_100_type)
    input_dict_dir = os.path.join(subreg_input_dir, '_input_lists.npy') 
    topo_nc = xr.open_dataset(topo_nc_dir)
    
    topo_gebco = topo_nc['gebco_avg'].values
    
    #   to build the IBOUND array and topo lists, use always GEBCO in the offshore part - this can be sometimes messed
    #   up in the coastal DEM case since it seems to have positive values even in the offshore part. So to tackle this
    #   the final ibound array will be a collage of two arrays - inland part taken from either merit or coastal dem 
    #   (if thats the scenario) and the offshore always GEBCO. In case the scenario specifies gebco inland and offshore
    #   then skip this step and just load all from the GEBCO ibound array.
    input_dict = np.load(input_dict_dir, allow_pickle = True).item()  
    #input_dict = np.load(input_dict_dir).item()     
    thk_lst = input_dict['lst_cs_thk']
    width_lst = input_dict['lst_cs_width'] 
    anchor_dist_lst = input_dict['lst_anchor_dist']
    anchor_depth_lst = input_dict['lst_anchor_depth']
    wtd_lst = fix_lst(input_dict['lst_wtd'])
    topo_lst = fix_lst(input_dict['lst_gebco_merit_avg_100m'])
    topo_gebco_500m_lst = fix_lst(input_dict['lst_topo_gebco_avg_500m'])
    soil_type_lst = fix_lst(input_dict['lst_soil_type'])
    soil_thk_lst = fix_lst(input_dict['lst_soil_thk'])
    nasa_lst = fix_lst(input_dict['lst_cs_nasa'])
    offshore_lst = fix_lst(input_dict['lst_cs_offshore'])
    pcr_lst = input_dict['lst_cs_pcr_rch']
    watergap_lst = input_dict['lst_cs_watergap_rch']
    p_min_et_lst = input_dict['lst_cs_p_min_et']
    k_soil_lst = fix_lst(input_dict['lst_cs_k_soil'])
    drn_rate_lst = fix_lst(input_dict['lst_cs_drn_rate'])
    glhymps_top_lay_lst = fix_lst(input_dict['lst_cs_glhymps_top_lay'])
    glhymps_bot_lay_lst = fix_lst(input_dict['lst_cs_glhymps_bot_lay'])
    
    if not os.path.exists(os.path.join(subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND_IBOUND.npy')):
        new_ibound_dict = rivbas.save_model_input_files_cst_type_topo_v2('', thk_lst, width_lst, anchor_dist_lst, anchor_depth_lst, wtd_lst, topo_lst, topo_gebco_500m_lst,\
                                                                      soil_type_lst, soil_thk_lst, nasa_lst, offshore_lst, pcr_lst, watergap_lst, p_min_et_lst, k_soil_lst, drn_rate_lst,\
                                                                      glhymps_top_lay_lst, glhymps_bot_lay_lst, 100., 100., 10., coscat_id, subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND')
    else:
        try:
            os.remove(os.path.join(subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND_IBOUND.npy'))
        except FileNotFoundError:
            pass
        new_ibound_dict = rivbas.save_model_input_files_cst_type_topo_v2('', thk_lst, width_lst, anchor_dist_lst, anchor_depth_lst, wtd_lst, topo_lst, topo_gebco_500m_lst,\
                                                                      soil_type_lst, soil_thk_lst, nasa_lst, offshore_lst, pcr_lst, watergap_lst, p_min_et_lst, k_soil_lst, drn_rate_lst,\
                                                                      glhymps_top_lay_lst, glhymps_bot_lay_lst, 100., 100., 10., coscat_id, subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND')    
        
