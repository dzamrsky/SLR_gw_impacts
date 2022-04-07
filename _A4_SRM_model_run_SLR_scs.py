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
#matplotlib.use('agg')
import matplotlib.pyplot as plt
#plt.ioff()
#from scipy.stats import norm
#   the id loop is the only input into the script, the rest will be found in the lookup table
#id_loop = int(sys.argv[1])
#param_combo = int(sys.argv[2])
id_loop = 1986
param_combo = 1

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
init_dir = r'f:\_cartesius_backup\_A3_backup\coscat_files'

"""
master_csv_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models/_coscat_reg_main_tb.csv'
model_sc_csv_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models/_model_scenarios_v2.csv'
input_dir = r'/projects/0/qt16165/_dzamrsky/_A4_MODEL_input_files'
cs_regs_shp_dir = r'/projects/0/qt16165/_dzamrsky/_A4_input_data/cs_reg_id.shp'
model_out_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_sens_analysis_v2'
init_dir = r'/projects/0/qt16165/_dzamrsky/A3_final_results/coscat_files'
swat_exe_dir = r'/home/dzamrsky/_1413_test/temp/swtv4'
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
    if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_delta', 'avg_sim_MERGED.nc')):
        init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_delta', 'avg_sim_MERGED.nc')
    else:
        if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_all', 'avg_sim_MERGED.nc')):
            init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_all', 'avg_sim_MERGED.nc')
        else:
            if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim_MERGED.nc')):
                init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim_MERGED.nc')
            else:
                if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_other', 'avg_sim_MERGED.nc')):
                    init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_other', 'avg_sim_MERGED.nc')
                else:
                    init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry', 'avg_sim_MERGED.nc')

elif cs_type == 2:
    if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_henry', 'avg_sim_MERGED.nc')):
        init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry', 'avg_sim_MERGED.nc')
    else:
        if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim_MERGED.nc')):
            init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim_MERGED.nc')
        else:
            init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_all', 'avg_sim_MERGED.nc')

elif cs_type == 3:
    if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_other', 'avg_sim_MERGED.nc')):
        init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_other', 'avg_sim_MERGED.nc')
    else:
        if os.path.exists(os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim_MERGED.nc')):
            init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_henry_other', 'avg_sim_MERGED.nc')
        else:
            init_dict_dir = os.path.join(init_dir, '_' + id_cs_str + '_all', 'avg_sim_MERGED.nc')

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
    os.remove(os.path.join(subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND_IBOUND.npy'))
    new_ibound_dict = rivbas.save_model_input_files_cst_type_topo_v2('', thk_lst, width_lst, anchor_dist_lst, anchor_depth_lst, wtd_lst, topo_lst, topo_gebco_500m_lst,\
                                                                  soil_type_lst, soil_thk_lst, nasa_lst, offshore_lst, pcr_lst, watergap_lst, p_min_et_lst, k_soil_lst, drn_rate_lst,\
                                                                  glhymps_top_lay_lst, glhymps_bot_lay_lst, 100., 100., 10., coscat_id, subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND')    
    
    
dict_in = np.load(os.path.join(subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND_IBOUND.npy'), allow_pickle = True).item()      
#dict_in = np.load(os.path.join(subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND_IBOUND.npy')).item()  
ibound_arr = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['ibound_arr']
top_elev = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['top_elev']
bot_elev = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['bot_elev']
lay_elev = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['lay_elev']
cst_offset = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['cst_offset'] 
lay_idx = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['lay_idx'] 
col_idx = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['col_idx'] 
end_idx = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['end_idx'] 
top = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['top'] 
zbot = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['zbot'] 
col_width = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['col_width']
idx_start = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['idx_start']
idx_end = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['idx_end']
x_start = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['x_start']
x_end = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['x_end']
soil_thk = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['soil_thk']
soil_type = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['soil_type'] 
cst_idx = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['cst_idx']
wtd_elev = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['wtd_elev']
k_soil = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['k_soil']
drn_dens = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['drn_rate'] 
glhymps_top_lay = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['glhymps_top_lay']
glhymps_bot_lay = dict_in[('_%s_' + '%s_IBOUND') % (topo_100_type, int(del_col_res))]['glhymps_bot_lay']        
        
#   make sure to change the elevation DEM type, but if we do it like below it also maintains the same extent of the model domain
#if topo_100_type != 'gebco':
#    topo_res = topo_nc['gebco_%s_avg' % (topo_100_type)].values.tolist()[int((200 + x_start) * (1000. / del_col_res)) : int((200 + x_end) * (1000. / del_col_res))]  
#else:
#    topo_res = topo_nc['gebco_avg'].values.tolist()[int((200 + x_start) * (1000. / del_col_res)) : int((200 + x_end) * (1000. / del_col_res))]          
        
#   get other input
rch_csv = os.path.join(os.path.join(input_dir, id_cs_str), '_SUBREG_representative_models_RCH_stats.csv')
glhymps_csv = os.path.join(os.path.join(input_dir, id_cs_str), '_SUBREG_representative_models_GEOLOGY_stats.csv')
geo_csv = os.path.join(input_dir, 'coscat_geology_input.csv') 

#   get the lognormal mu and stdev values from the csv file for RCH and GEOLOGY
#glhymps_bot_param_ln = mfc.read_param_from_csv_subreg(coscat_id, subreg_id, glhymps_csv, 'glhymps_bot')
#glhymps_top_param_ln = mfc.read_param_from_csv_subreg(coscat_id, subreg_id, glhymps_csv, 'glhymps_top')

#   define the model name and all directories
scenario_name = '_SRM_'
#model_name_basis = '%s_%s_%s' % (id_cs_str + '_SRM_' + str_id_subreg, topo_type, str(int(del_col_res)) + 'm_' + str(int(del_lay_res)) + 'm')
model_name_basis = '%s_SRM_%s' % (id_cs_str, str_id_subreg) 

#   create folders for the specific model run
out_modelrun_dir = os.path.join(model_out_dir, model_name_basis)
out_dir = model_out_dir
summary_folder = os.path.join(model_out_dir, r'_summary_%s' % (model_name_basis)) 
geo_summary_csv_dir = os.path.join(summary_folder, '_geo_summary_%s.csv' % (model_name_basis))
res_coscat_summary_csv_dir = os.path.join(summary_folder, '_res_coscat_summary_%s.csv' % (model_name_basis))

#   define stress periods, sp_00 last 10 000 years and serves as stabilizing stress period with the absolute minimum
#   LGM sea level. In that period we take the COSCATs concentration distribution at LGM, interpolate to the shape of the 
#   SRMs model domain and run for 10 000 years to get an initial pre-SLR condition. 
equilibrium_names = ['BP_30000_to_20000', 'BP_20000_to_19000', 'BP_19000_to_18000', 'BP_18000_to_17000', 'BP_17000_to_16000',\
                     'BP_16000_to_15000', 'BP_15000_to_14000', 'BP_14000_to_13000', 'BP_13000_to_12000', 'BP_12000_to_11000',\
                     'BP_11000_to_10000', 'BP_10000_to_09000', 'BP_09000_to_08000', 'BP_08000_to_07000', 'BP_07000_to_06000',\
                     'BP_06000_to_05000', 'BP_05000_to_04000', 'BP_04000_to_03000', 'BP_03000_to_02000', 'BP_02000_to_01900',\
                     'BP_01900_to_01800', 'BP_01800_to_01700', 'BP_01700_to_01600', 'BP_01600_to_01500', 'BP_01500_to_01400',\
                     'BP_01400_to_01300', 'BP_01300_to_01200', 'BP_01200_to_01100', 'BP_01100_to_01000', 'BP_01000_to_00900',\
                     'BP_00900_to_00800', 'BP_00800_to_00700', 'BP_00700_to_00600', 'BP_00600_to_00500', 'BP_00500_to_00400',\
                     'BP_00400_to_00300', 'BP_00300_to_00200', 'BP_00200_to_00100', 'BP_00100_to_00000']

sea_levels = [-130., -129., -125., -126., -124., -110., -95., -75., -67., -57,\
              -45, -32., -17., -7., -3, -2., -1.5, -0.7, -0.2, 0.,\
              0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0.,\
              0., 0., 0., 0., 0., 0., 0., 0.]

sp_duration = [10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,\
               1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 100,\
               100, 100, 100, 100, 100, 100, 100, 100, 100, 100,\
               100, 100, 100, 100, 100, 100, 100, 100, 100]

yrs_len_lst = [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,\
               1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 100,\
               100, 100, 100, 100, 100, 100, 100, 100, 100, 100,\
               100, 100, 100, 100, 100, 100, 100, 100, 100]

perlen_lst = [365.25 * i for i in yrs_len_lst]                     # an array of the stress period lengths
nstp_lst = [int(i / 50) for i in yrs_len_lst]                                 # number of time steps in each stress period
nper = 1#len(perlen)


#   loop through the SRM folder and get lognormal distribution of the GW recharge values, per stress period
#gw_rch_files = [f for f in os.listdir(subreg_input_dir) if f.endswith('_gw_rch.nc')]
gw_rch_files = [f for f in os.listdir(subreg_input_dir) if f.endswith('_gw_rch_clay.nc')]
topo_files = [f for f in os.listdir(subreg_input_dir) if f.endswith('pt_dist_100_topo.nc')]

#   create the gw_rch input lists, take from 9th element onwards because first stress period is 10k years and the groundwater
#   recharge values were extracted for each 1000 yers over that time period.
gw_rch_ln_vals = rivbas.get_gw_rch_lognorm(subreg_input_dir, gw_rch_files, topo_files, x_start, x_end)    
gw_rch_mu = gw_rch_ln_vals[0][9:]
gw_rch_std = gw_rch_ln_vals[1][9:]  
    

"""   -----------------------------------------------------------------------------------   """
"""                 Get the geology parameters and combinations of those                    """ 
"""   -----------------------------------------------------------------------------------   """

#   get the geology parameters
geo_params = mfc.read_geo_params(coscat_id, geo_csv)
#   derive the input for the synthetic geology profile generation
mud_pct, sand_pct, y_val = geo_params[-3], geo_params[-4], geo_params[-5]
qs_cat, sm_cat, y_cat = geo_params[0], geo_params[1], geo_params[2]
mud_shlf_pct, mud_slp_pct = geo_params[-1], geo_params[-2]
#   get the sed flux category
if qs_cat == 1:
    sed_flux = 'low'
elif qs_cat == 2:
    sed_flux = 'medium'
else:
    sed_flux = 'high'
#   get the sed type category (always going to be 'small' and the clay will be filled in based on the mud %)
sed_type = 'small'

#   specify the parameter values for the parameters that will be varied
#sand_lay_n = 3#[2, 5]             #   number of glhymps1 + glhymps 2 layers (one low/high stand combinaiton)
inland_aqf_scenario = 3         #   always will be 3 
inland_aqf_lrs = []
variation_pct = 25.             #   how much will the total % of the thickness vary for a sand layer 
tot_lay_thk = 0.                #   also has to be set to zero first
#p_fact = 0.75
off_lay_start = 0.
#clay_cap_thk = 20.

#   if the mud % is > 50%:
if mud_shlf_pct > 50.:
    clay_cap_shelf = True       #   insert the clay capping layer on top of the continental shelf/slope, the same reworking
else:
    clay_cap_shelf = False

if mud_slp_pct > 50.:    
    clay_cap_slope = True       #   parameter lay_pres_y1 will apply to these layers as well
else:
    clay_cap_slope = False       #   parameter lay_pres_y1 will apply to these layers as well
    
clay_cap_shelf_thk = clay_cap_thk                       
clay_cap_slope_thk = clay_cap_thk
lay_pres_y1 = round(y_val, 1)          #   parameter lay_pres_y1 will apply to these layers as well
const_geo_hk_vals = False

#   define the initial conditions of the model and the starting concentration and heads arrays
init_cond = 'salt'
sconc_arr = None
strt_arr = None

#   random seed is always going to be 1 in the case of the parameter testing runs
rand_seed = 1
random.seed(rand_seed)  

#   in case it is the first random seed create the overall csv file and also the folder to store the geology pictures
#if param_combo == 1:

#   create the summary folder and subfolders
if not os.path.exists(summary_folder):
    os.makedirs(summary_folder, exist_ok = True)
    os.makedirs(os.path.join(summary_folder, '_geology'), exist_ok = True)      
#   create the overall summary csv files
if not os.path.exists(geo_summary_csv_dir):
    #   the goeological summary
    f_geo = open(geo_summary_csv_dir,'w')
    geo_headers = 'id_seed, sed_flux, low_high_lay_N, mud_pct, pres_pct, mud_cap_shelf, mud_cap_slope, mud_cap_shelf_thk_m, mud_cap_slope_thk_m, p_fact, mud_lays_strt,\
                   high_lay_1, low_lay_1, high_lay_2, low_lay_2, high_lay_3, low_lay_3, high_lay_4, low_lay_4, high_lay_5, low_lay_5'
    f_geo.write(" ".join(geo_headers.split())) 
    f_geo.close()      
if not os.path.exists(res_coscat_summary_csv_dir):                 
    #   the results summary
    f_res_coscat = open(res_coscat_summary_csv_dir,'w')
    res_coscat_headers = 'id_seed, eq_0_t_yrs, eq_0_t_runtime_secs, eq_0_frsh_inl_vol, eq_0_frsh_inl_pct, eq_0_frsh_shelf_vol, eq_0_frsh_shelf_pct,\
                            eq_10_t_yrs, eq_10_t_runtime_secs, eq_10_frsh_inl_vol, eq_10_frsh_inl_pct, eq_10_frsh_shelf_vol, eq_10_frsh_shelf_pct,\
                            eq_10_DSP_t_yrs, eq_10_DSP_t_runtime_secs, eq_10_DSP_frsh_inl_vol, eq_10_DSP_frsh_inl_pct, eq_10_DSP_frsh_shelf_vol, eq_10_DSP_frsh_shelf_pct'
    f_res_coscat.write(" ".join(res_coscat_headers.split())) 
    f_res_coscat.close()           

#       loop through the number of sand layers and calculate a thickness of each sand layer
#       define the total fractions of the model domain based on mud percentage
pct_glh_1 = int(round(mud_pct))
pct_glh_2 = 100 - pct_glh_1

#   split the percentages into the number of layers
def chunk_lst(tot_sum, n_lays):
    dividers = sorted(random.sample(range(1, tot_sum), n_lays - 1))
    return [a - b for a, b in zip(dividers + [tot_sum], [0] + dividers)]

#   get the individual layer fractions
glh_1_layers = chunk_lst(pct_glh_1, sand_lay_n)
glh_2_layers = chunk_lst(pct_glh_2, sand_lay_n)

#   reset the list
inland_aqf_lrs = []
#   combine those into the final list
for a in range(sand_lay_n):
    inland_aqf_lrs.append([int(glh_1_layers[a]), 'glhymps_1'])
    inland_aqf_lrs.append([int(glh_2_layers[a]), 'glhymps_2'])      

#   sometimes a layer can have a negative thickness, if that is the case then remove 1% of each layer above and check if the 
#   previously negative thickness layer is now positive, if not then repeat. Also include the case when the layer has 0 thickness
for b in range(len(inland_aqf_lrs)):
    if inland_aqf_lrs[b][0] <= 0:
        while inland_aqf_lrs[b][0] <= 0:
            for c in range(b):
                inland_aqf_lrs[c][0] -= 1
                inland_aqf_lrs[b][0] += 1

#   get the pct of the layers into the table and string
lay_thk_pct = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
for h in range(len(inland_aqf_lrs)):
    lay_thk_pct[h] = inland_aqf_lrs[h][0]
lay_thk_str = ''
for g in range(len(lay_thk_pct)):
    lay_thk_str += str(lay_thk_pct[g]) + ','

def linreg(X, Y):
    """
    return a,b in solution to y = ax + b such that root mean square distance between trend line and original points is minimized
    """
    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in zip(X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y
    det = Sxx * N - Sx * Sx
    return (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det

""" ----------------------------------------------------------------------------------------------- """
""" CREATE THE INPUT BEFORE RUNNING THE MODEL - SO IT DOESNT NEED TO BE MODIFIED IN THE LOOP ITSELF """
""" ----------------------------------------------------------------------------------------------- """
#   create names
model_name_basis = '%s_SRM_%s' % (id_cs_str, str_id_subreg) 
if param_combo >= 10:
    param_combo_str = str(param_combo)
else:
    param_combo_str = '0' + str(param_combo)
sc_dir = os.path.join(out_dir, model_name_basis, 'sc_' + param_combo_str)
os.makedirs(sc_dir, exist_ok = True)

if not os.path.exists(os.path.join(sc_dir, '_INIT_dict.npy')):    
    #   create the model object and the top and bottom elevation lists 
    model = cs_model(None, top_elev, None, None, None, None, None, None, None, None, None, None, None, None, None, None,\
                     None, None, None, False)
    
    #   trum the top of the IBOUND array if there are only inactive layers
    trm_idx = 0
    for i in range(ibound_arr.shape[0]):
        if len([g for g in ibound_arr[i, 0, :] if g == 1]) == 0:
            trm_idx += 1
        else:
            break
    ibound_arr = ibound_arr[trm_idx:, :, :]
    
    model.ibound_arr = ibound_arr
    model.top_elev = top_elev
    model.bot_elev = bot_elev
    model.wtd_elev = wtd_elev
    model.lay_elev = lay_elev[trm_idx:]
    model.cst_offset = cst_offset
    model.lay_idx = lay_idx
    model.col_idx = col_idx
    model.end_idx = end_idx
    model.top = math.ceil(np.nanmax(top_elev) / del_lay) * del_lay  # top
    model.zbot = zbot
    model.col_width = col_width
    model.idx_start = idx_start
    model.idx_end = idx_end
    model.x_start = x_start
    model.x_end = x_end
    model.soil_thk = soil_thk
    model.soil_type = soil_type
    model.cst_idx = cst_idx
    model.delv = del_lay
    model.k_soil = k_soil
    model.drn_dens = drn_dens
    model.glhymps_top_lay = glhymps_top_lay
    model.glhymps_bot_lay = glhymps_bot_lay
    model.coscat_id = coscat_id
    model.del_col_res = del_col_res
    model.del_lay_res = del_lay_res
    model.del_col = 100.
    model.del_lay = 10.
    
    try:
        if len(model.soil_thk) < len(model.soil_type):
            model.soil_thk = [model.soil_thk] * len(model.soil_type)
    except TypeError:
        model.soil_thk = [model.soil_thk] * len(model.soil_type)
    
    #   check if there is ocean/sea on the left side of the model domain
    new_model = mfc.trim_model_domain(model, model.top_elev, model.ibound_arr)
    model.x_start = new_model[0]
    model.ibound_arr = new_model[1]
    model.ncol = new_model[2]
    trim_idx_st = new_model[-3] 
    trim_idx = new_model[-2]

    #   if yes then assign or cut the model properties so they match the new dimensions
    if new_model[-1] != 0 and trim_idx != -1:
        model.x_start = new_model[0]
        model.ibound_arr = new_model[1]
        model.ncol = new_model[2]
        model.top_elev = model.top_elev[trim_idx_st: trim_idx]
        model.bot_elev = model.bot_elev[trim_idx_st: trim_idx]
        model.wtd_elev = model.wtd_elev[trim_idx_st: trim_idx]
        model.soil_thk = model.soil_thk[trim_idx_st: trim_idx]
        model.soil_type = model.soil_type[trim_idx_st: trim_idx]
        model.k_soil = model.k_soil[trim_idx_st: trim_idx] 
        model.drn_dens = model.drn_dens[trim_idx_st: trim_idx] 
        model.glhymps_top_lay = model.glhymps_top_lay[trim_idx_st: trim_idx] 
        model.glhymps_bot_lay = model.glhymps_bot_lay[trim_idx_st: trim_idx]                 
    elif new_model[-1] == 0 and trim_idx_st == 0 and trim_idx != -1:
        model.ibound_arr = new_model[1]
        model.ncol = new_model[2]                
        model.top_elev = model.top_elev[new_model[-1]:trim_idx]
        model.bot_elev = model.bot_elev[new_model[-1]:trim_idx]
        model.wtd_elev = model.wtd_elev[new_model[-1]:trim_idx]
        model.soil_thk = model.soil_thk[new_model[-1]:trim_idx]
        model.soil_type = model.soil_type[new_model[-1]:trim_idx]
        model.k_soil = model.k_soil[new_model[-1]:trim_idx] 
        model.drn_dens = model.drn_dens[new_model[-1]:trim_idx] 
        model.glhymps_top_lay = model.glhymps_top_lay[new_model[-1]:trim_idx] 
        model.glhymps_bot_lay = model.glhymps_bot_lay[new_model[-1]:trim_idx]  
    elif new_model[-1] == 0 and trim_idx_st != 0 and trim_idx != -1:
        model.ibound_arr = new_model[1]
        model.ncol = new_model[2]                
        model.top_elev = model.top_elev[trim_idx_st: trim_idx]
        model.bot_elev = model.bot_elev[trim_idx_st: trim_idx]
        model.wtd_elev = model.wtd_elev[trim_idx_st: trim_idx]
        model.soil_thk = model.soil_thk[trim_idx_st: trim_idx]
        model.soil_type = model.soil_type[trim_idx_st: trim_idx]
        model.k_soil = model.k_soil[trim_idx_st: trim_idx] 
        model.drn_dens = model.drn_dens[trim_idx_st: trim_idx] 
        model.glhymps_top_lay = model.glhymps_top_lay[trim_idx_st: trim_idx] 
        model.glhymps_bot_lay = model.glhymps_bot_lay[trim_idx_st: trim_idx]                
    elif new_model[-1] != 0 and trim_idx == -1:
        model.top_elev = model.top_elev[new_model[-3]:]
        model.bot_elev = model.bot_elev[new_model[-3]:]
        model.wtd_elev = model.wtd_elev[new_model[-3]:]
        model.soil_thk = model.soil_thk[new_model[-3]:]
        model.soil_type = model.soil_type[new_model[-3]:]
        model.k_soil = model.k_soil[new_model[-3]:] 
        model.drn_dens = model.drn_dens[new_model[-3]:] 
        model.glhymps_top_lay = model.glhymps_top_lay[new_model[-3]:] 
        model.glhymps_bot_lay = model.glhymps_bot_lay[new_model[-3]:]                 
    elif new_model[-1] == 0 and trim_idx == -1:
        pass

    """
    #   if yes then assign or cut the model properties so they match the new dimensions
    if new_model[-1] != 0 and trim_idx != -1:
        model.x_start = new_model[0]
        model.ibound_arr = new_model[1]
        model.ncol = new_model[2]
        model.top_elev = model.top_elev[new_model[-3]:trim_idx]
        model.bot_elev = model.bot_elev[new_model[-3]:trim_idx]
        model.wtd_elev = model.wtd_elev[new_model[-3]:trim_idx]
        model.soil_thk = model.soil_thk[new_model[-3]:trim_idx]
        model.soil_type = model.soil_type[new_model[-3]:trim_idx]
        model.k_soil = model.k_soil[new_model[-3]:trim_idx] 
        model.drn_dens = model.drn_dens[new_model[-3]:trim_idx] 
        model.glhymps_top_lay = model.glhymps_top_lay[new_model[-3]:trim_idx] 
        model.glhymps_bot_lay = model.glhymps_bot_lay[new_model[-3]:trim_idx]                 
    elif new_model[-1] == 0 and trim_idx != -1:
        model.ibound_arr = new_model[1]
        model.ncol = new_model[2]                
        model.top_elev = model.top_elev[new_model[-1]:trim_idx]
        model.bot_elev = model.bot_elev[new_model[-1]:trim_idx]
        model.wtd_elev = model.wtd_elev[new_model[-1]:trim_idx]
        model.soil_thk = model.soil_thk[new_model[-1]:trim_idx]
        model.soil_type = model.soil_type[new_model[-1]:trim_idx]
        model.k_soil = model.k_soil[new_model[-1]:trim_idx] 
        model.drn_dens = model.drn_dens[new_model[-1]:trim_idx] 
        model.glhymps_top_lay = model.glhymps_top_lay[new_model[-1]:trim_idx] 
        model.glhymps_bot_lay = model.glhymps_bot_lay[new_model[-1]:trim_idx]                 
    elif new_model[-1] != 0 and trim_idx == -1:
        model.top_elev = model.top_elev[new_model[-3]:]
        model.bot_elev = model.bot_elev[new_model[-3]:]
        model.wtd_elev = model.wtd_elev[new_model[-3]:]
        model.soil_thk = model.soil_thk[new_model[-3]:]
        model.soil_type = model.soil_type[new_model[-3]:]
        model.k_soil = model.k_soil[new_model[-3]:] 
        model.drn_dens = model.drn_dens[new_model[-3]:] 
        model.glhymps_top_lay = model.glhymps_top_lay[new_model[-3]:] 
        model.glhymps_bot_lay = model.glhymps_bot_lay[new_model[-3]:]                 
    elif new_model[-1] == 0 and trim_idx == -1:
        pass
    """
    
    model.top = math.ceil(np.nanmax(model.top_elev) / del_lay) * del_lay  # top
    model.zbot = model.top - model.ibound_arr.shape[0] * model.del_lay
    
    #   find the shelf edge and indexes of cont shelf and coastline
    if init_cond != 'fresh':
        cont_shelf_edge = mgt.find_shelf_break(model, model.top_elev, True) #   try to find the shelf break     
    else:
        cont_shelf_edge = mgt_frsh.find_shelf_break_fresh(model, model.top_elev, True) #   try to find the shelf break     
    model.x_start = round(model.x_start + model.cst_offset, 2)
    model.x_end = model.ibound_arr.shape[-1] / 10. + model.x_start
    idx_cst_0 = int(abs(round(model.x_start + model.cst_offset, 2)) * 10)
    idx_shelf_edge_0 = int((cont_shelf_edge[0] - model.x_start + model.cst_offset) * 10)
    if idx_shelf_edge_0 > model.ibound_arr.shape[-1]:
        idx_shelf_edge_0 = model.ibound_arr.shape[-1] - 1
    
    #   max number of layers allowed is 999...
    if model.ibound_arr.shape[0] > 999:
        model.nlay = 999
        model.ibound_arr = model.ibound_arr[:999, :, :]
 
    #   create the HK array, only for "saline" initial condition, also make sure there are no nan values
    figname = '_COSCAT_%s' % (id_cs_str + '_SRM_' + str_id_subreg) + '_GeoRandSc_%s' % str(param_combo)  
    model.botm = model.lay_elev[1:] #   assign this to create hk_arr
    hk_arr = mgt.create_geology_profile(model, rand_seed, inland_aqf_lrs, p_fact, mud_pct, sed_type, lay_pres_y1, clay_cap_shelf_thk, clay_cap_slope_thk,\
                                   off_lay_start, sed_flux, sc_dir, os.path.join(summary_folder, '_geology'), const_geo_hk_vals, clay_cap_shelf, clay_cap_slope, figname)  
    model.hk_arr = hk_arr
    model.vk_arr = hk_arr * 0.1
    model.laytyp = laytyp_val
    
    hk_top_mean = np.nanmean(model.hk_vals_top, axis=0)
    hk_top_std = np.nanstd(model.hk_vals_top, axis=0)
    hk_bot_mean = np.nanmean(model.hk_vals_bot, axis=0)
    hk_bot_std = np.nanstd(model.hk_vals_bot, axis=0)    
    
    glh_1_val = round(hk_bot_mean, 4)
    glh_2_val = round(hk_top_mean, 4)
    
    for v in range(model.ibound_arr.shape[0]):
        for w in range(model.ibound_arr.shape[-1]):
            if model.ibound_arr[v, 0, w] == 1 and math.isnan(model.hk_arr[v, 0, w]):
                print(v, w)
    
    #   check last time that dimensions are ok, ibound_arr has priority
    if model.ncol != model.ibound_arr.shape[-1]:
        model.ncol = model.ibound_arr.shape[-1]
        model.delr = model.delr[:model.ibound_arr.shape[-1]]        
        
    #   calculate the x_st and x_end taking into account the coastal offset
    x_st = round(model.x_start + model.cst_offset, 2)
    x_end = model.ibound_arr.shape[-1] / 10. + x_st
    #y_top = math.ceil(max(model.top_elev) / 10.) * 10.
    y_top = model.top
    y_bot = model.zbot
    #   create lists with 
    x_cells = np.linspace(x_st + 0.05, x_end - 0.05, model.ibound_arr.shape[-1])
    x_cells = [round(k, 2) for k in x_cells]
    y_cells = np.linspace(y_top - 5., y_bot - 5., model.ibound_arr.shape[0] + 1)
    y_cells = [round(k, 2) for k in y_cells]      
    
    sc_arr = np.ones((model.ibound_arr.shape[0], model.ibound_arr.shape[1], model.ibound_arr.shape[-1])) * -1
    st_arr = np.ones((model.ibound_arr.shape[0], model.ibound_arr.shape[1], model.ibound_arr.shape[-1])) * -1
    #   for each column in the model domain, try to find the right column from the averaged profile

    #   open the COSCAT results and create initial head and concentration plots
    if os.path.isfile(init_dict_dir):
        last_nc = xr.open_dataset(init_dict_dir)
        last_nc.sel()
        conc_vals_all = last_nc['solute concentration']
        head_vals_all = last_nc['heads']

        #   select the concentration and heads at the last LGM - which corresponds to eq_21               
        eq_21_ts = last_nc['time'].values.tolist()[-1] #  set it to 0 first, in case the model type was fresh in A3, it will take the first concentration profile as initial condition
        for f in range(len(last_nc['eq'].values.tolist())):
            if last_nc['eq'].values.tolist()[f] == 'eq_21':
                eq_21_ts = last_nc['time'].values.tolist()[f]
                #print(last_nc['time'].values.tolist()[f])
        last_conc = last_nc.sel(time = eq_21_ts)['solute concentration']
        last_head = last_nc.sel(time = eq_21_ts)['heads']
  
        for a in range(model.ibound_arr.shape[0]):
            ibound_act_lay_idxs = [k for k, x in enumerate(model.ibound_arr[a, 0, :].tolist()) if x == 1]    
            if len(ibound_act_lay_idxs) > 0:
                try:
                    y_coord = y_cells[a]
                    #   get active cells in the column
                    ibound_act_lay_idxs = [k for k, x in enumerate(model.ibound_arr[a, 0, :].tolist()) if x == 1]    
                    y_coord_coscat = min(last_conc['y'].values.tolist(), key = lambda x : abs(x - y_coord))
                    ibound_act_lay_idxs_coscat = [k for k, x in enumerate(last_conc.sel(y = y_coord_coscat).values.tolist()) if x < 100.] 
                    #   select the y coordinates of the concentration array for the SRM and COSCAT arrays
                    conc_val_coscat = last_conc.sel(y = y_coord_coscat).values[ibound_act_lay_idxs_coscat[0] : ibound_act_lay_idxs_coscat[0] + len(ibound_act_lay_idxs_coscat)]
                    head_val_coscat = last_head.sel(y = y_coord_coscat).values[ibound_act_lay_idxs_coscat[0] : ibound_act_lay_idxs_coscat[0] + len(ibound_act_lay_idxs_coscat)]
                    #   resize the concentration into the size of the SRM column 
                    size = len(ibound_act_lay_idxs_coscat)
                    xloc = np.arange(size)
                    newsize = len(ibound_act_lay_idxs)
                    new_xloc = np.linspace(0, size, newsize)
                    srm_conc_vals = np.interp(new_xloc, xloc, conc_val_coscat)                                
                    srm_head_vals = np.interp(new_xloc, xloc, head_val_coscat)                                     
                    for b in range(len(ibound_act_lay_idxs)):
                        sc_arr[a, 0, ibound_act_lay_idxs[b]] = srm_conc_vals[b]
                        st_arr[a, 0, ibound_act_lay_idxs[b]] = srm_head_vals[b]                                    
                        
                    #   old version
                    #conc_val = last_conc.sel(x = x_coord).values[ibound_act_lay_idxs[b]]
                    #head_val = last_head.sel(x = x_coord).values[ibound_act_lay_idxs[b]]       
                except (IndexError, KeyError): # if the average COSCAT array is smaller than the actual HYBAS IBOUND
                    #   if the x_coord is further from coastline than 20km then assign fresh water concentration and head value of the topography
                    for b in range(len(ibound_act_lay_idxs)):
                        if y_coord < last_conc.y.values[0]:
                            conc_val = last_conc.sel(y = last_conc.y.values[0]).values[ibound_act_lay_idxs[b]]
                            head_val = last_head.sel(y = last_head.y.values[0]).values[ibound_act_lay_idxs[b]]
                        else:
                            conc_val = 999.99
                            head_val = -999.99
                        sc_arr[a, 0, ibound_act_lay_idxs[b]] = conc_val
                        st_arr[a, 0, ibound_act_lay_idxs[b]] = head_val                            
                    
                    try:
                        if len(last_conc.sel(y = y_coord_coscat).values.tolist()) > 1:
                            conc_val_last = [k for k in last_conc.sel(y = y_coord_coscat).values.tolist() if k < 100.][-1]
                            head_val_last = [k for k in last_head.sel(y = y_coord_coscat).values.tolist() if k != -999.99][-1]                                        
                        else:
                            conc_val_last = [k for k in last_conc.sel(y = y_coord_coscat).values.tolist()[0] if k < 100.][-1]
                            head_val_last = [k for k in last_head.sel(y = y_coord_coscat).values.tolist()[0] if k != -999.99][-1]
                    #   in cases the new ibound array stretches further than the COSCAT ibound
                    except (IndexError, KeyError):
                        if y_coord < -0.:
                            if len(last_conc.sel(y = y_coord_coscat).values.tolist()) > 1:
                                conc_val_last = [k for k in last_conc.sel(y = last_conc.y.values).values.tolist() if k < 100.][-1]
                                head_val_last = [k for k in last_head.sel(y = last_head.y.values).values.tolist() if k != -999.99][-1]        
                            else:
                                conc_val_last = [k for k in last_conc.sel(y = last_conc.y.values[0]).values.tolist()[0] if k < 100.][-1]
                                head_val_last = [k for k in last_head.sel(y = last_head.y.values[0]).values.tolist()[0] if k != -999.99][-1]  
                        else:
                            if len(last_conc.sel(y = y_coord_coscat).values.tolist()) > 1:
                                conc_val_last = [k for k in last_conc.sel(y = last_conc.y.values[-1]).values.tolist() if k < 100.][-1]
                                head_val_last = [k for k in last_head.sel(y = last_head.y.values[-1]).values.tolist() if k != -999.99][-1]    
                            else:
                                conc_val_last = [k for k in last_conc.sel(y = last_conc.y.values[-1]).values.tolist()[0] if k < 100.][-1]
                                head_val_last = [k for k in last_head.sel(y = last_head.y.values[-1]).values.tolist()[0] if k != -999.99][-1]  
                            
                    sc_arr[a, 0, ibound_act_lay_idxs[b]] = conc_val
                    st_arr[a, 0, ibound_act_lay_idxs[b]] = head_val    
        
                #   fill in random values in place of nan values
                for j in range(len(ibound_act_lay_idxs)):
                    if sc_arr[a, 0, ibound_act_lay_idxs[j]] == 999.99:
                        try:
                            mean_conc_val = statistics.mean([k for k in sc_arr[a, 0, ibound_act_lay_idxs[0] : ibound_act_lay_idxs[-1]].tolist() if k > 0.0 and i < 100.0])
                            std_conc_val = statistics.stdev([k for k in sc_arr[a, 0, ibound_act_lay_idxs[0] : ibound_act_lay_idxs[-1]].tolist() if k > 0.0 and i < 100.0])
                        except StatisticsError:
                            mean_conc_val, std_conc_val = 0.0, 0.0
                        sc_arr[a, 0, ibound_act_lay_idxs[j]] = round(abs(np.random.normal(mean_conc_val, std_conc_val)), 2)
                        #sc_arr[a, 0, ibound_act_lay_idxs[j]] = conc_val_last
                 
                for j in range(len(ibound_act_lay_idxs)):
                    if st_arr[a, 0, ibound_act_lay_idxs[j]] <= -999.:
                        try:
                            mean_head_val = statistics.mean([k for k in st_arr[a, 0, ibound_act_lay_idxs[0] : ibound_act_lay_idxs[-1]].tolist() if k > -990.])
                            std_head_val = statistics.stdev([k for k in st_arr[a, 0, ibound_act_lay_idxs[0] : ibound_act_lay_idxs[-1]].tolist() if k > -990.])
                        except StatisticsError:
                            try:
                                mean_head_val, std_head_val = [k for k in st_arr[a, 0, ibound_act_lay_idxs[0] : ibound_act_lay_idxs[-1]].tolist() if k > -990.][0], 0.0
                            except IndexError:
                                mean_head_val, std_head_val = 0.0, 0.0
                        st_arr[a, 0, ibound_act_lay_idxs[j]] = round(np.random.normal(mean_head_val, std_head_val), 2)
                        #st_arr[a, 0, ibound_act_lay_idxs[j]] = head_val_last

    else:
        #   if there is no COSCAT file from A3, just create an initial concentration - saline offshore/fresh inland
        sc_arr[:, 0, model.cst_idx:] = 35.
        sc_arr[:, 0, :model.cst_idx] = 0.
        #   head array just make all heads to be 0.
        for u in range(len(model.top_elev)):
            st_arr[:, 0, u] = model.top_elev[u]
          
    #   define the 100m and resolution x coordinates based on which the interpolation will be done 
    #x_lst_100 = np.arange(model.x_start, round(model.x_end, 1), del_col / 1000.)
    #x_lst_res = np.arange(model.x_start, round(model.x_end, 1), del_col_res / 1000.)                              

    x_lst_100 = np.arange(model.x_start, model.x_start + (model.ncol * 0.1), del_col / 1000.)
    x_lst_res = np.arange(model.x_start, model.x_start + (model.ncol * 0.1), del_col_res / 1000.)     

    if len(model.top_elev) - len(x_lst_res) == 1:
        x_lst_100 = np.arange(model.x_start, model.x_start + (model.ncol * 0.1) - 0.1, del_col / 1000.)
        x_lst_res = np.arange(model.x_start, model.x_start + (model.ncol * 0.1) - 0.1, del_col_res / 1000.)     

    y_lst_10 = np.arange(model.top, model.top - model.ibound_arr.shape[0] * del_lay, - del_lay)
    y_lst_res = np.arange(model.top, y_lst_10[-1] - del_lay_res, - del_lay_res)  
    
    #   final check if there are nan values in either array
    if math.isnan(np.sum(sc_arr)):                    
        sc_arr = np.nan_to_num(sc_arr, copy=True, nan = 0.0, posinf = None, neginf = None)
    if math.isnan(np.sum(st_arr)):                    
        st_arr = np.nan_to_num(st_arr, copy=True, nan = 0.0, posinf = None, neginf = None)

    #   adapt the inland part based on the topo_res input, read in the lists with updated x_start and x_end positions
    if topo_100_type != 'gebco':
        topo_res = topo_nc['gebco_%s_avg' % (topo_100_type)].values.tolist()[max(0, int((200 + model.x_start + model.cst_offset) * (1000. / model.del_col_res)) + 1) : min(int(400000. / del_col_res), int((200 + model.x_end + model.cst_offset) * (1000. / model.del_col_res)) + 1)]  
    else:
        topo_res = topo_nc['gebco_avg'].values.tolist()[max(0, int((200 + model.x_start + model.cst_offset) * (1000. / model.del_col_res)) + 1) : min(int(400000. / del_col_res), int((200 + model.x_end + model.cst_offset) * (1000. / model.del_col_res)) + 1)]   

    new_ibound_arr = model.ibound_arr
    new_sconc_arr = sc_arr
    new_strt_arr = st_arr
    new_hk_arr = model.hk_arr
    new_bot_elev = np.array(model.bot_elev)

    #   adapt the top_elev and ibound array in the inland domain
    for c in range(int(abs(model.x_start * 1000 / model.del_col_res))):
        new_top = math.ceil(np.nanmax(top_elev) / del_lay_res) * del_lay_res
        #new_top = math.ceil(np.nanmax(topo_res) / del_lay_res) * del_lay_res
        top_trim = int((model.top - new_top) / del_lay_res)
        model.top = new_top
        y_lst_res_new = y_lst_res
        #   2) if the top_trim is negative it means we need to add rows on top of the IBOUND array
        if top_trim < 0:
            for t in range(abs(top_trim)):
                new_ibound_arr = np.vstack((np.array([model.ibound_arr.shape[-1] * [0]]), model.ibound_arr[:, 0, :]))
                new_ibound_arr = np.expand_dims(new_ibound_arr, axis = 1)                       
                new_sconc_arr = np.vstack((np.array([model.ibound_arr.shape[-1] * [0]]), new_sconc_arr[:, 0, :]))
                new_sconc_arr = np.expand_dims(new_sconc_arr, axis = 1)                        
                new_strt_arr = np.vstack((np.array([model.ibound_arr.shape[-1] * [0]]), new_strt_arr[:, 0, :]))
                new_strt_arr = np.expand_dims(new_strt_arr, axis = 1)                               
                new_hk_arr = np.vstack((np.array([model.ibound_arr.shape[-1] * [0]]), new_hk_arr[:, 0, :]))
                new_hk_arr = np.expand_dims(new_hk_arr, axis = 1)                            
                y_lst_res_new = np.hstack((y_lst_res_new[0] + del_lay_res, y_lst_res))
        elif top_trim > 0:
            new_ibound_arr = model.ibound_arr[top_trim:, :, :]
            new_sconc_arr = new_sconc_arr[top_trim:, :, :]
            new_strt_arr = new_strt_arr[top_trim:, :, :]
            new_hk_arr = new_hk_arr[top_trim:, :, :]
            y_lst_res_new = y_lst_res[top_trim:]
        else:
            new_ibound_arr = model.ibound_arr
            new_sconc_arr = sc_arr
            new_strt_arr = st_arr
            new_hk_arr = model.hk_arr
            new_bot_elev = np.array(model.bot_elev)
            #new_ibound_arr = model.ibound_arr
            #new_sconc_arr = new_sconc_arr
            #new_strt_arr = new_strt_arr
            #new_hk_arr = new_hk_arr  
    
    #   if the del_col and del_col_res and del_lay and del_lay_res are not the same change the resolution of the model inputs            
    if del_col != del_col_res:
        new_ibound_lst, new_sconc_lst, new_strt_lst, new_hk_lst = [], [], [], []
        #   create - row by row - the new array with new columns resolution
        for y in range(model.ibound_arr.shape[0]):
            if model.del_col_res < 100.:
                new_ibound_row = np.repeat(new_ibound_arr[y, 0, :].tolist(), 100. / model.del_col_res)
                new_sconc_row = np.repeat(new_sconc_arr[y, 0, :].tolist(),100. / model.del_col_res)
                new_strt_row = np.repeat(new_strt_arr[y, 0, :].tolist(),100. / model.del_col_res)
                new_hk_row = np.repeat(new_hk_arr[y, 0, :].tolist(),100. / model.del_col_res)                
            else:                
                new_ibound_row = np.interp(x_lst_res, x_lst_100, model.ibound_arr[y, 0, :]).tolist()
                new_sconc_row = np.interp(x_lst_res, x_lst_100, sc_arr[y, 0, :]).tolist()
                new_strt_row = np.interp(x_lst_res, x_lst_100, st_arr[y, 0, :]).tolist()
                new_hk_row = np.interp(x_lst_res, x_lst_100, model.hk_arr[y, 0, :]).tolist()
            new_sconc_lst.append(new_sconc_row)
            new_strt_lst.append(new_strt_row)
            new_ibound_lst.append(new_ibound_row)
            new_hk_lst.append(new_hk_row)
        new_top_elev = np.interp(x_lst_res, x_lst_100, model.top_elev)
        new_bot_elev = np.interp(x_lst_res, x_lst_100, model.bot_elev)
        new_ibound_arr = np.array(new_ibound_lst)
        new_ibound_arr = np.expand_dims(new_ibound_arr, axis = 1)
        #new_ibound_arr = new_ibound_arr.astype(int)              
        new_hk_arr = np.array(new_hk_lst)
        new_hk_arr = np.expand_dims(new_hk_arr, axis = 1)
        new_sconc_arr = np.array(new_sconc_lst)
        new_sconc_arr = np.expand_dims(new_sconc_arr, axis = 1)
        new_strt_arr = np.array(new_strt_lst)
        new_strt_arr = np.expand_dims(new_strt_arr, axis = 1)  
        model.top_elev = new_top_elev.tolist()
        model.bot_elev = new_bot_elev.tolist()
        model.ibound_arr = new_ibound_arr
        
    if del_lay != del_lay_res:
        new_ibound_lst, new_sconc_lst, new_strt_lst, new_hk_lst = [], [], [], []
        if model.del_lay_res < 10:
            y_lst_res_new = np.linspace(model.top, model.top - new_ibound_arr.shape[0] *  (10. / model.del_lay_res) * model.del_lay_res, int(new_ibound_arr.shape[0] *  (10. / model.del_lay_res)), endpoint = False)
        else:
            y_lst_res_new = y_lst_res
        if del_col == del_col_res:
            new_ibound_arr = model.ibound_arr
            new_sconc_arr = sc_arr
            new_strt_arr = st_arr
            new_hk_arr = model.hk_arr
            new_bot_elev = np.array(model.bot_elev)
        for j in range(new_ibound_arr.shape[-1]):
            #   need to invert all the lists because of np.interp gives bullshit results if xp and x are not increasing.. 
            if model.del_lay_res < 10:
                new_ibound_col = np.repeat(new_ibound_arr[:, 0, j].tolist(),10. / model.del_lay_res)
                new_sconc_col = np.repeat(new_sconc_arr[:, 0, j].tolist(),10. / model.del_lay_res)
                new_strt_col = np.repeat(new_strt_arr[:, 0, j].tolist(),10. / model.del_lay_res)
                new_hk_col = np.repeat(new_hk_arr[:, 0, j].tolist(),10. / model.del_lay_res)
            else:
                new_ibound_col = np.interp(y_lst_res[::-1], y_lst_10[::-1], new_ibound_arr[:, 0, j].tolist()[::-1]).tolist()[::-1]
                new_sconc_col = np.interp(y_lst_res[::-1], y_lst_10[::-1], new_sconc_arr[:, 0, j].tolist()[::-1]).tolist()[::-1]
                new_strt_col = np.interp(y_lst_res[::-1], y_lst_10[::-1], new_strt_arr[:, 0, j].tolist()[::-1]).tolist()[::-1]
                new_hk_col = np.interp(y_lst_res[::-1], y_lst_10[::-1], new_hk_arr[:, 0, j].tolist()[::-1]).tolist()[::-1]
            new_sconc_lst.append(new_sconc_col)
            new_strt_lst.append(new_strt_col)      
            new_ibound_lst.append(new_ibound_col)             
            new_hk_lst.append(new_hk_col)
        new_ibound_arr = np.array(new_ibound_lst)
        new_ibound_arr = np.swapaxes(new_ibound_arr, 0, 1)
        new_ibound_arr = np.expand_dims(new_ibound_arr, axis = 1)
        #new_ibound_arr = new_ibound_arr.astype(int)                 
        new_hk_arr = np.array(new_hk_lst)
        new_hk_arr = np.swapaxes(new_hk_arr, 0, 1)
        new_hk_arr = np.expand_dims(new_hk_arr, axis = 1)      
        new_sconc_arr = np.array(new_sconc_lst)
        new_sconc_arr = np.swapaxes(new_sconc_arr, 0, 1)
        new_sconc_arr = np.expand_dims(new_sconc_arr, axis = 1)
        new_strt_arr = np.array(new_strt_lst)
        new_strt_arr = np.swapaxes(new_strt_arr, 0, 1)
        new_strt_arr = np.expand_dims(new_strt_arr, axis = 1)       
        model.ibound_arr = new_ibound_arr
    
    #   transform the IBOUND - inland part, so it fits the elevation of the resolution topography
    if del_col != del_col_res or del_lay != del_lay_res:
        #   1) check the inland indexes and only adapt the topography from the inland boundary until the coastline
        idx_cst_res = int(abs(x_lst_100[0] * 1000.) / del_col_res)  #   coast index in the res topo list

        #   for the layers above 0m current sea level, go column by column and adapt the active cells
        for i in range(model.ibound_arr.shape[-1]):
            topo_res_col = model.top_elev[i]
            #for j in range(int(y_lst_res_new[0] / del_lay_res)):
            for j in range(model.ibound_arr.shape[0]):
                mid_cell_elev = y_lst_res_new[j] - (del_lay_res / 2.)
                if topo_res_col > mid_cell_elev and mid_cell_elev > model.bot_elev[i]:
                    #print(i, j, topo_res[j], mid_cell_elev)
                    new_ibound_arr[j, 0, i] = 1

                    closest_vals = new_sconc_arr[max(j - 2, 0) : min(j + 2, model.ibound_arr.shape[0]), 0, max(i - 2, 0) : min(i + 2, model.ibound_arr.shape[-1])]
                    closest_vals[closest_vals < 0] = np.nan
                    closest_val = np.nanmean(closest_vals)
                    new_sconc_arr[j, 0, i] = closest_val

                    closest_vals = new_strt_arr[max(j - 2, 0) : min(j + 2, model.ibound_arr.shape[0]), 0, max(i - 2, 0) : min(i + 2, model.ibound_arr.shape[-1])]
                    closest_vals[np.abs(closest_vals) <= -999.] = np.nan
                    closest_val = np.nanmean(closest_vals)
                    new_strt_arr[j, 0, i] = closest_val

                    #closest_vals = new_hk_arr[max(j - 2, 0) : min(j + 2, model.ibound_arr.shape[0]), 0, max(i - 2, 0) : min(i + 2, model.ibound_arr.shape[-1])]
                    #closest_val = np.nanmean(closest_vals)
                    #new_hk_arr[j, 0, i] = closest_val
                    
                else:
                    new_ibound_arr[j, 0, i] = 0
                    new_sconc_arr[j, 0, i] = 0.0
                    new_strt_arr[j, 0, i] = -1
                    #new_hk_arr[j, 0, i] = np.nan
        
        #   check the max amount of layers
        if new_ibound_arr.shape[0] > 999:
            #   also cut all the arrays and lists so that it corresponds to the new bottom
            y_lst_res_new = y_lst_res_new[:999]               
            new_end_col = model.bot_elev.index(np.array(model.bot_elev)[np.array(model.bot_elev) > y_lst_res_new[-1]].min())  
            new_ibound_arr = new_ibound_arr[:999, :, :new_end_col]
            new_sconc_arr = new_sconc_arr[:999, :, :new_end_col]
            new_strt_arr = new_strt_arr[:999, :, :new_end_col]
            new_hk_arr = new_hk_arr[:999, :, :new_end_col]
            x_lst_res = x_lst_res[:new_end_col]
            #   remove negative values from conc array
            new_sconc_arr[new_sconc_arr < 0] = 0.0
            new_strt_arr[new_strt_arr <= -999] = 0.0
            model.ncol = new_ibound_arr.shape[-1]
            model.nlay = new_ibound_arr.shape[0]
            model.top_elev = model.top_elev[:new_end_col]
            new_bot_elev = new_bot_elev[:new_end_col]
            model.wtd_elev = model.wtd_elev[:new_end_col]

        #   assign the new arrays to the old ones
        model.ibound_arr = new_ibound_arr
        model.hk_arr = new_hk_arr
        sc_arr = new_sconc_arr
        st_arr = new_strt_arr
        model.bot_elev = new_bot_elev.tolist()
        y_lst_res = y_lst_res_new
        #model.top_elev = topo_res[:idx_cst_res] + model.top_elev[idx_cst_res:]
        
        #   change the top and botm lists so they match the new resolution
        model.top = y_lst_res[0]
        model.botm = y_lst_res[1:].tolist()
        model.botm.append(y_lst_res[-1] - del_lay_res)
        model.lay_elev = y_lst_res.tolist()
    
        #   there might be some values that are non integer in the IBOUND array
        new_ibound_arr[new_ibound_arr > 0] = 1

    #   now create the DIS package again
    model.dis_input(nrow, delc, del_col, del_lay, 1, [perlen_lst[0]], [nstp_lst[0]], laycbd, max_depth, cs_points_dist)
    if model.ibound_arr.shape[0] > 999:
        model.botm = model.botm[:999]
        model.zbot = model.botm[-1]
    
    #   check that the shape of the IBOUND array is the same as the end_idx and end_lay_idx - if not trim all the arrays and lists so that they have the same size, this happens when 
    #   the model domain stretches further than the maximum depth specified.
    if model.idx_end != model.ibound_arr.shape[-1]:
        model.ibound_arr = model.ibound_arr[:model.nlay_end_idx, :, :model.idx_end]
        model.hk_arr = model.hk_arr[:model.nlay_end_idx, :, :model.idx_end]
        model.vk_arr = model.vk_arr[:model.nlay_end_idx, :, :model.idx_end]
        sc_arr = sc_arr[:model.nlay_end_idx, :, :model.idx_end]
        st_arr = st_arr[:model.nlay_end_idx, :, :model.idx_end]
        x_lst_res = x_lst_res[:model.idx_end]
        model.top_elev = model.top_elev[:model.idx_end]
        model.bot_elev = model.bot_elev[:model.idx_end]
        model.x_end = model.idx_end / (1000. / del_col_res) + model.x_start
    
        model.ncol = model.ibound_arr.shape[-1]
        model.delr = model.ncol * [del_col_res]
        model.del_col = del_col_res
        model.ibound_arr = model.ibound_arr.astype(int)
        model.lay_elev = y_lst_res.tolist()
        model.top = model.lay_elev[0]
        model.zbot = model.lay_elev[-1] + del_lay

    if model.nlay_end_idx != model.ibound_arr.shape[0]:
        model.ibound_arr = model.ibound_arr[:model.nlay_end_idx, :, :]
        model.hk_arr = model.hk_arr[:model.nlay_end_idx, :, :]
        model.vk_arr = model.vk_arr[:model.nlay_end_idx, :, :]
        sc_arr = sc_arr[:model.nlay_end_idx, :, :]
        st_arr = st_arr[:model.nlay_end_idx, :, :]
        y_lst_res = y_lst_res[:model.nlay_end_idx]
        model.lay_elev = model.lay_elev[:model.nlay_end_idx]
        model.x_end = model.idx_end / (1000. / del_col_res) + model.x_start
    
        model.ncol = model.ibound_arr.shape[-1]
        model.delr = model.ncol * [del_col_res]
        model.del_col = del_col_res
        model.ibound_arr = model.ibound_arr.astype(int)
        model.lay_elev = y_lst_res.tolist()
        model.top = model.lay_elev[0]
        model.zbot = model.lay_elev[-1] + del_lay
        
    #   loop through each column, check where top elevation is lower than sea level 
    #   first find the column which will act as the coast, go from the offshore boundary and once the top elev is higher than the sea level stop
    cst_col = len(model.top_elev)
    for f in reversed(model.top_elev):
        if f < -130.:
            cst_col -= 1
        else:
            break
    for g in range(cst_col - 1, len(model.top_elev)):
    #for g in range(model.ibound_arr.shape[-1]):
        if model.top_elev[g] <= sea_levels[0]:
            #   get the uppermost active layer of the column
            try:
                top_lay = model.ibound_arr[:, 0, g].tolist().index(1)
                sc_arr[top_lay, 0, g] = 35.0
            except (IndexError, ValueError):
                pass
            
    #   check for nan values in active cells, set negative and too large values to nan first 
    sc_arr[np.abs(sc_arr) >= 999.] = np.nan
    sc_arr[sc_arr < -0.5] = np.nan   
    for a in range(model.ibound_arr.shape[0]):
        for b in range(model.ibound_arr.shape[-1]):
            if (math.isnan(sc_arr[a, 0, b]) or sc_arr[a, 0, b] > 100) and model.ibound_arr[a, 0, b] == 1:
                closest_vals = sc_arr[max(a - 2, 0) : min(a + 2, model.ibound_arr.shape[0]), 0, max(b - 2, 0) : min(b + 2, model.ibound_arr.shape[-1])]
                closest_val = np.nanmean(closest_vals)
                sc_arr[a, 0, b] = closest_val
    #   check for nan values in active cells 
    st_arr[np.abs(st_arr) <= -999.] = np.nan
    st_arr[st_arr == -1.0] = np.nan
    #   inland part - set heads to elevation
    for a in range(0, cst_col):
        for b in range(st_arr.shape[0]):
            if model.ibound_arr[b, 0, a] == 1:
                st_arr[b, 0, a] = model.top_elev[a]
    
    for a in range(model.ibound_arr.shape[0]):
        for b in range(model.ibound_arr.shape[-1]):
            if (st_arr[a, 0, b] <= -999. or math.isnan(st_arr[a, 0, b])) and model.ibound_arr[a, 0, b] == 1:
                closest_vals = st_arr[max(a - 2, 0) : min(a + 2, model.ibound_arr.shape[0]), 0, max(b - 2, 0) : min(b + 2, model.ibound_arr.shape[-1])]
                closest_vals[closest_vals <= -999.] = np.nan
                closest_val = np.nanmean(closest_vals)
                del closest_vals
                st_arr[a, 0, b] = closest_val                            
    
    #   assign the starting sconc and strt arrays
    sc_arr[np.abs(sc_arr) >= 999.] = np.nan
    sc_arr[sc_arr < -0.5] = np.nan          
    sc_arr[sc_arr > 35.0] = 35.0
    
    model.sconc_arr = sc_arr
    model.plot_conc_profile_yrs(sc_dir, x_lst_res, sc_arr, sea_levels[0], 0, 0, zoom = [-10, 10], conc_diff = False)   
 
    #   one last change of the strt array, now check that the heads are not lower than elevation in the column
    for i in range(st_arr.shape[-1]):
        ibound_act_lay_idxs = [k for k, x in enumerate(st_arr[:,0,i].tolist()) if x < model.top_elev[i]]    
        for lay in ibound_act_lay_idxs:
            st_arr[lay, 0, i] = model.top_elev[i]
    
    st_arr[np.abs(st_arr) <= -999.] = np.nan
    #st_arr[st_arr == -1.0] = np.nan
    #   also set heads that are lower than -130m to -130m
    st_arr[st_arr < sea_levels[0]] = sea_levels[0]
    model.strt_arr = st_arr
    model.plot_head_profile_yrs(sc_dir, x_lst_res, st_arr, sea_levels[0], 0, 0, None, None, zoom = [-10, 10], head_diff = False)   

    #   one more final check for nan values
    if math.isnan(np.sum(sc_arr)):                    
        sc_arr = np.nan_to_num(sc_arr, copy=True, nan = 0.0, posinf = None, neginf = None)
    if math.isnan(np.sum(st_arr)):                    
        st_arr = np.nan_to_num(st_arr, copy=True, nan = 0.0, posinf = None, neginf = None)            
    
    #   make sure that there are no nan values in the HK array after the refinement/downscaling
    for a in range(model.hk_arr.shape[0]):
        #   get active cells
        col_act_cells = [c for c, x in enumerate(model.ibound_arr[a, 0, :].tolist()) if x == 1]   
        #   check that at those coordinates the HK vals are not nan
        for idx in col_act_cells:
            hk_cell = model.hk_arr[a, 0, idx]
            if math.isnan(hk_cell):
                #print(a, 0, idx) 
                hk_val = round(np.nanmean(model.hk_arr[a - 1 : a + 1, :, idx - 2 : idx + 2]), 4)
                if math.isnan(hk_val):
                    hk_val = round(np.nanmean(model.hk_arr), 4)
                model.hk_arr[a, 0, idx] = hk_val

    #   check the HK array and the location of potential GHB cells - make sure that there are no large differences between neighbouring cells - causes non-convergence
    #       1) the firs inland column
    for a in range(model.ibound_arr.shape[0]):
        if model.ibound_arr[a, 0, 0] == 1:            
            neighb_arr = model.hk_arr[max(0, a-2) : min(a + 2,model.ibound_arr.shape[0]), :, :2]
            nanmean = np.nanmean(neighb_arr)
            hk_val = model.hk_arr[a, 0, 0]
            #   only change the value if there are at least 2 orders of magnitude difference
            if hk_val / nanmean > 100:
                model.hk_arr[a, 0, 0] = nanmean
    #       2) cells in the top of each column
    for b in range(1, model.ibound_arr.shape[-1]):
        try:
            ibound_act_lay_idx = [k for k, x in enumerate(model.ibound_arr[:, 0, b].tolist()) if x == 1][0] 
            neighb_arr = model.hk_arr[max(0, ibound_act_lay_idx - 1) : min(ibound_act_lay_idx + 1, model.ibound_arr.shape[0]), :, max(0, b - 1) : min(model.ibound_arr.shape[-1], b + 1)]
            neighb_arr[neighb_arr == 0] = np.nan
            nanmean = np.nanmean(neighb_arr)
            hk_val = model.hk_arr[ibound_act_lay_idx, 0, b]        
            if hk_val / nanmean > 10 or hk_val / nanmean < 0.1:
                model.hk_arr[ibound_act_lay_idx, 0, b] = nanmean
        except IndexError:
            pass
    #       3) in the cells that will be DRN check that the conductivity is same in the column (top 3 cells)
    for b in range(model.ibound_arr.shape[-1]):
        if model.top_elev[b] > sea_levels[0]:
            ibound_act_lay_idx = [k for k, x in enumerate(model.ibound_arr[:, 0, b].tolist()) if x == 1]
            if len(ibound_act_lay_idx) > 0:
                neighb_col = model.hk_arr[ibound_act_lay_idx[0] + 1 : min(ibound_act_lay_idx[0] + 5, model.ibound_arr.shape[0]), :, b]
                neighb_col[neighb_col == 0] = np.nan
                nanmean = np.nanmean(neighb_col)
                hk_val = model.hk_arr[ibound_act_lay_idx[0], 0, b]                
                if hk_val / nanmean > 10 or hk_val / nanmean < 0.1:
                    model.hk_arr[ibound_act_lay_idx[0], 0, b] = nanmean

    #       4) look for differences in cell values next to each other - smoothen them out if necessary
    #       also, look for too thin layers, go column by column
    """
    for b in range(model.ibound_arr.shape[-1]):
        col_vals = model.hk_arr[:, 0 , b].tolist()
        #   go through the list and if the neighboring cells hk ratio is higher than 10 (?) adapt the value
        for k in range(len(col_vals) - 1):
            #   has to be non equal to 0.
            if col_vals[k] != 0. and col_vals[k + 1] != 0.:
                curr_val = col_vals[k]
                next_val = col_vals[k + 1]
                ratio = round(curr_val / next_val, 2)
                if ratio > 10.:
                    col_vals[k + 1] = curr_val / 100.
                elif ratio < 0.01:
                    col_vals[k + 1] = curr_val * 10.
                    
        #   now change the values in the array
        for k in range(len(col_vals)):
            model.hk_arr[k, 0 , b] = col_vals[k]
    """
    
    #       5) check individual cells that are surrounded by cells with values that are too high/low
    for b in range(1, model.ibound_arr.shape[-1]):
        for c in range(1, model.ibound_arr.shape[0]):
            cell_val = model.hk_arr[c, 0 ,b]
            if cell_val != 0.:
                #   check the surrounding values left, right, top, bottom
                if b > 0:
                    left_cell = model.hk_arr[c, 0 ,b - 1]
                else:
                    left_cell = 0.
                if b < model.ibound_arr.shape[-1] - 1:
                    right_cell = model.hk_arr[c, 0 ,b + 1]
                else:
                    right_cell = 0.
                if c > 0:
                    top_cell = model.hk_arr[c - 1, 0 ,b]
                else:
                    top_cell = 0.                    
                if c < model.ibound_arr.shape[0] - 1:
                    bot_cell = model.hk_arr[c + 1, 0 ,b]
                else:
                    bot_cell = 0.              
                #   check the ratios
                ratio_lst = [left_cell / cell_val, right_cell / cell_val, top_cell / cell_val, bot_cell / cell_val]
                high_ratio = [i for i in ratio_lst if i > 100.]
                low_ratio = [i for i in ratio_lst if i < .0001]
                
                if len(high_ratio) > 1 and cell_val > 0.0001:
                    model.hk_arr[c, 0 ,b] = model.hk_arr[c, 0 ,b] * 10.
                    #print(b, c, model.hk_arr[c, 0 ,b], ratio_lst)
                    
                if len(low_ratio) > 1 and cell_val > 0.0001:
                    #print(b, c, model.hk_arr[c, 0 ,b], ratio_lst)  
                    model.hk_arr[c, 0 ,b] = model.hk_arr[c, 0 ,b] / 1000.
    
    model.vk_arr = model.hk_arr * 0.1

    #   make sure there are no nan values in the arrays
    model.hk_arr = np.nan_to_num(model.hk_arr, copy = True, nan = 0.0, posinf = None, neginf = None) 
    model.vk_arr = np.nan_to_num(model.vk_arr, copy = True, nan = 0.0, posinf = None, neginf = None) 
    model.sconc_arr = np.nan_to_num(model.sconc_arr, copy = True, nan = 0.0, posinf = None, neginf = None) 
    model.strt_arr = np.nan_to_num(model.strt_arr, copy = True, nan = 0.0, posinf = None, neginf = None) 

    if model.botm[0] != model.lay_elev[1]:
        model.botm = model.lay_elev[1:]
        model.botm.append(model.lay_elev[-1] - del_lay_res)
        
    #   save all the input arrays and lists into a numpy dictionary 
    out_dict_name = os.path.join(sc_dir, '_INIT_dict.npy')
    dict_out = dict()
    dict_out['dis_info'] = {'ibound_arr': model.ibound_arr,\
                            'top_elev': model.top_elev,\
                            'bot_elev': model.bot_elev,\
                            'lay_elev': model.lay_elev,\
                            'hk_arr': model.hk_arr,\
                            'vk_arr': model.vk_arr,\
                            'botm': model.botm,\
                            'top': model.top,\
                            'zbot': model.botm[-1],\
                            'sconc_arr': model.sconc_arr,\
                            'strt_arr': model.strt_arr,\
                            'nlay': model.nlay,\
                            'nrow': model.nrow,\
                            'nper': model.nper,\
                            'delr': model.delr,\
                            'delc': model.delc,\
                            'laycbd': model.laycbd,\
                            'del_lay': model.del_lay,\
                            'del_col': model.del_col,\
                            'ncol': model.ncol,\
                            'idx_end': model.idx_end,\
                            'nlay_end_val': model.nlay_end_val,\
                            'nlay_end_idx': model.nlay_end_idx,\
                            'end_idx': model.end_idx,\
                            'x0': 0,\
                            'x1': model.ibound_arr.shape[-1] * model.del_col,\
                            'wtd_elev': model.wtd_elev,\
                            'cst_offset': model.cst_offset,\
                            'lay_idx': model.lay_idx,\
                            'col_idx': model.col_idx,\
                            'idx_start': model.idx_start,\
                            'x_start': model.x_start,\
                            'x_end': model.x_end,\
                            'soil_thk': model.soil_thk,\
                            'soil_type': model.soil_type,\
                            'cst_idx': model.cst_idx,\
                            'delv': model.delv,\
                            'k_soil': model.k_soil,\
                            'drn_dens': model.drn_dens,\
                            'glhymps_top_lay': model.glhymps_top_lay,\
                            'glhymps_bot_lay': model.glhymps_bot_lay,\
                            'coscat_id': model.coscat_id,\
                            'del_col_res': model.del_col_res,\
                            'del_lay_res': model.del_lay_res,\
                            'x_lst_res' : x_lst_res,\
                            'y_lst_res' : y_lst_res} 
    #   save the dictionary
    np.save(out_dict_name, dict_out)         

    #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
    fig = plt.figure(figsize = (20, 12))
    ax1 = plt.subplot2grid((2, 2), (0, 1))  #   overall concentration profiles
    ax2 = plt.subplot2grid((2, 2), (1, 1))  #   zoomed in concentration profile
    ax1.set_position([0.135, 0.55, 0.85, 0.4]) # [left, bottom, width, height]
    ax2.set_position([0.135, 0.1, 0.85, 0.4])   
    
    #   define the extent of the figures
    x_start_new = model.x_start #-(self.cst_idx - self.idx_start) / 2.
    x_end_new = model.x_end
    y_max_new = model.top
    y_min_new = model.zbot
    
    #   add the location of the aquitard layers    
    mapimg = (model.hk_arr == 1e-04)
    vmin = 0.0
    vmax = 1.0
    #   create output name for the individual figure
    fig_out_dir = os.path.join(sc_dir, '_IBOUND.png')   

    im1 = ax1.imshow(model.ibound_arr[:, 0, :], aspect = 'auto', interpolation = None,
                   extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
    ax1.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                   extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
    
    ax2.imshow(model.ibound_arr[:, 0, :], aspect = 'auto', interpolation = None,
                   extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
    ax2.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                   extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
    
    ax1.set_xlim([x_start_new, x_end_new])
    ax1.set_ylim([y_min_new, y_max_new])
       
    ax2.set_xlim([max(-10, model.x_start), 10])
    ax2.set_ylim([-200, model.top])
    
    #   set the gridlines and constant lines in the plot
    #   add constant lines with elevation = 0m asl. and coastline 
    ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
    ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
    ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
    ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
 
    ax1.axvline(x = 0. + model.cst_offset, linewidth = 2, color = 'k', zorder = 2)
    ax2.axvline(x = 0. + model.cst_offset, linewidth = 2, color = 'k', zorder = 2)    
   
    ax1.axhline(y = 0, linewidth = 2, color = 'k', zorder = 2)    
    ax2.axhline(y = 0, linewidth = 2, color = 'k', zorder = 2)    

    x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new / 5.) * 5., math.ceil(x_end_new), 5.0), nbins = None)    
    x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new / 1.) * 1., math.ceil(x_end_new), 1.0), nbins = None)    
    ax1.xaxis.set_major_locator(x_major_locator)
    ax1.xaxis.set_minor_locator(x_minor_locator)
 
    y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(y_min_new / 500.) * 500., math.ceil(y_max_new), 500.0), nbins = None)    
    y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(y_min_new / 100.) * 100. , math.ceil(y_max_new), 100.0), nbins = None)    
    ax1.yaxis.set_major_locator(y_major_locator)
    ax1.yaxis.set_minor_locator(y_minor_locator)
    
    ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
    ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  

    ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
    ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)        
            
    ax1.set_xlabel('distance from coast (km)', fontsize = 18)
    ax1.set_ylabel('elevation (m asl.)', fontsize = 18)
    ax2.set_xlabel('distance from coast (km)', fontsize = 18)
    ax2.set_ylabel('elevation (m asl.)', fontsize = 18)    

    if len(x_lst_res) > len(model.top_elev):
        x_lst_res = x_lst_res[:len(model.top_elev)]

    ax1.plot(x_lst_res, model.top_elev, c = 'lime', linewidth = 1.0)
    ax2.plot(x_lst_res, model.top_elev, c = 'lime', linewidth = 1.0)
    
    if os.path.exists(fig_out_dir):
        print('Figure   ' + fig_out_dir + '   already exists.')
    else: 
        #   save the figure
        plt.savefig(fig_out_dir, dpi = 300)
        plt.close(fig)
        #del fig
    
    
""" ----------------------------------------------------------------------------------------------- """
"""                 RUN THE ACTUAL MODEL WITHOUT CREATING IBOUND ARRAY IN EACH STEP.                """
""" ----------------------------------------------------------------------------------------------- """
#   if the ibound_array exists then proceed to build the model
if os.path.exists(os.path.join(sc_dir, '_INIT_dict.npy')): 

    #   loop through all the equilibrium periods
    for i in range(len(sea_levels)):
        #   get the stress periods sea level, name and duration
        sea_level = sea_levels[i]
        eq_name = equilibrium_names[i]
        sp_dur = sp_duration[i]
        print (sea_level, eq_name, sp_dur)
        rch_lst_sp_all = []
        yrs_len = yrs_len_lst[i]
        perlen = [perlen_lst[i]]
        nstp = [nstp_lst[i]]
        nper = 1
        
        #   save the output only in middle and end of the stress periods
        if yrs_len == 1000:
            th_btn_time_step_val = [int((yrs_len / 50) + 1)]
        else:
            th_btn_time_step_val = [int((yrs_len / 10) + 1)]
        
        #   define the model name and check if the model folder already exists, if no create it
        modelname = '%s_%s_%s' % (id_cs_str + '_SRM_' + str_id_subreg, topo_type, del_col_res_str + 'm_' + str(int(del_lay_res)) + 'm') 
        topo_sc_dir = os.path.join(out_dir, model_name_basis, 'sc_' + param_combo_str, eq_name)
        
        print (modelname, topo_sc_dir)
        if not os.path.exists(topo_sc_dir):
            os.makedirs(topo_sc_dir, exist_ok = True)

        #   in the first stress period create an overall csv file where we will store the results throughout the whole model run
        if i == 0:
            res_summary_csv_dir = os.path.join(out_dir, model_name_basis, '_res_coscat_summary_%s_%s.csv' % (coscat_id, param_combo))
            if not os.path.exists(res_summary_csv_dir):
                f_res = open(res_summary_csv_dir,'w')
                res_headers = 'eq_name, t_sp_end, frsh_inl_vol, frsh_inl_pct, frsh_shelf_vol, frsh_shelf_pct'                
                f_res.write(" ".join(res_headers.split())) 
                f_res.close()   

        #   define the directories of all the dictionaries
        conc_dict_save_dir = os.path.join(topo_sc_dir, '_conc.npy')
        head_dict_save_dir = os.path.join(topo_sc_dir, '_head.npy')
        cbc_dict_save_dir = os.path.join(topo_sc_dir, '_cbc.npy')    
        ghb_check_dict_save_dir = os.path.join(topo_sc_dir, '_ghb_check.npy')                 
        list_dict_save_dir = os.path.join(topo_sc_dir, '_list.npy')           
        nc_folder = os.path.join(topo_sc_dir, '_nc')            
        #   check if the previous equilibrium was reached, based on the converged.txt
        converged_txt_dir = os.path.join(topo_sc_dir, 'converged.txt')
        csv_dir_ts = os.path.join(topo_sc_dir, '_time_steps_DH_DC.csv')   
        csv_dir_abs_cell_change = os.path.join(topo_sc_dir, '_abs_cell_change.csv')   
        
        #   if the text file exists it means that the convergence was reached and we can go to the next equilibrium stage
        #   load the last sconc and strt arrays as starting arrays for now (might be changed later)
        if os.path.exists(converged_txt_dir):
            print("Model   " + modelname + "   already exists, moving to next EQ.")
            #   read in the head and concentration from the previous simulation (the last time step obviously)
            nc_lst = []
            for filename in os.walk(nc_folder):
                for g in range(len(filename[-1])):
                    if filename[-1][g].endswith(".nc"):
                        print(filename[-1][g].split('_')[-1].split('.nc')[0])
                        nc_lst.append(int(filename[-1][g].split('_')[-1].split('.nc')[0]))
            #   choose the last timestep from the sorted list
            last_ts = sorted(nc_lst)[-1]
            last_nc_name = equilibrium_names[i] + '_' + str(last_ts) + '.nc'
            #   load the nc file and get the head and concentration arrays
            last_nc = xr.open_dataset(os.path.join(os.path.join(out_dir, model_name_basis, 'sc_' + param_combo_str, equilibrium_names[i]), '_nc', last_nc_name))
            sconc_arr = np.expand_dims(np.array(last_nc['solute concentration']), axis = 1)
            strt_arr = np.expand_dims(np.array(last_nc['heads']), axis = 1)
            #rch_lst_sp = last_nc['rch_rate'].values
            sconc_arr[np.isnan(sconc_arr)] = 0.  
            strt_arr[np.isnan(strt_arr)] = 0.            
            #model.sconc_arr = sconc_arr 
            #model.strt_arr = strt_arr
            actual_time = 0     #   model starts at 0 years
            sp = 1
            continue
        #   if the converged.txt doesnt exist, check if there are already the dictionaries created, if yes load the last
        #   time step sconc and strt arrays as starting arrays, also change th sp and actual time values 
        else:
            #if os.path.exists(conc_dict_save_dir) and os.path.exists(head_dict_save_dir):
            if os.path.exists(nc_folder):    
                #   loop through the folder and find all the files ending with .nc
                nc_lst = []
                for filename in os.walk(nc_folder):
                    for g in range(len(filename[-1])):
                        if filename[-1][g].endswith(".nc"):
                            print(filename[-1][g].split('_')[-1].split('.nc')[0])
                            nc_lst.append(int(filename[-1][g].split('_')[-1].split('.nc')[0]))
                #   choose the last timestep from the sorted list
                try:
                    last_ts = sorted(nc_lst)[-1]
                    last_nc_name = eq_name + '_' + str(last_ts) + '.nc'
                    #   load the nc file and get the head and concentration arrays
                    last_nc = xr.open_dataset(os.path.join(nc_folder, last_nc_name))
                    sconc_arr = np.expand_dims(np.array(last_nc['solute concentration']), axis = 1)
                    strt_arr = np.expand_dims(np.array(last_nc['heads']), axis = 1)
                    sconc_arr[np.isnan(sconc_arr)] = 0.  
                    strt_arr[np.isnan(strt_arr)] = 0.                    
                    #rch_lst_sp = last_nc['rch_rate'].values         
                    #   assign new values for the sp and the time parameter
                    sp = int(last_ts / yrs_len) + 1
                    actual_time = last_ts
                    
                #   in case the _nc folder is empty (the simulation could have ended at the first SP without any output to be read)
                #   so just pass this and reset the sp and time counters while the starting arrays will be the same as the last converged ones..
                except IndexError:
                    sp = 1
                    actual_time = 0
                    pass
                
                time_start_check = False        
                
            else:
                actual_time = 0     #   model starts at 0 years
                sp = 1                                        
                f = open(csv_dir_ts,'w')
                f_headers = 'SP, time_start, time_end, max_diff_head, min_diff_head, max_diff_conc, min_diff_conc, fresh_cell_diff_(end_start),\
                            brackish_cell_diff_(end_start), saline_cell_diff_(end_start), frsh_inl_vol,\
                            frsh_inl_pct, frsh_shelf_vol, frsh_shelf_pct, , prev_tot_cells, tot_cells, ghb_cells_check_fail, ghb_cells_check_succ,\
                            ghb_cell_dif_abs_max,  riv_in, rch_in, ghb_in, tot_in, riv_out, drn_out, ghb_out, tot_out, run_time (s)'
                f.write(" ".join(f_headers.split())) 
                f.close()     
                if eq_name != 'BP_30000_to_20000' and sp == 1:
                    sconc_arr[np.isnan(sconc_arr)] = 0.
                    strt_arr[np.isnan(strt_arr)] = 0.
                pass          

        time_start_check = True         
        time_start = time.time()
        #   start with converged as False..
        converged = False
                
        while converged is False:# and sp < 2:
            #   create the model object and the top and bottom elevation lists 
            model = cs_model(None, top_elev, None, None, None, None, None, None, None, None, None, None, None, None, None, None,\
                             None, None, None, False)

            model_dict_in = np.load(os.path.join(sc_dir, '_INIT_dict.npy'), allow_pickle = True).item()
            model.ibound_arr = model_dict_in['dis_info']['ibound_arr']                                   
            model.botm = model_dict_in['dis_info']['botm']                                     
            model.top = model_dict_in['dis_info']['top']                           
            model.lay_elev = model_dict_in['dis_info']['lay_elev'] 
            model.top_elev = model_dict_in['dis_info']['top_elev'] 
            model.bot_elev = model_dict_in['dis_info']['bot_elev'] 
            model.hk_arr = model_dict_in['dis_info']['hk_arr']
            model.vk_arr = model_dict_in['dis_info']['vk_arr']
            model.zbot = model_dict_in['dis_info']['zbot'] 
            model.nlay = model_dict_in['dis_info']['nlay'] 
            model.nrow = model_dict_in['dis_info']['nrow'] 
            model.nper = model_dict_in['dis_info']['nper'] 
            model.delr = model_dict_in['dis_info']['delr']
            model.delc = model_dict_in['dis_info']['delc']
            model.laycbd = model_dict_in['dis_info']['laycbd'] 
            model.del_lay = model_dict_in['dis_info']['del_lay_res']  #model_dict_in['dis_info']['del_lay'] 
            model.del_col = model_dict_in['dis_info']['del_col_res']  #model_dict_in['dis_info']['del_col'] 
            model.ncol = model_dict_in['dis_info']['ncol'] 
            model.idx_end = model_dict_in['dis_info']['idx_end'] 
            model.nlay_end_val = model_dict_in['dis_info']['nlay_end_val'] 
            model.nlay_end_idx = model_dict_in['dis_info']['nlay_end_idx'] 
            model.end_idx = model_dict_in['dis_info']['end_idx'] 
            model.x0 = model_dict_in['dis_info']['x0']
            model.x1 = model_dict_in['dis_info']['x1']    
            model.wtd_elev = model_dict_in['dis_info']['wtd_elev'] 
            model.cst_offset = model_dict_in['dis_info']['cst_offset']
            model.lay_idx = model_dict_in['dis_info']['lay_idx']
            model.col_idx = model_dict_in['dis_info']['col_idx'] 
            model.idx_start = model_dict_in['dis_info']['idx_start'] 
            model.x_start = model_dict_in['dis_info']['x_start']
            model.x_end = model_dict_in['dis_info']['x_end']
            model.soil_thk = model_dict_in['dis_info']['soil_thk'] 
            model.soil_type = model_dict_in['dis_info']['soil_type'] 
            model.cst_idx = model_dict_in['dis_info']['cst_idx']
            model.delv = model_dict_in['dis_info']['delv']
            model.k_soil = model_dict_in['dis_info']['k_soil'] 
            model.drn_dens = model_dict_in['dis_info']['drn_dens'] 
            model.glhymps_top_lay = model_dict_in['dis_info']['glhymps_top_lay']
            model.glhymps_bot_lay = model_dict_in['dis_info']['glhymps_bot_lay'] 
            model.coscat_id = model_dict_in['dis_info']['coscat_id']
            model.del_col_res = model_dict_in['dis_info']['del_col_res'] 
            model.del_lay_res = model_dict_in['dis_info']['del_lay_res'] 
            model.perlen = perlen
            model.nstp = nstp
            model.laytyp = laytyp_val
            x_lst_res = model_dict_in['dis_info']['x_lst_res'] 
            y_lst_res = model_dict_in['dis_info']['y_lst_res'] 

            #   make sure to change the elevation DEM type, but if we do it like below it also maintains the same extent of the model domain
            if topo_100_type != 'gebco':
                topo_lst_sc = [round(i, 2) for i in topo_nc['gebco_' + topo_100_type + '_avg'].values.tolist()]
            else:
                topo_lst_sc = [round(i, 2) for i in topo_nc['gebco_avg'].values.tolist()]
            #   calculate new index 
            mid_idx = int(2000 + model.cst_offset * 10)
            topo_lst_new = topo_lst_sc[mid_idx - model.cst_idx : mid_idx - model.cst_idx + len(model.top_elev)]
            
            #   now connect the two lists into one
            topo_final = topo_lst_new[: model.cst_idx] + model.top_elev[model.cst_idx :]
            model.top_elev = topo_final

            #plt.imshow(model.hk_arr[:,0,:])len()

            #   check the HK array for cells with largely different HK values - non convergence occurs when high and low permeable cells are too close to each other
            #   1) all cells with HK larger than 30 m/d will be limited to 30 m/d
            
            model.hk_arr[model.hk_arr > 20.] = 20.
            """
            #   2) go through the HK array and try to find cells that have a HK ratio difference larger than a defined limit
            hk_diff = 100.
            hk_new = 10.
            
            for r in range(model.hk_arr.shape[0]):
                for s in range(model.hk_arr.shape[-1]):
                    #   only check cells that have higher HK than 1 m/d
                    if model.hk_arr[r, 0, s] > 10.:
                        #   select the surrounding cells - 1 row up and below, 1 column left and right - if possible
                        surr_cells = model.hk_arr[max(0, r - 1) : min(r + 1, model.hk_arr.shape[0]) + 1, 0, max(0, s - 1) : min(s + 1, model.hk_arr.shape[-1]) + 1]
                        #   if there is a cell that has HK value 100 times lower than the i,j cell change the value of the HK cell to hk_new
                        for a in range(surr_cells.shape[0]):
                            for b in range(surr_cells.shape[-1]):
                                if surr_cells[a, b] < model.hk_arr[r, 0, s] / hk_diff and surr_cells[a, b] != 0.:
                                    #print(model.hk_arr[r - 1 + a, 0, s - 1 + b], model.hk_arr[r, 0, s] / 10.)
                                    model.hk_arr[r - 1 + a, 0, s - 1 + b] = hk_new#model.hk_arr[r, 0, s] / 10.
            """
            model.vk_arr = model.hk_arr * 0.1
                       
            #   only if it is the start of the first stress period load the pre-defined sconc and strt
            if eq_name == 'BP_30000_to_20000' and sp == 1:
                model.sconc_arr = model_dict_in['dis_info']['sconc_arr']
                model.strt_arr = model_dict_in['dis_info']['strt_arr']
            else:
                model.sconc_arr = sconc_arr
                model.strt_arr = strt_arr

            #   all the cells that represent the sea floor will be set to saline before starting simulation of each SP
            cst_col = len(model.top_elev)
            for f in reversed(model.top_elev):
                if f < sea_level:
                    cst_col -= 1
                else:
                    break
            for a in range(cst_col - 1, model.sconc_arr.shape[-1]):
                ibound_act_lay_idxs = [k for k, x in enumerate(model.ibound_arr[:, 0, a].tolist()) if x == 1]
                #   select the top elevation to see if its lower than sea level, if so turn the cell to saline
                if len(ibound_act_lay_idxs) > 0:
                    if model.top_elev[a] < sea_level:
                        model.sconc_arr[ibound_act_lay_idxs[0], 0, a] = 35.0
                
            if eq_name == 'BP_30000_to_20000' and sp == 1:
                model.strt_arr[:,:,cst_col:] = sea_level
                   
            #   find the shelf edge
            if init_cond != 'fresh':
                cont_shelf_edge = mgt.find_shelf_break(model, model.top_elev, True) #   try to find the shelf break     
            else:
                cont_shelf_edge = mgt_frsh.find_shelf_break_fresh(model, model.top_elev, True) #   try to find the shelf break     
            
            #   make sure the position of the coastline is correct - in the INIT dictionary it is 
            if eq_name == 'BP_30000_to_20000':
                idx_cst_0 = int(abs(round(model.x_start + model.cst_offset, 2)) * 10)
                idx_shelf_edge_0 = int((cont_shelf_edge[0] - model.x_start + model.cst_offset) * 10)                
                if idx_shelf_edge_0 > model.ibound_arr.shape[-1]:
                    idx_shelf_edge_0 = model.ibound_arr.shape[-1] - 1
            else:
                #   open the text file to find the position of the coast and cont shelf at the EQ 0  
                txt = open(os.path.join(os.path.join(out_dir, model_name_basis, 'sc_' + param_combo_str, 'BP_30000_to_20000'), 'converged.txt'),'r')
                lines = txt.readlines()
                idx_cst_0 = int(lines[0].strip('\n'))
                idx_shelf_edge_0 = int(lines[1])
                if idx_shelf_edge_0 > model.ibound_arr.shape[-1]:
                    idx_shelf_edge_0 = model.ibound_arr.shape[-1] - 1
                    
            #   the hk and vk arrays are already created so no need to create them again, just write the parameter values in to the csv file    
            if time_start_check:
                #   write into the geo_summary_csv file
                if i == 0:                  
                    f_geo = open(geo_summary_csv_dir,'a')
                    f_geo.write("\n")
                    f_geo.write(str(param_combo) + ',' + sed_type + ',' + str(sand_lay_n) + ',' + str(mud_pct) + ',' + str(y_val) + ',' + str(clay_cap_shelf)\
                                + ',' + str(clay_cap_slope) + ',' + str(clay_cap_shelf_thk) + ',' + str(clay_cap_slope_thk) + ',' + str(p_fact) + ',' + str(off_lay_start) + ',' + lay_thk_str) 
                    f_geo.close()       
                    #   create a list where we are going to store the results for each EQ period and use it to write into the final res_summary_csv
                    eq_res_lst = []

            #   specify input for the GHB package
            cond_factor = 1.
            clay_incr = True # means that there will not be any highering up of conductance if the calculated value is too low
                             # this usually happened in the clay layers where the conductivity is low and therefore the conductance too
            cond_limit = 1.            
            clay_cond = 1.
            #   check if the last column offshore has only one active cell, if it is more then also assign GHB to those cells
            last_col_act_cells = [i for i, x in enumerate(model.ibound_arr[:, 0, -1].tolist()) if x == 1]   
            if len(last_col_act_cells) > 1:
                last_col_ghb = True
            else:
                last_col_ghb = False
            
            #   set the time start check
            if actual_time != 0:
                time_start_check = False

            x_lst_100 = np.arange(model.x_start, model.x_end, 100. / 1000.)
            if sp ==1:
                if rch_type == 'paleo':
                    gw_rch_fig_out = os.path.join(topo_sc_dir, '_rch_p_min_et_input.png')
                    gw_rch_lst_100 = mfc.generate_ln_rand_lst(param_combo, gw_rch_mu[i], gw_rch_std[i], len(x_lst_100),\
                                                       'GW RCH input values per column', 'RCH rate in mm/d', gw_rch_fig_out, plot = True) 
                else:
                    gw_rch_fig_out = os.path.join(topo_sc_dir, '_rch_p_min_et_input.png')
                    gw_rch_lst_100 = mfc.generate_ln_rand_lst(param_combo, gw_rch_mu[-1], gw_rch_std[-1], len(x_lst_100),\
                                                       'GW RCH input values per column', 'RCH rate in mm/d', gw_rch_fig_out, plot = True)     
                
            if 'gw_rch_lst_100' not in locals():
                gw_rch_fig_out = os.path.join(topo_sc_dir, '_rch_p_min_et_input.png')
                gw_rch_lst_100 = mfc.generate_ln_rand_lst(param_combo, gw_rch_mu[-1], gw_rch_std[-1], len(x_lst_100),\
                                                   'GW RCH input values per column', 'RCH rate in mm/d', gw_rch_fig_out, plot = True)                     
            
            #   create the recharge input..without the recharge adjusting factor 
            gw_rch_lst_m_d = []
            if del_col != del_col_res:
                gw_rch_lst = np.interp(x_lst_res, x_lst_100, gw_rch_lst_100).tolist()
            else:
                gw_rch_lst = gw_rch_lst_100
            gw_rch_lst_m_d = [round(i * 0.001 / 365., 6) for i in gw_rch_lst]
            gw_rch_lst_m_d = gw_rch_lst_m_d[:model.ibound_arr.shape[-1]]                

            #   create input for all but the DIS and BAS packages
            #model.ghb_input(sea_level, 1., clay_cond, cond_limit, clay_incr, last_col_ghb, True, True, True, wtd_inl_lvl = True)  
            
            #   create the GHB package input here
            itype = flopy.mt3d.Mt3dSsm.itype_dict()
            model.ghb_arr = model.ibound_arr * 1
            model.ghb_input_lst = []
            model.ssmdata = []
            #       1) the inland column and the offshore last column where the CONDUCTANCE is calculated as Hk * del_lay * LD 
            for b in range(model.ibound_arr.shape[0]):
                if model.ibound_arr[b, 0, 0] == 1:
                    model.ghb_input_lst.append([b, 0, 0, model.top_elev[0], model.hk_arr[b, 0, 0] * model.del_lay_res * ld_rate])
                    model.ghb_arr[b, 0, 0] = -1
                    model.ssmdata.append([b, 0, 0, 0.0, itype['GHB']])            
            
            for c in range(model.ibound_arr.shape[0]):
                if model.ibound_arr[c, 0, model.ibound_arr.shape[-1] - 1] == 1:
                    model.ghb_input_lst.append([c, 0, model.ibound_arr.shape[-1] - 1, max(sea_level, model.top_elev[-1]), min(100., model.hk_arr[c, 0, model.ibound_arr.shape[-1] - 1] * model.del_lay_res * ld_rate)])
                    model.ghb_arr[c, 0, model.ibound_arr.shape[-1] - 1] = -1
                    model.ssmdata.append([c, 0, model.ibound_arr.shape[-1] - 1, 35.0, itype['GHB']])      

            #       2) the offshore columns that represent the sea floor, where the CONDUCTANCE is calculated as Hk * del_col * LD          
            for a in range(cst_col, len(model.top_elev)):
                #   get the index of the first non-zero ibound cell
                try:
                    top_cell_lay_idx = model.ibound_arr[:, 0, a].tolist().index(1)
                    topo_val = model.top_elev[a]
                    if topo_val < sea_level:
                        #model.ghb_input_lst.append([top_cell_lay_idx, 0, a, sea_level, model.vk_arr[top_cell_lay_idx, 0, a] * model.del_col_res * ld_rate])   #   sea level is 0.0m (present condition)
                        model.ghb_input_lst.append([top_cell_lay_idx, 0, a, sea_level, min(100., model.vk_arr[top_cell_lay_idx, 0, a] * model.del_col_res * ld_rate)])
                        model.ghb_arr[top_cell_lay_idx, 0, a] = -1
                        model.ssmdata.append([top_cell_lay_idx, 0, a, 35.0, itype['GHB']])
                except ValueError:
                    pass            
                
            #   write the final output dictionary, inlcude each stress period
            model.ghb_arr_in = {}
            for d in range(len(model.perlen)):
                model.ghb_arr_in[d] = model.ghb_input_lst                     
            model.rch_input_from_lst_vk_limit(gw_rch_lst_m_d, sea_level)
            
            #model.drn_input(sea_level, 0.0)
            #model.drn_input_wtd(sea_level)
            
            #  drainage is assigned only to cells that receive recharge - cells with elev above sea level
            drn_input_lst = []
            #for i in range(self.ncol):
            #for i in range(self.ibound_arr.shape[-1]):
            for i in range(len(model.top_elev)):
                ibound_col_lst = model.ibound_arr[:, 0, i].tolist()
                #   check the 1st column with ibound_val = 1 (active cell)
                try:
                    drn_lay = ibound_col_lst.index(1)
                    ##  now check if the elevation is below sea level, if so assign the cell to ghb list
                    if model.top_elev[i] >= sea_level:
                        cond_cell = ((model.hk_arr[drn_lay, 0, i] / 10.) * model.del_col)# / model.del_lay
                        #drn_input_lst.append([drn_lay, 0, i, model.top_elev[i], model.del_col_res * ld_rate / model.del_lay_res])
                        drn_input_lst.append([drn_lay, 0, i, model.top_elev[i], min(10., cond_cell)])
                        #drn_input_lst.append([drn_lay, 0, i, model.top_elev[i], model.vk_arr[drn_lay, 0, i] / model.del_lay_res])
                    else:
                        pass
                except ValueError:
                    pass
            #   write the final output dictionary, inlcude each stress period
            model.drn_arr_in = {}
            for c in range(len(model.perlen)):
                model.drn_arr_in[c] = drn_input_lst            
                
            model.pcg_input(hclose, rclose)
            model.oc_input(1)#(th_oc_time_step_val)
            model.btn_input(porosity, dt0, th_btn_time_step_val, ifmtcn, nprs, chkmas, nprobs, nprmas)
            model.adv_input(mixelm)
            model.dsp_input(al, trpt, trpv, dmcoef)
            model.gcg_input(iter1, mxiter, isolve, cclose)
            model.vdf_input(iwtable, densemin, densemax, denseref, denseslp, firstdt) 
            model.ssm_input()
            
            #   write all the input files
            model.write_all_input(modelname, mf_exe_dir, mt3d_exe_dir, swat_exe_dir, topo_sc_dir)

            #   plot the HK array, just as a double check that it stays the same throughout the simulation
            x_lst_res = np.arange(model.x_start, model.x_start + (model.hk_arr.shape[-1] * 0.1), del_col_res / 1000.)    
            model.plot_HK_arr_res(topo_sc_dir, x_lst_res, model.hk_arr, sea_level, 0, 0, zoom = [-10., 10.])
 
            #   run the model
            model.run_model()
                  
            #   read the model output and create output pictures and video
            model.read_output()
        
            #   check if the plotting directories exist, if yes then delete them
            dest_dir_conc = os.path.join(topo_sc_dir, '_conc_pngs')    
            if not os.path.exists(dest_dir_conc):
                os.makedirs(dest_dir_conc, exist_ok = True)        
            dest_dir_conc_flow = os.path.join(topo_sc_dir, '_conc_flow_pngs')    
            if not os.path.exists(dest_dir_conc_flow):
                os.makedirs(dest_dir_conc_flow, exist_ok = True)                                 
            dest_dir_head_fresh = os.path.join(topo_sc_dir, '_head_fresh_pngs')    
            if not os.path.exists(dest_dir_head_fresh):
                os.makedirs(dest_dir_head_fresh, exist_ok = True)                                   
            nc_dir = os.path.join(topo_sc_dir, '_nc')    
            if not os.path.exists(nc_dir):
                os.makedirs(nc_dir)         
                
            #   conc/head profiles to be read every stepth year and define the convergence criteria
            if yrs_len == 1000:
                step = 500.      
            elif yrs_len == 100:
                step = 50
                
            conc_limit = 0.05       
            head_limit = 0.05                    
            #   read in the number of time steps and transform the time from days to years
            conc_ts_lst, heads_ts_lst = model.time_steps, model.times_heads
            conc_ts_lst_yrs = [round(k / 365.25, 0) for k in conc_ts_lst]
            heads_ts_lst_yrs = [round(k / 365.25, 0) for k in heads_ts_lst]

            x_coord_lst = np.linspace(model.x_start + (model.del_col / 1000.) / 2, model.x_end - (model.del_col / 1000.) / 2, model.ncol)
            x_coord_lst = [round(k, 3) for k in x_coord_lst]
            y_coord_lst = np.linspace(model.top - 5., model.zbot + 5., model.nlay)
            y_coord_lst = [round(k, 3) for k in y_coord_lst]
            cst_idx_0 = int(abs(model.x_start) * 1000. / model.del_col_res)

            #   for each stepth year print a concentration, head profile, if it is the 1st SP then start with the 
            #   profiles at TS = 0yrs..
            if sp == 1:
                cbc_dict = dict()
                ghb_check_dict = dict()
                list_dict = dict()
                conc_arr_prev = model.sconc_arr.astype(dtype = np.float64)
                head_arr_prev = model.strt_arr.astype(dtype = np.float64)        
                zero_arr = model.ibound_arr.astype(dtype = np.float64) * 0.0
                cbc_dict[0] = {'q_ghb': zero_arr, 'q_drn': zero_arr, 'q_rch': zero_arr, 'q_riv' : zero_arr}                       
                ghb_check_dict[0] = {'ghb_check': zero_arr, 'ghb_head_diff_lst': [-1], 'ghb_head_diff_max': [-1], 'succ_cnt': -1, 'fail_cnt': -1}                      
                head_fresh_arr_prev = mfc.pt_hd_to_frsh_hd(head_arr_prev, conc_arr_prev, denseref, denseslp, model.top, model.botm)
                #   for the concentration, heads and cbc create a netcdf file (to save memory)
                xa_sum = xr.Dataset(data_vars = {'solute concentration' : (('y', 'x'), conc_arr_prev[:, 0, :]),
                                                 'heads' : (('y', 'x'),head_arr_prev[:, 0, :]),
                                                 'cbc Q_right' : (('y', 'x'), zero_arr[:, 0, :]),
                                                 'cbc Q_bottom' : (('y', 'x'), zero_arr[:, 0, :])},
                                    coords = {'x' : x_coord_lst,
                                              'y' : y_coord_lst})
                xa_sum = xa_sum.assign_coords(time = 0)
                xa_name = eq_name + '_0.nc'
                xa_sum.to_netcdf(os.path.join(nc_dir, xa_name))                
            else:
                cbc_dict = np.load(cbc_dict_save_dir, allow_pickle = True).item()  
                ghb_check_dict = np.load(ghb_check_dict_save_dir, allow_pickle = True).item()  
                list_dict = np.load(list_dict_save_dir, allow_pickle = True).item()
                #   the recharge adjusting factor in case there are huge head differences in the model = injecting recharge into low permeable systems
                #   if the drainage is larger than recharge then no need to decrease recharge anymore
                #   second condition is that all the heads have to be lower than topography + 0.1m 
                rch_adj_factor_lst = []
                for y in range(model.ncol):
                    if max(strt_arr[:,0,y]) > model.top_elev[y] + 1.0:
                        rch_adj_factor_lst.append(round(1 - (list_dict[list(list_dict.keys())[-1]]['DRAINS_OUT'] / list_dict[list(list_dict.keys())[-1]]['RECHARGE_IN']), 3))
                    else:
                        rch_adj_factor_lst.append(1.)

            if yrs_len == 1000:
                ts_to_plt = conc_ts_lst_yrs[::10]
                #   create a list of time steps to create plots of differences in concentrations and heads 
                ts_diff_to_plt = conc_ts_lst_yrs[::10]
            else:
                ts_to_plt = conc_ts_lst_yrs[::5]
                #   create a list of time steps to create plots of differences in concentrations and heads 
                ts_diff_to_plt = conc_ts_lst_yrs[::5]
                
                #   define a counter for counting the number of time steps where the fresh, brackish, saline are constant
            cnt_volumes_const = 0
            cell_diff_pos_lst, cell_diff_neg_lst = [], []
            
            #   loop through the list and plot the heads and concentration profiles for each of the chosen time steps
            for z in range(1, len(ts_to_plt)):
                #   select the index of the time step to be plotted and get the concentration and heads arrays
                ts_plt_conc = min(conc_ts_lst_yrs, key = lambda x:abs(x - ts_to_plt[z]))
                ts_plt_head = min(heads_ts_lst_yrs, key = lambda x:abs(x - ts_to_plt[z]))
                conc_arr = model.ucnobj.get_data(totim = model.time_steps[conc_ts_lst_yrs.index(ts_plt_conc)]).astype(dtype = np.float64)
                head_arr = model.hdsobj.get_data(totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)]).astype(dtype = np.float64)
                
                print('analyzing results for time = ', str(ts_plt_conc))
                
                #   calculate the freshwater head array
                head_fresh_arr = mfc.pt_hd_to_frsh_hd(head_arr, conc_arr, denseref, denseslp, model.top, model.botm)                    
                qx_in = model.cbbobj.get_data(text='flow right face', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]
                qz_in = model.cbbobj.get_data(text = 'flow lower face', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]                          
                #   create plots for given years
                real_time_yrs = int(actual_time + ts_plt_conc)       
                #   Plotting of 2D concentration profiles at given time
                #   find the location of the coastline for the zoom and plotting the coastline line in the plots
                for y in range(0, len(model.top_elev)):
                    if model.top_elev[y] > sea_level:
                        shift_cst = y / (1000. / del_col_res) + model.x_start
                        zoom_min = round(shift_cst - 10., 2)
                        zoom_max = round(shift_cst + 10., 2)
                        pass
                    else:
                        shift_cst = y / (1000. / del_col_res) + model.x_start 
                        zoom_min = round(shift_cst - 10., 2)
                        zoom_max = round(shift_cst + 10., 2)
                        break
                
                #   plot the arrays
                model.plot_head_profile_yrs(dest_dir_head_fresh, x_lst_res, head_fresh_arr, sea_level, real_time_yrs, shift_cst, None, None, zoom = [zoom_min, zoom_max], head_diff = False)
                model.plot_conc_profile_yrs(dest_dir_conc, x_lst_res, conc_arr, sea_level, real_time_yrs, shift_cst, zoom = [zoom_min, zoom_max], conc_diff = False)                    
                                
                #   get the arrays for the difference between the two consecutive time steps, if the z = 0 there is obviously no index before that so read
                #   in either the starting arrays (SP 0) or the last conc and head arrays from the previous SP
                ts_plt_conc_prev = min(conc_ts_lst_yrs, key = lambda x:abs(x - ts_to_plt[z - 1]))
                ts_plt_head_prev = min(heads_ts_lst_yrs, key = lambda x:abs(x - ts_to_plt[z - 1]))
                #   next calculate the difference arrays and also the number of cell changes between these two steps
                conc_arr_prev = model.ucnobj.get_data(totim = model.time_steps[conc_ts_lst_yrs.index(ts_plt_conc_prev)]).astype(dtype = np.float64)
                conc_arr_diff = conc_arr - conc_arr_prev
                conc_max_diff, conc_min_diff = np.nanmax(conc_arr_diff), np.nanmin(conc_arr_diff)   
                head_arr_prev = model.hdsobj.get_data(totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head_prev)]).astype(dtype = np.float64)
                head_fresh_arr_prev = mfc.pt_hd_to_frsh_hd(head_arr_prev, conc_arr_prev, denseref, denseslp, model.top, model.botm)      
                head_arr_diff = head_arr - head_arr_prev
                head_max_diff, head_min_diff = np.nanmax(head_arr_diff), np.nanmin(head_arr_diff)   
                head_fresh_arr_diff = head_fresh_arr - head_fresh_arr_prev

                #   read the other flow budget information                        
                qdrains = model.cbbobj.get_data(text = '          DRAINS', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]                                 
                qghb = model.cbbobj.get_data(text = ' HEAD DEP BOUNDS', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]
                qrch = model.cbbobj.get_data(text = '        RECHARGE', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]    
                rch_q = model.cbbobj.get_data(text = '        RECHARGE', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])
                drn_q = model.cbbobj.get_data(text = '          DRAINS', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])
                
                #   plot the HK array, just as a double check that it stays the same throughout the simulation
                if sp == 1 and z == 1:
                    if len(x_lst_res) != model.hk_arr.shape[-1]:
                        x_lst_res = np.arange(model.x_start, model.x_start + (model.hk_arr.shape[-1] * 0.1) - 0.1, del_col_res / 1000.)    
                    model.plot_HK_arr_res(topo_sc_dir, x_lst_res, model.hk_arr, sea_level, real_time_yrs, shift_cst, zoom = [-10., 10.])
                
                #   plot the respective arrays and overwrite the ts_plt_conc_prev and ts_plt_head_prev
                #   in non active cells assign the NaN value
                for x in range(conc_arr_diff.shape[0]):
                    for j in range(conc_arr_diff.shape[-1]):
                        if model.ibound_arr[x, 0, j] != 1:
                            conc_arr_diff[x, 0, j] = 999.999
                            head_arr_diff[x, 0, j] = 999.999
                            head_fresh_arr_diff[x, 0, j] = 999.999
                            
                #   get the number of fresh/brackish/saline cells in the model domain at TS_1 and TS_2 and check the change between them
                conc_ts_prev_fresh = ((-100.0 <= conc_arr_prev) & (conc_arr_prev <= 0.5)).sum()
                conc_ts_prev_brackish = ((0.5 < conc_arr_prev) & (conc_arr_prev <= 20.)).sum()
                conc_ts_prev_saline = ((conc_arr_prev > 20.0) & (conc_arr_prev < 100.)).sum()
                conc_ts_fresh = ((-100.0 <= conc_arr) & (conc_arr <= 0.5)).sum()
                conc_ts_brackish = ((0.5 < conc_arr) & (conc_arr <= 20.)).sum()
                conc_ts_saline = ((conc_arr > 20.0) & (conc_arr < 100.)).sum()
                
                prev_tot_cells = conc_ts_prev_fresh + conc_ts_prev_brackish + conc_ts_prev_saline
                tot_cells = conc_ts_fresh + conc_ts_brackish + conc_ts_saline
                
                #   change the values in arrays according to ranges to mark fresh, brackish and saline cells and see where the differences occur
                arr_curr, arr_prev = conc_arr * 1.0, conc_arr_prev * 1.0
                arr_curr[(-100.0 <= arr_curr) & (arr_curr <= 0.5)] = -1
                arr_curr[(0.5 < arr_curr) & (arr_curr <= 20.)] = 0
                arr_curr[(20. < arr_curr) & (arr_curr <= 100.)] = 1
                arr_prev[(-100.0 <= arr_prev) & (arr_prev <= 0.5)] = -1
                arr_prev[(0.5 < arr_prev) & (arr_prev <= 20.)] = 0
                arr_prev[(20. < arr_prev) & (arr_prev <= 100.)] = 1                    
                #   create the difference array and count the number of 1 or -1 cells
                arr_diff = arr_curr - arr_prev
                unique, counts = np.unique(arr_diff, return_counts = True)
                #   append to the list of cells that are different
                if 1 in unique:
                    pos_lst = list(np.where(arr_diff == 1.))
                    for l in range(len(pos_lst[0])):
                        cell_diff_pos_lst.append([pos_lst[0][l], pos_lst[1][l], pos_lst[2][l]])
                if -1 in unique:
                    neg_lst = np.where(arr_diff == -1.)
                    for l in range(len(neg_lst[0])):
                        cell_diff_neg_lst.append([neg_lst[0][l], neg_lst[1][l], neg_lst[2][l]])                    
                del arr_curr, arr_prev
                
                #   get the differences in number of cells of each category
                fresh_change, brack_change, saline_change = conc_ts_fresh - conc_ts_prev_fresh, conc_ts_brackish - conc_ts_prev_brackish, conc_ts_saline - conc_ts_prev_saline
                #   check if all the changes are 0, if yest increase the counter
                if fresh_change == 0 and brack_change == 0 and saline_change == 0:
                    cnt_volumes_const += 1
                #   append the results of the EQ to the final result list 
                if idx_shelf_edge_0 > model.ibound_arr.shape[-1]:
                    idx_shelf_edge_0 = model.ibound_arr.shape[-1] - 1
                    
                #   find the coastal index - taking current sea level as position of the coastline
                idx_shelf_edge_0 = int((cont_shelf_edge[0] - model.x_start + model.cst_offset) * (1000. / del_col_res))
                #volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx_0, idx_shelf_edge_0)
                volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx, idx_shelf_edge_0)
                
                #   write the results into the csv
                f_res = open(res_summary_csv_dir,'a')
                f_res.write("\n")
                f_res.write(eq_name + ',' + str(real_time_yrs) + ',' + str(volumes[1]) + ',' + str(round(100. * (volumes[1] / volumes[0]), 1))  + ',' +\
                            str(volumes[3]) + ',' +  str(round(100. * (volumes[3] / volumes[2]), 1))) 
                f_res.close()   
                list_dir = os.path.join(topo_sc_dir, modelname + '.list')
                swt_list = flopy.utils.mflistfile.SwtListBudget(list_dir)
                budget = swt_list.get_budget()
                #   check the GHB condition
                check_ghb = model.check_ghb_head(heads_ts_lst_yrs.index(ts_plt_head), 0.0001)
                #   write the results into the output CSV file
                f = open(csv_dir_ts,'a')
                f.write("\n")
                f.write(str(sp) + ',' + str(real_time_yrs - step) + ',' + str(real_time_yrs) + ',' + str(head_max_diff) + ',' + str(head_min_diff)\
                        + ',' + str(conc_max_diff) + ',' + str(conc_min_diff) + ',' + str(conc_ts_fresh - conc_ts_prev_fresh) + ',' +\
                        str(conc_ts_brackish - conc_ts_prev_brackish) + ',' + str(conc_ts_saline - conc_ts_prev_saline) + ',' +\
                        str(volumes[1]) + ',' + str(round(100. * (volumes[1] / volumes[0]), 1))  + ',' +\
                        str(volumes[3]) + ',' +  str(round(100. * (volumes[3] / volumes[2]), 1)) + ',' +\
                        str(prev_tot_cells) + ',' +  str(tot_cells) + ',' + str(model.ghb_fail) + ',' + str(model.ghb_succ) + ',' +\
                        str(round(model.ghb_head_check_max_diff, 5)) + ',' + str(budget[0][0][7])  + ',' + str(budget[0][0][6])  + ',' +\
                        str(budget[0][0][9]) + ',' + str(budget[0][0][12]) + ',' + str(budget[0][0][13]) + ',' + str(budget[0][0][16])  + ',' +\
                        str(model.run_time))
                f.close() 

                #if z % 5 == 0:
                    
                #   for the concentration, heads and cbc create a netcdf file (to save memory)
                xa_sum = xr.Dataset(data_vars = {'solute concentration' : (('y', 'x'), conc_arr[:, 0, :]),
                                                 'heads' : (('y', 'x'),head_arr[:, 0, :]),
                                                 'cbc Q_right' : (('y', 'x'), qx_in[:, 0, :]),
                                                 'cbc Q_bottom' : (('y', 'x'), qz_in[:, 0, :]),
                                                 'rch_rate' : (('x'), gw_rch_lst_m_d[:])},
                                    coords = {'x' : x_coord_lst,
                                              'y' : y_coord_lst})
                #xa_sum.expand_dims({'time':real_time_yrs})
                xa_sum = xa_sum.assign_coords(time = real_time_yrs)
                xa_name = eq_name + '_' + str(real_time_yrs) + '.nc'
                xa_sum.to_netcdf(os.path.join(nc_dir, xa_name))                         
                
                #   append the conc and head arrays at certain times to the output dictionaries
                ghb_check_dict[real_time_yrs] = {'ghb_check': model.ghb_head_check, 'ghb_head_diff_lst': model.ghb_head_diff_lst,\
                                                 'ghb_head_diff_max': model.ghb_head_check_max_diff, 'succ_cnt': model.ghb_succ, 'fail_cnt': model.ghb_fail}                    
                cbc_dict[real_time_yrs] = {'q_ghb': qghb, 'q_drn': qdrains, 'q_rch': qrch}                
                #   store the flow budget from the list file
                list_dict[real_time_yrs] = {'STORAGE_IN': budget[0][0][3], 'CONSTANT_HEAD_IN': budget[0][0][4], 'DRAINS_IN': budget[0][0][5],\
                                            'HEAD_DEP_BOUNDS_IN': budget[0][0][6], 'RECHARGE_IN': budget[0][0][7], 'DCDT_IN': budget[0][0][8],\
                                            'TOTAL_IN': budget[0][0][9], 'STORAGE_OUT': budget[0][0][10], 'CONSTANT_HEAD_OUT': budget[0][0][11],\
                                            'DRAINS_OUT': budget[0][0][12], 'HEAD_DEP_BOUNDS_OUT': budget[0][0][13], 'RECHARGE_OUT': budget[0][0][14],\
                                            'DCDT_OUT': budget[0][0][15], 'TOTAL_OUT': budget[0][0][16], 'IN-OUT': budget[0][0][17],\
                                            'PERCENT_DISCREPANCY': budget[0][0][18]}      
                """
                rch_adj_factor_lst = []
                for y in range(model.ncol):
                    if max(head_arr[:,0,y]) > model.top_elev[y] + 1.0:
                        rch_adj_factor_lst.append(round(1 - (list_dict[real_time_yrs]['DRAINS_OUT'] / list_dict[real_time_yrs]['RECHARGE_IN']), 3))
                    else:
                        rch_adj_factor_lst.append(1.)                        
                rch_adj_factor = round(1 - (list_dict[real_time_yrs]['DRAINS_OUT'] / list_dict[real_time_yrs]['RECHARGE_IN']), 3)
                """
                
                model.plot_conc_flow_vectors_profile_yrs(dest_dir_conc_flow, conc_arr, sea_level, real_time_yrs, shift_cst, qx_in, qz_in, zoom = [zoom_min, zoom_max], conc_diff = False)

            #   assign the prev arrays to be the last arrays of this SP and add the duration of the SP to the actual time counter
            conc_arr_prev, head_arr_prev = conc_arr, head_arr
            actual_time += (perlen[0] / 365.25)

            #   save the dictionary with the ibound arrays
            #   put the dictionaries in the right order (time wise)
            cbc_in_ordered = collections.OrderedDict(sorted(cbc_dict.items()))     
            ghb_check_in_ordered = collections.OrderedDict(sorted(ghb_check_dict.items()))          
            list_in_ordered = collections.OrderedDict(sorted(list_dict.items()))   
            np.save(cbc_dict_save_dir, cbc_in_ordered)        
            np.save(ghb_check_dict_save_dir, ghb_check_in_ordered)
            np.save(list_dict_save_dir, list_in_ordered)

            #   check if the model converged, convergence is defined as when the difference between concentrations between two SPs is less than one cell
            if cnt_volumes_const < int(perlen[0] / (step * 365.25)) - 1:
                #   now check the list and if it is only one cell changing back and forth 
                if cell_diff_pos_lst != [] or cell_diff_neg_lst != []:
                    
                    #   check which list is largest and then loop through the smaller one
                    if len(cell_diff_pos_lst) >= len(cell_diff_neg_lst):
                        lst_to_loop = [list(x) for x in set(tuple(x) for x in cell_diff_pos_lst)]
                        lst_no_loop = [list(x) for x in set(tuple(x) for x in cell_diff_neg_lst)]
                    else:
                        lst_to_loop = [list(x) for x in set(tuple(x) for x in cell_diff_neg_lst)]
                        lst_no_loop = [list(x) for x in set(tuple(x) for x in cell_diff_pos_lst)]
                    
                    #   loop through the list and for each cell check if it occurs in the other list, first set up counters to count how many cells are 
                    #   really changing and not just oscillating
                    non_loop_cells = len(lst_no_loop)
                    loop_cells = len(lst_to_loop)
                    #   loop through the list and if it is also in the smaller list decrease the total counter
                    for g in range(len(lst_to_loop)):
                        if lst_to_loop[g] in lst_no_loop:
                            non_loop_cells -= 1
                            loop_cells -= 1
                    #   if the counter itself is lower or equal to 1 it is considered a stable situation (change of 1 cell per 1000 years)
                    if loop_cells + non_loop_cells <= 1:
                        converged = True
                    else:
                        converged = False
                else:
                    loop_cells, non_loop_cells = -1, -1
                    converged = True
            else:
                converged = True

            #   then check if the time specific for the eq already expired, if yes then change the converged to True
            if not converged:
                #   only in case the time is reached
                if sp_dur == -1:
                    pass
                elif actual_time < sp_dur:
                    pass
                else:
                    converged = True

            #   check if the model converged or not
            if not converged:
                sconc_arr = model.ucnobj.get_data(totim = model.time_steps[conc_ts_lst_yrs.index(ts_plt_conc)]).astype(dtype = np.float64)
                strt_arr = model.hdsobj.get_data(totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)]).astype(dtype = np.float64)
                                
                sp += 1  
                gc.collect()
                #break
                
            else:
                #   save the dictionary with the ibound arrays
                cbc_dict_save_dir = os.path.join(topo_sc_dir, '_cbc.npy')    
                ghb_check_dict_save_dir = os.path.join(topo_sc_dir, '_ghb_check.npy')                 
                list_dict_save_dir = os.path.join(topo_sc_dir, '_list.npy')   
                #   put the dictionaries in the right order (time wise)
                cbc_in_ordered = collections.OrderedDict(sorted(cbc_dict.items()))     
                ghb_check_in_ordered = collections.OrderedDict(sorted(ghb_check_dict.items()))          
                list_in_ordered = collections.OrderedDict(sorted(list_dict.items()))   

                strt_arr = model.hdsobj.get_data(totim = model.times_heads[-1])
                sconc_arr = model.ucnobj.get_data(totim = model.time_steps[-1])                                    
                
                np.save(cbc_dict_save_dir, cbc_in_ordered)        
                np.save(ghb_check_dict_save_dir, ghb_check_in_ordered)
                np.save(list_dict_save_dir, list_in_ordered)

                del cbc_dict, ghb_check_dict, cbc_in_ordered, ghb_check_in_ordered,\
                list_dict_save_dir, list_in_ordered
                
                time_end = time.time()           
                run_time = time_end - time_start                   

                #volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx_0, idx_shelf_edge_0)
                #volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx_0, idx_shelf_edge_0)
                volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx, idx_shelf_edge_0)
                
                #   create converged text file
                txt = open(os.path.join(topo_sc_dir, 'converged.txt'),'w')
                txt.write(str(idx_cst_0) + '\n')
                txt.write(str(idx_shelf_edge_0))                
                txt.close()
                converged = True
                #break
            
            try:
                #   cleanup
                del model.hdsobj, model.ucnobj, model.cbbobj
                os.remove(os.path.join(model.foldername, model.modelname + '.hds'))
                os.remove(os.path.join(model.foldername, 'MT3D001.UCN'))
                os.remove(os.path.join(model.foldername, 'MT3D.CNF'))            
                os.remove(os.path.join(model.foldername, model.modelname + '.cbc'))
                del model

            #   except when the model doesnt converge
            except(IndexError, OSError):
                print('Model has not converged..')
                #   in that case open the LIST find the cells that have the largest head change in the 1st inner iteraion
                with open(os.path.join(topo_sc_dir, modelname + '.list'), 'r') as file:
                    for num, line in enumerate(file, 1):
                        if 'MAXIMUM HEAD CHANGE FOR EACH ITERATION (1 INDICATES THE FIRST INNER ITERATION)' in line:
                            print('found at line:', num)
                pass
                #sys.exit()

    """ ----------------------------------------------------------------------------------------------- """
    """                                 RUN THE SEA LEVEL RISE SCENARIOS                                """
    """ ----------------------------------------------------------------------------------------------- """
    
    equilibrium_names = ['AP_00000_to_00100', 'AP_00100_to_00200', 'AP_00200_to_00300', 'AP_00300_to_00400', 'AP_00400_to_00500']
    slr_names = ['RCP_26_', 'RCP_45_', 'RCP_85_']
    sea_levels_lst = [[0.4, 0.7, 0.8, 0.8, 0.8], [0.5, 1., 1.5, 1.5, 1.5], [0.7, 2.1, 3.7, 3.7, 3.7]]
    sp_duration = [100, 100, 100, 100, 100]
    yrs_len_lst = [100, 100, 100, 100, 100]
    
    perlen_lst = [365.25 * i for i in yrs_len_lst]                     # an array of the stress period lengths
    nstp_lst = [int(i / 50) for i in yrs_len_lst]                                 # number of time steps in each stress period
    nper = 1#len(perlen)
    
    #   loop through all the RCP scenarions (3)
    for a in range(len(slr_names)):
        sea_levels = sea_levels_lst[a]
        slr_name = slr_names[a]
        
        #   loop through the individual sea level stages
        for b in range(len(sea_levels)):
        
            #   get the stress periods sea level, name and duration
            sea_level = sea_levels[b]
            eq_name = slr_name + equilibrium_names[b]
            sp_dur = sp_duration[b]
            print (sea_level, eq_name, sp_dur)
            rch_lst_sp_all = []
            yrs_len = yrs_len_lst[b]
            perlen = [perlen_lst[b]]
            nstp = [nstp_lst[b]]
            nper = 1
        
            #   save the output only in middle and end of the stress periods
            if yrs_len == 1000:
                th_btn_time_step_val = [int((yrs_len / 50) + 1)]
            else:
                th_btn_time_step_val = [int((yrs_len / 10) + 1)]
            
            #   define the model name and check if the model folder already exists, if no create it
            modelname = '%s_%s_%s' % (id_cs_str + '_SRM_' + str_id_subreg, topo_type, del_col_res_str + 'm_' + str(int(del_lay_res)) + 'm') 
            topo_sc_dir = os.path.join(out_dir, model_name_basis, 'sc_' + param_combo_str, eq_name)
            
            print (modelname, topo_sc_dir)
            if not os.path.exists(topo_sc_dir):
                os.makedirs(topo_sc_dir, exist_ok = True)

            #   define the directories of all the dictionaries
            conc_dict_save_dir = os.path.join(topo_sc_dir, '_conc.npy')
            head_dict_save_dir = os.path.join(topo_sc_dir, '_head.npy')
            cbc_dict_save_dir = os.path.join(topo_sc_dir, '_cbc.npy')    
            ghb_check_dict_save_dir = os.path.join(topo_sc_dir, '_ghb_check.npy')                 
            list_dict_save_dir = os.path.join(topo_sc_dir, '_list.npy')           
            nc_folder = os.path.join(topo_sc_dir, '_nc')            
            #   check if the previous equilibrium was reached, based on the converged.txt
            converged_txt_dir = os.path.join(topo_sc_dir, 'converged.txt')
            csv_dir_ts = os.path.join(topo_sc_dir, '_time_steps_DH_DC.csv')   
            csv_dir_abs_cell_change = os.path.join(topo_sc_dir, '_abs_cell_change.csv')   
        
            #   if the text file exists it means that the convergence was reached and we can go to the next equilibrium stage
            #   load the last sconc and strt arrays as starting arrays for now (might be changed later)
            if os.path.exists(converged_txt_dir):
                print("Model   " + modelname + "   already exists, moving to next EQ.")
                #   read in the head and concentration from the previous simulation (the last time step obviously)
                nc_lst = []
                for filename in os.walk(nc_folder):
                    for g in range(len(filename[-1])):
                        if filename[-1][g].endswith(".nc"):
                            print(filename[-1][g].split('_')[-1].split('.nc')[0])
                            nc_lst.append(int(filename[-1][g].split('_')[-1].split('.nc')[0]))
                #   choose the last timestep from the sorted list
                last_ts = sorted(nc_lst)[-1]
                last_nc_name = eq_name + '_' + str(last_ts) + '.nc'
                #   load the nc file and get the head and concentration arrays
                last_nc = xr.open_dataset(os.path.join(os.path.join(out_dir, model_name_basis, 'sc_' + param_combo_str, eq_name), '_nc', last_nc_name))
                sconc_arr = np.expand_dims(np.array(last_nc['solute concentration']), axis = 1)
                strt_arr = np.expand_dims(np.array(last_nc['heads']), axis = 1)
                #rch_lst_sp = last_nc['rch_rate'].values
                sconc_arr[np.isnan(sconc_arr)] = 0.  
                strt_arr[np.isnan(strt_arr)] = 0.            
                #model.sconc_arr = sconc_arr 
                #model.strt_arr = strt_arr
                actual_time = 0     #   model starts at 0 years
                sp = 1
                continue
            #   if the converged.txt doesnt exist, check if there are already the dictionaries created, if yes load the last
            #   time step sconc and strt arrays as starting arrays, also change th sp and actual time values 
            else:
                #if os.path.exists(conc_dict_save_dir) and os.path.exists(head_dict_save_dir):
                if os.path.exists(nc_folder):    
                    #   loop through the folder and find all the files ending with .nc
                    nc_lst = []
                    for filename in os.walk(nc_folder):
                        for g in range(len(filename[-1])):
                            if filename[-1][g].endswith(".nc"):
                                print(filename[-1][g].split('_')[-1].split('.nc')[0])
                                nc_lst.append(int(filename[-1][g].split('_')[-1].split('.nc')[0]))
                    #   choose the last timestep from the sorted list
                    try:
                        last_ts = sorted(nc_lst)[-1]
                        last_nc_name = eq_name + '_' + str(last_ts) + '.nc'
                        #   load the nc file and get the head and concentration arrays
                        last_nc = xr.open_dataset(os.path.join(nc_folder, last_nc_name))
                        sconc_arr = np.expand_dims(np.array(last_nc['solute concentration']), axis = 1)
                        strt_arr = np.expand_dims(np.array(last_nc['heads']), axis = 1)
                        sconc_arr[np.isnan(sconc_arr)] = 0.  
                        strt_arr[np.isnan(strt_arr)] = 0.                    
                        #rch_lst_sp = last_nc['rch_rate'].values                
                        #   assign new values for the sp and the time parameter
                        sp = int(last_ts / yrs_len) + 1
                        actual_time = last_ts
                        
                    #   in case the _nc folder is empty (the simulation could have ended at the first SP without any output to be read)
                    #   so just pass this and reset the sp and time counters while the starting arrays will be the same as the last converged ones..
                    except IndexError:
                        sp = 1
                        actual_time = 0
                        pass
                    
                    time_start_check = False        
                    
                else:
                    actual_time = 0     #   model starts at 0 years
                    sp = 1                                        
                    f = open(csv_dir_ts,'w')
                    f_headers = 'SP, time_start, time_end, max_diff_head, min_diff_head, max_diff_conc, min_diff_conc, fresh_cell_diff_(end_start),\
                                brackish_cell_diff_(end_start), saline_cell_diff_(end_start), frsh_inl_vol,\
                                frsh_inl_pct, frsh_shelf_vol, frsh_shelf_pct, , prev_tot_cells, tot_cells, ghb_cells_check_fail, ghb_cells_check_succ,\
                                ghb_cell_dif_abs_max,  riv_in, rch_in, ghb_in, tot_in, riv_out, drn_out, ghb_out, tot_out, run_time (s)'
                    f.write(" ".join(f_headers.split())) 
                    f.close()     
                    if eq_name != 'BP_30000_to_20000' and sp == 1:
                        sconc_arr[np.isnan(sconc_arr)] = 0.
                        strt_arr[np.isnan(strt_arr)] = 0.
                    pass          

            time_start_check = True         
            time_start = time.time()
            #   start with converged as False..
            converged = False
                    
            while converged is False:# and sp < 2:
                #   create the model object and the top and bottom elevation lists 
                model = cs_model(None, top_elev, None, None, None, None, None, None, None, None, None, None, None, None, None, None,\
                                 None, None, None, False)
    
                model_dict_in = np.load(os.path.join(sc_dir, '_INIT_dict.npy'), allow_pickle = True).item()
                model.ibound_arr = model_dict_in['dis_info']['ibound_arr']                                   
                model.botm = model_dict_in['dis_info']['botm']                                     
                model.top = model_dict_in['dis_info']['top']                           
                model.lay_elev = model_dict_in['dis_info']['lay_elev'] 
                model.top_elev = model_dict_in['dis_info']['top_elev'] 
                model.bot_elev = model_dict_in['dis_info']['bot_elev'] 
                model.hk_arr = model_dict_in['dis_info']['hk_arr']
                model.vk_arr = model_dict_in['dis_info']['vk_arr']
                model.zbot = model_dict_in['dis_info']['zbot'] 
                model.nlay = model_dict_in['dis_info']['nlay'] 
                model.nrow = model_dict_in['dis_info']['nrow'] 
                model.nper = model_dict_in['dis_info']['nper'] 
                model.delr = model_dict_in['dis_info']['delr']
                model.delc = model_dict_in['dis_info']['delc']
                model.laycbd = model_dict_in['dis_info']['laycbd'] 
                model.del_lay = model_dict_in['dis_info']['del_lay_res']  #model_dict_in['dis_info']['del_lay'] 
                model.del_col = model_dict_in['dis_info']['del_col_res']  #model_dict_in['dis_info']['del_col'] 
                model.ncol = model_dict_in['dis_info']['ncol'] 
                model.idx_end = model_dict_in['dis_info']['idx_end'] 
                model.nlay_end_val = model_dict_in['dis_info']['nlay_end_val'] 
                model.nlay_end_idx = model_dict_in['dis_info']['nlay_end_idx'] 
                model.end_idx = model_dict_in['dis_info']['end_idx'] 
                model.x0 = model_dict_in['dis_info']['x0']
                model.x1 = model_dict_in['dis_info']['x1']    
                model.wtd_elev = model_dict_in['dis_info']['wtd_elev'] 
                model.cst_offset = model_dict_in['dis_info']['cst_offset']
                model.lay_idx = model_dict_in['dis_info']['lay_idx']
                model.col_idx = model_dict_in['dis_info']['col_idx'] 
                model.idx_start = model_dict_in['dis_info']['idx_start'] 
                model.x_start = model_dict_in['dis_info']['x_start']
                model.x_end = model_dict_in['dis_info']['x_end']
                model.soil_thk = model_dict_in['dis_info']['soil_thk'] 
                model.soil_type = model_dict_in['dis_info']['soil_type'] 
                model.cst_idx = model_dict_in['dis_info']['cst_idx']
                model.delv = model_dict_in['dis_info']['delv']
                model.k_soil = model_dict_in['dis_info']['k_soil'] 
                model.drn_dens = model_dict_in['dis_info']['drn_dens'] 
                model.glhymps_top_lay = model_dict_in['dis_info']['glhymps_top_lay']
                model.glhymps_bot_lay = model_dict_in['dis_info']['glhymps_bot_lay'] 
                model.coscat_id = model_dict_in['dis_info']['coscat_id']
                model.del_col_res = model_dict_in['dis_info']['del_col_res'] 
                model.del_lay_res = model_dict_in['dis_info']['del_lay_res'] 
                model.perlen = perlen
                model.nstp = nstp
                model.laytyp = laytyp_val
                x_lst_res = model_dict_in['dis_info']['x_lst_res'] 
                y_lst_res = model_dict_in['dis_info']['y_lst_res'] 

                #   make sure to change the elevation DEM type, but if we do it like below it also maintains the same extent of the model domain
                if topo_100_type != 'gebco':
                    topo_lst_sc = [round(i, 2) for i in topo_nc['gebco_' + topo_100_type + '_avg'].values.tolist()]
                else:
                    topo_lst_sc = [round(i, 2) for i in topo_nc['gebco_avg'].values.tolist()]
                #   calculate new index 
                mid_idx = int(2000 + model.cst_offset * 10)
                topo_lst_new = topo_lst_sc[mid_idx - model.cst_idx : mid_idx - model.cst_idx + len(model.top_elev)]
                
                #   now connect the two lists into one
                topo_final = topo_lst_new[: model.cst_idx] + model.top_elev[model.cst_idx :]
                model.top_elev = topo_final
    
                #   check the HK array for cells with largely different HK values - non convergence occurs when high and low permeable cells are too close to each other
                #   1) all cells with HK larger than 30 m/d will be limited to 30 m/d
                model.hk_arr[model.hk_arr > 30.] = 30.
                model.vk_arr = model.hk_arr * 0.1
                #   only if it is the start of the first stress period load the pre-defined sconc and strt
                model.sconc_arr = sconc_arr
                model.strt_arr = strt_arr

                #   all the cells that represent the sea floor will be set to saline before starting simulation of each SP
                cst_col = len(model.top_elev)
                for f in reversed(model.top_elev):
                    if f < sea_level:
                        cst_col -= 1
                    else:
                        break
                for x in range(cst_col - 1, model.sconc_arr.shape[-1]):
                    ibound_act_lay_idxs = [k for k, x in enumerate(model.ibound_arr[:, 0, x].tolist()) if x == 1]
                    #   select the top elevation to see if its lower than sea level, if so turn the cell to saline
                    if len(ibound_act_lay_idxs) > 0:
                        if model.top_elev[x] < sea_level:
                            model.sconc_arr[ibound_act_lay_idxs[0], 0, x] = 35.0
                    
                if eq_name == 'BP_30000_to_20000' and sp == 1:
                    model.strt_arr[:,:,cst_col:] = sea_level
                       
                #   find the shelf edge
                if init_cond != 'fresh':
                    cont_shelf_edge = mgt.find_shelf_break(model, model.top_elev, True) #   try to find the shelf break     
                else:
                    cont_shelf_edge = mgt_frsh.find_shelf_break_fresh(model, model.top_elev, True) #   try to find the shelf break     
                
                #   make sure the position of the coastline is correct - in the INIT dictionary it is 
                if eq_name == 'BP_30000_to_20000':
                    idx_cst_0 = int(abs(round(model.x_start + model.cst_offset, 2)) * 10)
                    idx_shelf_edge_0 = int((cont_shelf_edge[0] - model.x_start + model.cst_offset) * 10)                
                    if idx_shelf_edge_0 > model.ibound_arr.shape[-1]:
                        idx_shelf_edge_0 = model.ibound_arr.shape[-1] - 1
                else:
                    #   open the text file to find the position of the coast and cont shelf at the EQ 0  
                    txt = open(os.path.join(os.path.join(out_dir, model_name_basis, 'sc_' + param_combo_str, 'BP_30000_to_20000'), 'converged.txt'),'r')
                    lines = txt.readlines()
                    idx_cst_0 = int(lines[0].strip('\n'))
                    idx_shelf_edge_0 = int(lines[1])
                    if idx_shelf_edge_0 > model.ibound_arr.shape[-1]:
                        idx_shelf_edge_0 = model.ibound_arr.shape[-1] - 1
                        

                #   specify input for the GHB package
                cond_factor = 1.
                clay_incr = True # means that there will not be any highering up of conductance if the calculated value is too low
                                 # this usually happened in the clay layers where the conductivity is low and therefore the conductance too
                cond_limit = 1.            
                clay_cond = 1.
                #   check if the last column offshore has only one active cell, if it is more then also assign GHB to those cells
                last_col_act_cells = [i for i, x in enumerate(model.ibound_arr[:, 0, -1].tolist()) if x == 1]   
                if len(last_col_act_cells) > 1:
                    last_col_ghb = True
                else:
                    last_col_ghb = False
                
                #   set the time start check
                if actual_time != 0:
                    time_start_check = False
    
                x_lst_100 = np.arange(model.x_start, model.x_end, 100. / 1000.)
                if sp ==1:
                    if rch_type == 'paleo':
                        gw_rch_fig_out = os.path.join(topo_sc_dir, '_rch_p_min_et_input.png')
                        gw_rch_lst_100 = mfc.generate_ln_rand_lst(param_combo, gw_rch_mu[-1], gw_rch_std[-1], len(x_lst_100),\
                                                           'GW RCH input values per column', 'RCH rate in mm/d', gw_rch_fig_out, plot = True) 
                    else:
                        gw_rch_fig_out = os.path.join(topo_sc_dir, '_rch_p_min_et_input.png')
                        gw_rch_lst_100 = mfc.generate_ln_rand_lst(param_combo, gw_rch_mu[-1], gw_rch_std[-1], len(x_lst_100),\
                                                           'GW RCH input values per column', 'RCH rate in mm/d', gw_rch_fig_out, plot = True)     
                    
                if 'gw_rch_lst_100' not in locals():
                    gw_rch_fig_out = os.path.join(topo_sc_dir, '_rch_p_min_et_input.png')
                    gw_rch_lst_100 = mfc.generate_ln_rand_lst(param_combo, gw_rch_mu[-1], gw_rch_std[-1], len(x_lst_100),\
                                                       'GW RCH input values per column', 'RCH rate in mm/d', gw_rch_fig_out, plot = True)                     
                
                #   create the recharge input..without the recharge adjusting factor 
                gw_rch_lst_m_d = []
                if del_col != del_col_res:
                    gw_rch_lst = np.interp(x_lst_res, x_lst_100, gw_rch_lst_100).tolist()
                else:
                    gw_rch_lst = gw_rch_lst_100
                gw_rch_lst_m_d = [round(i * 0.001 / 365., 6) for i in gw_rch_lst]
                gw_rch_lst_m_d = gw_rch_lst_m_d[:model.ibound_arr.shape[-1]]                

                #   create input for all but the DIS and BAS packages
                #model.ghb_input(sea_level, 1., clay_cond, cond_limit, clay_incr, last_col_ghb, True, True, True, wtd_inl_lvl = True)  
                
                #   create the GHB package input here
                itype = flopy.mt3d.Mt3dSsm.itype_dict()
                model.ghb_arr = model.ibound_arr * 1
                model.ghb_input_lst = []
                model.ssmdata = []
                #       1) the inland column and the offshore last column where the CONDUCTANCE is calculated as Hk * del_lay * LD 
                for d in range(model.ibound_arr.shape[0]):
                    if model.ibound_arr[d, 0, 0] == 1:
                        model.ghb_input_lst.append([d, 0, 0, model.top_elev[0], model.hk_arr[d, 0, 0] * model.del_lay_res * ld_rate])
                        model.ghb_arr[d, 0, 0] = -1
                        model.ssmdata.append([d, 0, 0, 0.0, itype['GHB']])            
                
                for c in range(model.ibound_arr.shape[0]):
                    if model.ibound_arr[c, 0, model.ibound_arr.shape[-1] - 1] == 1:
                        model.ghb_input_lst.append([c, 0, model.ibound_arr.shape[-1] - 1, max(sea_level, model.top_elev[-1]), min(100., model.hk_arr[c, 0, model.ibound_arr.shape[-1] - 1] * model.del_lay_res * ld_rate)])
                        model.ghb_arr[c, 0, model.ibound_arr.shape[-1] - 1] = -1
                        model.ssmdata.append([c, 0, model.ibound_arr.shape[-1] - 1, 35.0, itype['GHB']])      

                #       2) the offshore columns that represent the sea floor, where the CONDUCTANCE is calculated as Hk * del_col * LD          
                for f in range(cst_col, len(model.top_elev)):
                    #   get the index of the first non-zero ibound cell
                    try:
                        top_cell_lay_idx = model.ibound_arr[:, 0, f].tolist().index(1)
                        topo_val = model.top_elev[f]
                        if topo_val < sea_level:
                            #model.ghb_input_lst.append([top_cell_lay_idx, 0, f, sea_level, model.vk_arr[top_cell_lay_idx, 0, f] * model.del_col_res * ld_rate])   #   sea level is 0.0m (present condition)
                            model.ghb_input_lst.append([top_cell_lay_idx, 0, f, sea_level, min(100., model.vk_arr[top_cell_lay_idx, 0, f] * model.del_col_res * ld_rate)])
                            model.ghb_arr[top_cell_lay_idx, 0, f] = -1
                            model.ssmdata.append([top_cell_lay_idx, 0, f, 35.0, itype['GHB']])
                    except ValueError:
                        pass            
                
                #   write the final output dictionary, inlcude each stress period
                model.ghb_arr_in = {}
                for d in range(len(model.perlen)):
                    model.ghb_arr_in[d] = model.ghb_input_lst                     
                model.rch_input_from_lst_vk_limit(gw_rch_lst_m_d, sea_level)
                
                #model.drn_input(sea_level, 0.0)
                #model.drn_input_wtd(sea_level)
                
                #  drainage is assigned only to cells that receive recharge - cells with elev above sea level
                drn_input_lst = []
                #for i in range(self.ncol):
                #for i in range(self.ibound_arr.shape[-1]):
                for i in range(len(model.top_elev)):
                    ibound_col_lst = model.ibound_arr[:, 0, i].tolist()
                    #   check the 1st column with ibound_val = 1 (active cell)
                    try:
                        drn_lay = ibound_col_lst.index(1)
                        ##  now check if the elevation is below sea level, if so assign the cell to ghb list
                        if model.top_elev[i] >= sea_level:
                            cond_cell = ((model.hk_arr[drn_lay, 0, i] / 10.) * model.del_col)# / model.del_lay
                            #drn_input_lst.append([drn_lay, 0, i, model.top_elev[i], model.del_col_res * ld_rate / model.del_lay_res])
                            drn_input_lst.append([drn_lay, 0, i, model.top_elev[i], cond_cell])
                            #drn_input_lst.append([drn_lay, 0, i, model.top_elev[i], model.vk_arr[drn_lay, 0, i] / model.del_lay_res])
                        else:
                            pass
                    except ValueError:
                        pass
                #   write the final output dictionary, inlcude each stress period
                model.drn_arr_in = {}
                for c in range(len(model.perlen)):
                    model.drn_arr_in[c] = drn_input_lst            
                    
                model.pcg_input(hclose, rclose)
                model.oc_input(1)#(th_oc_time_step_val)
                model.btn_input(porosity, dt0, th_btn_time_step_val, ifmtcn, nprs, chkmas, nprobs, nprmas)
                model.adv_input(mixelm)
                model.dsp_input(al, trpt, trpv, dmcoef)
                model.gcg_input(iter1, mxiter, isolve, cclose)
                model.vdf_input(iwtable, densemin, densemax, denseref, denseslp, firstdt) 
                model.ssm_input()
                
                #   write all the input files
                model.write_all_input(modelname, mf_exe_dir, mt3d_exe_dir, swat_exe_dir, topo_sc_dir)
                #   run the model
                model.run_model()
                      
                #   read the model output and create output pictures and video
                model.read_output()
            
                #   check if the plotting directories exist, if yes then delete them
                dest_dir_conc = os.path.join(topo_sc_dir, '_conc_pngs')    
                if not os.path.exists(dest_dir_conc):
                    os.makedirs(dest_dir_conc, exist_ok = True)        
                dest_dir_conc_flow = os.path.join(topo_sc_dir, '_conc_flow_pngs')    
                if not os.path.exists(dest_dir_conc_flow):
                    os.makedirs(dest_dir_conc_flow, exist_ok = True)                                 
                dest_dir_head_fresh = os.path.join(topo_sc_dir, '_head_fresh_pngs')    
                if not os.path.exists(dest_dir_head_fresh):
                    os.makedirs(dest_dir_head_fresh, exist_ok = True)                                   
                nc_dir = os.path.join(topo_sc_dir, '_nc')    
                if not os.path.exists(nc_dir):
                    os.makedirs(nc_dir)         
                    
                #   conc/head profiles to be read every stepth year and define the convergence criteria
                if yrs_len == 1000:
                    step = 500.      
                elif yrs_len == 100:
                    step = 50
                    
                conc_limit = 0.05       
                head_limit = 0.05                    
                #   read in the number of time steps and transform the time from days to years
                conc_ts_lst, heads_ts_lst = model.time_steps, model.times_heads
                conc_ts_lst_yrs = [round(k / 365.25, 0) for k in conc_ts_lst]
                heads_ts_lst_yrs = [round(k / 365.25, 0) for k in heads_ts_lst]
    
                x_coord_lst = np.linspace(model.x_start + (model.del_col / 1000.) / 2, model.x_end - (model.del_col / 1000.) / 2, model.ncol)
                x_coord_lst = [round(k, 3) for k in x_coord_lst]
                y_coord_lst = np.linspace(model.top - 5., model.zbot + 5., model.nlay)
                y_coord_lst = [round(k, 3) for k in y_coord_lst]
                cst_idx_0 = int(abs(model.x_start) * 1000. / model.del_col_res)

                #   for each stepth year print a concentration, head profile, if it is the 1st SP then start with the 
                #   profiles at TS = 0yrs..
                if sp == 1:
                    cbc_dict = dict()
                    ghb_check_dict = dict()
                    list_dict = dict()
                    conc_arr_prev = model.sconc_arr.astype(dtype = np.float64)
                    head_arr_prev = model.strt_arr.astype(dtype = np.float64)        
                    zero_arr = model.ibound_arr.astype(dtype = np.float64) * 0.0
                    cbc_dict[0] = {'q_ghb': zero_arr, 'q_drn': zero_arr, 'q_rch': zero_arr, 'q_riv' : zero_arr}                       
                    ghb_check_dict[0] = {'ghb_check': zero_arr, 'ghb_head_diff_lst': [-1], 'ghb_head_diff_max': [-1], 'succ_cnt': -1, 'fail_cnt': -1}                      
                    head_fresh_arr_prev = mfc.pt_hd_to_frsh_hd(head_arr_prev, conc_arr_prev, denseref, denseslp, model.top, model.botm)
                    #   for the concentration, heads and cbc create a netcdf file (to save memory)
                    xa_sum = xr.Dataset(data_vars = {'solute concentration' : (('y', 'x'), conc_arr_prev[:, 0, :]),
                                                     'heads' : (('y', 'x'),head_arr_prev[:, 0, :]),
                                                     'cbc Q_right' : (('y', 'x'), zero_arr[:, 0, :]),
                                                     'cbc Q_bottom' : (('y', 'x'), zero_arr[:, 0, :])},
                                        coords = {'x' : x_coord_lst,
                                                  'y' : y_coord_lst})
                    xa_sum = xa_sum.assign_coords(time = 0)
                    xa_name = eq_name + '_0.nc'
                    xa_sum.to_netcdf(os.path.join(nc_dir, xa_name))                
                else:
                    cbc_dict = np.load(cbc_dict_save_dir, allow_pickle = True).item()  
                    ghb_check_dict = np.load(ghb_check_dict_save_dir, allow_pickle = True).item()  
                    list_dict = np.load(list_dict_save_dir, allow_pickle = True).item()
                    #   the recharge adjusting factor in case there are huge head differences in the model = injecting recharge into low permeable systems
                    #   if the drainage is larger than recharge then no need to decrease recharge anymore
                    #   second condition is that all the heads have to be lower than topography + 0.1m 
                    rch_adj_factor_lst = []
                    for y in range(model.ncol):
                        if max(strt_arr[:,0,y]) > model.top_elev[y] + 1.0:
                            rch_adj_factor_lst.append(round(1 - (list_dict[list(list_dict.keys())[-1]]['DRAINS_OUT'] / list_dict[list(list_dict.keys())[-1]]['RECHARGE_IN']), 3))
                        else:
                            rch_adj_factor_lst.append(1.)

                if yrs_len == 1000:
                    ts_to_plt = conc_ts_lst_yrs[::10]
                    #   create a list of time steps to create plots of differences in concentrations and heads 
                    ts_diff_to_plt = conc_ts_lst_yrs[::10]
                else:
                    ts_to_plt = conc_ts_lst_yrs[::5]
                    #   create a list of time steps to create plots of differences in concentrations and heads 
                    ts_diff_to_plt = conc_ts_lst_yrs[::5]
                    
                    #   define a counter for counting the number of time steps where the fresh, brackish, saline are constant
                cnt_volumes_const = 0
                cell_diff_pos_lst, cell_diff_neg_lst = [], []
            
                #   loop through the list and plot the heads and concentration profiles for each of the chosen time steps
                for z in range(1, len(ts_to_plt)):
                    #   select the index of the time step to be plotted and get the concentration and heads arrays
                    ts_plt_conc = min(conc_ts_lst_yrs, key = lambda x:abs(x - ts_to_plt[z]))
                    ts_plt_head = min(heads_ts_lst_yrs, key = lambda x:abs(x - ts_to_plt[z]))
                    conc_arr = model.ucnobj.get_data(totim = model.time_steps[conc_ts_lst_yrs.index(ts_plt_conc)]).astype(dtype = np.float64)
                    head_arr = model.hdsobj.get_data(totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)]).astype(dtype = np.float64)
                    
                    print('analyzing results for time = ', str(ts_plt_conc))
                    
                    #   calculate the freshwater head array
                    head_fresh_arr = mfc.pt_hd_to_frsh_hd(head_arr, conc_arr, denseref, denseslp, model.top, model.botm)                    
                    qx_in = model.cbbobj.get_data(text='flow right face', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]
                    qz_in = model.cbbobj.get_data(text = 'flow lower face', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]                          
                    #   create plots for given years
                    real_time_yrs = int(actual_time + ts_plt_conc)       
                    #   Plotting of 2D concentration profiles at given time
                    #   find the location of the coastline for the zoom and plotting the coastline line in the plots
                    for y in range(0, len(model.top_elev)):
                        if model.top_elev[y] > sea_level:
                            shift_cst = y / (1000. / del_col_res) + model.x_start
                            zoom_min = round(shift_cst - 10., 2)
                            zoom_max = round(shift_cst + 10., 2)
                            pass
                        else:
                            shift_cst = y / (1000. / del_col_res) + model.x_start 
                            zoom_min = round(shift_cst - 10., 2)
                            zoom_max = round(shift_cst + 10., 2)
                            break
                    
                    #   plot the arrays
                    model.plot_head_profile_yrs(dest_dir_head_fresh, x_lst_res, head_fresh_arr, sea_level, real_time_yrs, shift_cst, None, None, zoom = [zoom_min, zoom_max], head_diff = False)
                    model.plot_conc_profile_yrs(dest_dir_conc, x_lst_res, conc_arr, sea_level, real_time_yrs, shift_cst, zoom = [zoom_min, zoom_max], conc_diff = False)                    
                                    
                    #   get the arrays for the difference between the two consecutive time steps, if the z = 0 there is obviously no index before that so read
                    #   in either the starting arrays (SP 0) or the last conc and head arrays from the previous SP
                    ts_plt_conc_prev = min(conc_ts_lst_yrs, key = lambda x:abs(x - ts_to_plt[z - 1]))
                    ts_plt_head_prev = min(heads_ts_lst_yrs, key = lambda x:abs(x - ts_to_plt[z - 1]))
                    #   next calculate the difference arrays and also the number of cell changes between these two steps
                    conc_arr_prev = model.ucnobj.get_data(totim = model.time_steps[conc_ts_lst_yrs.index(ts_plt_conc_prev)]).astype(dtype = np.float64)
                    conc_arr_diff = conc_arr - conc_arr_prev
                    conc_max_diff, conc_min_diff = np.nanmax(conc_arr_diff), np.nanmin(conc_arr_diff)   
                    head_arr_prev = model.hdsobj.get_data(totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head_prev)]).astype(dtype = np.float64)
                    head_fresh_arr_prev = mfc.pt_hd_to_frsh_hd(head_arr_prev, conc_arr_prev, denseref, denseslp, model.top, model.botm)      
                    head_arr_diff = head_arr - head_arr_prev
                    head_max_diff, head_min_diff = np.nanmax(head_arr_diff), np.nanmin(head_arr_diff)   
                    head_fresh_arr_diff = head_fresh_arr - head_fresh_arr_prev
    
                    #   read the other flow budget information 
                    if len(drn_input_lst) != 0:                       
                        qdrains = model.cbbobj.get_data(text = '          DRAINS', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]
                        drn_q = model.cbbobj.get_data(text = '          DRAINS', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])                                 
                    qghb = model.cbbobj.get_data(text = ' HEAD DEP BOUNDS', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]
                    qrch = model.cbbobj.get_data(text = '        RECHARGE', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])[0]    
                    rch_q = model.cbbobj.get_data(text = '        RECHARGE', totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)])
                    
                    
                    #   plot the HK array, just as a double check that it stays the same throughout the simulation
                    if sp == 1 and z == 1:
                        model.plot_HK_arr_res(topo_sc_dir, x_lst_res, model.hk_arr, sea_level, real_time_yrs, shift_cst, zoom = [-10., 10.])
                    
                    #   plot the respective arrays and overwrite the ts_plt_conc_prev and ts_plt_head_prev
                    #   in non active cells assign the NaN value
                    for x in range(conc_arr_diff.shape[0]):
                        for j in range(conc_arr_diff.shape[-1]):
                            if model.ibound_arr[x, 0, j] != 1:
                                conc_arr_diff[x, 0, j] = 999.999
                                head_arr_diff[x, 0, j] = 999.999
                                head_fresh_arr_diff[x, 0, j] = 999.999
                                
                    #   get the number of fresh/brackish/saline cells in the model domain at TS_1 and TS_2 and check the change between them
                    conc_ts_prev_fresh = ((-100.0 <= conc_arr_prev) & (conc_arr_prev <= 0.5)).sum()
                    conc_ts_prev_brackish = ((0.5 < conc_arr_prev) & (conc_arr_prev <= 20.)).sum()
                    conc_ts_prev_saline = ((conc_arr_prev > 20.0) & (conc_arr_prev < 100.)).sum()
                    conc_ts_fresh = ((-100.0 <= conc_arr) & (conc_arr <= 0.5)).sum()
                    conc_ts_brackish = ((0.5 < conc_arr) & (conc_arr <= 20.)).sum()
                    conc_ts_saline = ((conc_arr > 20.0) & (conc_arr < 100.)).sum()
                    
                    prev_tot_cells = conc_ts_prev_fresh + conc_ts_prev_brackish + conc_ts_prev_saline
                    tot_cells = conc_ts_fresh + conc_ts_brackish + conc_ts_saline
                    
                    #   change the values in arrays according to ranges to mark fresh, brackish and saline cells and see where the differences occur
                    arr_curr, arr_prev = conc_arr * 1.0, conc_arr_prev * 1.0
                    arr_curr[(-100.0 <= arr_curr) & (arr_curr <= 0.5)] = -1
                    arr_curr[(0.5 < arr_curr) & (arr_curr <= 20.)] = 0
                    arr_curr[(20. < arr_curr) & (arr_curr <= 100.)] = 1
                    arr_prev[(-100.0 <= arr_prev) & (arr_prev <= 0.5)] = -1
                    arr_prev[(0.5 < arr_prev) & (arr_prev <= 20.)] = 0
                    arr_prev[(20. < arr_prev) & (arr_prev <= 100.)] = 1                    
                    #   create the difference array and count the number of 1 or -1 cells
                    arr_diff = arr_curr - arr_prev
                    unique, counts = np.unique(arr_diff, return_counts = True)
                    #   append to the list of cells that are different
                    if 1 in unique:
                        pos_lst = list(np.where(arr_diff == 1.))
                        for l in range(len(pos_lst[0])):
                            cell_diff_pos_lst.append([pos_lst[0][l], pos_lst[1][l], pos_lst[2][l]])
                    if -1 in unique:
                        neg_lst = np.where(arr_diff == -1.)
                        for l in range(len(neg_lst[0])):
                            cell_diff_neg_lst.append([neg_lst[0][l], neg_lst[1][l], neg_lst[2][l]])                    
                    del arr_curr, arr_prev
                    
                    #   get the differences in number of cells of each category
                    fresh_change, brack_change, saline_change = conc_ts_fresh - conc_ts_prev_fresh, conc_ts_brackish - conc_ts_prev_brackish, conc_ts_saline - conc_ts_prev_saline
                    #   check if all the changes are 0, if yest increase the counter
                    if fresh_change == 0 and brack_change == 0 and saline_change == 0:
                        cnt_volumes_const += 1
                    #   append the results of the EQ to the final result list 
                    if idx_shelf_edge_0 > model.ibound_arr.shape[-1]:
                        idx_shelf_edge_0 = model.ibound_arr.shape[-1] - 1
                        
                    #   find the coastal index - taking current sea level as position of the coastline
                    idx_shelf_edge_0 = int((cont_shelf_edge[0] - model.x_start + model.cst_offset) * (1000. / del_col_res))
                    #volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx_0, idx_shelf_edge_0)
                    volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx, idx_shelf_edge_0)
                    
                    #   write the results into the csv
                    f_res = open(res_summary_csv_dir,'a')
                    f_res.write("\n")
                    f_res.write(eq_name + ',' + str(real_time_yrs) + ',' + str(volumes[1]) + ',' + str(round(100. * (volumes[1] / volumes[0]), 1))  + ',' +\
                                str(volumes[3]) + ',' +  str(round(100. * (volumes[3] / volumes[2]), 1))) 
                    f_res.close()   
                    list_dir = os.path.join(topo_sc_dir, modelname + '.list')
                    swt_list = flopy.utils.mflistfile.SwtListBudget(list_dir)
                    budget = swt_list.get_budget()
                    #   check the GHB condition
                    check_ghb = model.check_ghb_head(heads_ts_lst_yrs.index(ts_plt_head), 0.0001)
                    #   write the results into the output CSV file
                    f = open(csv_dir_ts,'a')
                    f.write("\n")
                    f.write(str(sp) + ',' + str(real_time_yrs - step) + ',' + str(real_time_yrs) + ',' + str(head_max_diff) + ',' + str(head_min_diff)\
                            + ',' + str(conc_max_diff) + ',' + str(conc_min_diff) + ',' + str(conc_ts_fresh - conc_ts_prev_fresh) + ',' +\
                            str(conc_ts_brackish - conc_ts_prev_brackish) + ',' + str(conc_ts_saline - conc_ts_prev_saline) + ',' +\
                            str(volumes[1]) + ',' + str(round(100. * (volumes[1] / volumes[0]), 1))  + ',' +\
                            str(volumes[3]) + ',' +  str(round(100. * (volumes[3] / volumes[2]), 1)) + ',' +\
                            str(prev_tot_cells) + ',' +  str(tot_cells) + ',' + str(model.ghb_fail) + ',' + str(model.ghb_succ) + ',' +\
                            str(round(model.ghb_head_check_max_diff, 5)) + ',' + str(budget[0][0][7])  + ',' + str(budget[0][0][6])  + ',' +\
                            str(budget[0][0][9]) + ',' + str(budget[0][0][12]) + ',' + str(budget[0][0][13]) + ',' + str(budget[0][0][16])  + ',' +\
                            str(model.run_time))
                    f.close() 
    
                    #if z % 5 == 0:
                        
                    #   for the concentration, heads and cbc create a netcdf file (to save memory)
                    xa_sum = xr.Dataset(data_vars = {'solute concentration' : (('y', 'x'), conc_arr[:, 0, :]),
                                                     'heads' : (('y', 'x'),head_arr[:, 0, :]),
                                                     'cbc Q_right' : (('y', 'x'), qx_in[:, 0, :]),
                                                     'cbc Q_bottom' : (('y', 'x'), qz_in[:, 0, :]),
                                                     'rch_rate' : (('x'), gw_rch_lst_m_d[:])},
                                        coords = {'x' : x_coord_lst,
                                                  'y' : y_coord_lst})
                    #xa_sum.expand_dims({'time':real_time_yrs})
                    xa_sum = xa_sum.assign_coords(time = real_time_yrs)
                    xa_name = eq_name + '_' + str(real_time_yrs) + '.nc'
                    xa_sum.to_netcdf(os.path.join(nc_dir, xa_name))                         
                    
                    #   append the conc and head arrays at certain times to the output dictionaries
                    if len(drn_input_lst) != 0:   
                        ghb_check_dict[real_time_yrs] = {'ghb_check': model.ghb_head_check, 'ghb_head_diff_lst': model.ghb_head_diff_lst,\
                                                         'ghb_head_diff_max': model.ghb_head_check_max_diff, 'succ_cnt': model.ghb_succ, 'fail_cnt': model.ghb_fail}                    
                        cbc_dict[real_time_yrs] = {'q_ghb': qghb, 'q_drn': qdrains, 'q_rch': qrch}                
                        #   store the flow budget from the list file
                        list_dict[real_time_yrs] = {'STORAGE_IN': budget[0][0][3], 'CONSTANT_HEAD_IN': budget[0][0][4], 'DRAINS_IN': budget[0][0][5],\
                                                    'HEAD_DEP_BOUNDS_IN': budget[0][0][6], 'RECHARGE_IN': budget[0][0][7], 'DCDT_IN': budget[0][0][8],\
                                                    'TOTAL_IN': budget[0][0][9], 'STORAGE_OUT': budget[0][0][10], 'CONSTANT_HEAD_OUT': budget[0][0][11],\
                                                    'DRAINS_OUT': budget[0][0][12], 'HEAD_DEP_BOUNDS_OUT': budget[0][0][13], 'RECHARGE_OUT': budget[0][0][14],\
                                                    'DCDT_OUT': budget[0][0][15], 'TOTAL_OUT': budget[0][0][16], 'IN-OUT': budget[0][0][17],\
                                                    'PERCENT_DISCREPANCY': budget[0][0][18]}      
                    else:
                        ghb_check_dict[real_time_yrs] = {'ghb_check': model.ghb_head_check, 'ghb_head_diff_lst': model.ghb_head_diff_lst,\
                                                         'ghb_head_diff_max': model.ghb_head_check_max_diff, 'succ_cnt': model.ghb_succ, 'fail_cnt': model.ghb_fail}                    
                        cbc_dict[real_time_yrs] = {'q_ghb': qghb, 'q_drn': -1, 'q_rch': qrch}                
                        #   store the flow budget from the list file
                        list_dict[real_time_yrs] = {'STORAGE_IN': budget[0][0][3], 'CONSTANT_HEAD_IN': budget[0][0][4], 'DRAINS_IN': -1,\
                                                    'HEAD_DEP_BOUNDS_IN': budget[0][0][5], 'RECHARGE_IN': budget[0][0][6], 'DCDT_IN': budget[0][0][7],\
                                                    'TOTAL_IN': budget[0][0][8], 'STORAGE_OUT': budget[0][0][9], 'CONSTANT_HEAD_OUT': budget[0][0][10],\
                                                    'DRAINS_OUT': -1, 'HEAD_DEP_BOUNDS_OUT': budget[0][0][11], 'RECHARGE_OUT': budget[0][0][12],\
                                                    'DCDT_OUT': budget[0][0][13], 'TOTAL_OUT': budget[0][0][14], 'IN-OUT': budget[0][0][15],\
                                                    'PERCENT_DISCREPANCY': budget[0][0][16]}                            
                    """
                    rch_adj_factor_lst = []
                    for y in range(model.ncol):
                        if max(head_arr[:,0,y]) > model.top_elev[y] + 1.0:
                            rch_adj_factor_lst.append(round(1 - (list_dict[real_time_yrs]['DRAINS_OUT'] / list_dict[real_time_yrs]['RECHARGE_IN']), 3))
                        else:
                            rch_adj_factor_lst.append(1.)                        
                    rch_adj_factor = round(1 - (list_dict[real_time_yrs]['DRAINS_OUT'] / list_dict[real_time_yrs]['RECHARGE_IN']), 3)
                    """
                    
                    model.plot_conc_flow_vectors_profile_yrs(dest_dir_conc_flow, conc_arr, sea_level, real_time_yrs, shift_cst, qx_in, qz_in, zoom = [zoom_min, zoom_max], conc_diff = False)
    
                #   assign the prev arrays to be the last arrays of this SP and add the duration of the SP to the actual time counter
                conc_arr_prev, head_arr_prev = conc_arr, head_arr
                actual_time += (perlen[0] / 365.25)
    
                #   save the dictionary with the ibound arrays
                #   put the dictionaries in the right order (time wise)
                cbc_in_ordered = collections.OrderedDict(sorted(cbc_dict.items()))     
                ghb_check_in_ordered = collections.OrderedDict(sorted(ghb_check_dict.items()))          
                list_in_ordered = collections.OrderedDict(sorted(list_dict.items()))   
                np.save(cbc_dict_save_dir, cbc_in_ordered)        
                np.save(ghb_check_dict_save_dir, ghb_check_in_ordered)
                np.save(list_dict_save_dir, list_in_ordered)
    
                #   check if the model converged, convergence is defined as when the difference between concentrations between two SPs is less than one cell
                if cnt_volumes_const < int(perlen[0] / (step * 365.25)) - 1:
                    #   now check the list and if it is only one cell changing back and forth 
                    if cell_diff_pos_lst != [] or cell_diff_neg_lst != []:
                        
                        #   check which list is largest and then loop through the smaller one
                        if len(cell_diff_pos_lst) >= len(cell_diff_neg_lst):
                            lst_to_loop = [list(x) for x in set(tuple(x) for x in cell_diff_pos_lst)]
                            lst_no_loop = [list(x) for x in set(tuple(x) for x in cell_diff_neg_lst)]
                        else:
                            lst_to_loop = [list(x) for x in set(tuple(x) for x in cell_diff_neg_lst)]
                            lst_no_loop = [list(x) for x in set(tuple(x) for x in cell_diff_pos_lst)]
                        
                        #   loop through the list and for each cell check if it occurs in the other list, first set up counters to count how many cells are 
                        #   really changing and not just oscillating
                        non_loop_cells = len(lst_no_loop)
                        loop_cells = len(lst_to_loop)
                        #   loop through the list and if it is also in the smaller list decrease the total counter
                        for g in range(len(lst_to_loop)):
                            if lst_to_loop[g] in lst_no_loop:
                                non_loop_cells -= 1
                                loop_cells -= 1
                        #   if the counter itself is lower or equal to 1 it is considered a stable situation (change of 1 cell per 1000 years)
                        if loop_cells + non_loop_cells <= 1:
                            converged = True
                        else:
                            converged = False
                    else:
                        loop_cells, non_loop_cells = -1, -1
                        converged = True
                else:
                    converged = True
    
                #   then check if the time specific for the eq already expired, if yes then change the converged to True
                if not converged:
                    #   only in case the time is reached
                    if sp_dur == -1:
                        pass
                    elif actual_time < sp_dur:
                        pass
                    else:
                        converged = True
    
                #   check if the model converged or not
                if not converged:
                    sconc_arr = model.ucnobj.get_data(totim = model.time_steps[conc_ts_lst_yrs.index(ts_plt_conc)]).astype(dtype = np.float64)
                    strt_arr = model.hdsobj.get_data(totim = model.times_heads[heads_ts_lst_yrs.index(ts_plt_head)]).astype(dtype = np.float64)
                                    
                    sp += 1  
                    gc.collect()
                    #break
                    
                else:
                    #   save the dictionary with the ibound arrays
                    cbc_dict_save_dir = os.path.join(topo_sc_dir, '_cbc.npy')    
                    ghb_check_dict_save_dir = os.path.join(topo_sc_dir, '_ghb_check.npy')                 
                    list_dict_save_dir = os.path.join(topo_sc_dir, '_list.npy')   
                    #   put the dictionaries in the right order (time wise)
                    cbc_in_ordered = collections.OrderedDict(sorted(cbc_dict.items()))     
                    ghb_check_in_ordered = collections.OrderedDict(sorted(ghb_check_dict.items()))          
                    list_in_ordered = collections.OrderedDict(sorted(list_dict.items()))   
    
                    strt_arr = model.hdsobj.get_data(totim = model.times_heads[-1])
                    sconc_arr = model.ucnobj.get_data(totim = model.time_steps[-1])                                    
                    
                    np.save(cbc_dict_save_dir, cbc_in_ordered)        
                    np.save(ghb_check_dict_save_dir, ghb_check_in_ordered)
                    np.save(list_dict_save_dir, list_in_ordered)
    
                    del cbc_dict, ghb_check_dict, cbc_in_ordered, ghb_check_in_ordered,\
                    list_dict_save_dir, list_in_ordered
                    
                    time_end = time.time()           
                    run_time = time_end - time_start                   
    
                    #volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx_0, idx_shelf_edge_0)
                    #volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx_0, idx_shelf_edge_0)
                    volumes = mfc.get_frsh_vol_pct(model.ibound_arr, conc_arr, cst_idx, idx_shelf_edge_0)
                    
                    #   create converged text file
                    txt = open(os.path.join(topo_sc_dir, 'converged.txt'),'w')
                    txt.write(str(idx_cst_0) + '\n')
                    txt.write(str(idx_shelf_edge_0))                
                    txt.close()
                    converged = True
                    #break
                
                try:
                    #   cleanup
                    del model.hdsobj, model.ucnobj, model.cbbobj
                    os.remove(os.path.join(model.foldername, model.modelname + '.hds'))
                    os.remove(os.path.join(model.foldername, 'MT3D001.UCN'))
                    os.remove(os.path.join(model.foldername, 'MT3D.CNF'))            
                    os.remove(os.path.join(model.foldername, model.modelname + '.cbc'))
                    del model
    
                #   except when the model doesnt converge
                except(IndexError, OSError):
                    print('Model has not converged..')
                    #   in that case open the LIST find the cells that have the largest head change in the 1st inner iteraion
                    with open(os.path.join(topo_sc_dir, modelname + '.list'), 'r') as file:
                        for num, line in enumerate(file, 1):
                            if 'MAXIMUM HEAD CHANGE FOR EACH ITERATION (1 INDICATES THE FIRST INNER ITERATION)' in line:
                                print('found at line:', num)
                    pass
                    #sys.exit()

if eq_name == 'RCP_85_AP_00400_to_00500':
    status = df_in.loc[df_in['id_loop'] == id_loop, 'sc_' + str(param_combo)].values[0]
    if status == -1:
        df_in.loc[df_in['id_loop'] == id_loop, 'sc_' + str(param_combo)] = 1
        df_in.update(df_in)
        df_in.to_csv(master_csv_dir, sep = ',', encoding = 'utf-8', index = False)

elif eq_name != 'BP_30000_to_20000':
    pass

else:
    status = df_in.loc[df_in['id_loop'] == id_loop, 'sc_' + str(param_combo)].values[0]
    if status == -1:
        df_in.loc[df_in['id_loop'] == id_loop, 'sc_' + str(param_combo)] = -9
        df_in.update(df_in)
        df_in.to_csv(master_csv_dir, sep = ',', encoding = 'utf-8', index = False)    