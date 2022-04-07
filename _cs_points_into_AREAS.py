# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 13:10:58 2019

@author: daniel
"""

"""
This script loops through a shapefile/postgresql table with coastal points and looks for continous areas of same coastal types.
If there are few points of different coastal type within a given area those can be neglected (to be specified in script below).

The main idea behind this script is to split the COSCAT regions into subregions based on location of coastal points and their
coastal type. In such way we can maintain the ARP model methdology but gain more detailed averaged profiles fitting better each
smaller region. The final ID number will be composite of the COSCAT_ID number and the sub_region ID that will start at 1 until
all sub-regions in the COSCAT region are assigned. 

"""

#import os
import numpy as np
import geopandas as gpd
from itertools import groupby
from operator import itemgetter

#   read in the original raw shapefile with coastal points and coastal types
cs_pts_types_shp_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\cs_coastal_types_noduplicates.shp'
cs_pts_types_raw = gpd.read_file(cs_pts_types_shp_dir)
cs_pts_types_raw['cst_type'] = np.nan
cs_pts_types_raw['cst_reg_id'] = np.nan
cs_pts_types_raw['cst_reg_id_str'] = np.nan


id_lst = list(set(list(cs_pts_types_raw['id_cs'])))
coscat_lst = list(set(list(cs_pts_types_raw['coscat'])))


"""
Function to split coastal point of each COSCAT region into sub-regions based on coastal types and length of the coastal stretch.

coscat_id = 807
cs_pts_tb = cs_pts_types_raw
min_len = 2
max_gap = 3
max_diff_type = 5

"""

def create_cs_zones(coscat_id, cs_pts_tb, min_len, max_gap, max_diff_type):

    #   select rows from the table to be processed
    rows_coscat = cs_pts_tb.loc[cs_pts_tb['coscat'] == coscat_id]
    #rows_tb = rows_coscat.loc[(rows_coscat['id_cs'] >= st_id) & (rows_coscat['id_cs'] <= end_id)]
    id_cs_loop = sorted(list(set(rows_coscat['id_cs'].values.tolist())))
    
    #   define the right string from the coscat id number
    if coscat_id < 10:
        id_cs_str = '000' + str(coscat_id)
    elif coscat_id >= 10 and coscat_id < 100:
        id_cs_str = '00' + str(coscat_id)
    elif coscat_id >= 100 and coscat_id < 1000:
        id_cs_str = '0' + str(coscat_id)
    else:
        id_cs_str = str(coscat_id)      

    #   find the missing id numbers
    def missing_elements(L):
        start, end = L[0], L[-1]
        return sorted(set(range(start, end + 1)).difference(L))
    missing_gid = missing_elements(id_cs_loop)
    
    #   get the break id numbers - process the above created list to find gaps longer than max_gap
    gap_lst = []
    for k, g in groupby(enumerate(missing_gid), lambda ix : ix[0] - ix[1]):
        group = list(map(itemgetter(1), g))
        gap_lst.append((group))
    gap_lst = [i for i in gap_lst if len(i) >= max_gap]

    #   split the id_cs_loop list into sublists using the indexes from the gap_lst
    id_cs_loop_2 = []
    id_cs_st_idx = 0
    for j in range(len(gap_lst)):
        id_cs_loop_2.append(id_cs_loop[id_cs_st_idx : id_cs_loop.index(gap_lst[j][0] - 1) + 1])
        id_cs_st_idx = id_cs_loop.index(gap_lst[j][0] - 1) + 1
    id_cs_loop_2.append(id_cs_loop[id_cs_st_idx :])
    
    #   remove sublists that are smaller than min_len
    id_cs_loop_2 = [i for i in id_cs_loop_2 if len(i) >= min_len]
    
    #   loop through each sublist of the list created above
    id_cs_type_lst = []
    reg_cnt = 1
    for m in range(len(id_cs_loop_2)):
        chnk_lst = id_cs_loop_2[m]
        #   now loop through the chunk list
        cst_type_lst = []
        for n in range(len(chnk_lst)):
            cst_type = rows_coscat.loc[rows_coscat['id_cs'] == chnk_lst[n], 'class_fin'].values.tolist()[0]
            cst_type_lst.append(cst_type)
        #   loop through the list and check that the coastal indexes differences are with the limit
        prev_cst_type = cst_type_lst[0]
        cst_type_lst_fin = [[chnk_lst[0], prev_cst_type, reg_cnt]]
        for p in range(1, len(cst_type_lst)):
            
            #   if the type is the same then assign the same regional ID number
            if cst_type_lst[p] == prev_cst_type:
                cst_type_lst_fin.append([chnk_lst[p], prev_cst_type, reg_cnt])
            #   if not, check if there is the same coastal type in the next coastal points (max_diff_type + 1)
            elif cst_type_lst[p] != prev_cst_type and prev_cst_type in cst_type_lst[p : p + max_diff_type + 1]:
                cst_type_lst_fin.append([chnk_lst[p], prev_cst_type, reg_cnt])

            #   if even that condition is not satisfied, then go to the next regional id 
            else:
                reg_cnt += 1
                cst_type_lst_fin.append([chnk_lst[p], cst_type_lst[p], reg_cnt])           
                prev_cst_type = cst_type_lst[p]
        
        #   increase the reg_cnt at the end of each chunk
        reg_cnt += 1
        id_cs_type_lst.append(cst_type_lst_fin)

    #   insert the column into the new file
    id_cs_type_lst_to_shp = [item for sublist in id_cs_type_lst for item in sublist]
    
    #   an additional check, if there are sub-regions with less than min_len coastal points, remove them
    reg_id_delete = []
    for z in range(1, reg_cnt):
        #   count number of occurrences
        cs_pt_cnt = len([f for f in id_cs_type_lst_to_shp if f[2] == z])  
        if cs_pt_cnt < min_len:
            reg_id_delete.append(z)
    
    id_cs_type_lst_to_shp = [i for i in id_cs_type_lst_to_shp if i[2] not in reg_id_delete]
    reg_old = list(set([i[-1] for i in id_cs_type_lst_to_shp]))
    reg_new = list(np.arange(1, len(reg_old) + 1))
    
    for g in range(len(reg_old)):
        for u in range(len(id_cs_type_lst_to_shp)):
            if id_cs_type_lst_to_shp[u][2] == reg_old[g]:
                id_cs_type_lst_to_shp[u][2] = reg_new[g]
  
    for r in range(len(id_cs_type_lst_to_shp)):
        cs_pts_types_raw.at[cs_pts_types_raw['id_cs'] == id_cs_type_lst_to_shp[r][0], 'cst_reg_id'] = int(id_cs_type_lst_to_shp[r][2])
        cs_pts_types_raw.at[cs_pts_types_raw['id_cs'] == id_cs_type_lst_to_shp[r][0], 'cst_type'] = int(id_cs_type_lst_to_shp[r][1])
        cs_pts_types_raw.at[cs_pts_types_raw['id_cs'] == id_cs_type_lst_to_shp[r][0], 'cst_reg_id_str'] = id_cs_str + '_' + str(int(id_cs_type_lst_to_shp[r][2]))
    

    
min_len = 5                 #   minimum number of points to be considered for a sub-region
max_gap = 2                 #   maximum gap between id_cs - how many blank points can there be in sub-region set of points
max_diff_type = 2           #   how many points of different coastal type can be skipped and blended into the majority coastal type
cs_pts_tb = cs_pts_types_raw
#coscat_id = 12    
    
#   loop through all the COSCAT regions and create the sub-regions
#all_coscat_cs_reg_lst = []
for coscat_id in coscat_lst:
    print(coscat_id)
    create_cs_zones(coscat_id, cs_pts_tb, min_len, max_gap, max_diff_type)

    
    
# create the GeoDatFrame
gdf = gpd.GeoDataFrame(cs_pts_types_raw, geometry = cs_pts_types_raw.geometry)
gdf = gdf.fillna(-1)
gdf['cst_reg_id'] = gdf['cst_reg_id'].astype('int64')
gdf['cst_type'] = gdf['cst_type'].astype('int64')
# save the GeoDataFrame
gdf.to_file(driver = 'ESRI Shapefile', filename = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\cs_reg_id.shp')

    


"""
Function to split coastal point of each COSCAT region into sub-regions based on coastal types and length of the coastal stretch.

coscat_id = 807
cs_pts_tb = cs_pts_types_raw
min_len = 2
max_gap = 3
max_diff_type = 5
check_file = in_file

"""

def create_cs_zones_v2(coscat_id, cs_pts_tb, min_len, max_gap, max_diff_type, check_file):

    #   select rows from the table to be processed
    rows_coscat = cs_pts_tb.loc[cs_pts_tb['coscat'] == coscat_id]
    #rows_tb = rows_coscat.loc[(rows_coscat['id_cs'] >= st_id) & (rows_coscat['id_cs'] <= end_id)]
    id_cs_loop = sorted(list(set(rows_coscat['id_cs'].values.tolist())))
    
    #   define the right string from the coscat id number
    if coscat_id < 10:
        id_cs_str = '000' + str(coscat_id)
    elif coscat_id >= 10 and coscat_id < 100:
        id_cs_str = '00' + str(coscat_id)
    elif coscat_id >= 100 and coscat_id < 1000:
        id_cs_str = '0' + str(coscat_id)
    else:
        id_cs_str = str(coscat_id)      

    #   find the missing id numbers
    def missing_elements(L):
        start, end = L[0], L[-1]
        return sorted(set(range(start, end + 1)).difference(L))
    missing_gid = missing_elements(id_cs_loop)
    
    #   get the break id numbers - process the above created list to find gaps longer than max_gap
    gap_lst = []
    for k, g in groupby(enumerate(missing_gid), lambda ix : ix[0] - ix[1]):
        group = list(map(itemgetter(1), g))
        gap_lst.append((group))
    gap_lst = [i for i in gap_lst if len(i) >= max_gap]

    #   split the id_cs_loop list into sublists using the indexes from the gap_lst
    id_cs_loop_2 = []
    id_cs_st_idx = 0
    for j in range(len(gap_lst)):
        id_cs_loop_2.append(id_cs_loop[id_cs_st_idx : id_cs_loop.index(gap_lst[j][0] - 1) + 1])
        id_cs_st_idx = id_cs_loop.index(gap_lst[j][0] - 1) + 1
    id_cs_loop_2.append(id_cs_loop[id_cs_st_idx :])
    
    #   remove sublists that are smaller than min_len
    id_cs_loop_2 = [i for i in id_cs_loop_2 if len(i) >= min_len]
    
    #   loop through each sublist of the list created above
    id_cs_type_lst = []
    #   set the counter to match the highest SRM id from previous version
    reg_cnt = max(check_file.loc[check_file['coscat_id'] == coscat_id]['cst_reg_id'].values.tolist()) + 1
    
    for m in range(len(id_cs_loop_2)):
        chnk_lst = id_cs_loop_2[m]
        #   now loop through the chunk list
        cst_type_lst = []
        for n in range(len(chnk_lst)):
            cst_type = rows_coscat.loc[rows_coscat['id_cs'] == chnk_lst[n], 'class_fin'].values.tolist()[0]
            cst_type_lst.append(cst_type)
        #   loop through the list and check that the coastal indexes differences are with the limit
        prev_cst_type = cst_type_lst[0]
        cst_type_lst_fin = [[chnk_lst[0], prev_cst_type, reg_cnt]]
        for p in range(1, len(cst_type_lst)):
            
            #   if the type is the same then assign the same regional ID number
            if cst_type_lst[p] == prev_cst_type:
                cst_type_lst_fin.append([chnk_lst[p], prev_cst_type, reg_cnt])
            #   if not, check if there is the same coastal type in the next coastal points (max_diff_type + 1)
            elif cst_type_lst[p] != prev_cst_type and prev_cst_type in cst_type_lst[p : p + max_diff_type + 1]:
                cst_type_lst_fin.append([chnk_lst[p], prev_cst_type, reg_cnt])

            #   if even that condition is not satisfied, then go to the next regional id 
            else:
                reg_cnt += 1
                cst_type_lst_fin.append([chnk_lst[p], cst_type_lst[p], reg_cnt])           
                prev_cst_type = cst_type_lst[p]
        
        #   increase the reg_cnt at the end of each chunk
        for u in range(len(cst_type_lst_fin)):
            if ((check_file['id_cs'] == cst_type_lst_fin[u][0]) & (check_file['cst_reg_id'] != -1)).any():
                exists = True
                break
            else:
                exists = False
                pass
        
        if not exists:
            reg_cnt += 1
            id_cs_type_lst.append(cst_type_lst_fin)

    #   insert the column into the new file
    id_cs_type_lst_to_shp = [item for sublist in id_cs_type_lst for item in sublist]
    
    #   an additional check, if there are sub-regions with less than min_len coastal points, remove them
    reg_id_delete = []
    for z in range(max(check_file.loc[check_file['coscat_id'] == coscat_id]['cst_reg_id'].values.tolist()) + 1, reg_cnt):
        #   count number of occurrences
        cs_pt_cnt = len([f for f in id_cs_type_lst_to_shp if f[2] == z])  
        if cs_pt_cnt < min_len:
            reg_id_delete.append(z)
    
    id_cs_type_lst_to_shp = [i for i in id_cs_type_lst_to_shp if i[2] not in reg_id_delete]
    reg_old = list(set([i[-1] for i in id_cs_type_lst_to_shp]))
    reg_new = list(np.arange(max(check_file.loc[check_file['coscat_id'] == coscat_id]['cst_reg_id'].values.tolist()) + 1,\
                             len(reg_old) + max(check_file.loc[check_file['coscat_id'] == coscat_id]['cst_reg_id'].values.tolist()) + 1))
    """
    for g in range(len(reg_old)):
        for u in range(len(id_cs_type_lst_to_shp)):
            if id_cs_type_lst_to_shp[u][2] == reg_old[g]:
                id_cs_type_lst_to_shp[u][2] = reg_new[g]
    """
    for u in range(len(id_cs_type_lst_to_shp)):
        id_cs_type_lst_to_shp[u][2] = reg_new[reg_old.index(id_cs_type_lst_to_shp[u][2])]    

    for r in range(len(id_cs_type_lst_to_shp)):
        cs_pts_types_raw.at[cs_pts_types_raw['id_cs'] == id_cs_type_lst_to_shp[r][0], 'cst_reg_id'] = int(id_cs_type_lst_to_shp[r][2])
        cs_pts_types_raw.at[cs_pts_types_raw['id_cs'] == id_cs_type_lst_to_shp[r][0], 'cst_type'] = int(id_cs_type_lst_to_shp[r][1])
        cs_pts_types_raw.at[cs_pts_types_raw['id_cs'] == id_cs_type_lst_to_shp[r][0], 'cst_reg_id_str'] = id_cs_str + '_' + str(int(id_cs_type_lst_to_shp[r][2]))


#   add new SRMs to the list - but check if they dont already exist in the previous one
first_version_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\cs_reg_id.shp'
in_file = gpd.read_file(first_version_dir)

min_len = 2                 #   minimum number of points to be considered for a sub-region
max_gap = 3                 #   maximum gap between id_cs - how many blank points can there be in sub-region set of points
max_diff_type = 5           #   how many points of different coastal type can be skipped and blended into the majority coastal type
cs_pts_tb = cs_pts_types_raw
#coscat_id = 12    
check_file = in_file
    
#   loop through all the COSCAT regions and create the sub-regions
#all_coscat_cs_reg_lst = []
for coscat_id in coscat_lst:
    print(coscat_id)
    create_cs_zones_v2(coscat_id, cs_pts_tb, min_len, max_gap, max_diff_type, in_file)

# create the GeoDatFrame
gdf = gpd.GeoDataFrame(cs_pts_types_raw, geometry = cs_pts_types_raw.geometry)
gdf = gdf.fillna(-1)
gdf['cst_reg_id'] = gdf['cst_reg_id'].astype('int64')
gdf['cst_type'] = gdf['cst_type'].astype('int64')
# save the GeoDataFrame
gdf.to_file(driver = 'ESRI Shapefile', filename = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\cs_reg_id_v2.shp')








    
