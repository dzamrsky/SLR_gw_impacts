# -*- coding: utf-8 -*-
"""
Created on Thu Dec 06 09:49:58 2018

@author: daniel
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import csv
import os.path
import tarfile
import shutil
import collections

"""
Function that checks if there is no ocean/sea on the left side of the model domain, if there is it finds the highest positioned
point of the inland area and clips the model accordingly (keeping the domain to the right side of that point). This simulates the
position of the water divide.

top_lst = model.top_elev
ib_arr = model.ibound_arr
model_obj = model

"""
def trim_model_domain(model_obj, top_lst, ib_arr):
    
    #   first check if there really is a sea on the left side of the model domain, if there actually is one then find the index of the 
    #   left coastline position. Trim the topo values in case there is also sea/ocean at the left side of the land (land looks like an island)
    
    left_side = top_lst[:25]
    neg_lst = [i for i in left_side if i < 0]
    if len(neg_lst) > 0:
        start = 0
        for f in range(100):
            cut_lst = top_lst[start : start + 10]
            if len([i for i in cut_lst if i < 0]) >= 5:
                start += 10
            else:
                break
        idx_neg = start#left_side.index([i for i in left_side if i < 0][-1])
        while top_lst[idx_neg] < -0.:
            idx_neg += 1
        new_st_idx = idx_neg
        
    else:
        idx_neg = 0
        while top_lst[idx_neg] < -0.:
            idx_neg += 1
        new_st_idx = 0
        
    left_cst_idx = idx_neg

    #   check that the landward boundary doesnt have any empty columns on the inland side of the model domain
    for z in range(ib_arr.shape[-1]):
        ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, z].tolist()) if x == 1]
        if ibound_act_lay_idxs != []:
            break#pass
        else:
            new_st_idx += 1
            pass#break
    ib_arr = ib_arr[:, :, new_st_idx:]
        
    #   check if there is any column that contains only inactive cells - if so then cut the IBOUND array to end before that column
    new_end_idx = -1
    for a in range(new_st_idx, model_obj.ibound_arr.shape[-1]):
        ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, a].tolist()) if x == 1]
        if ibound_act_lay_idxs != []:
            pass
        else:
            ib_arr = ib_arr[:, :, : a - new_st_idx]
            new_end_idx = a
            break

    #   second step also checks that the end of the model domain ends at the lowest point, start from the end of the model domain
    if new_end_idx != -1:
        offsh_indx = new_end_idx
    else:
        offsh_indx = len(top_lst)
    
    for b in range(offsh_indx - 11, 10, -10):
        
        top_right_lst = []
        for c in range(b, b + 10):
            top_right_lay = [j for j, x in enumerate(model_obj.ibound_arr[:, 0, c].tolist()) if x == 1][0]
            top_right_lst.append(top_right_lay)
        
        top_left_lst = []
        for d in range(b - 10, b):
            top_lay_left = [j for j, x in enumerate(model_obj.ibound_arr[:, 0, d].tolist()) if x == 1][0]  
            top_left_lst.append(top_lay_left)

        print(np.mean(top_left_lst), np.mean(top_right_lst))
        
        #   if the elevation on the right is higher than on the left, continue
        if np.mean(top_left_lst) >= np.mean(top_right_lst) + 1.0:
            pass
        else:
            new_end_idx = b            
            #new_end_lay = [j for j, x in enumerate(model_obj.ibound_arr[:, 0, b].tolist()) if x == 1][-1]
            #ib_arr = ib_arr[:new_end_lay, :, :b]
            break

    #   now check that the right side of the model domain is not actually even higher than the left side
    if top_lst[new_end_idx] > top_lst[0]:
        print(top_lst[new_end_idx], top_lst[0])
        new_end_idx = top_lst.index(min(top_lst))
        b = new_end_idx

    new_end_lay_1 = [j for j, x in enumerate(model_obj.ibound_arr[:, 0, new_end_idx].tolist()) if x == 1][-1]
    new_end_lay_2 = new_end_lay_1
    
    #for c in range(model_obj.cst_idx, new_end_idx):
    for c in range(int(model_obj.cst_idx / 2), new_end_idx):
        try:
            new_end_lay_3 = [j for j, x in enumerate(model_obj.ibound_arr[:, 0, c].tolist()) if x == 1][-1]
            if new_end_lay_3 > new_end_lay_2:
                new_end_lay_2 = new_end_lay_3
        except IndexError:
            pass
    
    new_end_lay = max(new_end_lay_1, new_end_lay_2)

    if b < new_end_idx:
        ib_arr = ib_arr[:new_end_lay, :, : b]
    else:
        ib_arr = ib_arr[:new_end_lay, :, :]
    """
    #   if the idx_neg is 0 then dont do anything otherwise adjust the model domain and 
    if idx_neg == 0:
        
        #   check for multiple landmasses
        import itertools
        splitted = [list(g) for i, g in itertools.groupby(top_lst,lambda x: x<0)]        
        
        #   sometimes there might be small chunks of below sea level areas, but those are not the real ocean yet
        #   check from the beginning of the splitted list and if there is a short (< 2km) below sea level area
        new_splitted = [x for x in splitted]
        
        g = 1
        while len(new_splitted[g]) <= 20 and new_splitted[g][0] < 0. or new_splitted[g][0] > 0.:
            new_splitted[0] = new_splitted[0] + splitted[g]
            del new_splitted[g]       
   
        if len(new_splitted) > 2 and new_splitted[0][0] > 0. and new_splitted[1][0] < 0.:
            new_end_idx = top_lst.index(min(new_splitted[1]))
            max_elev = max(top_lst[:new_end_idx])
            if top_lst.index(max_elev) > 0:
                max_elev_idx = top_lst.index(max_elev) - 1             
            else:
                max_elev_idx = top_lst.index(max_elev)
                
            new_x_start = model_obj.x_start
            new_ib_arr = ib_arr[:, :, :new_end_idx]
            new_ncol = new_ib_arr.shape[-1]

        else:        
            new_x_start = model_obj.x_start
            new_ib_arr = ib_arr
            new_ncol = ib_arr.shape[-1]
            max_elev = max(top_lst[:])
            max_elev_idx = top_lst.index(max_elev) - 1        
            idx_neg = 0
        
        return new_x_start, new_ib_arr, new_ncol, max_elev_idx, new_end_idx, idx_neg
    
    else:
    """
    #   find the position of the coastline on the right side of the model domain
    #cst_val  = next((x for x in top_lst[left_cst_idx:] if x < 0.0), None)
    #right_cst_idx = left_cst_idx + top_lst[left_cst_idx:].index(cst_val) - 1
    right_cst_idx = - 1 
    #   now find the maximum elevation between these two points 
    try:
        max_elev = max(top_lst[left_cst_idx:right_cst_idx])
        max_elev_idx = top_lst.index(max_elev) #- 1 - left_cst_idx
    except ValueError:
        max_elev = max(top_lst[left_cst_idx:])
        max_elev_idx = top_lst.index(max_elev) #- 1 - left_cst_idx       

    #   change the 
    new_x_start = -(max_elev_idx / 10.) #model_obj.x_start + (max_elev_idx / 10.)
    new_ib_arr = ib_arr[:, :, max_elev_idx : new_end_idx]
    new_end_idx = max_elev_idx + new_ib_arr.shape[-1]
    new_ncol = new_ib_arr.shape[-1]

    return new_x_start, new_ib_arr, new_ncol, max_elev_idx, new_end_idx, idx_neg


def find_lonely_cells(ibound_array):

    col_n = ibound_array.shape[-1]
    row_n = ibound_array.shape[0]

    #   go through the ibound_arr cell by cell
    for j in range(1, col_n - 1):
        for i in range(1, row_n - 1):
            if j == 0:
                left = 0
            else:
                left = ibound_array[i-1,0,j]
            
            if j == (col_n - 1):
                right = 0
            else:
                right = ibound_array[i+1,0,j]
                
            if i == 0:
                up = 0
            else:
                up = ibound_array[i,0,j-1]
                
            if i == (row_n -1):
                down = 0
            else:
                down = ibound_array[i,0,j+1]

            if (left == 0 and right == 0 and up == 0 and down == 0):
                ibound_array[i,0,j] = 0
 
    return ibound_array



"""
Function that reads a CSV file and returns the lognormal values of MU and STDEV depending
on what COSCAT region and which parameter we are looking for.

coscat_id = 16
subreg_id = 1
csv_dir = r'g:\Water_Nexus\_A4_GUM\_MODEL_input_files\_SUBREG_representative_models_GEOLOGY_stats.csv'
cst_type = 'henry_other'
param_type = 'glhymps_bot'
"""

def read_param_from_csv_subreg(coscat_id, subreg_id, csv_dir, param_type):
    #   read in the csv file with pandas and get the list of column names
    df = pd.read_csv(csv_dir, index_col = False, header = 0)
    col_names = list(df.columns.values)
    #   define the name of the column we are looking for based on the function parameters
    col_param_name = param_type  
    match_str = [s for s in col_names if col_param_name in s]
    #   find which string is the mu and which is the stdv
    mu_col = [r for r in match_str if 'mu' in r]
    std_col = [r for r in match_str if 'std' in r]
    #   find the values in the table
    mu_val = df.loc[(df[col_names[0]] == coscat_id) & (df[col_names[1]] == subreg_id), mu_col].values[0]
    std_val = df.loc[(df[col_names[0]] == coscat_id) & (df[col_names[1]] == subreg_id), std_col].values[0]
    #   return the final values
    return mu_val, std_val


def read_param_from_csv(coscat_id, csv_dir, param_type):
    #   read in the csv file with pandas and get the list of column names
    df = pd.read_csv(csv_dir, index_col = False, header = 0)
    col_names = list(df.columns.values)
    #   define the name of the column we are looking for based on the function parameters
    col_param_name = param_type  
    match_str = [s for s in col_names if col_param_name in s]
    #   find which string is the mu and which is the stdv
    mu_col = [r for r in match_str if 'mu' in r]
    std_col = [r for r in match_str if 'std' in r]
    #   find the values in the table
    mu_val = df.loc[df[[col_names[0]] == coscat_id], mu_col].values[0]
    std_val = df.loc[df[col_names[0]] == coscat_id, std_col].values[0]
    #   return the final values
    return mu_val, std_val

"""
Function that generates a set of randomly generated values based on a lognormal distribution.
The input variables are the mu_val and std_val create in the previous function. 

mu_in = mu_val
std_in = std_val
len_lst = 100
title = 'RCH input values per column'
y_label = 'RCH rate in mm/d'
"""

def generate_ln_rand_lst(rand_seed_val, mu_in, std_in, len_lst, title, y_label, fig_out_dir, plot = True):    
    #   set the random seed for this function
    np.random.seed(rand_seed_val)
    #   generate an array of values drawn from the normal distribution
    rand_vals = np.exp(np.random.normal(mu_in, std_in, size = len_lst))    
    #   if the plotting is set to True also give the output plots
    if plot is True:
        #   create an X value list
        x_col = np.linspace(0, len_lst - 1, len_lst)
        #   define window size, output and axes
        fig, ax = plt.subplots(figsize=[8,6])
        ax.bar(x_col, rand_vals)
        #   set plot title
        ax.set_title(title)
        #   set x-axis name
        ax.set_xlabel("MODFLOW column number")
        #   set y-axis name
        ax.set_ylabel(y_label)
        ax = plt.gca()
        ax.yaxis.grid() 
        #   save the figure
        plt.savefig(fig_out_dir, dpi = 300)
        plt.close()
        del fig        
    #   return the list of randomly generated values
    return rand_vals.tolist()
    

"""
Function that finds and provides the geological input parameters for each COSCAT region.

coscat_id = 1413
csv_dir = r'g:\Water_Nexus\_A2\_models_final\coscat_geology_input.csv'


"""

def read_geo_params(coscat_id, csv_dir):
    
    #   get the sand and mud % for the given coscat region
    df = pd.read_csv(csv_dir, index_col = False, header = 0)
    col_names = list(df.columns.values)
    sand_val = df.loc[df[col_names[0]] == coscat_id, 'Sand_gl'].values[0]
    mud_val = df.loc[df[col_names[0]] == coscat_id, 'Mud_gl'].values[0]    
    y_exact_val = df.loc[df[col_names[0]] == coscat_id, 'Y1'].values[0]        
    qs_val = df.loc[df[col_names[0]] == coscat_id, 'QS'].values[0]   
    sm_val = df.loc[df[col_names[0]] == coscat_id, 'SM'].values[0]    
    y_val = df.loc[df[col_names[0]] == coscat_id, 'Y'].values[0]  
    shlf_mud_val = df.loc[df[col_names[0]] == coscat_id, 'Mud_shlf'].values[0]
    slp_mud_val = df.loc[df[col_names[0]] == coscat_id, 'Mud_slp'].values[0]
    #   since the total sum of these two (sand and mud) doesnt amount to 100% 
    add_val = (100. - (sand_val + mud_val)) / 2.
    sand_pct, mud_pct = sand_val + add_val, mud_val + add_val
    #   calculate the % for the reworking parameter - basically multiply by 5 since the
    #   highest Y1 value is ~20 (100/20 = 5)
    y_val_out = y_exact_val * 5.

    #   return the final output values
    return qs_val, sm_val, y_val, y_val_out, sand_pct, mud_pct, shlf_mud_val, slp_mud_val

"""
Function that loads in the final output budget file and for each time step computes how much
recharge gets drained for the total model domain.

budget_in = r'g:\Water_Nexus\_A2\_models_testing_1413\coscat_1413_geo_scenario_1_rch_p_min_et\eq_1\_cbc.npy' 

rch_csv = r'g:\Water_Nexus\_A2\_models_final\_COSCAT_representative_models_RCH_stats.csv'
rch_param_p_min_et_ln = read_param_from_csv(1413, rch_csv, 'henry_other', 'cs_p_min_et_rch')

np.random.seed(1)
#   generate an array of values drawn from the normal distribution
rch_lst_in = np.exp(np.random.normal(rch_param_p_min_et_ln[0], rch_param_p_min_et_ln[1], size = 217))    



"""

def model_res_analysis_rch_drn(budget_in, rch_lst_in):
    
    #   load the budget dictionary
    dict_in = np.load(budget_in).item()
    #   loop through all the keys of the dictionary and get the RCH and DRN arrays
    for key in dict_in.keys():
        
        if key != 0:
            rch_arr = dict_in[key]['q_rch'][1][0]
            drn_arr = dict_in[key]['q_drn']
    
            drn_sum, rch_sum = 0., 0.
    
            for i in range(len(drn_arr)):
                rch_val = rch_arr[i]
                drn_val = drn_arr[i][1]
                rch_in = rch_lst_in[i] * 0.001
                
                drn_sum += drn_val   
                rch_sum += rch_val
                
                print(rch_in, rch_val, drn_val)
                
            print(rch_sum, drn_sum)



#   transform the point water head to fresh water head array
def pt_hd_to_frsh_hd(pt_hd_arr, conc_arr, densref, denseslp, top_elev, bot_lst):
    
    top_bot_lst = [top_elev] + bot_lst
    frsh_hd_arr = pt_hd_arr * 0.0
    
    for i in range(pt_hd_arr.shape[0]):
        for j in range(pt_hd_arr.shape[-1]):
            
            pt_head = pt_hd_arr[i, 0, j]
            conc = denseslp * conc_arr[i, 0, j] + densref
            z = (top_bot_lst[i] + top_bot_lst[i + 1]) / 2
            
            frsh_hd_arr[i, 0, j] = (conc / densref) * pt_head - ((conc - densref) / densref) * z

    return frsh_hd_arr



#   gets the volumes and % of freshwater in different zones of the model domain
"""    
arr_ibound = model.ibound_arr
arr_conc = conc_arr
idx_cst = cst_idx_0  #int(abs(model.x_start) * 10)
idx_shelf_edge = idx_shelf_edge_0  #int((cont_shelf_edge[0] - model.x_start) * 10)

"""
def get_frsh_vol_pct(arr_ibound, arr_conc, idx_cst, idx_shelf_edge):
    #   calculate the volume and % in the inland zone
    tot_inland_cells, fresh_inland_cells = 0, 0
    for a in range(idx_cst):
        ibound_act_lay_idxs = [i for i, x in enumerate(arr_ibound[:, 0, a].tolist()) if x == 1]
        #   loop through the active layers and check whether the concentration is fresh or not
        for x in range(len(ibound_act_lay_idxs)):
            if arr_conc[ibound_act_lay_idxs[x], 0, a] <= 0.1: 
                fresh_inland_cells += 1
        tot_inland_cells += len(ibound_act_lay_idxs)
        
    #   calculate the volume and % in the shelf zone
    tot_shelf_cells, fresh_shelf_cells = 0, 0
    for b in range(idx_cst, min(idx_shelf_edge, arr_conc.shape[-1])):
        ibound_act_lay_idxs = [i for i, x in enumerate(arr_ibound[:, 0, b].tolist()) if x == 1]
        #   loop through the active layers and check whether the concentration is fresh or not
        for y in range(len(ibound_act_lay_idxs)):
            if arr_conc[ibound_act_lay_idxs[y], 0, b] <= 0.1: 
                fresh_shelf_cells += 1
        tot_shelf_cells += len(ibound_act_lay_idxs)
        
    return tot_inland_cells, fresh_inland_cells, tot_shelf_cells, fresh_shelf_cells



#   define a function to plot the GHB boundary and top system summary (DRN and RCH) 
"""    
ghb_cell_head_lst = model.ghb_head_diff_lst
ghb_cell_head_max_diff = model.ghb_head_check_max_diff
ghb_in_lst = model.ghb_input_lst
rch_lst = rch_p_min_et_lst_m_d
q_rch = qrch
q_drn = qdrains
q_ghb = qghb
delcol = model.del_col
ncol = model.ncol

arr_hk = model.hk_arr
cond_lst = model.ghb_input_lst
index_coast = idx_cst 

top_el = model.top_elev
x_st = model.x_start
x_end = model.x_end
index_coast = model.rch_extent
"""

def plot_ghb_drn_rch_summary(index_coast, ghb_cell_head_lst, ghb_cell_head_max_diff, rch_lst, q_rch, q_drn, q_ghb, top_el, x_st, x_end, delcol, ncol, arr_hk, cond_lst, save_dir):
    
    #index_coast = int(abs(x_st) * 10.)
    fig = plt.figure(figsize = (20, 12))
    #   create a figure object and assign all the parts
    ax1 = plt.subplot2grid((3, 1), (0, 0))  #    
    ax2 = plt.subplot2grid((3, 1), (1, 0))  #
    ax3 = plt.subplot2grid((3, 1), (2, 0))  #   
    ax2b = ax2.twinx()  # instantiate a second axes that shares the same x-axis    
    ax3b = ax3.twinx()  # instantiate a second axes that shares the same x-axis        
    #   set the exact positions for each axis
    ax1.set_position([0.05, 0.675, 0.9, 0.25]) # [left, bottom, width, height]
    ax2.set_position([0.05, 0.35, 0.9, 0.25])
    ax3.set_position([0.05, 0.025, 0.9, 0.25])        
 
    #   create the first plot that will show the RCH rates and the GHB difference between assigned and modelled heads
    x_axis_pts = np.linspace(x_st + (delcol / (2 * 1000.)), x_end - (delcol / (2 * 1000.)), ncol)
    
    #   prepare the individual type of internal flow types for plotting
    qrch_to_plot = q_rch[1].tolist()[0]
    qdrn_to_plot = [round(i[1], 6) for i in q_drn.tolist()] + (ncol - len([round(i[1], 6) for i in q_drn.tolist()])) * [0.]
    qghb_to_plot = len([round(i[1], 6) for i in q_drn.tolist()]) * [0.] + [round(i[1], 6) for i in q_ghb.tolist()][:ncol - len([round(i[1], 6) for i in q_drn.tolist()])]
    
    tot_lst = qrch_to_plot + qdrn_to_plot + qghb_to_plot
    
    ax1.bar(x_axis_pts, qrch_to_plot, color = 'blue', width = 0.1, label = 'RCH flow in the cell')    
    ax1.bar(x_axis_pts, qdrn_to_plot, color = 'red', width = 0.1, label = 'DRN flow in the cell')  
    ax1.bar(x_axis_pts, qghb_to_plot, color = 'green', width = 0.1, label = 'GHB flow in the cell')  
    ax1.axis([x_st, x_end, 0., max(list(rch_lst[:index_coast + 1])) * 1000.])
    ax1.set_ylabel('Flow rate (m3/d)')  
    ax1.axhline(y = 0., linewidth = 1, color = 'k')
    ax1.set_xlim([x_st, x_end])
    ax1.set_ylim([math.floor(min(tot_lst) * 10.) / 10., math.ceil(max(tot_lst) * 10.) / 10.])
    ax1.set_title('Top system and boundary flow rates in the model domain', fontsize = 14, y = 1.01)
    ax1.legend(loc = 'best')
    
    rch_to_plot = list([i * 1000. for i in rch_lst[:index_coast]])
    rch_to_plot.extend([0] * (ncol - index_coast))
    ax2.bar(x_axis_pts, rch_to_plot, width = 0.1, label = 'RCH rate')
    #ax2.plot(list(x_axis_pts[:index_coast]), list([i * 1000. for i in rch_lst[:index_coast]]), 'co', markersize = 3, label = 'GEBCO elevation')
    ax2.axis([x_st, x_end, 0., max(list(rch_lst[:index_coast])) * 1000.])
    ax2.set_ylabel('Recharge rate (mm/d)')  

    ghb_to_plt_1 = [y[-1] for y in ghb_cell_head_lst if y[2] >= index_coast]
    ghb_to_plt_y = [0] * (ncol - len(ghb_to_plt_1)) + ghb_to_plt_1
    
    mask1 = np.asarray(ghb_to_plt_y) < 0.0001
    mask2 = np.asarray(ghb_to_plt_y) >= 0.0001
    if x_axis_pts[mask1] != []:
        ax2b.bar(x_axis_pts[mask1], np.asarray(ghb_to_plt_y)[mask1], width = 0.1, color = 'green', label = 'GHB head diff < 0.0001')
    if x_axis_pts[mask2] != []:        
        ax2b.bar(x_axis_pts[mask2], np.asarray(ghb_to_plt_y)[mask2], width = 0.1, color = 'red', label = 'GHB head diff > 0.0001')
    ax2b.set_ylim([0., min(.25, round(max(ghb_to_plt_y), 5))])
    ax2b.set_ylabel('GHB head difference (m) between assigned and simulated head')
    ax2.set_title('Recharge rate and GHB head difference (assigned vs. computed) in the model domain', fontsize = 14, y = 1.01)

    lines, labels = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc = 0)
    
    #   get the conductivity of the top active cell in each column
    hk_lst = []
    for  d in range(ncol):
        hk_lst.append(arr_hk[:, 0, d].tolist()[[i for i, x in enumerate(arr_hk[:, 0, d].tolist()) if not math.isnan(x)][0]])
    lst_cond_to_plot = len([round(i[1], 6) for i in q_drn.tolist()]) * [0.] + [i[-1] for i in cond_lst][:ncol - len([round(i[1], 6) for i in q_drn.tolist()])]
    
    ax3.axis([x_st, x_end, 0., math.ceil(max(hk_lst))])
    ax3.set_ylabel('Hydraulic conductivity (m/d)')  
    ax3b.plot(x_axis_pts, lst_cond_to_plot, linewidth = 0.1, color = 'blue', label = 'GHB conductance')    
    hk_mask1 = np.asarray(hk_lst) < 0.01
    hk_mask2 = np.asarray(hk_lst) >= 0.01  
    if x_axis_pts[hk_mask1] != []:
        ax3.bar(x_axis_pts[hk_mask1], np.asarray([math.ceil(max(hk_lst)) * 0.5] * ncol)[hk_mask1], width = 0.1, color = 'red', label = 'Hk, clay cells')
    if x_axis_pts[hk_mask2] != []:
        ax3.bar(x_axis_pts[hk_mask2], np.asarray(hk_lst)[hk_mask2], width = 0.1, color = 'green', label = 'Hk, non-clay cells')    
    ax3b.set_ylabel('Conductance (m2/d)')  
    ax3b.set_ylim([0., math.ceil(max(lst_cond_to_plot))])
    ax3.set_title('Horizontal hydraulic conductivity and GHB conductance in the model domain', fontsize = 14, y = 1.01)
    
    lines, labels = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3b.get_legend_handles_labels()
    ax3.legend(lines + lines2, labels + labels2, loc=0)
    
    if os.path.exists(save_dir):
        print('Figure   ' + save_dir + '   already exists.')
    else: 
        #   save the figure
        plt.savefig(save_dir, dpi = 300)
    plt.close()
    





"""
Function that wraps up at the end of each model simulation and creates the final output files as netcdf files. The input to the function
is just the directory pointing to the folder system with the model results. 

model_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_1'
summary_folder = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_1\summary_dir'

in_dir = model_dir
summary_dir = os.path.join(summary_folder, '_conc_results')
out_name = 'coscat_1103_param_combo_run_1'

char_dir_conc = os.path.join(char_time_dir, '_conc_results')

dict_ex = np.load(r'g:\Water_Nexus\_A2\_models_param_runs\_0403\all\coscat_0403_param_combo_run_1\eq_00\_conc.npy')

"""

def create_conc_tar(in_dir, summary_dir, char_dir_conc, out_name):
    
    #   initialize the simulation time
    simu_time = 0
    final_dict = dict()     
    
    #   create lists of folders where the _conc dictionaries should be found
    file_list = os.listdir(in_dir)
    file_list_no_char = [i for i in file_list if 'eq' in i and 'time' not in i]
    file_name_list_no_char = [i.split('_')[-1] for i in file_list_no_char]    
    file_name_list_no_char = sorted(file_name_list_no_char)
    #   make a special list for the char times results
    file_list_char = [i for i in file_list if 'eq' in i and 'time' in i]
    file_name_list_char = [i.split('_')[1].split('_')[0] for i in file_list_char]    
    file_name_list_char = sorted(file_name_list_char)
    
    #   loop through the char time list first 
    for file in file_name_list_char:
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(in_dir, 'eq_' + file + '_char_time')
        dict_path = os.path.join(dict_dir, '_conc.npy')
        if os.path.exists(dict_path):
            print('dict_dir :       ' + dict_path)
            final_dict_char_time = dict()  
            #   load the dictionary and individually read the data inside key by key
            conc_dict = np.load(dict_path)[()]
            #   read each key, change the value of the key according to the time stamp (real time since simulation start)
            for key in conc_dict:
                real_time = simu_time + key
                #data = xr.DataArray(np.round(conc_dict[key]['conc'], 3), dims = ('lay', 'row', 'col'), name = real_time)
                final_dict_char_time[real_time] = {'conc': np.round(conc_dict[key]['conc'], 3)}       
            fin_dict_char_time_dir = os.path.join(char_dir_conc, out_name + 'eq_' + file + '_char_time' + '.npy')   
            final_dict_char_time_ordered = collections.OrderedDict(sorted(final_dict_char_time.items()))    
            #   save the dictionary
            np.save(fin_dict_char_time_dir, final_dict_char_time_ordered)

    #   do the same for the list with no_char_time dictionaries
    for file in file_name_list_no_char:    
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(in_dir, 'eq_' + file)
        dict_path = os.path.join(dict_dir, '_conc.npy')
        if os.path.exists(dict_path):
            print('dict_dir :       ' + dict_path)    
            #   load the dictionary and individually read the data inside key by key
            conc_dict = np.load(dict_path)[()]
            #   read each key, change the value of the key according to the time stamp (real time since simulation start)
            for key in conc_dict:
                real_time = simu_time + key
                final_dict[real_time] = {'conc': np.round(conc_dict[key]['conc'], 3)}                
        simu_time += key
  
    #fin_dict_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.npy'    
    #fin_dict_tar_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.tar.xz'
    fin_dict_dir = os.path.join(summary_dir, out_name + '.npy')   
    final_dict_ordered = collections.OrderedDict(sorted(final_dict.items()))        
    np.save(fin_dict_dir, final_dict_ordered)
    #fin_dict_tar_dir =  os.path.join(summary_dir, out_name + '.tar.xz')

    #tar = tarfile.open(fin_dict_tar_dir, "w:xz")
    #tar.add(fin_dict_dir)
    #tar.close()
    #os.remove(fin_dict_dir)
    
"""  Same kind of function, but this time for the heads dictionary   """
def create_head_tar(in_dir, summary_dir, char_dir_head, out_name):

    #   initialize the simulation time
    simu_time = 0
    final_dict = dict()     
    
    #   create lists of folders where the _conc dictionaries should be found
    file_list = os.listdir(in_dir)
    file_list_no_char = [i for i in file_list if 'eq' in i and 'time' not in i]
    file_name_list_no_char = [i.split('_')[-1] for i in file_list_no_char]    
    file_name_list_no_char = sorted(file_name_list_no_char)
    #   make a special list for the char times results
    file_list_char = [i for i in file_list if 'eq' in i and 'time' in i]
    file_name_list_char = [i.split('_')[1].split('_')[0] for i in file_list_char]    
    file_name_list_char = sorted(file_name_list_char)
    
    #   loop through the char time list first 
    for file in file_name_list_char:
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(in_dir, 'eq_' + file + '_char_time')
        dict_path = os.path.join(dict_dir, '_head.npy')
        if os.path.exists(dict_path):
            print('dict_dir :       ' + dict_path)
            final_dict_char_time = dict()  
            #   load the dictionary and individually read the data inside key by key
            head_dict = np.load(dict_path)[()]
            #   read each key, change the value of the key according to the time stamp (real time since simulation start)
            for key in head_dict:
                real_time = simu_time + key
                #data = xr.DataArray(np.round(conc_dict[key]['conc'], 3), dims = ('lay', 'row', 'col'), name = real_time)
                final_dict_char_time[real_time] = {'head': np.round(head_dict[key]['head'], 3)}       
            fin_dict_char_time_dir = os.path.join(char_dir_head, out_name + 'eq_' + file + '_char_time' + '.npy')   
            final_dict_char_time_ordered = collections.OrderedDict(sorted(final_dict_char_time.items()))    
            #   save the dictionary
            np.save(fin_dict_char_time_dir, final_dict_char_time_ordered)

    #   do the same for the list with no_char_time dictionaries
    for file in file_name_list_no_char:    
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(in_dir, 'eq_' + file)
        dict_path = os.path.join(dict_dir, '_head.npy')
        if os.path.exists(dict_path):
            print('dict_dir :       ' + dict_path)    
            #   load the dictionary and individually read the data inside key by key
            head_dict = np.load(dict_path)[()]
            #   read each key, change the value of the key according to the time stamp (real time since simulation start)
            for key in head_dict:
                real_time = simu_time + key
                final_dict[real_time] = {'head': np.round(head_dict[key]['head'], 3)}                
        simu_time += key
  
    #fin_dict_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.npy'    
    #fin_dict_tar_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.tar.xz'
    fin_dict_dir = os.path.join(summary_dir, out_name + '.npy')   
    final_dict_ordered = collections.OrderedDict(sorted(final_dict.items()))        
    np.save(fin_dict_dir, final_dict_ordered)

    
"""  Same kind of function, but this time for the CBC dictionary   """    
def create_cbc_tar(in_dir, summary_dir, char_dir_cbc, out_name):

    #   initialize the simulation time
    simu_time = 0
    final_dict = dict()     
    
    #   create lists of folders where the _conc dictionaries should be found
    file_list = os.listdir(in_dir)
    file_list_no_char = [i for i in file_list if '_to_' in i]
    file_name_list_no_char = file_list_no_char#[i.split('_')[-1] for i in file_list_no_char] 

    #   do the same for the list with no_char_time dictionaries
    for file in file_name_list_no_char:    
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(in_dir, file)
        dict_path = os.path.join(dict_dir, '_cbc.npy')
        if os.path.exists(dict_path):
            print('dict_dir :       ' + dict_path)    
            #   load the dictionary and individually read the data inside key by key
            cbc_dict = np.load(dict_path, allow_pickle = True)[()]
            #   read each key, change the value of the key according to the time stamp (real time since simulation start)
            for key in cbc_dict:
                real_time = simu_time + key
                final_dict[real_time] = {'q_ghb': cbc_dict[key]['q_ghb'],\
                                         'q_drn': cbc_dict[key]['q_drn'],\
                                         'q_rch': cbc_dict[key]['q_rch']}
        simu_time += key
  
    #fin_dict_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.npy'    
    #fin_dict_tar_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.tar.xz'
    fin_dict_dir = os.path.join(summary_dir, out_name + '.npy')   
    final_dict_ordered = collections.OrderedDict(sorted(final_dict.items()))        
    np.save(fin_dict_dir, final_dict_ordered)    
    
    
"""  Same kind of function, but this time for the GHB check dictionary   """  
def create_ghb_check_tar(in_dir, summary_dir, char_dir_ghb_check, out_name):

    #   initialize the simulation time
    simu_time = 0
    final_dict = dict()     
    
    #   create lists of folders where the _conc dictionaries should be found
    file_list = os.listdir(in_dir)
    file_list_no_char = [i for i in file_list if '_to_' in i]
    file_name_list_no_char = file_list_no_char#[i.split('_')[-1] for i in file_list_no_char]    

    #   do the same for the list with no_char_time dictionaries
    for file in file_name_list_no_char:    
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(in_dir, file)
        dict_path = os.path.join(dict_dir, '_ghb_check.npy')
        if os.path.exists(dict_path):
            print('dict_dir :       ' + dict_path)    
            #   load the dictionary and individually read the data inside key by key
            ghb_check_dict = np.load(dict_path, allow_pickle = True)[()]
            #   read each key, change the value of the key according to the time stamp (real time since simulation start)
            for key in ghb_check_dict:
                real_time = simu_time + key
                final_dict[real_time] = {'ghb_check': ghb_check_dict[key]['ghb_check'],\
                                         'ghb_head_diff_lst': ghb_check_dict[key]['ghb_head_diff_lst'],\
                                         'ghb_head_diff_max': ghb_check_dict[key]['ghb_head_diff_max'],\
                                         'succ_cnt': ghb_check_dict[key]['succ_cnt'],\
                                         'fail_cnt': ghb_check_dict[key]['fail_cnt']}   
        simu_time += key
  
    #fin_dict_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.npy'    
    #fin_dict_tar_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.tar.xz'
    fin_dict_dir = os.path.join(summary_dir, out_name + '.npy')   
    final_dict_ordered = collections.OrderedDict(sorted(final_dict.items()))        
    np.save(fin_dict_dir, final_dict_ordered)    

"""  Same kind of function, but this time for the List dictionary   

in_dir = walk_dir
summary_dir = list_dict_dir
char_dir_list = list_dict_dir_ct
out_name = subdir + '_list_dict'

"""  
def create_list_tar(in_dir, summary_dir, char_dir_list, out_name):
    
    #   initialize the simulation time
    simu_time = 0
    final_dict = dict()     
    
    #   create lists of folders where the _conc dictionaries should be found
    file_list = os.listdir(in_dir)
    file_list_no_char = [i for i in file_list if '_to_' in i]
    file_name_list_no_char = file_list_no_char#[i.split('_')[-1] for i in file_list_no_char]    
    
    #   do the same for the list with no_char_time dictionaries
    for file in file_name_list_no_char:    
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(in_dir, file)
        dict_path = os.path.join(dict_dir, '_list.npy')
        if os.path.exists(dict_path):
            print('dict_dir :       ' + dict_path)    
            #   load the dictionary and individually read the data inside key by key
            list_dict = np.load(dict_path, allow_pickle = True)[()]
            #   read each key, change the value of the key according to the time stamp (real time since simulation start)
            for key in list_dict:
                real_time = simu_time + key
                final_dict[real_time] = {'STORAGE_IN': np.round(list_dict[key]['STORAGE_IN'], 4),\
                                         'CONSTANT_HEAD_IN': np.round(list_dict[key]['CONSTANT_HEAD_IN'], 4),\
                                         'DRAINS_IN': np.round(list_dict[key]['DRAINS_IN'], 4),\
                                         'HEAD_DEP_BOUNDS_IN': np.round(list_dict[key]['HEAD_DEP_BOUNDS_IN'], 4),\
                                         'DCDT_IN': np.round(list_dict[key]['DCDT_IN'], 4),\
                                         'TOTAL_IN': np.round(list_dict[key]['TOTAL_IN'], 4),\
                                         'STORAGE_OUT': np.round(list_dict[key]['STORAGE_OUT'], 4),\
                                         'CONSTANT_HEAD_OUT': np.round(list_dict[key]['CONSTANT_HEAD_OUT'], 4),\
                                         'HEAD_DEP_BOUNDS_OUT': np.round(list_dict[key]['HEAD_DEP_BOUNDS_OUT'], 4),\
                                         'RECHARGE_OUT': np.round(list_dict[key]['RECHARGE_OUT'], 4),\
                                         'DCDT_OUT': np.round(list_dict[key]['DCDT_OUT'], 4),\
                                         'TOTAL_OUT': np.round(list_dict[key]['TOTAL_OUT'], 4),\
                                         'IN-OUT': np.round(list_dict[key]['IN-OUT'], 4),\
                                         'PERCENT_DISCREPANCY': np.round(list_dict[key]['PERCENT_DISCREPANCY'], 4)}      
        simu_time += key
  
    #fin_dict_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.npy'    
    #fin_dict_tar_dir = r'g:\\Water_Nexus\\_A2\\_models_testing_1413\\cartesius_runs\\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1\\eq_31_char_time\\_conc_cmp.tar.xz'
    fin_dict_dir = os.path.join(summary_dir, out_name + '.npy')   
    final_dict_ordered = collections.OrderedDict(sorted(final_dict.items()))        
    np.save(fin_dict_dir, final_dict_ordered)    

    
"""
function that combines all the csv results into one final csv ouptut file

dir_in = walk_dir
dir_summary = csv_dir
id_sc = sc_id
folder_name_lst = eq_lst
"""
def create_summary_csv_srm(dir_in, dir_summary, folder_name_lst, id_sc):

    #   initialize the simulation time
    simu_time = 0

    #   add 0s to the scenario id to get a good naming structutre
    if id_sc < 10:
        id_sc_str = '00' + str(id_sc)
    elif id_sc < 100 and id_sc >= 10:
        id_sc_str = '0' + str(id_sc)
    else:
        id_sc_str = str(id_sc)

    #   create the output CSV file, the headers will match the headers of the usual individual csv file
    csv_dir_ts = os.path.join(dir_summary, 'sc_' + id_sc_str + '_summary_EQs_output.csv')
    f = open(csv_dir_ts,'w')
    f_headers = 'EQ, time_start, time_end, max_diff_head, min_diff_head, max_diff_conc, min_diff_conc, fresh_cell_diff_(end_start), brackish_cell_diff_(end_start),\
                saline_cell_diff_(end_start), frsh_inl_vol, frsh_inl_pct, frsh_shelf_vol, frsh_shelf_pct, prev_tot_cells, tot_cells, ghb_cells_check_fail, ghb_cells_check_succ,\
                ghb_cell_dif_abs_max, rch_in, ghb_in, tot_in, drn_out, ghb_out, tot_out, run_time (s)'
    f.write(" ".join(f_headers.split())) 
    f.write('\n')
    f.close()                    
    
    #   do the same for the list with no_char_time dictionaries
    for file in folder_name_lst:    
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(dir_in, file)
        dict_path = os.path.join(dict_dir, '_time_steps_DH_DC.csv')
        
        if os.path.exists(dict_path):
            #   open the csv file 
            with open(dict_path, 'r') as f_tmp:
                #   copy line by line, but adjust some of the values - EQ number and the total simulation time
                for cnt, line in enumerate(f_tmp):
                    #   don't copy the header
                    if cnt != 0:  
                        #   adjust the first three values in the line
                        end_line = line.split(',')[3:]
                        start_line = str(file) + ',' + str(int(float(line.split(',')[1])) + simu_time) + ',' + str(int(float(line.split(',')[2])) + simu_time)
                        new_line = start_line
                        for i in range(len(end_line)):
                            new_line += ',' + end_line[i]
                        #   write the line into the overall file
                        f = open(csv_dir_ts,'a')
                        #f.write("\n")
                        f.write(new_line)
                f.write("\n")
                f.close()                         
            simu_time = int(new_line.split(',')[2])   
           
    """
    #   also find the other summary file, which is only with the % of concentrations and just copy it to csv_results
    for dirpath, dirnames, filenames in os.walk(dir_in):
        for filename in [f for f in filenames if "res_coscat_summary" in f]:
            new_filename = 'scenario_' + id_sc_str + '_pct_fresh_summary.csv'
            csv_dir_old = os.path.join(dirpath, filename) 
            csv_dir_new = os.path.join(dir_summary, new_filename) 
            #   copy the file
            shutil.copy(csv_dir_old, csv_dir_new)
    """
    return csv_dir_ts



def create_summary_csv(dir_in, dir_summary, dir_summary_char_times, id_sc):

    #   initialize the simulation time
    simu_time = 0

    #   add 0s to the scenario id to get a good naming structutre
    if id_sc < 10:
        id_sc_str = '00' + str(id_sc)
    elif id_sc < 100 and id_sc >= 10:
        id_sc_str = '0' + str(id_sc)
    else:
        id_sc_str = str(id_sc)

    #   create the output CSV file, the headers will match the headers of the usual individual csv file
    csv_dir_ts = os.path.join(dir_summary, 'sc_' + id_sc_str + '_summary_EQs_output.csv')
    f = open(csv_dir_ts,'w')
    f_headers = 'EQ, time_start, time_end, max_diff_head, min_diff_head, max_diff_conc, min_diff_conc, fresh_cell_diff_(end_start), brackish_cell_diff_(end_start),\
                saline_cell_diff_(end_start), frsh_inl_vol, frsh_inl_pct, frsh_shelf_vol, frsh_shelf_pct, prev_tot_cells, tot_cells, ghb_cells_check_fail, ghb_cells_check_succ,\
                ghb_cell_dif_abs_max, rch_in, ghb_in, tot_in, drn_out, ghb_out, tot_out, run_time (s)'
    f.write(" ".join(f_headers.split())) 
    f.write('\n')
    f.close()                    
    
    #   create lists of folders where the _conc dictionaries should be found
    file_list = os.listdir(dir_in)
    file_list_no_char = [i for i in file_list if 'eq' in i and 'time' not in i]
    file_name_list_no_char = [i.split('_')[-1] for i in file_list_no_char]    
    file_name_list_no_char = sorted(file_name_list_no_char)
    #   make a special list for the char times results
    file_list_char = [i for i in file_list if 'eq' in i and 'time' in i]
    file_name_list_char = [i.split('_')[1].split('_')[0] for i in file_list_char]    
    file_name_list_char = sorted(file_name_list_char)
    
    #   loop through the char time list first 
    for file in file_name_list_char:
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(dir_in, 'eq_' + file + '_char_time')
        dict_path = os.path.join(dict_dir, '_time_steps_DH_DC.csv')
        if os.path.exists(dict_path):
            dst_dir = os.path.join(dir_summary_char_times, 'scenario_' + id_sc_str + '_eq' + dict_path.split('eq')[-1].replace('/', ''))
            print('dict_dir :       ' + dict_path) 
            print('dst_dir :       ' + dst_dir)    
            shutil.copyfile(dict_path, dst_dir)    
    
    #   do the same for the list with no_char_time dictionaries
    for file in file_name_list_no_char:    
        #   create the conc.npy filename and check if it exists
        dict_dir = os.path.join(dir_in, 'eq_' + file)
        dict_path = os.path.join(dict_dir, '_time_steps_DH_DC.csv')
        if os.path.exists(dict_path):
            
            #   open the csv file 
            with open(dict_path, 'r') as f_tmp:
                #   copy line by line, but adjust some of the values - EQ number and the total simulation time
                for cnt, line in enumerate(f_tmp):
                    #   don't copy the header
                    if cnt != 0:  
                        #   adjust the first three values in the line
                        end_line = line.split(',')[3:]
                        start_line = str(file) + ',' + str(int(float(line.split(',')[1])) + simu_time) + ',' + str(int(float(line.split(',')[2])) + simu_time)
                        new_line = start_line
                        for i in range(len(end_line)):
                            new_line += ',' + end_line[i]
                        if int(float(line.split(',')[1])) + simu_time == 0:
                            print(line)
                            print(new_line)

                        #   write the line into the overall file
                        f = open(csv_dir_ts,'a')
                        #f.write("\n")
                        f.write(new_line)
                f.write("\n")
                f.close()                         
                    
        simu_time = int(new_line.split(',')[2])              
    """
    #   also find the other summary file, which is only with the % of concentrations and just copy it to csv_results
    for dirpath, dirnames, filenames in os.walk(dir_in):
        for filename in [f for f in filenames if "res_coscat_summary" in f]:
            new_filename = 'scenario_' + id_sc_str + '_pct_fresh_summary.csv'
            csv_dir_old = os.path.join(dirpath, filename) 
            csv_dir_new = os.path.join(dir_summary, new_filename) 
            #   copy the file
            shutil.copy(csv_dir_old, csv_dir_new)
    """
    return csv_dir_ts

"""   Function that will tar a source folder to a destination    

src_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\_summary'
dst_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\_summ_test.tar.xz'
"""
def tar_dir(src_dir, dst_dir):

    tar = tarfile.open(dst_dir, "w:xz")
    tar.add(src_dir)
    tar.close()



"""
Function to rename and gather all the figures created during the model simulation into one folder

main_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\coscat_1413_geo_scenarios_GHB_incr_min_100_cond_plots_1'
fig_type_dir = '_ghb_rch_params_pngs'
splitting_str = '_ghb_rch_params_'
dst_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\_summary\_figures'
char_time_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\_summary\_char_time\_figures'

"""
def gather_figs(main_dir, fig_type_dir, splitting_str, char_time_dir, dst_dir):
    
    #   first loop through the model folder and find all the plots () and initialize the simulation time
    simu_time = 0
    dst_char_folder = os.path.join(char_time_dir, fig_type_dir)
    dst_folder = os.path.join(dst_dir, fig_type_dir)
    if os.path.exists(dst_char_folder):
        pass
    else:
        os.makedirs(dst_char_folder)
        
    if os.path.exists(dst_folder):
        pass
    else:
        os.makedirs(dst_folder)   
        
    #   first make a list of all the dictionaries located with that folder structure
    for dirpath, dirnames, filenames in os.walk(main_dir):
        #   if the dirpath contains the fig_type_dir then enter that folder
        if fig_type_dir in dirpath:
            
            time_end_lst = []
            
            for filename in os.walk(dirpath):
                for i in range(len(filename[-1])):
                    if filename[-1][i].endswith(".png"):
                        
                        #   if there is the string "diff" in the filename there are two different time stamps, if not theres only one
                        if 'diff' in filename[-1][i]:
                            #   get the old time stamps
                            str_split = filename[-1][i].split('diff_')[-1].split('_to_')
                            t_start_old = int(str_split[0])
                            t_end_old = int(str_split[-1].split('_yrs')[0])
                            #   create new ones, also add 0s to the string 
                            t_start_new = t_start_old + simu_time
                            t_end_new = t_end_old + simu_time
                            time_end_lst.append(t_end_old)
                            
                            if t_start_new < 1000:
                                t_start_str = '000' + str(t_start_new)
                            elif t_start_new >= 1000 and t_start_new < 10000:
                                t_start_str = '00' + str(t_start_new)
                            elif t_start_new >= 10000 and t_start_new < 100000:
                                t_start_str = '0' + str(t_start_new)
                            else:
                                t_start_str = str(t_start_new)
                                
                            if t_end_new < 1000:
                                t_end_str = '000' + str(t_end_new)
                            elif t_end_new >= 1000 and t_end_new < 10000:
                                t_end_str = '00' + str(t_end_new)
                            elif t_end_new >= 10000 and t_end_new < 100000:
                                t_end_str = '0' + str(t_end_new)
                            else:
                                t_end_str = str(t_end_new)                    

                            #   if there is char_time in the name folder than move those to the correspondant folder
                            if 'char_time' in dirpath:
                                eq_num = int(dirpath.split('eq')[-1].split('char_time/')[0].replace('_', ''))
                                new_name = filename[-1][i].split(str(t_start_old))[0] + t_start_str + '_to_' + t_end_str + '_eq_' + str(eq_num) + '.png'                                  
                                new_dir = os.path.join(dst_char_folder, new_name)
                                old_dir = os.path.join(dirpath, filename[-1][i])
                                #   copy the file
                                shutil.copy(old_dir, new_dir)                            
                            else:
                                eq_num = int(dirpath.split('eq')[-1].split('/')[0].replace('_', ''))
                                new_name = filename[-1][i].split(str(t_start_old))[0] + t_start_str + '_to_' + t_end_str + '_eq_' + str(eq_num) + '.png'                                  
                                new_dir = os.path.join(dst_folder, new_name)
                                old_dir = os.path.join(dirpath, filename[-1][i])
                                #   copy the file
                                shutil.copy(old_dir, new_dir)

                        #   if there is the string "diff" in the filename there are two different time stamps, if not theres only one
                        else:
                            #   get the old time stamps
                            #splitting_str = fig_type_dir.split('pngs')[0].split('_')[1].replace('_', '')
                            str_split = filename[-1][i].split(splitting_str)[-1].split('_yrs')
                            t_old = int(str_split[0])
                            #   create new ones, also add 0s to the string 
                            t_new = t_old + simu_time
                            time_end_lst.append(t_old)
                            
                            if t_new < 1000:
                                t_str = '000' + str(t_new)
                            elif t_new >= 1000 and t_new < 10000:
                                t_str = '00' + str(t_new)
                            elif t_new >= 10000 and t_new < 100000:
                                t_str = '0' + str(t_new)
                            else:
                                t_str = str(t_new)            

                            #   if there is char_time in the name folder than move those to the correspondant folder
                            if 'char_time' in dirpath:
                                eq_num = int(dirpath.split('eq')[-1].split('char_time/')[0].replace('_', ''))
                                new_name = filename[-1][i].split(str(t_old))[0] + t_str + '_eq_' + str(eq_num) + '.png'                                  
                                new_dir = os.path.join(dst_char_folder, new_name)
                                old_dir = os.path.join(dirpath, filename[-1][i])
                                #   copy the file
                                print(old_dir, new_dir)
                                shutil.copy(old_dir, new_dir)                            
                            else:
                                eq_num = int(dirpath.split('eq')[-1].split('/')[0].replace('_', ''))
                                new_name = filename[-1][i].split(str(t_old))[0] + t_str + '_eq_' + str(eq_num) + '.png'                                  
                                new_dir = os.path.join(dst_folder, new_name)
                                old_dir = os.path.join(dirpath, filename[-1][i])
                                #   copy the file
                                print(old_dir, new_dir)
                                shutil.copy(old_dir, new_dir)
                        
            if 'char_time' not in dirpath:   
                simu_time += max(time_end_lst)
                    
    





























