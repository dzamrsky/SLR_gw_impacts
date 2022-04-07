# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 11:52:22 2019

@author: daniel
"""

import numpy as np
import os, sys
import xarray as xr

import flopy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import pandas as pd
import math
import cv2

# Turn interactive plotting off
plt.ioff()

#   add the path with all additional tailor made functions to import the scripts
#sys.path.append(r'g:\Water_Nexus\_A2\_scripts_py3\_model_tools_scripts')
sys.path.append(r'/home/dzamrsky/_model_tools_scripts')

#   The main purpose of this script is to wrap up all the output files from different model realizations for each COSCAT region
#   and tar or zip them into one folder. Apart from that, a summary output figures and videos will also be created and inlcuded
#   in the final tar file. The goal is to have only one tar file per COSCAT + coastal type combination. This will provide an easy
#   way how to transfer and strore the final output files in a substantially smaller sized output file. 

#   There are 2 main variables to be specified, the directory to be tarred and the name of the ouptu file - specified as 
#   COSCAT id number and the name of the coastal type. Additionally, a suffix such as 'param_run', 'test_run' or any other can 
#   be added to the final file name for further clarity. On Cartesius these will be specified in tha bash file and fed into this
#   script. On Windows machine for testing purposes they are specified below.


"""
Function that creates an average array out of all the arrays created during all the model runs

in_dict_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1st_cartesius_run'
dict_key = 'conc'
time_key = 1000


"""

def create_avg_array(in_dict_dir, dict_key, time_key):
    
    #   create a list that will be filled with arrays
    arr_lst = []
    
    #   go through the directory and find all the numpy dictionaries 
    for dirpath, dirnames, filenames in os.walk(in_dict_dir):
        for filename in [f for f in filenames if f.endswith(".npy")]:
            dict_dir = os.path.join(dirpath, filename)    
            print(dict_dir)
    
            #   load the dictionary and find the closest time key to the one specified as a parameter
            dict_in = np.load(dict_dir)[()]
            all_keys = list(dict_in.keys())
            closest_time_key = min(all_keys, key = lambda x : abs(x - time_key))
            #   load in the dictionary and transform from 3D to 2D
            dict_val_arr = dict_in[closest_time_key][dict_key]
            dict_val_arr_2d = dict_val_arr[:, 0, :]
    
            #   append the array to the list
            arr_lst.append(dict_val_arr_2d)
    
    #   convert the list to array and calculate the mean array out of that   
    lst_to_arr = np.asarray(arr_lst)
    mean_arr = np.mean(lst_to_arr, axis = 0)

    
"""
Function that gather information from the model MODFLOW package DIS

dis_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\coscat_1413_geo_scenarios_GHB_incr_min_1000_cond_plots_1\eq_00\coscat_1413_geo_scenarios_GHB_incr_min_1000_cond_plots_1.dis'
nam_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\coscat_1413_geo_scenarios_GHB_incr_min_1000_cond_plots_1\eq_00\coscat_1413_geo_scenarios_GHB_incr_min_1000_cond_plots_1.nam_swt'

model_ws = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\coscat_0016_delta_param_combo_run\_param_combo_run_1\eq_00'
model_name = '_param_combo_run_1'



"""

def get_model_info(model_ws, model_name):
    
    #   load the model
    ml = flopy.modflow.Modflow.load(model_name + ".nam_swt", model_ws = model_ws, verbose = False, check = False, exe_name = "mfnwt")

    #   load the DRN package and the DIS package
    drn_sp = ml.drn.stress_period_data[0]
    tot_ncol = ml.dis.ncol
    #   loop through the stress period data and get the number of consecutive drainage cells (in case there is an island further offshore)
    tot_inl_cols = 0 #  because the count in modflow (1 based) is different than in Python (0 based)
    for row in drn_sp:
        if row[2] == tot_inl_cols:
            tot_inl_cols += 1
        else:
            break
    #   now calculate the number of offshore columns - that will be the end of the model
    offshore_cols = tot_ncol - tot_inl_cols
    
    #   assign the start and end and find the top and bottom as well
    x_start, x_end  = tot_inl_cols / (-10.), offshore_cols / (10.)
    top = ml.dis.top[0]
    botm_lst = ml.dis.botm[:, 0, 0].tolist()
    botm_lst = [top] + botm_lst
    bottom = botm_lst[-1]

    #   get the ibound array and the last active cell in the zoom area
    ibound_arr = ml.bas6.ibound[:, :, :]
	
    end_col = min(tot_inl_cols + 100, ibound_arr.shape[-1] - 1)
	
    #ibound_act_lay_idxs = [i for i, x in enumerate(ibound_arr[:, 0, tot_inl_cols + 100].tolist()) if x == 1]
    ibound_act_lay_idxs = [i for i, x in enumerate(ibound_arr[:, 0, end_col].tolist()) if x == 1]
    last_lay = ibound_act_lay_idxs[-1]
    last_bot_zoom = botm_lst[last_lay]

    #   also get the zoom limits
    #   check the extents for the zoomed in area
    if x_start > -10.: 
        x_st_zoom = x_start
        y_top_zoom = math.ceil((top / 100.0) * 100.0)
    else:
        x_st_zoom = -10.
        y_top_zoom = math.ceil((top / 100.0) * 100.0)
    if x_end < 10.:                
        x_end_zoom = x_end
        y_bot_zoom = bottom
    else:
        x_end_zoom = 10.
        y_bot_zoom = math.floor(last_bot_zoom / 100.) * 100.

    #   give those as ouptut
    return x_start, x_end, bottom, top, x_st_zoom, x_end_zoom, y_bot_zoom, y_top_zoom

    
    
"""
Function that creates a final output video out of a dictionary. Including time steps and % of fresh water in the coastal and shelf domains.

conc_dict_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\coscat_1413_geo_scenarios_GHB_incr_min_1000_cond_plots_1\_conc_results_all.npy'
cbc_dict_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\coscat_1413_geo_scenarios_GHB_incr_min_1000_cond_plots_1\_cbc_results_all.npy'
frsh_pct_csv = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\coscat_1413_geo_scenarios_GHB_incr_min_1000_cond_plots_1\_res_coscat_summary_1413_1.csv'
nth_ts = 1000
model_dir = r'g:\Water_Nexus\_A2\_models_testing_1413\cartesius_runs\coscat_1413_geo_scenarios_GHB_incr_min_1000_cond_plots_1'
model_name = 'coscat_1413_geo_scenarios_GHB_incr_min_1000_cond_plots_1'


conc_dict_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_1\coscat_1103_all_param_combo_run_1.npy'
frsh_pct_csv = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_1\scenario_001_pct_fresh_summary_no_char_time.csv'

conc_dict_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_5\_conc_coscat_1103_param_combo_run_5.npy'
cbc_dict_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_5\_cbc_coscat_1103_param_combo_run_5.npy'
frsh_pct_csv = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_5\scenario_005_pct_fresh_summary.csv'

nth_ts = 1000
model_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_1'
model_name = 'coscat_1103_param_combo_run_1'

eq_lst = ['eq_00', 'eq_01', 'eq_02', 'eq_03', 'eq_04', 'eq_05', 'eq_06', 'eq_07', 'eq_08', 'eq_09',\
          'eq_10', 'eq_11', 'eq_12', 'eq_13', 'eq_14', 'eq_15', 'eq_16', 'eq_17', 'eq_18', 'eq_19',\
          'eq_20', 'eq_21', 'eq_21_char_time', 'eq_22', 'eq_23', 'eq_24', 'eq_25', 'eq_26', 'eq_27',\
          'eq_28', 'eq_29', 'eq_30', 'eq_31', 'eq_31_char_time']
sea_lvl_lst = [0., -10., -28., -38., -41., -36., -42, -53., -62., -42., -58., -73., -90., -83., -70., -75.,\
              -80., -81., -78., -90., -105., -120., -120., -108., -96., -84., -72, -60, -48., -36, -24, -12., 0., 0.]

vid_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_1'
img_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1103\coscat_1103_param_combo_run_1\_img_to_vid'

"""

def conc_final_video(conc_dict_dir, cbc_dict_dir, frsh_pct_csv, nth_ts, model_dir, eq_lst, sea_lvl_lst, vid_dir, img_dir, model_name):

    #   create matplotlib objects necessary to export the video. The setup of the final video contains an overall concentration profile,
    #   a zoomed in profile into the current coastal zone (shifting with the sea level), a colorbar with the concentration values and also
    #   a graph that will show (together with written values) the % of fresh water in the coastal (as of 0m asl) and continental shelf 

    #   define the template for showing time
    time_template = 'Time (K years): = %.1f'

    #   load the concentration list and the first (starting) concentration profile
    conc_dict = np.load(conc_dict_dir)[()]
    csv_file = pd.read_csv(frsh_pct_csv)    
    tot_time = 0
    csv_row = 9
    time_sum = 0

    fresh_coast, fresh_shelf, time_lst, lst_sea_lvl = [0], [0], [0], [0.]

    #   calculate the total time span of the whole model simulation, needed for the axis extent
    for index, row in csv_file.iterrows():
        if 'char_time' not in row[0]:
            time_sum += 0.1 
        
    #   get the original (starting) position of the coast and sea level
    model_ws_st = os.path.join(model_dir, 'eq_00')
    dim_limits_st = get_model_info(model_ws_st, model_name)
    st_sea_lvl = 0.0
    st_cst_dist = dim_limits_st[0]   



    img_time_lst = [0]

    #   loop through the time steps and load each nth time step
    for key in range(0, len(conc_dict.keys())):
        time_key = list(conc_dict.keys())[key]
        
        #   check if it is the nth time step
        if time_key % nth_ts == 0:
            print(time_key)

            #   check which EQ it is and also add to the total time
            cur_row = csv_file.iloc[csv_row,:]
            cur_eq, cur_inl_pct, cur_shlf_pct = cur_row[0], cur_row[3], cur_row[5]
            
            if 'char_time' not in cur_eq:
        
                #   check which EQ it is and also add to the total time
                #cur_row = csv_file.iloc[csv_row,:]
                #cur_eq, cur_inl_pct, cur_shlf_pct = cur_row[0], cur_row[3], cur_row[5]
                
                fresh_coast.append(cur_inl_pct)
                fresh_shelf.append(cur_shlf_pct)
                time_lst.append(tot_time)
                
                #   also get the current sea level based on the EQ
                cur_sea_lvl = sea_lvl_lst[eq_lst.index(cur_eq)]
                lst_sea_lvl.append(cur_sea_lvl)
                
                #   create the directory name of the model path
                model_ws = os.path.join(model_dir, cur_eq)
        
                #   get the limits of the model domain
                dim_limits = get_model_info(model_ws, model_name)
                x_st, x_end, bot, top = dim_limits[0] ,dim_limits[1], dim_limits[2], dim_limits[3]
                x_st_zoom, x_end_zoom, bot_zoom, top_zoom = dim_limits[4] ,dim_limits[5], dim_limits[6], dim_limits[7]
        
                #   get the concentration profile at given time
                conc_key = conc_dict[time_key]['conc']
                conc_key[np.abs(conc_key) > 100.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0 
                if time_key == 0:
                    conc_key[np.abs(conc_key) == 0.0] = np.nan
                #conc_key[np.abs(conc_key) < 0.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0   
    
                #   define the constants for the colobrar
                cmap = plt.cm.jet
                cmaplist = [cmap(i) for i in range(cmap.N)]
                cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
                # define the bins and normalize
                bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]            
                norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
                cbar_title = 'Salinity (ppt)'
                
                fig = plt.Figure(figsize = (18, 12))
                ax1 = plt.subplot2grid((3, 3), (0, 1), colspan = 2, fig = fig)              #   overall concentration profiles
                ax2 = plt.subplot2grid((3, 3), (1, 0), colspan = 2, fig = fig)              #   zoomed in concentration profile
                ax3 = plt.subplot2grid((3, 3), (0, 0), fig = fig)                           #   color bar area 
                ax4 = plt.subplot2grid((3, 3), (2, 1), fig = fig)                           #   the legend area for the concentration plots
                ax5 = plt.subplot2grid((3, 3), (1, 2), fig = fig)              #   graph with % fresh area
                ax5b = ax5.twinx()                                                          #   instantiate a second axes that shares the same x-axis    
                ax6 = plt.subplot2grid((3, 3), (2, 0), fig = fig)                           #   area with the time stamp
                ax7 = plt.subplot2grid((3, 3), (2, 2), fig = fig)                           #   the legend area for fresh % graph
                #   specify the location of each of these figure parts
                ax1.set_position([0.175, 0.5, 0.8, 0.45])                        # [left, bottom, width, height]
                ax2.set_position([0.05, 0.15, 0.45, 0.3])
                ax3.set_position([0.05, 0.5, 0.025, 0.45])  
                ax4.set_position([0.2, 0.025, 0.275, 0.075])
                ax5.set_position([0.55, 0.15, 0.4, 0.3])  
                ax6.set_position([0.025, 0.05, 0.15, 0.075])   
                ax7.set_position([0.6, 0.025, 0.3, 0.075])
        
                #   draw the concentration profile for given time step
                im1 = ax1.imshow(conc_key[:, 0, :], aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5, animated = True)
                ax1.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', label = 'shifting sea level', linestyle = '--')
                ax1.axvline(x = 0, ymin = bot, ymax = top, lw = 3., color = 'black', label = 'shifting coastline')   
                #   also the zoomed in plot
                ax2.imshow(conc_key[:, 0, :], aspect = 'auto', interpolation = 'none', cmap = cmap, norm = norm,\
                           extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5)
                
                #   plot the original sea level position and coastline position
                ax1.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey', label = 'current sea level', linestyle = '--')
                ax1.axvline(x = x_st - st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey', label = 'current coastline')   
         
                #   set the limits to the concentration plots
                ax1.set_xlim([x_st, x_end])
                ax1.set_ylim([bot, top])
                ax2.set_xlim([x_st_zoom, x_end_zoom])
                ax2.set_ylim([bot_zoom, top_zoom])
        
                #   set the gridlines and constant lines in the plot add constant lines with elevation = 0m asl. and coastline 
                ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
                ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
                ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
                ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
         
                ax1.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)
                ax2.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)    
               
                #   plot the shifting sea level and coastline position
                ax2.axhline(y = 0, xmin = x_st, xmax = x_end, lw = 3., color = 'black', linestyle = '--')
                ax2.axvline(x = 0, ymin = bot, ymax = top, lw = 3., color = 'black')    
  
                #   plot the original sea level position and coastline position
                ax2.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey')
                ax2.axvline(x = x_st - st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey')    
            
                x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 5.) * 5., math.ceil(x_end), 5.0), nbins = None)    
                x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 1.) * 1., math.ceil(x_end), 1.0), nbins = None)    
                ax1.xaxis.set_major_locator(x_major_locator)
                ax1.xaxis.set_minor_locator(x_minor_locator)
             
                y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 500.) * 500., math.ceil(top), 500.0), nbins = None)    
                y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 100.) * 100. , math.ceil(top), 100.0), nbins = None)    
                ax1.yaxis.set_major_locator(y_major_locator)
                ax1.yaxis.set_minor_locator(y_minor_locator)
                
                ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
                ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
            
                x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 5.) * 5., math.ceil(x_end_zoom), 5.0), nbins = None)    
                x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 1.) * 1., math.ceil(x_end_zoom), 1.0), nbins = None)    
                ax2.xaxis.set_major_locator(x_major_locator_zoom)
                ax2.xaxis.set_minor_locator(x_minor_locator_zoom)
            
                ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
                ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)        
                     
                ax1.set_xlabel('distance from coast (km)', fontsize = 12)
                ax1.set_ylabel('elevation (m asl.)', fontsize = 12)
                ax2.set_xlabel('distance from coast (km)', fontsize = 12)
                ax2.set_ylabel('elevation (m asl.)', fontsize = 12)    
    
                #   plot the colorbar
                cbar = plt.colorbar(im1, cax = ax3, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
                #cbar.ax.set_title(cbar_title, fontsize = 14, y  = 1.025)
                cbar.ax.set_ylabel(cbar_title, rotation = 90)
                cbar.ax.yaxis.set_label_position('left')
                cbar.ax.tick_params(labelsize = 12)    
    
                #   add the time stamp
                txt = ax6.text(0.0, 0.1, time_template % (time_key / 1000.), fontsize = 14, fontweight = 'bold')
                txt.set_clip_on(False)
                ax6.axis('off')
        
                #   create the graph of fresh pct in coastal and shelf zones
                x_time = np.linspace(0, tot_time / 1000., len(fresh_coast))
                lns1 = ax5.plot(x_time, fresh_coast, color = 'blue', linewidth = 2, label = '% fresh in coastal zone')
                lns2 = ax5.plot(x_time, fresh_shelf, color = 'green', linewidth = 2, label = '% fresh in shelf zone')
                ax5.set_xlim([0, round(time_sum, 0)])
                ax5.set_ylim([-5, 105])
                ax5.set_xlabel('Time since start of simulation (K years)', fontsize = 10)
                ax5.set_ylabel('% of fresh water', fontsize = 10)
                
                lns3 = ax5b.plot(x_time, lst_sea_lvl, color = 'cyan', linewidth = 3, label = 'Sea level (m asl)')
                ax5b.set_ylim([-125., 5.])
                ax5b.set_ylabel('Sea level (m asl)')

                #   add the legend for the concentration plots
                h,l = ax1.get_legend_handles_labels() # get labels and handles from ax1
                ax4.legend(h, l, ncol = 2, title = "Concentration plots legend", loc = 10)   
                ax4.axis('off')

                #   add the legend for the fresh % pct
                lns = lns1+ lns2 + lns3
                labs = [l.get_label() for l in lns]
                #ax.legend(lns, labs, loc=0)                
                
                #h2,l2 = ax5.get_legend_handles_labels() # get labels and handles from ax1
                ax7.legend(lns, labs, ncol = 3, title = "Fresh % graph legend", loc = 10)   
                ax7.axis('off')

                figname = 'img_%s.png' % (str(tot_time))

                canvas = FigureCanvas(fig)
                canvas.print_figure(os.path.join(img_dir, figname))

                #   add to the total time and the row in the csv file
                tot_time += nth_ts
                csv_row += int(nth_ts / 100)
                img_time_lst.append(tot_time / 1000.)

                #plt.savefig(os.path.join(img_dir, figname))
                #plt.close(fig)
                
            else:
                key = key - 1
                tot_time -= nth_ts
                pass

    img_array = []
    os.chdir(img_dir)
    
    for i in range(len(img_time_lst) - 1):
        filename = 'img_' + str(int(img_time_lst[i] * 1000)) + '.png'
        print(filename)
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)    
    
    out = cv2.VideoWriter(os.path.join(vid_dir, 'conc_through_time.avi'), cv2.VideoWriter_fourcc(*'DIVX'), 5, size)
     
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()




"""
nc_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\_summary_0016_delta_param_combo_run\_nc\_param_combo_run_14'
nth_ts = 1000
vid_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\_summary_0016_delta_param_combo_run\_videos'
img_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\_summary_0016_delta_param_combo_run\_img_to_vid\_param_combo_run_14'
csv_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\_summary_0016_delta_param_combo_run\_csv_results\scenario_014_summary_EQs_output.csv'
model_dir = modeldir
modelname = model_name + '_' + str(14)



model_dir = walk_dir
nc_dir = walk_dir_nc
nth_ts = 1000
eq_lst = equilibrium_names
sea_lvl_lst = sea_levels
modelname = modelname
vid_dir = vid_dir
img_dir = walk_dir_img_to_vid
csv_dir = sum_csv_dir


(walk_dir, walk_dir_nc, 1000, eq_lst, sea_lvl_lst, subdir, vid_dir, walk_dir_img_to_vid, sum_csv_dir)

"""

def conc_final_video_from_nc_srm(model_dir, nc_dir, nth_ts, eq_lst, sea_lvl_lst, modelname, vid_dir, img_dir, csv_dir):

    #   create matplotlib objects necessary to export the video. The setup of the final video contains an overall concentration profile,
    #   a zoomed in profile into the current coastal zone (shifting with the sea level), a colorbar with the concentration values and also
    #   a graph that will show (together with written values) the % of fresh water in the coastal (as of 0m asl) and continental shelf 

    #   define the template for showing time
    time_template = 'Time (K years): = %.2f'

    #   load the concentration list and the first (starting) concentration profile
    tot_time = 0
    time_sum = 0

    fresh_coast, fresh_shelf, time_lst, lst_sea_lvl = [], [], [], []

    #   get a sorted list of all the netcdf output files
    ts_lst = []
    for path, subdirs, files in os.walk(nc_dir):
        for name in files:
            if not 'RCP' in name:
                print(os.path.join(path, name))  
                ts_lst.append(int(name.split('_')[-1].split('.')[0]))
    ts_lst = sorted(ts_lst)
            
    time_sum = ts_lst[-1] / 1000.
    
    #   get the original (starting) position of the coast and sea level
    st_sea_lvl = 0.0
    st_cst_dist = 0.0

    img_time_lst = []

    #   read in the csv file with the fresh percentages
    df_all = pd.read_csv(csv_dir)
    
    #   get the coast position and vertical extent for each EQ
    x_st_lst, top_zoom_lst, bot_zoom_lst  = [], [], []
    for eq_per in eq_lst:
        #   get the limits of the model domain
        try:
            dim_limits = get_model_info(os.path.join(model_dir, eq_per), modelname)    
            x_st, bot_zoom, top_zoom = dim_limits[0], dim_limits[6], dim_limits[7]            
            x_st_lst.append(x_st)
            top_zoom_lst.append(top_zoom)
            bot_zoom_lst.append(bot_zoom)
        except OSError:
            print('Model did not finish')
            break

    for eq in eq_lst:
        #   create empty list to collect all the different time steps values (will then sort the list based on the value)
        eq_ts = []    
        for dirpath, dirnames, filenames in os.walk(os.path.join(model_dir, eq, '_nc')):
            for filename in [f for f in filenames if f.endswith('.nc')]:
                #print(os.path.join(dirpath, filename), os.path.join(walk_dir_nc, filename))
                #print('-------  Copying all the .nc files for model run :     ' + subdir, filename)
                eq_ts.append(int(filename.split('_')[-1].split('.')[0]))
        #   sort the list because of time counting        
        if eq == 'BP_30000_to_20000':
            eq_ts = sorted(eq_ts)  # 1: to skip the time step 0
            tot_time = 0
            sim_time = 0
            for eq_t in eq_ts:
                if eq_t == 0:
                    #   open the netcdf file
                    nc_name = os.path.join(model_dir, '_INIT_dict.npy')
                    nc_in = np.load(nc_name, allow_pickle = True).item()
                    conc_key = nc_in['dis_info']['sconc_arr']       
                    ib_arr = nc_in['dis_info']['ibound_arr']     
                    for i in range(ib_arr.shape[0]):
                        for j in range(ib_arr.shape[-1]):
                            if ib_arr[i, 0, j] == 0:
                                conc_key[i, 0, j] = np.nan 
                    conc_key = conc_key[:,0,:]    
                    cur_inl_pct = df_all.loc[df_all[' time_end'] == 500][' frsh_inl_pct'].item()
                    cur_shlf_pct = df_all.loc[df_all[' time_end'] == 500][' frsh_shelf_pct'].item()   
                    #   check which EQ it is and also add to the total time                
                    fresh_coast.append(round(np.array(cur_inl_pct, ndmin = 1)[0], 2))
                    fresh_shelf.append(round(np.array(cur_shlf_pct, ndmin = 1)[0], 2))
                    time_lst.append(tot_time)
                    #   also get the current sea level based on the EQ
                    cur_eq_idx = eq_lst.index(eq)
                    cur_sea_lvl = sea_lvl_lst[cur_eq_idx]
                    lst_sea_lvl.append(cur_sea_lvl)     
                    tot_time += 500  
                    sim_time += 500  
                else:
                    #   open the netcdf file
                    nc_name = os.path.join(nc_dir, 'nc_' + str(tot_time) + '.nc')
                    nc_in = xr.open_dataset(nc_name)
                    #cur_eq = df_all.loc[df_all[' time_end'] == sim_time]['EQ'].item()
                    cur_inl_pct = df_all.loc[df_all[' time_end'] == sim_time][' frsh_inl_pct'].item()
                    cur_shlf_pct = df_all.loc[df_all[' time_end'] == sim_time][' frsh_shelf_pct'].item()                         
                    #   check which EQ it is and also add to the total time                
                    fresh_coast.append(round(np.array(cur_inl_pct, ndmin = 1)[0], 2))
                    fresh_shelf.append(round(np.array(cur_shlf_pct, ndmin = 1)[0], 2))
                    time_lst.append(tot_time)                    
                    #   also get the current sea level based on the EQ
                    cur_eq_idx = eq_lst.index(eq)
                    cur_sea_lvl = sea_lvl_lst[cur_eq_idx]
                    lst_sea_lvl.append(cur_sea_lvl)     
                    tot_time += 500  
                    sim_time += 500  
            #   if the model converged earlier then copy the last conc profile to fill the whole stress period
            if tot_time < 10500:
                for ts_fill in range(int((10500 - tot_time) / 500)):
                    fresh_coast.append(round(np.array(cur_inl_pct, ndmin = 1)[0], 2))
                    fresh_shelf.append(round(np.array(cur_shlf_pct, ndmin = 1)[0], 2))
                    lst_sea_lvl.append(cur_sea_lvl)     
                    time_lst.append(tot_time)                    
                    tot_time += 500          

        else:
            eq_ts_sorted = sorted(eq_ts)[1:]  # 1: to skip the time step 0
            time_add = eq_ts_sorted[-1] - eq_ts_sorted[-2]
            for ts in eq_ts_sorted:
                #   open the netcdf file
                nc_name = os.path.join(nc_dir, 'nc_' + str(tot_time) + '.nc')
                nc_in = xr.open_dataset(nc_name)
                #cur_eq = df_all.loc[df_all[' time_end'] == sim_time]['EQ'].item()
                cur_inl_pct = df_all.loc[df_all[' time_end'] == sim_time][' frsh_inl_pct'].item()
                cur_shlf_pct = df_all.loc[df_all[' time_end'] == sim_time][' frsh_shelf_pct'].item()                         
                #   check which EQ it is and also add to the total time                
                fresh_coast.append(round(np.array(cur_inl_pct, ndmin = 1)[0], 2))
                fresh_shelf.append(round(np.array(cur_shlf_pct, ndmin = 1)[0], 2))
                time_lst.append(tot_time)          
                cur_eq_idx = eq_lst.index(eq)
                cur_sea_lvl = sea_lvl_lst[cur_eq_idx]
                lst_sea_lvl.append(cur_sea_lvl)  
                tot_time += time_add   
                sim_time += time_add
            if eq == 'BP_03000_to_02000':
                tot_time -= 450
                sim_time -= 450

    #   create all the plots
    nc_name = os.path.join(nc_dir, 'nc_' + str(500) + '.nc')
    nc_in = xr.open_dataset(nc_name)
    x_lst =  nc_in.coords['x'].values
    y_lst =  nc_in.coords['y'].values
    x_st = x_lst[0]
    x_end = x_lst[-1]
    bot = y_lst[-1]
    top = y_lst[0]

    for a in range(len(ts_lst)):
        ts = ts_lst[a]
        nc_name = os.path.join(nc_dir, 'nc_' + str(ts) + '.nc')
        nc_in = xr.open_dataset(nc_name)
        #   get the concentration profile at given time
        conc_key = nc_in['solute concentration'].values
        conc_key[np.abs(conc_key) > 100.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0 
        if ts == 0:
            conc_key[np.abs(conc_key) == 0.0] = np.nan

        cur_sea_lvl = lst_sea_lvl[a]
        lst_lvl = lst_sea_lvl[:a + 1]
        fresh_cst = fresh_coast[:a + 1]
        fresh_shlf = fresh_shelf[:a + 1]

        #x_cur_cst = ((x_lst[0] - 0.05) - x_st_lst[cur_eq_idx])
        x_st_zoom = max(-10., x_st)
        x_end_zoom = x_st_zoom + 10.
        bot_zoom = max(bot_zoom_lst)#[cur_eq_idx]
        top_zoom = max(top_zoom_lst)#[cur_eq_idx]   
        
        #   define the constants for the colobrar
        cmap = plt.cm.jet
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
        # define the bins and normalize
        bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]            
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        cbar_title = 'Salinity (ppt)'
        
        fig = plt.Figure(figsize = (18, 12))
        ax1 = plt.subplot2grid((3, 3), (0, 1), colspan = 2, fig = fig)              #   overall concentration profiles
        ax2 = plt.subplot2grid((3, 3), (1, 0), colspan = 2, fig = fig)              #   zoomed in concentration profile
        ax3 = plt.subplot2grid((3, 3), (0, 0), fig = fig)                           #   color bar area 
        ax4 = plt.subplot2grid((3, 3), (2, 1), fig = fig)                           #   the legend area for the concentration plots
        ax5 = plt.subplot2grid((3, 3), (1, 2), fig = fig)              #   graph with % fresh area
        ax5b = ax5.twinx()                                                          #   instantiate a second axes that shares the same x-axis    
        ax6 = plt.subplot2grid((3, 3), (2, 0), fig = fig)                           #   area with the time stamp
        ax7 = plt.subplot2grid((3, 3), (2, 2), fig = fig)                           #   the legend area for fresh % graph
        #   specify the location of each of these figure parts
        ax1.set_position([0.175, 0.5, 0.8, 0.45])                        # [left, bottom, width, height]
        ax2.set_position([0.05, 0.15, 0.45, 0.3])
        ax3.set_position([0.05, 0.5, 0.025, 0.45])  
        ax4.set_position([0.2, 0.025, 0.275, 0.075])
        ax5.set_position([0.55, 0.15, 0.4, 0.3])  
        ax6.set_position([0.025, 0.05, 0.15, 0.075])   
        ax7.set_position([0.6, 0.025, 0.3, 0.075])

        #   draw the concentration profile for given time step
        im1 = ax1.imshow(conc_key[:, :], aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                   extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5, animated = True)
        ax1.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', label = 'shifting sea level', linestyle = '--')
        #ax1.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black', label = 'shifting coastline')   
        #   also the zoomed in plot
        ax2.imshow(conc_key[:, :], aspect = 'auto', interpolation = 'none', cmap = cmap, norm = norm,\
                   extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5)
        
        #   plot the original sea level position and coastline position
        ax1.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey', label = 'current sea level', linestyle = '--')
        ax1.axvline(x = x_st - st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey', label = 'current coastline')   
 
        #   set the limits to the concentration plots
        ax1.set_xlim([x_st, x_end])
        ax1.set_ylim([bot, top])
        ax2.set_xlim([x_st_zoom, x_end_zoom])
        ax2.set_ylim([bot_zoom, top_zoom])

        #   set the gridlines and constant lines in the plot add constant lines with elevation = 0m asl. and coastline 
        ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
 
        ax1.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)
        ax2.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)    
       
        #   plot the shifting sea level and coastline position
        #ax2.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', linestyle = '--')
        #ax2.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black')    
  
        #   plot the original sea level position and coastline position
        ax2.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', linestyle = '--')
        ax2.axvline(x = st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'black', linestyle = '--')    
    
        x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 5.) * 5., math.ceil(x_end), 5.0), nbins = None)    
        x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 1.) * 1., math.ceil(x_end), 1.0), nbins = None)    
        ax1.xaxis.set_major_locator(x_major_locator)
        ax1.xaxis.set_minor_locator(x_minor_locator)
     
        y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 500.) * 500., math.ceil(top), 500.0), nbins = None)    
        y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 100.) * 100. , math.ceil(top), 100.0), nbins = None)    
        ax1.yaxis.set_major_locator(y_major_locator)
        ax1.yaxis.set_minor_locator(y_minor_locator)
        
        ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
        ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
    
        x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 5.) * 5., math.ceil(x_end_zoom), 5.0), nbins = None)    
        x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 1.) * 1., math.ceil(x_end_zoom), 1.0), nbins = None)    
        ax2.xaxis.set_major_locator(x_major_locator_zoom)
        ax2.xaxis.set_minor_locator(x_minor_locator_zoom)
    
        ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
        ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)        
             
        ax1.set_xlabel('distance from coast (km)', fontsize = 12)
        ax1.set_ylabel('elevation (m asl.)', fontsize = 12)
        ax2.set_xlabel('distance from coast (km)', fontsize = 12)
        ax2.set_ylabel('elevation (m asl.)', fontsize = 12)    

        #   plot the colorbar
        cbar = plt.colorbar(im1, cax = ax3, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
        #cbar.ax.set_title(cbar_title, fontsize = 14, y  = 1.025)
        cbar.ax.set_ylabel(cbar_title, rotation = 90)
        cbar.ax.yaxis.set_label_position('left')
        cbar.ax.tick_params(labelsize = 12)    

        #   add the time stamp
        txt = ax6.text(0.0, 0.1, time_template % (ts / 1000.), fontsize = 14, fontweight = 'bold')
        txt.set_clip_on(False)
        ax6.axis('off')

        #   create the graph of fresh pct in coastal and shelf zones
        try:
            x_time = np.array([b / 1000 for b in ts_lst[:a + 1]])#np.linspace(0, tot_time / 1000., len(fresh_coast))
        except IndexError:
            x_time = np.array([b / 1000 for b in ts_lst])#np.linspace(0, tot_time / 1000., len(fresh_coast))            
        lns1 = ax5.plot(x_time, fresh_cst, color = 'blue', linewidth = 2, label = '% fresh in coastal zone')
        lns2 = ax5.plot(x_time, fresh_shlf, color = 'green', linewidth = 2, label = '% fresh in shelf zone')
        ax5.set_xlim([0, round(time_sum, 0)])
        ax5.set_ylim([-5, 105])
        ax5.set_xlabel('Time since start of simulation (K years)', fontsize = 10)
        ax5.set_ylabel('% of fresh water', fontsize = 10)
        
        lns3 = ax5b.plot(x_time, lst_lvl, color = 'cyan', linewidth = 3, label = 'Sea level (m asl)')
        ax5b.set_ylim([-131., 5.])
        ax5b.set_ylabel('Sea level (m asl)')

        #   add the legend for the concentration plots
        h,l = ax1.get_legend_handles_labels() # get labels and handles from ax1
        ax4.legend(h, l, ncol = 2, title = "Concentration plots legend", loc = 10)   
        ax4.axis('off')

        #   add the legend for the fresh % pct
        lns = lns1+ lns2 + lns3
        labs = [l.get_label() for l in lns]
        #ax.legend(lns, labs, loc=0)                
        
        #h2,l2 = ax5.get_legend_handles_labels() # get labels and handles from ax1
        ax7.legend(lns, labs, ncol = 3, title = "Fresh % graph legend", loc = 10)   
        ax7.axis('off')

        figname = 'img_%s.png' % (str(ts))

        canvas = FigureCanvas(fig)
        canvas.print_figure(os.path.join(img_dir, figname))

        #   add to the total time and the row in the csv file
        img_time_lst.append(ts)

        #plt.savefig(os.path.join(img_dir, figname))
        #plt.close(fig)
            
        #else:
        #    key = key - 1
        #    tot_time -= nth_ts
        #    pass

    img_array = []
    os.chdir(img_dir)
    
    for i in range(len(img_time_lst)):
        filename = 'img_' + str(int(img_time_lst[i])) + '.png'
        print(filename)
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)    
    
    out = cv2.VideoWriter(os.path.join(vid_dir, model_dir.split('\\')[-1] + '.avi'), cv2.VideoWriter_fourcc(*'DIVX'), 5, size)
    
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()


def conc_final_video_from_nc(model_dir, nc_dir, nth_ts, eq_lst, sea_lvl_lst, modelname, vid_dir, img_dir, csv_dir):

    #   create matplotlib objects necessary to export the video. The setup of the final video contains an overall concentration profile,
    #   a zoomed in profile into the current coastal zone (shifting with the sea level), a colorbar with the concentration values and also
    #   a graph that will show (together with written values) the % of fresh water in the coastal (as of 0m asl) and continental shelf 

    #   define the template for showing time
    time_template = 'Time (K years): = %.1f'

    #   load the concentration list and the first (starting) concentration profile
    tot_time = 0
    csv_row = 9
    time_sum = 0

    fresh_coast, fresh_shelf, time_lst, lst_sea_lvl = [0], [0], [0], [0.]

    #   get a sorted list of all the netcdf output files
    ts_lst = []
    for path, subdirs, files in os.walk(nc_dir):
        for name in files:
            print(os.path.join(path, name))  
            ts_lst.append(int(name.split('_')[-1].split('.')[0]))
    ts_lst = sorted(ts_lst)
            
    time_sum = ts_lst[-1] / 1000.
    
    #   get the original (starting) position of the coast and sea level
    st_sea_lvl = 0.0
    st_cst_dist = 0.0

    img_time_lst = [0]

    #   read in the csv file with the fresh percentages
    df_all = pd.read_csv(csv_dir)
    
    #   get the coast position and vertical extent for each EQ
    x_st_lst, top_zoom_lst, bot_zoom_lst  = [], [], []
    for eq_per in eq_lst:
        if 'char_time' not in eq_per:
            #   get the limits of the model domain
            try:
                dim_limits = get_model_info(os.path.join(model_dir, eq_per), modelname)    
                x_st, bot_zoom, top_zoom = dim_limits[0], dim_limits[6], dim_limits[7]            
                x_st_lst.append(x_st)
                top_zoom_lst.append(top_zoom)
                bot_zoom_lst.append(bot_zoom)
            except OSError:
                print('Model did not finish')
                break

    #   loop through the time steps and load each nth time step
    for key in range(0, len(ts_lst)):
        time_key = ts_lst[key]
        
        #   check if it is the nth time step
        if time_key % nth_ts == 0:
            print(time_key)

            #   open the netcdf file
            nc_name = os.path.join(nc_dir, 'nc_' + str(time_key) + '.nc')
            nc_in = xr.open_dataset(nc_name)

            #   check which EQ it is and also add to the total time
            if time_key == 0:
                cur_eq = 0
                cur_inl_pct = 0
                cur_shlf_pct = 0 
            else:               
                cur_eq = df_all.loc[df_all[' time_end'] == time_key]['EQ'].item()
                cur_inl_pct = df_all.loc[df_all[' time_end'] == time_key][' frsh_inl_pct'].item()
                cur_shlf_pct = df_all.loc[df_all[' time_end'] == time_key][' frsh_shelf_pct'].item()     

            #   check which EQ it is and also add to the total time                
            fresh_coast.append(round(np.array(cur_inl_pct, ndmin = 1)[0], 2))
            fresh_shelf.append(round(np.array(cur_shlf_pct, ndmin = 1)[0], 2))
            time_lst.append(tot_time)
            
            #   also get the current sea level based on the EQ
            if cur_eq < 10:
                eq_str = '0' + str(cur_eq)
            else:
                eq_str = str(cur_eq)
                
            cur_eq_idx = eq_lst.index('eq_' + eq_str)

            cur_sea_lvl = sea_lvl_lst[cur_eq_idx]
            lst_sea_lvl.append(cur_sea_lvl)
            
            #   we can only read in the lists once
            if time_key == 0:
                x_lst =  nc_in.coords['x'].values
                y_lst =  nc_in.coords['y'].values
                x_st = x_lst[0]
                x_end = x_lst[-1]
                bot = y_lst[-1]
                top = y_lst[0]
            
            x_cur_cst = ((x_lst[0] - 0.05) - x_st_lst[cur_eq_idx])
            x_st_zoom = x_cur_cst - 10.
            x_end_zoom = x_cur_cst + 10.
            bot_zoom = bot_zoom_lst[cur_eq_idx]
            top_zoom = top_zoom_lst[cur_eq_idx]   
            
            #   get the concentration profile at given time
            conc_key = nc_in['solute concentration'].values
            conc_key[np.abs(conc_key) > 100.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0 
            if time_key == 0:
                conc_key[np.abs(conc_key) == 0.0] = np.nan
            #conc_key[np.abs(conc_key) < 0.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0   

            #   define the constants for the colobrar
            cmap = plt.cm.jet
            cmaplist = [cmap(i) for i in range(cmap.N)]
            cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
            # define the bins and normalize
            bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]            
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
            cbar_title = 'Salinity (ppt)'
            
            fig = plt.Figure(figsize = (18, 12))
            ax1 = plt.subplot2grid((3, 3), (0, 1), colspan = 2, fig = fig)              #   overall concentration profiles
            ax2 = plt.subplot2grid((3, 3), (1, 0), colspan = 2, fig = fig)              #   zoomed in concentration profile
            ax3 = plt.subplot2grid((3, 3), (0, 0), fig = fig)                           #   color bar area 
            ax4 = plt.subplot2grid((3, 3), (2, 1), fig = fig)                           #   the legend area for the concentration plots
            ax5 = plt.subplot2grid((3, 3), (1, 2), fig = fig)              #   graph with % fresh area
            ax5b = ax5.twinx()                                                          #   instantiate a second axes that shares the same x-axis    
            ax6 = plt.subplot2grid((3, 3), (2, 0), fig = fig)                           #   area with the time stamp
            ax7 = plt.subplot2grid((3, 3), (2, 2), fig = fig)                           #   the legend area for fresh % graph
            #   specify the location of each of these figure parts
            ax1.set_position([0.175, 0.5, 0.8, 0.45])                        # [left, bottom, width, height]
            ax2.set_position([0.05, 0.15, 0.45, 0.3])
            ax3.set_position([0.05, 0.5, 0.025, 0.45])  
            ax4.set_position([0.2, 0.025, 0.275, 0.075])
            ax5.set_position([0.55, 0.15, 0.4, 0.3])  
            ax6.set_position([0.025, 0.05, 0.15, 0.075])   
            ax7.set_position([0.6, 0.025, 0.3, 0.075])
    
            #   draw the concentration profile for given time step
            im1 = ax1.imshow(conc_key[:, :], aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                       extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5, animated = True)
            ax1.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', label = 'shifting sea level', linestyle = '--')
            ax1.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black', label = 'shifting coastline')   
            #   also the zoomed in plot
            ax2.imshow(conc_key[:, :], aspect = 'auto', interpolation = 'none', cmap = cmap, norm = norm,\
                       extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5)
            
            #   plot the original sea level position and coastline position
            ax1.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey', label = 'current sea level', linestyle = '--')
            ax1.axvline(x = x_st - st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey', label = 'current coastline')   
     
            #   set the limits to the concentration plots
            ax1.set_xlim([x_st, x_end])
            ax1.set_ylim([bot, top])
            ax2.set_xlim([x_st_zoom, x_end_zoom])
            ax2.set_ylim([bot_zoom, top_zoom])
    
            #   set the gridlines and constant lines in the plot add constant lines with elevation = 0m asl. and coastline 
            ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
            ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
            ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
            ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
     
            ax1.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)
            ax2.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)    
           
            #   plot the shifting sea level and coastline position
            ax2.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', linestyle = '--')
            ax2.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black')    
  
            #   plot the original sea level position and coastline position
            ax2.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey')
            ax2.axvline(x = st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey')    
        
            x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 5.) * 5., math.ceil(x_end), 5.0), nbins = None)    
            x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 1.) * 1., math.ceil(x_end), 1.0), nbins = None)    
            ax1.xaxis.set_major_locator(x_major_locator)
            ax1.xaxis.set_minor_locator(x_minor_locator)
         
            y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 500.) * 500., math.ceil(top), 500.0), nbins = None)    
            y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 100.) * 100. , math.ceil(top), 100.0), nbins = None)    
            ax1.yaxis.set_major_locator(y_major_locator)
            ax1.yaxis.set_minor_locator(y_minor_locator)
            
            ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
            ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
        
            x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 5.) * 5., math.ceil(x_end_zoom), 5.0), nbins = None)    
            x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 1.) * 1., math.ceil(x_end_zoom), 1.0), nbins = None)    
            ax2.xaxis.set_major_locator(x_major_locator_zoom)
            ax2.xaxis.set_minor_locator(x_minor_locator_zoom)
        
            ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
            ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)        
                 
            ax1.set_xlabel('distance from coast (km)', fontsize = 12)
            ax1.set_ylabel('elevation (m asl.)', fontsize = 12)
            ax2.set_xlabel('distance from coast (km)', fontsize = 12)
            ax2.set_ylabel('elevation (m asl.)', fontsize = 12)    

            #   plot the colorbar
            cbar = plt.colorbar(im1, cax = ax3, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
            #cbar.ax.set_title(cbar_title, fontsize = 14, y  = 1.025)
            cbar.ax.set_ylabel(cbar_title, rotation = 90)
            cbar.ax.yaxis.set_label_position('left')
            cbar.ax.tick_params(labelsize = 12)    

            #   add the time stamp
            txt = ax6.text(0.0, 0.1, time_template % (time_key / 1000.), fontsize = 14, fontweight = 'bold')
            txt.set_clip_on(False)
            ax6.axis('off')
    
            #   create the graph of fresh pct in coastal and shelf zones
            x_time = np.linspace(0, tot_time / 1000., len(fresh_coast))
            lns1 = ax5.plot(x_time, fresh_coast, color = 'blue', linewidth = 2, label = '% fresh in coastal zone')
            lns2 = ax5.plot(x_time, fresh_shelf, color = 'green', linewidth = 2, label = '% fresh in shelf zone')
            ax5.set_xlim([0, round(time_sum, 0)])
            ax5.set_ylim([-5, 105])
            ax5.set_xlabel('Time since start of simulation (K years)', fontsize = 10)
            ax5.set_ylabel('% of fresh water', fontsize = 10)
            
            lns3 = ax5b.plot(x_time, lst_sea_lvl, color = 'cyan', linewidth = 3, label = 'Sea level (m asl)')
            ax5b.set_ylim([-125., 5.])
            ax5b.set_ylabel('Sea level (m asl)')

            #   add the legend for the concentration plots
            h,l = ax1.get_legend_handles_labels() # get labels and handles from ax1
            ax4.legend(h, l, ncol = 2, title = "Concentration plots legend", loc = 10)   
            ax4.axis('off')

            #   add the legend for the fresh % pct
            lns = lns1+ lns2 + lns3
            labs = [l.get_label() for l in lns]
            #ax.legend(lns, labs, loc=0)                
            
            #h2,l2 = ax5.get_legend_handles_labels() # get labels and handles from ax1
            ax7.legend(lns, labs, ncol = 3, title = "Fresh % graph legend", loc = 10)   
            ax7.axis('off')

            figname = 'img_%s.png' % (str(tot_time))

            canvas = FigureCanvas(fig)
            canvas.print_figure(os.path.join(img_dir, figname))

            #   add to the total time and the row in the csv file
            tot_time += nth_ts
            csv_row += int(nth_ts / 100)
            img_time_lst.append(tot_time / 1000.)

            #plt.savefig(os.path.join(img_dir, figname))
            #plt.close(fig)
                
            #else:
            #    key = key - 1
            #    tot_time -= nth_ts
            #    pass

    img_array = []
    os.chdir(img_dir)
    
    for i in range(len(img_time_lst) - 1):
        filename = 'img_' + str(int(img_time_lst[i] * 1000)) + '.png'
        print(filename)
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)    
    
    out = cv2.VideoWriter(os.path.join(vid_dir, '_sim_' + modelname + '.avi'), cv2.VideoWriter_fourcc(*'DIVX'), 5, size)
     
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()




"""

nth_ts = 1000
nc_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\coscat_1413_henry_other_param_combo_run\_summary\_nc'


nc_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\_WRAP_coscat_1413_henry_other_param_combo_run\_nc\_param_combo_run_14'
nth_ts = 1000
eq_lst = 

vid_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\coscat_1413_henry_other_param_combo_run\_summary'
img_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\coscat_1413_henry_other_param_combo_run\_summary\_img_to_vid'

"""



def avg_conc_final_video_from_nc(nc_dir, nth_ts, eq_lst, sea_lvl_lst, model_name, vid_dir, img_dir):

    #   create matplotlib objects necessary to export the video. The setup of the final video contains an overall concentration profile,
    #   a zoomed in profile into the current coastal zone (shifting with the sea level), a colorbar with the concentration values and also
    #   a graph that will show (together with written values) the % of fresh water in the coastal (as of 0m asl) and continental shelf 

    #   define the template for showing time
    time_template = 'Time (K years): = %.1f'

    #   load the concentration list and the first (starting) concentration profile
    tot_time = 0
    csv_row = 9
    time_sum = 0

    fresh_coast, fresh_shelf, time_lst, lst_sea_lvl = [0], [0], [0], [0.]

    #   get a sorted list of all the netcdf output files
    ts_lst = []
    for path, subdirs, files in os.walk(nc_dir):
        for name in files:
            print(os.path.join(path, name))  
            ts_lst.append(int(name.split('_')[-1].split('.')[0]))
    ts_lst = sorted(ts_lst)
            
    time_sum = ts_lst[-1] / 1000.
    
    #   get the original (starting) position of the coast and sea level
    st_sea_lvl = 0.0
    st_cst_dist = 0.0

    img_time_lst = [0]

    #   loop through the time steps and load each nth time step
    for key in range(0, len(ts_lst)):
        time_key = ts_lst[key]
        
        #   check if it is the nth time step
        if time_key % nth_ts == 0:
            print(time_key)

            #   open the netcdf file
            nc_name = os.path.join(vid_dir, '_nc', '_average_' + str(time_key) + '.nc')
            nc_in = xr.open_dataset(nc_name)

            #   check which EQ it is and also add to the total time
            cur_eq, cur_inl_pct, cur_shlf_pct = nc_in['eq'].values, nc_in['fresh_inl'].values, nc_in['fresh_shelf'].values
            cur_eq = np.array(cur_eq, ndmin = 1)[0]
            
            if 'char_time' not in cur_eq:
        
                #   check which EQ it is and also add to the total time                
                fresh_coast.append(round(np.array(cur_inl_pct, ndmin = 1)[0], 2))
                fresh_shelf.append(round(np.array(cur_shlf_pct, ndmin = 1)[0], 2))
                time_lst.append(tot_time)
                
                #   also get the current sea level based on the EQ
                cur_sea_lvl = sea_lvl_lst[eq_lst.index(cur_eq)]
                lst_sea_lvl.append(cur_sea_lvl)
                
                x_st = nc_in.coords['x'].values[0]
                x_end = nc_in.coords['x'].values[-1]
                bot = nc_in.coords['y'].values[-1]
                top = nc_in.coords['y'].values[0]
                x_cur_cst = np.array(nc_in['x_coast'].values, ndmin = 1)[0]
                x_st_zoom = x_cur_cst - 10.
                x_end_zoom = x_cur_cst + 10.
                bot_zoom = np.array(nc_in['bot_zoom'].values, ndmin = 1)[0]
                top_zoom = np.array(nc_in['top_zoom'].values, ndmin = 1)[0]
                
                #   get the concentration profile at given time
                conc_key = nc_in['solute concentration'].values
                conc_key[np.abs(conc_key) > 100.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0 
                if time_key == 0:
                    conc_key[np.abs(conc_key) == 0.0] = np.nan
                #conc_key[np.abs(conc_key) < 0.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0   
    
                #   define the constants for the colobrar
                cmap = plt.cm.jet
                cmaplist = [cmap(i) for i in range(cmap.N)]
                cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
                # define the bins and normalize
                bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]            
                norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
                cbar_title = 'Salinity (ppt)'
                
                fig = plt.Figure(figsize = (18, 12))
                ax1 = plt.subplot2grid((3, 3), (0, 1), colspan = 2, fig = fig)              #   overall concentration profiles
                ax2 = plt.subplot2grid((3, 3), (1, 0), colspan = 2, fig = fig)              #   zoomed in concentration profile
                ax3 = plt.subplot2grid((3, 3), (0, 0), fig = fig)                           #   color bar area 
                ax4 = plt.subplot2grid((3, 3), (2, 1), fig = fig)                           #   the legend area for the concentration plots
                ax5 = plt.subplot2grid((3, 3), (1, 2), fig = fig)              #   graph with % fresh area
                ax5b = ax5.twinx()                                                          #   instantiate a second axes that shares the same x-axis    
                ax6 = plt.subplot2grid((3, 3), (2, 0), fig = fig)                           #   area with the time stamp
                ax7 = plt.subplot2grid((3, 3), (2, 2), fig = fig)                           #   the legend area for fresh % graph
                #   specify the location of each of these figure parts
                ax1.set_position([0.175, 0.5, 0.8, 0.45])                        # [left, bottom, width, height]
                ax2.set_position([0.05, 0.15, 0.45, 0.3])
                ax3.set_position([0.05, 0.5, 0.025, 0.45])  
                ax4.set_position([0.2, 0.025, 0.275, 0.075])
                ax5.set_position([0.55, 0.15, 0.4, 0.3])  
                ax6.set_position([0.025, 0.05, 0.15, 0.075])   
                ax7.set_position([0.6, 0.025, 0.3, 0.075])
        
                #   draw the concentration profile for given time step
                im1 = ax1.imshow(conc_key[:, :], aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5, animated = True)
                ax1.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', label = 'shifting sea level', linestyle = '--')
                ax1.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black', label = 'shifting coastline')   
                #   also the zoomed in plot
                ax2.imshow(conc_key[:, :], aspect = 'auto', interpolation = 'none', cmap = cmap, norm = norm,\
                           extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5)
                
                #   plot the original sea level position and coastline position
                ax1.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey', label = 'current sea level', linestyle = '--')
                ax1.axvline(x = x_st - st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey', label = 'current coastline')   
         
                #   set the limits to the concentration plots
                ax1.set_xlim([x_st, x_end])
                ax1.set_ylim([bot, top])
                ax2.set_xlim([x_st_zoom, x_end_zoom])
                ax2.set_ylim([bot_zoom, top_zoom])
        
                #   set the gridlines and constant lines in the plot add constant lines with elevation = 0m asl. and coastline 
                ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
                ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
                ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
                ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
         
                ax1.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)
                ax2.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)    
               
                #   plot the shifting sea level and coastline position
                ax2.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', linestyle = '--')
                ax2.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black')    
  
                #   plot the original sea level position and coastline position
                ax2.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey')
                ax2.axvline(x = st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey')    
            
                x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 5.) * 5., math.ceil(x_end), 5.0), nbins = None)    
                x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 1.) * 1., math.ceil(x_end), 1.0), nbins = None)    
                ax1.xaxis.set_major_locator(x_major_locator)
                ax1.xaxis.set_minor_locator(x_minor_locator)
             
                y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 500.) * 500., math.ceil(top), 500.0), nbins = None)    
                y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 100.) * 100. , math.ceil(top), 100.0), nbins = None)    
                ax1.yaxis.set_major_locator(y_major_locator)
                ax1.yaxis.set_minor_locator(y_minor_locator)
                
                ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
                ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
            
                x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 5.) * 5., math.ceil(x_end_zoom), 5.0), nbins = None)    
                x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 1.) * 1., math.ceil(x_end_zoom), 1.0), nbins = None)    
                ax2.xaxis.set_major_locator(x_major_locator_zoom)
                ax2.xaxis.set_minor_locator(x_minor_locator_zoom)
            
                ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
                ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)        
                     
                ax1.set_xlabel('distance from coast (km)', fontsize = 12)
                ax1.set_ylabel('elevation (m asl.)', fontsize = 12)
                ax2.set_xlabel('distance from coast (km)', fontsize = 12)
                ax2.set_ylabel('elevation (m asl.)', fontsize = 12)    
    
                #   plot the colorbar
                cbar = plt.colorbar(im1, cax = ax3, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
                #cbar.ax.set_title(cbar_title, fontsize = 14, y  = 1.025)
                cbar.ax.set_ylabel(cbar_title, rotation = 90)
                cbar.ax.yaxis.set_label_position('left')
                cbar.ax.tick_params(labelsize = 12)    
    
                #   add the time stamp
                txt = ax6.text(0.0, 0.1, time_template % (time_key / 1000.), fontsize = 14, fontweight = 'bold')
                txt.set_clip_on(False)
                ax6.axis('off')
        
                #   create the graph of fresh pct in coastal and shelf zones
                x_time = np.linspace(0, tot_time / 1000., len(fresh_coast))
                lns1 = ax5.plot(x_time, fresh_coast, color = 'blue', linewidth = 2, label = '% fresh in coastal zone')
                lns2 = ax5.plot(x_time, fresh_shelf, color = 'green', linewidth = 2, label = '% fresh in shelf zone')
                ax5.set_xlim([0, round(time_sum, 0)])
                ax5.set_ylim([-5, 105])
                ax5.set_xlabel('Time since start of simulation (K years)', fontsize = 10)
                ax5.set_ylabel('% of fresh water', fontsize = 10)
                
                lns3 = ax5b.plot(x_time, lst_sea_lvl, color = 'cyan', linewidth = 3, label = 'Sea level (m asl)')
                ax5b.set_ylim([-125., 5.])
                ax5b.set_ylabel('Sea level (m asl)')

                #   add the legend for the concentration plots
                h,l = ax1.get_legend_handles_labels() # get labels and handles from ax1
                ax4.legend(h, l, ncol = 2, title = "Concentration plots legend", loc = 10)   
                ax4.axis('off')

                #   add the legend for the fresh % pct
                lns = lns1+ lns2 + lns3
                labs = [l.get_label() for l in lns]
                #ax.legend(lns, labs, loc=0)                
                
                #h2,l2 = ax5.get_legend_handles_labels() # get labels and handles from ax1
                ax7.legend(lns, labs, ncol = 3, title = "Fresh % graph legend", loc = 10)   
                ax7.axis('off')

                figname = 'img_%s.png' % (str(tot_time))

                canvas = FigureCanvas(fig)
                canvas.print_figure(os.path.join(img_dir, figname))

                #   add to the total time and the row in the csv file
                tot_time += nth_ts
                csv_row += int(nth_ts / 100)
                img_time_lst.append(tot_time / 1000.)

                #plt.savefig(os.path.join(img_dir, figname))
                #plt.close(fig)
                
            else:
                key = key - 1
                tot_time -= nth_ts
                pass

    img_array = []
    os.chdir(img_dir)
    
    for i in range(len(img_time_lst) - 1):
        filename = 'img_' + str(int(img_time_lst[i] * 1000)) + '.png'
        print(filename)
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)    
    
    out = cv2.VideoWriter(os.path.join(vid_dir, 'conc_through_time.avi'), cv2.VideoWriter_fourcc(*'DIVX'), 5, size)
     
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()






"""
    Function to create an averaged .netcdf file form all the individual files created as output by the models. 

in_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\coscat_1413_henry_other_param_combo_run'
out_dir_name = '_summary'
out_dir_name = 
ts_to_read = 500

eq_names = ['eq_00', 'eq_01', 'eq_02', 'eq_03', 'eq_04', 'eq_05', 'eq_06', 'eq_07', 'eq_08', 'eq_09',\
            'eq_10', 'eq_11', 'eq_12', 'eq_13', 'eq_14', 'eq_15', 'eq_16', 'eq_17', 'eq_18', 'eq_19',\
            'eq_20', 'eq_21', 'eq_21_char_time', 'eq_22', 'eq_23', 'eq_24', 'eq_25', 'eq_26', 'eq_27',\
            'eq_28', 'eq_29', 'eq_30', 'eq_31', 'eq_31_char_time']
sp_dur = [-1, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,\
          5000, 5000, 5000, 5000, 5000, 20000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 20000]


eq_names = ['eq_00', 'eq_01', 'eq_02', 'eq_03', 'eq_04', 'eq_05', 'eq_06', 'eq_07', 'eq_08', 'eq_09',\
          'eq_10', 'eq_11', 'eq_12', 'eq_13', 'eq_14', 'eq_15', 'eq_16', 'eq_17', 'eq_18', 'eq_19',\
          'eq_20', 'eq_21', 'eq_21_char_time', 'eq_22', 'eq_23', 'eq_24', 'eq_25', 'eq_26', 'eq_27',\
          'eq_28', 'eq_29', 'eq_30', 'eq_31', 'eq_31_char_time']
sp_dur = [-1, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,\
               5000, 5000, 5000, 5000, 5000, 20000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 20000]


#   create average netcdf file
in_dir = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\coscat_0016_delta_param_combo_run'
out_dir_name = r'g:\Water_Nexus\_A2\_models_param_runs\_1413\_summary_0016_delta_param_combo_run\avg_sim\_nc'
ts_to_read = 500

"""

def create_avg_netcdf(in_dir, out_dir_name, ts_to_read, eq_names, sp_dur):
    
    #   find all subfolders before creating the summary folder
    subfolders = next(os.walk(in_dir))[1]
    
    if '_summary' in subfolders:
        subfolders.remove('_summary')

    sum_dir = out_dir_name
        
    #   first check if the output directory exists, if not create it
    if not os.path.exists(out_dir_name):
        os.makedirs(sum_dir, exist_ok = True)        
    
    #if not os.path.exists(os.path.join(in_dir, out_dir_name)):
    #    sum_dir = os.path.join(in_dir, out_dir_name, '_nc')
    #    os.makedirs(sum_dir, exist_ok = True)

    #   set the time step variable
    total_time = 0
    
    #   define a function that opens a csv file and reads in the value of fresh water % for a specific time step
    def read_frsh_pct_from_csv(csv_dir, ts_val):
        #   open the csv file
        csv_file = pd.read_csv(csv_dir)                  
        frsh_inl_ts = csv_file.loc[csv_file[' time_end'] == 1000][' frsh_inl_pct'].item()
        frsh_shelf_ts = csv_file.loc[csv_file[' time_end'] == 1000][' frsh_shelf_pct'].item()     
        return frsh_inl_ts, frsh_shelf_ts

    #   loop through the equilibrium names and find the netcdf files for each time step
    for a in range(len(eq_names)):
        
        if 'char_time' not in eq_names[a]:
        
            total_time_sp = 0
            #   if the SP duration is equal to -1 it means that it can run for unlimited amount of time till it reaches equilibrium,
            #   some model runs will take less time to 'converge' if that is the case then just repeat the last profile for these model
            #   runs assuming it doesnt change anymore since steady state is achieved.
            #if sp_dur[a] == -1:

            #   create the directory name of the model path
            model_ws = os.path.join(in_dir, subfolders[0], eq_names[a])
            #   get the limits of the model domain
            dim_limits = get_model_info(model_ws, subfolders[0])
            x_st, bot_zoom, top_zoom = dim_limits[0], dim_limits[6], dim_limits[7]                
            #   set a counter that will check that there are still existing files with the given time_step
            file_exist = 1
            
            #   loop through an almost infinite loop and check if there are any existing files with the given time step
            for b in range(1000000):
            
                if file_exist > 0:
                    file_exist = 0
                    #   specify the name of the netcdf file to look for in the model folder directory
                    nc_name = eq_names[a] + '_' + str(total_time_sp) + '.nc'
                    print(nc_name)

                    conc_arr_lst, head_arr_lst, frsh_inl_lst, frsh_shelf_lst = [], [], [], []
           
                    #   create a name of the file for each of the subfolders, and check if it exists
                    for subf in subfolders:
                        nc_file = os.path.join(in_dir, subf, eq_names[a], '_nc', nc_name)
                        csv_file = os.path.join(in_dir, subf, eq_names[a], '_time_steps_DH_DC.csv')
            
                        if os.path.exists(nc_file):
                            print(nc_file)
                            file_exist += 1
                            
                            #   load the nc file and get the head and concentration arrays
                            last_nc = xr.open_dataset(nc_file)
                            conc_arr = np.expand_dims(np.array(last_nc['solute concentration'].values), axis = 1)
                            head_arr = np.expand_dims(np.array(last_nc['heads'].values), axis = 1)  
                            
                            frsh_pct = read_frsh_pct_from_csv(csv_file, total_time_sp)
                            x_lst = list(last_nc.coords['x'].values)
                            y_lst = list(last_nc.coords['y'].values)
                            
                        #   if not go into the subfolder and find the last .nc file (with the highest time step)
                        else:
                            eq_folder = os.path.join(in_dir, subf, eq_names[a], '_nc')
                            file_lst = os.listdir(eq_folder)
                            os.chdir(eq_folder)
                            newest = max(file_lst, key = lambda x: os.stat(x).st_mtime)
            
                            #   load the nc file and get the head and concentration arrays
                            last_nc = xr.open_dataset(newest)
                            conc_arr = np.expand_dims(np.array(last_nc['solute concentration'].values), axis = 1)
                            head_arr = np.expand_dims(np.array(last_nc['heads'].values), axis = 1)
                            x_lst = list(last_nc.coords['x'].values)
                            y_lst = list(last_nc.coords['y'].values)

                            frsh_pct = read_frsh_pct_from_csv(csv_file, int(newest.split('_')[-1].split('.')[0]))                                 
            
                        conc_arr_lst.append(conc_arr)
                        head_arr_lst.append(head_arr)
                        frsh_inl_lst.append(frsh_pct[0])
                        frsh_shelf_lst.append(frsh_pct[1])
                        
                    conc_arr_1 = conc_arr_lst[0]
                    for c in range(1, len(conc_arr_lst)):
                        conc_arr_1 = np.hstack((conc_arr_1, conc_arr_lst[c]))

                    head_arr_1 = head_arr_lst[0]
                    for c in range(1, len(head_arr_lst)):
                        head_arr_1 = np.hstack((head_arr_1, head_arr_lst[c]))
            
                    avg_conc = np.expand_dims(np.mean(conc_arr_1, axis = 1), axis = 1)
                    avg_head = np.expand_dims(np.mean(head_arr_1, axis = 1), axis = 1)
                    avg_frsh_inl = sum(frsh_inl_lst) / len(frsh_inl_lst)
                    avg_frsh_shelf = sum(frsh_shelf_lst) / len(frsh_shelf_lst)

                    #   for the concentration, heads and cbc create a netcdf file (to save memory)
                    xa_sum = xr.Dataset(data_vars = {'solute concentration' : (('y', 'x'), avg_conc[:, 0, :]),
                                                     'heads' : (('y', 'x'), avg_head[:, 0, :]),
                                                     'fresh_inl' : (avg_frsh_inl),
                                                     'fresh_shelf' : (avg_frsh_shelf),
                                                     'eq' : (eq_names[a]),
                                                     'x_coast' : ((x_lst[0] - 0.05) - x_st),
                                                     'bot_zoom' : (bot_zoom),
                                                     'top_zoom' : (top_zoom)},
                                        coords = {'x' : x_lst,
                                                  'y' : y_lst})
                    xa_sum = xa_sum.assign_coords(time = 0)
                    xa_name = '_average_' + str(total_time) + '.nc'
                    xa_sum.to_netcdf(os.path.join(sum_dir, xa_name))      
                
                    total_time_sp += ts_to_read         
                    total_time += ts_to_read
                

"""
all_frsh_csv_dir  = all_res_csv_dir
nth_th = 1000
eq_lst = equilibrium_names
eq_slr_lst = equilibrium_names_slr
avg_nc_dir = simu_name
img_dir = avg_img_dir
out_vid_dir = avg_dir
"""   
            
def plot_avg_conc_profile_SLR_scs(all_frsh_csv_dir, nth_ts, eq_lst, eq_slr_lst, avg_nc_dir, img_dir, out_vid_dir):
    
    #   a zoomed in profile into the current coastal zone (shifting with the sea level), a colorbar with the concentration values and also
    #   a graph that will show (together with written values) the % of fresh water in the coastal (as of 0m asl) and continental shelf 
    sea_levels = [-130., -129., -125., -126., -124., -110., -95., -75., -67., -57,\
                  -45, -32., -17., -7., -3, -2., -1.5, -0.7, -0.2, 0.,\
                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0.,\
                  0., 0., 0., 0., 0., 0., 0., 0.]
    slr_names = ['RCP_26_', 'RCP_45_', 'RCP_85_']
    sea_levels_slr = [[0.4, 0.7, 0.8, 0.8, 0.8], [0.5, 1., 1.5, 1.5, 1.5], [0.7, 2.1, 3.7, 3.7, 3.7]]
    
    #   first check if the output directory exists, if not create it
    os.makedirs(img_dir, exist_ok = True)        

    # make directories for the 9 combinations of DEM and RCP
    os.makedirs(os.path.join(img_dir, 'merit'), exist_ok = True)        
    os.makedirs(os.path.join(img_dir, 'coastal'), exist_ok = True)             
    os.makedirs(os.path.join(img_dir, 'gebco'), exist_ok = True)       

    #   select the different DEM and SLR scenarios
    sc_lst = ['merit', 'coastal', 'gebco']

    #   loop through the list
    for i in range(len(sc_lst)):
        #   define the current dem-rcp combo scenario list and its name
        sc_name = sc_lst[i]
        dir_out = os.path.join(img_dir, sc_name)

        for d in range(len(slr_names)):
            rcp = slr_names[d]
            sea_lvl_lst = sea_levels + sea_levels_slr[d]
            
            #   define the template for showing time
            time_template = 'Time (K years): = %.1f'
        
            #   load the concentration list and the first (starting) concentration profile
            tot_time = 0
            csv_row = 9
            time_sum = 0
        
            #fresh_coast, fresh_shelf, time_lst, lst_sea_lvl = [0], [0], [0], [0.]
            #inl_min_lst, inl_max_lst, shlf_min_lst, shlf_max_lst = [0], [0], [0], [0]
    
            fresh_coast, fresh_shelf, time_lst, lst_sea_lvl = [], [], [], []
            inl_min_lst, inl_max_lst, shlf_min_lst, shlf_max_lst = [], [], [], []
        
            #   get a sorted list of all the netcdf output files
            ts_lst = []
            for path, subdirs, files in os.walk(os.path.join(avg_nc_dir, sc_name)):
                for name in files:
                    print(os.path.join(path, name))  
                    ts_lst.append(int(name.split('_')[-1].split('.')[0]))
            ts_lst = sorted(ts_lst)
                    
            time_sum = ts_lst[-1] / 1000.
            
            #   get the original (starting) position of the coast and sea level
            st_sea_lvl = 0.0
            st_cst_dist = 0.0
        
            img_time_lst = []
        
            df = pd.read_csv(all_frsh_csv_dir, index_col = False, header = 0)
            df_hdr = []
            for a in df.columns:
                df_hdr.append(a)
        
            #   loop through the time steps and load each nth time step
            for key in range(len(ts_lst)):
                time_key = ts_lst[key]
                x_time = [a/1000 for a in ts_lst[: key+1]] #np.linspace(0, time_key / 1000., len(fresh_coast))
                print(time_key)
    
                """
                if key > 30000:
                    #   open the netcdf file
                    nc_name = os.path.join(avg_nc_dir, sc_name, '_average_' + rcp + str(time_key) + '.nc')                
                else:
                    #   open the netcdf file
                    nc_name = os.path.join(avg_nc_dir, sc_name, '_average_' + str(time_key) + '.nc')
                nc_in = xr.open_dataset(nc_name)
                """
                
                #   check which EQ it is and also add to the total time
                if time_key == 0:
                    nc_name = os.path.join(avg_nc_dir, sc_name, '_average_500.nc')
                    nc_in = xr.open_dataset(nc_name)                
                    cur_eq, cur_inl_pct, cur_shlf_pct = nc_in['eq'].values, nc_in['fresh_inl'].values, nc_in['fresh_shelf'].values
                    cur_eq = np.array(cur_eq, ndmin = 1)[0]
                    #   also get the current sea level based on the EQ
                    cur_sea_lvl = sea_levels[eq_lst.index(cur_eq)]
                    lst_sea_lvl.append(cur_sea_lvl)                  
                    
                elif time_key > 30000:
                    nc_name = os.path.join(avg_nc_dir, sc_name, '_average_' + rcp + str(time_key) + '.nc')  
                    nc_in = xr.open_dataset(nc_name)                
                    cur_eq, cur_inl_pct, cur_shlf_pct = nc_in['eq'].values, nc_in['fresh_inl'].values, nc_in['fresh_shelf'].values
                    cur_eq = np.array(cur_eq, ndmin = 1)[0]
                    #   also get the current sea level based on the EQ
                    cur_sea_lvl = sea_lvl_lst[eq_slr_lst.index('AP' + cur_eq.split('AP')[-1])]
                    lst_sea_lvl.append(cur_sea_lvl)      
                    
                else:
                    nc_name = os.path.join(avg_nc_dir, sc_name, '_average_' + str(time_key) + '.nc')  
                    nc_in = xr.open_dataset(nc_name)   
                    cur_eq, cur_inl_pct, cur_shlf_pct = nc_in['eq'].values, nc_in['fresh_inl'].values, nc_in['fresh_shelf'].values
                    cur_eq = np.array(cur_eq, ndmin = 1)[0]
                    cur_sea_lvl = sea_levels[eq_lst.index(cur_eq)]
                    lst_sea_lvl.append(cur_sea_lvl)  
                    
                #   check which EQ it is and also add to the total time                
                fresh_coast.append(round(np.array(cur_inl_pct, ndmin = 1)[0], 2))
                fresh_shelf.append(round(np.array(cur_shlf_pct, ndmin = 1)[0], 2))
                time_lst.append(tot_time)
                
                #   read in the csv file and find all the % of fresh water for the corresponding time step
                inl_lst, shlf_lst = [], []
                ts_vals = df[(df['time_end'] == time_key)]
                try:
                    for hdr in df_hdr:
                        if 'inl' in hdr:
                            inl_lst.append(ts_vals[hdr].values[0])
                        elif 'shlf' in hdr:
                            shlf_lst.append(ts_vals[hdr].values[0])
                #   happens for ts = 0 years (initial state)
                except IndexError:
                    for hdr in df_hdr:
                        if 'inl' in hdr:
                            inl_lst.append(0.)
                        elif 'shlf' in hdr:
                            shlf_lst.append(0.)
                    
                #   select a min and max values for each of the inland and shelf areas
                inl_min_lst.append(min(inl_lst))
                inl_max_lst.append(max(inl_lst))
                shlf_min_lst.append(min(shlf_lst))
                shlf_max_lst.append(max(shlf_lst))
                
                x_st = nc_in.coords['x'].values[0]
                x_end = nc_in.coords['x'].values[-1]
                bot = nc_in.coords['y'].values[-1]
                top = nc_in.coords['y'].values[0]
                x_cur_cst = np.array(nc_in['x_coast'].values, ndmin = 1)[0]
                x_st_zoom = max(-10., x_st)#x_cur_cst - 10.
                x_end_zoom = min(10., x_end) # x_cur_cst + 10.
                top_zoom = np.array(nc_in['top_zoom'].values, ndmin = 1)[0]            
                x_coord_sel_end = min(nc_in['x'].values.tolist(), key=lambda x:abs(x - x_end_zoom))
                bot_zoom = top_zoom - nc_in.sel(x = x_coord_sel_end)['heads'].values.tolist().index(next(x for x in nc_in.sel(x = x_coord_sel_end)['heads'].values.tolist()[::-1] if x > -999.)) * 10 - 10
                #bot_zoom = np.array(nc_in['bot_zoom'].values, ndmin = 1)[0]
                
                #   get the concentration profile at given time
                conc_key = nc_in['solute concentration'].values
                conc_key[np.abs(conc_key) > 100.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0 
                if time_key == 0:
                    conc_key[np.abs(conc_key) == 0.0] = np.nan
                #conc_key[np.abs(conc_key) < 0.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0   
    
                #   define the constants for the colobrar
                cmap = plt.cm.jet
                cmaplist = [cmap(i) for i in range(cmap.N)]
                cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
                # define the bins and normalize
                bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]            
                norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
                cbar_title = 'Salinity (ppt)'
                
                fig = plt.Figure(figsize = (18, 12))
                ax1 = plt.subplot2grid((3, 3), (0, 1), colspan = 2, fig = fig)              #   overall concentration profiles
                ax2 = plt.subplot2grid((3, 3), (1, 0), colspan = 2, fig = fig)              #   zoomed in concentration profile
                ax3 = plt.subplot2grid((3, 3), (0, 0), fig = fig)                           #   color bar area 
                ax4 = plt.subplot2grid((3, 3), (2, 1), fig = fig)                           #   the legend area for the concentration plots
                ax5 = plt.subplot2grid((3, 3), (1, 2), fig = fig)              #   graph with % fresh area
                ax5b = ax5.twinx()                                                          #   instantiate a second axes that shares the same x-axis    
                ax6 = plt.subplot2grid((3, 3), (2, 0), fig = fig)                           #   area with the time stamp
                ax7 = plt.subplot2grid((3, 3), (2, 2), fig = fig)                           #   the legend area for fresh % graph
                #   specify the location of each of these figure parts
                ax1.set_position([0.175, 0.5, 0.8, 0.45])                        # [left, bottom, width, height]
                ax2.set_position([0.05, 0.15, 0.45, 0.3])
                ax3.set_position([0.05, 0.5, 0.025, 0.45])  
                ax4.set_position([0.2, 0.025, 0.275, 0.075])
                ax5.set_position([0.55, 0.15, 0.4, 0.3])  
                ax6.set_position([0.025, 0.05, 0.15, 0.075])   
                ax7.set_position([0.6, 0.025, 0.3, 0.075])
        
                #   draw the concentration profile for given time step
                im1 = ax1.imshow(conc_key[:, :], aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                           extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5, animated = True)
                ax1.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', label = 'shifting sea level', linestyle = '--')
                ax1.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black', label = 'shifting coastline')   
                #   also the zoomed in plot
                ax2.imshow(conc_key[:, :], aspect = 'auto', interpolation = 'none', cmap = cmap, norm = norm,\
                           extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5)
                
                #   plot the original sea level position and coastline position
                ax1.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey', label = 'current sea level', linestyle = '--')
                ax1.axvline(x = x_st - st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey', label = 'current coastline')   
         
                #   set the limits to the concentration plots
                ax1.set_xlim([x_st, x_end])
                ax1.set_ylim([bot, top])
                ax2.set_xlim([x_st_zoom, x_end_zoom])
                ax2.set_ylim([bot_zoom, top_zoom])
        
                #   set the gridlines and constant lines in the plot add constant lines with elevation = 0m asl. and coastline 
                ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
                ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
                ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
                ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
         
                ax1.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)
                ax2.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)    
               
                #   plot the shifting sea level and coastline position
                ax2.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', linestyle = '--')
                ax2.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black')    
      
                #   plot the original sea level position and coastline position
                ax2.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey')
                ax2.axvline(x = x_st - st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey')    
            
                x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 5.) * 5., math.ceil(x_end), 5.0), nbins = None)    
                x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 1.) * 1., math.ceil(x_end), 1.0), nbins = None)    
                ax1.xaxis.set_major_locator(x_major_locator)
                ax1.xaxis.set_minor_locator(x_minor_locator)
             
                y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 500.) * 500., math.ceil(top), 500.0), nbins = None)    
                y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 100.) * 100. , math.ceil(top), 100.0), nbins = None)    
                ax1.yaxis.set_major_locator(y_major_locator)
                ax1.yaxis.set_minor_locator(y_minor_locator)
                
                ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
                ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
            
                x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 5.) * 5., math.ceil(x_end_zoom), 5.0), nbins = None)    
                x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 1.) * 1., math.ceil(x_end_zoom), 1.0), nbins = None)    
                ax2.xaxis.set_major_locator(x_major_locator_zoom)
                ax2.xaxis.set_minor_locator(x_minor_locator_zoom)
            
                ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
                ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)        
                     
                ax1.set_xlabel('distance from coast (km)', fontsize = 12)
                ax1.set_ylabel('elevation (m asl.)', fontsize = 12)
                ax2.set_xlabel('distance from coast (km)', fontsize = 12)
                ax2.set_ylabel('elevation (m asl.)', fontsize = 12)    
    
                #   plot the colorbar
                cbar = plt.colorbar(im1, cax = ax3, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
                #cbar.ax.set_title(cbar_title, fontsize = 14, y  = 1.025)
                cbar.ax.set_ylabel(cbar_title, rotation = 90)
                cbar.ax.yaxis.set_label_position('left')
                cbar.ax.tick_params(labelsize = 12)    
    
                #   add the time stamp
                txt = ax6.text(0.0, 0.1, time_template % (time_key / 1000.), fontsize = 14, fontweight = 'bold')
                txt.set_clip_on(False)
                ax6.axis('off')
        
                #   create the graph of fresh pct in coastal and shelf zones
    
                lns1 = ax5.plot(x_time, fresh_coast, color = 'blue', linewidth = 2, label = '% fresh in coastal zone')
                lns2 = ax5.plot(x_time, fresh_shelf, color = 'green', linewidth = 2, label = '% fresh in shelf zone')
                
                ax5.fill_between(x_time, inl_min_lst, inl_max_lst, facecolor='blue', alpha = 0.5)
                ax5.fill_between(x_time, shlf_min_lst, shlf_max_lst, facecolor='green', alpha = 0.5)
                
                ax5.set_xlim([0, round(time_sum, 0)])
                ax5.set_ylim([-5, 105])
                ax5.set_xlabel('Time since start of simulation (K years)', fontsize = 10)
                ax5.set_ylabel('% of fresh water', fontsize = 10)
                
                lns3 = ax5b.plot(x_time, lst_sea_lvl, color = 'cyan', linewidth = 3, label = 'Sea level (m asl)')
                ax5b.set_ylim([-135., 5.])
                ax5b.set_ylabel('Sea level (m asl)')
    
                #   add the legend for the concentration plots
                h,l = ax1.get_legend_handles_labels() # get labels and handles from ax1
                ax4.legend(h, l, ncol = 2, title = "Concentration plots legend", loc = 10)   
                ax4.axis('off')
    
                #   add the legend for the fresh % pct
                lns = lns1+ lns2 + lns3
                labs = [l.get_label() for l in lns]
                #ax.legend(lns, labs, loc=0)                
                
                #h2,l2 = ax5.get_legend_handles_labels() # get labels and handles from ax1
                ax7.legend(lns, labs, ncol = 3, title = "Fresh % graph legend", loc = 10)   
                ax7.axis('off')
    
                figname = '%simg_%s.png' % (rcp, str(time_key))
    
                canvas = FigureCanvas(fig)
                canvas.print_figure(os.path.join(dir_out, figname))
    
                #   add to the total time and the row in the csv file
                img_time_lst.append(time_key / 1000.)
    
                #plt.savefig(os.path.join(img_dir, figname))
                #plt.close(fig)
            
            img_array = []
            os.chdir(dir_out)
            
            for i in range(len(img_time_lst) - 1):
                filename = rcp + 'img_' + str(int(img_time_lst[i] * 1000)) + '.png'
                print(filename)
                img = cv2.imread(filename)
                height, width, layers = img.shape
                size = (width,height)
                img_array.append(img)    
            
            out = cv2.VideoWriter(os.path.join(out_vid_dir, sc_name + '_' + rcp + 'AVG_conc_through_time.avi'), cv2.VideoWriter_fourcc(*'DIVX'), 5, size)
             
            for i in range(len(img_array)):
                out.write(img_array[i])
            out.release()





def plot_avg_conc_profile(all_frsh_csv_dir, nth_ts, eq_lst, avg_nc_dir, img_dir, out_vid_dir):
    
    #   a zoomed in profile into the current coastal zone (shifting with the sea level), a colorbar with the concentration values and also
    #   a graph that will show (together with written values) the % of fresh water in the coastal (as of 0m asl) and continental shelf 
    sea_levels_26 = [-130., -129., -125., -126., -124., -110., -95., -75., -67., -57,\
                  -45, -32., -17., -7., -3, -2., -1.5, -0.7, -0.2, 0.,\
                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0.,\
                  0., 0., 0., 0., 0., 0., 0., 0., 0.4, 0.7, 0.8, 0.8, 0.8]
    
    sea_levels_45 = [-130., -129., -125., -126., -124., -110., -95., -75., -67., -57,\
                  -45, -32., -17., -7., -3, -2., -1.5, -0.7, -0.2, 0.,\
                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0.,\
                  0., 0., 0., 0., 0., 0., 0., 0., 0.5, 1., 1.5, 1.5, 1.5]
    sea_levels_85 = [-130., -129., -125., -126., -124., -110., -95., -75., -67., -57,\
                  -45, -32., -17., -7., -3, -2., -1.5, -0.7, -0.2, 0.,\
                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0.,\
                  0., 0., 0., 0., 0., 0., 0., 0., 0.7, 2.1, 3.7, 3.7, 3.7]

    #   first check if the output directory exists, if not create it
    os.makedirs(img_dir, exist_ok = True)        

    # make directories for the 9 combinations of DEM and RCP
    os.makedirs(os.path.join(img_dir, 'merit_26'), exist_ok = True)       
    os.makedirs(os.path.join(img_dir, 'merit_45'), exist_ok = True)       
    os.makedirs(os.path.join(img_dir, 'merit_85'), exist_ok = True)       
    os.makedirs(os.path.join(img_dir, 'coastal_26'), exist_ok = True)       
    os.makedirs(os.path.join(img_dir, 'coastal_45'), exist_ok = True)       
    os.makedirs(os.path.join(img_dir, 'coastal_85'), exist_ok = True)       
    os.makedirs(os.path.join(img_dir, 'gebco_26'), exist_ok = True)       
    os.makedirs(os.path.join(img_dir, 'gebco_45'), exist_ok = True)       
    os.makedirs(os.path.join(img_dir, 'gebco_85'), exist_ok = True)       

    #   select the different DEM and SLR scenarios
    sc_lst = ['merit_26', 'merit_45', 'merit_85', 'coastal_26',\
              'coastal_45', 'coastal_85', 'gebco_26', 'gebco_45', 'gebco_85']

    #   loop through the list
    for i in range(len(sc_lst)):
        #   define the current dem-rcp combo scenario list and its name
        sc_name = sc_lst[i]
        rcp = sc_name.split('_')[-1]
        dir_out = os.path.join(img_dir, sc_name)

        if rcp == '26':
            sea_lvl_lst = sea_levels_26
        elif rcp == '45':
            sea_lvl_lst = sea_levels_45
        elif rcp == '85':
            sea_lvl_lst = sea_levels_85
            
        #   define the template for showing time
        time_template = 'Time (K years): = %.1f'
    
        #   load the concentration list and the first (starting) concentration profile
        tot_time = 0
        csv_row = 9
        time_sum = 0
    
        #fresh_coast, fresh_shelf, time_lst, lst_sea_lvl = [0], [0], [0], [0.]
        #inl_min_lst, inl_max_lst, shlf_min_lst, shlf_max_lst = [0], [0], [0], [0]

        fresh_coast, fresh_shelf, time_lst, lst_sea_lvl = [], [], [], []
        inl_min_lst, inl_max_lst, shlf_min_lst, shlf_max_lst = [], [], [], []
    
        #   get a sorted list of all the netcdf output files
        ts_lst = []
        for path, subdirs, files in os.walk(os.path.join(avg_nc_dir, sc_name)):
            for name in files:
                print(os.path.join(path, name))  
                ts_lst.append(int(name.split('_')[-1].split('.')[0]))
        ts_lst = sorted(ts_lst)
                
        time_sum = ts_lst[-1] / 1000.
        
        #   get the original (starting) position of the coast and sea level
        st_sea_lvl = 0.0
        st_cst_dist = 0.0
    
        img_time_lst = []
    
        df = pd.read_csv(all_frsh_csv_dir, index_col = False, header = 0)
        df_hdr = []
        for a in df.columns:
            df_hdr.append(a)
    
        #   loop through the time steps and load each nth time step
        for key in range(len(ts_lst)):
            time_key = ts_lst[key]
            x_time = [a/1000 for a in ts_lst[: key+1]] #np.linspace(0, time_key / 1000., len(fresh_coast))
            print(time_key)

            #   open the netcdf file
            nc_name = os.path.join(avg_nc_dir, sc_name, '_average_' + str(time_key) + '.nc')
            nc_in = xr.open_dataset(nc_name)

            #   check which EQ it is and also add to the total time
            if time_key == 0:
                nc_name = os.path.join(avg_nc_dir, sc_name, '_average_500.nc')
                nc_in = xr.open_dataset(nc_name)                
                cur_eq, cur_inl_pct, cur_shlf_pct = nc_in['eq'].values, nc_in['fresh_inl'].values, nc_in['fresh_shelf'].values
                cur_eq = np.array(cur_eq, ndmin = 1)[0]
            else:
                cur_eq, cur_inl_pct, cur_shlf_pct = nc_in['eq'].values, nc_in['fresh_inl'].values, nc_in['fresh_shelf'].values
                cur_eq = np.array(cur_eq, ndmin = 1)[0]
            
            #   check which EQ it is and also add to the total time                
            fresh_coast.append(round(np.array(cur_inl_pct, ndmin = 1)[0], 2))
            fresh_shelf.append(round(np.array(cur_shlf_pct, ndmin = 1)[0], 2))
            time_lst.append(tot_time)
            
            #   read in the csv file and find all the % of fresh water for the corresponding time step
            inl_lst, shlf_lst = [], []
            ts_vals = df[(df['time_end'] == time_key)]
            try:
                for hdr in df_hdr:
                    if 'inl' in hdr:
                        inl_lst.append(ts_vals[hdr].values[0])
                    elif 'shlf' in hdr:
                        shlf_lst.append(ts_vals[hdr].values[0])
            #   happens for ts = 0 years (initial state)
            except IndexError:
                for hdr in df_hdr:
                    if 'inl' in hdr:
                        inl_lst.append(0.)
                    elif 'shlf' in hdr:
                        shlf_lst.append(0.)
                
            #   select a min and max values for each of the inland and shelf areas
            inl_min_lst.append(min(inl_lst))
            inl_max_lst.append(max(inl_lst))
            shlf_min_lst.append(min(shlf_lst))
            shlf_max_lst.append(max(shlf_lst))
            
            #   also get the current sea level based on the EQ
            cur_sea_lvl = sea_lvl_lst[eq_lst.index(cur_eq)]
            lst_sea_lvl.append(cur_sea_lvl)
            
            x_st = nc_in.coords['x'].values[0]
            x_end = nc_in.coords['x'].values[-1]
            bot = nc_in.coords['y'].values[-1]
            top = nc_in.coords['y'].values[0]
            x_cur_cst = np.array(nc_in['x_coast'].values, ndmin = 1)[0]
            x_st_zoom = max(-10., x_st)#x_cur_cst - 10.
            x_end_zoom = min(10., x_end) # x_cur_cst + 10.
            top_zoom = np.array(nc_in['top_zoom'].values, ndmin = 1)[0]            
            x_coord_sel_end = min(nc_in['x'].values.tolist(), key=lambda x:abs(x - x_end_zoom))
            bot_zoom = top_zoom - nc_in.sel(x = x_coord_sel_end)['heads'].values.tolist().index(next(x for x in nc_in.sel(x = x_coord_sel_end)['heads'].values.tolist()[::-1] if x > -999.)) * 10 - 10
            #bot_zoom = np.array(nc_in['bot_zoom'].values, ndmin = 1)[0]
            
            #   get the concentration profile at given time
            conc_key = nc_in['solute concentration'].values
            conc_key[np.abs(conc_key) > 100.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0 
            if time_key == 0:
                conc_key[np.abs(conc_key) == 0.0] = np.nan
            #conc_key[np.abs(conc_key) < 0.0] = np.nan  #   35.5 because sometimes the concentrations are > 35.0   

            #   define the constants for the colobrar
            cmap = plt.cm.jet
            cmaplist = [cmap(i) for i in range(cmap.N)]
            cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
            # define the bins and normalize
            bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]            
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
            cbar_title = 'Salinity (ppt)'
            
            fig = plt.Figure(figsize = (18, 12))
            ax1 = plt.subplot2grid((3, 3), (0, 1), colspan = 2, fig = fig)              #   overall concentration profiles
            ax2 = plt.subplot2grid((3, 3), (1, 0), colspan = 2, fig = fig)              #   zoomed in concentration profile
            ax3 = plt.subplot2grid((3, 3), (0, 0), fig = fig)                           #   color bar area 
            ax4 = plt.subplot2grid((3, 3), (2, 1), fig = fig)                           #   the legend area for the concentration plots
            ax5 = plt.subplot2grid((3, 3), (1, 2), fig = fig)              #   graph with % fresh area
            ax5b = ax5.twinx()                                                          #   instantiate a second axes that shares the same x-axis    
            ax6 = plt.subplot2grid((3, 3), (2, 0), fig = fig)                           #   area with the time stamp
            ax7 = plt.subplot2grid((3, 3), (2, 2), fig = fig)                           #   the legend area for fresh % graph
            #   specify the location of each of these figure parts
            ax1.set_position([0.175, 0.5, 0.8, 0.45])                        # [left, bottom, width, height]
            ax2.set_position([0.05, 0.15, 0.45, 0.3])
            ax3.set_position([0.05, 0.5, 0.025, 0.45])  
            ax4.set_position([0.2, 0.025, 0.275, 0.075])
            ax5.set_position([0.55, 0.15, 0.4, 0.3])  
            ax6.set_position([0.025, 0.05, 0.15, 0.075])   
            ax7.set_position([0.6, 0.025, 0.3, 0.075])
    
            #   draw the concentration profile for given time step
            im1 = ax1.imshow(conc_key[:, :], aspect='auto', interpolation='none', cmap = cmap, norm = norm,\
                       extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5, animated = True)
            ax1.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', label = 'shifting sea level', linestyle = '--')
            ax1.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black', label = 'shifting coastline')   
            #   also the zoomed in plot
            ax2.imshow(conc_key[:, :], aspect = 'auto', interpolation = 'none', cmap = cmap, norm = norm,\
                       extent = (x_st, x_end, bot, top), vmin = 0, vmax = 35.5)
            
            #   plot the original sea level position and coastline position
            ax1.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey', label = 'current sea level', linestyle = '--')
            ax1.axvline(x = x_st - st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey', label = 'current coastline')   
     
            #   set the limits to the concentration plots
            ax1.set_xlim([x_st, x_end])
            ax1.set_ylim([bot, top])
            ax2.set_xlim([x_st_zoom, x_end_zoom])
            ax2.set_ylim([bot_zoom, top_zoom])
    
            #   set the gridlines and constant lines in the plot add constant lines with elevation = 0m asl. and coastline 
            ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
            ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
            ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
            ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
     
            ax1.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)
            ax2.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)    
           
            #   plot the shifting sea level and coastline position
            ax2.axhline(y = cur_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'black', linestyle = '--')
            ax2.axvline(x = x_cur_cst, ymin = bot, ymax = top, lw = 3., color = 'black')    
  
            #   plot the original sea level position and coastline position
            ax2.axhline(y = st_sea_lvl, xmin = x_st, xmax = x_end, lw = 3., color = 'grey')
            ax2.axvline(x = x_st - st_cst_dist, ymin = bot, ymax = top, lw = 3., color = 'grey')    
        
            x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 5.) * 5., math.ceil(x_end), 5.0), nbins = None)    
            x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st / 1.) * 1., math.ceil(x_end), 1.0), nbins = None)    
            ax1.xaxis.set_major_locator(x_major_locator)
            ax1.xaxis.set_minor_locator(x_minor_locator)
         
            y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 500.) * 500., math.ceil(top), 500.0), nbins = None)    
            y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(bot / 100.) * 100. , math.ceil(top), 100.0), nbins = None)    
            ax1.yaxis.set_major_locator(y_major_locator)
            ax1.yaxis.set_minor_locator(y_minor_locator)
            
            ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
            ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
        
            x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 5.) * 5., math.ceil(x_end_zoom), 5.0), nbins = None)    
            x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_st_zoom / 1.) * 1., math.ceil(x_end_zoom), 1.0), nbins = None)    
            ax2.xaxis.set_major_locator(x_major_locator_zoom)
            ax2.xaxis.set_minor_locator(x_minor_locator_zoom)
        
            ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
            ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)        
                 
            ax1.set_xlabel('distance from coast (km)', fontsize = 12)
            ax1.set_ylabel('elevation (m asl.)', fontsize = 12)
            ax2.set_xlabel('distance from coast (km)', fontsize = 12)
            ax2.set_ylabel('elevation (m asl.)', fontsize = 12)    

            #   plot the colorbar
            cbar = plt.colorbar(im1, cax = ax3, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
            #cbar.ax.set_title(cbar_title, fontsize = 14, y  = 1.025)
            cbar.ax.set_ylabel(cbar_title, rotation = 90)
            cbar.ax.yaxis.set_label_position('left')
            cbar.ax.tick_params(labelsize = 12)    

            #   add the time stamp
            txt = ax6.text(0.0, 0.1, time_template % (time_key / 1000.), fontsize = 14, fontweight = 'bold')
            txt.set_clip_on(False)
            ax6.axis('off')
    
            #   create the graph of fresh pct in coastal and shelf zones

            lns1 = ax5.plot(x_time, fresh_coast, color = 'blue', linewidth = 2, label = '% fresh in coastal zone')
            lns2 = ax5.plot(x_time, fresh_shelf, color = 'green', linewidth = 2, label = '% fresh in shelf zone')
            
            ax5.fill_between(x_time, inl_min_lst, inl_max_lst, facecolor='blue', alpha = 0.5)
            ax5.fill_between(x_time, shlf_min_lst, shlf_max_lst, facecolor='green', alpha = 0.5)
            
            ax5.set_xlim([0, round(time_sum, 0)])
            ax5.set_ylim([-5, 105])
            ax5.set_xlabel('Time since start of simulation (K years)', fontsize = 10)
            ax5.set_ylabel('% of fresh water', fontsize = 10)
            
            lns3 = ax5b.plot(x_time, lst_sea_lvl, color = 'cyan', linewidth = 3, label = 'Sea level (m asl)')
            ax5b.set_ylim([-135., 5.])
            ax5b.set_ylabel('Sea level (m asl)')

            #   add the legend for the concentration plots
            h,l = ax1.get_legend_handles_labels() # get labels and handles from ax1
            ax4.legend(h, l, ncol = 2, title = "Concentration plots legend", loc = 10)   
            ax4.axis('off')

            #   add the legend for the fresh % pct
            lns = lns1+ lns2 + lns3
            labs = [l.get_label() for l in lns]
            #ax.legend(lns, labs, loc=0)                
            
            #h2,l2 = ax5.get_legend_handles_labels() # get labels and handles from ax1
            ax7.legend(lns, labs, ncol = 3, title = "Fresh % graph legend", loc = 10)   
            ax7.axis('off')

            figname = 'img_%s.png' % (str(time_key))

            canvas = FigureCanvas(fig)
            canvas.print_figure(os.path.join(dir_out, figname))

            #   add to the total time and the row in the csv file
            img_time_lst.append(time_key / 1000.)

            #plt.savefig(os.path.join(img_dir, figname))
            #plt.close(fig)
        
        img_array = []
        os.chdir(dir_out)
        
        for i in range(len(img_time_lst) - 1):
            filename = 'img_' + str(int(img_time_lst[i] * 1000)) + '.png'
            print(filename)
            img = cv2.imread(filename)
            height, width, layers = img.shape
            size = (width,height)
            img_array.append(img)    
        
        out = cv2.VideoWriter(os.path.join(out_vid_dir, sc_name + '_AVG_conc_through_time.avi'), cv2.VideoWriter_fourcc(*'DIVX'), 5, size)
         
        for i in range(len(img_array)):
            out.write(img_array[i])
        out.release()

"""
sum_dir = summary_dir
simu_name = model_name_basis
eq_lst = equilibrium_names
eq_slr_lst = equilibrium_names_slr
"""
def create_fin_coscat_run_csv_files_SLR_scs(coscat_id_str, main_dir, sum_dir, eq_lst, eq_slr_lst, sp_duration, simu_name):

    #   define the summary dir and the model simulation dir
    #sum_dir = os.path.join(main_dir, '_summary_' + coscat_id_str + '_' + cst_type + simu_name)
    #model_simu_dir = os.path.join(main_dir, 'coscat_' + coscat_id_str + '_' + cst_type + '_' + simu_name)

    #   define the output name of the final csv file
    frsh_vol_csv = os.path.join(sum_dir, '_frsh_vol_summary_' + 'COSCAT_' + simu_name + '.csv')
    tot_summary_inl_csv = os.path.join(sum_dir, '_frsh_INL_' + 'COSCAT_' + simu_name + '.csv')
    tot_summary_shlf_csv = os.path.join(sum_dir, '_frsh_SHLF_' + 'COSCAT_' + simu_name + '.csv')
    
    #   find all csv files in the folder 
    eq_0_ts_lst = []
    csv_dir = os.path.join(sum_dir, '_csv_results')
    #print(csv_dir)
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            #print(filename)
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            #simu_nr_int = int(simu_nr_str)        
            #   get all the time steps for eq_0
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)
            col_names = list(df.columns.values)
            eq_0_rows = df.loc[df['EQ'] == 'BP_30000_to_20000']

            #   get the last value of the time_end column, only if the model simulation converged, if not append a -1
            try:
                if len(df.loc[df['EQ'] == 'BP_00100_to_00000'][' time_end'].values) == 2:
                    eq_0_last_ts = eq_0_rows.tail(1).values[0][2] 
                    eq_0_ts_lst.append([eq_0_last_ts, simu_nr_str])
            except IndexError:
                #print('The csv file is empty (most probably)')
                continue

    #   get the index of the tuple with the maximum time step
    max_idx = [i for i, tupl in enumerate(eq_0_ts_lst) if tupl[0] == max([x[0] for x in eq_0_ts_lst])][0]
    max_ts_0 = eq_0_ts_lst[max_idx][0]

    #   column containing the time steps and EQ numbers, taken form the simulation with the longest EQ_0
    ts_csv_dir = os.path.join(csv_dir, 'sc_' + eq_0_ts_lst[max_idx][1] + '_summary_EQs_output.csv')
    ts_csv = pd.read_csv(ts_csv_dir, index_col = False, header = 0)
    ts_csv_eq_col = list(ts_csv['EQ'].values)
    ts_csv_ts_col = list(ts_csv[' time_end'].values)
    
    #   create a list of headers 
    frsh_headers = ['eq']
    for a in range(len(eq_0_ts_lst)):
        frsh_headers.append(eq_0_ts_lst[a][1] + '_frsh_inl')
        frsh_headers.append(eq_0_ts_lst[a][1] + '_frsh_shlf')
    frsh_headers.append('avg_frsh_inl')
    frsh_headers.append('avg_frsh_shlf')
    
    df_frsh_pct = pd.DataFrame(columns = frsh_headers)
    dict_frsh_pct = {}
    dict_frsh_pct['eq'] = ts_csv_eq_col
    frsh_inl_arr = np.zeros((len(eq_0_ts_lst), len(ts_csv_eq_col)))
    frsh_shlf_arr = np.zeros((len(eq_0_ts_lst), len(ts_csv_eq_col)))
    idx_lay = 0
    
    frsh_inl_all_lst, frsh_shlf_all_lst, headers_lst = [], [], []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)
            
            if len(df.loc[df['EQ'] == 'BP_00100_to_00000'][' time_end'].values) == 2:
                if len(df.loc[df['EQ'] == 'BP_30000_to_20000'][' time_end'].values) == 20:
                    inl_lst = df[' frsh_inl_pct'].values.tolist()
                    shlf_lst = df[' frsh_shelf_pct'].values.tolist()
                    
                else:
                    inl_lst_eq_0 = df.loc[df['EQ'] == 'BP_30000_to_20000'][' frsh_inl_pct'].values.tolist()
                    shlf_lst_eq_0 = df.loc[df['EQ'] == 'BP_30000_to_20000'][' frsh_shelf_pct'].values.tolist()
                    for f in range(20 - len(inl_lst_eq_0)):
                        inl_lst_eq_0.append(inl_lst_eq_0[-1])
                        shlf_lst_eq_0.append(shlf_lst_eq_0[-1])
                    inl_lst_rest = df.loc[df['EQ'] != 'BP_30000_to_20000'][' frsh_inl_pct'].values.tolist()
                    shlf_lst_rest = df.loc[df['EQ'] != 'BP_30000_to_20000'][' frsh_shelf_pct'].values.tolist()                    
                    
                    inl_lst = inl_lst_eq_0 + inl_lst_rest
                    shlf_lst = shlf_lst_eq_0 + shlf_lst_rest
                    
                #   set the column headers
                hd_frsh_inl = simu_nr_str + '_frsh_inl'
                hd_frsh_shlf = simu_nr_str + '_frsh_shlf'
                
                inl_lst.insert(0, inl_lst[0])
                shlf_lst.insert(0, shlf_lst[0])
                
                headers_lst.append([hd_frsh_inl, hd_frsh_shlf])                            
                frsh_inl_all_lst.append(inl_lst)
                frsh_shlf_all_lst.append(shlf_lst)
        
    inl_frsh_arr = np.array(frsh_inl_all_lst)
    shlf_frsh_arr = np.array(frsh_shlf_all_lst)        
        
    avg_inl_frsh_lst = list(np.mean(inl_frsh_arr, axis = 0))
    avg_shlf_frsh_lst = list(np.mean(shlf_frsh_arr, axis = 0))    
    
    avg_inl_frsh_lst = [round(i, 2) for i in avg_inl_frsh_lst]
    avg_shlf_frsh_lst = [round(i, 2) for i in avg_shlf_frsh_lst]

    eq_lst_to_dict, ts_lst_to_dict = [], []     
    for g in range(int(10000 / 500) + 1):
        eq_lst_to_dict.append('BP_30000_to_20000')
        ts_lst_to_dict.append(500. * g)
    time = ts_lst_to_dict[-1] + 500
    for h in range(1, len(eq_lst)):
        eq_name = eq_lst[h]
        ts_dur = sp_duration[h]    
        ts_step = ts_dur / 2
        for i in range(1, int(ts_dur / (ts_dur / 2) + 1)):
            eq_lst_to_dict.append(eq_name)
            ts_lst_to_dict.append(time)        
            time += ts_step
        if eq_name == 'BP_03000_to_02000':
            time -= 450

    ts_lst_to_dict = np.arange(0, 28500, 500).tolist() + np.arange(28050, 30050, 50).tolist()
    #ts_lst_to_dict = sorted(ts_lst_to_dict)[1:]
    
    ts_last = ts_lst_to_dict[-1]
    for slr_sc in ['RCP_26_', 'RCP_45_', 'RCP_85_']:
        eq_lst_to_dict = eq_lst_to_dict + [slr_sc + i for i in eq_slr_lst]
        ts_lst_to_dict = ts_lst_to_dict + np.linspace(ts_last + 50, ts_last + 500, 10).tolist()

    dict_frsh_pct['eq'] = eq_lst_to_dict
    dict_frsh_pct['time_end'] = ts_lst_to_dict                
            
    #   fill the dictionary
    for f in range(len(headers_lst)):
        dict_frsh_pct[headers_lst[f][0]] = frsh_inl_all_lst[f]
        dict_frsh_pct[headers_lst[f][1]] = frsh_shlf_all_lst[f]
    
    dict_frsh_pct['avg_frsh_inl'] = avg_inl_frsh_lst
    dict_frsh_pct['avg_frsh_shlf'] = avg_shlf_frsh_lst    
        
    df_out = pd.DataFrame.from_dict(dict_frsh_pct, orient = 'index').T   
    df_out.to_csv(frsh_vol_csv)

    """         create the summary csv files, with % fresh at certain times and the convergence times         """
    #geo_csv_dir = os.path.join(sum_dir, '_geo_summary_' + coscat_id_str + '.csv')
    #geo_df = pd.read_csv(geo_csv_dir, index_col = False, header = 0)    
    #geo_headers = geo_df.columns.values
    
    """
    sum_headers = ['coscat_id', 'simu_id', 'eq_0_t_end', 'eq_0_t_end_inl_frsh', 'eq_0_t_end_shlf_frsh',\
                   'eq_21_t_end_inl_frsh', 'eq_21_t_end_shlf_frsh', 'eq_31_t_end_inl_frsh', 'eq_31_t_end_shlf_frsh',\
                   'eq_21_CT_t_end', 'eq_21_CT_t_end_inl_frsh', 'eq_21_CT_t_end_shlf_frsh', 'eq_31_CT_t_end', 'eq_31_CT_t_end_inl_frsh', 'eq_31_CT_t_end_shlf_frsh']
    """
    
    sum_headers_inl = ['coscat_id', 'srm_id', 'sc_id', 'eq_0_t_end', 'eq_20000_BP_inl_frsh', 'eq_19000_BP_inl_frsh',\
                       'eq_18000_BP_inl_frsh', 'eq_17000_BP_inl_frsh','eq_16000_BP_inl_frsh', 'eq_15000_BP_inl_frsh',\
                       'eq_14000_BP_inl_frsh', 'eq_13000_BP_inl_frsh','eq_12000_BP_inl_frsh', 'eq_11000_BP_inl_frsh',\
                       'eq_10000_BP_inl_frsh', 'eq_09000_BP_inl_frsh','eq_08000_BP_inl_frsh', 'eq_07000_BP_inl_frsh',\
                       'eq_06000_BP_inl_frsh', 'eq_05000_BP_inl_frsh','eq_04000_BP_inl_frsh', 'eq_03000_BP_inl_frsh',\
                       'eq_02000_BP_inl_frsh', 'eq_01900_BP_inl_frsh','eq_01800_BP_inl_frsh', 'eq_01700_BP_inl_frsh',\
                       'eq_01600_BP_inl_frsh', 'eq_01500_BP_inl_frsh','eq_01400_BP_inl_frsh', 'eq_01300_BP_inl_frsh',\
                       'eq_01200_BP_inl_frsh', 'eq_01100_BP_inl_frsh','eq_01000_BP_inl_frsh', 'eq_00900_BP_inl_frsh',\
                       'eq_00800_BP_inl_frsh', 'eq_00700_BP_inl_frsh','eq_00600_BP_inl_frsh', 'eq_00500_BP_inl_frsh',\
                       'eq_00400_BP_inl_frsh', 'eq_00300_BP_inl_frsh','eq_00200_BP_inl_frsh', 'eq_00100_BP_inl_frsh',\
                       'eq_00000_BP_inl_frsh', 'RCP_26_050_AP_inl_frsh', 'RCP_26_100_AP_inl_frsh', 'RCP_26_150_AP_inl_frsh',\
                       'RCP_26_200_AP_inl_frsh', 'RCP_26_250_AP_inl_frsh', 'RCP_26_300_AP_inl_frsh', 'RCP_26_350_AP_inl_frsh',\
                       'RCP_26_400_AP_inl_frsh', 'RCP_26_450_AP_inl_frsh', 'RCP_26_500_AP_inl_frsh', 'RCP_45_050_AP_inl_frsh',\
                       'RCP_45_100_AP_inl_frsh', 'RCP_45_150_AP_inl_frsh', 'RCP_45_200_AP_inl_frsh', 'RCP_45_250_AP_inl_frsh',\
                       'RCP_45_300_AP_inl_frsh', 'RCP_45_350_AP_inl_frsh', 'RCP_45_400_AP_inl_frsh', 'RCP_45_450_AP_inl_frsh', 'RCP_45_500_AP_inl_frsh',\
                       'RCP_85_050_AP_inl_frsh', 'RCP_85_100_AP_inl_frsh', 'RCP_85_150_AP_inl_frsh', 'RCP_85_200_AP_inl_frsh', 'RCP_85_250_AP_inl_frsh',\
                       'RCP_85_300_AP_inl_frsh', 'RCP_85_350_AP_inl_frsh', 'RCP_85_400_AP_inl_frsh', 'RCP_85_450_AP_inl_frsh', 'RCP_85_500_AP_inl_frsh']                      
                       
    all_sc_lst = []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            #print(simu_nr_str)
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)    
            if len(df.loc[df['EQ'] == 'BP_00100_to_00000'][' time_end'].values) == 2:
                try:                      
                    sc_lst = []
                    #   read the specific rows
                    eq_0_t_end = df[(df['EQ'] == eq_lst[0])][' time_end'].values[-1]
                    eq_20000_BP_inl_frsh = df[(df['EQ'] == eq_lst[0])][' frsh_inl_pct'].values[-1]
                    eq_19000_BP_inl_frsh = df[(df['EQ'] == eq_lst[1])][' frsh_inl_pct'].values[-1]                
                    eq_18000_BP_inl_frsh = df[(df['EQ'] == eq_lst[2])][' frsh_inl_pct'].values[-1]
                    eq_17000_BP_inl_frsh = df[(df['EQ'] == eq_lst[3])][' frsh_inl_pct'].values[-1]     
                    eq_16000_BP_inl_frsh = df[(df['EQ'] == eq_lst[4])][' frsh_inl_pct'].values[-1]     
                    eq_15000_BP_inl_frsh = df[(df['EQ'] == eq_lst[5])][' frsh_inl_pct'].values[-1]                 
                    eq_14000_BP_inl_frsh = df[(df['EQ'] == eq_lst[6])][' frsh_inl_pct'].values[-1]
                    eq_13000_BP_inl_frsh = df[(df['EQ'] == eq_lst[7])][' frsh_inl_pct'].values[-1]     
                    eq_12000_BP_inl_frsh = df[(df['EQ'] == eq_lst[8])][' frsh_inl_pct'].values[-1]     
                    eq_11000_BP_inl_frsh = df[(df['EQ'] == eq_lst[9])][' frsh_inl_pct'].values[-1]                    
                    eq_10000_BP_inl_frsh = df[(df['EQ'] == eq_lst[10])][' frsh_inl_pct'].values[-1]
                    eq_09000_BP_inl_frsh = df[(df['EQ'] == eq_lst[11])][' frsh_inl_pct'].values[-1]     
                    eq_08000_BP_inl_frsh = df[(df['EQ'] == eq_lst[12])][' frsh_inl_pct'].values[-1]     
                    eq_07000_BP_inl_frsh = df[(df['EQ'] == eq_lst[13])][' frsh_inl_pct'].values[-1]                            
                    eq_06000_BP_inl_frsh = df[(df['EQ'] == eq_lst[14])][' frsh_inl_pct'].values[-1]
                    eq_05000_BP_inl_frsh = df[(df['EQ'] == eq_lst[15])][' frsh_inl_pct'].values[-1]     
                    eq_04000_BP_inl_frsh = df[(df['EQ'] == eq_lst[16])][' frsh_inl_pct'].values[-1]     
                    eq_03000_BP_inl_frsh = df[(df['EQ'] == eq_lst[17])][' frsh_inl_pct'].values[-1]                
                    eq_02000_BP_inl_frsh = df[(df['EQ'] == eq_lst[18])][' frsh_inl_pct'].values[-1]   
                    eq_01900_BP_inl_frsh = df[(df['EQ'] == eq_lst[19])][' frsh_inl_pct'].values[-1]     
                    eq_01800_BP_inl_frsh = df[(df['EQ'] == eq_lst[20])][' frsh_inl_pct'].values[-1]     
                    eq_01700_BP_inl_frsh = df[(df['EQ'] == eq_lst[21])][' frsh_inl_pct'].values[-1]                            
                    eq_01600_BP_inl_frsh = df[(df['EQ'] == eq_lst[22])][' frsh_inl_pct'].values[-1]
                    eq_01500_BP_inl_frsh = df[(df['EQ'] == eq_lst[23])][' frsh_inl_pct'].values[-1]     
                    eq_01400_BP_inl_frsh = df[(df['EQ'] == eq_lst[24])][' frsh_inl_pct'].values[-1]     
                    eq_01300_BP_inl_frsh = df[(df['EQ'] == eq_lst[25])][' frsh_inl_pct'].values[-1]                
                    eq_01200_BP_inl_frsh = df[(df['EQ'] == eq_lst[26])][' frsh_inl_pct'].values[-1]   
                    eq_01100_BP_inl_frsh = df[(df['EQ'] == eq_lst[27])][' frsh_inl_pct'].values[-1]                
                    eq_01000_BP_inl_frsh = df[(df['EQ'] == eq_lst[28])][' frsh_inl_pct'].values[-1]   
                    eq_00900_BP_inl_frsh = df[(df['EQ'] == eq_lst[29])][' frsh_inl_pct'].values[-1]     
                    eq_00800_BP_inl_frsh = df[(df['EQ'] == eq_lst[30])][' frsh_inl_pct'].values[-1]     
                    eq_00700_BP_inl_frsh = df[(df['EQ'] == eq_lst[31])][' frsh_inl_pct'].values[-1]                            
                    eq_00600_BP_inl_frsh = df[(df['EQ'] == eq_lst[32])][' frsh_inl_pct'].values[-1]
                    eq_00500_BP_inl_frsh = df[(df['EQ'] == eq_lst[33])][' frsh_inl_pct'].values[-1]     
                    eq_00400_BP_inl_frsh = df[(df['EQ'] == eq_lst[34])][' frsh_inl_pct'].values[-1]     
                    eq_00300_BP_inl_frsh = df[(df['EQ'] == eq_lst[35])][' frsh_inl_pct'].values[-1]                
                    eq_00200_BP_inl_frsh = df[(df['EQ'] == eq_lst[36])][' frsh_inl_pct'].values[-1]   
                    eq_00100_BP_inl_frsh = df[(df['EQ'] == eq_lst[37])][' frsh_inl_pct'].values[-1]                
                    eq_00000_BP_inl_frsh = df[(df['EQ'] == eq_lst[38])][' frsh_inl_pct'].values[-1]
                    
                    RCP_26_050_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[0])][' frsh_inl_pct'].values[0]     
                    RCP_26_100_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[0])][' frsh_inl_pct'].values[-1]     
                    RCP_26_150_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[1])][' frsh_inl_pct'].values[0]                
                    RCP_26_200_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[1])][' frsh_inl_pct'].values[-1]   
                    RCP_26_250_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[2])][' frsh_inl_pct'].values[0]                     
                    RCP_26_300_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[2])][' frsh_inl_pct'].values[-1]     
                    RCP_26_350_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[3])][' frsh_inl_pct'].values[0]     
                    RCP_26_400_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[3])][' frsh_inl_pct'].values[-1]                
                    RCP_26_450_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[4])][' frsh_inl_pct'].values[0]   
                    RCP_26_500_AP_inl_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[4])][' frsh_inl_pct'].values[-1]                              
                    
                    RCP_45_050_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[0])][' frsh_inl_pct'].values[0]     
                    RCP_45_100_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[0])][' frsh_inl_pct'].values[-1]     
                    RCP_45_150_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[1])][' frsh_inl_pct'].values[0]                
                    RCP_45_200_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[1])][' frsh_inl_pct'].values[-1]   
                    RCP_45_250_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[2])][' frsh_inl_pct'].values[0]                     
                    RCP_45_300_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[2])][' frsh_inl_pct'].values[-1]     
                    RCP_45_350_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[3])][' frsh_inl_pct'].values[0]     
                    RCP_45_400_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[3])][' frsh_inl_pct'].values[-1]                
                    RCP_45_450_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[4])][' frsh_inl_pct'].values[0]   
                    RCP_45_500_AP_inl_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[4])][' frsh_inl_pct'].values[-1]                            
                    
                    RCP_85_050_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[0])][' frsh_inl_pct'].values[0]     
                    RCP_85_100_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[0])][' frsh_inl_pct'].values[-1]     
                    RCP_85_150_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[1])][' frsh_inl_pct'].values[0]                
                    RCP_85_200_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[1])][' frsh_inl_pct'].values[-1]   
                    RCP_85_250_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[2])][' frsh_inl_pct'].values[0]                     
                    RCP_85_300_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[2])][' frsh_inl_pct'].values[-1]     
                    RCP_85_350_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[3])][' frsh_inl_pct'].values[0]     
                    RCP_85_400_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[3])][' frsh_inl_pct'].values[-1]                
                    RCP_85_450_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[4])][' frsh_inl_pct'].values[0]   
                    RCP_85_500_AP_inl_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[4])][' frsh_inl_pct'].values[-1]                         
                    
                    sc_lst.append([coscat_id_str, simu_name.split('SRM_')[-1], simu_nr_str, eq_0_t_end, eq_20000_BP_inl_frsh,\
                                   eq_19000_BP_inl_frsh, eq_18000_BP_inl_frsh, eq_17000_BP_inl_frsh, eq_16000_BP_inl_frsh,\
                                   eq_15000_BP_inl_frsh, eq_14000_BP_inl_frsh, eq_13000_BP_inl_frsh, eq_12000_BP_inl_frsh,\
                                   eq_11000_BP_inl_frsh, eq_10000_BP_inl_frsh, eq_09000_BP_inl_frsh, eq_08000_BP_inl_frsh,\
                                   eq_07000_BP_inl_frsh, eq_06000_BP_inl_frsh, eq_05000_BP_inl_frsh, eq_04000_BP_inl_frsh,\
                                   eq_03000_BP_inl_frsh, eq_02000_BP_inl_frsh, eq_01900_BP_inl_frsh, eq_01800_BP_inl_frsh,\
                                   eq_01700_BP_inl_frsh, eq_01600_BP_inl_frsh, eq_01500_BP_inl_frsh, eq_01400_BP_inl_frsh,\
                                   eq_01300_BP_inl_frsh, eq_01200_BP_inl_frsh, eq_01100_BP_inl_frsh, eq_01000_BP_inl_frsh,\
                                   eq_00900_BP_inl_frsh, eq_00800_BP_inl_frsh, eq_00700_BP_inl_frsh, eq_00600_BP_inl_frsh,\
                                   eq_00500_BP_inl_frsh, eq_00400_BP_inl_frsh, eq_00300_BP_inl_frsh, eq_00200_BP_inl_frsh,\
                                   eq_00100_BP_inl_frsh, eq_00000_BP_inl_frsh,\
                                   RCP_26_050_AP_inl_frsh, RCP_26_100_AP_inl_frsh, RCP_26_150_AP_inl_frsh, RCP_26_200_AP_inl_frsh,\
                                   RCP_26_250_AP_inl_frsh, RCP_26_300_AP_inl_frsh, RCP_26_350_AP_inl_frsh, RCP_26_400_AP_inl_frsh,\
                                   RCP_26_450_AP_inl_frsh, RCP_26_500_AP_inl_frsh,\
                                   RCP_45_050_AP_inl_frsh, RCP_45_100_AP_inl_frsh, RCP_45_150_AP_inl_frsh, RCP_45_200_AP_inl_frsh,\
                                   RCP_45_250_AP_inl_frsh, RCP_45_300_AP_inl_frsh, RCP_45_350_AP_inl_frsh, RCP_45_400_AP_inl_frsh,\
                                   RCP_45_450_AP_inl_frsh, RCP_45_500_AP_inl_frsh,\
                                   RCP_85_050_AP_inl_frsh, RCP_85_100_AP_inl_frsh, RCP_85_150_AP_inl_frsh, RCP_85_200_AP_inl_frsh,\
                                   RCP_85_250_AP_inl_frsh, RCP_85_300_AP_inl_frsh, RCP_85_350_AP_inl_frsh, RCP_85_400_AP_inl_frsh,\
                                   RCP_85_450_AP_inl_frsh, RCP_85_500_AP_inl_frsh])
            
                    all_sc_lst.append(sc_lst)
                except FileNotFoundError:
                    #print('Model simulation did not terminate')
                    continue

    all_sc_lst = [i[0] for i in all_sc_lst]
    df_sum = pd.DataFrame(all_sc_lst, columns = sum_headers_inl)
    df_sum.to_csv(tot_summary_inl_csv)            

    sum_headers_shlf = ['coscat_id', 'srm_id', 'sc_id', 'eq_0_t_end', 'eq_20000_BP_shlf_frsh', 'eq_19000_BP_shlf_frsh',\
                       'eq_18000_BP_shlf_frsh', 'eq_17000_BP_shlf_frsh','eq_16000_BP_shlf_frsh', 'eq_15000_BP_shlf_frsh',\
                       'eq_14000_BP_shlf_frsh', 'eq_13000_BP_shlf_frsh','eq_12000_BP_shlf_frsh', 'eq_11000_BP_shlf_frsh',\
                       'eq_10000_BP_shlf_frsh', 'eq_09000_BP_shlf_frsh','eq_08000_BP_shlf_frsh', 'eq_07000_BP_shlf_frsh',\
                       'eq_06000_BP_shlf_frsh', 'eq_05000_BP_shlf_frsh','eq_04000_BP_shlf_frsh', 'eq_03000_BP_shlf_frsh',\
                       'eq_02000_BP_shlf_frsh', 'eq_01900_BP_shlf_frsh','eq_01800_BP_shlf_frsh', 'eq_01700_BP_shlf_frsh',\
                       'eq_01600_BP_shlf_frsh', 'eq_01500_BP_shlf_frsh','eq_01400_BP_shlf_frsh', 'eq_01300_BP_shlf_frsh',\
                       'eq_01200_BP_shlf_frsh', 'eq_01100_BP_shlf_frsh','eq_01000_BP_shlf_frsh', 'eq_00900_BP_shlf_frsh',\
                       'eq_00800_BP_shlf_frsh', 'eq_00700_BP_shlf_frsh','eq_00600_BP_shlf_frsh', 'eq_00500_BP_shlf_frsh',\
                       'eq_00400_BP_shlf_frsh', 'eq_00300_BP_shlf_frsh','eq_00200_BP_shlf_frsh', 'eq_00100_BP_shlf_frsh',\
                       'eq_00000_BP_shlf_frsh', 'RCP_26_050_AP_shlf_frsh', 'RCP_26_100_AP_shlf_frsh', 'RCP_26_150_AP_shlf_frsh',\
                       'RCP_26_200_AP_shlf_frsh', 'RCP_26_250_AP_shlf_frsh', 'RCP_26_300_AP_shlf_frsh', 'RCP_26_350_AP_shlf_frsh',\
                       'RCP_26_400_AP_shlf_frsh', 'RCP_26_450_AP_shlf_frsh', 'RCP_26_500_AP_shlf_frsh', 'RCP_45_050_AP_shlf_frsh',\
                       'RCP_45_100_AP_shlf_frsh', 'RCP_45_150_AP_shlf_frsh', 'RCP_45_200_AP_shlf_frsh', 'RCP_45_250_AP_shlf_frsh',\
                       'RCP_45_300_AP_shlf_frsh', 'RCP_45_350_AP_shlf_frsh', 'RCP_45_400_AP_shlf_frsh', 'RCP_45_450_AP_shlf_frsh', 'RCP_45_500_AP_shlf_frsh',\
                       'RCP_85_050_AP_shlf_frsh', 'RCP_85_100_AP_shlf_frsh', 'RCP_85_150_AP_shlf_frsh', 'RCP_85_200_AP_shlf_frsh', 'RCP_85_250_AP_shlf_frsh',\
                       'RCP_85_300_AP_shlf_frsh', 'RCP_85_350_AP_shlf_frsh', 'RCP_85_400_AP_shlf_frsh', 'RCP_85_450_AP_shlf_frsh', 'RCP_85_500_AP_shlf_frsh']   

    all_sc_lst = []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            #print(simu_nr_str)
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)    
            if len(df.loc[df['EQ'] == 'BP_00100_to_00000'][' time_end'].values) == 2:
                try:                      
                    sc_lst = []
                    #   read the specific rows
                    eq_0_t_end = df[(df['EQ'] == eq_lst[0])][' time_end'].values[-1]
                    eq_20000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[0])][' frsh_shelf_pct'].values[-1]
                    eq_19000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[1])][' frsh_shelf_pct'].values[-1]                
                    eq_18000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[2])][' frsh_shelf_pct'].values[-1]
                    eq_17000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[3])][' frsh_shelf_pct'].values[-1]     
                    eq_16000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[4])][' frsh_shelf_pct'].values[-1]     
                    eq_15000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[5])][' frsh_shelf_pct'].values[-1]                 
                    eq_14000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[6])][' frsh_shelf_pct'].values[-1]
                    eq_13000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[7])][' frsh_shelf_pct'].values[-1]     
                    eq_12000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[8])][' frsh_shelf_pct'].values[-1]     
                    eq_11000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[9])][' frsh_shelf_pct'].values[-1]                    
                    eq_10000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[10])][' frsh_shelf_pct'].values[-1]
                    eq_09000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[11])][' frsh_shelf_pct'].values[-1]     
                    eq_08000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[12])][' frsh_shelf_pct'].values[-1]     
                    eq_07000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[13])][' frsh_shelf_pct'].values[-1]                            
                    eq_06000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[14])][' frsh_shelf_pct'].values[-1]
                    eq_05000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[15])][' frsh_shelf_pct'].values[-1]     
                    eq_04000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[16])][' frsh_shelf_pct'].values[-1]     
                    eq_03000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[17])][' frsh_shelf_pct'].values[-1]                
                    eq_02000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[18])][' frsh_shelf_pct'].values[-1]   
                    eq_01900_BP_shlf_frsh = df[(df['EQ'] == eq_lst[19])][' frsh_shelf_pct'].values[-1]     
                    eq_01800_BP_shlf_frsh = df[(df['EQ'] == eq_lst[20])][' frsh_shelf_pct'].values[-1]     
                    eq_01700_BP_shlf_frsh = df[(df['EQ'] == eq_lst[21])][' frsh_shelf_pct'].values[-1]                            
                    eq_01600_BP_shlf_frsh = df[(df['EQ'] == eq_lst[22])][' frsh_shelf_pct'].values[-1]
                    eq_01500_BP_shlf_frsh = df[(df['EQ'] == eq_lst[23])][' frsh_shelf_pct'].values[-1]     
                    eq_01400_BP_shlf_frsh = df[(df['EQ'] == eq_lst[24])][' frsh_shelf_pct'].values[-1]     
                    eq_01300_BP_shlf_frsh = df[(df['EQ'] == eq_lst[25])][' frsh_shelf_pct'].values[-1]                
                    eq_01200_BP_shlf_frsh = df[(df['EQ'] == eq_lst[26])][' frsh_shelf_pct'].values[-1]   
                    eq_01100_BP_shlf_frsh = df[(df['EQ'] == eq_lst[27])][' frsh_shelf_pct'].values[-1]                
                    eq_01000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[28])][' frsh_shelf_pct'].values[-1]   
                    eq_00900_BP_shlf_frsh = df[(df['EQ'] == eq_lst[29])][' frsh_shelf_pct'].values[-1]     
                    eq_00800_BP_shlf_frsh = df[(df['EQ'] == eq_lst[30])][' frsh_shelf_pct'].values[-1]     
                    eq_00700_BP_shlf_frsh = df[(df['EQ'] == eq_lst[31])][' frsh_shelf_pct'].values[-1]                            
                    eq_00600_BP_shlf_frsh = df[(df['EQ'] == eq_lst[32])][' frsh_shelf_pct'].values[-1]
                    eq_00500_BP_shlf_frsh = df[(df['EQ'] == eq_lst[33])][' frsh_shelf_pct'].values[-1]     
                    eq_00400_BP_shlf_frsh = df[(df['EQ'] == eq_lst[34])][' frsh_shelf_pct'].values[-1]     
                    eq_00300_BP_shlf_frsh = df[(df['EQ'] == eq_lst[35])][' frsh_shelf_pct'].values[-1]                
                    eq_00200_BP_shlf_frsh = df[(df['EQ'] == eq_lst[36])][' frsh_shelf_pct'].values[-1]   
                    eq_00100_BP_shlf_frsh = df[(df['EQ'] == eq_lst[37])][' frsh_shelf_pct'].values[-1]                
                    eq_00000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[38])][' frsh_shelf_pct'].values[-1]   

                    RCP_26_050_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[0]     
                    RCP_26_100_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_150_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[0]                
                    RCP_26_200_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[-1]   
                    RCP_26_250_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[0]                     
                    RCP_26_300_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_350_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[0]     
                    RCP_26_400_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[-1]                
                    RCP_26_450_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[0]   
                    RCP_26_500_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[-1]                              
                    
                    RCP_45_050_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[0]     
                    RCP_45_100_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_150_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[0]                
                    RCP_45_200_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[-1]   
                    RCP_45_250_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[0]                     
                    RCP_45_300_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_350_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[0]     
                    RCP_45_400_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[-1]                
                    RCP_45_450_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[0]   
                    RCP_45_500_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[-1]                            
                    
                    RCP_85_050_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[0]     
                    RCP_85_100_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_150_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[0]                
                    RCP_85_200_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[-1]   
                    RCP_85_250_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[0]                     
                    RCP_85_300_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_350_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[0]     
                    RCP_85_400_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[-1]                
                    RCP_85_450_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[0]   
                    RCP_85_500_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[-1]                         
                    
                    sc_lst.append([coscat_id_str, simu_name.split('SRM_')[-1], simu_nr_str, eq_0_t_end, eq_20000_BP_shlf_frsh,\
                                   eq_19000_BP_shlf_frsh, eq_18000_BP_shlf_frsh, eq_17000_BP_shlf_frsh, eq_16000_BP_shlf_frsh,\
                                   eq_15000_BP_shlf_frsh, eq_14000_BP_shlf_frsh, eq_13000_BP_shlf_frsh, eq_12000_BP_shlf_frsh,\
                                   eq_11000_BP_shlf_frsh, eq_10000_BP_shlf_frsh, eq_09000_BP_shlf_frsh, eq_08000_BP_shlf_frsh,\
                                   eq_07000_BP_shlf_frsh, eq_06000_BP_shlf_frsh, eq_05000_BP_shlf_frsh, eq_04000_BP_shlf_frsh,\
                                   eq_03000_BP_shlf_frsh, eq_02000_BP_shlf_frsh, eq_01900_BP_shlf_frsh, eq_01800_BP_shlf_frsh,\
                                   eq_01700_BP_shlf_frsh, eq_01600_BP_shlf_frsh, eq_01500_BP_shlf_frsh, eq_01400_BP_shlf_frsh,\
                                   eq_01300_BP_shlf_frsh, eq_01200_BP_shlf_frsh, eq_01100_BP_shlf_frsh, eq_01000_BP_shlf_frsh,\
                                   eq_00900_BP_shlf_frsh, eq_00800_BP_shlf_frsh, eq_00700_BP_shlf_frsh, eq_00600_BP_shlf_frsh,\
                                   eq_00500_BP_shlf_frsh, eq_00400_BP_shlf_frsh, eq_00300_BP_shlf_frsh, eq_00200_BP_shlf_frsh,\
                                   eq_00100_BP_shlf_frsh, eq_00000_BP_shlf_frsh,\
                                   RCP_26_050_AP_shlf_frsh, RCP_26_100_AP_shlf_frsh, RCP_26_150_AP_shlf_frsh, RCP_26_200_AP_shlf_frsh,\
                                   RCP_26_250_AP_shlf_frsh, RCP_26_300_AP_shlf_frsh, RCP_26_350_AP_shlf_frsh, RCP_26_400_AP_shlf_frsh,\
                                   RCP_26_450_AP_shlf_frsh, RCP_26_500_AP_shlf_frsh,\
                                   RCP_45_050_AP_shlf_frsh, RCP_45_100_AP_shlf_frsh, RCP_45_150_AP_shlf_frsh, RCP_45_200_AP_shlf_frsh,\
                                   RCP_45_250_AP_shlf_frsh, RCP_45_300_AP_shlf_frsh, RCP_45_350_AP_shlf_frsh, RCP_45_400_AP_shlf_frsh,\
                                   RCP_45_450_AP_shlf_frsh, RCP_45_500_AP_shlf_frsh,\
                                   RCP_85_050_AP_shlf_frsh, RCP_85_100_AP_shlf_frsh, RCP_85_150_AP_shlf_frsh, RCP_85_200_AP_shlf_frsh,\
                                   RCP_85_250_AP_shlf_frsh, RCP_85_300_AP_shlf_frsh, RCP_85_350_AP_shlf_frsh, RCP_85_400_AP_shlf_frsh,\
                                   RCP_85_450_AP_shlf_frsh, RCP_85_500_AP_shlf_frsh])

                    all_sc_lst.append(sc_lst)
                except FileNotFoundError:
                    #print('Model simulation did not terminate')
                    continue

    all_sc_lst = [i[0] for i in all_sc_lst]
    df_sum = pd.DataFrame(all_sc_lst, columns = sum_headers_shlf)
    df_sum.to_csv(tot_summary_shlf_csv)            

    #   select the different DEM and SLR scenarios
    merit_sc = [1, 2, 3, 4, 5, 6, 7, 8]  # gold
    coastal_sc = [9, 10, 11, 12, 13, 14, 15, 16] # greenyellow
    gebco_sc = [17, 18, 19, 20, 21, 22, 23, 24]

    def get_frsh_lst(in_csv, sc_lst, rcp):
        df = pd.read_csv(in_csv, index_col = False, header = 0)    
        out_lst = []
        for sc in sc_lst:
            try:
                lst = df.loc[:, df.columns.str.contains('_BP_')].loc[df['sc_id'] == sc].values.tolist()[0]
                lst_rcp = df.loc[:, df.columns.str.contains(rcp)].loc[df['sc_id'] == sc].values.tolist()[0]
                out_lst.append(lst + lst_rcp)
            except IndexError:
                pass
        return out_lst

    merit_26 = get_frsh_lst(tot_summary_inl_csv, merit_sc, 'RCP_26')
    merit_45 = get_frsh_lst(tot_summary_inl_csv, merit_sc, 'RCP_45')
    merit_85 = get_frsh_lst(tot_summary_inl_csv, merit_sc, 'RCP_85')
    coastal_26 = get_frsh_lst(tot_summary_inl_csv, coastal_sc, 'RCP_26')
    coastal_45 = get_frsh_lst(tot_summary_inl_csv, coastal_sc, 'RCP_45')
    coastal_85 = get_frsh_lst(tot_summary_inl_csv, coastal_sc, 'RCP_85')
    gebco_26 = get_frsh_lst(tot_summary_inl_csv, gebco_sc, 'RCP_26')
    gebco_45 = get_frsh_lst(tot_summary_inl_csv, gebco_sc, 'RCP_45')
    gebco_85 = get_frsh_lst(tot_summary_inl_csv, gebco_sc, 'RCP_85')

    merit_26_shlf = get_frsh_lst(tot_summary_shlf_csv, merit_sc, 'RCP_26')
    merit_45_shlf = get_frsh_lst(tot_summary_shlf_csv, merit_sc, 'RCP_45')
    merit_85_shlf = get_frsh_lst(tot_summary_shlf_csv, merit_sc, 'RCP_85')
    coastal_26_shlf = get_frsh_lst(tot_summary_shlf_csv, coastal_sc, 'RCP_26')
    coastal_45_shlf = get_frsh_lst(tot_summary_shlf_csv, coastal_sc, 'RCP_45')
    coastal_85_shlf = get_frsh_lst(tot_summary_shlf_csv, coastal_sc, 'RCP_85')
    gebco_26_shlf = get_frsh_lst(tot_summary_shlf_csv, gebco_sc, 'RCP_26')
    gebco_45_shlf = get_frsh_lst(tot_summary_shlf_csv, gebco_sc, 'RCP_45')
    gebco_85_shlf = get_frsh_lst(tot_summary_shlf_csv, gebco_sc, 'RCP_85')

    x_lst = [-20000, -19000, -18000, -17000, -16000, -15000, -14000, -13000, -12000, -11000, -10000, -9000, -8000, -7000,\
             -6000, -5000, -4000, -3000, -2000, -1900, -1800, -1700, -1600, -1500, -1400, -1300, -1200, -1100, -1000, -900,\
             -800, -700, -600, -500, -400, -300, -200, -100, 0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

    x_lst = [-38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25,\
             -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9,\
             -8, -7, -6, -5, -4, -3, -2, -1, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]

    #   make a graph with the fresh water concentrations in time
    fig = plt.figure(figsize = (9, 8))
    ax1 = plt.subplot2grid((5, 2), (1, 0), rowspan = 2)
    ax2 = plt.subplot2grid((5, 2), (1, 1))
    ax3 = plt.subplot2grid((5, 2), (2, 1))    
    ax4 = plt.subplot2grid((5, 2), (3, 1))    
    ax5 = plt.subplot2grid((5, 2), (4, 1))   
    
    ax1.set_position([0.025, 0.4, 0.05, 0.3]) # [left, bottom, width, height]
    ax2.set_position([0.09, 0.6, 0.85, 0.375]) 
    ax2.set_facecolor('whitesmoke')
    ax3.set_position([0.09, 0.22, 0.85, 0.375])     
    ax3.set_facecolor('whitesmoke')
    ax4.set_position([0.2, 0.125, 0.6, 0.05]) # [left, bottom, width, height]    
    ax5.set_position([0.09, 0.025, 0.85, 0.075])     
    
    x_axis_label = 'Time (ka BP)'
    y_axis_label = '% fresh groundwater (< 1 g/l TDS)'# compared to grid resolution 100m / 10m'
    ax1.text(0.15, 0.5, y_axis_label, horizontalalignment = 'center', verticalalignment = 'center', transform = ax1.transAxes, fontsize = 9, rotation = 'vertical')
    ax1.axis('off')
    ax4.text(0.5, 0.6, x_axis_label, horizontalalignment = 'center', verticalalignment = 'center', transform = ax4.transAxes, fontsize = 9)
    ax4.axis('off')

    """
    for i in range(len(merit_26)):
        ax2.plot(x_lst, merit_26[i], color = 'orange', linewidth = 0.5, alpha = 0.4)
    for i in range(len(merit_45)):
        ax2.plot(x_lst, merit_45[i], color = 'red', linewidth = 0.5, alpha = 0.4)
    for i in range(len(merit_85)):
        ax2.plot(x_lst, merit_85[i], color = 'darkred', linewidth = 0.5, alpha = 0.4)
    """
    if len(merit_26) > 0:
        ax2.plot(x_lst, np.mean(np.array(merit_26), axis = 0).tolist(), color = 'greenyellow', linewidth = 1.5, label = 'Merit_2.6')
    if len(merit_45) > 0:
        ax2.plot(x_lst, np.mean(np.array(merit_45), axis = 0).tolist(), color = 'gold', linewidth = 1.5, label = 'Merit_4.5')
    if len(merit_85) > 0:
        ax2.plot(x_lst, np.mean(np.array(merit_85), axis = 0).tolist(), color = 'lightskyblue', linewidth = 1.5, label = 'Merit_8.5')    
    if len(coastal_26) > 0:
        ax2.plot(x_lst, np.mean(np.array(coastal_26), axis = 0).tolist(), color = 'limegreen', linewidth = 1.5, label = 'CoastalDEM_2.6')
    if len(coastal_45) > 0:
        ax2.plot(x_lst, np.mean(np.array(coastal_45), axis = 0).tolist(), color = 'darkorange', linewidth = 1.5, label = 'CoastalDEM_4.5')
    if len(coastal_85) > 0:
        ax2.plot(x_lst, np.mean(np.array(coastal_85), axis = 0).tolist(), color = 'dodgerblue', linewidth = 1.5, label = 'CoastalDEM_8.5')    
    if len(gebco_26) > 0:
        ax2.plot(x_lst, np.mean(np.array(gebco_26), axis = 0).tolist(), color = 'forestgreen', linewidth = 1.5, label = 'GEBCO_2.6')
    if len(gebco_45) > 0: 
        ax2.plot(x_lst, np.mean(np.array(gebco_45), axis = 0).tolist(), color = 'red', linewidth = 1.5, label = 'GEBCO_4.5')
    if len(gebco_85) > 0:
        ax2.plot(x_lst, np.mean(np.array(gebco_85), axis = 0).tolist(), color = 'navy', linewidth = 1.5, label = 'GEBCO_8.5')    

    x_major_locator = matplotlib.ticker.FixedLocator([-35, -30, -25, -20, -15, -10, -5, 0, 5], nbins = None)   
    ax2.xaxis.set_major_locator(x_major_locator)
    ax2.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([-38, -37, -36, -34, -33, -32, -31, -29, -28, -27, -26,\
                                                                -24, -23, -22, -21, -19, -18, -17, -16, -14, -13, -12, -11, -9,\
                                                                -8, -7, -6, -4, -3, -2, -1, 1, 2, 3, 4], nbins = None))
    ax2.xaxis.set_ticklabels([])
    ax2.set_xlim([-38, 5])
    ax2.grid(which = 'major', color = 'grey', linestyle = '-', linewidth = 1.25, zorder = 1)    
    ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = 0.5, zorder = 1)    
    ax2.grid(axis='y', color = 'grey', alpha = 0.5, linewidth = 0.5)

    if len(merit_26_shlf) > 0:
        ax3.plot(x_lst, np.mean(np.array(merit_26_shlf), axis = 0).tolist(), color = 'greenyellow', linewidth = 1.5)
    if len(merit_45_shlf) > 0:
        ax3.plot(x_lst, np.mean(np.array(merit_45_shlf), axis = 0).tolist(), color = 'gold', linewidth = 1.5)
    if len(merit_85_shlf) > 0:
        ax3.plot(x_lst, np.mean(np.array(merit_85_shlf), axis = 0).tolist(), color = 'lightskyblue', linewidth = 1.5)  
    if len(coastal_26_shlf) > 0:
        ax3.plot(x_lst, np.mean(np.array(coastal_26_shlf), axis = 0).tolist(), color = 'limegreen', linewidth = 1.5)
    if len(coastal_45_shlf) > 0:
        ax3.plot(x_lst, np.mean(np.array(coastal_45_shlf), axis = 0).tolist(), color = 'darkorange', linewidth = 1.5)
    if len(coastal_85_shlf) > 0:
        ax3.plot(x_lst, np.mean(np.array(coastal_85_shlf), axis = 0).tolist(), color = 'dodgerblue', linewidth = 1.5)    
    if len(gebco_26_shlf) > 0:
        ax3.plot(x_lst, np.mean(np.array(gebco_26_shlf), axis = 0).tolist(), color = 'forestgreen', linewidth = 1.5)
    if len(gebco_45_shlf) > 0:
        ax3.plot(x_lst, np.mean(np.array(gebco_45_shlf), axis = 0).tolist(), color = 'red', linewidth = 1.5)
    if len(gebco_85_shlf) > 0:
        ax3.plot(x_lst, np.mean(np.array(gebco_85_shlf), axis = 0).tolist(), color = 'navy', linewidth = 1.5)    

    ax3.xaxis.set_major_locator(x_major_locator)
    ax3.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([-38, -37, -36, -34, -33, -32, -31, -29, -28, -27, -26,\
                                                                -24, -23, -22, -21, -19, -18, -17, -16, -14, -13, -12, -11, -9,\
                                                                -8, -7, -6, -4, -3, -2, -1, 1, 2, 3, 4], nbins = None))
    ax3.set_xticklabels(['-17', '-12', '-7', '-2', '-1.5', '-1', '-0.5', '0', '0.5'], fontsize = 7)
    #ax2.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([5., 9., 12., 14.], nbins = None))
    ax3.set_xlim([-38, 5])
    ax3.grid(which = 'major', color = 'grey', linestyle = '-', linewidth = 1.25, zorder = 1)    
    ax3.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = 0.5, zorder = 1)    

    h, l = ax2.get_legend_handles_labels() # get labels and handles from ax1
    ax5.legend(h, l, loc = 10, ncol = 3) 
    ax5.axis('off')
    
    plot_name = os.path.join(sum_dir, '_frsh_pct_' + 'COSCAT_' + simu_name + '.png')
    plt.savefig(plot_name, dpi = 600)
    plt.close(fig)    


                
    
"""

sum_dir = summary_dir
simu_name = model_name_basis

"""
            
            
def create_fin_coscat_run_csv_files(coscat_id_str, main_dir, sum_dir, eq_lst, sp_duration, simu_name):

    #   define the summary dir and the model simulation dir
    #sum_dir = os.path.join(main_dir, '_summary_' + coscat_id_str + '_' + cst_type + simu_name)
    #model_simu_dir = os.path.join(main_dir, 'coscat_' + coscat_id_str + '_' + cst_type + '_' + simu_name)

    #   define the output name of the final csv file
    frsh_vol_csv = os.path.join(sum_dir, '_frsh_vol_summary_' + 'COSCAT_' + simu_name + '.csv')
    tot_summary_inl_csv = os.path.join(sum_dir, '_frsh_INL_' + 'COSCAT_' + simu_name + '.csv')
    tot_summary_shlf_csv = os.path.join(sum_dir, '_frsh_SHLF_' + 'COSCAT_' + simu_name + '.csv')
    
    #   find all csv files in the folder 
    eq_0_ts_lst = []
    csv_dir = os.path.join(sum_dir, '_csv_results')
    #print(csv_dir)
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            #print(filename)
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            #simu_nr_int = int(simu_nr_str)        
            #   get all the time steps for eq_0
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)
            col_names = list(df.columns.values)
            eq_0_rows = df.loc[df['EQ'] == 'BP_30000_to_20000']

            #   get the last value of the time_end column, only if the model simulation converged, if not append a -1
            try:
                if len(df.loc[df['EQ'] == 'AP_00400_to_00500'][' time_end'].values) == 2:
                    eq_0_last_ts = eq_0_rows.tail(1).values[0][2] 
                    eq_0_ts_lst.append([eq_0_last_ts, simu_nr_str])
            except IndexError:
                #print('The csv file is empty (most probably)')
                continue

    #   get the index of the tuple with the maximum time step
    max_idx = [i for i, tupl in enumerate(eq_0_ts_lst) if tupl[0] == max([x[0] for x in eq_0_ts_lst])][0]
    max_ts_0 = eq_0_ts_lst[max_idx][0]

    #   column containing the time steps and EQ numbers, taken form the simulation with the longest EQ_0
    ts_csv_dir = os.path.join(csv_dir, 'sc_' + eq_0_ts_lst[max_idx][1] + '_summary_EQs_output.csv')
    ts_csv = pd.read_csv(ts_csv_dir, index_col = False, header = 0)
    ts_csv_eq_col = list(ts_csv['EQ'].values)
    ts_csv_ts_col = list(ts_csv[' time_end'].values)
    
    #   create a list of headers 
    frsh_headers = ['eq']
    for a in range(len(eq_0_ts_lst)):
        frsh_headers.append(eq_0_ts_lst[a][1] + '_frsh_inl')
        frsh_headers.append(eq_0_ts_lst[a][1] + '_frsh_shlf')
    frsh_headers.append('avg_frsh_inl')
    frsh_headers.append('avg_frsh_shlf')
    
    df_frsh_pct = pd.DataFrame(columns = frsh_headers)
    dict_frsh_pct = {}
    dict_frsh_pct['eq'] = ts_csv_eq_col
    frsh_inl_arr = np.zeros((len(eq_0_ts_lst), len(ts_csv_eq_col)))
    frsh_shlf_arr = np.zeros((len(eq_0_ts_lst), len(ts_csv_eq_col)))
    idx_lay = 0
    
    frsh_inl_all_lst, frsh_shlf_all_lst, headers_lst = [], [], []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)
            
            if len(df.loc[df['EQ'] == 'AP_00400_to_00500'][' time_end'].values) == 2:
                if len(df.loc[df['EQ'] == 'BP_30000_to_20000'][' time_end'].values) == 20:
                    inl_lst = df[' frsh_inl_pct'].values.tolist()
                    shlf_lst = df[' frsh_shelf_pct'].values.tolist()
                    
                else:
                    inl_lst_eq_0 = df.loc[df['EQ'] == 'BP_30000_to_20000'][' frsh_inl_pct'].values.tolist()
                    shlf_lst_eq_0 = df.loc[df['EQ'] == 'BP_30000_to_20000'][' frsh_shelf_pct'].values.tolist()
                    for f in range(20 - len(inl_lst_eq_0)):
                        inl_lst_eq_0.append(inl_lst_eq_0[-1])
                        shlf_lst_eq_0.append(shlf_lst_eq_0[-1])
                    inl_lst_rest = df.loc[df['EQ'] != 'BP_30000_to_20000'][' frsh_inl_pct'].values.tolist()
                    shlf_lst_rest = df.loc[df['EQ'] != 'BP_30000_to_20000'][' frsh_shelf_pct'].values.tolist()                    
                    
                    inl_lst = inl_lst_eq_0 + inl_lst_rest
                    shlf_lst = shlf_lst_eq_0 + shlf_lst_rest
                    
                #   set the column headers
                hd_frsh_inl = simu_nr_str + '_frsh_inl'
                hd_frsh_shlf = simu_nr_str + '_frsh_shlf'
                
                inl_lst.insert(0, inl_lst[0])
                shlf_lst.insert(0, shlf_lst[0])
                
                headers_lst.append([hd_frsh_inl, hd_frsh_shlf])                            
                frsh_inl_all_lst.append(inl_lst)
                frsh_shlf_all_lst.append(shlf_lst)
        
    inl_frsh_arr = np.array(frsh_inl_all_lst)
    shlf_frsh_arr = np.array(frsh_shlf_all_lst)        
        
    avg_inl_frsh_lst = list(np.mean(inl_frsh_arr, axis = 0))
    avg_shlf_frsh_lst = list(np.mean(shlf_frsh_arr, axis = 0))    
    
    avg_inl_frsh_lst = [round(i, 2) for i in avg_inl_frsh_lst]
    avg_shlf_frsh_lst = [round(i, 2) for i in avg_shlf_frsh_lst]

    eq_lst_to_dict, ts_lst_to_dict = [], []     
    for g in range(int(max_ts_0 / 500) + 1):
        eq_lst_to_dict.append('BP_30000_to_20000')
        ts_lst_to_dict.append(500. * g)
    time = ts_lst_to_dict[-1] + 500
    for h in range(1, len(eq_lst)):
        eq_name = eq_lst[h]
        ts_dur = sp_duration[h]    
        ts_step = ts_dur / 2
        for i in range(1, int(ts_dur / (ts_dur / 2) + 1)):
            eq_lst_to_dict.append(eq_name)
            ts_lst_to_dict.append(time)        
            time += ts_step
        if eq_name == 'BP_03000_to_02000':
            time -= 450

    dict_frsh_pct['eq'] = eq_lst_to_dict
    dict_frsh_pct['time_end'] = ts_lst_to_dict                
            
    #   fill the dictionary
    for f in range(len(headers_lst)):
        dict_frsh_pct[headers_lst[f][0]] = frsh_inl_all_lst[f]
        dict_frsh_pct[headers_lst[f][1]] = frsh_shlf_all_lst[f]
    
    dict_frsh_pct['avg_frsh_inl'] = avg_inl_frsh_lst
    dict_frsh_pct['avg_frsh_shlf'] = avg_shlf_frsh_lst    
        
    df_out = pd.DataFrame.from_dict(dict_frsh_pct, orient = 'index').T   
    df_out.to_csv(frsh_vol_csv)

    """         create the summary csv files, with % fresh at certain times and the convergence times         """
    #geo_csv_dir = os.path.join(sum_dir, '_geo_summary_' + coscat_id_str + '.csv')
    #geo_df = pd.read_csv(geo_csv_dir, index_col = False, header = 0)    
    #geo_headers = geo_df.columns.values
    
    """
    sum_headers = ['coscat_id', 'simu_id', 'eq_0_t_end', 'eq_0_t_end_inl_frsh', 'eq_0_t_end_shlf_frsh',\
                   'eq_21_t_end_inl_frsh', 'eq_21_t_end_shlf_frsh', 'eq_31_t_end_inl_frsh', 'eq_31_t_end_shlf_frsh',\
                   'eq_21_CT_t_end', 'eq_21_CT_t_end_inl_frsh', 'eq_21_CT_t_end_shlf_frsh', 'eq_31_CT_t_end', 'eq_31_CT_t_end_inl_frsh', 'eq_31_CT_t_end_shlf_frsh']
    """
    
    sum_headers_inl = ['coscat_id', 'srm_id', 'sc_id', 'eq_0_t_end', 'eq_20000_BP_inl_frsh', 'eq_19000_BP_inl_frsh',\
                       'eq_18000_BP_inl_frsh', 'eq_17000_BP_inl_frsh','eq_16000_BP_inl_frsh', 'eq_15000_BP_inl_frsh',\
                       'eq_14000_BP_inl_frsh', 'eq_13000_BP_inl_frsh','eq_12000_BP_inl_frsh', 'eq_11000_BP_inl_frsh',\
                       'eq_10000_BP_inl_frsh', 'eq_09000_BP_inl_frsh','eq_08000_BP_inl_frsh', 'eq_07000_BP_inl_frsh',\
                       'eq_06000_BP_inl_frsh', 'eq_05000_BP_inl_frsh','eq_04000_BP_inl_frsh', 'eq_03000_BP_inl_frsh',\
                       'eq_02000_BP_inl_frsh', 'eq_01900_BP_inl_frsh','eq_01800_BP_inl_frsh', 'eq_01700_BP_inl_frsh',\
                       'eq_01600_BP_inl_frsh', 'eq_01500_BP_inl_frsh','eq_01400_BP_inl_frsh', 'eq_01300_BP_inl_frsh',\
                       'eq_01200_BP_inl_frsh', 'eq_01100_BP_inl_frsh','eq_01000_BP_inl_frsh', 'eq_00900_BP_inl_frsh',\
                       'eq_00800_BP_inl_frsh', 'eq_00700_BP_inl_frsh','eq_00600_BP_inl_frsh', 'eq_00500_BP_inl_frsh',\
                       'eq_00400_BP_inl_frsh', 'eq_00300_BP_inl_frsh','eq_00200_BP_inl_frsh', 'eq_00100_BP_inl_frsh',\
                       'eq_00000_BP_inl_frsh', 'eq_00100_AP_inl_frsh','eq_00200_AP_inl_frsh', 'eq_00300_AP_inl_frsh',\
                       'eq_00400_BP_inl_frsh', 'eq_00500_AP_inl_frsh']

    all_sc_lst = []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            #print(simu_nr_str)
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)    
            if len(df.loc[df['EQ'] == 'AP_00400_to_00500'][' time_end'].values) == 2:
                try:                      
                    sc_lst = []
                    #   read the specific rows
                    eq_0_t_end = df[(df['EQ'] == eq_lst[0])][' time_end'].values[-1]
                    eq_20000_BP_inl_frsh = df[(df['EQ'] == eq_lst[0])][' frsh_inl_pct'].values[-1]
                    eq_19000_BP_inl_frsh = df[(df['EQ'] == eq_lst[1])][' frsh_inl_pct'].values[-1]                
                    eq_18000_BP_inl_frsh = df[(df['EQ'] == eq_lst[2])][' frsh_inl_pct'].values[-1]
                    eq_17000_BP_inl_frsh = df[(df['EQ'] == eq_lst[3])][' frsh_inl_pct'].values[-1]     
                    eq_16000_BP_inl_frsh = df[(df['EQ'] == eq_lst[4])][' frsh_inl_pct'].values[-1]     
                    eq_15000_BP_inl_frsh = df[(df['EQ'] == eq_lst[5])][' frsh_inl_pct'].values[-1]                 
                    eq_14000_BP_inl_frsh = df[(df['EQ'] == eq_lst[6])][' frsh_inl_pct'].values[-1]
                    eq_13000_BP_inl_frsh = df[(df['EQ'] == eq_lst[7])][' frsh_inl_pct'].values[-1]     
                    eq_12000_BP_inl_frsh = df[(df['EQ'] == eq_lst[8])][' frsh_inl_pct'].values[-1]     
                    eq_11000_BP_inl_frsh = df[(df['EQ'] == eq_lst[9])][' frsh_inl_pct'].values[-1]                    
                    eq_10000_BP_inl_frsh = df[(df['EQ'] == eq_lst[10])][' frsh_inl_pct'].values[-1]
                    eq_09000_BP_inl_frsh = df[(df['EQ'] == eq_lst[11])][' frsh_inl_pct'].values[-1]     
                    eq_08000_BP_inl_frsh = df[(df['EQ'] == eq_lst[12])][' frsh_inl_pct'].values[-1]     
                    eq_07000_BP_inl_frsh = df[(df['EQ'] == eq_lst[13])][' frsh_inl_pct'].values[-1]                            
                    eq_06000_BP_inl_frsh = df[(df['EQ'] == eq_lst[14])][' frsh_inl_pct'].values[-1]
                    eq_05000_BP_inl_frsh = df[(df['EQ'] == eq_lst[15])][' frsh_inl_pct'].values[-1]     
                    eq_04000_BP_inl_frsh = df[(df['EQ'] == eq_lst[16])][' frsh_inl_pct'].values[-1]     
                    eq_03000_BP_inl_frsh = df[(df['EQ'] == eq_lst[17])][' frsh_inl_pct'].values[-1]                
                    eq_02000_BP_inl_frsh = df[(df['EQ'] == eq_lst[18])][' frsh_inl_pct'].values[-1]   
                    eq_01900_BP_inl_frsh = df[(df['EQ'] == eq_lst[19])][' frsh_inl_pct'].values[-1]     
                    eq_01800_BP_inl_frsh = df[(df['EQ'] == eq_lst[20])][' frsh_inl_pct'].values[-1]     
                    eq_01700_BP_inl_frsh = df[(df['EQ'] == eq_lst[21])][' frsh_inl_pct'].values[-1]                            
                    eq_01600_BP_inl_frsh = df[(df['EQ'] == eq_lst[22])][' frsh_inl_pct'].values[-1]
                    eq_01500_BP_inl_frsh = df[(df['EQ'] == eq_lst[23])][' frsh_inl_pct'].values[-1]     
                    eq_01400_BP_inl_frsh = df[(df['EQ'] == eq_lst[24])][' frsh_inl_pct'].values[-1]     
                    eq_01300_BP_inl_frsh = df[(df['EQ'] == eq_lst[25])][' frsh_inl_pct'].values[-1]                
                    eq_01200_BP_inl_frsh = df[(df['EQ'] == eq_lst[26])][' frsh_inl_pct'].values[-1]   
                    eq_01100_BP_inl_frsh = df[(df['EQ'] == eq_lst[27])][' frsh_inl_pct'].values[-1]                
                    eq_01000_BP_inl_frsh = df[(df['EQ'] == eq_lst[28])][' frsh_inl_pct'].values[-1]   
                    eq_00900_BP_inl_frsh = df[(df['EQ'] == eq_lst[29])][' frsh_inl_pct'].values[-1]     
                    eq_00800_BP_inl_frsh = df[(df['EQ'] == eq_lst[30])][' frsh_inl_pct'].values[-1]     
                    eq_00700_BP_inl_frsh = df[(df['EQ'] == eq_lst[31])][' frsh_inl_pct'].values[-1]                            
                    eq_00600_BP_inl_frsh = df[(df['EQ'] == eq_lst[32])][' frsh_inl_pct'].values[-1]
                    eq_00500_BP_inl_frsh = df[(df['EQ'] == eq_lst[33])][' frsh_inl_pct'].values[-1]     
                    eq_00400_BP_inl_frsh = df[(df['EQ'] == eq_lst[34])][' frsh_inl_pct'].values[-1]     
                    eq_00300_BP_inl_frsh = df[(df['EQ'] == eq_lst[35])][' frsh_inl_pct'].values[-1]                
                    eq_00200_BP_inl_frsh = df[(df['EQ'] == eq_lst[36])][' frsh_inl_pct'].values[-1]   
                    eq_00100_BP_inl_frsh = df[(df['EQ'] == eq_lst[37])][' frsh_inl_pct'].values[-1]                
                    eq_00000_BP_inl_frsh = df[(df['EQ'] == eq_lst[38])][' frsh_inl_pct'].values[-1]   
                    eq_00100_AP_inl_frsh = df[(df['EQ'] == eq_lst[39])][' frsh_inl_pct'].values[-1]     
                    eq_00200_AP_inl_frsh = df[(df['EQ'] == eq_lst[40])][' frsh_inl_pct'].values[-1]     
                    eq_00300_AP_inl_frsh = df[(df['EQ'] == eq_lst[41])][' frsh_inl_pct'].values[-1]                
                    eq_00400_AP_inl_frsh = df[(df['EQ'] == eq_lst[42])][' frsh_inl_pct'].values[-1]   
                    eq_00500_AP_inl_frsh = df[(df['EQ'] == eq_lst[43])][' frsh_inl_pct'].values[-1]                
    
                    sc_lst.append([coscat_id_str, simu_name.split('SRM_')[-1], simu_nr_str, eq_0_t_end, eq_20000_BP_inl_frsh,\
                                   eq_19000_BP_inl_frsh, eq_18000_BP_inl_frsh, eq_17000_BP_inl_frsh, eq_16000_BP_inl_frsh,\
                                   eq_15000_BP_inl_frsh, eq_14000_BP_inl_frsh, eq_13000_BP_inl_frsh, eq_12000_BP_inl_frsh,\
                                   eq_11000_BP_inl_frsh, eq_10000_BP_inl_frsh, eq_09000_BP_inl_frsh, eq_08000_BP_inl_frsh,\
                                   eq_07000_BP_inl_frsh, eq_06000_BP_inl_frsh, eq_05000_BP_inl_frsh, eq_04000_BP_inl_frsh,\
                                   eq_03000_BP_inl_frsh, eq_02000_BP_inl_frsh, eq_01900_BP_inl_frsh, eq_01800_BP_inl_frsh,\
                                   eq_01700_BP_inl_frsh, eq_01600_BP_inl_frsh, eq_01500_BP_inl_frsh, eq_01400_BP_inl_frsh,\
                                   eq_01300_BP_inl_frsh, eq_01200_BP_inl_frsh, eq_01100_BP_inl_frsh, eq_01000_BP_inl_frsh,\
                                   eq_00900_BP_inl_frsh, eq_00800_BP_inl_frsh, eq_00700_BP_inl_frsh, eq_00600_BP_inl_frsh,\
                                   eq_00500_BP_inl_frsh, eq_00400_BP_inl_frsh, eq_00300_BP_inl_frsh, eq_00200_BP_inl_frsh,\
                                   eq_00100_BP_inl_frsh, eq_00000_BP_inl_frsh, eq_00100_AP_inl_frsh, eq_00200_AP_inl_frsh,\
                                   eq_00300_AP_inl_frsh, eq_00400_AP_inl_frsh, eq_00500_AP_inl_frsh])
            
                    all_sc_lst.append(sc_lst)
                except FileNotFoundError:
                    #print('Model simulation did not terminate')
                    continue

    all_sc_lst = [i[0] for i in all_sc_lst]
    df_sum = pd.DataFrame(all_sc_lst, columns = sum_headers_inl)
    df_sum.to_csv(tot_summary_inl_csv)            

 
    sum_headers_shlf = ['coscat_id', 'srm_id', 'sc_id', 'eq_0_t_end', 'eq_20000_BP_shlf_frsh', 'eq_19000_BP_shlf_frsh',\
                       'eq_18000_BP_shlf_frsh', 'eq_17000_BP_shlf_frsh','eq_16000_BP_shlf_frsh', 'eq_15000_BP_shlf_frsh',\
                       'eq_14000_BP_shlf_frsh', 'eq_13000_BP_shlf_frsh','eq_12000_BP_shlf_frsh', 'eq_11000_BP_shlf_frsh',\
                       'eq_10000_BP_shlf_frsh', 'eq_09000_BP_shlf_frsh','eq_08000_BP_shlf_frsh', 'eq_07000_BP_shlf_frsh',\
                       'eq_06000_BP_shlf_frsh', 'eq_05000_BP_shlf_frsh','eq_04000_BP_shlf_frsh', 'eq_03000_BP_shlf_frsh',\
                       'eq_02000_BP_shlf_frsh', 'eq_01900_BP_shlf_frsh','eq_01800_BP_shlf_frsh', 'eq_01700_BP_shlf_frsh',\
                       'eq_01600_BP_shlf_frsh', 'eq_01500_BP_shlf_frsh','eq_01400_BP_shlf_frsh', 'eq_01300_BP_shlf_frsh',\
                       'eq_01200_BP_shlf_frsh', 'eq_01100_BP_shlf_frsh','eq_01000_BP_shlf_frsh', 'eq_00900_BP_shlf_frsh',\
                       'eq_00800_BP_shlf_frsh', 'eq_00700_BP_shlf_frsh','eq_00600_BP_shlf_frsh', 'eq_00500_BP_shlf_frsh',\
                       'eq_00400_BP_shlf_frsh', 'eq_00300_BP_shlf_frsh','eq_00200_BP_shlf_frsh', 'eq_00100_BP_shlf_frsh',\
                       'eq_00000_BP_shlf_frsh', 'eq_00100_AP_shlf_frsh','eq_00200_AP_shlf_frsh', 'eq_00300_AP_shlf_frsh',\
                       'eq_00400_BP_shlf_frsh', 'eq_00500_AP_shlf_frsh']

    all_sc_lst = []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            #print(simu_nr_str)
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)    
            if len(df.loc[df['EQ'] == 'AP_00400_to_00500'][' time_end'].values) == 2:
                try:                      
                    sc_lst = []
                    #   read the specific rows
                    eq_0_t_end = df[(df['EQ'] == eq_lst[0])][' time_end'].values[-1]
                    eq_20000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[0])][' frsh_shelf_pct'].values[-1]
                    eq_19000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[1])][' frsh_shelf_pct'].values[-1]                
                    eq_18000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[2])][' frsh_shelf_pct'].values[-1]
                    eq_17000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[3])][' frsh_shelf_pct'].values[-1]     
                    eq_16000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[4])][' frsh_shelf_pct'].values[-1]     
                    eq_15000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[5])][' frsh_shelf_pct'].values[-1]                 
                    eq_14000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[6])][' frsh_shelf_pct'].values[-1]
                    eq_13000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[7])][' frsh_shelf_pct'].values[-1]     
                    eq_12000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[8])][' frsh_shelf_pct'].values[-1]     
                    eq_11000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[9])][' frsh_shelf_pct'].values[-1]                    
                    eq_10000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[10])][' frsh_shelf_pct'].values[-1]
                    eq_09000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[11])][' frsh_shelf_pct'].values[-1]     
                    eq_08000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[12])][' frsh_shelf_pct'].values[-1]     
                    eq_07000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[13])][' frsh_shelf_pct'].values[-1]                            
                    eq_06000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[14])][' frsh_shelf_pct'].values[-1]
                    eq_05000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[15])][' frsh_shelf_pct'].values[-1]     
                    eq_04000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[16])][' frsh_shelf_pct'].values[-1]     
                    eq_03000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[17])][' frsh_shelf_pct'].values[-1]                
                    eq_02000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[18])][' frsh_shelf_pct'].values[-1]   
                    eq_01900_BP_shlf_frsh = df[(df['EQ'] == eq_lst[19])][' frsh_shelf_pct'].values[-1]     
                    eq_01800_BP_shlf_frsh = df[(df['EQ'] == eq_lst[20])][' frsh_shelf_pct'].values[-1]     
                    eq_01700_BP_shlf_frsh = df[(df['EQ'] == eq_lst[21])][' frsh_shelf_pct'].values[-1]                            
                    eq_01600_BP_shlf_frsh = df[(df['EQ'] == eq_lst[22])][' frsh_shelf_pct'].values[-1]
                    eq_01500_BP_shlf_frsh = df[(df['EQ'] == eq_lst[23])][' frsh_shelf_pct'].values[-1]     
                    eq_01400_BP_shlf_frsh = df[(df['EQ'] == eq_lst[24])][' frsh_shelf_pct'].values[-1]     
                    eq_01300_BP_shlf_frsh = df[(df['EQ'] == eq_lst[25])][' frsh_shelf_pct'].values[-1]                
                    eq_01200_BP_shlf_frsh = df[(df['EQ'] == eq_lst[26])][' frsh_shelf_pct'].values[-1]   
                    eq_01100_BP_shlf_frsh = df[(df['EQ'] == eq_lst[27])][' frsh_shelf_pct'].values[-1]                
                    eq_01000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[28])][' frsh_shelf_pct'].values[-1]   
                    eq_00900_BP_shlf_frsh = df[(df['EQ'] == eq_lst[29])][' frsh_shelf_pct'].values[-1]     
                    eq_00800_BP_shlf_frsh = df[(df['EQ'] == eq_lst[30])][' frsh_shelf_pct'].values[-1]     
                    eq_00700_BP_shlf_frsh = df[(df['EQ'] == eq_lst[31])][' frsh_shelf_pct'].values[-1]                            
                    eq_00600_BP_shlf_frsh = df[(df['EQ'] == eq_lst[32])][' frsh_shelf_pct'].values[-1]
                    eq_00500_BP_shlf_frsh = df[(df['EQ'] == eq_lst[33])][' frsh_shelf_pct'].values[-1]     
                    eq_00400_BP_shlf_frsh = df[(df['EQ'] == eq_lst[34])][' frsh_shelf_pct'].values[-1]     
                    eq_00300_BP_shlf_frsh = df[(df['EQ'] == eq_lst[35])][' frsh_shelf_pct'].values[-1]                
                    eq_00200_BP_shlf_frsh = df[(df['EQ'] == eq_lst[36])][' frsh_shelf_pct'].values[-1]   
                    eq_00100_BP_shlf_frsh = df[(df['EQ'] == eq_lst[37])][' frsh_shelf_pct'].values[-1]                
                    eq_00000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[38])][' frsh_shelf_pct'].values[-1]   
                    eq_00100_AP_shlf_frsh = df[(df['EQ'] == eq_lst[39])][' frsh_shelf_pct'].values[-1]     
                    eq_00200_AP_shlf_frsh = df[(df['EQ'] == eq_lst[40])][' frsh_shelf_pct'].values[-1]     
                    eq_00300_AP_shlf_frsh = df[(df['EQ'] == eq_lst[41])][' frsh_shelf_pct'].values[-1]                
                    eq_00400_AP_shlf_frsh = df[(df['EQ'] == eq_lst[42])][' frsh_shelf_pct'].values[-1]   
                    eq_00500_AP_shlf_frsh = df[(df['EQ'] == eq_lst[43])][' frsh_shelf_pct'].values[-1]                
    
                    sc_lst.append([coscat_id_str, simu_name.split('SRM_')[-1], simu_nr_str, eq_0_t_end, eq_20000_BP_shlf_frsh,\
                                   eq_19000_BP_shlf_frsh, eq_18000_BP_shlf_frsh, eq_17000_BP_shlf_frsh, eq_16000_BP_shlf_frsh,\
                                   eq_15000_BP_shlf_frsh, eq_14000_BP_shlf_frsh, eq_13000_BP_shlf_frsh, eq_12000_BP_shlf_frsh,\
                                   eq_11000_BP_shlf_frsh, eq_10000_BP_shlf_frsh, eq_09000_BP_shlf_frsh, eq_08000_BP_shlf_frsh,\
                                   eq_07000_BP_shlf_frsh, eq_06000_BP_shlf_frsh, eq_05000_BP_shlf_frsh, eq_04000_BP_shlf_frsh,\
                                   eq_03000_BP_shlf_frsh, eq_02000_BP_shlf_frsh, eq_01900_BP_shlf_frsh, eq_01800_BP_shlf_frsh,\
                                   eq_01700_BP_shlf_frsh, eq_01600_BP_shlf_frsh, eq_01500_BP_shlf_frsh, eq_01400_BP_shlf_frsh,\
                                   eq_01300_BP_shlf_frsh, eq_01200_BP_shlf_frsh, eq_01100_BP_shlf_frsh, eq_01000_BP_shlf_frsh,\
                                   eq_00900_BP_shlf_frsh, eq_00800_BP_shlf_frsh, eq_00700_BP_shlf_frsh, eq_00600_BP_shlf_frsh,\
                                   eq_00500_BP_shlf_frsh, eq_00400_BP_shlf_frsh, eq_00300_BP_shlf_frsh, eq_00200_BP_shlf_frsh,\
                                   eq_00100_BP_shlf_frsh, eq_00000_BP_shlf_frsh, eq_00100_AP_shlf_frsh, eq_00200_AP_shlf_frsh,\
                                   eq_00300_AP_shlf_frsh, eq_00400_AP_shlf_frsh, eq_00500_AP_shlf_frsh])
            
                    all_sc_lst.append(sc_lst)
                except FileNotFoundError:
                    #print('Model simulation did not terminate')
                    continue

    all_sc_lst = [i[0] for i in all_sc_lst]
    df_sum = pd.DataFrame(all_sc_lst, columns = sum_headers_shlf)
    df_sum.to_csv(tot_summary_shlf_csv)            

    #   select the different DEM and SLR scenarios
    merit_26_sc = [1, 4, 7, 10, 13, 16, 19, 22]  # gold
    merit_45_sc = [2, 5, 8, 11, 14, 17, 20, 23]  # darkorange
    merit_85_sc = [3, 6, 9, 12, 15, 18, 21, 24]  # red
    coastal_26_sc = [25, 28, 31, 34, 37, 40, 43, 46] # greenyellow
    coastal_45_sc = [26, 29, 32, 35, 38, 41, 44, 47] # limegreen
    coastal_85_sc = [27, 30, 33, 36, 39, 42, 45, 48] # forestgreen
    gebco_26_sc = [49, 52, 55, 58, 61, 64, 67, 70]
    gebco_45_sc = [50, 53, 56, 59, 62, 65, 68, 71]
    gebco_85_sc = [51, 54, 57, 60, 63, 66, 69, 72]

    def get_frsh_lst(in_csv, sc_lst):
        df = pd.read_csv(in_csv, index_col = False, header = 0)    
        out_lst = []
        for sc in sc_lst:
            try:
                lst = df.loc[df['sc_id'] == sc].values.tolist()[0]
                out_lst.append(lst[len(lst) - 44 :])
            except IndexError:
                pass
        return out_lst

    merit_26 = get_frsh_lst(tot_summary_inl_csv, merit_26_sc)
    merit_45 = get_frsh_lst(tot_summary_inl_csv, merit_45_sc)
    merit_85 = get_frsh_lst(tot_summary_inl_csv, merit_85_sc)
    coastal_26 = get_frsh_lst(tot_summary_inl_csv, coastal_26_sc)
    coastal_45 = get_frsh_lst(tot_summary_inl_csv, coastal_45_sc)
    coastal_85 = get_frsh_lst(tot_summary_inl_csv, coastal_85_sc)
    gebco_26 = get_frsh_lst(tot_summary_inl_csv, gebco_26_sc)
    gebco_45 = get_frsh_lst(tot_summary_inl_csv, gebco_45_sc)
    gebco_85 = get_frsh_lst(tot_summary_inl_csv, gebco_85_sc)

    merit_26_shlf = get_frsh_lst(tot_summary_shlf_csv, merit_26_sc)
    merit_45_shlf = get_frsh_lst(tot_summary_shlf_csv, merit_45_sc)
    merit_85_shlf = get_frsh_lst(tot_summary_shlf_csv, merit_85_sc)
    coastal_26_shlf = get_frsh_lst(tot_summary_shlf_csv, coastal_26_sc)
    coastal_45_shlf = get_frsh_lst(tot_summary_shlf_csv, coastal_45_sc)
    coastal_85_shlf = get_frsh_lst(tot_summary_shlf_csv, coastal_85_sc)
    gebco_26_shlf = get_frsh_lst(tot_summary_shlf_csv, gebco_26_sc)
    gebco_45_shlf = get_frsh_lst(tot_summary_shlf_csv, gebco_45_sc)
    gebco_85_shlf = get_frsh_lst(tot_summary_shlf_csv, gebco_85_sc)


    x_lst = [-20000, -19000, -18000, -17000, -16000, -15000, -14000, -13000, -12000, -11000, -10000, -9000, -8000, -7000,\
             -6000, -5000, -4000, -3000, -2000, -1900, -1800, -1700, -1600, -1500, -1400, -1300, -1200, -1100, -1000, -900,\
             -800, -700, -600, -500, -400, -300, -200, -100, 0, 100, 200, 300, 400, 500]

    x_lst = [-38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25,\
             -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9,\
             -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]

    #   make a graph with the fresh water concentrations in time
    fig = plt.figure(figsize = (9, 8))
    ax1 = plt.subplot2grid((5, 2), (1, 0), rowspan = 2)
    ax2 = plt.subplot2grid((5, 2), (1, 1))
    ax3 = plt.subplot2grid((5, 2), (2, 1))    
    ax4 = plt.subplot2grid((5, 2), (3, 1))    
    ax5 = plt.subplot2grid((5, 2), (4, 1))   
    
    ax1.set_position([0.025, 0.4, 0.05, 0.3]) # [left, bottom, width, height]
    ax2.set_position([0.09, 0.6, 0.85, 0.375]) 
    ax2.set_facecolor('whitesmoke')
    ax3.set_position([0.09, 0.22, 0.85, 0.375])     
    ax3.set_facecolor('whitesmoke')
    ax4.set_position([0.2, 0.125, 0.6, 0.05]) # [left, bottom, width, height]    
    ax5.set_position([0.09, 0.025, 0.85, 0.075])     
    
    x_axis_label = 'Time (ka BP)'
    y_axis_label = '% fresh groundwater (< 1 g/l TDS)'# compared to grid resolution 100m / 10m'
    ax1.text(0.15, 0.5, y_axis_label, horizontalalignment = 'center', verticalalignment = 'center', transform = ax1.transAxes, fontsize = 9, rotation = 'vertical')
    ax1.axis('off')
    ax4.text(0.5, 0.6, x_axis_label, horizontalalignment = 'center', verticalalignment = 'center', transform = ax4.transAxes, fontsize = 9)
    ax4.axis('off')

    """
    for i in range(len(merit_26)):
        ax2.plot(x_lst, merit_26[i], color = 'orange', linewidth = 0.5, alpha = 0.4)
    for i in range(len(merit_45)):
        ax2.plot(x_lst, merit_45[i], color = 'red', linewidth = 0.5, alpha = 0.4)
    for i in range(len(merit_85)):
        ax2.plot(x_lst, merit_85[i], color = 'darkred', linewidth = 0.5, alpha = 0.4)
    """
    ax2.plot(x_lst, np.mean(np.array(merit_26), axis = 0).tolist(), color = 'greenyellow', linewidth = 1.5, label = 'Merit_2.6')
    ax2.plot(x_lst, np.mean(np.array(merit_45), axis = 0).tolist(), color = 'gold', linewidth = 1.5, label = 'Merit_4.5')
    ax2.plot(x_lst, np.mean(np.array(merit_85), axis = 0).tolist(), color = 'lightskyblue', linewidth = 1.5, label = 'Merit_8.5')    
    ax2.plot(x_lst, np.mean(np.array(coastal_26), axis = 0).tolist(), color = 'limegreen', linewidth = 1.5, label = 'CoastalDEM_2.6')
    ax2.plot(x_lst, np.mean(np.array(coastal_45), axis = 0).tolist(), color = 'darkorange', linewidth = 1.5, label = 'CoastalDEM_4.5')
    ax2.plot(x_lst, np.mean(np.array(coastal_85), axis = 0).tolist(), color = 'dodgerblue', linewidth = 1.5, label = 'CoastalDEM_8.5')    
    ax2.plot(x_lst, np.mean(np.array(gebco_26), axis = 0).tolist(), color = 'forestgreen', linewidth = 1.5, label = 'GEBCO_2.6')
    ax2.plot(x_lst, np.mean(np.array(gebco_45), axis = 0).tolist(), color = 'red', linewidth = 1.5, label = 'GEBCO_4.5')
    ax2.plot(x_lst, np.mean(np.array(gebco_85), axis = 0).tolist(), color = 'navy', linewidth = 1.5, label = 'GEBCO_8.5')    

    x_major_locator = matplotlib.ticker.FixedLocator([-35, -30, -25, -20, -15, -10, -5, 0, 5], nbins = None)   
    ax2.xaxis.set_major_locator(x_major_locator)
    ax2.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([-38, -37, -36, -34, -33, -32, -31, -29, -28, -27, -26,\
                                                                -24, -23, -22, -21, -19, -18, -17, -16, -14, -13, -12, -11, -9,\
                                                                -8, -7, -6, -4, -3, -2, -1, 1, 2, 3, 4], nbins = None))
    ax2.xaxis.set_ticklabels([])
    ax2.set_xlim([-38, 5])
    ax2.grid(which = 'major', color = 'grey', linestyle = '-', linewidth = 1.25, zorder = 1)    
    ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = 0.5, zorder = 1)    
    ax2.grid(axis='y', color = 'grey', alpha = 0.5, linewidth = 0.5)

    ax3.plot(x_lst, np.mean(np.array(merit_26_shlf), axis = 0).tolist(), color = 'greenyellow', linewidth = 1.5)
    ax3.plot(x_lst, np.mean(np.array(merit_45_shlf), axis = 0).tolist(), color = 'gold', linewidth = 1.5)
    ax3.plot(x_lst, np.mean(np.array(merit_85_shlf), axis = 0).tolist(), color = 'lightskyblue', linewidth = 1.5)    
    ax3.plot(x_lst, np.mean(np.array(coastal_26_shlf), axis = 0).tolist(), color = 'limegreen', linewidth = 1.5)
    ax3.plot(x_lst, np.mean(np.array(coastal_45_shlf), axis = 0).tolist(), color = 'darkorange', linewidth = 1.5)
    ax3.plot(x_lst, np.mean(np.array(coastal_85_shlf), axis = 0).tolist(), color = 'dodgerblue', linewidth = 1.5)    
    ax3.plot(x_lst, np.mean(np.array(gebco_26_shlf), axis = 0).tolist(), color = 'forestgreen', linewidth = 1.5)
    ax3.plot(x_lst, np.mean(np.array(gebco_45_shlf), axis = 0).tolist(), color = 'red', linewidth = 1.5)
    ax3.plot(x_lst, np.mean(np.array(gebco_85_shlf), axis = 0).tolist(), color = 'navy', linewidth = 1.5)    

    ax3.xaxis.set_major_locator(x_major_locator)
    ax3.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([-38, -37, -36, -34, -33, -32, -31, -29, -28, -27, -26,\
                                                                -24, -23, -22, -21, -19, -18, -17, -16, -14, -13, -12, -11, -9,\
                                                                -8, -7, -6, -4, -3, -2, -1, 1, 2, 3, 4], nbins = None))
    ax3.set_xticklabels(['-17', '-12', '-7', '-2', '-1.5', '-1', '-0.5', '0', '0.5'], fontsize = 7)
    #ax2.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([5., 9., 12., 14.], nbins = None))
    ax3.set_xlim([-38, 5])
    ax3.grid(which = 'major', color = 'grey', linestyle = '-', linewidth = 1.25, zorder = 1)    
    ax3.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = 0.5, zorder = 1)    

    h, l = ax2.get_legend_handles_labels() # get labels and handles from ax1
    ax5.legend(h, l, loc = 10, ncol = 3) 
    ax5.axis('off')
    
    plot_name = os.path.join(sum_dir, '_frsh_pct_' + 'COSCAT_' + simu_name + '.png')
    plt.savefig(plot_name, dpi = 600)
    plt.close(fig)    
    
"""
Same function as above, but fitted for the case with freshwater starting condition

t_cur = time_present
sum_dir = summary_dir
simu_name = model_name
"""       
def create_fin_coscat_run_csv_files_FRESH(t_cur, coscat_id, cst_type, main_dir, sum_dir, eq_lst, sp_duration, simu_name):

    #   adjust the coscat_id if necessary by adding 0s in front of the number
    if coscat_id < 10:
        coscat_id_str = '000' + str(coscat_id)
    elif coscat_id >= 10 and coscat_id < 100:
        coscat_id_str = '00' + str(coscat_id)
    elif coscat_id >= 100 and coscat_id < 1000:
        coscat_id_str = '0' + str(coscat_id)
    else:
        coscat_id_str = str(coscat_id)

    #   define the summary dir and the model simulation dir
    #sum_dir = os.path.join(main_dir, '_summary_' + coscat_id_str + '_' + cst_type + simu_name)
    #model_simu_dir = os.path.join(main_dir, 'coscat_' + coscat_id_str + '_' + cst_type + '_' + simu_name)

    #   define the output name of the final csv file
    frsh_vol_csv = os.path.join(sum_dir, '_frsh_vol_summary_' + coscat_id_str + '_' + cst_type + simu_name + '.csv')
    tot_summary_csv = os.path.join(sum_dir, '_total_summary_' + coscat_id_str + '_' + cst_type + simu_name + '.csv')
    
    #   find all csv files in the folder 
    eq_0_ts_lst = []
    csv_dir = os.path.join(sum_dir, '_csv_results')
    print(csv_dir)
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            print(filename)
            simu_nr_str = filename.split('scenario_')[-1].split('_')[0]
            #simu_nr_int = int(simu_nr_str)        
            #   get all the time steps for eq_0
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)
            col_names = list(df.columns.values)
            eq_0_rows = df.loc[df['EQ'] == 0]

            #   get the last value of the time_end column
            try:
                eq_0_last_ts = eq_0_rows.tail(1).values[0][2] 
                eq_0_ts_lst.append([eq_0_last_ts, simu_nr_str])
            except IndexError:
                print('The csv file is empty (most probably)')
                continue

    #   get the index of the tuple with the maximum time step
    print(eq_0_ts_lst)
    max_idx = [i for i, tupl in enumerate(eq_0_ts_lst) if tupl[0] == max([x[0] for x in eq_0_ts_lst])][0]
    max_ts_0 = eq_0_ts_lst[max_idx][0]

    #   column containing the time steps and EQ numbers, taken form the simulation with the longest EQ_0
    ts_csv_dir = os.path.join(csv_dir, 'scenario_' + eq_0_ts_lst[max_idx][1] + '_summary_EQs_output.csv')
    ts_csv = pd.read_csv(ts_csv_dir, index_col = False, header = 0)
    ts_csv_eq_col = list(ts_csv['EQ'].values)
    ts_csv_ts_col = list(ts_csv[' time_end'].values)
    
    #   create a list of headers 
    frsh_headers = ['time_yrs', 'eq']
    for a in range(len(eq_0_ts_lst)):
        frsh_headers.append(eq_0_ts_lst[a][1] + '_frsh_inl')
        frsh_headers.append(eq_0_ts_lst[a][1] + '_frsh_shlf')
    frsh_headers.append('avg_frsh_inl')
    frsh_headers.append('avg_frsh_shlf')
    
    df_frsh_pct = pd.DataFrame(columns = frsh_headers)
    dict_frsh_pct = {}
    dict_frsh_pct['eq'] = ts_csv_eq_col
    dict_frsh_pct['time_yrs'] = ts_csv_ts_col    
    frsh_inl_arr = np.zeros((len(eq_0_ts_lst), len(ts_csv_eq_col)))
    frsh_shlf_arr = np.zeros((len(eq_0_ts_lst), len(ts_csv_eq_col)))
    idx_lay = 0
    
    frsh_inl_all_lst, frsh_shlf_all_lst, headers_lst = [], [], []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('scenario_')[-1].split('_')[0]
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)
            
            #   create empty lists to store the fractions of fresh water per time step
            inl_lst, shlf_lst = [], []
            
            #   set the column headers
            hd_frsh_inl = simu_nr_str + '_frsh_inl'
            hd_frsh_shlf = simu_nr_str + '_frsh_shlf'
            ts_step = 100.   #   step in years
            tot_time = 100.
            
            headers_lst.append([hd_frsh_inl, hd_frsh_shlf])
            
            #    first read all the output for the EQ = 0, if the time step doesnt exist just repeat the last one in the model simulation
            for b in range(1, int(max_ts_0 / ts_step) + 1):
                try:
                    cur_row = df[(df['EQ'] == 0) & (df[' time_end'] == b * ts_step)]
                    frsh_inl = cur_row[' frsh_inl_pct'].values[0]
                    frsh_shlf = cur_row[' frsh_shelf_pct'].values[0]         
                    inl_lst.append(frsh_inl)
                    shlf_lst.append(frsh_shlf)
                    tot_time += ts_step
                except IndexError:
                    try:
                        inl_lst.append(frsh_inl)
                        shlf_lst.append(frsh_shlf)
                    except UnboundLocalError:
                        inl_lst.append(-1)
                        shlf_lst.append(-1)
                
            for c in range(1, len(eq_lst)):
                eq_num = int(eq_lst[c].split('_')[-1])
                ts_dur = sp_duration[c]
            
                for d in range(1, int(ts_dur / ts_step) + 1):
                    #print(d, tot_time, tot_time + (d * ts_step))
                    #tot_time += ts_step
                    try:
                        cur_row = df[(df['EQ'] == eq_num) & (df[' time_end'] == tot_time)]# - (max_ts_0 - eq_len))]
                        frsh_inl = cur_row[' frsh_inl_pct'].values[0]
                        frsh_shlf = cur_row[' frsh_shelf_pct'].values[0]         
                        inl_lst.append(frsh_inl)
                        shlf_lst.append(frsh_shlf)
                        tot_time += ts_step
                    except IndexError:
                        try:
                            inl_lst.append(frsh_inl)
                            shlf_lst.append(frsh_shlf)
                        except UnboundLocalError:
                            inl_lst.append(-1)
                            shlf_lst.append(-1)	
                        
            frsh_inl_all_lst.append(inl_lst)
            frsh_shlf_all_lst.append(shlf_lst)
            
            #for d in range(len(frsh_inl_all_lst)):
            #   print(len(frsh_inl_all_lst[d]))
        
    inl_frsh_arr = np.array(frsh_inl_all_lst)
    shlf_frsh_arr = np.array(frsh_shlf_all_lst)        
        
    avg_inl_frsh_lst = list(np.mean(inl_frsh_arr, axis = 0))
    avg_shlf_frsh_lst = list(np.mean(shlf_frsh_arr, axis = 0))    
    
    avg_inl_frsh_lst = [round(i, 2) for i in avg_inl_frsh_lst]
    avg_shlf_frsh_lst = [round(i, 2) for i in avg_shlf_frsh_lst]

    eq_lst_to_dict, ts_lst_to_dict = [], []     
    for g in range(int(max_ts_0 / ts_step)):
        eq_lst_to_dict.append(0)
        ts_lst_to_dict.append(100. + ts_step * g)
    time = ts_lst_to_dict[-1] + ts_step
    for h in range(1, len(eq_lst)):
        eq_num = int(eq_lst[h].split('_')[-1])
        ts_dur = sp_duration[h]    
        for i in range(1, int(ts_dur / ts_step) + 1):
            eq_lst_to_dict.append(eq_num)
            ts_lst_to_dict.append(time)        
            time += ts_step

    dict_frsh_pct['eq'] = eq_lst_to_dict
    dict_frsh_pct['time_end'] = ts_lst_to_dict                
            
    #   fill the dictionary
    for f in range(len(headers_lst)):
        dict_frsh_pct[headers_lst[f][0]] = frsh_inl_all_lst[f]
        dict_frsh_pct[headers_lst[f][1]] = frsh_shlf_all_lst[f]
    
    dict_frsh_pct['avg_frsh_inl'] = avg_inl_frsh_lst
    dict_frsh_pct['avg_frsh_shlf'] = avg_shlf_frsh_lst    
        
    df_out = pd.DataFrame.from_dict(dict_frsh_pct, orient = 'index').T   
    df_out.to_csv(frsh_vol_csv)

    """         create the summary csv files, with % fresh at certain times and the convergence times         """
    #geo_csv_dir = os.path.join(sum_dir, '_geo_summary_' + coscat_id_str + '.csv')
    #geo_df = pd.read_csv(geo_csv_dir, index_col = False, header = 0)    
    #geo_headers = geo_df.columns.values
    sum_headers = ['coscat_id', 'cst_type', 'simu_id', 'eq_0_t_end', 't_cur_min_2000_inl_frsh', 't_cur_min_2000_shlf_frsh',\
                   't_cur_inl_frsh', 't_cur_shlf_frsh', 't_cur_plus_2000_inl_frsh', 't_cur_plus_2000_shlf_frsh',\
                   't_cur_plus_10000_inl_frsh', 't_cur_plus_10000_shlf_frsh', 't_cur_plus_20000_inl_frsh', 't_cur_plus_20000_shlf_frsh',
                   't_end_inl_frsh', 't_end_shlf_frsh']
    all_sc_lst = []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('scenario_')[-1].split('_')[0]
            print(simu_nr_str)
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)    
            try:                                      
                sc_lst = []
                #   read the specific rows
                eq_0_t_end = df[(df['EQ'] == 0)][' time_end'].values[-1]
                
                if t_cur - 2000 > 0:
                    t_cur_min_2000_inl_frsh = df[(df[' time_end'] == t_cur - 2000)][' frsh_inl_pct'].values[-1]
                    t_cur_min_2000_shlf_frsh = df[(df[' time_end'] == t_cur - 2000)][' frsh_shelf_pct'].values[-1]
                else:
                    t_cur_min_2000_inl_frsh = df[(df[' time_end'] == 100)][' frsh_inl_pct'].values[-1]
                    t_cur_min_2000_shlf_frsh = df[(df[' time_end'] == 100)][' frsh_shelf_pct'].values[-1]                
                
                t_cur_inl_frsh = df[(df[' time_end'] == t_cur)][' frsh_inl_pct'].values[-1]
                t_cur_shlf_frsh = df[(df[' time_end'] == t_cur)][' frsh_shelf_pct'].values[-1]                

                if t_cur + 2000 > eq_0_t_end:
                    t_cur_plus_2000_inl_frsh = df[(df[' time_end'] == eq_0_t_end)][' frsh_inl_pct'].values[-1]
                    t_cur_plus_2000_shlf_frsh = df[(df[' time_end'] == eq_0_t_end)][' frsh_shelf_pct'].values[-1]                           
                else:
                    t_cur_plus_2000_inl_frsh = df[(df[' time_end'] == t_cur + 2000)][' frsh_inl_pct'].values[-1]
                    t_cur_plus_2000_shlf_frsh = df[(df[' time_end'] == t_cur + 2000)][' frsh_shelf_pct'].values[-1]                        
                
                if t_cur + 10000 > eq_0_t_end:
                    t_cur_plus_10000_inl_frsh = df[(df[' time_end'] == eq_0_t_end)][' frsh_inl_pct'].values[-1]
                    t_cur_plus_10000_shlf_frsh = df[(df[' time_end'] == eq_0_t_end)][' frsh_shelf_pct'].values[-1]                           
                else:
                    t_cur_plus_10000_inl_frsh = df[(df[' time_end'] == t_cur + 10000)][' frsh_inl_pct'].values[-1]
                    t_cur_plus_10000_shlf_frsh = df[(df[' time_end'] == t_cur + 10000)][' frsh_shelf_pct'].values[-1]                  
                
                if t_cur + 20000 > eq_0_t_end:
                    t_cur_plus_20000_inl_frsh = df[(df[' time_end'] == eq_0_t_end)][' frsh_inl_pct'].values[-1]
                    t_cur_plus_20000_shlf_frsh = df[(df[' time_end'] == eq_0_t_end)][' frsh_shelf_pct'].values[-1]                           
                else:
                    t_cur_plus_20000_inl_frsh = df[(df[' time_end'] == t_cur + 20000)][' frsh_inl_pct'].values[-1]
                    t_cur_plus_20000_shlf_frsh = df[(df[' time_end'] == t_cur + 20000)][' frsh_shelf_pct'].values[-1]                       
                
                t_end_inl_frsh = df[(df[' time_end'] == eq_0_t_end)][' frsh_inl_pct'].values[-1]
                t_end_shlf_frsh = df[(df[' time_end'] == eq_0_t_end)][' frsh_shelf_pct'].values[-1] 
                
                sc_lst.append([coscat_id, cst_type, int(simu_nr_str), eq_0_t_end, t_cur_min_2000_inl_frsh, t_cur_min_2000_shlf_frsh,\
                               t_cur_inl_frsh, t_cur_shlf_frsh, t_cur_plus_2000_inl_frsh, t_cur_plus_2000_shlf_frsh,\
                               t_cur_plus_10000_inl_frsh, t_cur_plus_10000_shlf_frsh, t_cur_plus_20000_inl_frsh, t_cur_plus_20000_shlf_frsh,\
                               t_end_inl_frsh, t_end_shlf_frsh])
                all_sc_lst.append(sc_lst)
                
            except (IndexError,FileNotFoundError):
                print('Model simulation did not terminate')
                continue

    all_sc_lst = [i[0] for i in all_sc_lst]
    df_sum = pd.DataFrame(all_sc_lst, columns = sum_headers)
    df_sum.to_csv(tot_summary_csv)       




"""
in_dir = main_dir
summary_dir = summary_dir
out_dir_name = simu_name
all_res_csv_dir = all_res_csv_dir
ts_to_read = 500
eq_names = eq_lst
sp_dur = sp_duration
model_name = modelname
frsh_inl_csv = frsh_csv_dir
"""
            
    
def create_avg_netcdf_from_csv(in_dir, summary_dir, out_dir_name, all_res_csv_dir, frsh_inl_csv, ts_to_read, eq_names, sp_dur, model_name):
    
    #   find all subfolders before creating the summary folder
    subfolders = next(os.walk(in_dir))[1]
    
    if '_summary' in subfolders:
        subfolders.remove('_summary')

    sum_dir = out_dir_name
        
    #   first check if the output directory exists, if not create it
    if not os.path.exists(out_dir_name):
        os.makedirs(sum_dir, exist_ok = True)        

    # make directories for the 9 combinations of DEM and RCP
    os.makedirs(os.path.join(sum_dir, 'merit_26'), exist_ok = True)       
    os.makedirs(os.path.join(sum_dir, 'merit_45'), exist_ok = True)       
    os.makedirs(os.path.join(sum_dir, 'merit_85'), exist_ok = True)       
    os.makedirs(os.path.join(sum_dir, 'coastal_26'), exist_ok = True)       
    os.makedirs(os.path.join(sum_dir, 'coastal_45'), exist_ok = True)       
    os.makedirs(os.path.join(sum_dir, 'coastal_85'), exist_ok = True)       
    os.makedirs(os.path.join(sum_dir, 'gebco_26'), exist_ok = True)       
    os.makedirs(os.path.join(sum_dir, 'gebco_45'), exist_ok = True)       
    os.makedirs(os.path.join(sum_dir, 'gebco_85'), exist_ok = True)       

    #   select the different DEM and SLR scenarios
    sc_lst = [[[1, 4, 7, 10, 13, 16, 19, 22], 'merit_26'],\
              [[2, 5, 8, 11, 14, 17, 20, 23], 'merit_45'],\
              [[3, 6, 9, 12, 15, 18, 21, 24], 'merit_85'],\
              [[25, 28, 31, 34, 37, 40, 43, 46], 'coastal_26'],\
              [[26, 29, 32, 35, 38, 41, 44, 47], 'coastal_45'],\
              [[27, 30, 33, 36, 39, 42, 45, 48], 'coastal_85'],\
              [[49, 52, 55, 58, 61, 64, 67, 70], 'gebco_26'],\
              [[50, 53, 56, 59, 62, 65, 68, 71], 'gebco_45'],\
              [[51, 54, 57, 60, 63, 66, 69, 72], 'gebco_85']]

    #   loop through the list
    for i in range(len(sc_lst)):
        #   define the current dem-rcp combo scenario list and its name
        scs = sc_lst[i][0]
        dir_out = os.path.join(sum_dir, sc_lst[i][1])

        #   read the overall csv file and then go time step by time step    
        frsh_file = pd.read_csv(frsh_inl_csv)
        csv_file = pd.read_csv(all_res_csv_dir)    

        #   deal with the first stress period (this is the one where the simulation sometimes converges earlier than 10 000 years)
        #       1) the initial concentration and head profiles are the same for all the models, so just find the first simulation that converged
        #          to do that get all the scenario ids that converged and match it with the scs for the current dem-rcp combination
        all_sc_lst = frsh_file['sc_id'].values.tolist()
        sc_converged_lst = list(set(scs).intersection(all_sc_lst))
        if sc_converged_lst[0] < 10:
            sc_first = '0' + str(sc_converged_lst[0])
        else:
            sc_first = str(sc_converged_lst[0])

        #   create the directory name of the model path, read from the first subfolder and first EQ, the dimensions stay the same anyway
        model_ws = os.path.join(in_dir, 'sc_' + sc_first, eq_names[0])
        #   get the limits of the model domain
        dim_limits = get_model_info(model_ws, model_name + '_' + sc_lst[i][1].split('_')[0] + '_100m_10m')   
        x_st, bot_zoom, top_zoom = dim_limits[0], dim_limits[6], dim_limits[7]    
        
        #   the conc profile at ts = 0 is read manually, because it is the same for all the profiles anyway (all salt or all fresh)
        nc_name_0 = 'nc_' + str(0) + '.nc'
        nc_file_0 = os.path.join(summary_dir, '_nc', 'sc_' + sc_first, nc_name_0)
    
        #   load the nc file and get the head and concentration arrays
        nc_0 = xr.open_dataset(nc_file_0)
        conc_arr_0 = np.expand_dims(np.array(nc_0['solute concentration'].values), axis = 1)
        head_arr_0 = np.expand_dims(np.array(nc_0['heads'].values), axis = 1)      
        x_lst = list(nc_0.coords['x'].values)
        y_lst = list(nc_0.coords['y'].values)
    
        #   for the concentration, heads and cbc create a netcdf file (to save memory)
        xa_sum_0 = xr.Dataset(data_vars = {'solute concentration' : (('y', 'x'), conc_arr_0[:, 0, :]),
                                         'heads' : (('y', 'x'), head_arr_0[:, 0, :]),
                                         'fresh_inl' : (0.),
                                         'fresh_shelf' : (0.),
                                         'eq' : (eq_names[0]),
                                         'x_coast' : ((x_lst[0] - 0.05) - x_st),
                                         'bot_zoom' : (bot_zoom),
                                         'top_zoom' : (top_zoom)},
                               coords = {'x' : x_lst,
                                         'y' : y_lst})
        xa_sum_0 = xa_sum_0.assign_coords(time = 0)
        xa_name_0 = '_average_' + str(0) + '.nc'
        xa_sum_0.to_netcdf(os.path.join(dir_out, xa_name_0))      

        #   get a sorted list of all the netcdf output files
        ts_lst = []
        for path, subdirs, files in os.walk(os.path.join(summary_dir, '_nc', 'sc_' + sc_first)):
            for name in files:
                #print(os.path.join(path, name))  
                ts_lst.append(int(name.split('_')[-1].split('.')[0]))
                
        ts_lst = np.arange(0, 28500, 500).tolist() + np.arange(28050, 30550, 50).tolist()
        ts_lst = sorted(ts_lst)[1:]

        #   for each time step create an averaged 
        for ts in ts_lst:
            #   specify the name of the netcdf file to look for in the model folder directory
            conc_arr_lst, head_arr_lst, frsh_inl_lst, frsh_shelf_lst = [], [], [], []

            for sc in scs:
                if sc < 10:
                    sc_str = '0' + str(sc)
                else:
                    sc_str = str(sc)
                #   find the row where the time_end is equal to the total time step
                try:
                    frsh_inl_ts = csv_file.loc[csv_file['time_end'] == ts]['0' + sc_str + '_frsh_inl'].item()
                    frsh_shelf_ts = csv_file.loc[csv_file['time_end'] == ts]['0' + sc_str + '_frsh_shlf'].item()    
                    eq_cur = csv_file.loc[csv_file['time_end'] == ts]['eq'].item()    
                    #   define the path of the given netcdf file
                    nc_file = os.path.join(summary_dir, '_nc', 'sc_' + sc_str, 'nc_' + str(ts) + '.nc')
                    #   load the nc file and get the head and concentration arrays
                    last_nc = xr.open_dataset(nc_file)
                    conc_arr = np.expand_dims(np.array(last_nc['solute concentration'].values), axis = 1)
                    head_arr = np.expand_dims(np.array(last_nc['heads'].values), axis = 1)  
                    conc_arr_lst.append(conc_arr)
                    head_arr_lst.append(head_arr)
                    frsh_inl_lst.append(frsh_inl_ts)
                    frsh_shelf_lst.append(frsh_shelf_ts)      
                    x_lst = list(last_nc.coords['x'].values)
                    y_lst = list(last_nc.coords['y'].values)
                #   keyerror occures when a scenario did not converge - the values are not written in the final CSV file..
                except KeyError:
                    pass
  
            conc_arr_1 = conc_arr_lst[0]
            for c in range(1, len(conc_arr_lst)):
                conc_arr_1 = np.hstack((conc_arr_1, conc_arr_lst[c]))

            head_arr_1 = head_arr_lst[0]
            for c in range(1, len(head_arr_lst)):
                head_arr_1 = np.hstack((head_arr_1, head_arr_lst[c]))
    
            avg_conc = np.expand_dims(np.mean(conc_arr_1, axis = 1), axis = 1)
            avg_head = np.expand_dims(np.mean(head_arr_1, axis = 1), axis = 1)
            avg_frsh_inl = sum(frsh_inl_lst) / len(frsh_inl_lst)
            avg_frsh_shelf = sum(frsh_shelf_lst) / len(frsh_shelf_lst)

            #   for the concentration, heads and cbc create a netcdf file (to save memory)
            xa_sum = xr.Dataset(data_vars = {'solute concentration' : (('y', 'x'), avg_conc[:, 0, :]),
                                             'heads' : (('y', 'x'), avg_head[:, 0, :]),
                                             'fresh_inl' : (avg_frsh_inl),
                                             'fresh_shelf' : (avg_frsh_shelf),
                                             'eq' : (eq_cur),
                                             'x_coast' : ((x_lst[0] - 0.05) - x_st),
                                             'bot_zoom' : (bot_zoom),
                                             'top_zoom' : (top_zoom)},
                                coords = {'x' : x_lst,
                                          'y' : y_lst})
            xa_sum = xa_sum.assign_coords(time = 0)
            xa_name = '_average_' + str(ts) + '.nc'
            xa_sum.to_netcdf(os.path.join(dir_out, xa_name))      
        

"""
in_dir = main_dir
summary_dir = summary_dir
out_dir_name = simu_name
all_res_csv_dir = all_res_csv_dir
ts_to_read = 500
eq_names = equilibrium_names
eq_names_slr = equilibrium_names_slr
sp_dur = sp_duration
model_name = modelname
frsh_inl_csv = frsh_csv_dir
"""
            
    
def create_avg_netcdf_from_csv_SLR_scs(in_dir, summary_dir, out_dir_name, all_res_csv_dir, frsh_inl_csv, ts_to_read, eq_names, eq_names_slr, sp_dur, model_name):
    
    #   find all subfolders before creating the summary folder
    subfolders = next(os.walk(in_dir))[1]
    
    if '_summary' in subfolders:
        subfolders.remove('_summary')

    sum_dir = out_dir_name
        
    #   first check if the output directory exists, if not create it
    if not os.path.exists(out_dir_name):
        os.makedirs(sum_dir, exist_ok = True)        

    # make directories for the 9 combinations of DEM and RCP
    os.makedirs(os.path.join(sum_dir, 'merit'), exist_ok = True)         
    os.makedirs(os.path.join(sum_dir, 'coastal'), exist_ok = True)          
    os.makedirs(os.path.join(sum_dir, 'gebco'), exist_ok = True)           

    #   select the different DEM and SLR scenarios
    sc_lst = [[[1, 2, 3, 4, 5, 6, 7, 8], 'merit'],\
              [[9, 10, 11, 12, 13, 14, 15, 16], 'coastal'],\
              [[17, 18, 19, 20, 21, 22, 23, 24], 'gebco']]

    rcp_lst = ['RCP_26', 'RCP_45', 'RCP_85']

    #   loop through the list
    for i in range(len(sc_lst)):
        try:
            #   define the current dem-rcp combo scenario list and its name
            scs = sc_lst[i][0]
            dir_out = os.path.join(sum_dir, sc_lst[i][1])
            rcp = rcp_lst[i]
    
            #   read the overall csv file and then go time step by time step    
            frsh_file = pd.read_csv(frsh_inl_csv)
            csv_file = pd.read_csv(all_res_csv_dir)    
    
            #   deal with the first stress period (this is the one where the simulation sometimes converges earlier than 10 000 years)
            #       1) the initial concentration and head profiles are the same for all the models, so just find the first simulation that converged
            #          to do that get all the scenario ids that converged and match it with the scs for the current dem-rcp combination
            all_sc_lst = frsh_file['sc_id'].values.tolist()
            sc_converged_lst = list(set(scs).intersection(all_sc_lst))
            if sc_converged_lst[0] < 10:
                sc_first = '0' + str(sc_converged_lst[0])
            else:
                sc_first = str(sc_converged_lst[0])
    
            #   create the directory name of the model path, read from the first subfolder and first EQ, the dimensions stay the same anyway
            model_ws = os.path.join(in_dir, 'sc_' + sc_first, eq_names[0])
            
            #   check that the file exists
            model_name_test = model_name + '_' + sc_lst[i][1].split('_')[0] + '_050m_2m_CST_ZOOM'  + '.nam_swt'
            """
            if not os.path.exists(os.path.join(model_ws, model_name_test + '')):
                #import tarfile
                #   untar the folder
                #folder_name = os.path.join(in_dir, 'sc_' + sc_first + '.tgz')
                #folder_name = r'g:\Water_Nexus\_A4_models\_SLR_models\0003_SRM_004\sc_01.tgz'
                #tar = tarfile.open(folder_name, "r:gz")
                #tar.extractall(path = r'g:\Water_Nexus\_A4_models\_SLR_models\0003_SRM_004')
                #tar.close()

                #   get the limits of the model domain
                dim_limits = get_model_info(model_ws, model_name + '_' + sc_lst[i][1].split('_')[0] + '_050m_2m_CST_ZOOM')   
                x_st, bot_zoom, top_zoom = dim_limits[0], dim_limits[6], dim_limits[7]    

            else:
                #   get the limits of the model domain
                dim_limits = get_model_info(model_ws, model_name + '_' + sc_lst[i][1].split('_')[0] + '_050m_2m_CST_ZOOM')   
                x_st, bot_zoom, top_zoom = dim_limits[0], dim_limits[6], dim_limits[7]    
            """
            #   the conc profile at ts = 0 is read manually, because it is the same for all the profiles anyway (all salt or all fresh)
            nc_name_0 = 'nc_' + str(0) + '.nc'
            nc_file_0 = os.path.join(summary_dir, '_nc', 'sc_' + sc_first, nc_name_0)
        
            #   load the nc file and get the head and concentration arrays
            nc_0 = xr.open_dataset(nc_file_0)
            conc_arr_0 = np.expand_dims(np.array(nc_0['solute concentration'].values), axis = 1)
            head_arr_0 = np.expand_dims(np.array(nc_0['heads'].values), axis = 1)      
            x_lst = list(nc_0.coords['x'].values)
            y_lst = list(nc_0.coords['y'].values)
            
            x_st = x_lst[0]
            bot_zoom = y_lst[-1]
            top_zoom = y_lst[0]
        
            #   for the concentration, heads and cbc create a netcdf file (to save memory)
            xa_sum_0 = xr.Dataset(data_vars = {'solute concentration' : (('y', 'x'), conc_arr_0[:, 0, :]),
                                             'heads' : (('y', 'x'), head_arr_0[:, 0, :]),
                                             'fresh_inl' : (0.),
                                             'fresh_shelf' : (0.),
                                             'eq' : (eq_names[0]),
                                             'x_coast' : ((x_lst[0] - 0.05) - x_st),
                                             'bot_zoom' : (bot_zoom),
                                             'top_zoom' : (top_zoom)},
                                   coords = {'x' : x_lst,
                                             'y' : y_lst})
            xa_sum_0 = xa_sum_0.assign_coords(time = 0)
            xa_name_0 = '_average_' + str(0) + '.nc'
            xa_sum_0.to_netcdf(os.path.join(dir_out, xa_name_0))      
    
            #   get a sorted list of all the netcdf output files
            ts_lst = []
            for path, subdirs, files in os.walk(os.path.join(summary_dir, '_nc', 'sc_' + sc_first)):
                for name in files:
                    #print(os.path.join(path, name))  
                    ts_lst.append(int(name.split('_')[-1].split('.')[0]))
                    
            ts_lst = np.arange(0, 5500, 500).tolist()
            ts_lst = sorted(ts_lst)[1:]
    
            #   for each time step create an averaged 
            for ts in ts_lst:
                #   specify the name of the netcdf file to look for in the model folder directory
                conc_arr_lst, head_arr_lst, frsh_inl_lst, frsh_shelf_lst = [], [], [], []
    
                for sc in scs:
                    if sc < 10:
                        sc_str = '0' + str(sc)
                    else:
                        sc_str = str(sc)
                    #   find the row where the time_end is equal to the total time step
                    try:
                        frsh_inl_ts = csv_file.loc[csv_file['time_end'] == ts]['0' + sc_str + '_frsh_inl'].item()
                        frsh_shelf_ts = csv_file.loc[csv_file['time_end'] == ts]['0' + sc_str + '_frsh_shlf'].item()    
                        eq_cur = csv_file.loc[csv_file['time_end'] == ts]['eq'].item()    
                        #   define the path of the given netcdf file
                        nc_file = os.path.join(summary_dir, '_nc', 'sc_' + sc_str, 'nc_' + str(ts) + '.nc')
                        if os.path.exists(nc_file):
                            #   load the nc file and get the head and concentration arrays
                            last_nc = xr.open_dataset(nc_file)
                            conc_arr = np.expand_dims(np.array(last_nc['salinity'].values), axis = 1)
                            head_arr = np.expand_dims(np.array(last_nc['heads'].values), axis = 1)  
                            conc_arr_lst.append(conc_arr)
                            head_arr_lst.append(head_arr)
                            frsh_inl_lst.append(frsh_inl_ts)
                            frsh_shelf_lst.append(frsh_shelf_ts)      
                            x_lst = list(last_nc.coords['x'].values)
                            y_lst = list(last_nc.coords['y'].values)
                    #   keyerror occures when a scenario did not converge - the values are not written in the final CSV file..
                    except KeyError:
                        pass
      
                if len(conc_arr_lst) != 0:
                    conc_arr_1 = conc_arr_lst[0]
                    for c in range(1, len(conc_arr_lst)):
                        conc_arr_1 = np.hstack((conc_arr_1, conc_arr_lst[c]))
        
                    head_arr_1 = head_arr_lst[0]
                    for c in range(1, len(head_arr_lst)):
                        head_arr_1 = np.hstack((head_arr_1, head_arr_lst[c]))
        
                avg_conc = np.expand_dims(np.mean(conc_arr_1, axis = 1), axis = 1)
                avg_head = np.expand_dims(np.mean(head_arr_1, axis = 1), axis = 1)
                avg_frsh_inl = sum(frsh_inl_lst) / len(frsh_inl_lst)
                avg_frsh_shelf = sum(frsh_shelf_lst) / len(frsh_shelf_lst)
    
                #   for the concentration, heads and cbc create a netcdf file (to save memory)
                xa_sum = xr.Dataset(data_vars = {'salinity' : (('y', 'x'), avg_conc[:, 0, :]),
                                                 'heads' : (('y', 'x'), avg_head[:, 0, :]),
                                                 'fresh_inl' : (avg_frsh_inl),
                                                 'fresh_shelf' : (avg_frsh_shelf),
                                                 'eq' : (eq_cur),
                                                 'x_coast' : ((x_lst[0] - 0.05) - x_st),
                                                 'bot_zoom' : (bot_zoom),
                                                 'top_zoom' : (top_zoom)},
                                    coords = {'x' : x_lst,
                                              'y' : y_lst})
                xa_sum = xa_sum.assign_coords(time = 0)
                xa_name = '_average_' + str(ts) + '.nc'
                xa_sum.to_netcdf(os.path.join(dir_out, xa_name))      
            
            #   now create average files for the RCP scenarios
            ts_slr_lst = np.arange(5005, 5205, 5).tolist()
            rcp_lst = ['RCP_26', 'RCP_45', 'RCP_85']
    
            eq_names_slr_to_nc = []
            for h in range(len(eq_names_slr)):
                for k in range(2):
                    eq_names_slr_to_nc.append(eq_names_slr[h])
    
            for rcp in rcp_lst:
                #   for each time step create an averaged 
                for ts in ts_slr_lst:
                    #   specify the name of the netcdf file to look for in the model folder directory
                    conc_arr_lst, head_arr_lst, frsh_inl_lst, frsh_shelf_lst = [], [], [], []
        
                    for sc in scs:
                        if sc < 10:
                            sc_str = '0' + str(sc)
                        else:
                            sc_str = str(sc)
                        #   find the row where the time_end is equal to the total time step
                        try:
                            frsh_inl_ts = csv_file[(csv_file['time_end'] == ts) & (csv_file['eq'] == rcp + '_' + eq_names_slr_to_nc[ts_slr_lst.index(ts)])]['0' + sc_str + '_frsh_inl'].item()          
                            frsh_shelf_ts = csv_file[(csv_file['time_end'] == ts) & (csv_file['eq'] == rcp + '_' + eq_names_slr_to_nc[ts_slr_lst.index(ts)])]['0' + sc_str + '_frsh_shlf'].item()      
                            eq_cur = rcp + '_' + eq_names_slr_to_nc[ts_slr_lst.index(ts)]
                            #   define the path of the given netcdf file
                            nc_file = os.path.join(summary_dir, '_nc', 'sc_' + sc_str, 'nc_' + rcp + '_' + str(ts) + '.nc')
                            #   load the nc file and get the head and concentration arrays
                            if os.path.exists(nc_file):
                                last_nc = xr.open_dataset(nc_file)
                                conc_arr = np.expand_dims(np.array(last_nc['solute concentration'].values), axis = 1)
                                head_arr = np.expand_dims(np.array(last_nc['heads'].values), axis = 1)  
                                conc_arr_lst.append(conc_arr)
                                head_arr_lst.append(head_arr)
                                frsh_inl_lst.append(frsh_inl_ts)
                                frsh_shelf_lst.append(frsh_shelf_ts)      
                                x_lst = list(last_nc.coords['x'].values)
                                y_lst = list(last_nc.coords['y'].values)
                        #   keyerror occures when a scenario did not converge - the values are not written in the final CSV file..
                        except KeyError:
                            pass
          
                    if len(conc_arr_lst) != 0:
                        conc_arr_1 = conc_arr_lst[0]
                        for c in range(1, len(conc_arr_lst)):
                            conc_arr_1 = np.hstack((conc_arr_1, conc_arr_lst[c]))
            
                        head_arr_1 = head_arr_lst[0]
                        for c in range(1, len(head_arr_lst)):
                            head_arr_1 = np.hstack((head_arr_1, head_arr_lst[c]))
            
                    avg_conc = np.expand_dims(np.mean(conc_arr_1, axis = 1), axis = 1)
                    avg_head = np.expand_dims(np.mean(head_arr_1, axis = 1), axis = 1)
                    avg_frsh_inl = sum(frsh_inl_lst) / len(frsh_inl_lst)
                    avg_frsh_shelf = sum(frsh_shelf_lst) / len(frsh_shelf_lst)
        
                    #   for the concentration, heads and cbc create a netcdf file (to save memory)
                    xa_sum = xr.Dataset(data_vars = {'salinity' : (('y', 'x'), avg_conc[:, 0, :]),
                                                     'heads' : (('y', 'x'), avg_head[:, 0, :]),
                                                     'fresh_inl' : (avg_frsh_inl),
                                                     'fresh_shelf' : (avg_frsh_shelf),
                                                     'eq' : (eq_cur),
                                                     'x_coast' : ((x_lst[0] - 0.05) - x_st),
                                                     'bot_zoom' : (bot_zoom),
                                                     'top_zoom' : (top_zoom)},
                                        coords = {'x' : x_lst,
                                                  'y' : y_lst})
                    xa_sum = xa_sum.assign_coords(time = 0)
                    xa_name = '_average_' + rcp + '_' + str(ts) + '.nc'
                    xa_sum.to_netcdf(os.path.join(dir_out, xa_name))     

        except IndexError:
            print('No simulations for this DEM scenario converged..')
    
    
    
    
    

"""
sum_dir = summary_dir
simu_name = model_name_basis
eq_lst = equilibrium_names
eq_slr_lst = equilibrium_names_slr
"""
def create_fin_coscat_run_csv_files_SLR_scs_PBL(coscat_id_str, main_dir, sum_dir, eq_lst, eq_slr_lst, sp_duration, simu_name):

    #   define the summary dir and the model simulation dir
    #sum_dir = os.path.join(main_dir, '_summary_' + coscat_id_str + '_' + cst_type + simu_name)
    #model_simu_dir = os.path.join(main_dir, 'coscat_' + coscat_id_str + '_' + cst_type + '_' + simu_name)

    #   define the output name of the final csv file
    frsh_vol_csv = os.path.join(sum_dir, '_frsh_vol_summary_' + 'COSCAT_' + simu_name + '.csv')
    tot_summary_inl_csv = os.path.join(sum_dir, '_frsh_INL_' + 'COSCAT_' + simu_name + '.csv')
    tot_summary_shlf_csv = os.path.join(sum_dir, '_frsh_SHLF_' + 'COSCAT_' + simu_name + '.csv')
    
    #   find all csv files in the folder 
    eq_0_ts_lst = []
    csv_dir = os.path.join(sum_dir, '_csv_results')
    #print(csv_dir)
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            #print(filename)
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            #simu_nr_int = int(simu_nr_str)        
            #   get all the time steps for eq_0
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)
            col_names = list(df.columns.values)
            eq_0_rows = df.loc[df['EQ'] == '3000BC_to_2000BC']

            #   get the last value of the time_end column, only if the model simulation converged, if not append a -1
            try:
                if len(df.loc[df['EQ'] == '1000AD_to_2000AD'][' time_end'].values) == 2:
                    eq_0_last_ts = eq_0_rows.tail(1).values[0][2] 
                    eq_0_ts_lst.append([eq_0_last_ts, simu_nr_str])
            except IndexError:
                #print('The csv file is empty (most probably)')
                continue

    #   get the index of the tuple with the maximum time step
    max_idx = [i for i, tupl in enumerate(eq_0_ts_lst) if tupl[0] == max([x[0] for x in eq_0_ts_lst])][0]
    max_ts_0 = eq_0_ts_lst[max_idx][0]

    #   column containing the time steps and EQ numbers, taken form the simulation with the longest EQ_0
    ts_csv_dir = os.path.join(csv_dir, 'sc_' + eq_0_ts_lst[max_idx][1] + '_summary_EQs_output.csv')
    ts_csv = pd.read_csv(ts_csv_dir, index_col = False, header = 0)
    ts_csv_eq_col = list(ts_csv['EQ'].values)
    ts_csv_ts_col = list(ts_csv[' time_end'].values)
    
    #   create a list of headers 
    frsh_headers = ['eq']
    for a in range(len(eq_0_ts_lst)):
        frsh_headers.append(eq_0_ts_lst[a][1] + '_frsh_inl')
        frsh_headers.append(eq_0_ts_lst[a][1] + '_frsh_shlf')
    frsh_headers.append('avg_frsh_inl')
    frsh_headers.append('avg_frsh_shlf')
    
    df_frsh_pct = pd.DataFrame(columns = frsh_headers)
    dict_frsh_pct = {}
    dict_frsh_pct['eq'] = ts_csv_eq_col
    frsh_inl_arr = np.zeros((len(eq_0_ts_lst), len(ts_csv_eq_col)))
    frsh_shlf_arr = np.zeros((len(eq_0_ts_lst), len(ts_csv_eq_col)))
    idx_lay = 0
    
    frsh_inl_all_lst, frsh_shlf_all_lst, headers_lst = [], [], []
    max_len = 0
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)                  
            if len(df.loc[df['EQ'] == '1000AD_to_2000AD'][' time_end'].values) == 2:
                if len(df.loc[df['EQ'] == '3000BC_to_2000BC'][' time_end'].values) == 2:
                    inl_lst = df[' frsh_inl_pct'].values.tolist()
                    shlf_lst = df[' frsh_shelf_pct'].values.tolist()                    
                    if len(inl_lst) > max_len:
                        max_len = len(inl_lst)
                    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)
            
            if len(df.loc[df['EQ'] == '1000AD_to_2000AD'][' time_end'].values) == 2:
                if len(df.loc[df['EQ'] == '3000BC_to_2000BC'][' time_end'].values) == 2:
                    inl_lst = df[' frsh_inl_pct'].values.tolist()
                    shlf_lst = df[' frsh_shelf_pct'].values.tolist()                    
                    
                if len(inl_lst) < max_len:
                    for f in range(max_len - len(inl_lst)):
                        inl_lst.append(np.nan)
                        shlf_lst.append(np.nan)
        
                #   set the column headers
                hd_frsh_inl = simu_nr_str + '_frsh_inl'
                hd_frsh_shlf = simu_nr_str + '_frsh_shlf'
                
                inl_lst.insert(0, inl_lst[0])
                shlf_lst.insert(0, shlf_lst[0])
                
                headers_lst.append([hd_frsh_inl, hd_frsh_shlf])                            
                frsh_inl_all_lst.append(inl_lst)
                frsh_shlf_all_lst.append(shlf_lst)
        
    inl_frsh_arr = np.array(frsh_inl_all_lst)
    shlf_frsh_arr = np.array(frsh_shlf_all_lst)   
        
    avg_inl_frsh_lst = list(np.nanmean(inl_frsh_arr, axis = 0))
    avg_shlf_frsh_lst = list(np.nanmean(shlf_frsh_arr, axis = 0))    
    
    avg_inl_frsh_lst = [round(i, 2) for i in avg_inl_frsh_lst]
    avg_shlf_frsh_lst = [round(i, 2) for i in avg_shlf_frsh_lst]


    
    """
    for g in range(int(5000 / 500) + 1):
        eq_lst_to_dict.append('3000BC_to_2000BC')
        ts_lst_to_dict.append(500. * g)
    time = ts_lst_to_dict[-1] + 500
    
    for h in range(1, len(eq_lst)):
        eq_name = eq_lst[h]
        ts_dur = sp_duration[h]    
        ts_step = ts_dur / 2
        for i in range(1, int(ts_dur / (ts_dur / 2) + 1)):
            eq_lst_to_dict.append(eq_name)
            ts_lst_to_dict.append(time)        
            time += ts_step
        if eq_name == '3000BC_to_2000BC':
            time -= 450
    """


    eq_lst_col = ['3000BC_to_2000BC']
    for e in range(len(eq_lst)):
       for g in range(2):
           eq_lst_col.append(eq_lst[e])
    
    eq_lst_to_dict, ts_lst_to_dict = [i for i in eq_lst_col], []
    
    ts_lst_to_dict = np.arange(0, 5500, 500).tolist() + np.arange(5005, 5205, 5).tolist()
    #ts_lst_to_dict = sorted(ts_lst_to_dict)[1:]
    eq_slr_lst_repeat = list(np.repeat(eq_slr_lst, 2))
    ts_last = ts_lst_to_dict[-1]
    for slr_sc in ['RCP_26_', 'RCP_45_', 'RCP_85_']:     
        eq_lst_to_dict = eq_lst_to_dict + [slr_sc + i for i in eq_slr_lst_repeat]
        ts_lst_to_dict = ts_lst_to_dict + np.arange(5005, 5205, 5).tolist()

    dict_frsh_pct['eq'] = eq_lst_to_dict
    dict_frsh_pct['time_end'] = ts_lst_to_dict                
            
    #   fill the dictionary
    for f in range(len(headers_lst)):
        dict_frsh_pct[headers_lst[f][0]] = frsh_inl_all_lst[f]
        dict_frsh_pct[headers_lst[f][1]] = frsh_shlf_all_lst[f]
    
    dict_frsh_pct['avg_frsh_inl'] = avg_inl_frsh_lst
    dict_frsh_pct['avg_frsh_shlf'] = avg_shlf_frsh_lst    
        
    df_out = pd.DataFrame.from_dict(dict_frsh_pct, orient = 'index').T   
    df_out.to_csv(frsh_vol_csv)

    """         create the summary csv files, with % fresh at certain times and the convergence times         """
    #geo_csv_dir = os.path.join(sum_dir, '_geo_summary_' + coscat_id_str + '.csv')
    #geo_df = pd.read_csv(geo_csv_dir, index_col = False, header = 0)    
    #geo_headers = geo_df.columns.values
    
    """
    sum_headers = ['coscat_id', 'simu_id', 'eq_0_t_end', 'eq_0_t_end_inl_frsh', 'eq_0_t_end_shlf_frsh',\
                   'eq_21_t_end_inl_frsh', 'eq_21_t_end_shlf_frsh', 'eq_31_t_end_inl_frsh', 'eq_31_t_end_shlf_frsh',\
                   'eq_21_CT_t_end', 'eq_21_CT_t_end_inl_frsh', 'eq_21_CT_t_end_shlf_frsh', 'eq_31_CT_t_end', 'eq_31_CT_t_end_inl_frsh', 'eq_31_CT_t_end_shlf_frsh']
    """
    sum_headers_inl = ['coscat_id', 'srm_id', 'sc_id', 'eq_0_t_end', 'eq_5000_BP_inl_frsh', 'eq_4500_BP_inl_frsh',\
                       'eq_4000_BP_inl_frsh', 'eq_3500_BP_inl_frsh','eq_3000_BP_inl_frsh', 'eq_2500_BP_inl_frsh',\
                       'eq_2000_BP_inl_frsh', 'eq_1500_BP_inl_frsh','eq_1000_BP_inl_frsh', 'eq_0500_BP_inl_frsh',\
                       'eq_0000_BP_inl_frsh', 'RCP_26_0005_AP_inl_frsh', 'RCP_26_0010_AP_inl_frsh', 'RCP_26_0015_AP_inl_frsh',\
                       'RCP_26_0020_AP_inl_frsh', 'RCP_26_0025_AP_inl_frsh', 'RCP_26_0030_AP_inl_frsh', 'RCP_26_0035_AP_inl_frsh',\
                       'RCP_26_0040_AP_inl_frsh', 'RCP_26_0045_AP_inl_frsh', 'RCP_26_0050_AP_inl_frsh', 'RCP_26_0055_AP_inl_frsh',\
                       'RCP_26_0060_AP_inl_frsh', 'RCP_26_0065_AP_inl_frsh', 'RCP_26_0070_AP_inl_frsh', 'RCP_26_0075_AP_inl_frsh',\
                       'RCP_26_0080_AP_inl_frsh', 'RCP_26_0085_AP_inl_frsh', 'RCP_26_0090_AP_inl_frsh', 'RCP_26_0095_AP_inl_frsh',\
                       'RCP_26_0100_AP_inl_frsh', 'RCP_26_0105_AP_inl_frsh', 'RCP_26_0110_AP_inl_frsh', 'RCP_26_0115_AP_inl_frsh',\
                       'RCP_26_0120_AP_inl_frsh', 'RCP_26_0125_AP_inl_frsh', 'RCP_26_0130_AP_inl_frsh', 'RCP_26_0135_AP_inl_frsh',\
                       'RCP_26_0140_AP_inl_frsh', 'RCP_26_0145_AP_inl_frsh', 'RCP_26_0150_AP_inl_frsh', 'RCP_26_0155_AP_inl_frsh',\
                       'RCP_26_0160_AP_inl_frsh', 'RCP_26_0165_AP_inl_frsh', 'RCP_26_0170_AP_inl_frsh', 'RCP_26_0175_AP_inl_frsh',\
                       'RCP_26_0180_AP_inl_frsh', 'RCP_26_0185_AP_inl_frsh', 'RCP_26_0190_AP_inl_frsh', 'RCP_26_0195_AP_inl_frsh', 'RCP_26_0200_AP_inl_frsh',\
                       'RCP_45_0005_AP_inl_frsh', 'RCP_45_0010_AP_inl_frsh', 'RCP_45_0015_AP_inl_frsh',\
                       'RCP_45_0020_AP_inl_frsh', 'RCP_45_0025_AP_inl_frsh', 'RCP_45_0030_AP_inl_frsh', 'RCP_45_0035_AP_inl_frsh',\
                       'RCP_45_0040_AP_inl_frsh', 'RCP_45_0045_AP_inl_frsh', 'RCP_45_0050_AP_inl_frsh', 'RCP_45_0055_AP_inl_frsh',\
                       'RCP_45_0060_AP_inl_frsh', 'RCP_45_0065_AP_inl_frsh', 'RCP_45_0070_AP_inl_frsh', 'RCP_45_0075_AP_inl_frsh',\
                       'RCP_45_0080_AP_inl_frsh', 'RCP_45_0085_AP_inl_frsh', 'RCP_45_0090_AP_inl_frsh', 'RCP_45_0095_AP_inl_frsh',\
                       'RCP_45_0100_AP_inl_frsh', 'RCP_45_0105_AP_inl_frsh', 'RCP_45_0110_AP_inl_frsh', 'RCP_45_0115_AP_inl_frsh',\
                       'RCP_45_0120_AP_inl_frsh', 'RCP_45_0125_AP_inl_frsh', 'RCP_45_0130_AP_inl_frsh', 'RCP_45_0135_AP_inl_frsh',\
                       'RCP_45_0140_AP_inl_frsh', 'RCP_45_0145_AP_inl_frsh', 'RCP_45_0150_AP_inl_frsh', 'RCP_45_0155_AP_inl_frsh',\
                       'RCP_45_0160_AP_inl_frsh', 'RCP_45_0165_AP_inl_frsh', 'RCP_45_0170_AP_inl_frsh', 'RCP_45_0175_AP_inl_frsh',\
                       'RCP_45_0180_AP_inl_frsh', 'RCP_45_0185_AP_inl_frsh', 'RCP_45_0190_AP_inl_frsh', 'RCP_45_0195_AP_inl_frsh', 'RCP_45_0200_AP_inl_frsh',\
                       'RCP_85_0005_AP_inl_frsh', 'RCP_85_0010_AP_inl_frsh', 'RCP_85_0015_AP_inl_frsh',\
                       'RCP_85_0020_AP_inl_frsh', 'RCP_85_0025_AP_inl_frsh', 'RCP_85_0030_AP_inl_frsh', 'RCP_85_0035_AP_inl_frsh',\
                       'RCP_85_0040_AP_inl_frsh', 'RCP_85_0045_AP_inl_frsh', 'RCP_85_0050_AP_inl_frsh', 'RCP_85_0055_AP_inl_frsh',\
                       'RCP_85_0060_AP_inl_frsh', 'RCP_85_0065_AP_inl_frsh', 'RCP_85_0070_AP_inl_frsh', 'RCP_85_0075_AP_inl_frsh',\
                       'RCP_85_0080_AP_inl_frsh', 'RCP_85_0085_AP_inl_frsh', 'RCP_85_0090_AP_inl_frsh', 'RCP_85_0095_AP_inl_frsh',\
                       'RCP_85_0100_AP_inl_frsh', 'RCP_85_0105_AP_inl_frsh', 'RCP_85_0110_AP_inl_frsh', 'RCP_85_0115_AP_inl_frsh',\
                       'RCP_85_0120_AP_inl_frsh', 'RCP_85_0125_AP_inl_frsh', 'RCP_85_0130_AP_inl_frsh', 'RCP_85_0135_AP_inl_frsh',\
                       'RCP_85_0140_AP_inl_frsh', 'RCP_85_0145_AP_inl_frsh', 'RCP_85_0150_AP_inl_frsh', 'RCP_85_0155_AP_inl_frsh',\
                       'RCP_85_0160_AP_inl_frsh', 'RCP_85_0165_AP_inl_frsh', 'RCP_85_0170_AP_inl_frsh', 'RCP_85_0175_AP_inl_frsh',\
                       'RCP_85_0180_AP_inl_frsh', 'RCP_85_0185_AP_inl_frsh', 'RCP_85_0190_AP_inl_frsh', 'RCP_85_0195_AP_inl_frsh', 'RCP_85_0200_AP_inl_frsh']
    all_sc_lst = []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            #print(simu_nr_str)
            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)    
            if len(df.loc[df['EQ'] == '3000BC_to_2000BC'][' time_end'].values) == 2:
                
                inl_lst = df[' frsh_inl_pct'].values.tolist()                  
                if len(inl_lst) < max_len:
                    for f in range(max_len - len(inl_lst)):
                        inl_lst.append(np.nan)

                try:                      
                    sc_lst = []
                    #   read the specific rows
                    eq_0_t_end = inl_lst[0]
                    eq_5000_BP_inl_frsh = inl_lst[0]
                    eq_4500_BP_inl_frsh = inl_lst[0]           
                    eq_4000_BP_inl_frsh = inl_lst[1]
                    eq_3500_BP_inl_frsh = inl_lst[2]
                    eq_3000_BP_inl_frsh = inl_lst[3]   
                    eq_2500_BP_inl_frsh = inl_lst[4]              
                    eq_2000_BP_inl_frsh = inl_lst[5]
                    eq_1500_BP_inl_frsh = inl_lst[6]
                    eq_1000_BP_inl_frsh = inl_lst[7] 
                    eq_0500_BP_inl_frsh = inl_lst[8]
                    eq_0000_BP_inl_frsh = inl_lst[9]    
                                        
                    RCP_26_005_AP_inl_frsh = inl_lst[10]
                    RCP_26_010_AP_inl_frsh = inl_lst[11]
                    RCP_26_015_AP_inl_frsh = inl_lst[12]       
                    RCP_26_020_AP_inl_frsh = inl_lst[13]  
                    RCP_26_025_AP_inl_frsh = inl_lst[14]                
                    RCP_26_030_AP_inl_frsh = inl_lst[15]
                    RCP_26_035_AP_inl_frsh = inl_lst[16]
                    RCP_26_040_AP_inl_frsh = inl_lst[17]     
                    RCP_26_045_AP_inl_frsh = inl_lst[18]
                    RCP_26_050_AP_inl_frsh = inl_lst[19]                 
                    RCP_26_055_AP_inl_frsh = inl_lst[20]
                    RCP_26_060_AP_inl_frsh = inl_lst[21]    
                    RCP_26_065_AP_inl_frsh = inl_lst[22]  
                    RCP_26_070_AP_inl_frsh = inl_lst[23]
                    RCP_26_075_AP_inl_frsh = inl_lst[24]                 
                    RCP_26_080_AP_inl_frsh = inl_lst[25]   
                    RCP_26_085_AP_inl_frsh = inl_lst[26]   
                    RCP_26_090_AP_inl_frsh = inl_lst[27]            
                    RCP_26_095_AP_inl_frsh = inl_lst[28]
                    RCP_26_100_AP_inl_frsh = inl_lst[29]                          
                    RCP_26_105_AP_inl_frsh = inl_lst[30]
                    RCP_26_110_AP_inl_frsh = inl_lst[31]
                    RCP_26_115_AP_inl_frsh = inl_lst[32]              
                    RCP_26_120_AP_inl_frsh = inl_lst[33]
                    RCP_26_125_AP_inl_frsh = inl_lst[34]          
                    RCP_26_130_AP_inl_frsh = inl_lst[35]  
                    RCP_26_135_AP_inl_frsh = inl_lst[36]
                    RCP_26_140_AP_inl_frsh = inl_lst[37]           
                    RCP_26_145_AP_inl_frsh = inl_lst[38]  
                    RCP_26_150_AP_inl_frsh = inl_lst[39]                  
                    RCP_26_155_AP_inl_frsh = inl_lst[40]
                    RCP_26_160_AP_inl_frsh = inl_lst[41]
                    RCP_26_165_AP_inl_frsh = inl_lst[42]      
                    RCP_26_170_AP_inl_frsh = inl_lst[43]
                    RCP_26_175_AP_inl_frsh = inl_lst[44]               
                    RCP_26_180_AP_inl_frsh = inl_lst[45] 
                    RCP_26_185_AP_inl_frsh = inl_lst[46]     
                    RCP_26_190_AP_inl_frsh = inl_lst[47]              
                    RCP_26_195_AP_inl_frsh = inl_lst[48]
                    RCP_26_200_AP_inl_frsh = inl_lst[49]

                    RCP_45_005_AP_inl_frsh = inl_lst[50]
                    RCP_45_010_AP_inl_frsh = inl_lst[51]
                    RCP_45_015_AP_inl_frsh = inl_lst[52]       
                    RCP_45_020_AP_inl_frsh = inl_lst[53]  
                    RCP_45_025_AP_inl_frsh = inl_lst[54]                
                    RCP_45_030_AP_inl_frsh = inl_lst[55]
                    RCP_45_035_AP_inl_frsh = inl_lst[56]
                    RCP_45_040_AP_inl_frsh = inl_lst[57]     
                    RCP_45_045_AP_inl_frsh = inl_lst[58]
                    RCP_45_050_AP_inl_frsh = inl_lst[59]                 
                    RCP_45_055_AP_inl_frsh = inl_lst[60]
                    RCP_45_060_AP_inl_frsh = inl_lst[61]    
                    RCP_45_065_AP_inl_frsh = inl_lst[62]  
                    RCP_45_070_AP_inl_frsh = inl_lst[63]
                    RCP_45_075_AP_inl_frsh = inl_lst[64]                 
                    RCP_45_080_AP_inl_frsh = inl_lst[65]   
                    RCP_45_085_AP_inl_frsh = inl_lst[66]   
                    RCP_45_090_AP_inl_frsh = inl_lst[67]            
                    RCP_45_095_AP_inl_frsh = inl_lst[68]
                    RCP_45_100_AP_inl_frsh = inl_lst[69]                          
                    RCP_45_105_AP_inl_frsh = inl_lst[70]
                    RCP_45_110_AP_inl_frsh = inl_lst[71]
                    RCP_45_115_AP_inl_frsh = inl_lst[72]              
                    RCP_45_120_AP_inl_frsh = inl_lst[73]
                    RCP_45_125_AP_inl_frsh = inl_lst[74]          
                    RCP_45_130_AP_inl_frsh = inl_lst[75]  
                    RCP_45_135_AP_inl_frsh = inl_lst[76]
                    RCP_45_140_AP_inl_frsh = inl_lst[77]           
                    RCP_45_145_AP_inl_frsh = inl_lst[78]  
                    RCP_45_150_AP_inl_frsh = inl_lst[79]                  
                    RCP_45_155_AP_inl_frsh = inl_lst[80]
                    RCP_45_160_AP_inl_frsh = inl_lst[81]
                    RCP_45_165_AP_inl_frsh = inl_lst[82]      
                    RCP_45_170_AP_inl_frsh = inl_lst[83]
                    RCP_45_175_AP_inl_frsh = inl_lst[84]               
                    RCP_45_180_AP_inl_frsh = inl_lst[85] 
                    RCP_45_185_AP_inl_frsh = inl_lst[86]     
                    RCP_45_190_AP_inl_frsh = inl_lst[87]              
                    RCP_45_195_AP_inl_frsh = inl_lst[88]
                    RCP_45_200_AP_inl_frsh = inl_lst[89]
                    
                    RCP_85_005_AP_inl_frsh = inl_lst[90]
                    RCP_85_010_AP_inl_frsh = inl_lst[91]
                    RCP_85_015_AP_inl_frsh = inl_lst[92]       
                    RCP_85_020_AP_inl_frsh = inl_lst[93]  
                    RCP_85_025_AP_inl_frsh = inl_lst[94]                
                    RCP_85_030_AP_inl_frsh = inl_lst[95]
                    RCP_85_035_AP_inl_frsh = inl_lst[96]
                    RCP_85_040_AP_inl_frsh = inl_lst[97]     
                    RCP_85_045_AP_inl_frsh = inl_lst[98]
                    RCP_85_050_AP_inl_frsh = inl_lst[99]                 
                    RCP_85_055_AP_inl_frsh = inl_lst[100]
                    RCP_85_060_AP_inl_frsh = inl_lst[101]    
                    RCP_85_065_AP_inl_frsh = inl_lst[102]  
                    RCP_85_070_AP_inl_frsh = inl_lst[103]
                    RCP_85_075_AP_inl_frsh = inl_lst[104]                 
                    RCP_85_080_AP_inl_frsh = inl_lst[105]   
                    RCP_85_085_AP_inl_frsh = inl_lst[106]   
                    RCP_85_090_AP_inl_frsh = inl_lst[107]            
                    RCP_85_095_AP_inl_frsh = inl_lst[108]
                    RCP_85_100_AP_inl_frsh = inl_lst[109]                          
                    RCP_85_105_AP_inl_frsh = inl_lst[110]
                    RCP_85_110_AP_inl_frsh = inl_lst[111]
                    RCP_85_115_AP_inl_frsh = inl_lst[112]              
                    RCP_85_120_AP_inl_frsh = inl_lst[113]
                    RCP_85_125_AP_inl_frsh = inl_lst[114]          
                    RCP_85_130_AP_inl_frsh = inl_lst[115]  
                    RCP_85_135_AP_inl_frsh = inl_lst[116]
                    RCP_85_140_AP_inl_frsh = inl_lst[117]           
                    RCP_85_145_AP_inl_frsh = inl_lst[118]  
                    RCP_85_150_AP_inl_frsh = inl_lst[119]                  
                    RCP_85_155_AP_inl_frsh = inl_lst[120]
                    RCP_85_160_AP_inl_frsh = inl_lst[121]
                    RCP_85_165_AP_inl_frsh = inl_lst[122]      
                    RCP_85_170_AP_inl_frsh = inl_lst[123]
                    RCP_85_175_AP_inl_frsh = inl_lst[124]               
                    RCP_85_180_AP_inl_frsh = inl_lst[125] 
                    RCP_85_185_AP_inl_frsh = inl_lst[126]     
                    RCP_85_190_AP_inl_frsh = inl_lst[127]              
                    RCP_85_195_AP_inl_frsh = inl_lst[128]
                    RCP_85_200_AP_inl_frsh = inl_lst[129]
                    
                    sc_lst.append([coscat_id_str, simu_name.split('SRM_')[-1], simu_nr_str, eq_0_t_end, eq_5000_BP_inl_frsh,\
                                   eq_4500_BP_inl_frsh, eq_4000_BP_inl_frsh, eq_3500_BP_inl_frsh, eq_3000_BP_inl_frsh,\
                                   eq_2500_BP_inl_frsh, eq_2000_BP_inl_frsh, eq_1500_BP_inl_frsh, eq_1000_BP_inl_frsh,\
                                   eq_0500_BP_inl_frsh, eq_0000_BP_inl_frsh,\
                                   RCP_26_005_AP_inl_frsh, RCP_26_010_AP_inl_frsh, RCP_26_015_AP_inl_frsh, RCP_26_020_AP_inl_frsh,\
                                   RCP_26_025_AP_inl_frsh, RCP_26_030_AP_inl_frsh, RCP_26_035_AP_inl_frsh, RCP_26_040_AP_inl_frsh,\
                                   RCP_26_045_AP_inl_frsh, RCP_26_050_AP_inl_frsh, RCP_26_055_AP_inl_frsh, RCP_26_060_AP_inl_frsh,\
                                   RCP_26_065_AP_inl_frsh, RCP_26_070_AP_inl_frsh, RCP_26_075_AP_inl_frsh, RCP_26_080_AP_inl_frsh,\
                                   RCP_26_085_AP_inl_frsh, RCP_26_090_AP_inl_frsh, RCP_26_095_AP_inl_frsh, RCP_26_100_AP_inl_frsh,\
                                   RCP_26_105_AP_inl_frsh, RCP_26_110_AP_inl_frsh, RCP_26_115_AP_inl_frsh, RCP_26_120_AP_inl_frsh,\
                                   RCP_26_125_AP_inl_frsh, RCP_26_130_AP_inl_frsh, RCP_26_135_AP_inl_frsh, RCP_26_140_AP_inl_frsh,\
                                   RCP_26_145_AP_inl_frsh, RCP_26_150_AP_inl_frsh, RCP_26_155_AP_inl_frsh, RCP_26_160_AP_inl_frsh,\
                                   RCP_26_165_AP_inl_frsh, RCP_26_170_AP_inl_frsh, RCP_26_175_AP_inl_frsh, RCP_26_180_AP_inl_frsh,\
                                   RCP_26_185_AP_inl_frsh, RCP_26_190_AP_inl_frsh, RCP_26_195_AP_inl_frsh, RCP_26_200_AP_inl_frsh,\
                                   RCP_45_005_AP_inl_frsh, RCP_45_010_AP_inl_frsh, RCP_45_015_AP_inl_frsh, RCP_45_020_AP_inl_frsh,\
                                   RCP_45_025_AP_inl_frsh, RCP_45_030_AP_inl_frsh, RCP_45_035_AP_inl_frsh, RCP_45_040_AP_inl_frsh,\
                                   RCP_45_045_AP_inl_frsh, RCP_45_050_AP_inl_frsh, RCP_45_055_AP_inl_frsh, RCP_45_060_AP_inl_frsh,\
                                   RCP_45_065_AP_inl_frsh, RCP_45_070_AP_inl_frsh, RCP_45_075_AP_inl_frsh, RCP_45_080_AP_inl_frsh,\
                                   RCP_45_085_AP_inl_frsh, RCP_45_090_AP_inl_frsh, RCP_45_095_AP_inl_frsh, RCP_45_100_AP_inl_frsh,\
                                   RCP_45_105_AP_inl_frsh, RCP_45_110_AP_inl_frsh, RCP_45_115_AP_inl_frsh, RCP_45_120_AP_inl_frsh,\
                                   RCP_45_125_AP_inl_frsh, RCP_45_130_AP_inl_frsh, RCP_45_135_AP_inl_frsh, RCP_45_140_AP_inl_frsh,\
                                   RCP_45_145_AP_inl_frsh, RCP_45_150_AP_inl_frsh, RCP_45_155_AP_inl_frsh, RCP_45_160_AP_inl_frsh,\
                                   RCP_45_165_AP_inl_frsh, RCP_45_170_AP_inl_frsh, RCP_45_175_AP_inl_frsh, RCP_45_180_AP_inl_frsh,\
                                   RCP_45_185_AP_inl_frsh, RCP_45_190_AP_inl_frsh, RCP_45_195_AP_inl_frsh, RCP_45_200_AP_inl_frsh,\
                                   RCP_85_005_AP_inl_frsh, RCP_85_010_AP_inl_frsh, RCP_85_015_AP_inl_frsh, RCP_85_020_AP_inl_frsh,\
                                   RCP_85_025_AP_inl_frsh, RCP_85_030_AP_inl_frsh, RCP_85_035_AP_inl_frsh, RCP_85_040_AP_inl_frsh,\
                                   RCP_85_045_AP_inl_frsh, RCP_85_050_AP_inl_frsh, RCP_85_055_AP_inl_frsh, RCP_85_060_AP_inl_frsh,\
                                   RCP_85_065_AP_inl_frsh, RCP_85_070_AP_inl_frsh, RCP_85_075_AP_inl_frsh, RCP_85_080_AP_inl_frsh,\
                                   RCP_85_085_AP_inl_frsh, RCP_85_090_AP_inl_frsh, RCP_85_095_AP_inl_frsh, RCP_85_100_AP_inl_frsh,\
                                   RCP_85_105_AP_inl_frsh, RCP_85_110_AP_inl_frsh, RCP_85_115_AP_inl_frsh, RCP_85_120_AP_inl_frsh,\
                                   RCP_85_125_AP_inl_frsh, RCP_85_130_AP_inl_frsh, RCP_85_135_AP_inl_frsh, RCP_85_140_AP_inl_frsh,\
                                   RCP_85_145_AP_inl_frsh, RCP_85_150_AP_inl_frsh, RCP_85_155_AP_inl_frsh, RCP_85_160_AP_inl_frsh,\
                                   RCP_85_165_AP_inl_frsh, RCP_85_170_AP_inl_frsh, RCP_85_175_AP_inl_frsh, RCP_85_180_AP_inl_frsh,\
                                   RCP_85_185_AP_inl_frsh, RCP_85_190_AP_inl_frsh, RCP_85_195_AP_inl_frsh, RCP_85_200_AP_inl_frsh])
            
                    all_sc_lst.append(sc_lst)
                except FileNotFoundError:
                    #print('Model simulation did not terminate')
                    continue

    all_sc_lst = [i[0] for i in all_sc_lst]
    df_sum = pd.DataFrame(all_sc_lst, columns = sum_headers_inl)
    df_sum.to_csv(tot_summary_inl_csv)            

    sum_headers_shlf = ['coscat_id', 'srm_id', 'sc_id', 'eq_0_t_end', 'eq_5000_BP_shlf_frsh', 'eq_4500_BP_shlf_frsh',\
                       'eq_4000_BP_shlf_frsh', 'eq_3500_BP_shlf_frsh','eq_3000_BP_shlf_frsh', 'eq_2500_BP_shlf_frsh',\
                       'eq_2000_BP_shlf_frsh', 'eq_1500_BP_shlf_frsh','eq_1000_BP_shlf_frsh', 'eq_0500_BP_shlf_frsh',\
                       'eq_0000_BP_shlf_frsh', 'RCP_26_0005_AP_shlf_frsh', 'RCP_26_0010_AP_shlf_frsh', 'RCP_26_0015_AP_shlf_frsh',\
                       'RCP_26_0020_AP_shlf_frsh', 'RCP_26_0025_AP_shlf_frsh', 'RCP_26_0030_AP_shlf_frsh', 'RCP_26_0035_AP_shlf_frsh',\
                       'RCP_26_0040_AP_shlf_frsh', 'RCP_26_0045_AP_shlf_frsh', 'RCP_26_0050_AP_shlf_frsh', 'RCP_26_0055_AP_shlf_frsh',\
                       'RCP_26_0060_AP_shlf_frsh', 'RCP_26_0065_AP_shlf_frsh', 'RCP_26_0070_AP_shlf_frsh', 'RCP_26_0075_AP_shlf_frsh',\
                       'RCP_26_0080_AP_shlf_frsh', 'RCP_26_0085_AP_shlf_frsh', 'RCP_26_0090_AP_shlf_frsh', 'RCP_26_0095_AP_shlf_frsh',\
                       'RCP_26_0100_AP_shlf_frsh', 'RCP_26_0105_AP_shlf_frsh', 'RCP_26_0110_AP_shlf_frsh', 'RCP_26_0115_AP_shlf_frsh',\
                       'RCP_26_0120_AP_shlf_frsh', 'RCP_26_0125_AP_shlf_frsh', 'RCP_26_0130_AP_shlf_frsh', 'RCP_26_0135_AP_shlf_frsh',\
                       'RCP_26_0140_AP_shlf_frsh', 'RCP_26_0145_AP_shlf_frsh', 'RCP_26_0150_AP_shlf_frsh', 'RCP_26_0155_AP_shlf_frsh',\
                       'RCP_26_0160_AP_shlf_frsh', 'RCP_26_0165_AP_shlf_frsh', 'RCP_26_0170_AP_shlf_frsh', 'RCP_26_0175_AP_shlf_frsh',\
                       'RCP_26_0180_AP_shlf_frsh', 'RCP_26_0185_AP_shlf_frsh', 'RCP_26_0190_AP_shlf_frsh', 'RCP_26_0195_AP_shlf_frsh', 'RCP_26_0200_AP_shlf_frsh',\
                       'RCP_45_0005_AP_shlf_frsh', 'RCP_45_0010_AP_shlf_frsh', 'RCP_45_0015_AP_shlf_frsh',\
                       'RCP_45_0020_AP_shlf_frsh', 'RCP_45_0025_AP_shlf_frsh', 'RCP_45_0030_AP_shlf_frsh', 'RCP_45_0035_AP_shlf_frsh',\
                       'RCP_45_0040_AP_shlf_frsh', 'RCP_45_0045_AP_shlf_frsh', 'RCP_45_0050_AP_shlf_frsh', 'RCP_45_0055_AP_shlf_frsh',\
                       'RCP_45_0060_AP_shlf_frsh', 'RCP_45_0065_AP_shlf_frsh', 'RCP_45_0070_AP_shlf_frsh', 'RCP_45_0075_AP_shlf_frsh',\
                       'RCP_45_0080_AP_shlf_frsh', 'RCP_45_0085_AP_shlf_frsh', 'RCP_45_0090_AP_shlf_frsh', 'RCP_45_0095_AP_shlf_frsh',\
                       'RCP_45_0100_AP_shlf_frsh', 'RCP_45_0105_AP_shlf_frsh', 'RCP_45_0110_AP_shlf_frsh', 'RCP_45_0115_AP_shlf_frsh',\
                       'RCP_45_0120_AP_shlf_frsh', 'RCP_45_0125_AP_shlf_frsh', 'RCP_45_0130_AP_shlf_frsh', 'RCP_45_0135_AP_shlf_frsh',\
                       'RCP_45_0140_AP_shlf_frsh', 'RCP_45_0145_AP_shlf_frsh', 'RCP_45_0150_AP_shlf_frsh', 'RCP_45_0155_AP_shlf_frsh',\
                       'RCP_45_0160_AP_shlf_frsh', 'RCP_45_0165_AP_shlf_frsh', 'RCP_45_0170_AP_shlf_frsh', 'RCP_45_0175_AP_shlf_frsh',\
                       'RCP_45_0180_AP_shlf_frsh', 'RCP_45_0185_AP_shlf_frsh', 'RCP_45_0190_AP_shlf_frsh', 'RCP_45_0195_AP_shlf_frsh', 'RCP_45_0200_AP_shlf_frsh',\
                       'RCP_85_0005_AP_shlf_frsh', 'RCP_85_0010_AP_shlf_frsh', 'RCP_85_0015_AP_shlf_frsh',\
                       'RCP_85_0020_AP_shlf_frsh', 'RCP_85_0025_AP_shlf_frsh', 'RCP_85_0030_AP_shlf_frsh', 'RCP_85_0035_AP_shlf_frsh',\
                       'RCP_85_0040_AP_shlf_frsh', 'RCP_85_0045_AP_shlf_frsh', 'RCP_85_0050_AP_shlf_frsh', 'RCP_85_0055_AP_shlf_frsh',\
                       'RCP_85_0060_AP_shlf_frsh', 'RCP_85_0065_AP_shlf_frsh', 'RCP_85_0070_AP_shlf_frsh', 'RCP_85_0075_AP_shlf_frsh',\
                       'RCP_85_0080_AP_shlf_frsh', 'RCP_85_0085_AP_shlf_frsh', 'RCP_85_0090_AP_shlf_frsh', 'RCP_85_0095_AP_shlf_frsh',\
                       'RCP_85_0100_AP_shlf_frsh', 'RCP_85_0105_AP_shlf_frsh', 'RCP_85_0110_AP_shlf_frsh', 'RCP_85_0115_AP_shlf_frsh',\
                       'RCP_85_0120_AP_shlf_frsh', 'RCP_85_0125_AP_shlf_frsh', 'RCP_85_0130_AP_shlf_frsh', 'RCP_85_0135_AP_shlf_frsh',\
                       'RCP_85_0140_AP_shlf_frsh', 'RCP_85_0145_AP_shlf_frsh', 'RCP_85_0150_AP_shlf_frsh', 'RCP_85_0155_AP_shlf_frsh',\
                       'RCP_85_0160_AP_shlf_frsh', 'RCP_85_0165_AP_shlf_frsh', 'RCP_85_0170_AP_shlf_frsh', 'RCP_85_0175_AP_shlf_frsh',\
                       'RCP_85_0180_AP_shlf_frsh', 'RCP_85_0185_AP_shlf_frsh', 'RCP_85_0190_AP_shlf_frsh', 'RCP_85_0195_AP_shlf_frsh', 'RCP_85_0200_AP_shlf_frsh']
    all_sc_lst = []
    
    #   loop again and write in the values into the right columns this time
    for dirpath, dirnames, filenames in os.walk(csv_dir):
        for filename in [f for f in filenames if f.endswith('.csv')]:
            simu_nr_str = filename.split('sc_')[-1].split('_')[0]
            #print(simu_nr_str)

            df = pd.read_csv(os.path.join(csv_dir, filename), index_col = False, header = 0)    
            if len(df.loc[df['EQ'] == '3000BC_to_2000BC'][' time_end'].values) == 2:
                
                shlf_lst = df[' frsh_shelf_pct'].values.tolist()                    
                if len(shlf_lst) < max_len:
                    for f in range(max_len - len(shlf_lst)):
                        shlf_lst.append(np.nan)                

                try:                                       
                    sc_lst = []

                    eq_0_t_end = shlf_lst[0]
                    eq_5000_BP_shlf_frsh = shlf_lst[0]
                    eq_4500_BP_shlf_frsh = shlf_lst[0]           
                    eq_4000_BP_shlf_frsh = shlf_lst[1]
                    eq_3500_BP_shlf_frsh = shlf_lst[2]
                    eq_3000_BP_shlf_frsh = shlf_lst[3]   
                    eq_2500_BP_shlf_frsh = shlf_lst[4]              
                    eq_2000_BP_shlf_frsh = shlf_lst[5]
                    eq_1500_BP_shlf_frsh = shlf_lst[6]
                    eq_1000_BP_shlf_frsh = shlf_lst[7] 
                    eq_0500_BP_shlf_frsh = shlf_lst[8]
                    eq_0000_BP_shlf_frsh = shlf_lst[9]    
                                        
                    RCP_26_005_AP_shlf_frsh = shlf_lst[10]
                    RCP_26_010_AP_shlf_frsh = shlf_lst[11]
                    RCP_26_015_AP_shlf_frsh = shlf_lst[12]       
                    RCP_26_020_AP_shlf_frsh = shlf_lst[13]  
                    RCP_26_025_AP_shlf_frsh = shlf_lst[14]                
                    RCP_26_030_AP_shlf_frsh = shlf_lst[15]
                    RCP_26_035_AP_shlf_frsh = shlf_lst[16]
                    RCP_26_040_AP_shlf_frsh = shlf_lst[17]     
                    RCP_26_045_AP_shlf_frsh = shlf_lst[18]
                    RCP_26_050_AP_shlf_frsh = shlf_lst[19]                 
                    RCP_26_055_AP_shlf_frsh = shlf_lst[20]
                    RCP_26_060_AP_shlf_frsh = shlf_lst[21]    
                    RCP_26_065_AP_shlf_frsh = shlf_lst[22]  
                    RCP_26_070_AP_shlf_frsh = shlf_lst[23]
                    RCP_26_075_AP_shlf_frsh = shlf_lst[24]                 
                    RCP_26_080_AP_shlf_frsh = shlf_lst[25]   
                    RCP_26_085_AP_shlf_frsh = shlf_lst[26]   
                    RCP_26_090_AP_shlf_frsh = shlf_lst[27]            
                    RCP_26_095_AP_shlf_frsh = shlf_lst[28]
                    RCP_26_100_AP_shlf_frsh = shlf_lst[29]                          
                    RCP_26_105_AP_shlf_frsh = shlf_lst[30]
                    RCP_26_110_AP_shlf_frsh = shlf_lst[31]
                    RCP_26_115_AP_shlf_frsh = shlf_lst[32]              
                    RCP_26_120_AP_shlf_frsh = shlf_lst[33]
                    RCP_26_125_AP_shlf_frsh = shlf_lst[34]          
                    RCP_26_130_AP_shlf_frsh = shlf_lst[35]  
                    RCP_26_135_AP_shlf_frsh = shlf_lst[36]
                    RCP_26_140_AP_shlf_frsh = shlf_lst[37]           
                    RCP_26_145_AP_shlf_frsh = shlf_lst[38]  
                    RCP_26_150_AP_shlf_frsh = shlf_lst[39]                  
                    RCP_26_155_AP_shlf_frsh = shlf_lst[40]
                    RCP_26_160_AP_shlf_frsh = shlf_lst[41]
                    RCP_26_165_AP_shlf_frsh = shlf_lst[42]      
                    RCP_26_170_AP_shlf_frsh = shlf_lst[43]
                    RCP_26_175_AP_shlf_frsh = shlf_lst[44]               
                    RCP_26_180_AP_shlf_frsh = shlf_lst[45] 
                    RCP_26_185_AP_shlf_frsh = shlf_lst[46]     
                    RCP_26_190_AP_shlf_frsh = shlf_lst[47]              
                    RCP_26_195_AP_shlf_frsh = shlf_lst[48]
                    RCP_26_200_AP_shlf_frsh = shlf_lst[49]

                    RCP_45_005_AP_shlf_frsh = shlf_lst[50]
                    RCP_45_010_AP_shlf_frsh = shlf_lst[51]
                    RCP_45_015_AP_shlf_frsh = shlf_lst[52]       
                    RCP_45_020_AP_shlf_frsh = shlf_lst[53]  
                    RCP_45_025_AP_shlf_frsh = shlf_lst[54]                
                    RCP_45_030_AP_shlf_frsh = shlf_lst[55]
                    RCP_45_035_AP_shlf_frsh = shlf_lst[56]
                    RCP_45_040_AP_shlf_frsh = shlf_lst[57]     
                    RCP_45_045_AP_shlf_frsh = shlf_lst[58]
                    RCP_45_050_AP_shlf_frsh = shlf_lst[59]                 
                    RCP_45_055_AP_shlf_frsh = shlf_lst[60]
                    RCP_45_060_AP_shlf_frsh = shlf_lst[61]    
                    RCP_45_065_AP_shlf_frsh = shlf_lst[62]  
                    RCP_45_070_AP_shlf_frsh = shlf_lst[63]
                    RCP_45_075_AP_shlf_frsh = shlf_lst[64]                 
                    RCP_45_080_AP_shlf_frsh = shlf_lst[65]   
                    RCP_45_085_AP_shlf_frsh = shlf_lst[66]   
                    RCP_45_090_AP_shlf_frsh = shlf_lst[67]            
                    RCP_45_095_AP_shlf_frsh = shlf_lst[68]
                    RCP_45_100_AP_shlf_frsh = shlf_lst[69]                          
                    RCP_45_105_AP_shlf_frsh = shlf_lst[70]
                    RCP_45_110_AP_shlf_frsh = shlf_lst[71]
                    RCP_45_115_AP_shlf_frsh = shlf_lst[72]              
                    RCP_45_120_AP_shlf_frsh = shlf_lst[73]
                    RCP_45_125_AP_shlf_frsh = shlf_lst[74]          
                    RCP_45_130_AP_shlf_frsh = shlf_lst[75]  
                    RCP_45_135_AP_shlf_frsh = shlf_lst[76]
                    RCP_45_140_AP_shlf_frsh = shlf_lst[77]           
                    RCP_45_145_AP_shlf_frsh = shlf_lst[78]  
                    RCP_45_150_AP_shlf_frsh = shlf_lst[79]                  
                    RCP_45_155_AP_shlf_frsh = shlf_lst[80]
                    RCP_45_160_AP_shlf_frsh = shlf_lst[81]
                    RCP_45_165_AP_shlf_frsh = shlf_lst[82]      
                    RCP_45_170_AP_shlf_frsh = shlf_lst[83]
                    RCP_45_175_AP_shlf_frsh = shlf_lst[84]               
                    RCP_45_180_AP_shlf_frsh = shlf_lst[85] 
                    RCP_45_185_AP_shlf_frsh = shlf_lst[86]     
                    RCP_45_190_AP_shlf_frsh = shlf_lst[87]              
                    RCP_45_195_AP_shlf_frsh = shlf_lst[88]
                    RCP_45_200_AP_shlf_frsh = shlf_lst[89]
                    
                    RCP_85_005_AP_shlf_frsh = shlf_lst[90]
                    RCP_85_010_AP_shlf_frsh = shlf_lst[91]
                    RCP_85_015_AP_shlf_frsh = shlf_lst[92]       
                    RCP_85_020_AP_shlf_frsh = shlf_lst[93]  
                    RCP_85_025_AP_shlf_frsh = shlf_lst[94]                
                    RCP_85_030_AP_shlf_frsh = shlf_lst[95]
                    RCP_85_035_AP_shlf_frsh = shlf_lst[96]
                    RCP_85_040_AP_shlf_frsh = shlf_lst[97]     
                    RCP_85_045_AP_shlf_frsh = shlf_lst[98]
                    RCP_85_050_AP_shlf_frsh = shlf_lst[99]                 
                    RCP_85_055_AP_shlf_frsh = shlf_lst[100]
                    RCP_85_060_AP_shlf_frsh = shlf_lst[101]    
                    RCP_85_065_AP_shlf_frsh = shlf_lst[102]  
                    RCP_85_070_AP_shlf_frsh = shlf_lst[103]
                    RCP_85_075_AP_shlf_frsh = shlf_lst[104]                 
                    RCP_85_080_AP_shlf_frsh = shlf_lst[105]   
                    RCP_85_085_AP_shlf_frsh = shlf_lst[106]   
                    RCP_85_090_AP_shlf_frsh = shlf_lst[107]            
                    RCP_85_095_AP_shlf_frsh = shlf_lst[108]
                    RCP_85_100_AP_shlf_frsh = shlf_lst[109]                          
                    RCP_85_105_AP_shlf_frsh = shlf_lst[110]
                    RCP_85_110_AP_shlf_frsh = shlf_lst[111]
                    RCP_85_115_AP_shlf_frsh = shlf_lst[112]              
                    RCP_85_120_AP_shlf_frsh = shlf_lst[113]
                    RCP_85_125_AP_shlf_frsh = shlf_lst[114]          
                    RCP_85_130_AP_shlf_frsh = shlf_lst[115]  
                    RCP_85_135_AP_shlf_frsh = shlf_lst[116]
                    RCP_85_140_AP_shlf_frsh = shlf_lst[117]           
                    RCP_85_145_AP_shlf_frsh = shlf_lst[118]  
                    RCP_85_150_AP_shlf_frsh = shlf_lst[119]                  
                    RCP_85_155_AP_shlf_frsh = shlf_lst[120]
                    RCP_85_160_AP_shlf_frsh = shlf_lst[121]
                    RCP_85_165_AP_shlf_frsh = shlf_lst[122]      
                    RCP_85_170_AP_shlf_frsh = shlf_lst[123]
                    RCP_85_175_AP_shlf_frsh = shlf_lst[124]               
                    RCP_85_180_AP_shlf_frsh = shlf_lst[125] 
                    RCP_85_185_AP_shlf_frsh = shlf_lst[126]     
                    RCP_85_190_AP_shlf_frsh = shlf_lst[127]              
                    RCP_85_195_AP_shlf_frsh = shlf_lst[128]
                    RCP_85_200_AP_shlf_frsh = shlf_lst[129]
                    
                    """
                    eq_0_t_end = df[(df['EQ'] == eq_lst[0])][' time_end'].values[-1]
                    eq_5000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[0])][' frsh_shelf_pct'].values[0]
                    eq_4500_BP_shlf_frsh = df[(df['EQ'] == eq_lst[0])][' frsh_shelf_pct'].values[0]                
                    eq_4000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[0])][' frsh_shelf_pct'].values[-1]
                    eq_3500_BP_shlf_frsh = df[(df['EQ'] == eq_lst[1])][' frsh_shelf_pct'].values[0]     
                    eq_3000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[1])][' frsh_shelf_pct'].values[-1]     
                    eq_2500_BP_shlf_frsh = df[(df['EQ'] == eq_lst[2])][' frsh_shelf_pct'].values[0]                 
                    eq_2000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[2])][' frsh_shelf_pct'].values[-1]
                    eq_1500_BP_shlf_frsh = df[(df['EQ'] == eq_lst[3])][' frsh_shelf_pct'].values[0]     
                    eq_1000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[3])][' frsh_shelf_pct'].values[-1]     
                    eq_0500_BP_shlf_frsh = df[(df['EQ'] == eq_lst[4])][' frsh_shelf_pct'].values[0] 
                    eq_0000_BP_shlf_frsh = df[(df['EQ'] == eq_lst[4])][' frsh_shelf_pct'].values[-1]     
                    
                    RCP_26_005_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[0]     
                    RCP_26_010_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_015_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[0]                
                    RCP_26_020_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[-1]   
                    RCP_26_025_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[0]                     
                    RCP_26_030_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_035_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[0]     
                    RCP_26_040_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[-1]                
                    RCP_26_045_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[0]   
                    RCP_26_050_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[-1]                     
                    RCP_26_055_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[5])][' frsh_shelf_pct'].values[0]     
                    RCP_26_060_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[5])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_065_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[6])][' frsh_shelf_pct'].values[0]                
                    RCP_26_070_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[6])][' frsh_shelf_pct'].values[-1]   
                    RCP_26_075_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[7])][' frsh_shelf_pct'].values[0]                     
                    RCP_26_080_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[7])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_085_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[8])][' frsh_shelf_pct'].values[0]     
                    RCP_26_090_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[8])][' frsh_shelf_pct'].values[-1]                
                    RCP_26_095_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[9])][' frsh_shelf_pct'].values[0]   
                    RCP_26_100_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[9])][' frsh_shelf_pct'].values[-1]                              
                    RCP_26_105_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[10])][' frsh_shelf_pct'].values[0]     
                    RCP_26_110_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[10])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_115_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[11])][' frsh_shelf_pct'].values[0]                
                    RCP_26_120_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[11])][' frsh_shelf_pct'].values[-1]   
                    RCP_26_125_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[12])][' frsh_shelf_pct'].values[0]                     
                    RCP_26_130_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[12])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_135_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[13])][' frsh_shelf_pct'].values[0]     
                    RCP_26_140_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[13])][' frsh_shelf_pct'].values[-1]                
                    RCP_26_145_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[14])][' frsh_shelf_pct'].values[0]   
                    RCP_26_150_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[14])][' frsh_shelf_pct'].values[-1]                     
                    RCP_26_155_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[15])][' frsh_shelf_pct'].values[0]     
                    RCP_26_160_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[15])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_165_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[16])][' frsh_shelf_pct'].values[0]                
                    RCP_26_170_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[16])][' frsh_shelf_pct'].values[-1]   
                    RCP_26_175_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[17])][' frsh_shelf_pct'].values[0]                     
                    RCP_26_180_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[17])][' frsh_shelf_pct'].values[-1]     
                    RCP_26_185_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[18])][' frsh_shelf_pct'].values[0]     
                    RCP_26_190_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[18])][' frsh_shelf_pct'].values[-1]                
                    RCP_26_195_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[19])][' frsh_shelf_pct'].values[0]   
                    RCP_26_200_AP_shlf_frsh = df[(df['EQ'] == 'RCP_26_' + eq_slr_lst[19])][' frsh_shelf_pct'].values[-1]    

                    RCP_45_005_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[0]     
                    RCP_45_010_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_015_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[0]                
                    RCP_45_020_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[-1]   
                    RCP_45_025_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[0]                     
                    RCP_45_030_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_035_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[0]     
                    RCP_45_040_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[-1]                
                    RCP_45_045_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[0]   
                    RCP_45_050_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[-1]                     
                    RCP_45_055_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[5])][' frsh_shelf_pct'].values[0]     
                    RCP_45_060_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[5])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_065_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[6])][' frsh_shelf_pct'].values[0]                
                    RCP_45_070_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[6])][' frsh_shelf_pct'].values[-1]   
                    RCP_45_075_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[7])][' frsh_shelf_pct'].values[0]                     
                    RCP_45_080_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[7])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_085_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[8])][' frsh_shelf_pct'].values[0]     
                    RCP_45_090_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[8])][' frsh_shelf_pct'].values[-1]                
                    RCP_45_095_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[9])][' frsh_shelf_pct'].values[0]   
                    RCP_45_100_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[9])][' frsh_shelf_pct'].values[-1]                              
                    RCP_45_105_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[10])][' frsh_shelf_pct'].values[0]     
                    RCP_45_110_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[10])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_115_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[11])][' frsh_shelf_pct'].values[0]                
                    RCP_45_120_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[11])][' frsh_shelf_pct'].values[-1]   
                    RCP_45_125_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[12])][' frsh_shelf_pct'].values[0]                     
                    RCP_45_130_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[12])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_135_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[13])][' frsh_shelf_pct'].values[0]     
                    RCP_45_140_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[13])][' frsh_shelf_pct'].values[-1]                
                    RCP_45_145_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[14])][' frsh_shelf_pct'].values[0]   
                    RCP_45_150_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[14])][' frsh_shelf_pct'].values[-1]                     
                    RCP_45_155_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[15])][' frsh_shelf_pct'].values[0]     
                    RCP_45_160_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[15])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_165_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[16])][' frsh_shelf_pct'].values[0]                
                    RCP_45_170_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[16])][' frsh_shelf_pct'].values[-1]   
                    RCP_45_175_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[17])][' frsh_shelf_pct'].values[0]                     
                    RCP_45_180_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[17])][' frsh_shelf_pct'].values[-1]     
                    RCP_45_185_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[18])][' frsh_shelf_pct'].values[0]     
                    RCP_45_190_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[18])][' frsh_shelf_pct'].values[-1]                
                    RCP_45_195_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[19])][' frsh_shelf_pct'].values[0]   
                    RCP_45_200_AP_shlf_frsh = df[(df['EQ'] == 'RCP_45_' + eq_slr_lst[19])][' frsh_shelf_pct'].values[-1] 
                    
                    RCP_85_005_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[0]     
                    RCP_85_010_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[0])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_015_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[0]                
                    RCP_85_020_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[1])][' frsh_shelf_pct'].values[-1]   
                    RCP_85_025_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[0]                     
                    RCP_85_030_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[2])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_035_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[0]     
                    RCP_85_040_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[3])][' frsh_shelf_pct'].values[-1]                
                    RCP_85_045_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[0]   
                    RCP_85_050_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[4])][' frsh_shelf_pct'].values[-1]                     
                    RCP_85_055_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[5])][' frsh_shelf_pct'].values[0]     
                    RCP_85_060_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[5])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_065_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[6])][' frsh_shelf_pct'].values[0]                
                    RCP_85_070_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[6])][' frsh_shelf_pct'].values[-1]   
                    RCP_85_075_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[7])][' frsh_shelf_pct'].values[0]                     
                    RCP_85_080_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[7])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_085_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[8])][' frsh_shelf_pct'].values[0]     
                    RCP_85_090_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[8])][' frsh_shelf_pct'].values[-1]                
                    RCP_85_095_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[9])][' frsh_shelf_pct'].values[0]   
                    RCP_85_100_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[9])][' frsh_shelf_pct'].values[-1]                              
                    RCP_85_105_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[10])][' frsh_shelf_pct'].values[0]     
                    RCP_85_110_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[10])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_115_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[11])][' frsh_shelf_pct'].values[0]                
                    RCP_85_120_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[11])][' frsh_shelf_pct'].values[-1]   
                    RCP_85_125_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[12])][' frsh_shelf_pct'].values[0]                     
                    RCP_85_130_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[12])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_135_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[13])][' frsh_shelf_pct'].values[0]     
                    RCP_85_140_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[13])][' frsh_shelf_pct'].values[-1]                
                    RCP_85_145_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[14])][' frsh_shelf_pct'].values[0]   
                    RCP_85_150_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[14])][' frsh_shelf_pct'].values[-1]                     
                    RCP_85_155_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[15])][' frsh_shelf_pct'].values[0]     
                    RCP_85_160_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[15])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_165_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[16])][' frsh_shelf_pct'].values[0]                
                    RCP_85_170_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[16])][' frsh_shelf_pct'].values[-1]   
                    RCP_85_175_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[17])][' frsh_shelf_pct'].values[0]                     
                    RCP_85_180_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[17])][' frsh_shelf_pct'].values[-1]     
                    RCP_85_185_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[18])][' frsh_shelf_pct'].values[0]     
                    RCP_85_190_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[18])][' frsh_shelf_pct'].values[-1]                
                    RCP_85_195_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[19])][' frsh_shelf_pct'].values[0]   
                    RCP_85_200_AP_shlf_frsh = df[(df['EQ'] == 'RCP_85_' + eq_slr_lst[19])][' frsh_shelf_pct'].values[-1] 
                    """
                    
                    sc_lst.append([coscat_id_str, simu_name.split('SRM_')[-1], simu_nr_str, eq_0_t_end, eq_5000_BP_shlf_frsh,\
                                   eq_4500_BP_shlf_frsh, eq_4000_BP_shlf_frsh, eq_3500_BP_shlf_frsh, eq_3000_BP_shlf_frsh,\
                                   eq_2500_BP_shlf_frsh, eq_2000_BP_shlf_frsh, eq_1500_BP_shlf_frsh, eq_1000_BP_shlf_frsh,\
                                   eq_0500_BP_shlf_frsh, eq_0000_BP_shlf_frsh,\
                                   RCP_26_005_AP_shlf_frsh, RCP_26_010_AP_shlf_frsh, RCP_26_015_AP_shlf_frsh, RCP_26_020_AP_shlf_frsh,\
                                   RCP_26_025_AP_shlf_frsh, RCP_26_030_AP_shlf_frsh, RCP_26_035_AP_shlf_frsh, RCP_26_040_AP_shlf_frsh,\
                                   RCP_26_045_AP_shlf_frsh, RCP_26_050_AP_shlf_frsh, RCP_26_055_AP_shlf_frsh, RCP_26_060_AP_shlf_frsh,\
                                   RCP_26_065_AP_shlf_frsh, RCP_26_070_AP_shlf_frsh, RCP_26_075_AP_shlf_frsh, RCP_26_080_AP_shlf_frsh,\
                                   RCP_26_085_AP_shlf_frsh, RCP_26_090_AP_shlf_frsh, RCP_26_095_AP_shlf_frsh, RCP_26_100_AP_shlf_frsh,\
                                   RCP_26_105_AP_shlf_frsh, RCP_26_110_AP_shlf_frsh, RCP_26_115_AP_shlf_frsh, RCP_26_120_AP_shlf_frsh,\
                                   RCP_26_125_AP_shlf_frsh, RCP_26_130_AP_shlf_frsh, RCP_26_135_AP_shlf_frsh, RCP_26_140_AP_shlf_frsh,\
                                   RCP_26_145_AP_shlf_frsh, RCP_26_150_AP_shlf_frsh, RCP_26_155_AP_shlf_frsh, RCP_26_160_AP_shlf_frsh,\
                                   RCP_26_165_AP_shlf_frsh, RCP_26_170_AP_shlf_frsh, RCP_26_175_AP_shlf_frsh, RCP_26_180_AP_shlf_frsh,\
                                   RCP_26_185_AP_shlf_frsh, RCP_26_190_AP_shlf_frsh, RCP_26_195_AP_shlf_frsh, RCP_26_200_AP_shlf_frsh,\
                                   RCP_45_005_AP_shlf_frsh, RCP_45_010_AP_shlf_frsh, RCP_45_015_AP_shlf_frsh, RCP_45_020_AP_shlf_frsh,\
                                   RCP_45_025_AP_shlf_frsh, RCP_45_030_AP_shlf_frsh, RCP_45_035_AP_shlf_frsh, RCP_45_040_AP_shlf_frsh,\
                                   RCP_45_045_AP_shlf_frsh, RCP_45_050_AP_shlf_frsh, RCP_45_055_AP_shlf_frsh, RCP_45_060_AP_shlf_frsh,\
                                   RCP_45_065_AP_shlf_frsh, RCP_45_070_AP_shlf_frsh, RCP_45_075_AP_shlf_frsh, RCP_45_080_AP_shlf_frsh,\
                                   RCP_45_085_AP_shlf_frsh, RCP_45_090_AP_shlf_frsh, RCP_45_095_AP_shlf_frsh, RCP_45_100_AP_shlf_frsh,\
                                   RCP_45_105_AP_shlf_frsh, RCP_45_110_AP_shlf_frsh, RCP_45_115_AP_shlf_frsh, RCP_45_120_AP_shlf_frsh,\
                                   RCP_45_125_AP_shlf_frsh, RCP_45_130_AP_shlf_frsh, RCP_45_135_AP_shlf_frsh, RCP_45_140_AP_shlf_frsh,\
                                   RCP_45_145_AP_shlf_frsh, RCP_45_150_AP_shlf_frsh, RCP_45_155_AP_shlf_frsh, RCP_45_160_AP_shlf_frsh,\
                                   RCP_45_165_AP_shlf_frsh, RCP_45_170_AP_shlf_frsh, RCP_45_175_AP_shlf_frsh, RCP_45_180_AP_shlf_frsh,\
                                   RCP_45_185_AP_shlf_frsh, RCP_45_190_AP_shlf_frsh, RCP_45_195_AP_shlf_frsh, RCP_45_200_AP_shlf_frsh,\
                                   RCP_85_005_AP_shlf_frsh, RCP_85_010_AP_shlf_frsh, RCP_85_015_AP_shlf_frsh, RCP_85_020_AP_shlf_frsh,\
                                   RCP_85_025_AP_shlf_frsh, RCP_85_030_AP_shlf_frsh, RCP_85_035_AP_shlf_frsh, RCP_85_040_AP_shlf_frsh,\
                                   RCP_85_045_AP_shlf_frsh, RCP_85_050_AP_shlf_frsh, RCP_85_055_AP_shlf_frsh, RCP_85_060_AP_shlf_frsh,\
                                   RCP_85_065_AP_shlf_frsh, RCP_85_070_AP_shlf_frsh, RCP_85_075_AP_shlf_frsh, RCP_85_080_AP_shlf_frsh,\
                                   RCP_85_085_AP_shlf_frsh, RCP_85_090_AP_shlf_frsh, RCP_85_095_AP_shlf_frsh, RCP_85_100_AP_shlf_frsh,\
                                   RCP_85_105_AP_shlf_frsh, RCP_85_110_AP_shlf_frsh, RCP_85_115_AP_shlf_frsh, RCP_85_120_AP_shlf_frsh,\
                                   RCP_85_125_AP_shlf_frsh, RCP_85_130_AP_shlf_frsh, RCP_85_135_AP_shlf_frsh, RCP_85_140_AP_shlf_frsh,\
                                   RCP_85_145_AP_shlf_frsh, RCP_85_150_AP_shlf_frsh, RCP_85_155_AP_shlf_frsh, RCP_85_160_AP_shlf_frsh,\
                                   RCP_85_165_AP_shlf_frsh, RCP_85_170_AP_shlf_frsh, RCP_85_175_AP_shlf_frsh, RCP_85_180_AP_shlf_frsh,\
                                   RCP_85_185_AP_shlf_frsh, RCP_85_190_AP_shlf_frsh, RCP_85_195_AP_shlf_frsh, RCP_85_200_AP_shlf_frsh])

                    all_sc_lst.append(sc_lst)
                except FileNotFoundError:
                    #print('Model simulation did not terminate')
                    continue

    all_sc_lst = [i[0] for i in all_sc_lst]
    df_sum = pd.DataFrame(all_sc_lst, columns = sum_headers_shlf)
    df_sum.to_csv(tot_summary_shlf_csv)            

    