# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 16:56:41 2018

@author: daniel
"""


import matplotlib.pyplot as plt
import math
import matplotlib
from matplotlib import rcParams
import numpy as np
import os
import matplotlib as mpl

"""
plot_arr = conc.get_data(totim = t_conc[1])
zoom_in = [-2000., 2000.]
sea_lvl = sea_level
real_time_yrs = 1000.

dest_dir = dest_dir_conc
xmin = x_start
xmax = x_foot_of_slope
ymax = y_start_top
ymin = y_foot_of_slope_bot   
delcol = del_col
dellay = del_lay
hk_array = hk_arr
aqt_hk_val = hk_aqt_val
"""


#   Plotting of 2D concentration profiles at given time
def plot_conc_profile_yrs(dest_dir, plot_arr, top_el, bot_el, sea_lvl, real_time_yrs, xmin, xmax, ymin, ymax, delcol, dellay,
                          hk_array, aqt_hk_val, zoom):

    #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
    fig = plt.figure(figsize = (20, 15))
    ax1 = plt.subplot2grid((2, 2), (0, 1))  #   overall concentration profiles
    ax2 = plt.subplot2grid((2, 2), (1, 1))  #   zoomed in concentration profile
    ax3 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)  #   color bar area 
    ax1.set_position([0.135, 0.55, 0.85, 0.4]) # [left, bottom, width, height]
    ax2.set_position([0.135, 0.1, 0.85, 0.4])
    ax3.set_position([0.04, 0.2, 0.025, 0.6])        
    
    #   adapt the time in years to fit better the folder naming structure
    if real_time_yrs < 10.:
        str_plot_time = '000000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 10. and real_time_yrs < 100.:
        str_plot_time = '00000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 100. and real_time_yrs < 1000.:
        str_plot_time = '0000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 100. and real_time_yrs < 10000.:
        str_plot_time = '000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 1000. and real_time_yrs < 10000.:
        str_plot_time = '00' + str(int(real_time_yrs)) 
    else:
        str_plot_time = '0' + str(int(real_time_yrs))     
    
    if zoom:
        xmin_zoom = zoom[0]
        xmax_zoom = zoom[1]
        idx_start_zoom = int(abs((xmin - xmin_zoom) / delcol))
        idx_end_zoom = int(abs((xmin - xmax_zoom) / delcol))
        ymax_zoom = math.ceil(top_el[idx_start_zoom] / dellay) * dellay
        ymin_zoom = math.floor(bot_el[idx_end_zoom] / dellay) * dellay
        
    #   add the location of the aquitard layers    
    mapimg = (hk_array == aqt_hk_val)
       
    #   read in the concentration profile from the dictionary, find the closest time step to the input tume_step_yrs
    
    vmin = 0.0
    vmax = 35.0
    title = 'Concentration profile at ' + str(real_time_yrs) + ' years since simulation start'
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    # define the bins and normalize
    bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    cbar_title = 'Concentration of \n salt (g/l)'
    plot_arr[np.abs(plot_arr) >= 999.] = np.nan
    #   create output name for the individual figure
    fig_out_dir = os.path.join(dest_dir, '_conc_' + str_plot_time + '_yrs.png')    
    
    im1 = ax1.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                   extent = (xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
    ax1.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                   extent = (xmin, xmax, ymin, ymax), alpha = 0.25)
    
    ax2.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                   extent = (xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
    ax2.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                   extent = (xmin, xmax, ymin, ymax), alpha = 0.25)
    
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
       
    ax2.set_xlim([xmin_zoom, xmax_zoom])
    ax2.set_ylim([ymin_zoom, ymax_zoom])
    
    #   set the gridlines and constant lines in the plot
    #   add constant lines with elevation = 0m asl. and coastline 
    ax1.axhline(y = 0., linewidth = 2, color = 'k', zorder = 2)
    ax1.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)
    ax2.axhline(y = 0., linewidth = 2, color = 'k', zorder = 2)
    ax2.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)    
    
    ax1.axhline(y = sea_lvl, linewidth = 3, color = 'red', zorder = 2)    
    ax2.axhline(y = sea_lvl, linewidth = 3, color = 'red', zorder = 2)    

    x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(xmin / 5.) * 5., math.ceil(xmax), 1000.0), nbins = None)    
    x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(xmin / 1.) * 1., math.ceil(xmax), 250.0), nbins = None)    
    ax1.xaxis.set_major_locator(x_major_locator)
    ax1.xaxis.set_minor_locator(x_minor_locator)
 
    y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(ymin / 500.) * 500., math.ceil(ymax), 250.0), nbins = None)    
    y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(ymin / 100.) * 100. , math.ceil(ymax), 50.0), nbins = None)    
    ax1.yaxis.set_major_locator(y_major_locator)
    ax1.yaxis.set_minor_locator(y_minor_locator)
    
    ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
    ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  

    x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(xmin_zoom / 5.) * 5., math.ceil(xmax_zoom), 1000.0), nbins = None)    
    x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(xmin_zoom / 1.) * 1., math.ceil(xmax_zoom), 250.0), nbins = None)    
    ax2.xaxis.set_major_locator(x_major_locator_zoom)
    ax2.xaxis.set_minor_locator(x_minor_locator_zoom)

    ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
    ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)        
            
    ax1.set_title(title, fontsize = 24, y = 1.025)
    ax1.set_xlabel('distance from coast (m)', fontsize = 18)
    ax1.set_ylabel('elevation (m asl.)', fontsize = 18)
    ax2.set_xlabel('distance from coast (m)', fontsize = 18)
    ax2.set_ylabel('elevation (m asl.)', fontsize = 18)    
    
    #   plot the colorbar
    cbar = plt.colorbar(im1, cax = ax3, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
    cbar.ax.set_title(cbar_title, fontsize = 14, y  = 1.025)
    cbar.ax.tick_params(labelsize = 12)    
    
    #   set the font 
    rcParams['font.family'] = 'Garamond'
    rcParams['axes.facecolor'] = 'white'
    rcParams['savefig.facecolor'] = 'white'    
    
    if os.path.exists(fig_out_dir):
        print('Figure   ' + fig_out_dir + '   already exists.')
    else: 
        #   save the figure
        plt.savefig(fig_out_dir, dpi = 300)
        plt.close()
        del fig




#   Plotting of 2D concentration profiles at given time
def plot_head_profile_yrs(dest_dir, plot_arr, top_el, bot_el, sea_lvl, real_time_yrs, xmin, xmax, ymin, ymax, delcol, dellay,
                          hk_array, aqt_hk_val, qx1, qz1, n_col, n_lay, dis_package, zoom):

    #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
    fig = plt.figure(figsize = (20, 20))
    ax1 = plt.subplot2grid((2, 2), (0, 1))  #   overall concentration profiles
    ax2 = plt.subplot2grid((2, 2), (1, 1))  #   zoomed in concentration profile
    ax3 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)  #   color bar area 
    ax1.set_position([0.135, 0.55, 0.85, 0.4]) # [left, bottom, width, height]
    ax2.set_position([0.135, 0.1, 0.85, 0.4])
    ax3.set_position([0.04, 0.2, 0.025, 0.6])        
    
    #   adapt the time in years to fit better the folder naming structure
    if real_time_yrs < 10.:
        str_plot_time = '000000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 10. and real_time_yrs < 100.:
        str_plot_time = '00000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 100. and real_time_yrs < 1000.:
        str_plot_time = '0000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 100. and real_time_yrs < 10000.:
        str_plot_time = '000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 1000. and real_time_yrs < 10000.:
        str_plot_time = '00' + str(int(real_time_yrs)) 
    else:
        str_plot_time = '0' + str(int(real_time_yrs)) 

    if zoom:
        xmin_zoom = zoom[0]
        xmax_zoom = zoom[1]
        idx_start_zoom = int(abs((xmin - xmin_zoom) / delcol))
        idx_end_zoom = int(abs((xmin - xmax_zoom) / delcol))
        ymax_zoom = math.ceil(top_el[idx_start_zoom] / dellay) * dellay
        ymin_zoom = math.floor(bot_el[idx_end_zoom] / dellay) * dellay
    
    #   add the location of the aquitard layers    
    mapimg = (hk_array == aqt_hk_val)
       
    #   read in the concentration profile from the dictionary, find the closest time step to the input tume_step_yrs
    title = 'Heads profile and flow pattern at ' + str(int(real_time_yrs)) + ' years since simulation start'

    plot_arr[np.abs(plot_arr) >= 999.] = np.nan

    vmax = math.ceil(np.nanmax(plot_arr))
    vmin = math.floor(np.nanmin(plot_arr))
    
    cmap = plt.cm.plasma
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
    #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)        
    
    bounds = np.linspace(vmin, vmax, 20)
    bounds = [round(x, 1) for x in bounds]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    cbar_title = 'Head elevation (m asl.)'
    
    qx1[np.abs(qx1) == 0.] = np.nan        
    qz1[np.abs(qz1) == 0.] = np.nan        
    
    qx1_avg = np.empty(qx1.shape, dtype = qx1.dtype)
    qz1_avg = np.empty(qz1.shape, dtype = qz1.dtype)
    qx1_avg[:, :, :] = 0.5 * (qx1[:, :, 0 : n_col] + qx1[:, :, : n_col])
    qx1_avg[:, :, 0] = 0.5 * qx1[:, :, 0]
    qz1_avg[:, :, :] = 0.5 * (qz1[0 : n_lay, :, :] + qz1[: n_lay, :, :])
    qz1_avg[0, :, :] = 0.5 * qz1[0, :, :]
    
    y, x, z = dis_package.get_node_coordinates()
    x = np.linspace(xmin, xmax, plot_arr.shape[-1])
    X, Z = np.meshgrid(x, z[:, 0, 0])
    iskip = 3

    #   create output name for the individual figure
    fig_out_dir = os.path.join(dest_dir, '_head_' + str_plot_time + '_yrs.png')           
        
    im1 = ax1.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                   extent = (xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
    ax1.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                   extent = (xmin, xmax, ymin, ymax), alpha = 0.25)
    
    ax2.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                   extent = (xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
    ax2.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                   extent = (xmin, xmax, ymin, ymax), alpha = 0.25)
        
    ax1.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip], qx1_avg[::iskip, 0, ::iskip], -qz1_avg[::iskip, 0, ::iskip],
       color = 'white', scale = 5, headwidth = 2, headlength = 2, headaxislength = 2, width = 0.0015)
    ax2.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip], qx1_avg[::iskip, 0, ::iskip], -qz1_avg[::iskip, 0, ::iskip],
       color = 'white', scale = 5, headwidth = 3, headlength = 2, headaxislength = 2, width = 0.0015)    
    
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
       
    ax2.set_xlim([xmin_zoom, xmax_zoom])
    ax2.set_ylim([ymin_zoom, ymax_zoom])
    
    #   set the gridlines and constant lines in the plot
    #   add constant lines with elevation = 0m asl. and coastline 
    ax1.axhline(y = 0., linewidth = 2, color = 'k', zorder = 2)
    ax1.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)
    ax2.axhline(y = 0., linewidth = 2, color = 'k', zorder = 2)
    ax2.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)    
    
    ax1.axhline(y = sea_lvl, linewidth = 3, color = 'red', zorder = 2)    
    ax2.axhline(y = sea_lvl, linewidth = 3, color = 'red', zorder = 2)    

    x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(xmin / 5.) * 5., math.ceil(xmax), 1000.0), nbins = None)    
    x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(xmin / 1.) * 1., math.ceil(xmax), 250.0), nbins = None)    
    ax1.xaxis.set_major_locator(x_major_locator)
    ax1.xaxis.set_minor_locator(x_minor_locator)
 
    y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(ymin / 500.) * 500., math.ceil(ymax), 250.0), nbins = None)    
    y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(ymin / 100.) * 100. , math.ceil(ymax), 50.0), nbins = None)    
    ax1.yaxis.set_major_locator(y_major_locator)
    ax1.yaxis.set_minor_locator(y_minor_locator)
    
    ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
    ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  

    x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(xmin_zoom / 5.) * 5., math.ceil(xmax_zoom), 1000.0), nbins = None)    
    x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(xmin_zoom / 1.) * 1., math.ceil(xmax_zoom), 250.0), nbins = None)    
    ax2.xaxis.set_major_locator(x_major_locator_zoom)
    ax2.xaxis.set_minor_locator(x_minor_locator_zoom)

    ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
    ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)        
            
    ax1.set_title(title, fontsize = 24, y = 1.025)
    ax1.set_xlabel('distance from coast (km)', fontsize = 18)
    ax1.set_ylabel('elevation (m asl.)', fontsize = 18)
    ax2.set_xlabel('distance from coast (km)', fontsize = 18)
    ax2.set_ylabel('elevation (m asl.)', fontsize = 18)    
    
    #   plot the colorbar
    cbar = plt.colorbar(im1, cax = ax3, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
    cbar.ax.set_title(cbar_title, fontsize = 14, y  = 1.025)
    cbar.ax.tick_params(labelsize = 12)    
    
    #   set the font 
    rcParams['font.family'] = 'Garamond'
    rcParams['axes.facecolor'] = 'white'
    rcParams['savefig.facecolor'] = 'white'    
    
    if os.path.exists(fig_out_dir):
        print('Figure   ' + fig_out_dir + '   already exists.')
    else: 
        #   save the figure
        plt.savefig(fig_out_dir, dpi = 300)
        plt.close()
        del fig


"""
cell_coord_lst = [layer, 0, col_idx]
top_elev_lst = top_elev
head_file = head_arr
hk_array = hk_arr
riv_dict_in = riv_dict
dis_pckg = model.dis
riv_cells = riv_q
drn_cells = drn_q
drn_dict = model.drn_arr_in
rch = rch_q
"""


def get_cell_info(cell_coord_lst, top_elev_lst, head_file, rch, hk_array, riv_dict_in, dis_pckg, riv_cells, drn_cells, drn_dict):
    #   first get the node from the DIS package
    
    node = dis_pckg.get_node((cell_coord_lst[0], 0, cell_coord_lst[-1]))[0] 
    coords = dis_pckg.get_lrc(node)

    elev_val = top_elev_lst[coords[0][-1] - 1]
    head_val = head_file[coords[0][0] - 1, 0, coords[0][-1] - 1]
    hk_val = hk_array[coords[0][0] - 1, 0, coords[0][-1] - 1]
    
    #   check if the node is in the river cells, if yes get all the infor
    for riv_cell in riv_cells[0]:
        if riv_cell[0] == node:
            riv_flow_val = riv_cell[-1]
            coords = dis_pckg.get_lrc(node)
            #   there is some weird shift in the node and riv cell coordinate counting, throws error for the last riv cell..
            #try:            
            riv_cell_info = [subl for subl in riv_dict_in[0] if subl[0] == coords[0][0] - 1 and subl[1] == coords[0][1] - 1 and subl[2] == coords[0][2] - 1][0][3:]
            #except IndexError:
            #    riv_cell_info = [subl for subl in riv_dict[0] if subl[0] == cell_coord_lst[0] and subl[1] == cell_coord_lst[1] and subl[2] == cell_coord_lst[2] - 1][0][3:]
            stage, cond, rbot = round(riv_cell_info[0], 3), round(riv_cell_info[1], 3), round(riv_cell_info[2], 3)            
            break
        else:
            riv_flow_val, stage, cond, rbot = '-', '-', '-', '-'

    rch_flow_val = rch[0][-1][0][coords[0][-1] - 1]
    
    #   check if the node is in the river cells 
    for drn_cell in drn_cells[0]:
        if drn_cell[0] == node:
            drn_flow_val = drn_cell[-1]
            drn_cell_info = [subl for subl in drn_dict[0] if subl[0] == coords[0][0] - 1 and subl[1] == coords[0][1] - 1 and subl[2] == coords[0][2] - 1][0][3:]
            drn_cond = drn_cell_info[-1]
            break
        else:
            drn_flow_val = '-'
            drn_cond = '-'
    
    if drn_flow_val != '-':
        drn_flow_val = round(drn_flow_val, 3)

    if riv_flow_val != '-':
        riv_flow_val = round(riv_flow_val, 3)        

    if rch_flow_val != '-':
        rch_flow_val = round(rch_flow_val, 3)     
        
    cell_text = [[coords[0][0] - 1, 0, coords[0][-1] - 1], hk_val, [drn_cond, cond], stage, rbot, elev_val, round(head_val, 3), riv_flow_val, drn_flow_val, rch_flow_val]
    return cell_text



"""
dest_dir = dest_dir_RIV
riv_id = g
plot_arr = head_arr
ibound_arr = model.ibound_arr
tb_data = tb_data
xmin = model.x_start
xmax = model.x_end 
ymin = zbot
ymax = top
y_min_zoom = top - (layer + 10) * 10.
hk_array = hk_arr
aqt_hk_val = 0.001
qx1 = qx_in
qz1 = qz_in
n_col = model.ncol
n_lay = model.nlay
dis_package = model.dis
"""
#   Plotting of 2D concentration profiles at given time
def plot_RIV_cells(dest_dir, riv_id, plot_arr, ibound_arr, tb_data, real_time_yrs, riv_start, riv_width, riv_buffer, xmin, xmax, ymin, ymax, y_min_zoom,\
                   del_col, del_lay, hk_array, aqt_hk_val, qx1, qz1, n_col, n_lay, dis_package):

    #   extract information for the river cells and their surrounding
    #       1) define the rows of interest
    rows = ['lay, row, col', 'Hk (m/d)', 'DRN/RIV COND (m2/d?)', 'RIV head elev (m)', 'RIV bot elev (m)', 'Elev (m asl.)', 'Head (m asl.)', 'RIV leakage (?)', 'DRN (?)', 'RCH (?)']
    
    #   adapt the time in years to fit better the folder naming structure
    if real_time_yrs < 10.:
        str_plot_time = '000000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 10. and real_time_yrs < 100.:
        str_plot_time = '00000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 100. and real_time_yrs < 1000.:
        str_plot_time = '0000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 100. and real_time_yrs < 10000.:
        str_plot_time = '000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 1000. and real_time_yrs < 10000.:
        str_plot_time = '00' + str(int(real_time_yrs)) 
    else:
        str_plot_time = '0' + str(int(real_time_yrs))     
    
    #   read in the concentration profile from the dictionary, find the closest time step to the input tume_step_yrs
    title = 'Heads profile and fluxes (RIV, RCH, DRN) for time ' + str(real_time_yrs) + ' years' 

    #   adapt the time in years to fit better the folder naming structure
    if real_time_yrs < 10.:
        str_plot_time = '000000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 10. and real_time_yrs < 100.:
        str_plot_time = '00000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 100. and real_time_yrs < 1000.:
        str_plot_time = '0000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 100. and real_time_yrs < 10000.:
        str_plot_time = '000' + str(int(real_time_yrs)) 
    elif real_time_yrs >= 1000. and real_time_yrs < 10000.:
        str_plot_time = '00' + str(int(real_time_yrs)) 
    else:
        str_plot_time = '0' + str(int(real_time_yrs)) 
    
    #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
    fig = plt.figure(figsize = (18, 7.5))
    
    #ax1 = plt.subplot2grid((3, 2), (0, 0), colspan = 2)  #   the title area - mainly with RIV parameters
    ax2 = plt.subplot2grid((3, 2), (0, 1))  #   overall concentration profiles
    ax3 = plt.subplot2grid((3, 2), (1, 1))  #   zoomed in concentration profile
    ax4 = plt.subplot2grid((3, 2), (0, 0), rowspan=2)  #   color bar area 
    ax5 = plt.subplot2grid((3, 2), (2, 0), colspan = 2)  #   table with head values and other output
    
    #ax1.set_position([0.05, 0.9, 0.9, 0.075]) # [left, bottom, width, height]
    ax2.set_position([0.125, 0.6, 0.825, 0.32])
    ax3.set_position([0.125, 0.25, 0.825, 0.32])        
    ax4.set_position([0.05, 0.4, 0.025, 0.35])
    ax5.set_position([0.05, 0.005, 0.9, 0.18])        
    
    
    #   add the location of the aquitard layers    
    #mapimg = (hk_array == aqt_hk_val)
    mapimg = (hk_array <= aqt_hk_val)
       
    #   read in the concentration profile from the dictionary, find the closest time step to the input tume_step_yrs
    #title = 'Heads profile and flow pattern at ' + str(int(real_time_yrs)) + ' years since simulation start'
    
    plot_arr[np.abs(plot_arr) >= 999.] = np.nan
    
    vmax = math.ceil(np.nanmax(plot_arr))
    vmin = math.floor(np.nanmin(plot_arr))
    
    cmap = plt.cm.plasma
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
    #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)        
    
    bounds = np.linspace(vmin, vmax, 20)
    bounds = [round(x, 1) for x in bounds]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    cbar_title = 'Head elevation (m asl.)'
    
    qx1[np.abs(qx1) == 0.] = np.nan        
    qz1[np.abs(qz1) == 0.] = np.nan        
    
    qx1_avg = np.empty(qx1.shape, dtype = qx1.dtype)
    qz1_avg = np.empty(qz1.shape, dtype = qz1.dtype)
    qx1_avg[:, :, :] = 0.5 * (qx1[:, :, 0 : n_col] + qx1[:, :, : n_col])
    qx1_avg[:, :, 0] = 0.5 * qx1[:, :, 0]
    qz1_avg[:, :, :] = 0.5 * (qz1[0 : n_lay, :, :] + qz1[: n_lay, :, :])
    qz1_avg[0, :, :] = 0.5 * qz1[0, :, :]
    
    y, x, z = dis_package.get_node_coordinates()
    x = np.linspace(xmin, xmax, plot_arr.shape[-1])
    X, Z = np.meshgrid(x, z[:, 0, 0])
    iskip = 3
    
    im2 = ax2.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                   extent = (xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
    ax2.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                   extent = (xmin, xmax, ymin, ymax), alpha = 0.25)
    ax2.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip], qx1_avg[::iskip, 0, ::iskip], -qz1_avg[::iskip, 0, ::iskip],
       color = 'white', scale = 10, headwidth = 2, headlength = 2, headaxislength = 2, width = 0.0015)
    
    ax2.set_xlim([xmin, 0.])
    ax2.set_ylim([y_min_zoom, ymax])
    
    #   plot river cells in the zoom in area and its surrounding
    riv_min_x, riv_max_x = riv_start - riv_buffer, riv_start + riv_width + riv_buffer
    
    ax3.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                   extent = (xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
    ax3.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                   extent = (xmin, xmax, ymin, ymax), alpha = 0.25)
    ax3.quiver(X[:, :], Z[:, :], qx1_avg[:, 0, :], -qz1_avg[:, 0, :],
       color = 'white', scale = 10, headwidth = 2, headlength = 2, headaxislength = 2, width = 0.0015)
    
    ax3.set_xlim([riv_min_x / 1000., riv_max_x / 1000.])
    #   find the min max 
    riv_max_col = int(abs(xmin * 1000. // del_col)) - int(abs(riv_max_x // del_col))
    riv_max_cell_lst = list(ibound_arr[:, 0, riv_max_col])
    indices = [i for i, x in enumerate(riv_max_cell_lst) if x == 1]
    y_min_zoom = ymax - indices[-1] * del_lay - del_lay
    ax3.set_ylim([y_min_zoom, ymax])
    
    #   plot the colorbar
    cbar = plt.colorbar(im2, cax = ax4, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)
    cbar.ax.set_title(cbar_title, fontsize = 12, rotation = 90, y = 0.25, x = -1.1)
    cbar.ax.tick_params(labelsize = 10)   
    cbar.ax.yaxis.set_ticks_position('left')
    
    ax2.set_title(title, fontsize = 24, y = 1.025)
    ax3.set_xlabel('distance from coast (km)', fontsize = 12)
    ax2.set_ylabel('elevation (m asl.)', fontsize = 12)
    ax3.set_ylabel('elevation (m asl.)', fontsize = 12)    
    
    ax5.axis("off")
    ax5.table(cellText = list(map(list, zip(*tb_data))), cellLoc = 'center', rowLoc = 'center',\
              rowLabels = rows, colWidths=[0.8 / len(tb_data) for x in tb_data], loc = 'center',\
              bbox = [0.07, 0., 0.915, 1.])

    #   set the font 
    rcParams['font.family'] = 'Garamond'
    rcParams['axes.facecolor'] = 'white'
    rcParams['savefig.facecolor'] = 'white'    

    fig_out_dir = os.path.join(dest_dir, '_RIV_' + str(riv_id) + '_' + str_plot_time + '_yrs.png')   
    
    if os.path.exists(fig_out_dir):
        print('Figure   ' + fig_out_dir + '   already exists.')
    else: 
        #   save the figure
        plt.savefig(fig_out_dir, dpi = 300)
        plt.close()
        del fig







































