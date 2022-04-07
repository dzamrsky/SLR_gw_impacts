# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:31:49 2019

@author: daniel
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 09:20:13 2018

@author: daniel

basemap help/tutorial online

-   https://matplotlib.org/basemap/api/basemap_api.html
-   https://media.readthedocs.org/pdf/basemaptutorial/latest/basemaptutorial.pdf




"""
import matplotlib.pyplot as plt
import numpy as np
import math
import os
#import matplotlib
#from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
from matplotlib.patches import PathPatch
from matplotlib.colors import Normalize
from matplotlib import colors
from matplotlib import colorbar
import matplotlib.cm as cmx
import matplotlib.patches as mpatches
#import matplotlib 
#import pandas as pd
from shapely.geometry import mapping, shape
import fiona
#from scipy.interpolate import UnivariateSpline
from scipy.signal import argrelextrema
import rdp
from collections import Counter
from scipy.stats import norm
from matplotlib import rcParams

import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.ioff()

import sys
sys.path.append(r'g:\Water_Nexus\_A3\_scripts_py3')
import ws_masterscript_py3 as ws_ms
import ws_datatools_py3 as ws               # connect to dbase, run sqls etc.
from modeltools_py3 import cs_model 

model_dir = r'g:\Water_Nexus\_A3\_model_final_v2'
plot_dir = r'g:\Water_Nexus\_A3\_model_final_v2\_plots'

ne_coastline = r'g:\_ORIGINAL_DATA\natural_earth\ne_10m_land'
coscat_polys = r'g:\_CREATED_DATA\_A2_data\coscat_analysis\coscat_dissolved'
f = open('g:\Water_Nexus\_A2\_figures\_representative_profiles_COSCATs\errors.txt','w')
# Read the original GLIM_su Shapefile
input = fiona.open('g:/_CREATED_DATA/_A2_data/GLIM_su_only.shp', 'r')

#table_db = "coastline_coscat"
table_coscat_polys = "coscat_polys"
#id_point = "32806"

del_lay = 10.
del_col = 100.
cs_points_dist = 500.
max_depth = -10000.

#  define all the inputs for the DIS package
nrow = 1                                    # number of rows in model domain (2D so = 1)
nper = 2                                    # number of model stress periods
delc = 1.0                                  # an array of spacings along a column, using the list created earlier (2D so = 1)
laycbd = 0                                  # 0 indicates no confining bed... -> this needs to be checked!!
perlen = [365.25 * 1000, 365.25 * 50]         # an array of the stress period lengths
nstp = [25, 50]                             # number of time steps in each stress period

#   get all glim classes from the original glim_su.shp
all_su_classes_id, all_su_classes_name = [], []
for element in input:
    all_su_classes_id.append(element['properties'].get('id_su'))
    all_su_classes_name.append(str(element['properties'].get('litho_2')))

#all_su_classes_id = sorted(list(set(all_su_classes_id)))
#all_su_classes_name = sorted(list(set(all_su_classes_name)))
all_su_classes_id = sorted(list(set(all_su_classes_id)), key=lambda x: (x is None, x))
all_su_classes_name = sorted(list(set(all_su_classes_name)), key=lambda x: (x is None, x))

#   create colormap for the GLIM classes
cmap = plt.get_cmap('hot')      


#   first create necessary dictionaries (with K values for soil and GLIM classes)
glim_col_dict = dict()
idx = 0
for glim_id in all_su_classes_id:
    glim_col_dict[glim_id] = dict()
    glim_col_dict[glim_id]['glim_su_name'] = all_su_classes_name[idx]
    glim_col_dict[glim_id]['glim_su_id'] = glim_id
    glim_col_dict[glim_id]['plot_cl'] = cmap(float(idx) / float(len(all_su_classes_id)))
    idx += 1

#   function that finds the extent for the map based on the input list of selected coastal points
def get_geoplot_extent(point_id_lst, name_db, host_db, username_db, pass_db, table_db):
    db_connect = ws.dbase_tools(name_db, username_db, host_db, pass_db)
    #   define the min max for both X and Y coordinates
    x_lst, y_lst, pt_class_lst = [], [], []
    #   loop through the list of coastal points and select X and Y coordinates for each of them
    for id_point in point_id_lst:
        #   create a sql to find the coordinates of the point
        sql_xy = "SELECT coscat_id, ST_X(geom), ST_Y(geom), class_fin FROM %s WHERE id_cs = '%s'" % (table_db, id_point)
        db_connect.cur.execute(sql_xy)
        point_coords = db_connect.cur.fetchall()[0]
        x_coord, y_coord = point_coords[1], point_coords[2]
        coscat_id = point_coords[0]
        pt_class = point_coords[3]
        #   append the coordinates to the overall lists
        x_lst.append(x_coord)
        y_lst.append(y_coord)
        pt_class_lst.append(pt_class)
    #   give back all the max and min coordinates
    min_x, max_x, min_y, max_y = min(x_lst), max(x_lst), min(y_lst), max(y_lst)
    print("x_min : " + str(min_x) + "\n" + "x_max : " + str(max_x) + "\n" + "y_min : " + str(min_y) + "\n" + "y_max : " + str(max_y)) 
    return min_x, max_x, min_y, max_y, x_lst, y_lst, int(coscat_id), pt_class_lst

#   get extent of the coscat region
def get_coscat_extent(coscat_id, name_db, host_db, username_db, pass_db, table_db):
    db_connect = ws.dbase_tools(name_db, username_db, host_db, pass_db)
    #   sql to extract the geometry for given coscat region
    sql_xy = "SELECT array_to_string(array_agg, ',') FROM \
                (SELECT array_agg( ST_x(geom)||' '||ST_y(geom)) FROM \
                (SELECT (ST_dumppoints(geom)).geom FROM %s WHERE coscat = '%s' \
                ) AS foo_1\
                ) AS foo_2" % (table_db, coscat_id)
    db_connect.cur.execute(sql_xy)
    point_coords = db_connect.cur.fetchall()[0]
    #   process the whole string into separate x and y lists to get the extent of the coscat region
    x_lst, y_lst = [], []
    str_split = point_coords[0].split(',')
    for string in str_split:
        x_lst.append(string.split(' ')[0])
        y_lst.append(string.split(' ')[1])
    #   output is the min max X and Y and the overall lists
    return float(min(x_lst)), float(max(x_lst)), float(min(y_lst)), float(max(y_lst)), x_lst, y_lst
    

#   fucntion that plots the map with coastal points and GLIM subclasses for given COSCAT region and
#   coastal type of the coastal points (islands, deltas, 2D henry profiles)
def plot_map(id_coscat, out_dir):
    #   create the overall figure
    fig = plt.figure(figsize=(20 , 20))
    ax = fig.add_subplot(111)
    m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = y_min - 1.0, urcrnrlat = y_max + 1.0,\
                llcrnrlon = x_min - 1.0, urcrnrlon = x_max + 1.0, lat_ts = 20, resolution = 'i')
    # draw parallels and meridians.
    m.readshapefile(ne_coastline, 'areas')
    m.readshapefile(clipped_name.split('.')[0], 'glim_su')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90.,91.,15.), labels = [True, True, False, False], dashes = [2,2])
    m.drawmeridians(np.arange(-180.,181.,15.), labels= [False, False, False, True], dashes = [2,2])
    m.drawcountries(linewidth = 2, linestyle = 'solid', color = 'k' ) 
    m.drawrivers(linewidth = 1, linestyle = 'solid', color = 'blue')
    #   add the polygons of the land masses / coastline
    patches = []
    for info, shape_2 in zip(m.areas_info, m.areas):
        patches.append( Polygon(np.array(shape_2), True) )
    ax.add_collection(PatchCollection(patches, facecolor = 'k', alpha = 0.1, edgecolor = 'k', linewidths = 1., zorder = 2))
    #   add the polygons of the land masses / coastline
    patches_id, patches_color = [], []
    idx_1 = 0
    for info, shape_3 in zip(m.glim_su_info, m.glim_su):
        #patches_glim.append( Polygon(np.array(shape_3), True) )
        patch = [Polygon(np.array(shape_3), True)]
        patch_id = m.glim_su_info[idx_1].get('id_su')
        #patch_color = glim_col_dict.keys()[patch_id]
        patch_color = glim_col_dict[patch_id].get('plot_cl')
        patches_color.append(patch_color)
        patches_id.append(patch_id)
        pc = PatchCollection(patch)
        pc.set_facecolor(patch_color)
        pc.set_edgecolor('none') 
        ax.add_collection(pc)
        idx_1 += 1
    patches_id_legend = list(set(patches_id))
    #   plot the coastal points
    x1, y1, x2, y2, x3, y3 = [], [], [], [], [], []
    for i in range(len(x_pt_lst)):
        pt_class = int(cst_pt_type_lst[i])
        if pt_class == 1:
            x1.append(x_pt_lst[i])
            y1.append(y_pt_lst[i])
        elif pt_class == 2:
            x2.append(x_pt_lst[i])
            y2.append(y_pt_lst[i])
        else:
            x3.append(x_pt_lst[i])
            y3.append(y_pt_lst[i])    
    ax.plot(x1, y1, marker = 'o', color = 'cyan', markersize = 5, markeredgewidth = .5, linestyle = 'None', label = 'Delta')
    ax.plot(x2, y2, marker = 'o', color = 'fuchsia', markersize = 5, markeredgewidth = .5, linestyle = 'None', label = 'Henry profile')
    ax.plot(x3, y3, marker = 'o', color = 'yellow', markersize = 5, markeredgewidth = .5, linestyle = 'None', label = 'Other')
    #   define title and plot the legend
    plt.title("COSCAT region nr.: " + str(id_coscat), fontsize = 20)
    handles = []
    for j in patches_id_legend:
        patch_name = glim_col_dict[j].get('glim_su_name')
        patch_col = glim_col_dict[j].get('plot_cl') 
        handles += [mpatches.Patch(color = patch_col, label = patch_name)]
    cs_pt_legend = plt.legend(title = 'Costal point classes', loc = 1, numpoints = 1)
    ax = plt.gca().add_artist(cs_pt_legend)
    plt.legend(title = 'GLIM class', handles = handles, loc = 4, fontsize = 'medium')
    plt.savefig(out_dir, dpi = 300, facecolor = 'w', edgecolor = 'w',
        orientation = 'portrait', papertype = None, format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
        frameon = None)
    plt.close()

 
#   histogram plotting and 
def plot_histogram(input_lst, title, x_axis_label, out_dir):
    # fixed bin size
    #bins = np.arange(0, int(math.ceil(max(input_lst) / 100.0)) * 100, 25) # fixed bin size
    fig1 = plt.figure(figsize = (10, 10))
    ax = fig1.add_subplot(1, 1, 1)
    bins = np.linspace(0, int(math.ceil(max(input_lst) / 100.0)) * 100, 25) # fixed number of bins    
    ax.set_xlim([0, int(math.ceil(max(input_lst) / 100.0)) * 100])
    ax.axvline(x = sum(input_lst)/float(len(input_lst)), linewidth = 2, color = 'red')    
    ax.hist(input_lst, bins = bins, alpha = 0.5)
    ax.set_title(title)
    ax.set_xlabel(x_axis_label)
    ax.set_ylabel('Count')
    plt.savefig(out_dir, dpi = 300, facecolor = 'w', edgecolor = 'w',
            orientation = 'portrait', papertype = None, format = None,
            transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
            frameon = None)
    plt.close()     
    
    
#   histogram plotting and 
def plot_histograms(input_lst_cs_thk, input_lst_cs_width, input_lst_cs_anch_dist, input_lst_anch_depth, out_dir):
    # fixed bin size
    #bins = np.arange(0, int(math.ceil(max(input_lst) / 100.0)) * 100, 25) # fixed bin size
    fig1 = plt.figure(figsize = (20, 20))
    ax1 = fig1.add_subplot(2, 2, 1)
    bins = np.linspace(0, int(math.ceil(max(input_lst_cs_thk) / 100.0)) * 100, 25) # fixed number of bins    
    ax1.set_xlim([0, int(math.ceil(max(input_lst_cs_thk) / 100.0)) * 100])
    ax1.axvline(x = sum(input_lst_cs_thk)/float(len(input_lst_cs_thk)), linewidth = 2, color = 'red')    
    ax1.hist(input_lst_cs_thk, bins = bins, alpha = 0.5)
    ax1.set_title('Coastal thickness')
    ax1.set_xlabel('Thickness (m)')
    ax1.set_ylabel('Count')
    
    ax2 = fig1.add_subplot(2, 2, 2)
    bins = np.linspace(0, int(math.ceil(max(input_lst_cs_width) / 100.0)) * 100, 25) # fixed number of bins    
    ax2.set_xlim([0, int(math.ceil(max(input_lst_cs_width) / 100.0)) * 100])
    ax2.axvline(x = sum(input_lst_cs_width)/float(len(input_lst_cs_width)), linewidth = 2, color = 'red')    
    ax2.hist(input_lst_cs_width, bins = bins, alpha = 0.5)
    ax2.set_title('Coastal plain width')
    ax2.set_xlabel('Width (km)')
    ax2.set_ylabel('Count')    
    
    ax3 = fig1.add_subplot(2, 2, 3)
    bins = np.linspace(0, int(math.ceil(max(input_lst_cs_anch_dist) / 100.0)) * 100, 25) # fixed number of bins    
    ax3.set_xlim([0, int(math.ceil(max(input_lst_cs_anch_dist) / 100.0)) * 100])
    ax3.axvline(x = sum(input_lst_cs_anch_dist)/float(len(input_lst_cs_anch_dist)), linewidth = 2, color = 'red')    
    ax3.hist(input_lst_cs_anch_dist, bins = bins, alpha = 0.5)
    ax3.set_title('Anchor point dist. from. coast')
    ax3.set_xlabel('Distance (km)')
    ax3.set_ylabel('Count')
    
    ax4 = fig1.add_subplot(2, 2, 4)
    input_lst_anch_depth_int = [int(i) for i in input_lst_anch_depth]
    bins = np.linspace(int(math.floor(min(input_lst_anch_depth_int) / 100.0)) * 100, int(math.ceil(max(input_lst_anch_depth_int) / 100.0)) * 100, 25) # fixed number of bins    
    ax4.set_xlim([int(math.floor(min(input_lst_anch_depth_int) / 100.0)) * 100, int(math.ceil(max(input_lst_anch_depth_int) / 100.0)) * 100])
    ax4.axvline(x = sum(input_lst_anch_depth)/float(len(input_lst_anch_depth_int)), linewidth = 2, color = 'red')    
    ax4.hist(input_lst_anch_depth, bins = bins, alpha = 0.5)
    ax4.set_title('Anchor point depth bsl.')
    ax4.set_xlabel('Depth (m)')
    ax4.set_ylabel('Count')      
    
    plt.savefig(out_dir, dpi = 300, facecolor = 'w', edgecolor = 'w',
            orientation = 'portrait', papertype = None, format = None,
            transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
            frameon = None)
    plt.close()       
      
    
    
#   plot representative profile for a coastal type for each COSCAT region
def plot_representative_profile(thk_lst, width_lst, anchor_dist_lst, anchor_depth_lst, wtd_lst, topo_lst, soil_type_lst, soil_thk_lst,\
                                nasa_lst, offshore_lst, plot_name, err_file, id_cosc):
    """
    thk_lst = lst_cs_thk_henry
    width_lst = lst_cs_width_henry
    anchor_dist_lst = lst_anchor_dist_henry
    anchor_depth_lst = lst_anchor_depth_henry
    wtd_lst = lst_wtd_henry
    topo_lst = lst_topo_henry
    soil_type_lst = lst_soil_type_henry
    soil_thk_lst = lst_soil_thk_henry
    nasa_lst = lst_cs_nasa_henry
    offshore_lst = lst_cs_offshore_henry
    plot_name = plot_name_henry
    err_file = f
    id_cosc = id_num_coscat        
    """
    
    #   clean the lists - if the values are = non_val then change to NaN not to screw up the calculation of averages
    arr_wtd = np.array(wtd_lst)
    arr_wtd[arr_wtd <= -999.0] = np.nan     # sometimes there are values of -9999.0 might be an error in the dataset, otherwise nonval = -999
    arr_topo = np.array(topo_lst)
    arr_soil_type = np.array(soil_type_lst)
    arr_soil_thk = np.array(soil_thk_lst)
    arr_soil_thk[arr_soil_thk == -32768.0] = np.nan
    arr_nasa = np.array(nasa_lst) 
    arr_nasa[arr_nasa == -1.0] = np.nan
    arr_offshore = np.array(offshore_lst) # -32768
    arr_offshore[arr_offshore == -32768.0] = np.nan
    
    avg_wtd = np.nanmean(arr_wtd, axis = 0)
    avg_topo = np.nanmean(arr_topo, axis = 0)
    avg_nasa = np.nanmean(arr_nasa, axis = 0)
    avg_offshore = np.nanmean(arr_offshore, axis = 0)
    
    try:
        #max_wtd = np.nanmax(arr_wtd, axis = 0)
        max_topo = np.nanmax(arr_topo, axis = 0)
        #max_nasa = np.nanmax(arr_nasa, axis = 0)
        #max_offshore = np.nanmax(arr_offshore, axis = 0)
        
        #min_wtd = np.nanmin(arr_wtd, axis = 0)
        min_topo = np.nanmin(arr_topo, axis = 0)
        #min_nasa = np.nanmin(arr_nasa, axis = 0)
        #min_offshore = np.nanmin(arr_offshore, axis = 0)    
    except ValueError:
        print('Some error.. who knows whats happening here..')
        
    #   combine the wtd and gebco values - inland take wtd and offshore the gebco
    lst_wtd_topo = []
    wtd_vals_round = [round(i, 1) for i in avg_wtd]
    for a in range(len(wtd_vals_round)):
        if avg_topo[a] >= 0.0:
            lst_wtd_topo.append(avg_topo[a] - wtd_vals_round[a])
        else:
            lst_wtd_topo.append(avg_topo[a])   
    arr_wtd_topo = np.array(lst_wtd_topo)        
    #avg_wtd_topo = np.nanmean(arr_wtd_topo, axis = 0)
    
    end_cst_plain = round(sum(width_lst) / float(len(width_lst)), 1)
    x_nasa = round(sum(anchor_dist_lst) / float(len(anchor_dist_lst)), 1)
    y_nasa = round(sum(anchor_depth_lst) / float(len(anchor_depth_lst)), 1)
    cst_thick_est_avg = round(sum(thk_lst) / float(len(thk_lst)), 1)
    #   cst_thick_est_stdev = round(cs_model_point.sed_thick_est_stdev, 0) 
    
    model = cs_model(33843, avg_topo.tolist(), avg_nasa.tolist(), x_nasa, y_nasa, cst_thick_est_avg, end_cst_plain,\
                     arr_soil_thk.tolist(), arr_soil_type.tolist(), avg_wtd.tolist(), avg_offshore.tolist(), False)
    model.get_top_bot_lst_find_cst(cst_look_up_idx, fos[1])
    model.get_top_bot_col_lst(del_col, cs_points_dist, model.cst_idx, smooth = True)
    model.modflow_bot_elev_list(del_lay, True)


    model.get_top_bot_col_lst(del_col, cs_points_dist, True)
    model.modflow_bot_elev_list(del_lay, True)
    
    try:
        model.modflow_col_list_constant(cs_points_dist, del_col)
        model.modflow_create_IBOUND_arr()
        model.find_topo_divide(10., -1000., del_col)         #   indicate the value of topographic difference to find the divide 
        model.dis_input(nrow, delc, del_col, del_lay, nper, perlen, nstp, laycbd, max_depth, cs_points_dist)
        
        model.top = model.lay_elev[model.lay_idx]            #   top elevation of the model domain - check if it should be like this!
        model.botm = model.lay_elev[model.lay_idx + 1:]   
        #model.zbot = model.botm[-1] 
        zbot = math.floor((model.bot_elev[-1] / 100.0) * 100.0)
        
        #cst_offset_idx = (np.abs(avg_topo - 0.0)).argmin()
        #cst_offset_idx = (np.abs(avg_topo[390:] - 0.0)).argmin()    
        cst_offset_val  = next((x for x in avg_topo[cst_look_up_idx:] if x < 0.0), None)
        cst_offset_idx = avg_topo[cst_look_up_idx:].tolist().index(cst_offset_val)
        cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
        
        print(model.x_start, model.x_end)
        
        ibound_arr = model.ibound_arr * 1.0
        ibound_arr[np.abs(ibound_arr) == 0.] = np.nan
        
        fig1 = plt.figure(figsize = (20, 10))
        ax = fig1.add_subplot(3, 1, 1)
        x_axis = np.linspace(-200.0, 200.0, 801)
        ax.plot(x_axis, arr_topo.T, color = "grey", alpha = 0.25)
        ax.plot(x_axis, arr_wtd_topo.T, linewidth = 2, color = "red")
        ax.plot(x_axis, avg_topo.T, linewidth = 2, color = "blue")
        ax.plot(x_axis, max_topo.T, linewidth = 1, color = "blue")
        ax.plot(x_axis, min_topo.T, linewidth = 1, color = "blue")
        #   round to closest 0.5 to fit the cross-section points..
        nasa_idx = int(2 * (200. - round(x_nasa * 2) / 2))
        ax.plot(0.0 - round(x_nasa * 2) / 2, avg_topo[nasa_idx] - avg_nasa[nasa_idx], 'ko')
        ax.plot(0.0, - cst_thick_est_avg, 'ko')
        ax.axvline(x = - end_cst_plain, linewidth = 2, color = 'k')    
        ax.axvline(x = 0, linewidth = 2, color = 'k')    
        ax.axhline(y = 0, linewidth = 1, color = 'k')    
        
        ax2 = fig1.add_subplot(3, 1, 2)
        x_axis = np.linspace(-200.0, 200.0, 801)
        ax2.plot(x_axis, arr_topo.T, color = "grey", alpha = 0.25)
        ax2.plot(x_axis, arr_wtd_topo.T, linewidth = 2, color = "red")
        ax2.plot(x_axis, avg_topo.T, linewidth = 2, color = "blue")
        ax2.plot(x_axis, max_topo.T, linewidth = 1, color = "blue")
        ax2.plot(x_axis, min_topo.T, linewidth = 1, color = "blue")
        #   round to closest 0.5 to fit the cross-section points..
        nasa_idx = int(2 * (200. - round(x_nasa * 2) / 2))
        ax2.plot(0.0 - round(x_nasa * 2) / 2, avg_topo[nasa_idx] - avg_nasa[nasa_idx], 'ko')
        ax2.plot(0.0, - cst_thick_est_avg, 'ko')
        ax2.axvline(x = - end_cst_plain, linewidth = 2, color = 'k')    
        ax2.axvline(x = 0, linewidth = 2, color = 'k')    
        ax2.axhline(y = 0, linewidth = 1, color = 'k') 
        if cst_offset > 0:
            ax2.set_xlim([model.x_start + cst_offset, model.x_end - cst_offset])
        else:
            ax2.set_xlim([model.x_start - cst_offset, model.x_end + cst_offset])
        ax2.set_ylim([min(model.top_elev), max(model.top_elev)])    
        
        ax3 = fig1.add_subplot(3, 1, 3)
        x_axis = np.linspace(-200.0, 200.0, 801)
        ax3.imshow(ibound_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = 'Set2',
                   extent = (model.x_start + cst_offset, model.x_end + cst_offset, zbot, model.top), vmin = 0., vmax = np.nanmax(ibound_arr)) # 
        ax3.plot(x_axis, arr_wtd_topo.T, linewidth = 2, color = "red")
        ax3.plot(x_axis, avg_topo.T, linewidth = 2, color = "blue")
        ax3.plot(x_axis, max_topo.T, linewidth = 1, color = "blue")
        ax3.plot(x_axis, min_topo.T, linewidth = 1, color = "blue")
        #   round to closest 0.5 to fit the cross-section points..
        nasa_idx = int(2 * (200. - round(x_nasa * 2) / 2))
        ax3.plot(0.0 - round(x_nasa * 2) / 2, avg_topo[nasa_idx] - avg_nasa[nasa_idx], 'ko')
        ax3.plot(0.0, - cst_thick_est_avg, 'ko')
        ax3.axvline(x = - end_cst_plain, linewidth = 2, color = 'k')    
        ax3.axvline(x = 0, linewidth = 2, color = 'k')    
        ax3.axhline(y = 0, linewidth = 1, color = 'k') 
        if cst_offset > 0:
            ax3.set_xlim([model.x_start + cst_offset, model.x_end - cst_offset])
        else:
            ax3.set_xlim([model.x_start - cst_offset, model.x_end + cst_offset])
        ax3.set_ylim([zbot, max(model.top_elev)])   

        plt.savefig(plot_name, dpi = 300, facecolor = 'w', edgecolor = 'w',
            orientation = 'portrait', papertype = None, format = None,
            transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
            frameon = None)
        plt.close()
        
    except IndexError:
        err_file.write("Cannot find coastline - all average elevation values are > 0m asl., id_coscat = " + str(id_cosc))
        #continue

#   Function that creates a global map with shaded COSCAT region 
def coscat_global_map(id_coscat, out_dir):
    #   create the overall figure
    fig = plt.figure(figsize=(10 , 10))
    ax = fig.add_subplot(111)
    m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -90.0, urcrnrlat = 90.0,\
                llcrnrlon = -180.0, urcrnrlon = 180.0, lat_ts = 20, resolution = 'i')
    # draw parallels and meridians.
    #m.readshapefile(ne_coastline, 'areas')
    m.readshapefile(coscat_polys, 'coscat')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90.,91.,15.), labels = [True, True, False, False], dashes = [2,2])
    m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], dashes = [2,2])
    patches   = []
    for info, shape_5 in zip(m.coscat_info, m.coscat):
        if info['COSCAT'] == id_coscat:
            patches.append( Polygon(np.array(shape_5), True) )
    ax.add_collection(PatchCollection(patches, facecolor= 'red', edgecolor='red', linewidths = 2, zorder=2))
    plt.title("COSCAT region nr.: " + str(id_coscat), fontsize = 20)
    plt.savefig(out_dir, dpi = 300, facecolor = 'w', edgecolor = 'w',
        orientation = 'portrait', papertype = None, format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
        frameon = None)
    plt.close()  



#   function that creates the .csv file with the different model scenarios
#def create_model_scenarios_table(id_coscat, input_lst, out_dir): 


#   combine two sets of topographical data - inland and offshore, with a buffer around the coastline to interpolate
#   the elevations if there is a mismatch between e.g. maximum topography inland and minimum offshore..
#   cst_buff_width = 10 means 10km
def combine_topo(inland_topo_lst, offshore_topo_lst, cst_buff_width = 10):
    
    #   first trim the input lists to only the buffer part around the coast
    in_topo_lst = inland_topo_lst[400 - 2 * cst_buff_width : 400 + 2 * cst_buff_width]
    off_topo_lst = offshore_topo_lst[400 - 2 * cst_buff_width : 400 + 2 * cst_buff_width]
    
    #   take the first value from in_topo_lst and last one from off_topo_lst as end values
    #   the rest will be based on % difference going from 100% inland to 100% offshore
    inter_lst = [in_topo_lst[0]]
    pct_diff = 100. / (float(cst_buff_width * 4) - 1)
    for i in range(1, cst_buff_width * 4 - 1):
        diff = in_topo_lst[i] - off_topo_lst[i]
        inter_val = in_topo_lst[i] - (i * pct_diff * diff) / 100.
        inter_lst.append(inter_val)
        #print diff, i, i * pct_diff, in_topo_lst[i], off_topo_lst[i], in_topo_lst[i] - (i * pct_diff * diff) / 100.
    inter_lst.append(off_topo_lst[-1])
    #   create the final output list
    fin_lst = inland_topo_lst[: 400 - 2 * cst_buff_width] + inter_lst + offshore_topo_lst[400 + 2 * cst_buff_width :]
    return fin_lst 




#   function that saves the individual model inputs into a specific folder
def save_model_input_files_cst_type_topo(cst_type, thk_lst, width_lst, anchor_dist_lst, anchor_depth_lst, wtd_lst, topo_lst,\
                                         soil_type_lst, soil_thk_lst, nasa_lst, offshore_lst, pcr_lst, watergap_lst, p_min_et_lst,\
                                         k_soil_lst, drn_rate_lst, riv_cond_01_lst, riv_cond_05_lst, riv_cond_10_lst, riv_width_lst,\
                                         riv_bot_elev_lst, riv_head_elev_lst, glhymps_top_lay_lst, glhymps_bot_lay_lst, id_coscat, out_dir):

    """
    thk_lst = lst_cs_thk_henry
    width_lst = lst_cs_width_henry
    anchor_dist_lst = lst_anchor_dist_henry
    anchor_depth_lst = lst_anchor_depth_henry
    wtd_lst = lst_wtd_henry
    topo_lst = lst_topo_henry
    soil_type_lst = lst_soil_type_henry
    soil_thk_lst = lst_soil_thk_henry
    nasa_lst = lst_cs_nasa_henry
    offshore_lst = lst_cs_offshore_henry
    err_file = f
    id_cosc = id_num_coscat   
    out_dir = coscat_dir    
    cst_type = 'henry'
    pcr_lst = lst_cs_pcr_rch_henry
    watergap_lst = lst_cs_watergap_rch_henry
    p_min_et_lst = lst_cs_p_min_et_henry
    k_soil_lst = lst_cs_k_soil_henry
    drn_rate_lst = lst_cs_drn_rate_henry
    riv_cond_01_lst = lst_cs_riv_cond_01_henry
    riv_cond_05_lst = lst_cs_riv_cond_05_henry
    riv_cond_10_lst = lst_cs_riv_cond_10_henry
    riv_width_lst = lst_cs_riv_width_henry
    riv_bot_elev_lst = lst_cs_riv_bot_elev_henry
    riv_head_elev_lst = lst_cs_riv_head_elev_henry
    glhymps_top_lay_lst = lst_cs_glhymps_top_lay_henry
    glhymps_bot_lay_lst = lst_cs_glhymps_bot_lay_henry
    id_coscat = id_num_coscat
    out_dir = coscat_dir      
    topo_diff = 10.
    bot_limit = -10000.
    ibound_act_cells_limit = 2
    offshore_dist_limit = 200000.

    thk_lst = lst_cs_thk_henry_other
    width_lst = lst_cs_width_henry_other
    anchor_dist_lst = lst_anchor_dist_henry_other
    anchor_depth_lst = lst_anchor_depth_henry_other
    wtd_lst = lst_wtd_henry_other
    topo_lst = lst_topo_henry_other
    soil_type_lst = lst_soil_type_henry_other
    soil_thk_lst = lst_soil_thk_henry_other
    nasa_lst = lst_cs_nasa_henry_other
    offshore_lst = lst_cs_offshore_henry_other
    err_file = f
    id_cosc = id_num_coscat   
    out_dir = coscat_dir 
    cst_type = 'henry_other'
    pcr_lst = lst_cs_pcr_rch_henry_other
    watergap_lst = lst_cs_watergap_rch_henry_other
    p_min_et_lst = lst_cs_p_min_et_henry_other
    k_soil_lst = lst_cs_k_soil_henry_other
    drn_rate_lst = lst_cs_drn_rate_henry_other
    riv_cond_01_lst = lst_cs_riv_cond_01_henry_other
    riv_cond_05_lst = lst_cs_riv_cond_05_henry_other
    riv_cond_10_lst = lst_cs_riv_cond_10_henry_other
    riv_width_lst = lst_cs_riv_width_henry_other
    riv_bot_elev_lst = lst_cs_riv_bot_elev_henry_other
    riv_head_elev_lst = lst_cs_riv_head_elev_henry_other
    glhymps_top_lay_lst = lst_cs_glhymps_top_lay_henry_other
    glhymps_bot_lay_lst = lst_cs_glhymps_bot_lay_henry_other
    id_coscat = id_num_coscat
    out_dir = coscat_dir    
    topo_diff = 10.
    bot_limit = -10000.
    ibound_act_cells_limit = 2
    offshore_dist_limit = 200000.
    
    
    thk_lst = lst_cs_thk_other
    width_lst = lst_cs_width_other
    anchor_dist_lst = lst_anchor_dist_other
    anchor_depth_lst = lst_anchor_depth_other
    wtd_lst = lst_wtd_other
    topo_lst = lst_topo_other
    soil_type_lst = lst_soil_type_other
    soil_thk_lst = lst_soil_thk_other
    nasa_lst = lst_cs_nasa_other
    offshore_lst = lst_cs_offshore_other
    err_file = f
    id_cosc = id_num_coscat   
    out_dir = coscat_dir   
    cst_type = 'other'
    pcr_lst = lst_cs_pcr_rch_other
    watergap_lst = lst_cs_watergap_rch_other
    p_min_et_lst = lst_cs_p_min_et_other
    k_soil_lst = lst_cs_k_soil_other
    drn_rate_lst = lst_cs_drn_rate_other
    riv_cond_01_lst = lst_cs_riv_cond_01_other
    riv_cond_05_lst = lst_cs_riv_cond_05_other
    riv_cond_10_lst = lst_cs_riv_cond_10_other
    riv_width_lst = lst_cs_riv_width_other
    riv_bot_elev_lst = lst_cs_riv_bot_elev_other
    riv_head_elev_lst = lst_cs_riv_head_elev_other
    glhymps_top_lay_lst = lst_cs_glhymps_top_lay_other
    glhymps_bot_lay_lst = lst_cs_glhymps_bot_lay_other
    id_coscat = id_num_coscat
    out_dir = coscat_dir   

    thk_lst = lst_cs_thk_delta
    width_lst = lst_cs_width_delta
    anchor_dist_lst = lst_anchor_dist_delta
    anchor_depth_lst = lst_anchor_depth_delta
    wtd_lst = lst_wtd_delta
    topo_lst = lst_topo_delta
    soil_type_lst = lst_soil_type_delta
    soil_thk_lst = lst_soil_thk_delta
    nasa_lst = lst_cs_nasa_delta
    offshore_lst = lst_cs_offshore_delta
    err_file = f
    id_cosc = id_num_coscat   
    out_dir = coscat_dir    
    cst_type = 'delta'
    pcr_lst = lst_cs_pcr_rch_delta
    watergap_lst = lst_cs_watergap_rch_delta
    p_min_et_lst = lst_cs_p_min_et_delta
    k_soil_lst = lst_cs_k_soil_delta
    drn_rate_lst = lst_cs_drn_rate_delta
    riv_cond_01_lst = lst_cs_riv_cond_01_delta
    riv_cond_05_lst = lst_cs_riv_cond_05_delta
    riv_cond_10_lst = lst_cs_riv_cond_10_delta
    riv_width_lst = lst_cs_riv_width_delta
    riv_bot_elev_lst = lst_cs_riv_bot_elev_delta
    riv_head_elev_lst = lst_cs_riv_head_elev_delta
    glhymps_top_lay_lst = lst_cs_glhymps_top_lay_delta
    glhymps_bot_lay_lst = lst_cs_glhymps_bot_lay_delta
    id_coscat = id_num_coscat
    #out_dir = coscat_dirs[0]   
    """

    last_ibound_col_n_layers = -55
    ls_elev = -55
    ls_elev_50 = -55   
    fos_pt_dist = -55
    
    #   clean the lists - if the values are = non_val then change to NaN not to screw up the calculation of averages
    arr_wtd = np.array(wtd_lst)
    arr_wtd[arr_wtd <= -999.0] = np.nan     # sometimes there are values of -9999.0 might be an error in the dataset, otherwise nonval = -999
    arr_topo = np.array(topo_lst)
    arr_soil_type = np.array(soil_type_lst)
    arr_soil_thk = np.array(soil_thk_lst)
    arr_soil_thk[arr_soil_thk <= -9999.0] = np.nan  # -32768.0
    arr_nasa = np.array(nasa_lst) 
    arr_nasa[arr_nasa == -1.0] = np.nan
    arr_offshore = np.array(offshore_lst) # -32768
    arr_offshore[arr_offshore == -32768.0] = np.nan
    
    arr_pcr_rch = np.array(pcr_lst)
    arr_pcr_rch[arr_pcr_rch == 0.0] = np.nan
    arr_watergap_rch = np.array(watergap_lst)
    arr_watergap_rch[arr_watergap_rch == 0.0] = np.nan
    arr_p_min_et = np.array(p_min_et_lst)
    arr_k_soil = np.array(k_soil_lst)
    arr_k_soil[arr_k_soil <= -9999.] = np.nan
    arr_drn_rate = np.array(drn_rate_lst)   
    arr_drn_rate[arr_drn_rate <= -9999.] = np.nan
    
    arr_riv_cond_01 = np.array(riv_cond_01_lst)   
    arr_riv_cond_05 = np.array(riv_cond_05_lst)   
    arr_riv_cond_10 = np.array(riv_cond_10_lst)   
    arr_riv_cond_01[arr_riv_cond_01 <= -9999.] = np.nan    
    arr_riv_cond_05[arr_riv_cond_05 <= -9999.] = np.nan    
    arr_riv_cond_10[arr_riv_cond_10 <= -9999.] = np.nan    
    arr_riv_width = np.array(riv_width_lst)   
    arr_riv_width[arr_riv_width <= -9999.] = np.nan       
    arr_riv_bot_elev = np.array(riv_bot_elev_lst)   
    arr_riv_bot_elev[arr_riv_bot_elev <= -9999.] = np.nan    
    arr_riv_head_elev = np.array(riv_head_elev_lst)   
    arr_riv_head_elev[arr_riv_head_elev <= -9999.] = np.nan   
    arr_glhymps_top_lay = np.array(glhymps_top_lay_lst)   
    arr_glhymps_bot_lay = np.array(glhymps_bot_lay_lst)   
    arr_glhymps_top_lay[arr_glhymps_top_lay <= -9999.] = np.nan   
    arr_glhymps_bot_lay[arr_glhymps_bot_lay <= -9999.] = np.nan   
        
    avg_wtd = np.nanmean(arr_wtd, axis = 0)
    avg_topo = np.nanmean(arr_topo, axis = 0)
    avg_nasa = np.nanmean(arr_nasa, axis = 0)
    avg_offshore = np.nanmean(arr_offshore, axis = 0)
    try:
        avg_pcr_rch = np.nanmean(arr_pcr_rch, axis = 0)
    except TypeError:
        arr_pcr_rch_data = np.array(arr_pcr_rch, dtype='float')
        avg_pcr_rch = np.nanmean(arr_pcr_rch_data, axis = 0)    
    avg_watergap_rch = np.nanmean(arr_watergap_rch, axis = 0)
    avg_p_min_et = np.nanmean(arr_p_min_et, axis = 0)
    try:
        avg_k_soil = np.nanmean(arr_k_soil, axis = 0)
    except ZeroDivisionError:
        arr_k_soil_data = np.array(arr_k_soil, dtype='float')
        avg_k_soil = np.nanmean(arr_k_soil_data, axis = 0)
        
    avg_drn_rate = np.nanmean(arr_drn_rate, axis = 0)
    avg_soil_thk = np.nanmean(arr_soil_thk, axis = 0)
    avg_soil_type = []
    
    avg_riv_cond_01 = np.nanmean(arr_riv_cond_01, axis = 0)    
    avg_riv_cond_05 = np.nanmean(arr_riv_cond_05, axis = 0)    
    avg_riv_cond_10 = np.nanmean(arr_riv_cond_10, axis = 0)    
    avg_riv_width = np.nanmean(arr_riv_width, axis = 0)    
    avg_riv_bot_elev = np.nanmean(arr_riv_bot_elev, axis = 0)    
    avg_riv_head_elev = np.nanmean(arr_riv_head_elev, axis = 0)        
    avg_glhymps_top_lay = np.nanmean(arr_glhymps_top_lay, axis = 0)        
    avg_glhymps_bot_lay = np.nanmean(arr_glhymps_bot_lay, axis = 0)            
    
    for i in range(avg_topo.shape[0]):
        if avg_topo[i] >= 0:
            col_lst = arr_soil_type[:, i].tolist()
            col_lst2 = [x for x in col_lst if x != -9999]
            cnt = Counter(col_lst2)
            avg_soil_type.append(cnt.most_common(1)[0][0])
        else:
            avg_soil_type.append(-1)

        
    #   combine the wtd and gebco values - inland take wtd and offshore the gebco
    lst_wtd_topo = []
    wtd_vals_round = [round(i, 1) for i in avg_wtd]
    for a in range(len(wtd_vals_round)):
        if avg_topo[a] >= 0.0:
            lst_wtd_topo.append(avg_topo[a] - wtd_vals_round[a])
        else:
            lst_wtd_topo.append(avg_topo[a])   
    arr_wtd_topo = np.array(lst_wtd_topo)        
    #avg_wtd_topo = np.nanmean(arr_wtd_topo, axis = 0)
    
    end_cst_plain = round(sum(width_lst) / float(len(width_lst)), 1)
    x_nasa = round(sum(anchor_dist_lst) / float(len(anchor_dist_lst)), 1)
    #y_nasa = round(sum(anchor_depth_lst) / float(len(anchor_depth_lst)), 1)
    cst_thick_est_avg = round(sum(thk_lst) / float(len(thk_lst)), 1)
        
    #   cst_thick_est_stdev = round(cs_model_point.sed_thick_est_stdev, 0) 
    
    #   list with all the topographical combinations inland_offshore
    scenarios = ['avg_avg']#['min_min', 'min_avg', 'min_max', 'avg_min', 'avg_avg', 'avg_max', 'max_min', 'max_avg', 'max_max']
    
    topo_dict = dict()

    for topo_sc in scenarios:
        topo_dict[topo_sc] = dict()
        topo_dict[topo_sc]['ibound_arr'] = []
        topo_dict[topo_sc]['top_elev'] = []
        topo_dict[topo_sc]['bot_elev'] = []
        topo_dict[topo_sc]['lay_elev'] = []
        topo_dict[topo_sc]['cst_offset'] = []
        topo_dict[topo_sc]['lay_idx'] = []
        topo_dict[topo_sc]['col_idx'] = []
        topo_dict[topo_sc]['end_idx'] = []
        topo_dict[topo_sc]['top'] = []
        topo_dict[topo_sc]['zbot'] = []
        topo_dict[topo_sc]['col_width'] = []
        topo_dict[topo_sc]['idx_start'] = []
        topo_dict[topo_sc]['idx_end'] = []
        topo_dict[topo_sc]['x_start'] = []
        topo_dict[topo_sc]['x_end'] = []
        topo_dict[topo_sc]['soil_thk'] = []
        topo_dict[topo_sc]['soil_type'] = []
        topo_dict[topo_sc]['cst_idx'] = []
        topo_dict[topo_sc]['wtd_elev'] = []
        topo_dict[topo_sc]['pcr_rch'] = []
        topo_dict[topo_sc]['watergap_rch'] = []
        topo_dict[topo_sc]['p_min_et'] = []
        topo_dict[topo_sc]['k_soil'] = []
        topo_dict[topo_sc]['drn_rate'] = []
        topo_dict[topo_sc]['riv_cond_01'] = []
        topo_dict[topo_sc]['riv_cond_05'] = []
        topo_dict[topo_sc]['riv_cond_10'] = []
        topo_dict[topo_sc]['riv_width'] = []        
        topo_dict[topo_sc]['riv_bot_elev'] = []
        topo_dict[topo_sc]['riv_head_elev'] = []        
        topo_dict[topo_sc]['glhymps_top_lay'] = []
        topo_dict[topo_sc]['glhymps_bot_lay'] = []
    
    #   create an overall figure where all the topographical models will be stored
    fig_all = plt.figure(figsize = (20, 10))
    
    in_lst = avg_topo.tolist()
    off_lst = avg_topo.tolist()
    final_topo = avg_topo.tolist()

    #   first check the topographical profiles both offshore and inland - if offshore all > 0 then dont model
    #   same case is for inland all < 0
    off_pos = [i for i in off_lst if i < 0]
    in_pos = [i for i in in_lst if i > 0]
        
    if off_pos == [] or len([x for x in off_pos if x > -5]) == len(off_pos):
        name = 'no_data'
        ax_all = fig_all.add_subplot(1, 1, 1)
        ax_all.text(0.15, 0.5, r'Cant build model, all topography offshore > 0')
        ax_all.set_title(name)
        print("Cannot find coastline - all average offshore elevation values are > 0m asl., id_coscat = ") 
        last_ibound_col_n_layers = -2
        ls_elev = -2
        ls_elev_50 = -2   
        fos_pt_dist = -2
        #return None, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist        

    
    elif in_pos == []:
        name = 'no_data'
        ax_all = fig_all.add_subplot(1, 1, 1)
        ax_all.text(0.15, 0.5, r'Cant build model, all topography offshore > 0')
        ax_all.set_title(name)
        print("Cannot find coastline - all average inland elevation values are < 0m asl., id_coscat = ") 
        last_ibound_col_n_layers = -2
        ls_elev = -2
        ls_elev_50 = -2   
        fos_pt_dist = -2
        #return None, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist    
    
    else:
        #   define the index from which we will look for the coastline 
        cst_look_up_idx = 200
        #   find the position of the coastline
        cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < 0.0), None)
        #   check if the coastal offset exists, if not it is probably beacuse the whole profile is above sea level
        if cst_offset_val is None:
            name = 'no_data'
            ax_all = fig_all.add_subplot(1, 1, 1)
            ax_all.text(0.15, 0.5, r'Cant build model, all topography offshore > 0')
            ax_all.set_title(name)
            print("Cannot find coastline - all average elevation values are > 0m asl., id_coscat = ")# + str(id_cosc))            
            last_ibound_col_n_layers = -3
            ls_elev = -3
            ls_elev_50 = -3   
            fos_pt_dist = -3
            #return None, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist    
        
        else:
            print("Pica") 
            cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
            cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
            cst_off_idx = cst_look_up_idx + cst_offset_idx
    
            if cst_off_idx < cst_look_up_idx:
                print("Pica 2")                 
                #   define the index from which we will look for the coastline 
                cst_look_up_idx = 300
                #   find the position of the coastline
                cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < 0.0), None)
                cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                cst_off_idx = cst_look_up_idx + cst_offset_idx
                
    
                if cst_off_idx < cst_look_up_idx:
                    print("Pica 3")                         
                    
                    #   define the index from which we will look for the coastline 
                    cst_look_up_idx = 350
                    #   find the position of the coastline
                    cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < 0.0), None)
                    cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                    cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                    cst_off_idx = cst_look_up_idx + cst_offset_idx                
                    last_ibound_col_n_layers = -3
                    #return None, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist    

                    if cst_off_idx < cst_look_up_idx:
                        print("Pica 3a")           
    
                        #   define the index from which we will look for the coastline 
                        cst_look_up_idx = 380
                        #   find the position of the coastline
                        cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < 0.0), None)
                        cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                        cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                        cst_off_idx = cst_look_up_idx + cst_offset_idx        
                        
                        if cst_off_idx < cst_look_up_idx:
                            print("Pica 3b the island is too small ty kundo")    
                            ls_elev = -3
                            ls_elev_50 = -3   
                            fos_pt_dist = -3
                            model_run = False   
                            
                        else:
                            model_run = True
                    else:
                        model_run = True
                else:
                    model_run = True
            else:
                model_run = True
            
            if model_run == True:
                print("Pica 4")    
                #   the Anchor point is defined as the point where it is first time indicated that the regolith
                #   thickness is equal or more than 50m, therefore put the anchor point at topo_elev - 50m
                anchor_topo_elev = final_topo[int(round((cst_off_idx - 2 * x_nasa) * 2) / 2)]
                y_nasa = anchor_topo_elev - 50.0    
            
                try: 
                    #   first try to find the shelf break and foot of continental slope points
                    print("Pica 5")    
                    offshore_dir = os.path.join(out_dir, cst_type)
                    if not os.path.exists(offshore_dir):
                        os.makedirs(offshore_dir)
                    fos = find_FOS(final_topo, offshore_dir, topo_sc, avg= False) 
                    print(fos) 
                    model = cs_model(1, final_topo, avg_nasa.tolist(), x_nasa, y_nasa, cst_thick_est_avg, end_cst_plain,\
                                     avg_soil_thk.tolist(), arr_soil_type.tolist(), arr_wtd_topo.tolist(), avg_offshore.tolist(),\
                                     avg_pcr_rch.tolist(), avg_watergap_rch.tolist(), avg_p_min_et.tolist(), avg_k_soil.tolist(),\
                                     avg_drn_rate.tolist(), avg_riv_cond_01.tolist(), avg_riv_cond_05.tolist(), avg_riv_cond_10.tolist(),\
                                     avg_riv_width.tolist(), avg_riv_bot_elev.tolist(), avg_riv_head_elev.tolist(),\
                                     avg_glhymps_top_lay.tolist(), avg_glhymps_bot_lay.tolist(), False)
                    
                    model.get_top_bot_lst_find_cst(cst_look_up_idx, fos[1], 0.0)
                    model.get_top_bot_col_lst(del_col, cs_points_dist, model.cst_idx, smooth = True)
                    model.get_top_bot_col_lst_top_sys_geo(del_col, cs_points_dist, model.cst_idx, smooth = True)
                    model.modflow_bot_elev_list(del_lay, True)

                except AttributeError:
                    print("Pica 6")   
                    name = 'no_data'
                    ax_all = fig_all.add_subplot(1, 1, 1)
                    ax_all.text(0.15, 0.5, r'Cant build model, all topography offshore > 0')
                    ax_all.set_title(name)                    
                    last_ibound_col_n_layers = -4
                    ls_elev = -4
                    ls_elev_50 = -4   
                    fos_pt_dist = -4
                    #return None, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist

                try:
                    print("Pica 7")  
                    model.modflow_col_list_constant(cs_points_dist, del_col)
                    model.modflow_create_IBOUND_arr()                    
                    try:
                        model.find_topo_divide(10., -10000., del_col, model.top_elev[model.idx_start], fos[1], 2, 200000.)         #   indicate the value of topographic difference to find the divide 
                    #   in case the island is small cut off the inland part with the start of the model.top_elev list
                    except IndexError:
                        model.find_topo_divide(10., -10000., del_col, model.top_elev[0], fos[1], 2, 200000.)  

                    model.dis_input(nrow, delc, del_col, del_lay, nper, perlen, nstp, laycbd, max_depth, cs_points_dist)
                    model.top = model.lay_elev[model.lay_idx]            #   top elevation of the model domain - check if it should be like this!
                    model.botm = model.lay_elev[model.lay_idx + 1:]   
                    #model.zbot = model.botm[-1] 
                    zbot = math.floor((model.bot_elev[-1] / 100.0) * 100.0)
                                 
                    ibound_arr = model.ibound_arr * 1.0
                    ibound_arr[np.abs(ibound_arr) == 0.] = np.nan
                    
                    str_0 = cst_type + "\\" + topo_sc + '.png'
                    plot_dir = os.path.join(out_dir, str_0)           
                    
                    x_start_new = model.x_start# + (model.cst_idx * del_col) / 1000.
                    x_shift = round(model.x_start + (model.cst_idx * del_col) / 1000., 1)
                    x_end_new = model.x_end # + (model.cst_idx * del_col) / 1000.
                           
                    #model_start_idx = final_topo.index(min(final_topo, key=lambda x:abs(x-model.top_elev[0])))   
                    #model_end_idx = final_topo.index(min(final_topo, key=lambda x:abs(x-model.top_elev[-1])))  
                           
                    model_end_idx = int(round(model.idx_end * del_col / 1000.)) * 2 + cst_off_idx + int(round(x_start_new * 2))       
                    #model_start_idx = 400 + int(round(x_start_new * 2))                           
                           
                    topo_lst_plot = final_topo[cst_off_idx + int(round(x_start_new * 2)): model_end_idx + 1]

                    #x_axis_plot2 = np.linspace(x_start_new, x_end_new, model.ibound_arr.shape[-1])# - model.col_idx)
                    #x_axis_plot = np.linspace(x_start_new, x_end_new, len(topo_lst_plot))

                    if cst_offset >= 0:
                        x_axis_plot2 = np.linspace(x_start_new - x_shift, x_end_new + x_shift, model.ibound_arr.shape[-1])# - model.col_idx)
                        x_axis_plot = np.linspace(x_start_new - x_shift, x_end_new + x_shift, len(topo_lst_plot))
                    else:
                        x_axis_plot2 = np.linspace(x_start_new - x_shift, x_end_new + x_shift, model.ibound_arr.shape[-1])# - model.col_idx)
                        x_axis_plot = np.linspace(x_start_new - x_shift, x_end_new + x_shift, len(topo_lst_plot))

                           
                    """
                    #model_end_idx = 401 + int((model.x_end + cst_offset) * 2)           
                    model_start_idx = 400 + int((model.x_start) * 2)   # CH  + cst_offset
                    if model.x_end == 200.0:
                        model_end_idx = model_start_idx + len(model.top_elev) / 5
                        #x_axis_plot2 = np.linspace(model.x_start, model.x_end, model.ibound_arr.shape[-1] - model.col_idx)
                    else:
                        model_end_idx = 400 + int(round((model.x_end) * 2))   # CH  + cst_offset
                        #x_axis_plot2 = np.linspace(model.x_start, model.x_end, 5 * (model_end_idx - model_start_idx))
                        #x_axis_plot2 = np.linspace(model.x_start, model.x_end, model.ibound_arr.shape[-1] - model.col_idx)
                    topo_lst_plot = final_topo[model_start_idx + (model.col_idx / 5) : model_end_idx]

                    if cst_offset > 0:
                        x_axis_plot2 = np.linspace(model.x_start - cst_offset, model.x_end - cst_offset, model.ibound_arr.shape[-1] - model.col_idx)
                        x_axis_plot = np.linspace(model.x_start - cst_offset, model.x_end - cst_offset, len(topo_lst_plot))
                    else:
                        x_axis_plot2 = np.linspace(model.x_start - cst_offset, model.x_end - cst_offset, model.ibound_arr.shape[-1] - model.col_idx)
                        x_axis_plot = np.linspace(model.x_start - cst_offset, model.x_end - cst_offset, len(topo_lst_plot))
                    """
                    
                    
                    fig1 = plt.figure(figsize = (20, 10))
                    ax = fig1.add_subplot(1, 1, 1)
                    
                    print("x_axis_plot, topo_lst_plot extent", x_axis_plot[0], x_axis_plot[-1], topo_lst_plot[0], topo_lst_plot[-1]) 
                    print("x_axis_plot2, model.top_elev", x_axis_plot2[0], x_axis_plot2[-1], model.top_elev[0], model.top_elev[-1])
                    
                
                    ax.imshow(ibound_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = 'Set2',
                               extent = (x_axis_plot[0], x_axis_plot[-1], model.botm[-1], model.top), vmin = 0., vmax = np.nanmax(ibound_arr))
                    
                    #ax.plot(x_axis_plot, topo_lst_plot, linewidth = 2, color = "blue")
                    ax.plot(x_axis_plot2, model.top_elev[:], linewidth = 2, color = "red")
                    ax.plot(x_axis_plot2, model.bot_elev[:], linewidth = 2, color = "green")                    
                    
                    if fos[0] is not None and cst_offset >= 0:
                        ax.plot(fos[0][0], fos[0][1], 'bo')
                    elif fos[0] is not None and cst_offset <= 0:
                        ax.plot(fos[0][0] + x_shift, fos[0][1], 'bo')
                        
                    ax.axhline(y = 0, linewidth = 1, color = 'k') 
                    ax.axvline(x = 0, linewidth = 2, color = 'k')
 
                    ax.plot(- x_nasa, y_nasa, 'ko')    #   CH
                    ax.plot(0.0, - cst_thick_est_avg, 'ko')   #   CH 
                   
                    if cst_offset >= 0:
                        ax.set_xlim([x_start_new - x_shift, x_end_new + x_shift])
                    else:
                        ax.set_xlim([x_start_new - x_shift, x_end_new + x_shift])
                    ax.set_ylim([math.floor(zbot / 100.) * 100, max(model.top_elev) + 50.0])                       
                    
                    
                    """
                    ax.plot(x_axis_plot2, model.top_elev[model.col_idx : model.ibound_arr.shape[-1]], linewidth = 1, color = "red")
                    ax.plot(x_axis_plot2, model.bot_elev[model.col_idx : model.ibound_arr.shape[-1]], linewidth = 2, color = "green")
                    ax.plot(- x_nasa, y_nasa, 'ko')    #   CH
                    ax.plot(0.0, - cst_thick_est_avg, 'ko')   #   CH 
                    ax.axvline(x = - end_cst_plain, linewidth = 2, color = 'k') 
                    ax.axvline(x = 0, linewidth = 2, color = 'k') 

                    if fos[0] is not None and cst_offset > 0:
                        ax.plot(fos[0][0] - cst_offset, fos[0][1], 'bo')
                    elif fos[0] is not None and cst_offset <= 0:
                        ax.plot(fos[0][0] - cst_offset, fos[0][1], 'bo')
                    ax.axhline(y = 0, linewidth = 1, color = 'k') 
                    if cst_offset > 0:
                        ax.set_xlim([model.x_start - cst_offset, model.x_end - cst_offset])
                    else:
                        ax.set_xlim([model.x_start - cst_offset, model.x_end - cst_offset])
                    ax.set_ylim([zbot - 500.0, max(model.top_elev) + 50.0])   
                    """
                    
                    print(plot_dir)    
                    plt.savefig(plot_dir, dpi = 300, facecolor = 'w', edgecolor = 'w',
                        orientation = 'portrait', papertype = None, format = None,
                        transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
                        frameon = None)     
                    plt.close()
    
                    
                    ax_all = fig_all.add_subplot(1, 1, 1)
                    ax_all.imshow(ibound_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = 'Set2',
                               extent = (x_axis_plot[0], x_axis_plot[-1], zbot, model.top), vmin = 0., vmax = np.nanmax(ibound_arr))
                    #ax_all.plot(x_axis_plot, topo_lst_plot, linewidth = 2, color = "blue")
                    #ax_all.plot(x_axis_plot2, model.top_elev[: 5 * (model_end_idx - model_start_idx)], linewidth = 1, color = "red")
                    ax_all.plot(x_axis_plot2, model.top_elev[:], linewidth = 2, color = "red")
                    ax_all.plot(- x_nasa, y_nasa, 'ko')
                    ax_all.plot(0.0, - cst_thick_est_avg, 'ko')
                    if fos[0] is not None and cst_offset > 0:
                        ax_all.plot(fos[0][0] - cst_offset, fos[0][1], 'bo')
                    elif fos[0] is not None and cst_offset <= 0:
                        ax_all.plot(fos[0][0] - cst_offset, fos[0][1], 'bo')
                    ax_all.axvline(x = - end_cst_plain, linewidth = 2, color = 'k')    
                    ax_all.axvline(x = 0, linewidth = 2, color = 'k')    
                    ax_all.axhline(y = 0, linewidth = 1, color = 'k') 
                    if cst_offset > 0:
                        ax_all.set_xlim([x_start_new, x_end_new])
                        #ax_all.set_xlim([x_start_new - cst_offset, x_end_new])
                    else:
                        ax_all.set_xlim([x_start_new, x_end_new])
                    ax_all.set_ylim([math.floor(zbot / 100.) * 100, max(model.top_elev) + 50.0])     
                    ax_all.set_title(topo_sc)
                    """
    
    
                   
                    ax_all = fig_all.add_subplot(3, 3, pl_pos)
                    ax_all.imshow(ibound_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = 'Set2',
                               extent = (x_axis_plot[0], x_axis_plot[-1], zbot, model.top), vmin = 0., vmax = np.nanmax(ibound_arr))
                    ax_all.plot(x_axis_plot, topo_lst_plot, linewidth = 2, color = "blue")
                    #ax_all.plot(x_axis_plot2, model.top_elev[: 5 * (model_end_idx - model_start_idx)], linewidth = 1, color = "red")
                    ax_all.plot(x_axis_plot2, model.top_elev[model.col_idx : model.ibound_arr.shape[-1]], linewidth = 1, color = "red")
                    ax_all.plot(- x_nasa, y_nasa, 'ko')
                    ax_all.plot(0.0, - cst_thick_est_avg, 'ko')
                    if fos[0] is not None and cst_offset > 0:
                        ax_all.plot(fos[0][0] - cst_offset, fos[0][1], 'bo')
                    elif fos[0] is not None and cst_offset <= 0:
                        ax_all.plot(fos[0][0] - cst_offset, fos[0][1], 'bo')
                    ax_all.axvline(x = - end_cst_plain, linewidth = 2, color = 'k')    
                    ax_all.axvline(x = 0, linewidth = 2, color = 'k')    
                    ax_all.axhline(y = 0, linewidth = 1, color = 'k') 
                    if cst_offset > 0:
                        ax_all.set_xlim([model.x_start - cst_offset, model.x_end - cst_offset])
                    else:
                        ax_all.set_xlim([model.x_start - cst_offset, model.x_end - cst_offset])
                    ax_all.set_ylim([zbot - 50.0, max(model.top_elev) + 50.0])   
                    ax_all.set_title(topo_sc)
                    pl_pos += 1
                    """
        
                    #   fill in the dictionary                    
                    topo_dict[topo_sc] = {'ibound_arr': model.ibound_arr,\
                                          'top_elev': model.top_elev,\
                                          'bot_elev': model.bot_elev,\
                                          'lay_elev': model.lay_elev,\
                                          'cst_offset': model.cst_offset,\
                                          'lay_idx': model.lay_idx,\
                                          'col_idx': model.col_idx,\
                                          'end_idx': model.end_idx,\
                                          'top': model.top,\
                                          'zbot': zbot,\
                                          'col_width': model.col_width,\
                                          'idx_start': model.idx_start,\
                                          'idx_end': model.idx_end,\
                                          'x_start': model.x_start - cst_offset,\
                                          'x_end': model.x_end - cst_offset,\
                                          #'x_start': model.x_start,\ 
                                          #'x_end': model.x_end,\ 
                                          'soil_thk': model.soil_thk,\
                                          'soil_type': model.soil_type,\
                                          'cst_idx': cst_off_idx,\
                                          'wtd_elev': model.wtd_elev,\
                                          'pcr_rch': model.pcr_rch,\
                                          'watergap_rch': model.watergap_rch,\
                                          'p_min_et': model.p_min_et,\
                                          'k_soil': model.k_soil,\
                                          'drn_rate': model.drn_dens,\
                                          'riv_cond_01': model.riv_cond_01,\
                                          'riv_cond_05': model.riv_cond_05,\
                                          'riv_cond_10': model.riv_cond_10,\
                                          'riv_bot_elev': model.riv_bot_elev,\
                                          'riv_head_elev': model.riv_head_elev,\
                                          'riv_width': model.riv_width,\
                                          'glhymps_top_lay': model.glhymps_top_lay,\
                                          'glhymps_bot_lay': model.glhymps_bot_lay}

                    #   get the counters for the final csv file
                    last_ibound_col_n_layers = len([i for i, x in enumerate(model.ibound_arr[:, 0, -1].tolist()) if x == 1])
                    #   check elevations offshore
                    #offshore_elev = model.top_elev[abs(int(x_start_new)):]
                    offshore_elev = model.top_elev[abs(int((model.x_start) * 10.) - 1):]

                    if last_ibound_col_n_layers < 5:
                        fos_pt_dist = x_end_new #len(offshore_elev) * 10 / 2.
                    else:
                        fos_pt_dist = 0
                    
                    ls_elev_bool = all(i >= -120. for i in offshore_elev)
                    if ls_elev_bool is False:
                        ls_elev = 0
                        ls_lst_below_low_stand = [i for i, x in enumerate(offshore_elev) if x < -120.]
                        if len(offshore_elev) / 2 < len(ls_lst_below_low_stand):
                            ls_elev_50 = 0
                        else:
                            ls_elev_50 = 1
                    else:
                        ls_elev = 1
                        ls_elev_50 = 1
                        

                except IndexError:
                    print("Pica 8")  
                    name = 'no_data'
                    ax_all = fig_all.add_subplot(1, 1, 1)
                    ax_all.text(0.15, 0.5, r'Cant build model, all topography offshore > 0')
                    ax_all.set_title(name)
                    last_ibound_col_n_layers = -5
                    ls_elev = -5
                    ls_elev_50 = -5   
                    fos_pt_dist = -5
                    #return None, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist

    #   plot the overview of all topographical profile
    plot_name_all = os.path.join(out_dir, cst_type + '_all_model_IBOUNDs.png')
    plt.savefig(plot_name_all, dpi = 300, facecolor = 'w', edgecolor = 'w',
        orientation = 'portrait', papertype = None, format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
        frameon = None)
    plt.close()

    #   save the dictionary with the ibound arrays
    dict_save_dir = os.path.join(out_dir, cst_type + '_IBOUND.npy')
    np.save(dict_save_dir, topo_dict)

    return topo_dict, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist

#   function that creates all the folders to store both input and model data for given coscat id
def create_directories(id_coscat, home_dir):
    #   check if the directory already exists
    coscat_dir = os.path.join(home_dir, str(id_coscat))
    if not os.path.exists(coscat_dir):
        os.makedirs(coscat_dir)
        input_file_dir = os.path.join(coscat_dir, '_input_files')
        model_file_dir = os.path.join(coscat_dir, '_model_files') 
        os.makedirs(input_file_dir)
        os.makedirs(model_file_dir)
    else:
        input_file_dir = os.path.join(coscat_dir, '_input_files')
        model_file_dir = os.path.join(coscat_dir, '_model_files')         
        print('Directories for COSCAT id: ' + str(id_coscat) + ' already exists.')
    return input_file_dir, model_file_dir

#   Function that finds the FOS (foot of the contintental slope), based on a paper
def find_FOS(topo_lst, out_dir, cst_type, avg = True):
    
    #topo_lst = lst_topo_henry
    #topo_lst = lst_topo_other
    #topo_lst = final_topo
    #topo_lst = final_topo
    #out_dir = offshore_dir
    #cst_type = topo_sc
    #avg= False
    
    #   clean the lists - if the values are = non_val then change to NaN not to screw up the calculation of averages
    arr_topo = np.array(topo_lst)
    #   distinguish between the cases when all the topographical profiles are supplied vs. when only one topo list is
    if avg:
        avg_topo_lst = np.nanmean(arr_topo, axis = 0)
    else:
        avg_topo_lst = topo_lst
        
    x_axis = np.linspace(0.0, 200.0, 401) 
    
    #   calculate the gradient and second derivative
    grad = np.gradient(avg_topo_lst[400:]) 
    #y_spl = UnivariateSpline(x_axis, avg_topo[400:], s = 0, k = 4)
    #y_spl_2d = y_spl.derivative(n = 2)
    #x_range = np.linspace(x_axis[0], x_axis[-1], 401)
    grad2 = np.gradient(grad) 
    
    loc_max = argrelextrema(grad2, np.greater)    
    loc_min = argrelextrema(grad2, np.less)    
    
    loc_minmax = np.concatenate((loc_max[0], loc_min[0])) # we still need to do this unfortunatly.
    loc_minmax.sort()
    
    loc_max_plt, loc_min_plt, loc_minmax_plt1, loc_minmax_plt2 = [], [], [], []

    for a in range(loc_max[0].shape[0]):
        loc_max_plt.append([x_axis[loc_max[0][a]], avg_topo_lst[400 + loc_max[0][a]]])
    
    for b in range(loc_min[0].shape[0]):
        loc_min_plt.append([x_axis[loc_min[0][b]], avg_topo_lst[400 + loc_min[0][b]]])
        
    for c in range(loc_minmax.shape[0]):
        #loc_minmax_plt1.append([avg_topo[loc_minmax[c]], avg_topo[loc_minmax[c]]])
        loc_minmax_plt1.append(avg_topo_lst[400 + loc_minmax[c]])
        loc_minmax_plt2.append(avg_topo_lst[400 + loc_minmax[c]])
    
    rdp_in1, rdp_in2 = [], []
    for i in range(len(loc_minmax_plt1)):
        rdp_in1.append([loc_minmax_plt1[i], loc_minmax_plt1[i]])
        rdp_in2.append([x_axis[loc_minmax[i]], loc_minmax_plt1[i]])
    rdp_points = rdp.rdp(rdp_in2, epsilon = 2)
    
    for j in range(len(rdp_points)):
        #   if it is the first point in the list start with coastal point - 0.0, 0.0
        if j == 0:
            slope_left = (rdp_points[j][1]) / (rdp_points[j][0])
            slope_right = (rdp_points[j + 1][1] - rdp_points[j][1]) / (rdp_points[j + 1][0] - rdp_points[j][0])
        #   if it is the last point than take the last point as 200.0, avg_topo[-1] 
        elif j == len(rdp_points) - 1:
            slope_left = (rdp_points[j - 1][1] - rdp_points[j][1]) / (rdp_points[j - 1][0] - rdp_points[j][0])
            slope_right = (rdp_points[j][1] - avg_topo_lst[-1]) / (200.0 - rdp_points[j][0])            
        #   otherwise just take the elements before and after the given point
        else:
            slope_left = (rdp_points[j - 1][1] - rdp_points[j][1]) / (rdp_points[j - 1][0] - rdp_points[j][0])
            slope_right = (rdp_points[j + 1][1] - rdp_points[j][1]) / (rdp_points[j + 1][0] - rdp_points[j][0])
            
        rdp_points[j].append(slope_left)
        rdp_points[j].append(slope_right)
        
        #print j, rdp_points[j], round(((slope_left) + (slope_right)) / 2., 2) # slope_left, slope_right

    #   reset the sb_point and fos_point
    sb_point, fos_point = None, None

    #   loop through the points to identify the shelf break (sb_point) sand the foot of the slope (fos_point)
    #   skipping the first point in the list because that one is almost at the coastline..
    for g in range(1, len(rdp_points) - 1):
        #   get all the different values form the list
        pt_dist, pt_depth, pt_slope_left, pt_slope_right = rdp_points[g][0], rdp_points[g][1], rdp_points[g][2], rdp_points[g][3]
        print(g, pt_dist, pt_depth, pt_slope_left, pt_slope_right)
        #   first find the sb_point, it should be above -200m bsl. and the slope on the right from the point
        #   should be higher than the slope on the left from the point
        if pt_depth > -200.0 and abs(pt_slope_left) < abs(pt_slope_right):
            sb_point = [pt_dist, pt_depth]
        #   the fos_point should be deeper than -200m bsl, the absolute slope gradient on the left should be higher
        #   than on the right, and the slope on the right is either positive or less than 25% than the absolute value
        #   of the slope on the left
        elif pt_depth < -1000.0 and abs(pt_slope_left) > abs(pt_slope_right):
            if pt_slope_right > 0 or round(abs(pt_slope_right)) <= round((75 * abs(pt_slope_left)) / 100.):
                #   check if it is just one peak (concave hull) in the profile, do this by checking the slope and 
                #   deptg(elev) of the following part of the topographical profile
                if rdp_points[g + 1][1] > pt_depth and abs(50 * rdp_points[g + 1][3]) / 100. > abs(pt_slope_left):
                    continue
                else:
                    fos_point = [pt_dist, pt_depth]
                    break
            else:
                pass
        
    try:
        #   initialize the figure
        fig = plt.figure(figsize = (20, 20))
        ax1 = fig.add_subplot(1, 1, 1)
        ax2 = ax1.twinx()
        ax1.plot(x_axis, avg_topo_lst[400:], linewidth = 2, color = "blue", label = 'Original topography/bathymetry')
        ax2.plot(x_axis, grad, linewidth = 0.75, color = "red", label = '1st order gradient (slope)')
        ax2.plot(x_axis, grad2, linewidth = 0.75, color = "green", label = '2nd order gradient')
        #   plot the vertical lines demarcating the continental shelf and slope
        try:
            ax1.axvline(x = sb_point[0], linewidth = 1, linestyle = ':', color = 'k', label = 'Location of continental shelf break')  
        except TypeError:
            pass
        try:
            ax1.axvline(x = fos_point[0], linewidth = 1, linestyle = '--', color = 'k',  label = 'Location of foot of continental slope')  
        except TypeError:
            pass
        #ax1.plot([x[0] for x in loc_max_plt], [x[1] for x in loc_max_plt], 'bo')
        #ax1.plot([x[0] for x in loc_min_plt], [x[1] for x in loc_min_plt], 'ro')
        ax1.plot([x[0] for x in rdp_points], [x[1] for x in rdp_points], linewidth = 2., color = "red", label = 'Douglas-Peucker (D-P) simplification')
        ax1.plot([x[0] for x in rdp_points], [x[1] for x in rdp_points], 'ko', label = 'D-P points')
        # ask matplotlib for the plotted objects and their labels
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines + lines2, labels + labels2, loc = 'best') 
        ax1.set_xlabel("Distance from coast (km)")
        ax1.set_ylabel("Depth below sea level (m)")
        ax2.set_ylabel("Gradient (%)")
        #ax2.plot(x_range, y_spl_2d(x_range), linewidth = 2, color = "green")
        #plt.show()
        if out_dir is not None:
            plot_name = os.path.join(out_dir, cst_type + '_FOS_point.png')
            plt.savefig(plot_name, dpi = 300, facecolor = 'w', edgecolor = 'w',
                orientation = 'portrait', papertype = None, format = None,
                transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
                frameon = None)
            plt.close()
            return sb_point, fos_point
        else:
            return sb_point, fos_point, rdp_points
    
    except UnboundLocalError:
        print ('Couldnt find shelf break / foot of continental slope.')
        return [None, None]


#   function that plots all the topographical profiles for given coastal type
def plot_all_topo_profiles(id_coscat, topo_lst, ate_lst, fos_pt, cst_type, out_dir):

    #topo_lst = lst_topo_henry
    #ate_lst = lst_cs_thk_henry
    #coscat_id = id_coscat
    #cst_type = 'henry'
    #out_dir = coscat_dir
    #fos_pt = find_FOS(topo_lst, None, cst_type, avg = True)
    
    plot_title = 'All topographical profiles of type ' + cst_type + ' for COSCAT region ' + str(id_coscat)    
    plot_dir = os.path.join(out_dir, '_' + cst_type + '_all_topo_profiles.png')    
    
    #   calculate the average topo_lst 
    arr_topo = np.array(topo_lst)
    avg_topo = np.nanmean(arr_topo, axis = 0).tolist()        
    avg_ate = sum(ate_lst)/len(ate_lst)    
    
    #   create a figure where all the plots will be stored
    fig = plt.figure(figsize = (15, 7.5))
    ax1 = plt.subplot2grid((1, 2), (0, 0))  #   overall figure
    ax2 = plt.subplot2grid((1, 2), (0, 1))  #   zoomed in coastal zone
    ax1.set_position([0.05, 0.05, 0.65, 0.9]) # [left, bottom, width, height]
    ax2.set_position([0.775, 0.05, 0.2, 0.9])
    
    x_axis = np.linspace(-200., 200., 801)

    ax1.axhline(y = 0.0, color = 'black', linestyle = '--', alpha = 0.5)    
    ax1.axvline(x = 0.0, color = 'black', linestyle = '--', alpha = 0.5)  
    
    ax2.axhline(y = 0.0, color = 'black', linestyle = '--', alpha = 0.5)    
    ax2.axvline(x = 0.0, color = 'black', linestyle = '--', alpha = 0.5)     

    ax1.plot(x_axis, topo_lst[0], color = 'grey', alpha = 0.2, label = 'topographical profile')
    ax1.plot([0.0], [avg_topo[400] - ate_lst[0]], marker = 'o', markersize = 2, color = 'grey', alpha = 0.5)
    ax2.plot(x_axis, topo_lst[0], color = 'grey', alpha = 0.2)
    ax2.plot([0.0], [avg_topo[400] - ate_lst[0]], marker = 'o', markersize = 6, color = 'grey', alpha = 0.5, label = 'ATE point', linestyle = 'None')  
    
    for x in range(1, len(topo_lst)):
        ax1.plot(x_axis, topo_lst[x], color = 'grey', alpha = 0.2)
        ax1.plot([0.0], [avg_topo[400] - ate_lst[x]], marker = 'o', markersize = 2, color = 'grey', alpha = 0.5)
        ax2.plot(x_axis, topo_lst[x], color = 'grey', alpha = 0.2)
        ax2.plot([0.0], [avg_topo[400] - ate_lst[x]], marker = 'o', markersize = 6, color = 'grey', alpha = 0.5)            

    ax1.plot(x_axis, avg_topo, color = 'red', linewidth = 2., label = 'Average topographical profile')
    ax1.plot([0.0], [avg_topo[400] - avg_ate], marker = 'o', markersize = 4, color = 'red', alpha = 0.75)
    ax2.plot(x_axis, avg_topo, color = 'red', linewidth = 2.)
    ax2.plot([0.0], [avg_topo[400] - avg_ate], marker = 'o', markersize = 12, color = 'red', alpha = 0.75, label = 'Average ATE point', linestyle = 'None')
    
    ax1.set_xlabel("Distance from coast (km)")
    ax1.set_ylabel("Elevation above sea level (m)")
    ax1.set_title(plot_title, fontsize = 18)

    # Create a Rectangle patch
    top_y = math.ceil(avg_topo[400] / 10.) * 100
    bot_y = top_y + abs(math.floor((avg_topo[400] - max(ate_lst)) / 100.) * 100)
    rect = patches.Rectangle((-2.5, top_y), 5., -bot_y, linewidth = 1, edgecolor = 'blue', facecolor = 'none', alpha = 0.4)
    ax1.add_patch(rect)
    
    ax2.set_xlim([-2.5, 2.5])
    ax2.set_ylim([-bot_y, top_y])

    if fos_pt[0] is not None:
        ax1.plot([fos_pt[0][0]], [fos_pt[0][1]], marker = 'o', markersize = 6, color = 'blue', alpha = 0.75, label = 'Shelf edge point', linestyle = 'None')

    if fos_pt[1] is not None:
        ax1.plot([fos_pt[1][0]], [fos_pt[1][1]], marker = 'o', markersize = 6, color = 'green', alpha = 0.75, label = 'Foot of cont. slope point', linestyle = 'None')

    
    rcParams['font.family'] = 'Garamond'
    rcParams['axes.facecolor'] = 'white'
    rcParams['savefig.facecolor'] = 'white'  
    rcParams['legend.numpoints'] = 1

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    
    ax1.legend(lines + lines2, labels + labels2, loc = 'best') 

    fig.savefig(plot_dir, dpi = 300, facecolor = 'w', edgecolor = 'w',
                orientation = 'portrait', papertype = None, format = None,
                transparent = False, bbox_inches = 'tight', pad_inches = 0.3,
                frameon = None)
    plt.close()




    """
    from scipy.spatial import ConvexHull
    hull = ConvexHull(rdp_points)

    plt.plot(rdp_points[hull.vertices][0], rdp_points[hull.vertices][1], 'r--', lw=2)
    plt.plot(rdp_points[hull.vertices[0],0], rdp_points[hull.vertices[0],1], 'ro')
    plt.show()
    """


def get_normal_dist(in_lst, plot_dir, x_axis_lbl = None, nan_val = None, no_negative_vals = True, plot = True):
  
    #in_lst = lst_cs_thk_delta
    #nan_val = None
    #plot_dir =  os.path.join(geometry_pic_dir, '_delta_cs_thk.png')
    #x_axis_lbl = 'Log normal distribution of coastal thickness'
    #   transform the list into an array and mark the Nan values if required

    mu_out, std_out, arr = -1., -1., -1.

    arr = np.array(in_lst)
    #   make sure that the type is float
    if arr.dtype != 'float32':   
        arr = arr.astype(np.float)
    if nan_val is not None:
        arr[arr == nan_val] = np.nan
        #   if required remove the negative values
        if no_negative_vals is True:        
            arr[arr < 0] = np.nan
        #   next, remove the Nan values
        arr = arr[np.logical_not(np.isnan(arr))]
        
    #   then flatten the array to 1D
    arr = arr.flatten()
    arr_ln = np.log(arr)
    # Fit a normal distribution to the data and round the mean and std values
    mu, std = norm.fit(arr_ln)    
    mu_out, std_out = round(mu, 4), round(std, 4)
    
    if plot:    
        try:
            # Plot the histogram.
            plt.hist(arr_ln, normed = True, bins = 'auto', color = 'g')
            # Plot the PDF.
            xmin, xmax = plt.xlim()
            x = np.linspace(xmin, xmax, 100)
            p = norm.pdf(x, mu_out, std_out)
            plt.plot(x, p, 'k', linewidth = 2, label = 'normal distribution PDF')
            plt.legend(loc = 0)        
            plt.xlabel("Natural log (ln) of " + x_axis_lbl)     
            #plt.text(-11, 0.4, "Fit results: mu = %s,  std = %s" % (mu_out, std_out))        
            title = "Fit results: mu = %s,  std = %s" % (mu_out, std_out)
            plt.title(title)
            rcParams['font.family'] = 'Garamond'
            rcParams['axes.facecolor'] = 'white'
            rcParams['savefig.facecolor'] = 'white'      
            plt.savefig(plot_dir)        
            plt.close()
        except ValueError:
            print('couldnt print figure')
            pass
    
    return mu_out, std_out, arr





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
coscat_ids = dbase_connect_world.get_ids_coscat("coscat_id", "cs_coastal_types")
#coscat_ids = dbase_connect.get_ids_coscat("coscat_id", "cs_coastal_types")


csv_dir_ts = os.path.join(model_dir, '_COSCAT_representative_models_summary_ATE_MEAN_STDEV.csv')  

headers = ['id_COSCAT', 'henry_cs_thk_mu', 'henry_cs_thk_std', 'other_cs_thk_mu', 'other_cs_thk_std', 'delta_cs_thk_mu', 'delta_cs_thk_std',\
'delta_other_cs_thk_mu', 'delta_other_cs_thk_std', 'henry_other_cs_thk_mu', 'henry_other_cs_thk_std', 'all_cs_thk_mu', 'all_cs_thk_std',\
'henry_cs_width_mu', 'henry_cs_width_std', 'other_cs_width_mu', 'other_cs_width_std', 'delta_cs_width_mu', 'delta_cs_width_std',\
'delta_other_cs_width_mu', 'delta_other_cs_width_std', 'henry_other_cs_width_mu', 'henry_other_cs_width_std', 'all_cs_width_mu', 'all_cs_width_std']
         
"""
f = open(csv_dir_ts,'w')
f.write('id_COSCAT, henry_cs_thk_mu, henry_cs_thk_std, other_cs_thk_mu, other_cs_thk_std, delta_cs_thk_mu, delta_cs_thk_std,\
         delta_other_cs_thk_mu, delta_other_cs_thk_std, henry_other_cs_thk_mu, henry_other_cs_thk_std, all_cs_thk_mu, all_cs_thk_std,\
         henry_cs_width_mu, henry_cs_width_std, other_cs_width_mu, other_cs_width_std,delta_cs_width_mu, delta_cs_width_std,\
         delta_other_cs_width_mu, delta_other_cs_width_std, henry_other_cs_width_mu, henry_other_cs_width_std, all_cs_width_mu, all_cs_width_std')
         #henry_anchor_dist_mu, henry_anchor_dist_std, other_anchor_dist_mu, other_anchor_dist_std,delta_anchor_dist_mu, delta_anchor_dist_std,\
         #delta_other_anchor_dist_mu, delta_other_anchor_dist_std, all_anchor_dist_mu, all_anchor_dist_std,\
         #henry_anchor_depth_mu, henry_anchor_depth_std, other_anchor_depth_mu, other_anchor_depth_std,delta_anchor_depth_mu, delta_anchor_depth_std,\
         #delta_other_anchor_depth_mu, delta_other_anchor_depth_std, all_anchor_depth_mu, all_anchor_depth_std') #write the header
f.close()
"""

rows_to_csv = []

#   loop through all the COSCAT ids
for coscat_id in coscat_ids:
    
    id_num_coscat = coscat_id[0]
    #   select all the cross-sections that are within the given COSCAT region
    cond_str = "coscat_id = " + str(id_num_coscat)
    
    coscat_ids = dbase_connect_world.get_data_condition("id_cs", cond_str, "cs_coastal_types")
    #coscat_ids = dbase_connect.get_data_condition("id_cs", cond_str, "cs_coastal_types")
    
    #if os.path.isfile(os.path.join(plot_dir, str(id_num_coscat) + ".png")):    
    
    #   add 0s in front of the name to have an ordered list of COSCAT ids in the folder 
    if id_num_coscat < 10:
        str_id_num_coscat = '000' + str(id_num_coscat)
    elif id_num_coscat >= 10 and id_num_coscat < 100:
        str_id_num_coscat = '00' + str(id_num_coscat)
    elif id_num_coscat >= 100 and id_num_coscat < 1000:
        str_id_num_coscat = '0' + str(id_num_coscat)    
    else:
        str_id_num_coscat = str(id_num_coscat)   
    
    #   define the model directory where the model files and figures will be stored
    coscat_dir = os.path.join(model_dir, str_id_num_coscat)
    rch_pic_dir = os.path.join(coscat_dir, '_RCH_hist')
    geometry_pic_dir = os.path.join(coscat_dir, '_GEOMETRY_hist')
    geology_pic_dir = os.path.join(coscat_dir, '_GEOLOGY_hist')    
    
    if not os.path.exists(coscat_dir):
        os.makedirs(coscat_dir)
        os.makedirs(rch_pic_dir)
        os.makedirs(geometry_pic_dir)
        os.makedirs(geology_pic_dir)
    else:
        print('Directories for COSCAT id: ' + str_id_num_coscat + ' already exists.')

    """
    #   create a map of the COSCAT region for better understanding 
    lst_in = [int(i[0]) for i in coscat_ids]
    test = get_geoplot_extent(lst_in, db_name_world, db_host, db_user, db_pass, table_db)
    test_coscat = get_coscat_extent(test[6], db_name_world, db_host, db_user, db_pass, table_coscat_polys)     
      
    #   based on the x_min, x_max, y_min and y_max from the COSCAT and coastal points extent select the larger block area
    x_min = min(test[0], test_coscat[1])
    x_max = max(test[1], test_coscat[0])
    y_min = min(test[2], test_coscat[2])
    y_max = max(test[3], test_coscat[3])  
    x_pt_lst, y_pt_lst, cst_pt_type_lst = test[4], test[5], test[7]        
    
    # clip the shapefile with the raster bounds 
    clipped = input.filter(bbox=((x_min , y_min, x_max, y_max)))
    # create the clipped shapefile with the same schema
    clipped_schema = input.schema.copy()
    clipped_name = 'g:/_CREATED_DATA/_A2_data/glim_su_' + str(id_num_coscat) + '.shp'
    #clipped_name = 'g:/_CREATED_DATA/_A2_data/glim_su_' + str(808) + '.shp'
    with fiona.collection(clipped_name, 'w', 'ESRI Shapefile', clipped_schema) as output:
        for elem in clipped:
               output.write({'properties': elem['properties'],'geometry': mapping(shape(elem['geometry']))})        
    
    map_dir = os.path.join(plot_dir, str(id_num_coscat) + "_map.png")
    map_plot = plot_map(int(id_num_coscat), map_dir)      
    
    gl_map_dir = os.path.join(plot_dir, str(id_num_coscat) + "_global_map.png")
    gl_map_plot = coscat_global_map(int(id_num_coscat), gl_map_dir)
    """
    
    #   define the model dimensions and empty lists to be filled 
    del_col = 100.
    del_lay = 10.
    cs_points_dist = 500.
    
    lst_cs_thk_henry_other, lst_cs_thk_delta, lst_cs_thk_henry, lst_cs_thk_other, lst_cs_thk_all, lst_cs_thk_delta_other = [], [], [], [], [], []
    lst_cs_width_henry_other, lst_cs_width_delta, lst_cs_width_henry, lst_cs_width_other, lst_cs_width_all, lst_cs_width_delta_other = [], [], [], [], [], []
    
    #   set the counters of different coastal types
    henry_other_n, henry_n, other_n, delta_n, delta_other_n, all_n = 0, 0, 0, 0, 0, 0        
    
    #   loop through the cross-sections and select required info from all the necessary tables
    for id_num in coscat_ids:
        id_cs = int(id_num[0])
        #   info about coastal plain width and the anchor point

        #cs_info = dbase_connect_world.get_data_condition("cst_plain_width, nasa_point_dist, nasa_point_depth, overall_avg",\
        #                                                 "id_cs = " + str(id_cs), "cs_sed_thick_est") 
        cs_info = dbase_connect.get_data_condition("cst_plain_width, nasa_point_dist, nasa_point_depth, overall_avg",\
                                                   "id_cs = " + str(id_cs), "cs_sed_thick_est_v2")
        #   also get the classification/delta of the coastal point
        cond_str_id_cs = "id_cs = '" + str(id_cs) + "'"
        
        #cs_class = dbase_connect_world.get_data_condition("class_fin, deltaid", cond_str_id_cs, "cs_coastal_types_all")
        cs_class = dbase_connect_world.get_data_condition("class_fin", cond_str_id_cs, "cs_coastal_types")
        #cs_class = dbase_connect.get_data_condition("class_fin", cond_str_id_cs, "cs_coastal_types")
        
        #delta = int(cs_class[0][1])
        coastal_class = cs_class[0][0]

        try:
            cs_plain_width, anchor_dist, anchor_depth, avg_depth = float(cs_info[0][0]), 200.0 - float(cs_info[0][1]), float(cs_info[0][2]), float(cs_info[0][3]) 
        except IndexError:
            continue
        
        #   if it is delta then and the anchor depth is more than 200km then assign 200km
        if coastal_class == 1:
            if anchor_dist < 0.:
                anchor_dist = 200.0 + anchor_dist

        cs_model_point = ws_ms.point(None, None, 400, 200000, id_cs = id_cs)
        
        #   happens when there are no data extracted for the id_cs
        try:
            cs_model_point.get_input_data(adjust_to_ocean = True)    

            lst_cs_thk_all.append(avg_depth)
            lst_cs_width_all.append(cs_plain_width)
    
            all_n += 1             
            
            #   first check that the type is different from delta (based on Sywitskis shapefile)
            if coastal_class == 1:
                print('delta')
                lst_cs_thk_delta.append(avg_depth)
                lst_cs_width_delta.append(cs_plain_width)
    
                lst_cs_thk_delta_other.append(avg_depth)
                lst_cs_width_delta_other.append(cs_plain_width)
                
                delta_n += 1
                delta_other_n += 1
                
            elif coastal_class == 2:
                print('henry')
                lst_cs_thk_henry.append(avg_depth)
                lst_cs_width_henry.append(cs_plain_width)
                
                lst_cs_thk_henry_other.append(avg_depth)
                lst_cs_width_henry_other.append(cs_plain_width)
    
                henry_n += 1    
                henry_other_n += 1                 
                
            elif coastal_class == 3:
                print('other')
                lst_cs_thk_other.append(avg_depth)
                lst_cs_width_other.append(cs_plain_width)              
    
                lst_cs_thk_delta_other.append(avg_depth)
                lst_cs_width_delta_other.append(cs_plain_width)  
    
                lst_cs_thk_henry_other.append(avg_depth)
                lst_cs_width_henry_other.append(cs_plain_width)
    
                other_n += 1
                delta_other_n += 1
                henry_other_n += 1
                
            else:
                print('other')
                lst_cs_thk_other.append(avg_depth)
                lst_cs_width_other.append(cs_plain_width)  
                
                lst_cs_thk_delta_other.append(avg_depth)
                lst_cs_width_delta_other.append(cs_plain_width)
    
                lst_cs_thk_henry_other.append(avg_depth)
                lst_cs_width_henry_other.append(cs_plain_width)
    
                other_n += 1
                delta_other_n += 1    
                henry_other_n += 1
            
            
        except IndexError:

            if coastal_class == 1:
                all_n += 1    
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
                
                #   first check that the type is different from delta (based on Sywitskis shapefile)
                if coastal_class == 1:
                    print('delta')
                    lst_cs_thk_delta.append(avg_depth)
                    lst_cs_width_delta.append(cs_plain_width)
   
                    lst_cs_thk_delta_other.append(avg_depth)
                    lst_cs_width_delta_other.append(cs_plain_width)
                    
                    delta_n += 1
                    delta_other_n += 1

        

    
    mu_delta_cs_thk, std_delta_cs_thk = -1, -1
    mu_delta_cs_width, std_delta_cs_width = -1, -1

    mu_henry_cs_thk, std_henry_cs_thk = -1, -1
    mu_henry_cs_width, std_henry_cs_width = -1, -1
    
    mu_other_cs_thk, std_other_cs_thk = -1, -1
    mu_other_cs_width, std_other_cs_width = -1, -1
    
    mu_henry_other_cs_thk, std_henry_other_cs_thk = -1, -1
    mu_henry_other_cs_width, std_henry_other_cs_width = -1, -1

    mu_delta_other_cs_thk, std_delta_other_cs_thk = -1, -1
    mu_delta_other_cs_width, std_delta_other_cs_width = -1, -1

    mu_all_cs_thk, std_all_cs_thk = -1, -1
    mu_all_cs_width, std_all_cs_width = -1, -1
    
    #   create plots for each coastal type 
    if lst_cs_thk_delta == []:
        delta_ls_elev = -1
        delta_ls_elev_50 = -1
        pass
    else:
        try:
            mu_delta_cs_thk, std_delta_cs_thk = norm.fit(lst_cs_thk_delta)
            mu_delta_cs_width, std_delta_cs_width = norm.fit(lst_cs_width_delta)
        except TypeError:
            pass

    #   create plots for each coastal type 
    if lst_cs_thk_henry_other == []:
        henry_other_ls_elev = -1
        henry_other_ls_elev_50 = -1
        pass
    else:
        try:
            mu_henry_other_cs_thk, std_henry_other_cs_thk = norm.fit(lst_cs_thk_henry_other)
            mu_henry_other_cs_width, std_henry_other_cs_width = norm.fit(lst_cs_width_henry_other)            
        except TypeError:
            pass
    
    if lst_cs_thk_henry == []:   
        henry_ls_elev = -1
        henry_ls_elev_50 = -1            
        pass
    else:
        try:         
            mu_henry_cs_thk, std_henry_cs_thk = norm.fit(lst_cs_thk_henry)
            mu_henry_cs_width, std_henry_cs_width = norm.fit(lst_cs_width_henry)         
        except TypeError:
            pass

    if lst_cs_thk_other == []:   
        other_ls_elev = -1
        other_ls_elev_50 = -1            
        pass
    else:          
        try:
            mu_other_cs_thk, std_other_cs_thk = norm.fit(lst_cs_thk_other)     
            mu_other_cs_width, std_other_cs_width = norm.fit(lst_cs_width_other)         
        except (TypeError, ZeroDivisionError):
            pass
 
    if lst_cs_thk_all == []:   
        all_ls_elev = -1
        all_ls_elev_50 = -1            
        pass
    else:          
        try:
            mu_all_cs_thk, std_all_cs_thk = norm.fit(lst_cs_thk_all)     
            mu_all_cs_width, std_all_cs_width = norm.fit(lst_cs_width_all)                     
        except TypeError:
            pass     
 
    if lst_cs_thk_delta_other == []:   
        delta_other_ls_elev = -1
        delta_other_ls_elev_50 = -1            
        pass
    else:          
        try:
            mu_delta_other_cs_thk, std_delta_other_cs_thk = norm.fit(lst_cs_thk_delta_other)     
            mu_delta_other_cs_width_mu, std_delta_other_cs_width_mu = norm.fit(lst_cs_width_delta_other)                         
        except TypeError:
            pass          

    rows_to_csv.append([int(id_num_coscat), round(mu_henry_cs_thk, 1), round(std_henry_cs_thk, 1), round(mu_other_cs_thk, 1), round(std_other_cs_thk, 1),\
                        round(mu_delta_cs_thk, 1), round(std_delta_cs_thk, 1), round(mu_delta_other_cs_thk, 1), round(std_delta_other_cs_thk, 1),\
                        round(mu_henry_other_cs_thk, 1), round(std_henry_other_cs_thk, 1), round(mu_all_cs_thk, 1), round(std_all_cs_thk, 1),\
                        round(mu_henry_cs_width, 1), round(std_henry_cs_width, 1), round(mu_other_cs_width, 1), round(std_other_cs_width, 1),\
                        round(mu_delta_cs_width, 1), round(std_delta_cs_width, 1), round(mu_delta_other_cs_width, 1), round(std_delta_other_cs_width, 1),\
                        round(mu_henry_other_cs_width, 1), round(std_henry_other_cs_width, 1), round(mu_all_cs_width, 1), round(std_all_cs_width, 1)])
    
import pandas as pd
    
df = pd.DataFrame(rows_to_csv, columns = headers)   
df.to_csv(csv_dir_ts)  

"""
f = open(csv_dir_ts, 'a')
f.write("\n")
f.write(str(id_num_coscat) + ',' + str(mu_henry_cs_thk) + ',' + str(std_henry_cs_thk) + ',' + str(mu_other_cs_thk) + ',' + str(std_other_cs_thk)\
        + ',' + str(mu_delta_cs_thk) + ',' + str(std_delta_cs_thk) + ',' + str(mu_delta_other_cs_thk) + ',' + str(std_delta_other_cs_thk)\
        + ',' + str(mu_henry_other_cs_thk) + ',' + str(std_henry_other_cs_thk) + ',' + str(mu_all_cs_thk) +  ',' + str(std_all_cs_thk)\
        + ',' + str(mu_henry_cs_width) + ',' + str(std_henry_cs_width) + ',' + str(mu_other_cs_width) + ',' + str(std_other_cs_width)\
        + ',' + str(mu_delta_cs_width) + ',' + str(std_delta_cs_width) + ',' + str(mu_delta_other_cs_width) + ',' + str(std_delta_other_cs_width)\
        + ',' + str(mu_henry_other_cs_width) + ',' + str(std_henry_other_cs_width) + ',' + str(mu_all_cs_width) +  ',' + str(std_all_cs_width))
f.close()    
"""










