# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 15:25:34 2019

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
import xarray as xr
import psycopg2
#import gdal
import utm
from time import time
import pandas as pd
import geopandas as gpd
import random 
from shapely.geometry import Point, MultiPoint
from shapely.ops import nearest_points
import matplotlib 
#matplotlib.use('agg')
import matplotlib.pyplot as plt
#plt.ioff()

from osgeo import gdal

sys_dir = r'g:\Water_Nexus\_A4_GUM\_scripts\_fc_scripts'
input_files_dir = r'g:\Water_Nexus\_A4_GUM\_input_data_COSCAT_REG'
plot_dir = r'g:\Water_Nexus\_A4_GUM\_MODEL_input_files\_plots'
model_dir = r'g:\Water_Nexus\_A4_GUM\_MODEL_input_files'

"""
sys_dir = r'home/dzamrsky/_A4/_scripts/_fc_scripts'
#input_files_dir = r'projects/0/qt16165/_dzamrsky/_A4_input_data'
plot_dir = r'projects/0/qt16165/_dzamrsky/_A4_MODEL_input_files/_plots'
model_dir = r'projects/0/qt16165/_dzamrsky/_A4_MODEL_input_files'

input_files_dir = r'/gpfs/scratch1/nodespecific/int1/dzamrsky/_A4_input_data'
"""

import sys
sys.path.append(sys_dir)
import ws_masterscript_py3 as ws_ms
import ws_datatools_py3 as ws               # connect to dbase, run sqls etc.
from modeltools_py3 import cs_model 
#import coastal_db_GIS_py3 as gis_tools



ne_coastline = os.path.join(input_files_dir, 'natural_earth', 'ne_10m_land')
coscat_polys = os.path.join(input_files_dir, 'coscat_dissolved')  
input = fiona.open(os.path.join(input_files_dir, 'GLIM_su_only.shp'), 'r')


#ne_coastline = r'g:\_ORIGINAL_DATA\natural_earth\ne_10m_land'
#coscat_polys = r'g:\_CREATED_DATA\_A2_data\coscat_analysis\coscat_dissolved'
#f = open('g:\Water_Nexus\_A2\_figures\_representative_profiles_COSCATs\errors.txt','w')
# Read the original GLIM_su Shapefile
#input = fiona.open('g:/_CREATED_DATA/_A2_data/GLIM_su_only.shp', 'r')


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
def save_model_input_files_cst_type_topo(cst_type, thk_lst, width_lst, anchor_dist_lst, anchor_depth_lst, wtd_lst, topo_lst, topo_gebco_500m_lst,\
                                         soil_type_lst, soil_thk_lst, nasa_lst, offshore_lst, pcr_lst, watergap_lst, p_min_et_lst,\
                                         k_soil_lst, drn_rate_lst,\
                                         glhymps_top_lay_lst, glhymps_bot_lay_lst, topo_res, del_col, del_lay, id_coscat, out_dir, out_name):

    """
    thk_lst = lst_cs_thk
    width_lst = lst_cs_width
    anchor_dist_lst = lst_anchor_dist
    anchor_depth_lst = lst_anchor_depth
    wtd_lst = lst_wtd
    topo_lst = lst_gebco_merit_avg_100m
    topo_gebco_500m_lst = lst_gebco_merit_avg_500m
    soil_type_lst = soil_type_lst
    soil_thk_lst = lst_soil_thk
    nasa_lst = lst_cs_nasa
    offshore_lst = lst_cs_offshore 
    pcr_lst = lst_cs_pcr_rch
    watergap_lst = lst_cs_watergap_rch
    p_min_et_lst = lst_cs_p_min_et
    k_soil_lst = lst_cs_k_soil
    drn_rate_lst = lst_cs_drn_rate
    glhymps_top_lay_lst = lst_cs_glhymps_top_lay
    glhymps_bot_lay_lst = lst_cs_glhymps_bot_lay
    id_coscat = id_cs
    out_dir = id_riv_bas_dir      
    topo_diff = 10.
    bot_limit = -10000.
    ibound_act_cells_limit = 2
    offshore_dist_limit = 200000.
    topo_res = 100.
    del_lay = 10.
    del_col = 100.
    out_name = '_gebco_merit_avg_100m'
 
    
    new_ibound_dict = rivbas.save_model_input_files_cst_type_topo('', thk_lst, width_lst, anchor_dist_lst, anchor_depth_lst, wtd_lst, topo_lst, topo_gebco_500m_lst,\
                                                              soil_type_lst, soil_thk_lst, nasa_lst, offshore_lst, pcr_lst, watergap_lst, p_min_et_lst, k_soil_lst, drn_rate_lst,\
                                                              glhymps_top_lay_lst, glhymps_bot_lay_lst, 100., 100., 10., coscat_id, subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND')

    
    
    
    id_coscat = coscat_id
    out_dir = subreg_input_dir
    out_name = '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND'
    
    for i in range(len(lst_gebco_merit_avg_100m[0])):
        if math.isnan(lst_gebco_merit_avg_100m[0][i]):
            print(i)
    
    for i in range(len(lst_gebco_merit_avg_500m[0])):
        if math.isnan(lst_gebco_merit_avg_500m[0][i]):
            print(i)
        
    
    
    """
    sea_level = 0.0

    last_ibound_col_n_layers = -55
    ls_elev = -55
    ls_elev_50 = -55   
    fos_pt_dist = -55
    
    #   clean the lists - if the values are = non_val then change to NaN not to screw up the calculation of averages
    arr_wtd = np.array(wtd_lst)
    arr_wtd[arr_wtd <= -999.0] = np.nan     # sometimes there are values of -9999.0 might be an error in the dataset, otherwise nonval = -999
    arr_topo = np.array(topo_lst)
    arr_topo_gebco_500m = np.array(topo_gebco_500m_lst)
    arr_topo[arr_topo == -9999.0] = np.nan
    arr_topo_gebco_500m[arr_topo_gebco_500m == -9999.0] = np.nan
    arr_soil_type = np.array(soil_type_lst)
    arr_soil_thk = np.array(soil_thk_lst)
    arr_soil_thk[arr_soil_thk <= -9999.0] = np.nan  # -32768.0
    arr_nasa = np.array(nasa_lst) 
    arr_nasa[arr_nasa == -1.0] = np.nan
    arr_offshore = np.array(offshore_lst) # -32768
    arr_offshore[arr_offshore == -32768.0] = np.nan
    
    arr_k_soil = np.array(k_soil_lst)
    arr_k_soil[arr_k_soil <= -9999.] = np.nan
    arr_drn_rate = np.array(drn_rate_lst)   
    arr_drn_rate[arr_drn_rate <= -9999.] = np.nan
    
    arr_glhymps_top_lay = np.array(glhymps_top_lay_lst)   
    arr_glhymps_bot_lay = np.array(glhymps_bot_lay_lst)   
    arr_glhymps_top_lay[arr_glhymps_top_lay <= -9999.] = np.nan   
    arr_glhymps_bot_lay[arr_glhymps_bot_lay <= -9999.] = np.nan   
    
    avg_wtd = np.nanmean(arr_wtd, axis = 0)
    avg_topo = np.nanmean(arr_topo, axis = 0)
    avg_topo_gebco_500m = np.nanmean(arr_topo_gebco_500m, axis = 0)    
    avg_nasa = np.nanmean(arr_nasa, axis = 0)
    avg_offshore = np.nanmean(arr_offshore, axis = 0)

    try:
        avg_k_soil = np.nanmean(arr_k_soil, axis = 0)
    except ZeroDivisionError:
        arr_k_soil_data = np.array(arr_k_soil, dtype='float')
        avg_k_soil = np.nanmean(arr_k_soil_data, axis = 0)
     
    avg_glhymps_top_lay = np.nanmean(arr_glhymps_top_lay, axis = 0)        
    avg_glhymps_bot_lay = np.nanmean(arr_glhymps_bot_lay, axis = 0)        
        
    avg_drn_rate = np.nanmean(arr_drn_rate, axis = 0)
    avg_soil_thk = np.nanmean(arr_soil_thk, axis = 0)
    avg_soil_type = []
    
    for i in range(avg_topo_gebco_500m.shape[0]):
        if avg_topo_gebco_500m[i] >= 0:
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
            lst_wtd_topo.append(avg_topo_gebco_500m[a] - wtd_vals_round[a])
        else:
            lst_wtd_topo.append(avg_topo_gebco_500m[a])   
    arr_wtd_topo = np.array(lst_wtd_topo)        
    #avg_wtd_topo = np.nanmean(arr_wtd_topo, axis = 0)
    
    #   refine the lists based on the resolution of the topography list 
    x_lst = np.arange(-200., 200 + 500. / 1000., 500. / 1000.)  # 500. because that is the default spacing of cross-section points
    x_refined_lst = np.arange(-200., 200 + topo_res / 1000., topo_res / 1000.)    
    
    avg_wtd_topo_refined = np.interp(x_refined_lst, x_lst, arr_wtd_topo).tolist()    
    avg_wtd_topo_refined = [round(i, 4) for i in avg_wtd_topo_refined]
    avg_nasa_refined = np.interp(x_refined_lst, x_lst, avg_nasa).tolist()    
    avg_nasa_refined = [round(i, 4) for i in avg_nasa_refined]    
    avg_offshore_refined = np.interp(x_refined_lst, x_lst, avg_offshore).tolist()    
    avg_offshore_refined = [round(i, 4) for i in avg_offshore_refined]   
    #avg_watergap_rch_refined = np.interp(x_refined_lst, x_lst, avg_watergap_rch).tolist()    
    #avg_watergap_rch_refined = [round(i, 1) for i in avg_watergap_rch_refined]  
    #avg_pcr_rch_refined = np.interp(x_refined_lst, x_lst, avg_pcr_rch).tolist()    
    #avg_pcr_rch_refined = [round(i, 1) for i in avg_pcr_rch_refined]      
    #avg_p_min_et_refined = np.interp(x_refined_lst, x_lst, avg_p_min_et).tolist()    
    #avg_p_min_et_refined = [round(i, 1) for i in avg_p_min_et_refined]   
    avg_k_soil_refined = np.interp(x_refined_lst, x_lst, avg_k_soil).tolist()    
    avg_k_soil_refined = [round(i, 7) for i in avg_k_soil_refined]   
    avg_drn_rate_refined = np.interp(x_refined_lst, x_lst, avg_drn_rate).tolist()    
    avg_drn_rate_refined = [round(i, 4) for i in avg_drn_rate_refined]   
    avg_soil_thk_refined = np.interp(x_refined_lst, x_lst, avg_soil_thk).tolist()    
    avg_soil_thk_refined = [round(i, 4) for i in avg_soil_thk_refined]   
    avg_soil_type_refined = np.interp(x_refined_lst, x_lst, avg_soil_type).tolist()    
    avg_soil_type_refined = [round(i, 4) for i in avg_soil_type_refined]   
    avg_glhymps_top_lay_refined = np.interp(x_refined_lst, x_lst, avg_glhymps_top_lay).tolist()  
    avg_glhymps_top_lay_refined = [round(i, 7) for i in avg_glhymps_top_lay_refined]
    avg_glhymps_bot_lay_refined = np.interp(x_refined_lst, x_lst, avg_glhymps_bot_lay).tolist()  
    avg_glhymps_bot_lay_refined = [round(i, 7) for i in avg_glhymps_bot_lay_refined]


    end_cst_plain = round(sum(width_lst) / float(len(width_lst)), 1)
    x_nasa = round(sum(anchor_dist_lst) / float(len(anchor_dist_lst)), 1)
    #y_nasa = round(sum(anchor_depth_lst) / float(len(anchor_depth_lst)), 1)
    cst_thick_est_avg = round(sum(thk_lst) / float(len(thk_lst)), 1)
        
    #   cst_thick_est_stdev = round(cs_model_point.sed_thick_est_stdev, 0) 
    topo_dict = dict()
    topo_dict['ibound_arr'] = []
    topo_dict['top_elev'] = []
    topo_dict['bot_elev'] = []
    topo_dict['lay_elev'] = []
    topo_dict['cst_offset'] = []
    topo_dict['lay_idx'] = []
    topo_dict['col_idx'] = []
    topo_dict['end_idx'] = []
    topo_dict['top'] = []
    topo_dict['zbot'] = []
    topo_dict['col_width'] = []
    topo_dict['idx_start'] = []
    topo_dict['idx_end'] = []
    topo_dict['x_start'] = []
    topo_dict['x_end'] = []
    topo_dict['soil_thk'] = []
    topo_dict['soil_type'] = []
    topo_dict['cst_idx'] = []
    topo_dict['wtd_elev'] = []
    topo_dict['k_soil'] = []
    topo_dict['drn_rate'] = []     
    topo_dict['glhymps_top_lay'] = []
    topo_dict['glhymps_bot_lay'] = []
            
    in_lst = avg_topo.tolist()
    off_lst = avg_topo.tolist()
    final_topo = avg_topo.tolist()
        
    #   first check the topographical profiles both offshore and inland - if offshore all > 0 then dont model
    #   same case is for inland all < 0
    off_pos = [i for i in off_lst if i < 0]
    in_pos = [i for i in in_lst if i > 0]    

    #   create an overall figure where all the topographical models will be stored
    fig_all = plt.figure(figsize = (20, 10))

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
        #cst_look_up_idx = 200
        #cst_look_up_idx =  int((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2)
        cst_look_up_idx =  int((((len(avg_topo) - 1) / 2) / 2))
        
        #   find the position of the coastline
        cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < sea_level), None)
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
            #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
            cst_offset = ((cst_look_up_idx + cst_offset_idx) - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000.
            cst_off_idx = cst_look_up_idx + cst_offset_idx
    
            if cst_off_idx < cst_look_up_idx:
                print("Pica 2")                 
                #   define the index from which we will look for the coastline 
                cst_look_up_idx = int((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) + int(((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) / 2)
                #   find the position of the coastline
                cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < sea_level), None)
                cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                cst_offset = ((cst_look_up_idx + cst_offset_idx) - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000.
                cst_off_idx = cst_look_up_idx + cst_offset_idx
                
                if cst_off_idx < cst_look_up_idx:
                    print("Pica 3")                         
                    
                    #   define the index from which we will look for the coastline 
                    cst_look_up_idx = int((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) + int(((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) / 4) * 3
                    #   find the position of the coastline
                    cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < sea_level), None)
                    cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                    #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                    cst_offset = ((cst_look_up_idx + cst_offset_idx) - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000. 
                    cst_off_idx = cst_look_up_idx + cst_offset_idx                
                    last_ibound_col_n_layers = -3
                    #return None, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist    

                    if cst_off_idx < cst_look_up_idx:
                        print("Pica 3a")           
    
                        #   define the index from which we will look for the coastline 
                        cst_look_up_idx = int((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) + int(((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) / 8) * 7
                        #   find the position of the coastline
                        cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < sea_level), None)
                        cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                        #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                        cst_offset = ((cst_look_up_idx + cst_offset_idx) - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000.
                        cst_off_idx = cst_look_up_idx + cst_offset_idx        

                        if cst_off_idx < cst_look_up_idx:
                            print("Pica 3a")           
        
                            #   define the index from which we will look for the coastline 
                            cst_look_up_idx = int((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) + int(((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) / 8) * 10
                            #   find the position of the coastline
                            cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < sea_level), None)
                            cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                            #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                            cst_offset = ((cst_look_up_idx + cst_offset_idx) - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000.
                            cst_off_idx = cst_look_up_idx + cst_offset_idx       

                            """
                            if cst_off_idx < cst_look_up_idx:
                                print("Pica 3b the island is too small ty kundo")    
                                ls_elev = -3
                                ls_elev_50 = -3   
                                fos_pt_dist = -3
                                model_run = False   
                            """ 
                            if cst_off_idx < cst_look_up_idx:
                                cst_look_up_idx = int(len(avg_topo) / 2)
                                print("Pica 3b the island is too small ty kundo")    
                                model_run = True   
                                
                            else:
                                model_run = True
                                
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
                anchor_topo_elev = final_topo[int(round((cst_off_idx - (200 - x_nasa) * (1000 / topo_res))))]
                #anchor_topo_elev = final_topo[int(round((cst_off_idx - x_nasa * (1000 / topo_res))))]
                y_nasa = anchor_topo_elev - 50.0    
            
                try: 
                    #   first try to find the shelf break and foot of continental slope points
                    print("Pica 5")    
                    offshore_dir = os.path.join(out_dir, '_topo_figs')
                    if not os.path.exists(offshore_dir):
                        os.makedirs(offshore_dir)
                    fos = find_FOS(final_topo, offshore_dir, 'avg_avg', avg = False) 
                    print(fos) 
                    model = cs_model(1, final_topo, avg_nasa_refined, x_nasa, y_nasa, cst_thick_est_avg, end_cst_plain,\
                                     avg_soil_thk_refined, avg_soil_type_refined, avg_wtd_topo_refined, avg_offshore_refined,\
                                     final_topo, final_topo, final_topo, avg_k_soil_refined,\
                                     avg_drn_rate_refined, avg_glhymps_top_lay_refined, avg_glhymps_bot_lay_refined, topo_res, False)
                    model.get_top_bot_lst_find_cst(cst_look_up_idx, fos[1], 0.0)
                    #model.get_top_bot_col_lst(del_col, cs_points_dist, model.cst_idx, smooth = True)
                    #model.get_top_bot_col_lst_top_sys_geo(del_col, cs_points_dist, model.cst_idx, smooth = True)
                    model.get_top_bot_col_lst(del_col, topo_res, model.cst_idx, smooth = True)
                    model.get_top_bot_col_lst_top_sys_geo(del_col, topo_res, model.cst_idx, smooth = True)                    
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
                        model.find_topo_divide_SRM_v2(10, -2000., del_col, del_lay, fos[1], 100., 20., 50000., 200000., 10000, 200000)
                        #model.find_topo_divide(10., -10000., del_col, model.top_elev[model.idx_start], fos[1], 2, 200000.)         #   indicate the value of topographic difference to find the divide 
                    #   in case the island is small cut off the inland part with the start of the model.top_elev list
                    except IndexError:
                        model.find_topo_divide_SRM(10, -2000., del_col, del_lay, fos[1], 20., 50000., 200000., 10000, 200000)
                        #model.find_topo_divide(10., -10000., del_col, model.top_elev[0], fos[1], 2, 200000.)                      
                    
                    #try:
                    #    model.find_topo_divide(10., -10000., del_col, model.top_elev[model.idx_start], fos[1], 2, 200000.)         #   indicate the value of topographic difference to find the divide 
                    #   in case the island is small cut off the inland part with the start of the model.top_elev list
                    #except IndexError:
                    #    model.find_topo_divide(10., -10000., del_col, model.top_elev[0], fos[1], 2, 200000.)  

                    model.dis_input(nrow, delc, del_col, del_lay, nper, perlen, nstp, laycbd, max_depth, cs_points_dist)
                    model.top = model.lay_elev[model.lay_idx]            #   top elevation of the model domain - check if it should be like this!
                    model.botm = model.lay_elev[model.lay_idx + 1:]   
                    #model.zbot = model.botm[-1] 
                    zbot = math.floor((model.bot_elev[-1] / 100.0) * 100.0)
                                 
                    ibound_arr = model.ibound_arr * 1.0
                    ibound_arr[np.abs(ibound_arr) == 0.] = np.nan
                    
                    str_0 = out_name + '.png'
                    plot_dir = os.path.join(out_dir, str_0)           
                    
                    x_start_new = model.x_start# + (model.cst_idx * del_col) / 1000.
                    x_shift = round(model.x_start + (model.cst_idx * del_col) / 1000., 1)
                    x_end_new = model.x_end # + (model.cst_idx * del_col) / 1000.
                           
                    #model_start_idx = final_topo.index(min(final_topo, key=lambda x:abs(x-model.top_elev[0])))   
                    #model_end_idx = final_topo.index(min(final_topo, key=lambda x:abs(x-model.top_elev[-1])))  
                           
                    model_end_idx = int(round(model.idx_end * del_col / 1000.)) * 2 + cst_off_idx + int(round(x_start_new * 2))       
                    #model_start_idx = 400 + int(round(x_start_new * 2))                           
                           
                    topo_lst_plot = final_topo[cst_off_idx + int(round(x_start_new * 2)): model_end_idx + 1]
                    #topo_lst_plot = final_topo[cst_off_idx + int(round(x_start_new * 2)): model_end_idx + 1]

                    #x_axis_plot2 = np.linspace(x_start_new, x_end_new, model.ibound_arr.shape[-1])# - model.col_idx)
                    #x_axis_plot = np.linspace(x_start_new, x_end_new, len(topo_lst_plot))

                    if cst_offset >= 0:
                        x_axis_plot2 = np.linspace(x_start_new - x_shift, x_end_new + x_shift, model.ibound_arr.shape[-1])# - model.col_idx)
                        x_axis_plot = np.linspace(x_start_new - x_shift, x_end_new + x_shift, len(topo_lst_plot))
                    else:
                        x_axis_plot2 = np.linspace(x_start_new - x_shift, x_end_new + x_shift, model.ibound_arr.shape[-1])# - model.col_idx)
                        x_axis_plot = np.linspace(x_start_new - x_shift, x_end_new + x_shift, len(topo_lst_plot))

                    top_elev_arr = np.array(model.top_elev)
                    mask = np.isnan(top_elev_arr)   
                    top_elev_arr[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), top_elev_arr[~mask])
                    
                    model.top_elev = top_elev_arr.tolist()
                           
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
                    ax_all.set_title(out_name)
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
                    topo_dict[out_name] = {'ibound_arr': model.ibound_arr,\
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
                                          'k_soil': model.k_soil,\
                                          'drn_rate': model.drn_dens,\
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
    plot_name_all = os.path.join(out_dir, '__' + out_name + '_all_model_IBOUNDs.png')
    plt.savefig(plot_name_all, dpi = 300, facecolor = 'w', edgecolor = 'w',
        orientation = 'portrait', papertype = None, format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
        frameon = None)
    plt.close()

    #   save the dictionary with the ibound arrays
    dict_save_dir = os.path.join(out_dir, out_name + '_IBOUND.npy')
    np.save(dict_save_dir, topo_dict)

    return topo_dict, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist, fos


""""
g_10 = lst_topo_gebco_10m
g_25 = lst_topo_gebco_25m
g_50 = lst_topo_gebco_50m
g_100 = lst_topo_gebco_100m
g_250 = lst_topo_gebco_250m
g_500 = lst_topo_gebco_500m

gm_10 = lst_gebco_merit_avg_10m
gm_25 = lst_gebco_merit_avg_25m
gm_50 = lst_gebco_merit_avg_50m
gm_100 = lst_gebco_merit_avg_100m
gm_250 = lst_gebco_merit_avg_250m
gm_500 = lst_gebco_merit_avg_500m

dir_out = id_riv_bas_dir

"""



def create_topo_nc(gm_10, gm_25, gm_50, gm_100, gm_250, gm_500, g_10, g_25, g_50, g_100, g_250, g_500,\
                   gcst_10, gcst_25, gcst_50, gcst_100, gcst_250, gcst_500, dir_out):
    
    #   this will just create a final netcdf files for each resolution with the same x coordinates
    nc_10_dir = os.path.join(dir_out, '_topo_10m_avg.nc')
    nc_25_dir = os.path.join(dir_out, '_topo_25m_avg.nc')
    nc_50_dir = os.path.join(dir_out, '_topo_50m_avg.nc')
    nc_100_dir = os.path.join(dir_out, '_topo_100m_avg.nc')
    nc_250_dir = os.path.join(dir_out, '_topo_250m_avg.nc')
    nc_500_dir = os.path.join(dir_out, '_topo_500m_avg.nc')
    
    #   make the x_coord lists
    x_coord_10 = np.arange(-200., 200 + 10. / 1000., 10. / 1000.)
    x_coord_25 = np.arange(-200., 200 + 25. / 1000., 25. / 1000.)
    x_coord_50 = np.arange(-200., 200 + 50. / 1000., 50. / 1000.)
    x_coord_100 = np.arange(-200., 200 + 100. / 1000., 100. / 1000.)
    x_coord_250 = np.arange(-200., 200 + 250. / 1000., 250. / 1000.)
    x_coord_500 = np.arange(-200., 200 + 500. / 1000., 500. / 1000.)
    
    #   create the average topo arrays to be stored in the netcdf files
    arr_topo_gm_10 = np.array(gm_10)
    arr_topo_gm_10[arr_topo_gm_10 == -9999.0] = np.nan
    avg_topo_gm_10 = np.nanmean(arr_topo_gm_10, axis = 0)
    arr_topo_g_10 = np.array(g_10)
    arr_topo_g_10[arr_topo_g_10 == -9999.0] = np.nan
    avg_topo_g_10 = np.nanmean(arr_topo_g_10, axis = 0)
    arr_topo_gcst_10 = np.array(gcst_10)
    arr_topo_gcst_10[arr_topo_gcst_10 == -9999.0] = np.nan
    avg_topo_gcst_10 = np.nanmean(arr_topo_gcst_10, axis = 0)

    arr_topo_gm_25 = np.array(gm_25)
    arr_topo_gm_25[arr_topo_gm_25 == -9999.0] = np.nan
    avg_topo_gm_25 = np.nanmean(arr_topo_gm_25, axis = 0)
    arr_topo_g_25 = np.array(g_25)
    arr_topo_g_25[arr_topo_g_25 == -9999.0] = np.nan
    avg_topo_g_25 = np.nanmean(arr_topo_g_25, axis = 0)
    arr_topo_gcst_25 = np.array(gcst_25)
    arr_topo_gcst_25[arr_topo_gcst_25 == -9999.0] = np.nan
    avg_topo_gcst_25 = np.nanmean(arr_topo_gcst_25, axis = 0)
    
    arr_topo_gm_50 = np.array(gm_50)
    arr_topo_gm_50[arr_topo_gm_50 == -9999.0] = np.nan
    avg_topo_gm_50 = np.nanmean(arr_topo_gm_50, axis = 0)
    arr_topo_g_50 = np.array(g_50)
    arr_topo_g_50[arr_topo_g_50 == -9999.0] = np.nan
    avg_topo_g_50 = np.nanmean(arr_topo_g_50, axis = 0)
    arr_topo_gcst_50 = np.array(gcst_50)
    arr_topo_gcst_50[arr_topo_gcst_50 == -9999.0] = np.nan
    avg_topo_gcst_50 = np.nanmean(arr_topo_gcst_50, axis = 0)    

    arr_topo_gm_100 = np.array(gm_100)
    arr_topo_gm_100[arr_topo_gm_100 == -9999.0] = np.nan
    avg_topo_gm_100 = np.nanmean(arr_topo_gm_100, axis = 0)
    arr_topo_g_100 = np.array(g_100)
    arr_topo_g_100[arr_topo_g_100 == -9999.0] = np.nan
    avg_topo_g_100 = np.nanmean(arr_topo_g_100, axis = 0)
    arr_topo_gcst_100 = np.array(gcst_100)
    arr_topo_gcst_100[arr_topo_gcst_100 == -9999.0] = np.nan
    avg_topo_gcst_100 = np.nanmean(arr_topo_gcst_100, axis = 0)

    arr_topo_gm_250 = np.array(gm_250)
    arr_topo_gm_250[arr_topo_gm_250 == -9999.0] = np.nan
    avg_topo_gm_250 = np.nanmean(arr_topo_gm_250, axis = 0)
    arr_topo_g_250 = np.array(g_250)
    arr_topo_g_250[arr_topo_g_250 == -9999.0] = np.nan
    avg_topo_g_250 = np.nanmean(arr_topo_g_250, axis = 0)
    arr_topo_gcst_250 = np.array(gcst_250)
    arr_topo_gcst_250[arr_topo_gcst_250 == -9999.0] = np.nan
    avg_topo_gcst_250 = np.nanmean(arr_topo_gcst_250, axis = 0)
    
    arr_topo_gm_500 = np.array(gm_500)
    arr_topo_gm_500[arr_topo_gm_500 == -9999.0] = np.nan
    avg_topo_gm_500 = np.nanmean(arr_topo_gm_500, axis = 0)
    arr_topo_g_500 = np.array(g_500)
    arr_topo_g_500[arr_topo_g_500 == -9999.0] = np.nan
    avg_topo_g_500 = np.nanmean(arr_topo_g_500, axis = 0)
    arr_topo_gcst_500 = np.array(gcst_500)
    arr_topo_gcst_500[arr_topo_gcst_500 == -9999.0] = np.nan
    avg_topo_gcst_500 = np.nanmean(arr_topo_gcst_500, axis = 0)
       
    nc_out_10 = xr.Dataset(data_vars = {'gebco_avg' : (('x'), avg_topo_g_10),
                                        'gebco_merit_avg' : (('x'), avg_topo_gm_10),
                                        'gebco_coastal_avg' : (('x'), avg_topo_gcst_10)},
                           coords = {'x' : x_coord_10})
    nc_out_10.to_netcdf(nc_10_dir)      
    
    nc_out_25 = xr.Dataset(data_vars = {'gebco_avg' : (('x'), avg_topo_g_25),
                                        'gebco_merit_avg' : (('x'), avg_topo_gm_25),
                                        'gebco_coastal_avg' : (('x'), avg_topo_gcst_25)},
                           coords = {'x' : x_coord_25})
    nc_out_25.to_netcdf(nc_25_dir)    

    nc_out_50 = xr.Dataset(data_vars = {'gebco_avg' : (('x'), avg_topo_g_50),
                                        'gebco_merit_avg' : (('x'), avg_topo_gm_50),
                                        'gebco_coastal_avg' : (('x'), avg_topo_gcst_50)},
                           coords = {'x' : x_coord_50})
    nc_out_50.to_netcdf(nc_50_dir)    

    nc_out_100 = xr.Dataset(data_vars = {'gebco_avg' : (('x'), avg_topo_g_100),
                                        'gebco_merit_avg' : (('x'), avg_topo_gm_100),
                                        'gebco_coastal_avg' : (('x'), avg_topo_gcst_100)},
                           coords = {'x' : x_coord_100})
    nc_out_100.to_netcdf(nc_100_dir)    

    nc_out_250 = xr.Dataset(data_vars = {'gebco_avg' : (('x'), avg_topo_g_250),
                                        'gebco_merit_avg' : (('x'), avg_topo_gm_250),
                                        'gebco_coastal_avg' : (('x'), avg_topo_gcst_250)},
                           coords = {'x' : x_coord_250})
    nc_out_250.to_netcdf(nc_250_dir)    

    nc_out_500 = xr.Dataset(data_vars = {'gebco_avg' : (('x'), avg_topo_g_500),
                                        'gebco_merit_avg' : (('x'), avg_topo_gm_500),
                                        'gebco_coastal_avg' : (('x'), avg_topo_gcst_500)},
                           coords = {'x' : x_coord_500})
    nc_out_500.to_netcdf(nc_500_dir)        
    



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
    #cst_type = 'avg_avg'
    #avg= True
    
    #   clean the lists - if the values are = non_val then change to NaN not to screw up the calculation of averages
    arr_topo = np.array(topo_lst)
    #   distinguish between the cases when all the topographical profiles are supplied vs. when only one topo list is
    if avg and len(arr_topo.shape) > 1:
        avg_topo_lst = np.nanmean(arr_topo, axis = 0)
    elif avg:
        avg_topo_lst = arr_topo
    else:
        avg_topo_lst = topo_lst
    #x_axis = np.linspace(0.0, 200.0, 401) 
    #x_axis = np.linspace(0.0, 200.0, len(topo_lst)) 
    x_axis = np.linspace(0.0, 200.0, int((len(topo_lst) / 2 + 1))) 
    
    #   calculate the gradient and second derivative
    #grad = np.gradient(avg_topo_lst[400:])
    grad = np.gradient(avg_topo_lst[int((len(topo_lst) - 1) / 2) :])
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
        #loc_max_plt.append([x_axis[loc_max[0][a]], avg_topo_lst[400 + loc_max[0][a]]])
        loc_max_plt.append([x_axis[loc_max[0][a]], avg_topo_lst[int((len(topo_lst) - 1) / 2) + loc_max[0][a]]])
    
    for b in range(loc_min[0].shape[0]):
        #loc_min_plt.append([x_axis[loc_min[0][b]], avg_topo_lst[400 + loc_min[0][b]]])
        loc_min_plt.append([x_axis[loc_min[0][b]], avg_topo_lst[int((len(topo_lst) - 1) / 2) + loc_min[0][b]]])
        
    for c in range(loc_minmax.shape[0]):
        #loc_minmax_plt1.append([avg_topo[loc_minmax[c]], avg_topo[loc_minmax[c]]])
        #loc_minmax_plt1.append(avg_topo_lst[400 + loc_minmax[c]])
        #loc_minmax_plt2.append(avg_topo_lst[400 + loc_minmax[c]])
        loc_minmax_plt1.append(avg_topo_lst[int((len(topo_lst) - 1) / 2) + loc_minmax[c]])
        loc_minmax_plt2.append(avg_topo_lst[int((len(topo_lst) - 1) / 2) + loc_minmax[c]])
    
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
        #ax1.plot(x_axis, avg_topo_lst[400:], linewidth = 2, color = "blue", label = 'Original topography/bathymetry')
        ax1.plot(x_axis, avg_topo_lst[int((len(topo_lst) - 1) / 2):], linewidth = 2, color = "blue", label = 'Original topography/bathymetry')
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
def plot_all_topo_profiles(topo_lst, ate_lst, fos_pt, topo_name, subreg_name, out_dir):

    #topo_lst = lst_gebco_merit_avg_500m
    #ate_lst = lst_cs_thk_henry
    #coscat_id = id_coscat
    #cst_type = 'henry'
    #out_dir = coscat_dir
    #fos_pt = find_FOS(topo_lst, None, cst_type, avg = True)
    
    plot_title = 'All topographical profiles (' + topo_name + ') for SUBREG ' + subreg_name    
    plot_dir = os.path.join(out_dir, '_' + topo_name + '_all_topo_profiles.png')    
    
    #   calculate the average topo_lst 
    arr_topo = np.array(topo_lst)
    avg_topo = np.nanmean(arr_topo, axis = 0)      
    avg_ate = sum(ate_lst)/len(ate_lst)    
    
    mask = np.isnan(avg_topo)   
    avg_topo[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), avg_topo[~mask])
    
    avg_topo = avg_topo.tolist()        
    
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



#   function to connect to the database
def connect_to_dtbase(db_name, db_user, db_host, db_pass):
    #   create connection string
    conn_string = str("dbname=%s user=%s host=%s password=%s") % (db_name, db_user, db_host, db_pass)
    try:
        conn = psycopg2.connect(conn_string)
        print("Successfully connected to database : " + str(db_name))
        #   set the cursor
        cur = conn.cursor()
        return conn, cur
    except:
        print("I am unable to connect to database : " + str(db_name))


#  Function to create a database table with specified columns
#  Intended for creation of values extracted from input data at cross-section
#  points location..
def sql_create_table(db_conn, db_cursor, tb_name, id_name, tb_cols):
    ##  first check that the cursor exists
    if db_cursor:
        ##  if yes then create the SQL and run it
        sql_command = "CREATE TABLE IF NOT EXISTS %s (\
                            %s serial PRIMARY KEY, %s);" % (tb_name, id_name, tb_cols)
        #   print sql_command
        db_cursor.execute(sql_command)
        db_conn.commit()
    else:
        print('Database cursor doesnt exist!')


#  Function to create a list of columns depending on distance between coastal
#  point and the points on the cross-section
#  col_count (equal to n_points), dist - distance between points on cross-section
def create_column_string(col_count, dist, col_type):
    col_string = ""
    for i in range(col_count + 1):
        dist_to_cs = i - (col_count/2)
        ##  for negative distance values
        if dist_to_cs < 0:
            col_string += 'dist_minus_' + str(int((abs(dist_to_cs) * dist))) + 'm ' + col_type + ','
        ##  for coastal point - distance is 0..
        elif dist_to_cs == 0:
            col_string += 'dist_' + str(int(dist_to_cs * dist)) + 'm ' + col_type + ','
        ##  for positive distance values
        else:
            col_string += 'dist_plus_' + str(int(dist_to_cs * dist)) + 'm ' + col_type + ','
    ##  remove the last comma in the string to get correct SQL command
    col_string = col_string[:-1]
    return col_string


#  Function for insert SQL into specific columns of the table
def sql_insert(db_conn, db_cursor, db_table, db_columns, ins_vals):
    ##  first check that the cursor exists
    if db_cursor:
        ##  if yes then perform the insert SQL
        sql_command = "INSERT INTO %s (%s) VALUES (%s)" % (db_table, db_columns, ins_vals)
        db_cursor.execute(sql_command)
        db_conn.commit()
    else:
        print('Database cursor doesnt exist!')


#  Function for inserting values in existing rows - updating the row
def sql_update(db_conn, db_cursor, db_table, db_columns, ins_vals, where_condition):
    ##  first check that the cursor exists
    #print db_cursor, db_table, db_columns, ins_vals, where_condition
    if db_cursor:
        ##  if yes then perform the insert SQL
        sql_command = "UPDATE %s SET %s = %s WHERE %s" % (db_table, db_columns, ins_vals, where_condition)
        #print sql_command
        db_cursor.execute(sql_command)
        db_conn.commit()
    else:
        print('Database cursor doesnt exist!')


#   function to calculate the coordinates from triangle information
def coord_from_triangle(x0, y0, x1, y1, b_len):
    c_len = math.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
    ##  abs value in case (x0-x1) < 0 - fixes the wrong angle the line has in case d_len < 0
    d_len = abs(x0 - x1)
    # print x0, y0, x1, y1, b_len, c_len, d_len
    delta = math.degrees(math.acos(d_len / c_len))
    omega = 90 - delta
    ##  calculate the coordinates for both directions from the point_5km coastal point
    ##  conditions below make sure that the cross-section is perpendicular
    if x0 > x1 and y0 > y1:
        C_coord_x_sea = y0 + math.sin(math.radians(omega)) * b_len
        C_coord_y_sea = x0 - math.cos(math.radians(omega)) * b_len
        C_coord_x_land = y0 - math.sin(math.radians(omega)) * b_len
        C_coord_y_land = x0 + math.cos(math.radians(omega)) * b_len
    elif x0 < x1 and y0 < y1:
        C_coord_x_sea = y0 + math.sin(math.radians(omega)) * b_len
        C_coord_y_sea = x0 - math.cos(math.radians(omega)) * b_len
        C_coord_x_land = y0 - math.sin(math.radians(omega)) * b_len
        C_coord_y_land = x0 + math.cos(math.radians(omega)) * b_len
    else:
        C_coord_x_sea = y0 - math.sin(math.radians(omega)) * b_len
        C_coord_y_sea = x0 - math.cos(math.radians(omega)) * b_len
        C_coord_x_land = y0 + math.sin(math.radians(omega)) * b_len
        C_coord_y_land = x0 + math.cos(math.radians(omega)) * b_len
    ##  return all the calculated coordinates
    return C_coord_x_land, C_coord_y_land, C_coord_x_sea, C_coord_y_sea


#   Function to get the closest raster pixel value to point coordinates
#   this fc is used in case the (coastal) point lies on a novalue pixel..
#   avg_type -> either average or maximal frequency value (for land_cover etc)
#   in case the raster provides data in classed values
def find_closest_neighbour(point_x_col, point_y_row, r_band, r_noval ,raster_x_size, raster_y_size, avg_type):
    ##  check if the col or row are not out of bounds of the raster extent
    if point_x_col >= raster_x_size:
        point_x_col_min_1 = point_x_col - 1
        point_x_col = point_x_col - raster_x_size
        point_x_col_plus_1 = point_x_col + 1
    elif point_x_col <= 0:
        point_x_col_plus_1 = point_x_col + 1
        point_x_col = point_x_col + raster_x_size - 1
        point_x_col_min_1 = point_x_col - 1
    elif point_y_row >= raster_y_size:
        point_y_row = point_y_row - raster_y_size
    elif point_y_row < 0:
        point_y_row = point_y_row + raster_y_size
    else:
        point_x_col_plus_1 = point_x_col + 1
        if point_x_col_plus_1 >= raster_x_size:
            point_x_col_plus_1 = point_x_col_plus_1 - raster_x_size
        point_x_col_min_1 = point_x_col - 1
    ##  get an average value from the surrounding pixels, omit the noval pixels
    non_noval_list = []
    non_noval_list.append(r_band.ReadAsArray(point_x_col_min_1, point_y_row - 1, 1, 1)[0][0])
    non_noval_list.append(r_band.ReadAsArray(point_x_col_min_1, point_y_row, 1, 1)[0][0])
    non_noval_list.append(r_band.ReadAsArray(point_x_col_min_1, point_y_row + 1, 1, 1)[0][0])
    non_noval_list.append(r_band.ReadAsArray(point_x_col, point_y_row - 1, 1, 1)[0][0])
    non_noval_list.append(r_band.ReadAsArray(point_x_col, point_y_row + 1, 1, 1)[0][0])
    non_noval_list.append(r_band.ReadAsArray(point_x_col_plus_1, point_y_row - 1, 1, 1)[0][0])
    non_noval_list.append(r_band.ReadAsArray(point_x_col_plus_1, point_y_row, 1, 1)[0][0])
    non_noval_list.append(r_band.ReadAsArray(point_x_col_plus_1, point_y_row + 1, 1, 1)[0][0])
    ##  remove the noval values from the list
    non_noval_list = [x for x in non_noval_list if x != r_noval]
    #print non_noval_list
    ##  calculate the average value or max_freq values
    if avg_type == 'avg':
        try:
            pixel_value = sum(non_noval_list) / float(len(non_noval_list))
        ##  in case the non_noval_list is empty
        except ZeroDivisionError:
            pixel_value = -9999
    elif avg_type == 'max_freq':
        try:
            pixel_value = Counter(non_noval_list).most_common(1)[0][0]
        ##  in case the non_noval_list is empty
        except IndexError:
            pixel_value = -9999
    ##  return the calculated pixel_value
    return pixel_value


#  Function to check the position to the equator (south or north)
#  if north ( > 10 000 000m) then substract it and move the zone_let from M to N
def equator_position(coords_1, coords_2, coords_3, coords_4, zone_num, zone_let):
    #print coords_1, coords_2, coords_3, coords_4, zone_num, zone_let
    if coords_1 > 10000000 and coords_3 > 10000000:
        zone_let_c_land, zone_let_c_sea = 'N', 'N'
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2 - 10000000, zone_num, zone_let_c_land)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4 - 10000000, zone_num, zone_let_c_sea)

    elif coords_1 > 10000000 and coords_3 < 10000000:
        zone_let_c_land, zone_let_c_sea = 'N', zone_let
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2 - 10000000, zone_num, zone_let_c_land)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4, zone_num, zone_let_c_sea)

    elif coords_1 < 10000000 and coords_3 > 10000000:
        zone_let_c_land, zone_let_c_sea = zone_let , 'N'
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2, zone_num, zone_let_c_land)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4 - 10000000, zone_num, zone_let_c_sea)

    else:
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2, zone_num, zone_let)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4, zone_num, zone_let)
    return c_coords_wgs_land, c_coords_wgs_sea


#  Function to check the position to the equator (south or north)
#  if north ( > 10 000 000m) then substract it and move the zone_let from M to N
def equator_position_angles(coords_1, coords_2, coords_3, coords_4, zone_num, zone_let):
    if coords_2 < 0 and coords_4 < 0:
        #print ('-2')
        zone_let_c_land, zone_let_c_sea = 'M', 'M'
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2 + 10000000, zone_num, zone_let_c_land)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4 + 10000000, zone_num, zone_let_c_sea)

    elif coords_2 < 0 and coords_4 > 0:
        #print ('-1')
        zone_let_c_land, zone_let_c_sea = 'M', zone_let
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2 + 10000000, zone_num, zone_let_c_land)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4, zone_num, zone_let_c_sea)

    elif coords_2 > 0 and coords_4 < 0:
        #print ('0')
        zone_let_c_land, zone_let_c_sea = zone_let, 'M'
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2, zone_num, zone_let)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4 + 10000000, zone_num, zone_let_c_sea)

    elif coords_2 > 10000000 and coords_4 > 10000000:
        #print ('1')
        zone_let_c_land, zone_let_c_sea = 'N', 'N'
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2 - 10000000, zone_num, zone_let_c_land)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4 - 10000000, zone_num, zone_let_c_sea)

    elif coords_2 > 10000000 and coords_4 < 10000000:
        #print ('2')
        zone_let_c_land, zone_let_c_sea = 'N', zone_let
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2 - 10000000, zone_num, zone_let_c_land)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4, zone_num, zone_let_c_sea)

    elif coords_2 < 10000000 and coords_4 > 10000000:
        #print ('3')
        zone_let_c_land, zone_let_c_sea = zone_let , 'N'
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2, zone_num, zone_let_c_land)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4 - 10000000, zone_num, zone_let_c_sea)

    else:
        #print ('4')
        c_coords_wgs_land = utm.to_latlon(coords_1, coords_2, zone_num, zone_let)
        c_coords_wgs_sea = utm.to_latlon(coords_3, coords_4, zone_num, zone_let)

    return c_coords_wgs_land, c_coords_wgs_sea


#   function to calculate the coordinates from triangle information
#   works for the non-perpendicular cross-sections, need to specify alfa
#   X0 is the coastal point and X1 is a cross-section point
def coord_from_triangle_alfa(x0, y0, x1, y1, b_len, alfa):
    ##  calculate the distance between the coastal point and the cross-section point
    c_len = math.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
    #print c_len
    ##  calculate the difference between y coordinates of two points
    a_1 = abs(y0 - y1)
    #print a_1
    ##  calculate angle alfa_1
    alfa_1 = math.degrees(math.asin(a_1 / c_len))
    #print alfa_1
    ##  calculate alfa_2 in radians
    alfa_2 = math.radians(90 - alfa_1 - alfa)
    #print alfa_2, 90 - alfa_1 - alfa
    ##  with alfa_2 calculate the distances b1 and b2
    b_1 = b_len * math.cos(alfa_2)
    b_2 = b_len * math.sin(alfa_2)
    #print b_1, b_2
    ##  with these two distances it is now possible to calculate the coordinates
    if x0 > x1 and y0 > y1:
        coord_x_sea = x0 + b_1
        coord_y_sea = y0 - b_2
        coord_x_land = x0 - b_1
        coord_y_land = y0 + b_2
        if alfa == 45:
            return coord_x_sea, coord_y_sea, coord_x_land, coord_y_land
        else:
            return coord_x_land, coord_y_land, coord_x_sea, coord_y_sea
    elif x0 < x1 and y0 < y1:
        coord_x_sea = x0 - b_1
        coord_y_sea = y0 + b_2
        coord_x_land = x0 + b_1
        coord_y_land = y0 - b_2
        if alfa == 135:
            return coord_x_sea, coord_y_sea, coord_x_land, coord_y_land
        else:
            return coord_x_land, coord_y_land, coord_x_sea, coord_y_sea
    #elif x0 < x1 and y0 > y1:
    #    coord_x_sea = x0 - b_1
    #    coord_y_sea = y0 + b_2
    #    coord_x_land = x0 + b_1
    #    coord_y_land = y0 - b_2
    else:
        coord_x_sea = x0 - b_1
        coord_y_sea = y0 - b_2
        coord_x_land = x0 + b_1
        coord_y_land = y0 + b_2
        if alfa == 135:
            return coord_x_sea, coord_y_sea, coord_x_land, coord_y_land
        else:
            return coord_x_land, coord_y_land, coord_x_sea, coord_y_sea

#   function to calculate the average value of a profile point
#   creates a perpendicular cross-section of the cross-section at given
#   point and gives a list of coordinates that can then be used for
#   extracting values from a raster
def avg_cs_point_coord(x0, y0, x1, y1, b_len):
    c_len = math.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
    ##  in case c_len = 0 change the value to prevent division by zero (this is for the
    ##  coastal point averages!)
    if c_len == 0:
        c_len = 0.0000001
    ##  abs value in case (x0-x1) < 0 - fixes the wrong angle the line has in case d_len < 0
    d_len = abs(x0 - x1)
    # print x0, y0, x1, y1, b_len, c_len, d_len
    delta = math.degrees(math.acos(d_len / c_len))
    omega = 90 - delta
    ##  calculate the coordinates for both directions from the point_5km coastal point
    ##  conditions below make sure that the cross-section is perpendicular
    if x0 > x1 and y0 > y1:
        C_coord_x_sea = y1 + math.sin(math.radians(omega)) * b_len
        C_coord_y_sea = x1 - math.cos(math.radians(omega)) * b_len
        C_coord_x_land = y1 - math.sin(math.radians(omega)) * b_len
        C_coord_y_land = x1 + math.cos(math.radians(omega)) * b_len
    elif x0 < x1 and y0 < y1:
        C_coord_x_sea = y1 + math.sin(math.radians(omega)) * b_len
        C_coord_y_sea = x1 - math.cos(math.radians(omega)) * b_len
        C_coord_x_land = y1 - math.sin(math.radians(omega)) * b_len
        C_coord_y_land = x1 + math.cos(math.radians(omega)) * b_len
    else:
        C_coord_x_sea = y1 - math.sin(math.radians(omega)) * b_len
        C_coord_y_sea = x1 - math.cos(math.radians(omega)) * b_len
        C_coord_x_land = y1 + math.sin(math.radians(omega)) * b_len
        C_coord_y_land = x1 + math.cos(math.radians(omega)) * b_len
    ##  return all the calculated coordinates
    return C_coord_x_land, C_coord_y_land, C_coord_x_sea, C_coord_y_sea


#   function to extract values from raster. In case the point coordinates are in
#   the edge of the raster (x > 180 or x < 180) read raster values from the other
#   'side' of the raster. E.g. if the raster column is 2 pixels higher than the
#   total of raster columns than read the column 2 of the raster = earth is round!
def read_raster_val(in_rb, in_gt, raster_x_size, raster_y_size, point_coord_x, point_coord_y):
    ##  calculate the col and row for the point coordinates
    px = int((point_coord_x - in_gt[0]) / in_gt[1]) #x pixel
    py = int((point_coord_y - in_gt[3]) / in_gt[5]) #y pixel
    #print point_coord_x, point_coord_y
    #print px, py
    #print('-----------------------------------------------------')
    ##  check if the col and row are within the raster extent range, if not
    ##  change them so they are - earth is round!
    if px >= raster_x_size:
        px = px - raster_x_size
    elif px < 0:
        px = px + raster_x_size
    elif py >= raster_y_size:
        py = py - raster_y_size
    elif py < 0:
        py = py + raster_y_size
    ##  get the raster value in the px and py
    #print px, py
    pixel_val = in_rb.ReadAsArray(px, py, 1, 1)[0][0]
    ##  return the pixel value
    return pixel_val




"""
Function that calculates lognormal distribution of groundwater recharge values based on input consisting of
list of files with gw rch profiles for each coastal point in the SRM,

gw_lst_in = gw_rch_files
srm_dir = subreg_input_dir
x_st = x_start
x_en = x_end

"""
def get_gw_rch_lognorm(srm_dir, gw_lst_in, topo_files, x_st, x_en):
    
    #   1) get list of time steps for which we will create gw_rch lognorm distribution, just use the first value in the list
    ts_lst = xr.open_dataset(os.path.join(srm_dir, gw_lst_in[0]))['time'].values.tolist()
    cs_id_lst = [f.split('_')[1] for f in topo_files]
    
    #   2) create an empty list where all the lognorm values for each cs will be stored
    main_lst = []
    
    #   3) loop through all the cs points of the SRM and all the time steps 
    for cs_pt_id in cs_id_lst:
        
    #   a) first check the coastal direction so the gw_rch lists are rightly positioned    
        topo_name = [f for f in topo_files if str(cs_pt_id) in f][0]
        topo_lst = xr.open_dataset(os.path.join(srm_dir, topo_name))['GEBCO_elev_AVG'].values.tolist()
        cs_topo = topo_lst[::5]
        
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
    
    #   b) create an empty list where we will store the lognorm values for each time step
        cs_lst = []
        gw_rch_file = [f for f in gw_lst_in if str(cs_pt_id) in f][0]
        gw_rch_nc = xr.open_dataset(os.path.join(srm_dir, gw_rch_file))
        for ts in ts_lst:
            gw_rch_ts = gw_rch_nc.sel(time = ts)['gw_rch_AVG'].values.tolist()
            if cont_direction == 'right':
                #   reverse the list of values so it seems like ocean is on the right hand side
                gw_rch_ts = list(reversed(gw_rch_ts))
            gw_rch_ts = gw_rch_ts[max(0, int(2000 + x_st * 10)) : int(2000 + x_en * 10)]
            gw_rch_ts = [i for i in gw_rch_ts if not math.isnan(i)]
            cs_lst.append(gw_rch_ts)
        main_lst.append(cs_lst)
    
    #   4) for each time step get the mean and stdev 
    mu_lst, std_lst = [], []
    for a in range(len(ts_lst)):
        ts_gw_rch_lst = []
        for b in range(len(main_lst)):
            ts_gw_rch_lst.append(main_lst[b][a])
    
        flat_list = [item for sublist in ts_gw_rch_lst for item in sublist]
        arr = np.array(flat_list)
        #   make sure that the type is float 
        arr = arr.astype(np.float)
        arr[arr < 0] = np.nan
        #   next, remove the Nan values
        arr = arr[np.logical_not(np.isnan(arr))]
        #   then flatten the array to 1D
        arr = arr.flatten()
        arr_ln = np.log(arr)
        
        # Fit a normal distribution to the data and round the mean and std values
        try:
            mu, std = norm.fit(arr_ln)    
            mu_lst.append(round(mu, 4))
            std_lst.append(round(std, 4))
        except RuntimeError:
            #   it fails when there are negative values in the log array
            arr_ln_pos = arr_ln[arr_ln >= 0]
            mu, std = norm.fit(arr_ln_pos)    
            mu_lst.append(round(mu, 4))
            std_lst.append(round(std, 4))
            pass
        
    #Plotting:
    #plt.plot(ts_lst, [np.exp(mu) for mu in mu_lst]) #mean curve.
    return mu_lst, std_lst 
    #return [np.exp(mu) for mu in mu_lst], [np.exp(std) for std in std_lst]    








def get_riv_input(cs_id, flo1k, cond_10, cond_05, cond_01, pcr_cond, pcr_stage, pcr_rivbot,\
                  gaia_width, gaia_depth, x_st, y_st, px_x, px_y, cs_plain_width, cont_direction):

    ##  first open the raster file via GDAL
    #raster_in = gdal.Open(riv_rast_dir)
    #gt = raster_in.GetGeoTransform()
    #rb = raster_in.GetRasterBand(1)
    #noval = rb.GetNoDataValue()
    #rx = raster_in.RasterXSize
    #ry = raster_in.RasterYSize

    flo1k_in = gdal.Open(flo1k)
    gt_flo = flo1k_in.GetGeoTransform()
    rb_flo = flo1k_in.GetRasterBand(1)
    noval_flo = rb_flo.GetNoDataValue()
    rx_flo = flo1k_in.RasterXSize
    ry_flo = flo1k_in.RasterYSize    

    cond_10_in = gdal.Open(cond_10)
    gt_cond_10 = cond_10_in.GetGeoTransform()
    rb_cond_10 = cond_10_in.GetRasterBand(1)
    noval_cond_10 = rb_cond_10.GetNoDataValue()
    rx_cond_10 = cond_10_in.RasterXSize
    ry_cond_10 = cond_10_in.RasterYSize    

    cond_05_in = gdal.Open(cond_05)
    gt_cond_05 = cond_05_in.GetGeoTransform()
    rb_cond_05 = cond_05_in.GetRasterBand(1)
    noval_cond_05 = rb_cond_05.GetNoDataValue()
    rx_cond_05 = cond_05_in.RasterXSize
    ry_cond_05 = cond_05_in.RasterYSize    
    
    cond_01_in = gdal.Open(cond_01)
    gt_cond_01 = cond_01_in.GetGeoTransform()
    rb_cond_01 = cond_01_in.GetRasterBand(1)
    noval_cond_01 = rb_cond_01.GetNoDataValue()
    rx_cond_01 = cond_01_in.RasterXSize
    ry_cond_01 = cond_01_in.RasterYSize    

    pcr_cond_in = gdal.Open(pcr_cond)
    gt_pcr_cond = pcr_cond_in.GetGeoTransform()
    rb_pcr_cond = pcr_cond_in.GetRasterBand(1)
    noval_pcr_cond = rb_pcr_cond.GetNoDataValue()
    rx_pcr_cond = pcr_cond_in.RasterXSize
    ry_pcr_cond = pcr_cond_in.RasterYSize    

    pcr_stage_in = gdal.Open(pcr_stage)
    gt_pcr_stage = pcr_stage_in.GetGeoTransform()
    rb_pcr_stage = pcr_stage_in.GetRasterBand(1)
    noval_pcr_stage = rb_pcr_stage.GetNoDataValue()
    rx_pcr_stage = pcr_stage_in.RasterXSize
    ry_pcr_stage = pcr_stage_in.RasterYSize    
    
    pcr_rivbot_in = gdal.Open(pcr_rivbot)
    gt_pcr_rivbot = pcr_rivbot_in.GetGeoTransform()
    rb_pcr_rivbot = pcr_rivbot_in.GetRasterBand(1)
    noval_pcr_rivbot = rb_pcr_rivbot.GetNoDataValue()
    rx_pcr_rivbot = pcr_rivbot_in.RasterXSize
    ry_pcr_rivbot = pcr_rivbot_in.RasterYSize        
    
    gaia_width_in = gdal.Open(gaia_width)
    gt_gaia_width = gaia_width_in.GetGeoTransform()
    rb_gaia_width = gaia_width_in.GetRasterBand(1)
    noval_gaia_width = rb_gaia_width.GetNoDataValue()
    rx_gaia_width = gaia_width_in.RasterXSize
    ry_gaia_width = gaia_width_in.RasterYSize       
    
    gaia_depth_in = gdal.Open(gaia_depth)
    gt_gaia_depth = gaia_depth_in.GetGeoTransform()
    rb_gaia_depth = gaia_depth_in.GetRasterBand(1)
    noval_gaia_depth = rb_gaia_depth.GetNoDataValue()
    rx_gaia_depth = gaia_depth_in.RasterXSize
    ry_gaia_depth = gaia_depth_in.RasterYSize   

    lst_to_csv = []

    n_points = 799
    max_dist = 399500
    #   define the number of points and distance from the cross-section point
    n_points_avg = 5
    avg_dist = 2500

    #db_name = "'_coastal_dbase_sed_thick_valid'"
    #db_name = "'cs_db_deltas'"
    db_name_world = "'cs_db_v1'"
    db_host = "'localhost'"
    db_user = "'postgres'"
    db_pass = "'postgres'"
    
    #table_db = "coastline_coscat"
    #table_db = "cs_coastal_types"
    coastline_points_raw_table = 'ne_10m_coastline_points_5km'#'ne_points_5km'#'ne_points_5km_test_raster_edge'##'ne_coatline_points_5km_test_avg'#'ne_coatline_points_5km_test_avg'
   
    #dbase_connect = ws.dbase_tools(db_name, db_user, db_host, db_pass)
    dbase_connect_world = ws.dbase_tools(db_name_world, db_user, db_host, db_pass)

    ##  define the database cursor - necessary to move around the dbase tables
    dtbase_cur_world = dbase_connect_world.cur

    ##   first run query to find all the 5km points
    sql_string_select_5km_points = "SELECT id_cs, ST_AsText((ST_Dump(%s.geom)).geom),\
                                    ST_X(geom), ST_Y(geom) AS the_POINT_geom FROM %s WHERE gid = %s"\
                                    % (coastline_points_raw_table, coastline_points_raw_table, cs_id)

    ##  run the sql and get all the coastline points from the raw table
    dtbase_cur_world.execute(sql_string_select_5km_points)
    point_5km = dtbase_cur_world.fetchall()[0]

    point_id = point_5km[0]
    point_geom = point_5km[1]
    point_x, point_y = point_5km[2], point_5km[3]

    ##  create the sql string to find closest node to point_5km
    sql_string_select_points = "SELECT gid, ST_AsText(geom), ST_X(geom), ST_Y(geom) FROM ne_coastline_nodes_test\
                                ORDER BY ne_coastline_nodes_test.geom <->\
                                ST_GeomFromText('%s') LIMIT 1;" % (point_geom)
    ##  open new database cursor
    dtbase_cur2 = dbase_connect_world.cur
    ##  run the sql and get results
    dtbase_cur2.execute(sql_string_select_points)
    closest_node = dtbase_cur2.fetchall()

    ##  assign the results of the query to variables below
    closest_node_x = closest_node[0][2]
    closest_node_y = closest_node[0][3]

    ##  calculate the distance between the point_5m and closest node
    ##  to do this we convert the lat/lon coordinates to UTM coordinates (using the utm package)
    point_5km_utm = utm.from_latlon(point_y, point_x)
    closest_node_utm = utm.from_latlon(closest_node_y, closest_node_x)
    zone_num = point_5km_utm[2]
    zone_let = point_5km_utm[3]
    
    ##  check that zone_num is between 1 and 60
    if zone_num > 60:
        zone_num = zone_num - 60

    def find_pix_val(rb, gt, rx, ry, noval, pt_x, pt_y):
        #px = int((pt_x - x_st) / px_x) #x pixel
        #py = int((pt_y - y_st) / px_y) #y pixel
        px = int((pt_x - gt[0]) / gt[1]) #x pixel
        py = int((pt_y - gt[3]) / gt[5]) #y pixel        
        
        ##  check if the col and row are within the raster extent range, if not
        ##  change them so they are - earth is round!
        if px >= rx:
            px = px - rx
        elif px < 0:
            px = px + rx
        elif py >= ry:
            py = py - ry
        elif py < 0:
            py = py + ry
        ##  get the raster value in the px and py
        #rb.ReadAsArray(px, py, 1, 1).astype(np.float)

        pixel_val = rb.ReadAsArray(px, py, 1, 1)[0][0]
    
        if pixel_val == noval or math.isnan(pixel_val):
            pixel_val = 0.0

        return pixel_val

    #pixel_val_cs_point = find_pix_val(point_x, point_y)
    #lst_to_csv.append([cs_id, 'cst', 0.0, point_x, point_y, pixel_val_cs_point])

    pixel_val_cs_point_flo = find_pix_val(rb_flo, gt_flo, rx_flo, ry_flo, noval_flo, point_x, point_y)
    pixel_val_cs_point_cond_10 = find_pix_val(rb_cond_10, gt_cond_10, rx_cond_10, ry_cond_10, noval_cond_10, point_x, point_y)
    pixel_val_cs_point_cond_05 = find_pix_val(rb_cond_05, gt_cond_05, rx_cond_05, ry_cond_05, noval_cond_05, point_x, point_y)
    pixel_val_cs_point_cond_01 = find_pix_val(rb_cond_01, gt_cond_01, rx_cond_01, ry_cond_01, noval_cond_01, point_x, point_y)
    pixel_val_cs_point_pcr_cond = find_pix_val(rb_pcr_cond, gt_pcr_cond, rx_pcr_cond, ry_pcr_cond, noval_pcr_cond, point_x, point_y)
    pixel_val_cs_point_pcr_stage = find_pix_val(rb_pcr_stage, gt_pcr_stage, rx_pcr_stage, ry_pcr_stage, noval_pcr_stage, point_x, point_y)    
    pixel_val_cs_point_pcr_rivbot = find_pix_val(rb_pcr_rivbot, gt_pcr_rivbot, rx_pcr_rivbot, ry_pcr_rivbot, noval_pcr_rivbot, point_x, point_y)  
    pixel_val_cs_point_gaia_width = find_pix_val(rb_gaia_width, gt_gaia_width, rx_gaia_width, ry_gaia_width, noval_gaia_width, point_x, point_y) 
    pixel_val_cs_point_gaia_depth = find_pix_val(rb_gaia_depth, gt_gaia_depth, rx_gaia_depth, ry_gaia_depth, noval_gaia_depth, point_x, point_y) 

    lst_to_csv.append([cs_id, 'cst', 0.0, point_x, point_y, pixel_val_cs_point_flo, pixel_val_cs_point_flo, pixel_val_cs_point_cond_01,\
                       pixel_val_cs_point_cond_05, pixel_val_cs_point_cond_10, pixel_val_cs_point_pcr_cond, pixel_val_cs_point_pcr_stage,\
                       pixel_val_cs_point_pcr_rivbot, pixel_val_cs_point_gaia_width, pixel_val_cs_point_gaia_depth])

    for a in range(n_points):
        try:
            ##  calculate the b_length
            b_length = max_dist - a * (max_dist/n_points)
            
            #   check if the distance is larger than the coastal plain(*1000 to km)
            if b_length > cs_plain_width * 1000.:
                pass
            else:            
                c_coords_utm = coord_from_triangle(point_5km_utm[0], point_5km_utm[1], closest_node_utm[0], closest_node_utm[1], b_length)
    
                #print ('c_coords :  '), c_coords_utm
                ##  transform the coordinates back to wgs84, position of the equator is taken into account
                c_coords_wgs = equator_position_angles(c_coords_utm[1], c_coords_utm[0], c_coords_utm[3], c_coords_utm[2], zone_num, zone_let)
    
                c_coord_wgs_x_land, c_coord_wgs_y_land, c_coord_wgs_x_sea, c_coord_wgs_y_sea = \
                c_coords_wgs[0][1], c_coords_wgs[0][0], c_coords_wgs[1][1], c_coords_wgs[1][0]
    
                ##  read values from the GEBCO raster
                #pixel_val_land = find_pix_val(c_coord_wgs_x_land, c_coord_wgs_y_land)
                #pixel_val_sea = find_pix_val(c_coord_wgs_x_sea, c_coord_wgs_y_sea)

                pixel_val_cs_point_flo_land = find_pix_val(rb_flo, gt_flo, rx_flo, ry_flo, noval_flo, c_coord_wgs_x_land, c_coord_wgs_y_land)
                pixel_val_cs_point_cond_10_land = find_pix_val(rb_cond_10, gt_cond_10, rx_cond_10, ry_cond_10, noval_cond_10, c_coord_wgs_x_land, c_coord_wgs_y_land)
                pixel_val_cs_point_cond_05_land = find_pix_val(rb_cond_05, gt_cond_05, rx_cond_05, ry_cond_05, noval_cond_05, c_coord_wgs_x_land, c_coord_wgs_y_land)
                pixel_val_cs_point_cond_01_land = find_pix_val(rb_cond_01, gt_cond_01, rx_cond_01, ry_cond_01, noval_cond_01, c_coord_wgs_x_land, c_coord_wgs_y_land)
                pixel_val_cs_point_pcr_cond_land = find_pix_val(rb_pcr_cond, gt_pcr_cond, rx_pcr_cond, ry_pcr_cond, noval_pcr_cond, c_coord_wgs_x_land, c_coord_wgs_y_land)
                pixel_val_cs_point_pcr_stage_land = find_pix_val(rb_pcr_stage, gt_pcr_stage, rx_pcr_stage, ry_pcr_stage, noval_pcr_stage, c_coord_wgs_x_land, c_coord_wgs_y_land)    
                pixel_val_cs_point_pcr_rivbot_land = find_pix_val(rb_pcr_rivbot, gt_pcr_rivbot, rx_pcr_rivbot, ry_pcr_rivbot, noval_pcr_rivbot, c_coord_wgs_x_land, c_coord_wgs_y_land)  
                pixel_val_cs_point_gaia_width_land = find_pix_val(rb_gaia_width, gt_gaia_width, rx_gaia_width, ry_gaia_width, noval_gaia_width, c_coord_wgs_x_land, c_coord_wgs_y_land) 
                pixel_val_cs_point_gaia_depth_land = find_pix_val(rb_gaia_depth, gt_gaia_depth, rx_gaia_depth, ry_gaia_depth, noval_gaia_depth, c_coord_wgs_x_land, c_coord_wgs_y_land) 
                
                pixel_val_cs_point_flo_sea = find_pix_val(rb_flo, gt_flo, rx_flo, ry_flo, noval_flo, c_coord_wgs_x_sea, c_coord_wgs_y_sea)
                pixel_val_cs_point_cond_10_sea = find_pix_val(rb_cond_10, gt_cond_10, rx_cond_10, ry_cond_10, noval_cond_10, c_coord_wgs_x_sea, c_coord_wgs_y_sea)
                pixel_val_cs_point_cond_05_sea = find_pix_val(rb_cond_05, gt_cond_05, rx_cond_05, ry_cond_05, noval_cond_05, c_coord_wgs_x_sea, c_coord_wgs_y_sea)
                pixel_val_cs_point_cond_01_sea = find_pix_val(rb_cond_01, gt_cond_01, rx_cond_01, ry_cond_01, noval_cond_01, c_coord_wgs_x_sea, c_coord_wgs_y_sea)
                pixel_val_cs_point_pcr_cond_sea = find_pix_val(rb_pcr_cond, gt_pcr_cond, rx_pcr_cond, ry_pcr_cond, noval_pcr_cond, c_coord_wgs_x_sea, c_coord_wgs_y_sea)
                pixel_val_cs_point_pcr_stage_sea = find_pix_val(rb_pcr_stage, gt_pcr_stage, rx_pcr_stage, ry_pcr_stage, noval_pcr_stage, c_coord_wgs_x_sea, c_coord_wgs_y_sea)    
                pixel_val_cs_point_pcr_rivbot_sea = find_pix_val(rb_pcr_rivbot, gt_pcr_rivbot, rx_pcr_rivbot, ry_pcr_rivbot, noval_pcr_rivbot, c_coord_wgs_x_sea, c_coord_wgs_y_sea)  
                pixel_val_cs_point_gaia_width_sea = find_pix_val(rb_gaia_width, gt_gaia_width, rx_gaia_width, ry_gaia_width, noval_gaia_width, c_coord_wgs_x_sea, c_coord_wgs_y_sea) 
                pixel_val_cs_point_gaia_depth_sea = find_pix_val(rb_gaia_depth, gt_gaia_depth, rx_gaia_depth, ry_gaia_depth, noval_gaia_depth, c_coord_wgs_x_sea, c_coord_wgs_y_sea) 

                ##  define sum and count of avg_values for final average value
                sum_avg_land, cnt_avg_land = 0, 0
                sum_avg_sea, cnt_avg_sea = 0, 0

                ##  calculate the average value for each cross-section point
                for b in range(n_points_avg):
                    avg_point_dist = avg_dist - b * (max_dist/n_points)

                    avg_point_coord_land = avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[1], c_coords_utm[0], avg_point_dist)
                    avg_point_coord_sea = avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[3], c_coords_utm[2], avg_point_dist)

                    ##  transform to wgs coordinates
                    avg_point_coord_wgs_land = equator_position_angles(avg_point_coord_land[1], avg_point_coord_land[0], avg_point_coord_land[3], avg_point_coord_land[2], zone_num, zone_let)
                    avg_point_coord_wgs_sea = equator_position_angles(avg_point_coord_sea[1], avg_point_coord_sea[0], avg_point_coord_sea[3], avg_point_coord_sea[2], zone_num, zone_let)

                    ##  get the pixel value for the cross-section avg point for the land pixel
                    #avg_pixel_val_land_1 = find_pix_val(avg_point_coord_wgs_land[0][1], avg_point_coord_wgs_land[0][0])
                    avg_pixel_val_land_1 = find_pix_val(rb_flo, gt_flo, rx_flo, ry_flo, noval_flo, avg_point_coord_wgs_land[0][1], avg_point_coord_wgs_land[0][0])
                    sum_avg_land += avg_pixel_val_land_1

                    #avg_pixel_val_land_2 = find_pix_val(avg_point_coord_wgs_land[1][1], avg_point_coord_wgs_land[1][0])
                    avg_pixel_val_land_2 = find_pix_val(rb_flo, gt_flo, rx_flo, ry_flo, noval_flo, avg_point_coord_wgs_land[1][1], avg_point_coord_wgs_land[1][0])
                    sum_avg_land += avg_pixel_val_land_2
                    cnt_avg_land += 2

                    ##  get the pixel value for the cross-section avg point for the sea pixel
                    #avg_pixel_val_sea_1 = find_pix_val(avg_point_coord_wgs_sea[0][1], avg_point_coord_wgs_sea[0][0])
                    avg_pixel_val_sea_1 = find_pix_val(rb_flo, gt_flo, rx_flo, ry_flo, noval_flo, avg_point_coord_wgs_sea[0][1], avg_point_coord_wgs_sea[0][0])
                    sum_avg_sea += avg_pixel_val_sea_1

                    #avg_pixel_val_sea_2 = find_pix_val(avg_point_coord_wgs_sea[1][1], avg_point_coord_wgs_sea[1][0])
                    avg_pixel_val_sea_2 = find_pix_val(rb_flo, gt_flo, rx_flo, ry_flo, noval_flo, avg_point_coord_wgs_sea[1][1], avg_point_coord_wgs_sea[1][0]) 
                    sum_avg_sea += avg_pixel_val_sea_2
                    cnt_avg_sea += 2

                #pixel_val_land_AVG = sum_avg_land / cnt_avg_land
                #pixel_val_sea_AVG = sum_avg_sea / cnt_avg_sea
                pixel_val_cs_point_flo_land_AVG = sum_avg_land / cnt_avg_land
                pixel_val_cs_point_flo_sea_AVG = sum_avg_sea / cnt_avg_sea

                if cont_direction == 'right':
                    #lst_to_csv.append([cs_id, 'land', b_length, c_coord_wgs_x_land, c_coord_wgs_y_land, pixel_val_land, pixel_val_land_AVG])
                    #lst_to_csv.append([cs_id, 'sea', b_length, c_coord_wgs_x_sea, c_coord_wgs_y_sea, pixel_val_sea, pixel_val_sea_AVG])                   
                    lst_to_csv.append([cs_id, 'land', b_length, c_coord_wgs_x_land, c_coord_wgs_y_land, pixel_val_cs_point_flo_land, pixel_val_cs_point_flo_land_AVG,\
                                       pixel_val_cs_point_cond_01_land,\
                                       pixel_val_cs_point_cond_05_land, pixel_val_cs_point_cond_10_land, pixel_val_cs_point_pcr_cond_land, pixel_val_cs_point_pcr_stage_land,\
                                       pixel_val_cs_point_pcr_rivbot_land, pixel_val_cs_point_gaia_width_land, pixel_val_cs_point_gaia_depth_land])
                    lst_to_csv.append([cs_id, 'sea', b_length, c_coord_wgs_x_sea, c_coord_wgs_y_sea, pixel_val_cs_point_flo_sea, pixel_val_cs_point_flo_sea_AVG,\
                                       pixel_val_cs_point_cond_01_sea,\
                                       pixel_val_cs_point_cond_05_sea, pixel_val_cs_point_cond_10_sea, pixel_val_cs_point_pcr_cond_sea, pixel_val_cs_point_pcr_stage_sea,\
                                       pixel_val_cs_point_pcr_rivbot_sea, pixel_val_cs_point_gaia_width_sea, pixel_val_cs_point_gaia_depth_sea])

                else:
                    #lst_to_csv.append([cs_id, 'sea', b_length, c_coord_wgs_x_land, c_coord_wgs_y_land, pixel_val_land, pixel_val_land_AVG])
                    #lst_to_csv.append([cs_id, 'land', b_length, c_coord_wgs_x_sea, c_coord_wgs_y_sea, pixel_val_sea, pixel_val_sea_AVG])                   

                    lst_to_csv.append([cs_id, 'sea', b_length, c_coord_wgs_x_land, c_coord_wgs_y_land, pixel_val_cs_point_flo_land, pixel_val_cs_point_flo_land_AVG,\
                                       pixel_val_cs_point_cond_01_land,\
                                       pixel_val_cs_point_cond_05_land, pixel_val_cs_point_cond_10_land, pixel_val_cs_point_pcr_cond_land, pixel_val_cs_point_pcr_stage_land,\
                                       pixel_val_cs_point_pcr_rivbot_land, pixel_val_cs_point_gaia_width_land, pixel_val_cs_point_gaia_depth_land])
                    lst_to_csv.append([cs_id, 'land', b_length, c_coord_wgs_x_sea, c_coord_wgs_y_sea, pixel_val_cs_point_flo_sea, pixel_val_cs_point_flo_sea_AVG,\
                                       pixel_val_cs_point_cond_01_sea,\
                                       pixel_val_cs_point_cond_05_sea, pixel_val_cs_point_cond_10_sea, pixel_val_cs_point_pcr_cond_sea, pixel_val_cs_point_pcr_stage_sea,\
                                       pixel_val_cs_point_pcr_rivbot_sea, pixel_val_cs_point_gaia_width_sea, pixel_val_cs_point_gaia_depth_sea])

        ##  ZeroDivisionError happens when point is located on top of the node (closest node)
        ##  this is the case in small islands where equidistant points are created..
        except ZeroDivisionError:
            continue
            print('Division by zero for id ' + str(point_id + 1))

    return lst_to_csv






#   get the river discharge difference between the start of the model domain and the coast
#x_dst_cst = cs_plain_width * 1000.   
#riv_csv = os.path.join(id_riv_bas_dir, '_flo1k_vals.csv')
#cs_id_lst = id_cs_lst

def get_riv_Q(riv_csv, cs_id_lst, x_dst_cst):
    
    df = pd.read_csv(riv_csv)
    
    #   get all the riv_Q at end of the model domain
    riv_q_upstr = list(df.loc[(df['direction'] == 'land') & (df['dist_to_cst'] == x_dst_cst), 'flo1k_val_AVG'].values)

    #   loop through individual cs profiles, sometimes the 'cst' point has nan value because the raster and coastal point mismatch
    #   avg - average of all values that are != 0
    #   absolute max of the list
    #   last non 0.0 value - last pixel before ocean
    riv_q_downstr_avg, riv_q_downstr_max, riv_q_downstr_last_non_nan = [], [], []
    for cs_id in cs_id_lst:
        #   select all the land points where the distance is lower than x_dst_cst
        riv_q_lst = list(df.loc[(df['direction'] == 'land') & (df['dist_to_cst'] < x_dst_cst) & (df['id_cs'] == cs_id), 'flo1k_val_AVG'].values)
        non_nan_lst = [x for x in riv_q_lst if x != 0.0]
        riv_q_downstr_avg.append(sum(non_nan_lst) / len(non_nan_lst))
        riv_q_downstr_max.append(max(non_nan_lst))
        riv_q_downstr_last_non_nan.append(non_nan_lst[-1])

    print(sum(riv_q_upstr) / len(riv_q_upstr))                                 #   average upstream Q
    print(sum(riv_q_downstr_avg) / len(riv_q_downstr_avg))                     #   average downstream Q
    print(sum(riv_q_downstr_max) / len(riv_q_downstr_max))                     #   average MAXIMUM Q (from all cross-sections)
    print(max(riv_q_downstr_max))                                              #   absolute maximum Q
   
    return round(sum(riv_q_upstr) / len(riv_q_upstr), 2), round(sum(riv_q_downstr_avg) / len(riv_q_downstr_avg), 2),\
           round(sum(riv_q_downstr_max) / len(riv_q_downstr_max), 2)
    


"""
riv_csv = os.path.join(id_riv_bas_dir, '_flo1k_vals.csv')
cs_id_lst = id_cs_lst
x_dst_cst =  max(lst_cs_width) * 1000.
"""

#   get frequency of RIV cell (based on GAIA width dataset) in each cross-section cell across all individual cross-sections -> ARP
def get_riv_info(riv_csv, cs_id_lst, x_dst_cst):

    df = pd.read_csv(riv_csv)
    
    #   create an empty list where we will store the counts/occurrences of river cells from individual cross-sections
    riv_occ_lst, riv_all_width_lst, riv_stage_pcr_lst, riv_bot_elev_pcr_lst = [], [], [], []

    #   get all the river width (GAIA) from the coast till the end of the model domain, loop by cross-section
    for cs_id in cs_id_lst:
        #   select all the land points where the distance is lower than x_dst_cst
        riv_width_land = list(df.loc[(df['direction'] == 'land') & (df['dist_to_cst'] < x_dst_cst) & (df['id_cs'] == cs_id), 'R_width_GAIA'].values)
        #riv_width_cst = list(df.loc[(df['direction'] == 'cst') & (df['dist_to_cst'] <= x_dst_cst) & (df['id_cs'] == cs_id), 'R_width_GAIA'].values)
        #riv_width_all = riv_width_land + riv_width_cst
        riv_stage_pcr = list(df.loc[(df['direction'] == 'land') & (df['dist_to_cst'] < x_dst_cst) & (df['id_cs'] == cs_id), 'R_head_elev_PCR'].values)
        riv_bot_elev_pcr = list(df.loc[(df['direction'] == 'land') & (df['dist_to_cst'] < x_dst_cst) & (df['id_cs'] == cs_id), 'R_riv_bot_elev_PCR'].values)
        
        riv_width_all = riv_width_land
        riv_width_lst, riv_width_val_lst, riv_stage_lst, riv_bot_elev_lst = [], [], [], []
        
        #   loop through the list and assign 1 or 0 if there is or istn a river cell
        for a in range(len(riv_width_all)):
            if riv_width_all[a] != 0.0:
                riv_width_lst.append(1)
                riv_width_val_lst.append(riv_width_all[a])
                riv_stage_lst.append(riv_stage_pcr[a])
                riv_bot_elev_lst.append(riv_bot_elev_pcr[a])
            else:
                riv_width_lst.append(0)
                riv_width_val_lst.append(0.0)
                riv_stage_lst.append(0.0)
                riv_bot_elev_lst.append(0.0)
                
        #   check if the length of the list corresponds to the span of the whole average profile
        if len(riv_width_lst) / 2 < x_dst_cst / 1000.:
            #print(cs_id, len(riv_width_lst) / 2, x_dst_cst / 1000.)
            #for b in range(int(x_dst_cst / 1000. - len(riv_width_lst) / 2) * 2):
            for b in range(int(x_dst_cst / 1000.) * 2 - len(riv_width_lst)):                
                riv_width_lst.insert(0, 0)
                riv_width_val_lst.insert(0, 0.0)
                riv_stage_lst.insert(0, 0.0)
                riv_bot_elev_lst.insert(0, 0.0)

        #   append the final list to the riv_occ_lst 
        riv_occ_lst.append(riv_width_lst)
        riv_all_width_lst.append(riv_width_val_lst)
        riv_stage_pcr_lst.append(riv_stage_lst)
        riv_bot_elev_pcr_lst.append(riv_bot_elev_lst)

    #for i in range(len(riv_occ_lst)):
    #    print(len(riv_occ_lst[i]))

    #   sum all the sublists in the riv_occ_lst and then calculate the % chance of riv cell occurring in the given cell
    riv_arr = np.array(riv_occ_lst)
    riv_lst = np.sum(riv_arr, axis = 0).tolist()
    riv_lst_out = [round(i / len(cs_id_lst) * 100, 2) for i in riv_lst]
    
    #   do the same with the actual widths, stage and bottom elevation - create overall average and also the average of all non-zero widths
    riv_width_arr = np.array(riv_all_width_lst)
    riv_width_avg_lst = np.sum(riv_width_arr, axis = 0).tolist()
    riv_width_avg_lst_out = [round(i / len(cs_id_lst), 2) for i in riv_width_avg_lst]    
    
    riv_stage_arr = np.array(riv_stage_pcr_lst)
    riv_stage_avg_lst = np.sum(riv_stage_arr, axis = 0).tolist()
    riv_stage_avg_lst_out = [round(i / len(cs_id_lst), 2) for i in riv_stage_avg_lst]    
    
    riv_bot_elev_arr = np.array(riv_bot_elev_pcr_lst)
    riv_bot_elev_avg_lst = np.sum(riv_bot_elev_arr, axis = 0).tolist()
    riv_bot_elev_avg_lst_out = [round(i / len(cs_id_lst), 2) for i in riv_bot_elev_avg_lst]        
    
    riv_width_avg_nonzero_lst_out, riv_stage_avg_nonzero_lst_out, riv_bot_elev_avg_nonzero_lst_out = [], [], []
    
    for c in range(riv_width_arr.shape[1]):
        try:
            avg_width_non_zero = riv_width_avg_lst[c] / len([i for i in list(riv_width_arr[:, c]) if i != 0.0])
            avg_stage_non_zero = riv_stage_avg_lst[c] / len([i for i in list(riv_stage_arr[:, c]) if i != 0.0])
            avg_bot_elev_non_zero = riv_bot_elev_avg_lst[c] / len([i for i in list(riv_bot_elev_arr[:, c]) if i != 0.0])
        except ZeroDivisionError:
            avg_width_non_zero = 0.0
            avg_stage_non_zero = 0.0
            avg_bot_elev_non_zero = 0.0     
            
        riv_width_avg_nonzero_lst_out.append(avg_width_non_zero)
        riv_stage_avg_nonzero_lst_out.append(avg_stage_non_zero)    
        riv_bot_elev_avg_nonzero_lst_out.append(avg_bot_elev_non_zero)    

    riv_occ_n = []
    for g in range(riv_arr.shape[0]):
        riv_occ_n.append(np.count_nonzero(riv_arr[g, :] == 1))

    #   create a final dictionary
    dict_out = {}
    dict_out['riv_avg_occurence'] = riv_occ_n                                   #   how many times rivers cross the cross-section (based on gaia width raster.)
    dict_out['riv_loc_pct'] = riv_lst_out                                       #   % chance of river crossing in the given cell
    dict_out['riv_width_avg'] = riv_width_avg_lst_out                           
    dict_out['riv_width_avg_non_zero'] = riv_width_avg_nonzero_lst_out
    dict_out['riv_stage_avg'] = riv_stage_avg_lst_out
    dict_out['riv_stage_avg_non_zero'] = riv_stage_avg_nonzero_lst_out
    dict_out['riv_bot_elev_avg'] = riv_bot_elev_avg_lst_out
    dict_out['riv_bot_elev_avg_non_zero'] = riv_bot_elev_avg_nonzero_lst_out

    return dict_out#riv_lst_out, riv_width_avg_lst_out, riv_width_avg_nonzero_lst_out

    
#   function that builds the MODFLOW RIV package input file, mainly the stress_period dictionary
"""
cs_id = id_cs
riv_dict = riv_dict_input
x_dst_cst = abs(model.x_start) * 1000.
del_col = 100.0
rand_n = 1

1, riv_dict_input, abs(model.x_start) * 1000., 100.

"""
def create_RIV_input(rand_n, riv_dict, x_dst_cst, del_col):    
    
    #   get the average number of river crossings 
    riv_n = int(round(sum(riv_dict['riv_avg_occurence']) / len(riv_dict['riv_avg_occurence']), 0))
    
    #   based on the probability per cross-section point list randomly select the cross-section points where rivers will be assigned
    np.random.seed(rand_n)  

    cs_pt_indexes = np.arange(0, len(riv_dict['riv_loc_pct']), 1)
    p_func_norm = [round((i / round(sum(riv_dict['riv_loc_pct']), 2)), 2) for i in riv_dict['riv_loc_pct']]       #   normalized probability function
    rand_idx_lst = np.random.choice(cs_pt_indexes, riv_n, p_func_norm).tolist()
    #   remove duplicates
    rand_idx_lst = list(set(rand_idx_lst))
    
    #for i in rand_idx_lst:
    #    print(p_func_norm[i])

    #   loop through the randomly selected indexes and get the river widths at these locations
    riv_widths = [round(riv_dict['riv_width_avg_non_zero'][i], 0) for i in rand_idx_lst]

    #   relate the width to the width of the model column and create the final list of model indexes where the river will be placed
    riv_col_cells = []  #   [start_index, number of columns after]
    for h in range(len(riv_widths)):
        if riv_widths[h] != 0.0:
            #   always at least 1 column when the width is not 0
            riv_col_cells.append([rand_idx_lst[h], 1 + int(round(riv_widths[h] / del_col, 0)), riv_widths[h]])

    #   create output list
    lst_out = [] #  model_col index, riv_stage, riv_bot_elev, pct_Q_avg
    riv_id = 1
    for k in range(len(riv_col_cells)):
        cells_to_lst = np.arange(int(riv_col_cells[k][0] * (500. / del_col)), int(riv_col_cells[k][0] * (500. / del_col)) + riv_col_cells[k][1], 1).tolist()
        for l in range(len(cells_to_lst)):
            lst_out.append([riv_id, cells_to_lst[l], round(riv_dict['riv_stage_avg_non_zero'][riv_col_cells[k][0]], 2),\
                            round(riv_dict['riv_bot_elev_avg_non_zero'][riv_col_cells[k][0]], 2), round(((riv_col_cells[k][-1] / len(cells_to_lst)) / sum(riv_widths)) * 100., 2)])
        riv_id += 1
    print(lst_out)        
    return lst_out    
        
    
    """
    #   relate the width to the width of the model column and create the final list of model indexes where the river will be placed
    riv_col_cells = []  #   [start_index, number of columns after]
    for h in range(len(riv_widths)):        
        if (round(riv_widths[h] / del_col, 0)) == 0:
            pass
        else:
            riv_col_cells.append([rand_idx_lst[h], int(round(riv_widths[h] / del_col, 0))])
            #print(round(riv_widths[h] / del_col, 0))
        
        
    #   create output list
    lst_out = [] #  model_col index, riv_stage, riv_bot_elev
    for k in range(len(riv_col_cells)):
        cells_to_lst = np.arange(int(riv_col_cells[k][0] * (500. / del_col)), int(riv_col_cells[k][0] * (500. / del_col)) + riv_col_cells[k][1], 1).tolist()
        for l in range(len(cells_to_lst)):
            lst_out.append([cells_to_lst[l], riv_dict['riv_stage_avg_non_zero'][riv_col_cells[k][0]], riv_dict['riv_bot_elev_avg_non_zero'][riv_col_cells[k][0]]])
    """


    














"""
Function that reads in the NC file(s) for a coscat region and creates a dictionary with average concentration profile.

id_cs = coscat_id
cs_dir = r'g:\Water_Nexus\_A4\_COSCAT_avg_NC_files'
dir_out = r'g:\Water_Nexus\_A4\_COSCAT_input_files'

"""
def create_avg_COSCAT_profile(id_cs, cs_dir, dir_out):

    #   define the right string from the coscat id number
    if id_cs < 10:
        id_cs_str = '000' + str(id_cs)
    elif id_cs >= 10 and id_cs < 100:
        id_cs_str = '00' + str(id_cs)
    elif id_cs >= 100 and id_cs < 1000:
        id_cs_str = '0' + str(id_cs)
    else:
        id_cs_str = str(id_cs)  

    #   define the directories 
    cs_in_dir = os.path.join(cs_dir, id_cs_str)
    cs_out_dir = os.path.join(dir_out, id_cs_str)
    os.makedirs(cs_out_dir, exist_ok = True)

    #   find all the netcdf files in the right directory
    x_coord_lst, y_coord_lst = [], []
    nc_arr_conc, nc_arr_head = [], []
    
    for file in os.listdir(cs_in_dir):
        if file.endswith('.nc'):
            nc_file = xr.open_dataset(os.path.join(cs_in_dir, file))
            x_coord_lst.append([nc_file.coords['x'].values.tolist()[0], nc_file.coords['x'].values.tolist()[-1]])
            y_coord_lst.append([nc_file.coords['y'].values.tolist()[0], nc_file.coords['y'].values.tolist()[-1]])
            
            eq_21_time = max(nc_file.where(nc_file.eq == 'eq_21').dropna(dim = 'time')['time'].values.tolist())
            eq_21_conc = nc_file.where(nc_file.time == eq_21_time).dropna(dim = 'time')['solute concentration']
            eq_21_head = nc_file.where(nc_file.time == eq_21_time).dropna(dim = 'time')['heads']         
            
            nc_arr_conc.append(eq_21_conc)
            nc_arr_head.append(eq_21_head)
            
            #print(os.path.join(cs_in_dir, file))
            #print(nc_file.dims['x'], nc_file.dims['y'])
            #print(nc_file.coords['x'].values.tolist(), nc_file.dims['y'])
    
    #   choose the max, min for both x and y axis to capture the whole extent of all netcdf files (different coastal types merged into one)
    min_x = min([item[0] for item in x_coord_lst])
    max_x = max([item[1] for item in x_coord_lst])
    max_y = max([item[0] for item in y_coord_lst])
    min_y = min([item[1] for item in y_coord_lst])    
    
    avg_arr_conc = np.zeros((int((abs(max_y) + abs(min_y)) / 10) + 1, int((abs(max_x) + abs(min_x)) * 10) + 1))   #(rows, columns)
    avg_arr_head = np.zeros((int((abs(max_y) + abs(min_y)) / 10) + 1, int((abs(max_x) + abs(min_x)) * 10) + 1))   #(rows, columns)
    x_coords = np.linspace(min_x, max_x, avg_arr_conc.shape[1])
    y_coords = np.linspace(max_y, min_y, avg_arr_conc.shape[0])
    
    #   make sure the values are correctly rounded
    x_coords = [round(i, 2) for i in x_coords]
    y_coords = [round(i, 2) for i in y_coords]
    
    #   go cell by cell and assign the average non-nan value to the avg_arr
    for row in range(avg_arr_conc.shape[0]):
        print(row)
        for col in range(avg_arr_conc.shape[1]):
            x_coord_val = x_coords[col]
            y_coord_val = y_coords[row]
            conc_vals, head_vals = [], []

            #   get the values at the coordinates from the concentration and head arrays, no matter the dimensions
            for z in range(len(nc_arr_conc)):
                try:
                    conc_vals.append(nc_arr_conc[z].sel(x = x_coord_val).sel(y = y_coord_val).values[0])
                    head_vals.append(nc_arr_head[z].sel(x = x_coord_val).sel(y = y_coord_val).values[0])
                except KeyError:
                    conc_vals.append(1.00000002e+30)
                    head_vals.append(-999.9)                    

            #   get the non-nan values from the lists
            conc_val = [i for i in conc_vals if i < 100.]
            head_val = [i for i in head_vals if i > -999.]
            
            try:
                head_avg = sum(head_val) / len(head_val)
            except ZeroDivisionError:
                head_avg = -999.99

            try:
                conc_avg = sum(conc_val) / len(conc_val)
            except ZeroDivisionError:
                conc_avg = 999.99
            
            avg_arr_conc[row, col] = round(conc_avg, 2)
            avg_arr_head[row, col] = round(head_avg, 2)

    #   it can happen that there are nan values within the total model domain since multiple coastal type with different
    #   dimensions can be merged. To tackle that loop through the columns and inerpolate values if this occurs.
    
    #   first check if there is a column where the last row has a non-nan value and is followed by a column with nan values
    mid_col = [i for i in range(len(avg_arr_conc[-1, :])) if avg_arr_conc[-1, :][i] != 999.99][-1] 
    
    for col in range(mid_col):
        conc_vals = avg_arr_conc[:, col]
        head_vals = avg_arr_head[:, col]

        #   remove nan values from the start and end of the list - in such way extract a list of presumably only non-nan values
        for a in range(len(conc_vals)):
            if conc_vals[a] == 999.99:
                pass
            else:
                break
        
        for b in reversed(range(len(conc_vals))):
            if conc_vals[b] == 999.99:
                pass
            else:
                break        
            
        conc_nonnan_lst = conc_vals[a : b].tolist()
        head_nonnan_lst = head_vals[a : b].tolist()

        if b - a > 2:
            conc_nan_idx = [i + a for i in range(len(conc_nonnan_lst)) if conc_nonnan_lst[i] == 999.99] 
            
            if conc_nan_idx != []:
                if b != avg_arr_conc.shape[0] - 1:                
                    conc_new_vals = np.linspace(conc_nonnan_lst[conc_nan_idx[0]- a - 1], conc_nonnan_lst[conc_nan_idx[-1]- a + 1], len(conc_nan_idx))
                    head_new_vals = np.linspace(head_nonnan_lst[conc_nan_idx[0]- a - 1], head_nonnan_lst[conc_nan_idx[-1]- a + 1], len(conc_nan_idx))        
                else:                
                    conc_new_vals = np.linspace(conc_nonnan_lst[conc_nan_idx[0]- a - 1], conc_vals[-1], len(conc_nan_idx))
                    head_new_vals = np.linspace(head_nonnan_lst[conc_nan_idx[0]- a - 1], head_vals[-1], len(conc_nan_idx))    
                
                for c in range(len(conc_nan_idx)):
                    avg_arr_conc[conc_nan_idx[c], col] = conc_new_vals[c]
                    avg_arr_head[conc_nan_idx[c], col] = head_new_vals[c]
        
    if mid_col < avg_arr_conc.shape[1]:
        
        last_conc_val = avg_arr_conc[-1, mid_col]
        last_head_val = avg_arr_head[-1, mid_col]
        
        for col in range(mid_col, avg_arr_conc.shape[1]):
        
            conc_vals = avg_arr_conc[:, col]
            head_vals = avg_arr_head[:, col]
    
            #   remove nan values from the start and end of the list - in such way extract a list of presumably only non-nan values
            for a in range(len(conc_vals)):
                if conc_vals[a] == 999.99:
                    pass
                else:
                    break

            conc_nonnan_lst = conc_vals[a:].tolist()
            head_nonnan_lst = head_vals[a:].tolist()
    
            conc_nan_idx = [i + a for i in range(len(conc_nonnan_lst)) if conc_nonnan_lst[i] == 999.99] 
            
            if conc_nan_idx != []:            
                conc_new_vals = np.linspace(conc_nonnan_lst[conc_nan_idx[0]- a - 1], last_conc_val, len(conc_nan_idx))
                head_new_vals = np.linspace(head_nonnan_lst[conc_nan_idx[0]- a - 1], last_head_val, len(conc_nan_idx))        
                
                for c in range(len(conc_nan_idx)):
                    avg_arr_conc[conc_nan_idx[c], col] = conc_new_vals[c]
                    avg_arr_head[conc_nan_idx[c], col] = head_new_vals[c]        

    #   for the concentration, heads and cbc create a netcdf file (to save memory)
    nc_out = xr.Dataset(data_vars = {'solute concentration' : (('y', 'x'), avg_arr_conc),
                                     'heads' : (('y', 'x'), avg_arr_head)},
                           coords = {'x' : x_coords,
                                     'y' : y_coords})

    nc_out_name = '_COSCAT_' + id_cs_str + '_avg.nc'
    nc_out.to_netcdf(os.path.join(cs_out_dir, nc_out_name))  


    """
    #   plot the arrays depending on how many coastal types there are
    if len(nc_arr_conc) == 1:
        fig = plt.figure(figsize = (15, 8))
        ax1 = plt.subplot2grid((2, 2), (0, 1))  #   coastal type concentration profile
        ax2 = plt.subplot2grid((2, 2), (1, 1))  #   average concentration profile
        ax3 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)  #   color bar area 
        ax1.set_position([0.135, 0.55, 0.85, 0.4]) # [left, bottom, width, height]
        ax2.set_position([0.135, 0.1, 0.85, 0.4])
        ax3.set_position([0.04, 0.2, 0.025, 0.6])               
        
        #   define the extent of the figures
        #x_start_new = self.x_start #-(self.cst_idx - self.idx_start) / 2.
        #x_end_new = self.x_end
        #y_max_new = math.ceil((self.top / 100.0) * 100.0)
        #y_min_new = min(self.bot_elev)            
        
    elif len(nc_arr_conc) == 2:
        fig = plt.figure(figsize = (15, 12))
        
    elif len(nc_arr_conc) == 3:
        fig = plt.figure(figsize = (15, 16))
    """
    











"""
Function to load data from a raster file. Added cross-section (mostly deltas) have no data input for e.g. soil info and recharge.
This funciton goes point by point (0.5km cross-section points) and reads in the data. Returns a list.

in_raster_dir = r'g:\_ORIGINAL_DATA\water_table_depth\wtd_world.tif'
id_cs = 32921
n_points = 400
max_dist = 200000
#   define the number of points and distance from the cross-section point
n_points_avg = 5
avg_dist = 2500

#   connect to the databases
db_name = "'cs_db_deltas'"
db_name_world = "'cs_db_v1'"


db_cursor = dtbase_cur
db_cursor_write = dbase_connect.cur


tb_name = 'cs_pcrglob_thick'
num_type = 'real'
avg = False
"""
def read_cs_raster_data(in_raster_dir, id_cs, n_points, max_dist, n_points_avg, avg_dist, db_name, db_name_world, tb_name, num_type = 'real', avg = False):
    
    #   connect to databases
    db_host = "'localhost'"
    db_user = "'postgres'"
    db_pass = "'postgres'"
    #table_db = "cs_coastal_types"
    dbase_connect = ws.dbase_tools(db_name, db_user, db_host, db_pass)
    dbase_connect_world = ws.dbase_tools(db_name_world, db_user, db_host, db_pass)
    db_cur = dbase_connect.cur
    db_cur_world = dbase_connect_world.cur
    
    #   
    print('Reading raster infor for ID_CS = ' + str(id_cs) + ' ...... from ' + in_raster_dir)
    
    ##  first open the raster file via GDAL
    raster_in = gdal.Open(in_raster_dir)
    ##  assign the geotransform, band and noval
    gt = raster_in.GetGeoTransform()
    rb = raster_in.GetRasterBand(1)
    noval = rb.GetNoDataValue()
    rx = raster_in.RasterXSize
    ry = raster_in.RasterYSize

    ##   first run query to find all the 5km points
    sql_string_select_5km_point = "SELECT id_cs, ST_AsText((ST_Dump(%s.geom)).geom),\
                                    ST_X(geom), ST_Y(geom) AS the_POINT_geom FROM %s WHERE gid = %s"\
                                    % ('ne_10m_coastline_points_5km', 'ne_10m_coastline_points_5km', str(id_cs))
    ##  run the sql and get all the coastline points from the raw table
    db_cur_world.execute(sql_string_select_5km_point)
    point_5km = db_cur_world.fetchall()[0]

    ##  create table in database that will contain the extracted values from raster
    ##  dist_columns runs the function to get a list of columns for the created table
    dist_columns = create_column_string(2 * n_points, 500, num_type)
    sql_create_table(dbase_connect_world.conn, db_cur_world, tb_name, 'id_cs', dist_columns)

    ##  get id and coordinates
    point_id = point_5km[0]
    point_geom = point_5km[1]
    point_x, point_y = point_5km[2], point_5km[3]

    #   check if the point already exists in the database
    dtbase_cur0 = dbase_connect_world.cur
    sql_check_if_exists = "SELECT id_cs FROM %s WHERE id_cs = %s" % (tb_name, point_id)
    dtbase_cur0 = dbase_connect_world.cur
    dtbase_cur0.execute(sql_check_if_exists)
    exists = dtbase_cur0.fetchall()    

    if exists:
        print ('Data for cs_id = ' + str(point_id) + ' already exist, moving to next..')
        sql_select_id = "SELECT * FROM %s WHERE id_cs = %s" % (tb_name, id_cs)
        dbase_connect_world.cur.execute(sql_select_id)
        return dbase_connect_world.cur.fetchall()[0]
        
    else:
        print ('Extracting data for cs_id = ' + str(point_id))
        ##  create the sql string to find closest node to point_5km
        sql_string_select_points = "SELECT gid, ST_AsText(geom), ST_X(geom), ST_Y(geom) FROM ne_coastline_nodes_test ORDER BY ne_coastline_nodes_test.geom <-> ST_GeomFromText('%s') LIMIT 1;" % (point_geom)
        ##  open new database cursor
        dtbase_cur2 = dbase_connect_world.cur
        ##  run the sql and get results
        dtbase_cur2.execute(sql_string_select_points)
        closest_node = dtbase_cur2.fetchall()
        ##  assign the results of the query to variables below
        #closest_node_id = closest_node[0][0]
        #closest_node_geom = closest_node[0][1]
        closest_node_x = closest_node[0][2]
        closest_node_y = closest_node[0][3]
        ##  calculate the distance between the point_5m and closest node
        ##  to do this we convert the lat/lon coordinates to UTM coordinates (using the utm package)
        point_5km_utm = utm.from_latlon(point_y, point_x)
        closest_node_utm = utm.from_latlon(closest_node_y, closest_node_x)
        zone_num = point_5km_utm[2]
        zone_let = point_5km_utm[3]
        
        ##  check that zone_num is between 1 and 60
        if zone_num > 60:
            zone_num = zone_num - 60
        ##  formula for distance from point_5km to closest node
        #dist_to_node = math.sqrt((point_5km_utm[0] - closest_node_utm[0]) ** 2 + (point_5km_utm[1] - closest_node_utm[1]) ** 2)
        ##  extract value on the coastline point
        pixel_val_cs_point = read_raster_val(rb, gt, rx, ry, point_x, point_y)
    
        ##  insert the cs id to the table of the raster info
        sql_insert(dbase_connect_world.conn, dtbase_cur2, tb_name, 'id_cs', str(point_id))
        dtbase_cur3 = dbase_connect_world.cur
    
        ##  check if the point location falls on a novalue pixel of the raster
        ##  if yes then get the average/max freq value from neighbouring pixels
        if pixel_val_cs_point == noval:
            ##  get the raster coordinates
            px_cs_point = int((point_x - gt[0]) / gt[1]) #x pixel
            py_cs_point = int((point_y - gt[3]) / gt[5]) #y pixel
            ##  find the nearest value for the point
            near_val_cs_point = find_closest_neighbour(px_cs_point, py_cs_point, rb, noval, rx, ry, 'avg')
            ##  insert the value at coastline point to the database (dist_0m)
            sql_update(dbase_connect_world.conn, dtbase_cur3, tb_name, 'dist_0m', near_val_cs_point, 'id_cs = ' + str(point_id))
        ##  if not insert to database the extracted value
        else:
            sql_update(dbase_connect_world.conn, dtbase_cur3, tb_name, 'dist_0m', pixel_val_cs_point, 'id_cs = ' + str(point_id))
    
        ##  write the point into the database tables - cs and cs_gebco2014
        ##  check if the point is already written to cs table
        dtbase_cur4 = dbase_connect.cur
        sql_string_check_cs_id = "SELECT id_cs FROM cs WHERE id_cs = %s" % (str(point_id))
        dtbase_cur4.execute(sql_string_check_cs_id)
        id_check = dtbase_cur4.fetchall()
        if id_check:
            #print ('cs id : ' + str(point_id) + ' exists already..')
            pass
        else:
            sql_insert(dbase_connect.conn, dtbase_cur4, 'cs', 'id, x_coord, y_coord', str(point_id) + ',' + str(point_x) + ',' + str(point_y))            
        
        ##  calculate the coordinates of the point on the cross-section perpendicular to line
        ##  closest_node_point_5km and passing through point_5km
        for a in range(n_points):
            try:
                ##  calculate the b_length
                b_length = max_dist - a * (max_dist/n_points)
                c_coords_utm = coord_from_triangle(point_5km_utm[0], point_5km_utm[1], closest_node_utm[0], closest_node_utm[1], b_length)
    
                #print ('c_coords :  '), c_coords_utm
                ##  transform the coordinates back to wgs84, position of the equator is taken into account
                c_coords_wgs = equator_position_angles(c_coords_utm[1], c_coords_utm[0], c_coords_utm[3], c_coords_utm[2], zone_num, zone_let)
    
                c_coord_wgs_x_land, c_coord_wgs_y_land, c_coord_wgs_x_sea, c_coord_wgs_y_sea = \
                c_coords_wgs[0][1], c_coords_wgs[0][0], c_coords_wgs[1][1], c_coords_wgs[1][0]
    
                ##  read values from the GEBCO raster
                pixel_val_land = read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_land, c_coord_wgs_y_land)
                pixel_val_sea = read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_sea, c_coord_wgs_y_sea)
    
                dtbase_cur5 = dbase_connect_world.cur
    
                ##  write the values to the database table
                #   create string so it matches the column name in format dist_plus/minus_XXXXm
                col_name_land = 'dist_plus_' + str(int(b_length)) + 'm'
                col_name_sea = 'dist_minus_' + str(int(b_length)) + 'm'
                #   write the values to the database
                sql_update(dbase_connect_world.conn, dtbase_cur5, tb_name, col_name_land, pixel_val_land, 'id_cs = ' + str(point_id))
                sql_update(dbase_connect_world.conn, dtbase_cur5, tb_name, col_name_sea, pixel_val_sea, 'id_cs = ' + str(point_id))
    
            ##  ZeroDivisionError happens when point is located on top of the node (closest node)
            ##  this is the case in small islands where equidistant points are created..
            except ZeroDivisionError:
                continue
                print('Division by zero for id ' + str(point_id))        
        
        sql_select_id = "SELECT * FROM %s WHERE id_cs = %s" % (tb_name, id_cs)
        dbase_connect_world.cur.execute(sql_select_id)
        return dbase_connect_world.cur.fetchall()[0]       
        
        
    
    
def read_cs_raster_data_no_sql(in_raster_dir, id_cs, n_points, max_dist, n_points_avg, avg_dist, point_5km_tb, closest_node_x, closest_node_y, num_type = 'real', avg = False):
    
    print('Reading raster infor for ID_CS = ' + str(id_cs) + ' ...... from ' + in_raster_dir)
    
    ##  first open the raster file via GDAL
    raster_in = gdal.Open(in_raster_dir)
    ##  assign the geotransform, band and noval
    gt = raster_in.GetGeoTransform()
    rb = raster_in.GetRasterBand(1)
    noval = rb.GetNoDataValue()
    rx = raster_in.RasterXSize
    ry = raster_in.RasterYSize

    pts_info = point_5km_tb.loc[point_5km_tb['id_cs'] == id_cs]
    point_x, point_y = pts_info.x.values.tolist()[0], pts_info.y.values.tolist()[0]
    point_id = pts_info.id_cs.values.tolist()[0]

    print ('Extracting data for cs_id = ' + str(point_id))
    ##  create the sql string to find closest node to point_5km   
    ##  calculate the distance between the point_5m and closest node
    ##  to do this we convert the lat/lon coordinates to UTM coordinates (using the utm package)
    point_5km_utm = utm.from_latlon(point_y, point_x)
    closest_node_utm = utm.from_latlon(closest_node_y, closest_node_x)
    zone_num = point_5km_utm[2]
    zone_let = point_5km_utm[3]
    
    ##  check that zone_num is between 1 and 60
    if zone_num > 60:
        zone_num = zone_num - 60
    ##  formula for distance from point_5km to closest node
    #dist_to_node = math.sqrt((point_5km_utm[0] - closest_node_utm[0]) ** 2 + (point_5km_utm[1] - closest_node_utm[1]) ** 2)
    ##  extract value on the coastline point
    pixel_val_cs_point = read_raster_val(rb, gt, rx, ry, point_x, point_y)

    x_dist = np.arange(-1 * max_dist, max_dist + 500, 500).tolist()  
    val_arr = np.zeros(len(x_dist)).tolist()

    ##  check if the point location falls on a novalue pixel of the raster
    ##  if yes then get the average/max freq value from neighbouring pixels
    if pixel_val_cs_point == noval:
        ##  get the raster coordinates
        px_cs_point = int((point_x - gt[0]) / gt[1]) #x pixel
        py_cs_point = int((point_y - gt[3]) / gt[5]) #y pixel
        ##  find the nearest value for the point
        near_val_cs_point = find_closest_neighbour(px_cs_point, py_cs_point, rb, noval, rx, ry, 'avg')
        idx_0 = x_dist.index(0)
        val_arr[idx_0] = near_val_cs_point
    ##  if not insert to database the extracted value
    else:
        idx_0 = x_dist.index(0)
        val_arr[idx_0] = pixel_val_cs_point

    ##  calculate the coordinates of the point on the cross-section perpendicular to line
    ##  closest_node_point_5km and passing through point_5km
    for a in range(n_points):
        try:
            ##  calculate the b_length
            b_length = max_dist - a * (max_dist/n_points)
            c_coords_utm = coord_from_triangle(point_5km_utm[0], point_5km_utm[1], closest_node_utm[0], closest_node_utm[1], b_length)

            #print ('c_coords :  '), c_coords_utm
            ##  transform the coordinates back to wgs84, position of the equator is taken into account
            c_coords_wgs = equator_position_angles(c_coords_utm[1], c_coords_utm[0], c_coords_utm[3], c_coords_utm[2], zone_num, zone_let)

            c_coord_wgs_x_land, c_coord_wgs_y_land, c_coord_wgs_x_sea, c_coord_wgs_y_sea = \
            c_coords_wgs[0][1], c_coords_wgs[0][0], c_coords_wgs[1][1], c_coords_wgs[1][0]

            ##  read values from the GEBCO raster
            pixel_val_land = read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_land, c_coord_wgs_y_land)
            pixel_val_sea = read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_sea, c_coord_wgs_y_sea)

            ##  write the values to the database table
            #   create string so it matches the column name in format dist_plus/minus_XXXXm
            #   write the values to the database
            
            idx_land = x_dist.index(int(b_length))
            idx_sea = x_dist.index(int(-1 * b_length))
            
            val_arr[idx_land] = pixel_val_land
            val_arr[idx_sea] = pixel_val_sea

        ##  ZeroDivisionError happens when point is located on top of the node (closest node)
        ##  this is the case in small islands where equidistant points are created..
        except ZeroDivisionError:
            continue
            print('Division by zero for id ' + str(point_id))        
    
    return val_arr
        
    
    
"""
Function that reads in RIV input data and additionally provides cross-section points as output.
These will be later used as input for creating the polygon shapefile representing the given 
regional extent. 

cs_id = id_cs




"""


"""
id_coscat = 808
id_coscat_str = id_cs_str
id_subreg = 1
sc_topo = 'gebco_merit_avg'
topo_res = 100
in_dir = input_dir
subreg_dir = subreg_input_dir
topo_dir = topo_nc_dir
topo_500m_dir = topo_500m_nc_dir
gis_data_dir =  r'g:\Water_Nexus\_A4_GUM\_GIS_data'
cs_regs_shp_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\cs_reg_id.shp'
out_name = '_merit_avg_100m'
"""



#   function that saves the individual model inputs into a specific folder
def create_model_shape(id_coscat, id_subreg, id_coscat_str, in_dir, cs_regs_shp_dir, gis_data_dir, topo_dir, sc_topo, topo_res, topo_500m_dir, subreg_dir, out_name):

    #   define all necessary directories and database table names
    offshore_sed_dir = r'g:\_ORIGINAL_DATA\seabed_lithology\seabed_lithology_v1.nc'
    wtd_dir = r'g:\_ORIGINAL_DATA\water_table_depth\wtd_world_3.tif'       
    soil_type_dir = r'g:\_CREATED_DATA\soilgrids\merge_all_LZW.tif'
    soil_thk_dir = r'g:\_ORIGINAL_DATA\soilgrids\absolute_depth_to_bedrock_BDTICM_M_1km_ll.tif'
    glhymps_top_dir = r'g:\_CREATED_DATA\_A2_data\GLHYMPS\GLHYMPS_2_upper_layer_NODATA.tif'
    glhymps_bot_dir = r'g:\_CREATED_DATA\_A2_data\GLHYMPS\GLHYMPS_1_bottom_layer_compressed.tif'
    k_soil_dir = r'g:\_ORIGINAL_DATA\Hydraul_Param_SoilGrids_Schaap_0\ks_100cm_mm_d.tiff'
    drn_rate_dir = r'g:\_ORIGINAL_DATA\LCS_D_global_drainage\LCS_Dd_global.tif'
    

    """
    offshore_sed_dir = os.path.join(gis_data_dir, 'seabed_lithology', 'seabed_lithology_v1.nc') 
    wtd_dir = os.path.join(gis_data_dir, 'water_table_depth', 'wtd_world_3.tif')     
    soil_type_dir = os.path.join(gis_data_dir, 'soilgrids', 'merge_all_LZW.tif')   
    soil_thk_dir = os.path.join(gis_data_dir, 'soilgrids', 'absolute_depth_to_bedrock_BDTICM_M_1km_ll.tif')   
    glhymps_top_dir = os.path.join(gis_data_dir, 'GLHYMPS', 'GLHYMPS_2_upper_layer_NODATA.tif')  
    glhymps_bot_dir = os.path.join(gis_data_dir, 'GLHYMPS', 'GLHYMPS_1_bottom_layer_compressed.tif')  
    k_soil_dir = os.path.join(gis_data_dir, 'Hydraul_Param_SoilGrids_Schaap_0', 'ks_100cm_mm_d.tiff')  
    drn_rate_dir = os.path.join(gis_data_dir, 'LCS_D_global_drainage', 'LCS_Dd_global.tif')
    """
    
    df_points_5km = pd.read_csv(os.path.join(gis_data_dir, 'ne_10m_coastline_points_5km.csv'))
    #   open the nodes CSV file and transform it into geopandas dataframe
    df_nodes = pd.read_csv(os.path.join(gis_data_dir, 'ne_coastline_nodes_test.csv')) 
    geometry_nodes = [Point(xy) for xy in zip(df_nodes.x, df_nodes.y)]
    df_nodes = df_nodes.drop(['x', 'y'], axis = 1)
    nodes_geom = MultiPoint(geometry_nodes)
        
    #   first get the other necessary attributes :
    cs_regs_raw = gpd.read_file(cs_regs_shp_dir)
    id_cs_lst = cs_regs_raw.loc[(cs_regs_raw['coscat'] == id_coscat) & (cs_regs_raw['cst_reg_id'] == id_subreg)]['id_cs'].values.tolist()
        
    #       thickness at coast     
    #gem_csv = pd.read_csv(os.path.join(in_dir, id_coscat_str, '_SUBREG_representative_models_GEOMETRY_stats.csv'))
    #cst_thk_mu = gem_csv.loc[gem_csv[' id_SUBREG'] == id_subreg][' cs_thk_mu'].values[0]   
    #cst_thk = round(math.exp(cst_thk_mu), 0)
    
    #       anchor point distance and depth
    df_sed_thk_v2 = pd.read_csv(os.path.join(gis_data_dir, 'cs_sed_thick_est_v2.csv'))
    cs_width_lst, anch_dist_lst, anch_depth_lst, avg_depth_lst = [], [], [], []    
    for id_cs in id_cs_lst:
        cs_info = df_sed_thk_v2.loc[df_sed_thk_v2['id_cs'] == id_cs]
        try:
            cs_width_lst.append(cs_info.cst_plain_width.values[0])
            anchor_dist = cs_info.nasa_point_dist.values[0]
            if 200. - anchor_dist < 0:
                anchor_dist = 200.
            else:
                anchor_dist = 200. - anchor_dist
            anch_dist_lst.append(anchor_dist)
            anch_depth_lst.append(cs_info.nasa_point_depth.values[0])
            avg_depth_lst.append(cs_info.overall_avg.values[0])        
        except IndexError:
            continue    
    
    #   get the average topo lists
    avg_topo = xr.open_dataset(topo_dir)[sc_topo].values.tolist()
    avg_topo_gebco_500m = xr.open_dataset(topo_500m_dir)['gebco_avg'].values.tolist()
    x_nasa = round(sum(anch_dist_lst) / float(len(anch_dist_lst)), 1)
    end_cst_plain = round(sum(cs_width_lst) / float(len(cs_width_lst)), 1)
    cst_thick_est_avg = round(sum(avg_depth_lst) / float(len(avg_depth_lst)), 1)

    lst_wtd, lst_soil_type, lst_soil_thk   = [], [], []
    lst_cs_nasa, lst_cs_offshore, lst_cs_k_soil = [], [], []
    lst_cs_drn_rate, lst_cs_glhymps_top_lay, lst_cs_glhymps_bot_lay = [], [], []
    
    for id_cs in id_cs_lst:
        cs_model_point = ws_ms.point(None, None, 400, 200000, id_cs = id_cs)

        pts_info = df_points_5km.loc[df_points_5km['id_cs'] == id_cs]
        point_lat, point_lon = pts_info.x.values.tolist()[0], pts_info.y.values.tolist()[0]
        #   get the closest node to the coastal point
        nearest_geoms = nearest_points(Point(point_lat, point_lon), nodes_geom)
        near_idx1 = nearest_geoms[1]
        closest_node_lat, closest_node_lon = near_idx1.x, near_idx1.y   

        #   happens when there are no data extracted for the id_cs
        try:        
            cs_model_point.get_input_data(adjust_to_ocean = True)    
            cs_wtd = cs_model_point.dta_wtd_depth 
            cs_soil_type = cs_model_point.dta_soil_type
            cs_soil_thk = cs_model_point.dta_soil_thk
            cs_offshore = cs_model_point.dta_offshore_sed
            cs_offshore_flt = [float(i) for i in cs_offshore]
            cs_nasa = cs_model_point.dta_thk_pel
            cs_nasa_flt = [float(i) for i in cs_nasa]    
            cs_k_soil = cs_model_point.dta_k_soil
            cs_drn_rate = cs_model_point.dta_drn_rate
            cs_glhymps_top_lay = cs_model_point.dta_glhymps_top_lay                   
            cs_glhymps_bot_lay = cs_model_point.dta_glhymps_bot_lay  
    
            lst_wtd.append(cs_wtd)
            lst_soil_type.append(cs_soil_type)
            lst_soil_thk.append(cs_soil_thk)
            lst_cs_nasa.append(cs_nasa_flt)
            lst_cs_offshore.append(cs_offshore_flt)                
            lst_cs_k_soil.append(cs_k_soil)
            lst_cs_drn_rate.append(cs_drn_rate)  
            lst_cs_glhymps_top_lay.append(cs_glhymps_top_lay)
            lst_cs_glhymps_bot_lay.append(cs_glhymps_bot_lay)                
    
            cont_direction = cs_model_point.land_pos

        except IndexError:            
            try:
                df_cs_topo = pd.read_csv(os.path.join(gis_data_dir, '_csv_input_files', 'cs_gebco_2014_avg_DELTAS.csv'))
                cs_topo = df_cs_topo.loc[df_cs_topo['id_cs'] == id_cs].values.tolist()[0][400:1201]
                del df_cs_topo
                df_cs_nasa = pd.read_csv(os.path.join(gis_data_dir,'_csv_input_files', 'cs_nasa_thick_avg_DELTAS.csv'))
                cs_nasa = df_cs_nasa.loc[df_cs_nasa['id_cs'] == id_cs].values.tolist()[0][400:1201]
                del df_cs_nasa
                df_cs_pcr = pd.read_csv(os.path.join(gis_data_dir,'_csv_input_files', 'cs_aq_thick_DELTAS.csv'))
                cs_pcr = df_cs_pcr.loc[df_cs_pcr['id_cs'] == id_cs].values.tolist()[0][400:1201]
                del df_cs_pcr
                df_cs_glim = pd.read_csv(os.path.join(gis_data_dir,'_csv_input_files', 'cs_glim_litho_DELTAS.csv'))
                cs_glim = df_cs_glim.loc[df_cs_glim['id_cs'] == id_cs].values.tolist()[0][400:1201]
                cs_nasa_flt = cs_nasa
                del df_cs_glim            
                
                cs_offshore = read_cs_raster_data_no_sql(offshore_sed_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)
                cs_offshore_flt  = [float(i) if i is not None else -32768.0 for i in cs_offshore]  #[float(i) for i in cs_offshore]  
                cs_wtd = read_cs_raster_data_no_sql(wtd_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_soil_thk = read_cs_raster_data_no_sql(soil_thk_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_soil_type = read_cs_raster_data_no_sql(soil_type_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'int', avg = False)[1:]
                cs_glhymps_top_lay = read_cs_raster_data_no_sql(glhymps_top_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_glhymps_bot_lay = read_cs_raster_data_no_sql(glhymps_bot_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_k_soil = read_cs_raster_data_no_sql(k_soil_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]
                cs_drn_rate = read_cs_raster_data_no_sql(drn_rate_dir, id_cs, 400, 200000, 2, 2500, df_points_5km, closest_node_lat, closest_node_lon, num_type = 'real', avg = False)[1:]

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
                    cs_soil_thk = list(reversed(cs_soil_thk))
                    cs_soil_type = list(reversed(cs_soil_type))
                    cs_glhymps_top_lay = list(reversed(cs_glhymps_top_lay)) 
                    cs_glhymps_bot_lay = list(reversed(cs_glhymps_bot_lay))
                    cs_k_soil = list(reversed(cs_k_soil))
                    cs_drn_rate = list(reversed(cs_drn_rate))
            
                else:
                    continue     
                
                lst_cs_nasa.append(cs_nasa_flt)
                lst_wtd.append(cs_wtd)
                lst_soil_type.append(cs_soil_type)
                lst_soil_thk.append(cs_soil_thk)
                lst_cs_nasa.append(cs_nasa_flt)
                lst_cs_offshore.append(cs_offshore_flt)                
                lst_cs_k_soil.append(cs_k_soil)
                lst_cs_drn_rate.append(cs_drn_rate)  
                lst_cs_glhymps_top_lay.append(cs_glhymps_top_lay)
                lst_cs_glhymps_bot_lay.append(cs_glhymps_bot_lay)
                
            except IndexError:
                pass
            
    #   clean the lists - if the values are = non_val then change to NaN not to screw up the calculation of averages
    arr_wtd = np.array(lst_wtd)
    arr_wtd[arr_wtd <= -999.0] = np.nan     # sometimes there are values of -9999.0 might be an error in the dataset, otherwise nonval = -999
    arr_soil_type = np.array(lst_soil_type)
    arr_soil_thk = np.array(lst_soil_thk)
    arr_soil_thk[arr_soil_thk <= -9999.0] = np.nan  # -32768.0
    arr_nasa = np.array(lst_cs_nasa) 
    arr_nasa[arr_nasa == -1.0] = np.nan
    arr_offshore = np.array(lst_cs_offshore) # -32768
    arr_offshore[arr_offshore == -32768.0] = np.nan
    arr_k_soil = np.array(lst_cs_k_soil)
    arr_k_soil[arr_k_soil <= -9999.] = np.nan
    arr_drn_rate = np.array(lst_cs_drn_rate)   
    arr_drn_rate[arr_drn_rate <= -9999.] = np.nan
    arr_glhymps_top_lay = np.array(lst_cs_glhymps_top_lay)   
    arr_glhymps_bot_lay = np.array(lst_cs_glhymps_bot_lay)   
    arr_glhymps_top_lay[arr_glhymps_top_lay <= -9999.] = np.nan   
    arr_glhymps_bot_lay[arr_glhymps_bot_lay <= -9999.] = np.nan   
    avg_wtd = np.nanmean(arr_wtd, axis = 0)
    avg_nasa = np.nanmean(arr_nasa, axis = 0)
    avg_offshore = np.nanmean(arr_offshore, axis = 0)

    try:
        avg_k_soil = np.nanmean(arr_k_soil, axis = 0)
    except ZeroDivisionError:
        arr_k_soil_data = np.array(arr_k_soil, dtype='float')
        avg_k_soil = np.nanmean(arr_k_soil_data, axis = 0)
     
    avg_glhymps_top_lay = np.nanmean(arr_glhymps_top_lay, axis = 0)        
    avg_glhymps_bot_lay = np.nanmean(arr_glhymps_bot_lay, axis = 0)        

    avg_drn_rate = np.nanmean(arr_drn_rate, axis = 0)
    avg_soil_thk = np.nanmean(arr_soil_thk, axis = 0)
    avg_soil_type = []
    
    for i in range(len(avg_topo_gebco_500m)):
        if avg_topo_gebco_500m[i] >= 0:
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
            lst_wtd_topo.append(avg_topo_gebco_500m[a] - wtd_vals_round[a])
        else:
            lst_wtd_topo.append(avg_topo_gebco_500m[a])   
    arr_wtd_topo = np.array(lst_wtd_topo)        
    
    #   refine the lists based on the resolution of the topography list 
    x_lst = np.arange(-200., 200 + 500. / 1000., 500. / 1000.)  # 500. because that is the default spacing of cross-section points
    x_refined_lst = np.arange(-200., 200 + topo_res / 1000., topo_res / 1000.)    
    
    avg_wtd_topo_refined = np.interp(x_refined_lst, x_lst, arr_wtd_topo).tolist()    
    avg_wtd_topo_refined = [round(i, 4) for i in avg_wtd_topo_refined]
    avg_nasa_refined = np.interp(x_refined_lst, x_lst, avg_nasa).tolist()    
    avg_nasa_refined = [round(i, 4) for i in avg_nasa_refined]    
    avg_offshore_refined = np.interp(x_refined_lst, x_lst, avg_offshore).tolist()    
    avg_offshore_refined = [round(i, 4) for i in avg_offshore_refined]   
    avg_k_soil_refined = np.interp(x_refined_lst, x_lst, avg_k_soil).tolist()    
    avg_k_soil_refined = [round(i, 4) for i in avg_k_soil_refined]   
    avg_drn_rate_refined = np.interp(x_refined_lst, x_lst, avg_drn_rate).tolist()    
    avg_drn_rate_refined = [round(i, 4) for i in avg_drn_rate_refined]   
    avg_soil_thk_refined = np.interp(x_refined_lst, x_lst, avg_soil_thk).tolist()    
    avg_soil_thk_refined = [round(i, 4) for i in avg_soil_thk_refined]   
    avg_soil_type_refined = np.interp(x_refined_lst, x_lst, avg_soil_type).tolist()    
    avg_soil_type_refined = [round(i, 4) for i in avg_soil_type_refined]   
    avg_glhymps_top_lay_refined = np.interp(x_refined_lst, x_lst, avg_glhymps_top_lay).tolist()  
    avg_glhymps_top_lay_refined = [round(i, 4) for i in avg_glhymps_top_lay_refined]
    avg_glhymps_bot_lay_refined = np.interp(x_refined_lst, x_lst, avg_glhymps_bot_lay).tolist()  
    avg_glhymps_bot_lay_refined = [round(i, 4) for i in avg_glhymps_bot_lay_refined]

    #   cst_thick_est_stdev = round(cs_model_point.sed_thick_est_stdev, 0) 
    topo_dict = dict()
    topo_dict['ibound_arr'] = []
    topo_dict['top_elev'] = []
    topo_dict['bot_elev'] = []
    topo_dict['lay_elev'] = []
    topo_dict['cst_offset'] = []
    topo_dict['lay_idx'] = []
    topo_dict['col_idx'] = []
    topo_dict['end_idx'] = []
    topo_dict['top'] = []
    topo_dict['zbot'] = []
    topo_dict['col_width'] = []
    topo_dict['idx_start'] = []
    topo_dict['idx_end'] = []
    topo_dict['x_start'] = []
    topo_dict['x_end'] = []
    topo_dict['soil_thk'] = []
    topo_dict['soil_type'] = []
    topo_dict['cst_idx'] = []
    topo_dict['wtd_elev'] = []
    topo_dict['k_soil'] = []
    topo_dict['drn_rate'] = []     
    topo_dict['glhymps_top_lay'] = []
    topo_dict['glhymps_bot_lay'] = []

    in_lst, off_lst, final_topo = avg_topo, avg_topo, avg_topo
        
    #   first check the topographical profiles both offshore and inland - if offshore all > 0 then dont model
    #   same case is for inland all < 0
    off_pos = [i for i in off_lst if i < 0]
    in_pos = [i for i in in_lst if i > 0]    

    #   create an overall figure where all the topographical models will be stored
    fig_all = plt.figure(figsize = (20, 10))

    if off_pos == [] or len([x for x in off_pos if x > -5]) == len(off_pos):
        name = 'no_data'
        ax_all = fig_all.add_subplot(1, 1, 1)
        ax_all.text(0.15, 0.5, r'Cant build model, all topography offshore > 0')
        ax_all.set_title(name)
        print("Cannot find coastline - all average offshore elevation values are > 0m asl., id_coscat = ") 

    elif in_pos == []:
        name = 'no_data'
        ax_all = fig_all.add_subplot(1, 1, 1)
        ax_all.text(0.15, 0.5, r'Cant build model, all topography offshore > 0')
        ax_all.set_title(name)
        print("Cannot find coastline - all average inland elevation values are < 0m asl., id_coscat = ") 
    
    else:
        #   define the index from which we will look for the coastline, just place it at the middle of the coastal profile 
        cst_look_up_idx =  int(((len(avg_topo) - 1) / 2))
        #   find the position of the coastline
        cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < 0.0), None)
        #   check if the coastal offset exists, if not it is probably beacuse the whole profile is above sea level
        if cst_offset_val is None:
            name = 'no_data'
            ax_all = fig_all.add_subplot(1, 1, 1)
            ax_all.text(0.15, 0.5, r'Cant build model, all topography offshore > 0')
            ax_all.set_title(name)
            print("Cannot find coastline - all average elevation values are > 0m asl.")            
        
        else:
            print("Pica") 
            cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
            #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
            cst_offset = ((cst_look_up_idx + cst_offset_idx) - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000.
            cst_off_idx = cst_look_up_idx + cst_offset_idx
    
            if cst_off_idx < cst_look_up_idx:
                print("Pica 2")                 
                #   define the index from which we will look for the coastline 
                cst_look_up_idx = int((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) + int(((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) / 2)
                #   find the position of the coastline
                cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < 0.0), None)
                cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                cst_offset = ((cst_look_up_idx + cst_offset_idx) - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000.
                cst_off_idx = cst_look_up_idx + cst_offset_idx
                
                if cst_off_idx < cst_look_up_idx:
                    print("Pica 3")                         
                    
                    #   define the index from which we will look for the coastline 
                    cst_look_up_idx = int((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) + int(((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) / 4) * 3
                    #   find the position of the coastline
                    cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < 0.0), None)
                    cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                    #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                    cst_offset = ((cst_look_up_idx + cst_offset_idx) - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000. 
                    cst_off_idx = cst_look_up_idx + cst_offset_idx                
                    last_ibound_col_n_layers = -3
                    #return None, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist    

                    if cst_off_idx < cst_look_up_idx:
                        print("Pica 3a")           
    
                        #   define the index from which we will look for the coastline 
                        cst_look_up_idx = int((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) + int(((len(avg_topo) - 1) * (500. / topo_res) / 10.  / 2) / 8) * 7
                        #   find the position of the coastline
                        cst_offset_val  = next((x for x in final_topo[cst_look_up_idx:] if x < 0.0), None)
                        cst_offset_idx = final_topo[cst_look_up_idx:].index(cst_offset_val) - 1
                        #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                        cst_offset = ((cst_look_up_idx + cst_offset_idx) - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000.
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
                anchor_topo_elev = final_topo[int(round((cst_off_idx - x_nasa * (1000 / topo_res))))]
                y_nasa = anchor_topo_elev - 50.0    
            
                try: 
                    #   first try to find the shelf break and foot of continental slope points
                    print("Pica 5")    
                    offshore_dir = os.path.join(subreg_dir, '_topo_figs')
                    if not os.path.exists(offshore_dir):
                        os.makedirs(offshore_dir)
                    fos = find_FOS(final_topo, offshore_dir, 'avg_avg', avg = False) 
                    print(fos) 
                    
                    model = cs_model(1, final_topo, avg_nasa_refined, x_nasa, y_nasa, cst_thick_est_avg, end_cst_plain,\
                                     avg_soil_thk_refined, avg_soil_type_refined, avg_wtd_topo_refined, avg_offshore_refined,\
                                     final_topo, final_topo, final_topo, avg_k_soil_refined,\
                                     avg_drn_rate_refined, avg_glhymps_top_lay_refined, avg_glhymps_bot_lay_refined, topo_res, False)
                    
                    model.get_top_bot_lst_find_cst(cst_look_up_idx, fos[1], 0.0)
                    #model.get_top_bot_col_lst(del_col, cs_points_dist, model.cst_idx, smooth = True)
                    #model.get_top_bot_col_lst_top_sys_geo(del_col, cs_points_dist, model.cst_idx, smooth = True)
                    model.get_top_bot_col_lst(del_col, topo_res, model.cst_idx, smooth = True)
                    model.get_top_bot_col_lst_top_sys_geo(del_col, topo_res, model.cst_idx, smooth = True)                    
                    model.modflow_bot_elev_list(del_lay, True)

                except AttributeError:
                    print("Pica 6")   
                    name = 'no_data'
                    ax_all = fig_all.add_subplot(1, 1, 1)
                    ax_all.text(0.15, 0.5, r'Cant build model, all topography offshore > 0')
                    ax_all.set_title(name)                    
                    #return None, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist

                try:
                    print("Pica 7")  
                    model.modflow_col_list_constant(cs_points_dist, del_col)
                    model.modflow_create_IBOUND_arr()                
                    
                    try:
                        model.find_topo_divide_SRM_v2(10, -2000., del_col, del_lay, fos[1], 100., 20., 50000., 100000., 10000, 50000)

                        #model.find_topo_divide(10., -10000., del_col, model.top_elev[model.idx_start], fos[1], 2, 200000.)         #   indicate the value of topographic difference to find the divide 
                    #   in case the island is small cut off the inland part with the start of the model.top_elev list
                    except IndexError:
                        #model.find_topo_divide(10., -10000., del_col, model.top_elev[0], fos[1], 2, 200000.)  
                        model.find_topo_divide_SRM(10, -2000., del_col, del_lay, fos[1], 20., 50000., 100000., 10000, 50000)

                    model.dis_input(nrow, delc, del_col, del_lay, nper, perlen, nstp, laycbd, max_depth, cs_points_dist)
                    model.top = model.lay_elev[model.lay_idx]            #   top elevation of the model domain - check if it should be like this!
                    model.botm = model.lay_elev[model.lay_idx + 1:]   
                    #model.zbot = model.botm[-1] 
                    zbot = math.floor((model.bot_elev[-1] / 100.0) * 100.0)
                                 
                    ibound_arr = model.ibound_arr * 1.0
                    ibound_arr[np.abs(ibound_arr) == 0.] = np.nan
                    
                    str_0 = '\\' + out_name + '.png'
                    plot_dir = os.path.join(subreg_dir, str_0)           
                    
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

                    top_elev_arr = np.array(model.top_elev)
                    mask = np.isnan(top_elev_arr)   
                    top_elev_arr[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), top_elev_arr[~mask])
                    
                    model.top_elev = top_elev_arr.tolist()
    
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
                    ax_all.set_title(out_name)
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
                    topo_dict[out_name] = {'ibound_arr': model.ibound_arr,\
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
                                          'k_soil': model.k_soil,\
                                          'drn_rate': model.drn_dens,\
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
    plot_name_all = os.path.join(subreg_dir, '__' + out_name + '_IBOUNDs.png')
    plt.savefig(plot_name_all, dpi = 300, facecolor = 'w', edgecolor = 'w',
        orientation = 'portrait', papertype = None, format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
        frameon = None)
    plt.close()

    #   save the dictionary with the ibound arrays
    dict_save_dir = os.path.join(subreg_dir, out_name + '_IBOUND.npy')
    np.save(dict_save_dir, topo_dict)

    return topo_dict, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist, fos    
    
    




#   function that saves the individual model inputs into a specific folder
def save_model_input_files_cst_type_topo_v2(cst_type, thk_lst, width_lst, anchor_dist_lst, anchor_depth_lst, wtd_lst, topo_lst, topo_gebco_500m_lst,\
                                         soil_type_lst, soil_thk_lst, nasa_lst, offshore_lst, pcr_lst, watergap_lst, p_min_et_lst,\
                                         k_soil_lst, drn_rate_lst,\
                                         glhymps_top_lay_lst, glhymps_bot_lay_lst, topo_res, del_col, del_lay, id_coscat, out_dir, out_name):

    """
    thk_lst = lst_thk
    width_lst = lst_cs_width
    anchor_dist_lst = lst_anchor_dist
    anchor_depth_lst = lst_anchor_depth
    wtd_lst = lst_wtd
    topo_lst = lst_gebco_merit_avg_100m
    topo_gebco_500m_lst = lst_gebco_merit_avg_500m
    soil_type_lst = soil_type_lst
    soil_thk_lst = lst_soil_thk
    nasa_lst = lst_cs_nasa
    offshore_lst = lst_cs_offshore 
    pcr_lst = lst_cs_pcr_rch
    watergap_lst = lst_cs_watergap_rch
    p_min_et_lst = lst_cs_p_min_et
    k_soil_lst = lst_cs_k_soil
    drn_rate_lst = lst_cs_drn_rate
    glhymps_top_lay_lst = lst_cs_glhymps_top_lay
    glhymps_bot_lay_lst = lst_cs_glhymps_bot_lay
    id_coscat = id_cs
    out_dir = id_riv_bas_dir      
    topo_diff = 10.
    bot_limit = -10000.
    ibound_act_cells_limit = 2
    offshore_dist_limit = 200000.
    topo_res = 100.
    del_lay = 10.
    del_col = 100.
    out_name = '_gebco_merit_avg_100m'
 
    
    new_ibound_dict = rivbas.save_model_input_files_cst_type_topo('', thk_lst, width_lst, anchor_dist_lst, anchor_depth_lst, wtd_lst, topo_lst, topo_gebco_500m_lst,\
                                                              soil_type_lst, soil_thk_lst, nasa_lst, offshore_lst, pcr_lst, watergap_lst, p_min_et_lst, k_soil_lst, drn_rate_lst,\
                                                              glhymps_top_lay_lst, glhymps_bot_lay_lst, 100., 100., 10., coscat_id, subreg_input_dir, '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND')

    
    
    
    id_coscat = coscat_id
    out_dir = subreg_input_dir
    out_name = '_' + topo_100_type + '_' + str(int(del_col_res)) + '_IBOUND'
    
    for i in range(len(lst_gebco_merit_avg_100m[0])):
        if math.isnan(lst_gebco_merit_avg_100m[0][i]):
            print(i)
    
    for i in range(len(lst_gebco_merit_avg_500m[0])):
        if math.isnan(lst_gebco_merit_avg_500m[0][i]):
            print(i)
        
    
    
    """
    sea_level = 0.0

    last_ibound_col_n_layers = -55
    ls_elev = -55
    ls_elev_50 = -55   
    fos_pt_dist = -55
    
    #   clean the lists - if the values are = non_val then change to NaN not to screw up the calculation of averages
    arr_wtd = np.array(wtd_lst)
    arr_wtd[arr_wtd <= -999.0] = np.nan     # sometimes there are values of -9999.0 might be an error in the dataset, otherwise nonval = -999
    arr_topo = np.array(topo_lst)
    arr_topo_gebco_500m = np.array(topo_gebco_500m_lst)
    arr_topo[arr_topo == -9999.0] = np.nan
    arr_topo_gebco_500m[arr_topo_gebco_500m == -9999.0] = np.nan
    arr_soil_type = np.array(soil_type_lst)
    arr_soil_thk = np.array(soil_thk_lst)
    arr_soil_thk[arr_soil_thk <= -9999.0] = np.nan  # -32768.0
    arr_nasa = np.array(nasa_lst) 
    arr_nasa[arr_nasa == -1.0] = np.nan
    arr_offshore = np.array(offshore_lst) # -32768
    arr_offshore[arr_offshore == -32768.0] = np.nan
    
    arr_k_soil = np.array(k_soil_lst)
    arr_k_soil[arr_k_soil <= -9999.] = np.nan
    arr_drn_rate = np.array(drn_rate_lst)   
    arr_drn_rate[arr_drn_rate <= -9999.] = np.nan
    
    arr_glhymps_top_lay = np.array(glhymps_top_lay_lst)   
    arr_glhymps_bot_lay = np.array(glhymps_bot_lay_lst)   
    arr_glhymps_top_lay[arr_glhymps_top_lay <= -9999.] = np.nan   
    arr_glhymps_bot_lay[arr_glhymps_bot_lay <= -9999.] = np.nan   
    
    avg_wtd = np.nanmean(arr_wtd, axis = 0)
    if len(arr_topo.shape) == 1:
        avg_topo = arr_topo
    else:
        avg_topo = np.nanmean(arr_topo, axis = 0)
    avg_topo_gebco_500m = np.nanmean(arr_topo_gebco_500m, axis = 0)    
    avg_nasa = np.nanmean(arr_nasa, axis = 0)
    avg_offshore = np.nanmean(arr_offshore, axis = 0)

    try:
        avg_k_soil = np.nanmean(arr_k_soil, axis = 0)
    except ZeroDivisionError:
        arr_k_soil_data = np.array(arr_k_soil, dtype='float')
        avg_k_soil = np.nanmean(arr_k_soil_data, axis = 0)
     
    avg_glhymps_top_lay = np.nanmean(arr_glhymps_top_lay, axis = 0)        
    avg_glhymps_bot_lay = np.nanmean(arr_glhymps_bot_lay, axis = 0)        
        
    avg_drn_rate = np.nanmean(arr_drn_rate, axis = 0)
    avg_soil_thk = np.nanmean(arr_soil_thk, axis = 0)
    avg_soil_type = []
    
    for i in range(avg_topo_gebco_500m.shape[0]):
        if avg_topo_gebco_500m[i] >= 0:
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
            lst_wtd_topo.append(avg_topo_gebco_500m[a] - wtd_vals_round[a])
        else:
            lst_wtd_topo.append(avg_topo_gebco_500m[a])   
    arr_wtd_topo = np.array(lst_wtd_topo)        
    #avg_wtd_topo = np.nanmean(arr_wtd_topo, axis = 0)
    
    #   refine the lists based on the resolution of the topography list 
    x_lst = np.arange(-200., 200 + 500. / 1000., 500. / 1000.)  # 500. because that is the default spacing of cross-section points
    x_refined_lst = np.arange(-200., 200 + topo_res / 1000., topo_res / 1000.)    
    
    try:
        avg_wtd_topo_refined = np.interp(x_refined_lst, x_lst, arr_wtd_topo).tolist()    
    except ValueError:
        if arr_wtd_topo.shape[0] < len(x_lst):
            for z in range(len(x_lst) - arr_wtd_topo.shape[0]):
                arr_wtd_topo = np.append(arr_wtd_topo, arr_wtd_topo[-1])
            avg_wtd_topo_refined = np.interp(x_refined_lst, x_lst, arr_wtd_topo).tolist()            
    avg_wtd_topo_refined = [round(i, 4) for i in avg_wtd_topo_refined]

    try:
        avg_nasa_refined = np.interp(x_refined_lst, x_lst, avg_nasa).tolist()    
    except ValueError:
        if avg_nasa.shape[0] < len(x_lst):
            for z in range(len(x_lst) - avg_nasa.shape[0]):
                avg_nasa = np.append(avg_nasa, avg_nasa[-1])
            avg_nasa_refined = np.interp(x_refined_lst, x_lst, avg_nasa).tolist()            
    avg_nasa_refined = [round(i, 4) for i in avg_nasa_refined]
    
    try:
        avg_offshore_refined = np.interp(x_refined_lst, x_lst, avg_offshore).tolist()    
    except ValueError:
        if avg_offshore.shape[0] < len(x_lst):
            for z in range(len(x_lst) - avg_offshore.shape[0]):
                avg_offshore = np.append(avg_offshore, avg_offshore[-1])
            avg_offshore_refined = np.interp(x_refined_lst, x_lst, avg_offshore).tolist()            
    avg_offshore_refined = [round(i, 4) for i in avg_offshore_refined]    
    
    try:
        avg_k_soil_refined = np.interp(x_refined_lst, x_lst, avg_k_soil).tolist()    
    except ValueError:
        if avg_k_soil.shape[0] < len(x_lst):
            for z in range(len(x_lst) - avg_k_soil.shape[0]):
                avg_k_soil = np.append(avg_k_soil, avg_k_soil[-1])
            avg_k_soil_refined = np.interp(x_refined_lst, x_lst, avg_k_soil).tolist()            
    avg_k_soil_refined = [round(i, 4) for i in avg_k_soil_refined]        
    
    try:
        avg_drn_rate_refined = np.interp(x_refined_lst, x_lst, avg_drn_rate).tolist()    
    except ValueError:
        if avg_drn_rate.shape[0] < len(x_lst):
            for z in range(len(x_lst) - avg_drn_rate.shape[0]):
                avg_drn_rate = np.append(avg_drn_rate, avg_drn_rate[-1])
            avg_drn_rate_refined = np.interp(x_refined_lst, x_lst, avg_drn_rate).tolist()            
    avg_drn_rate_refined = [round(i, 4) for i in avg_drn_rate_refined]     

    try:
        avg_soil_thk_refined = np.interp(x_refined_lst, x_lst, avg_soil_thk).tolist()    
    except ValueError:
        if avg_soil_thk.shape[0] < len(x_lst):
            for z in range(len(x_lst) - avg_soil_thk.shape[0]):
                avg_soil_thk = np.append(avg_soil_thk, avg_soil_thk[-1])
            avg_soil_thk_refined = np.interp(x_refined_lst, x_lst, avg_soil_thk).tolist()            
    avg_soil_thk_refined = [round(i, 4) for i in avg_soil_thk_refined]    
    
    try:
        avg_soil_type_refined = np.interp(x_refined_lst, x_lst, avg_soil_type).tolist()    
    except ValueError:
        if avg_soil_type.shape[0] < len(x_lst):
            for z in range(len(x_lst) - avg_soil_type.shape[0]):
                avg_soil_type = np.append(avg_soil_type, avg_soil_type[-1])
            avg_soil_type_refined = np.interp(x_refined_lst, x_lst, avg_soil_type).tolist()            
    avg_soil_type_refined = [round(i, 4) for i in avg_soil_type_refined]    
    
    try:
        avg_glhymps_top_lay_refined = np.interp(x_refined_lst, x_lst, avg_glhymps_top_lay).tolist()    
    except ValueError:
        if avg_glhymps_top_lay.shape[0] < len(x_lst):
            for z in range(len(x_lst) - avg_glhymps_top_lay.shape[0]):
                avg_glhymps_top_lay = np.append(avg_glhymps_top_lay, avg_glhymps_top_lay[-1])
            avg_glhymps_top_lay_refined = np.interp(x_refined_lst, x_lst, avg_glhymps_top_lay).tolist()            
    avg_glhymps_top_lay_refined = [round(i, 4) for i in avg_glhymps_top_lay_refined]    

    try:
        avg_glhymps_bot_lay_refined = np.interp(x_refined_lst, x_lst, avg_glhymps_bot_lay).tolist()    
    except ValueError:
        if avg_glhymps_bot_lay.shape[0] < len(x_lst):
            for z in range(len(x_lst) - avg_glhymps_bot_lay.shape[0]):
                avg_glhymps_bot_lay = np.append(avg_glhymps_bot_lay, avg_glhymps_bot_lay[-1])
            avg_glhymps_bot_lay_refined = np.interp(x_refined_lst, x_lst, avg_glhymps_bot_lay).tolist()            
    avg_glhymps_bot_lay_refined = [round(i, 4) for i in avg_glhymps_bot_lay_refined]    

    end_cst_plain = round(sum(width_lst) / float(len(width_lst)), 1)
    x_nasa = round(sum(anchor_dist_lst) / float(len(anchor_dist_lst)), 1)
    
    #   added this, if the x_nasa is outside of the coastal plain then 
    if x_nasa > end_cst_plain:
        x_nasa = end_cst_plain
    
    #y_nasa = round(sum(anchor_depth_lst) / float(len(anchor_depth_lst)), 1)
    cst_thick_est_avg = round(sum(thk_lst) / float(len(thk_lst)), 1)
        
    #   cst_thick_est_stdev = round(cs_model_point.sed_thick_est_stdev, 0) 
    topo_dict = dict()
    topo_dict['ibound_arr'] = []
    topo_dict['top_elev'] = []
    topo_dict['bot_elev'] = []
    topo_dict['lay_elev'] = []
    topo_dict['cst_offset'] = []
    topo_dict['lay_idx'] = []
    topo_dict['col_idx'] = []
    topo_dict['end_idx'] = []
    topo_dict['top'] = []
    topo_dict['zbot'] = []
    topo_dict['col_width'] = []
    topo_dict['idx_start'] = []
    topo_dict['idx_end'] = []
    topo_dict['x_start'] = []
    topo_dict['x_end'] = []
    topo_dict['soil_thk'] = []
    topo_dict['soil_type'] = []
    topo_dict['cst_idx'] = []
    topo_dict['wtd_elev'] = []
    topo_dict['k_soil'] = []
    topo_dict['drn_rate'] = []     
    topo_dict['glhymps_top_lay'] = []
    topo_dict['glhymps_bot_lay'] = []
            
    in_lst = avg_topo.tolist()
    off_lst = avg_topo.tolist()
    final_topo = avg_topo.tolist()
        
    #   first check the topographical profiles both offshore and inland - if offshore all > 0 then dont model
    #   same case is for inland all < 0
    off_pos = [i for i in off_lst if i < 0]
    in_pos = [i for i in in_lst if i > 0]    

    #   create an overall figure where all the topographical models will be stored
    fig_all = plt.figure(figsize = (20, 10))

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

    #   look for coastlin in the vicinity of the assumed coastline (middle of the cross-section)
    else:
        idx_range = 100
        left_lst, right_lst = [], []
        #   loop through the list that is cut out around the assumed coastline - look in the vicinity for the submerged value
        for h in range(0, idx_range):
            #   get first five values in each direction
            left_vals = final_topo[int((len(avg_topo) - 1) / 2) - h - 5 : int((len(avg_topo) - 1) / 2) - h]
            right_vals = final_topo[int((len(avg_topo) - 1) / 2) + h : int((len(avg_topo) - 1) / 2) + 5 + h]
            #   check the number of negative values in the lists
            neg_left = [i for i in left_vals if i < sea_level]
            neg_right = [i for i in right_vals if i < sea_level]
            #   append to the lists so we know how many negative values are there with increasing distance from the coast
            left_lst.append(len(neg_left))
            right_lst.append(len(neg_right))
        #   now loop through the overall lists and find the first value 5 - indicating that all the 5 values are below sea level
        left_coast_val  = next((x[0] for x in enumerate(left_lst) if x[1] == 5), None)
        right_coast_val  = next((x[0] for x in enumerate(right_lst) if x[1] == 5), None)
        #   now use these indexes (if they exist) to find the coastline offset value
        if left_coast_val is not None and right_coast_val is None:
            cst_offset_val  = final_topo[int((len(avg_topo) - 1) / 2) - left_coast_val]
        elif left_coast_val is None and right_coast_val is not None:
            cst_offset_val  = final_topo[int((len(avg_topo) - 1) / 2) + right_coast_val]
        #   if neither is None
        elif left_coast_val is not None and right_coast_val is not None:
            #   first check the position of first 0 on the left side (presumably the inland)
            left_coast_val  = next((x[0] for x in enumerate(left_lst) if x[1] == 0), None)
            cst_offset_val  = final_topo[int((len(avg_topo) - 1) / 2) - left_coast_val]
            
        #   if they dont exist then expand the search further from the coastline
        else:
            idx_range = 500
            left_lst, right_lst = [], []
            #   loop through the list that is cut out around the assumed coastline - look in the vicinity for the submerged value
            for h in range(0, idx_range):
                #   get first five values in each direction
                left_vals = final_topo[int((len(avg_topo) - 1) / 2) - h - 5 : int((len(avg_topo) - 1) / 2) - h]
                right_vals = final_topo[int((len(avg_topo) - 1) / 2) + h : int((len(avg_topo) - 1) / 2) + 5 + h]
                #   check the number of negative values in the lists
                neg_left = [i for i in left_vals if i < sea_level]
                neg_right = [i for i in right_vals if i < sea_level]
                #   append to the lists so we know how many negative values are there with increasing distance from the coast
                left_lst.append(len(neg_left))
                right_lst.append(len(neg_right))
            #   now loop through the overall lists and find the first value 5 - indicating that all the 5 values are below sea level
            left_coast_val  = next((x[0] for x in enumerate(left_lst) if x[1] == 5), None)
            right_coast_val  = next((x[0] for x in enumerate(right_lst) if x[1] == 5), None)
            #   now use these indexes (if they exist) to find the coastline offset value
            if left_coast_val is not None and right_coast_val is None:
                cst_offset_val  = final_topo[int((len(avg_topo) - 1) / 2) - left_coast_val]
            elif left_coast_val is None and right_coast_val is not None:
                cst_offset_val  = final_topo[int((len(avg_topo) - 1) / 2) + right_coast_val]
            #   if neither is None
            elif left_coast_val is not None and right_coast_val is not None:
                #   first check the position of first 0 on the left side (presumably the inland)
                left_coast_val  = next((x[0] for x in enumerate(left_lst) if x[1] == 0), None)
                cst_offset_val  = final_topo[int((len(avg_topo) - 1) / 2) - left_coast_val]
                
            #   if they dont exist then expand the search further from the coastline
            else:
                idx_range = 1000
                left_lst, right_lst = [], []
                #   loop through the list that is cut out around the assumed coastline - look in the vicinity for the submerged value
                for h in range(0, idx_range):
                    #   get first five values in each direction
                    left_vals = final_topo[int((len(avg_topo) - 1) / 2) - h - 5 : int((len(avg_topo) - 1) / 2) - h]
                    right_vals = final_topo[int((len(avg_topo) - 1) / 2) + h : int((len(avg_topo) - 1) / 2) + 5 + h]
                    #   check the number of negative values in the lists
                    neg_left = [i for i in left_vals if i < sea_level]
                    neg_right = [i for i in right_vals if i < sea_level]
                    #   append to the lists so we know how many negative values are there with increasing distance from the coast
                    left_lst.append(len(neg_left))
                    right_lst.append(len(neg_right))
                #   now loop through the overall lists and find the first value 5 - indicating that all the 5 values are below sea level
                left_coast_val  = next((x[0] for x in enumerate(left_lst) if x[1] == 5), None)
                right_coast_val  = next((x[0] for x in enumerate(right_lst) if x[1] == 5), None)
                #   now use these indexes (if they exist) to find the coastline offset value
                if left_coast_val is not None and right_coast_val is None:
                    cst_offset_val  = final_topo[int((len(avg_topo) - 1) / 2) - left_coast_val]
                elif left_coast_val is None and right_coast_val is not None:
                    cst_offset_val  = final_topo[int((len(avg_topo) - 1) / 2) + right_coast_val]
                #   if neither is None
                elif left_coast_val is not None and right_coast_val is not None:
                    #   first check the position of first 0 on the left side (presumably the inland)
                    left_coast_val  = next((x[0] for x in enumerate(left_lst) if x[1] == 0), None)
                    cst_offset_val  = final_topo[int((len(avg_topo) - 1) / 2) - left_coast_val]
                
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
        model_run = False

    else:
        print("Pica") 
        cst_off_idx = final_topo.index(cst_offset_val) - 1
        cst_offset = (cst_off_idx - int((len(avg_topo) - 1) * (500. / topo_res) / 10.)) * topo_res / 1000.
        model_run = True   
        
    #   decide where the coastline is based on the wtd dataset - where the first value is nan is where we put the coast
    wtd_cst = False
    """
    for g in range(cst_off_idx, len(avg_wtd_topo_refined)):
        if math.isnan(avg_wtd_topo_refined[g]):
            print(g)
            #cst_off_idx = g
            wtd_cst = True
            break
    """
    
    if model_run == True:
        print("Pica 4")    
        #   the Anchor point is defined as the point where it is first time indicated that the regolith
        #   thickness is equal or more than 50m, therefore put the anchor point at topo_elev - 50m
        #anchor_topo_elev = final_topo[int(round((cst_off_idx - (200 - x_nasa) * (1000 / topo_res))))]
        anchor_topo_elev = final_topo[int(round((cst_off_idx - x_nasa * (1000 / topo_res))))]
        y_nasa = anchor_topo_elev - 50.0    
    
        try: 
            #   first try to find the shelf break and foot of continental slope points
            print("Pica 5")    
            offshore_dir = os.path.join(out_dir, '_topo_figs')
            if not os.path.exists(offshore_dir):
                os.makedirs(offshore_dir)
            fos = find_FOS(final_topo, offshore_dir, 'avg_avg', avg = False) 
            print(fos) 
            model = cs_model(1, final_topo, avg_nasa_refined, x_nasa, y_nasa, cst_thick_est_avg, end_cst_plain,\
                             avg_soil_thk_refined, avg_soil_type_refined, avg_wtd_topo_refined, avg_offshore_refined,\
                             final_topo, final_topo, final_topo, avg_k_soil_refined,\
                             avg_drn_rate_refined, avg_glhymps_top_lay_refined, avg_glhymps_bot_lay_refined, topo_res, False)
            model.get_top_bot_lst_find_cst(cst_off_idx, fos[1], 0.0, wtd_cst)
            print('wtd_cst ' + str(wtd_cst), str(model.cst_idx))
            #model.get_top_bot_col_lst(del_col, cs_points_dist, model.cst_idx, smooth = True)
            #model.get_top_bot_col_lst_top_sys_geo(del_col, cs_points_dist, model.cst_idx, smooth = True)
            model.get_top_bot_col_lst(del_col, topo_res, model.cst_idx, smooth = True)
            model.get_top_bot_col_lst_top_sys_geo(del_col, topo_res, model.cst_idx, smooth = True)                    
            model.modflow_bot_elev_list(del_lay, True)
            print('CST_IDX          ', str(model.cst_idx))

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
                model.find_topo_divide_SRM_v2(10, -2000., del_col, del_lay, fos[1], 100., 20., 50000., 100000, 10000, 200000)
                #model.find_topo_divide(10., -10000., del_col, model.top_elev[model.idx_start], fos[1], 2, 200000.)         #   indicate the value of topographic difference to find the divide 
            #   in case the island is small cut off the inland part with the start of the model.top_elev list
            except IndexError:
                model.find_topo_divide_SRM(10, -2000., del_col, del_lay, fos[1], 20., 50000., 200000., 10000, 200000)
                #model.find_topo_divide(10., -10000., del_col, model.top_elev[0], fos[1], 2, 200000.)                      
            
            #try:
            #    model.find_topo_divide(10., -10000., del_col, model.top_elev[model.idx_start], fos[1], 2, 200000.)         #   indicate the value of topographic difference to find the divide 
            #   in case the island is small cut off the inland part with the start of the model.top_elev list
            #except IndexError:
            #    model.find_topo_divide(10., -10000., del_col, model.top_elev[0], fos[1], 2, 200000.)  

            model.dis_input(nrow, delc, del_col, del_lay, nper, perlen, nstp, laycbd, max_depth, cs_points_dist)
            #model.top = model.lay_elev[model.lay_idx]            #   top elevation of the model domain - check if it should be like this!
            model.top = model.lay_elev[0]      
            #model.botm = model.lay_elev[model.lay_idx + 1:]   
            model.botm = model.lay_elev[1:]   
            #model.zbot = model.botm[-1] 
            zbot = math.floor((model.bot_elev[-1] / 100.0) * 100.0)
                         
            ibound_arr = model.ibound_arr * 1.0
            ibound_arr[np.abs(ibound_arr) == 0.] = np.nan
            
            str_0 = out_name + '.png'
            plot_dir = os.path.join(out_dir, str_0)           
            
            x_start_new = model.x_start# + (model.cst_idx * del_col) / 1000.
            x_shift = round(model.x_start + (model.cst_idx * del_col) / 1000., 1)
            x_end_new = model.x_end # + (model.cst_idx * del_col) / 1000.
                   
            #model_start_idx = final_topo.index(min(final_topo, key=lambda x:abs(x-model.top_elev[0])))   
            #model_end_idx = final_topo.index(min(final_topo, key=lambda x:abs(x-model.top_elev[-1])))  
                   
            model_end_idx = int(round(model.idx_end * del_col / 1000.)) * 2 + cst_off_idx + int(round(x_start_new * 2))       
            #model_start_idx = 400 + int(round(x_start_new * 2))                           
                   
            topo_lst_plot = final_topo[cst_off_idx + int(round(x_start_new * 2)): model_end_idx + 1]
            #topo_lst_plot = final_topo[cst_off_idx + int(round(x_start_new * 2)): model_end_idx + 1]

            #x_axis_plot2 = np.linspace(x_start_new, x_end_new, model.ibound_arr.shape[-1])# - model.col_idx)
            #x_axis_plot = np.linspace(x_start_new, x_end_new, len(topo_lst_plot))

            if cst_offset >= 0:
                x_axis_plot2 = np.linspace(x_start_new - x_shift, x_end_new + x_shift, model.ibound_arr.shape[-1])# - model.col_idx)
                x_axis_plot = np.linspace(x_start_new - x_shift, x_end_new + x_shift, len(topo_lst_plot))
            else:
                x_axis_plot2 = np.linspace(x_start_new - x_shift, x_end_new + x_shift, model.ibound_arr.shape[-1])# - model.col_idx)
                x_axis_plot = np.linspace(x_start_new - x_shift, x_end_new + x_shift, len(topo_lst_plot))

            top_elev_arr = np.array(model.top_elev)
            mask = np.isnan(top_elev_arr)   
            top_elev_arr[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), top_elev_arr[~mask])
            
            model.top_elev = top_elev_arr.tolist()
                   
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
            
            #ax.plot(x_axis_plot2, [x1 - x2 for (x1, x2) in zip(model.top_elev, model.wtd_elev)], linewidth = 1, color = "blue")             
        
            
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
            ax_all.set_title(out_name)
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
            print('CST_IDX            ', str(model.cst_idx))
            #   fill in the dictionary                    
            topo_dict[out_name] = {'ibound_arr': model.ibound_arr,\
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
                                  #'cst_idx': cst_off_idx,\
                                  'cst_idx': model.cst_idx,\
                                  'wtd_elev': model.wtd_elev,\
                                  'k_soil': model.k_soil,\
                                  'drn_rate': model.drn_dens,\
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
    plot_name_all = os.path.join(out_dir, '__' + out_name + '_all_model_IBOUNDs.png')
    plt.savefig(plot_name_all, dpi = 300, facecolor = 'w', edgecolor = 'w',
        orientation = 'portrait', papertype = None, format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1,
        frameon = None)
    plt.close()

    #   save the dictionary with the ibound arrays
    dict_save_dir = os.path.join(out_dir, out_name + '_IBOUND.npy')
    np.save(dict_save_dir, topo_dict)

    return topo_dict, last_ibound_col_n_layers, ls_elev, ls_elev_50, fos_pt_dist, fos




















