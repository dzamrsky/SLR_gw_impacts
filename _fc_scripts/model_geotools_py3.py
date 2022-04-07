# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 09:39:55 2018

@author: daniel
"""

import random
import os
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cmx
import matplotlib 
from matplotlib import rcParams
from scipy.signal import argrelextrema
import rdp
import itertools

#save_dir = r'g:\Water_Nexus\_A2\_figures\_model_1341_geology_scenarios'

#   define all variables
#rand_seed = 1                   #   randomization seed, for reproductive purposes
laytyp_val = 0                  #   MODFLOW layer type, 0 = confined
sand_pct = 70                   #   sand percentage of the sediment volume (Maria)
clay_pct = 30                   #   clay (mud) percentage of the sediment volume (Maria)
n_aqt_in = 2                    #   number of aquitard layers in the inland part of the domain
aqt_in_x1 = -1.0                #   
n_aqt_off = 3                   #   number of aquitard layers in the offshore part of the domain
aqt_off_x0 = 1.0
soil_thk = 10.0                 #   
aq_vals_smooth_hor = False      #
top_soil = False                #   presence of top soil layer in the upper part of the inland domain
top_offshore = True             #   presence of top aquitard clay (mud) layer on top of the offshore part of the domain
rand_aqt_inland = True          #   
rand_aqt_offshore = True        #
end_at_cst = False              #
fos_point = None            
            
#model_obj = model  
#inland_aqf_scenario = 1
#inland_aqf_scenario = 2
#inland_aqf_scenario = 3
#   the list below specifies the thickness of each of the different layers in the bottom aquifer (in %)
#inland_aqf_lrs = [[20, 'glhymps_1'], [40, 'glhymps_2'], [10, 'glhymps_1'], [30, 'glhymps_2']]
#sed_type = 'medium'#'large' #  'small'
#sed_num_cells = 25. 
#off_lay_num = lay_num     #   make this equal to the number of layers inland
#sed_flux = 'medium' #'low' # 'high' 'medium' 


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


"""
line1 = (line_A_lst[j], line_B_lst[j])
line2 = line_thk_end
"""

def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1]) #Typo was here

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

def plot_Kh_arr_GLHMYPS_rand(model_obj, figname, cont_shelf_edge, end_col_idxs, end_lay_idxs, top_soil = False):
    #   calculate the offset of the coastline for plotting - this happens because the first offset is calculated from the
    #   raw topo data (from database with 0.5km spacing). When creating the top_elev interpolation the coast can be shifted by a bit
    cst_offset_val  = next((x for x in model_obj.top_elev if x < 0.0), None)
    cst_offset_idx = model_obj.top_elev.index(cst_offset_val) - 1
    model_obj.cst_offset_plot = round(model_obj.x_start + cst_offset_idx / (1000. /model_obj.del_col), 2) #  (1000. / model_obj.del_col) - to get the distance in km
    #   create a plot if desired
    plot_hk_arr = model_obj.hk_arr
    plot_hk_arr[np.abs(plot_hk_arr) == 0.] = np.nan
    #   find all the unique values in the KH array, remove the ones that are marked as no-value
    unique = np.unique(model_obj.hk_arr)
    unique_nan = [value for value in unique if not math.isnan(value)]
    #   define the colormap
    cmap = cmx.viridis
    #cmaplist = [cmap(i) for i in range(cmap.N)]
    # define the bins and normalize
    
    #unique_n = [0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 2.5, 5.0] + list(np.arange(10.0, round(max(unique_nan) / 10) * 10, 10.))
    unique_n = [0.0, 0.001, 0.1, 1.0, 5.0] + list(np.arange(10.0, round(max(unique_nan) / 10) * 10, 10.))
    norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)            

    # Concentration and Flow Plot
    fig = plt.figure(figsize = (20, 20))
    ax1 = plt.subplot2grid((2, 1), (0, 0))  #   overall concentration profiles
    ax2 = plt.subplot2grid((2, 1), (1, 0))  #   color bar area 

    ax1.set_position([0.05, 0.15, 0.9, 0.8])
    ax2.set_position([0.3, 0.025, 0.4, 0.05])        
        
    im1 = ax1.imshow(plot_hk_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                   extent = (model_obj.x_start, model_obj.x_end, model_obj.zbot, math.ceil((model_obj.top / 100.0) * 100.0)), vmin = 0., vmax = np.amax(plot_hk_arr))

    x_lines_topelev = np.linspace(model_obj.x_start + (model_obj.del_col / (2 * 1000.)), model_obj.x_end - (model_obj.del_col / (2 * 1000.)), model_obj.ncol)
    y_lines = np.linspace(model_obj.top, model_obj.zbot, model_obj.ibound_arr.shape[0] + 1)     
    
    lineA, = ax1.plot(x_lines_topelev, model_obj.top_elev, c = 'black', linewidth = 3, label = 'GEBCO elevation')
    ax1.set_xlim([model_obj.x_start, model_obj.x_end])

    ax1.set_ylim([min(model_obj.bot_elev), max(model_obj.top_elev)])
    ax1.axhline(y = 0., linewidth = 1, color = 'k')
    ax1.axvline(x = 0., linewidth = 1, color = 'k')
    ax1.axvline(x = cont_shelf_edge[0], linewidth = 1, color = 'k')    
    
    for x_line in end_col_idxs:
        ax1.axvline(x = x_lines_topelev.tolist()[x_line], linestyle='--', linewidth = .5, color = 'k')            
    for y_line in end_lay_idxs:
        ax1.axhline(y = y_lines.tolist()[y_line], linestyle='--', linewidth = .5, color = 'k')  
        
    ax1.set_title('Horizontal hydraulic conductivity distribution in the model_obj domain', fontsize = 26)
    ax1.set_xlabel('distance from coast (km)', fontsize = 14)
    ax1.set_ylabel('elevation (m asl.)', fontsize = 14)

    #   plot the colorbar
    cbar = plt.colorbar(im1, cax = ax2, cmap = cmap, norm = norm, spacing = 'uniform', ticks = unique_n, boundaries = unique_n, orientation='horizontal')
    ax2.set_xticklabels([str(e) for e in unique_n])   
    ax2.set_title("Horizontal hydraulic conductivity (m/d)", fontsize = 14)#, y  = 1.025)
    ax2.tick_params(labelsize = 12)    

    rcParams['font.family'] = 'Garamond'
    rcParams['axes.facecolor'] = 'white'
    rcParams['savefig.facecolor'] = 'white'  

    plt.show()

    #plt.savefig(os.path.join(save_dir, figname))
    #plt.close()
    #del fig


#   Function that finds the FOS (foot of the contintental slope), based on a paper
"""
model_obj = model
topo_lst = model_obj.top_elev
"""
def find_shelf_break(model_obj, topo_lst, avg = True):
    
    #topo_lst = model.top_elev
    #   trim the topo values in case there is also sea/ocean at the left side of the land (land looks like an island)
    idx_neg = 0
    while topo_lst[idx_neg] < 0.:
        idx_neg += 1

    topo_lst = topo_lst[idx_neg:]

    cst_val  = next((x for x in topo_lst if x < 0.0), None)
    cst_idx = topo_lst.index(cst_val) - 1

    #cst_val  = min(topo_lst, key=lambda x:abs(x - 0.))  
    #cst_idx = topo_lst.index(cst_val)        

    x_axis = np.linspace(0.0, 200.0, 401)
    x_axis = np.arange(0.0, 200.0, 0.1)
    #x_axis = np.arange(0.0, model.x_end, 0.1)

    off_topo_lst = topo_lst[cst_idx:]

    try:
        #   calculate the gradient and second derivative
        grad = np.gradient(off_topo_lst) 
        grad2 = np.gradient(grad) 
        
        loc_max = argrelextrema(grad2, np.greater)    
        loc_min = argrelextrema(grad2, np.less)    
        
        loc_minmax = np.concatenate((loc_max[0], loc_min[0])) # we still need to do this unfortunatly.
        loc_minmax.sort()
        
        loc_max_plt, loc_min_plt, loc_minmax_plt1, loc_minmax_plt2 = [], [], [], []
    
        for a in range(loc_max[0].shape[0]):
            loc_max_plt.append([x_axis[loc_max[0][a]], off_topo_lst[loc_max[0][a]]])
        
        for b in range(loc_min[0].shape[0]):
            loc_min_plt.append([x_axis[loc_min[0][b]], off_topo_lst[loc_min[0][b]]])
            
        for c in range(loc_minmax.shape[0]):
            #loc_minmax_plt1.append([avg_topo[loc_minmax[c]], avg_topo[loc_minmax[c]]])
            loc_minmax_plt1.append(off_topo_lst[loc_minmax[c]])
            loc_minmax_plt2.append(off_topo_lst[loc_minmax[c]])
        
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
                slope_right = (rdp_points[j][1] - off_topo_lst[-1]) / (model_obj.x_end - rdp_points[j][0])            
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
            #   first find the sb_point, it should be above -200m bsl. and the slope on the right from the point
            #   should be higher than the slope on the left from the point
            if pt_depth > -200.0 and abs(pt_slope_left) < abs(pt_slope_right):
                sb_point = [round(pt_dist, 1), pt_depth]
            #   the fos_point should be deeper than -200m bsl, the absolute slope gradient on the left should be higher
            #   than on the right, and the slope on the right is either positive or less than 25% than the absolute value
            #   of the slope on the left
            elif pt_depth < -1000.0 and abs(pt_slope_left) > abs(pt_slope_right):
                if pt_slope_right > 0 or abs(pt_slope_right) < (50 * abs(pt_slope_left)) / 100.:
                    #   check if it is just one peak (concave hull) in the profile, do this by checking the slope and 
                    #   deptg(elev) of the following part of the topographical profile
                    if rdp_points[g + 1][1] > pt_depth and abs(50 * rdp_points[g + 1][3]) / 100. > abs(pt_slope_left):
                        continue
                    else:
                        fos_point = [round(pt_dist, 1), pt_depth]
                        break
                else:
                    pass
                
    except (ValueError, ZeroDivisionError, IndexError):
        print('Not possible to find shelf break based on slope analysis - the shelf decline is too steep probably and matches the slope of the continental slope')
        #   instead, return the value where the bathymetry reaches the 125.m depth 
        shelf_val  = min(topo_lst, key=lambda x:abs(x-(-125.)))  
        
        #shelf_idx = topo_lst.index(shelf_val)        
        shelf_idx = len(topo_lst) - topo_lst[::-1].index(shelf_val)  
        sb_point = [(shelf_idx - cst_idx) / 10., shelf_val]
        
    if sb_point is None:
        print('Not possible to find shelf break based on slope analysis - the shelf decline is too steep probably and matches the slope of the continental slope')
        #   instead, return the value where the bathymetry reaches the 125.m depth 
        shelf_val  = min(topo_lst, key=lambda x:abs(x-(-125.)))     
        shelf_idx = topo_lst.index(shelf_val)        
        sb_point = [(shelf_idx - cst_idx) / 10., shelf_val]
        
    try:
        return sb_point
    
    except UnboundLocalError:
        print ('Couldnt find shelf break..')
        return [None]


"""
model_obj = model
inland_aqf_lrs = inland_aqf_lrs
sed_type = 'small'
#sed_flux = 'medium'
save_dir = topo_sc_dir
summary_save_dir = geo_summary_csv_dir

glh_1_mu, glh_1_std = glhymps_bot_param_ln[0][0], glhymps_bot_param_ln[1][0]
glh_2_mu, glh_2_std = glhymps_top_param_ln[0][0], glhymps_top_param_ln[1][0]
rand_seed_in = rand_seed

const_geo_hk_vals = False
glh_1_val = 0.1
glh_2_val = 10.
clay_val = 0.0001
rand_seed_in = 1

p_fact = 0.5
off_lay_thk_ratio = mud_pct
lay_pres_y1 = y_val           
clay_cap_shelf = True       #   insert the clay capping layer on top of the continental shelf/slope, the same reworking
clay_cap_slope = True       #   parameter lay_pres_y1 will apply to these layers as well
clay_cap_thk = 20.    


model, rand_seed, inland_aqf_lrs, p_fact, off_lay_thk_ratio, sed_type, lay_pres_y1, clay_cap_thk,\
                                                           off_lay_start, 'low', topo_sc_dir, const_geo_hk_vals, clay_cap_shelf, clay_cap_slope

mud_shlf_pct, mud_slp_pct
"""
 
#def create_geology_profile(model_obj, rand_seed_in, inland_aqf_lrs, p_fact, off_lay_thk_ratio, sed_type, lay_pres_y1, clay_cap_shelf_thk, clay_cap_slope_thk,\
#                           off_lay_start, sed_flux, save_dir, summary_save_dir, const_geo_hk_vals, clay_cap_shelf, clay_cap_slope, figname):           

"""
model_obj = model
rand_seed_in = rand_seed
inland_aqf_lrs
p_fact
off_lay_thk_ratio = mud_pct
sed_type
lay_pres_y1
clay_cap_shelf_thk
clay_cap_slope_thk
off_lay_start
sed_flux
save_dir = sc_dir
summary_save_dir = geo_summary_csv_dir
const_geo_hk_vals
clay_cap_shelf
clay_cap_slope
figname

"""


def create_geology_profile(model_obj, rand_seed_in, inland_aqf_lrs, p_fact, off_lay_thk_ratio, sed_type, lay_pres_y1, clay_cap_shelf_thk, clay_cap_slope_thk,\
                           off_lay_start, sed_flux, save_dir, summary_save_dir, const_geo_hk_vals, clay_cap_shelf, clay_cap_slope, figname):  
           
    random.seed(rand_seed_in)             
           
    """                    First prepare the datasets and lists                     """           
    cont_shelf_edge = find_shelf_break(model_obj, model_obj.top_elev, True) #   try to find the shelf break              
              
    #   find the coastal index - that will determine the limit between the inland and offshore domain
    #cst_offset_idx_clip = (model_obj.cst_idx - model_obj.idx_start) * 5 + 1  
     
    cst_offset_idx_clip =  int(abs(round(model_obj.x_start, 1) * 10) - 1)
    #cst_offset_val = model_obj.top_elev[cst_offset_idx_clip]
    #cst_offset_val  = next((x for x in model_obj.top_elev if (x < -10.0 and x < int(abs(model_obj.x_start) * 10)) or (x < 0.0 and x => int(abs(model_obj.x_start) * 10))), None)

    cst_offset_val = next((x for x in model_obj.top_elev[:int(abs(model_obj.x_start) * 10)] if x < -10.0), None)
    if not cst_offset_val:
        cst_offset_val = next((x for x in model_obj.top_elev[int(abs(model_obj.x_start) * 10):] if x < 0.0), None)
        cst_offset_idx = int(abs(model_obj.x_start) * 10) + model_obj.top_elev[int(abs(model_obj.x_start) * 10):].index(cst_offset_val) - 1
    else:
        cst_offset_idx = model_obj.top_elev.index(cst_offset_val) - 1
        
    model_obj.cst_offset_plot = round(model_obj.x_start + cst_offset_idx / (1000. /model_obj.del_col), 2) #  (1000. / model_obj.del_col) - to get the distance in km
          
    #cst_offset_idx = int(model_obj.x_start * -10.)  
       
    #   trim the soilgrids and replace all potential NaN values with the mean value of all the other thicknesses
    model_obj.soilgrids_thk = model_obj.soil_thk[: cst_offset_idx_clip]  
    soilgrids_thk_lst = [round(p / 100., 0) for p in model_obj.soilgrids_thk if p > -1]
    soilgrids_thk_mean = round(np.nanmean(soilgrids_thk_lst, axis=0), 0)
    soilgrids_thk_lst = [soilgrids_thk_mean if math.isnan(p) else p for p in soilgrids_thk_lst]
    
    #   get the mean and stdev values for glhymps 1 and 2
    #rand_val_glh_1 = np.exp(np.random.normal(glh_1_mu, glh_1_std, size = 1))[0]      
    #rand_val_glh_2 = np.exp(np.random.normal(glh_2_mu, glh_2_std, size = 1))[0] * 3600 * 24  
    
    #   trim and post-process the GLHYMPS and GLHYMPS 2.0 input lists
    model_obj.hk_vals_bot = model_obj.glhymps_bot_lay[: cst_offset_idx_clip]
    model_obj.hk_vals_top = model_obj.glhymps_top_lay[: cst_offset_idx_clip]
    
    #   check that the glhymps 2 layer is composed of non-nan values
    glh_2_check = [i for i in model_obj.hk_vals_top if not math.isnan(i)]
    
    if glh_2_check == []:
        model_obj.hk_vals_top = model_obj.k_soil[: cst_offset_idx_clip]
        model_obj.hk_vals_top = [round(i, 4) for i in model_obj.k_soil]
    else:
        model_obj.hk_vals_top = [round(i * 3600 *24, 4) for i in model_obj.glhymps_top_lay]
            
    #   recalculate the values based on the GLHYMPS 2 formula
    model_obj.hk_vals_bot = [round(i, 4) for i in model_obj.glhymps_bot_lay]

    #   replace the negative values - sometimes they are there for some reason, also the totaly unrealistic huge values, limit to 100. m/d
    for x in range(len(model_obj.hk_vals_top)):
        if model_obj.hk_vals_top[x] < 0. or model_obj.hk_vals_top[x] > 100.:
           model_obj.hk_vals_top[x] = np.mean([i for i in model_obj.hk_vals_top if i > 0 and i < 100.])
        try:
            if model_obj.hk_vals_bot[x] < 0. or model_obj.hk_vals_bot[x] > 100.:
               model_obj.hk_vals_bot[x] = np.mean([i for i in model_obj.hk_vals_bot if i > 0 and i < 100.])
        #   happens when the hk_vals_bot list is shorter than the top one, not sure why it happens sometimes
        #   if it happens then just skip the index in the loop, go to the next one 
        except IndexError:
            pass
    
    #   calculate the statistics for both top and bottom aquifer layers
    hk_top_mean = np.nanmean(model_obj.hk_vals_top, axis=0)
    hk_top_std = np.nanstd(model_obj.hk_vals_top, axis=0)
    hk_bot_mean = np.nanmean(model_obj.hk_vals_bot, axis=0)
    hk_bot_std = np.nanstd(model_obj.hk_vals_bot, axis=0)    

    #   limit the stdev value to 50% of the mean?
    if hk_bot_std > abs(hk_bot_mean):
        hk_bot_std = abs(hk_bot_mean)
    if hk_top_std > abs(hk_top_mean):
        hk_top_std = abs(hk_top_mean)

    glh_1_val = round(hk_bot_mean, 4)
    glh_2_val = round(hk_top_mean, 4)
    clay_val = 0.0001
    
    
    #   create the HK_ARR array that will be filled in the next steps
    #model_obj.hk_arr = np.zeros([model_obj.nlay, 1, model_obj.ncol])
    model_obj.hk_arr = np.zeros([model_obj.ibound_arr.shape[0], 1, model_obj.ibound_arr.shape[-1]])
    
    """                Next, create the inland part of the domain                   """       
    #   There are two different scenarios on how to create the stratigraphy in the inland part of the domain
    #       1) Based on the SOILGRIDS thickness (soilgrids_thk_lst) fill the upper part of the model_obj domnain
    #          with the GLHYMPS 2.0 values and the whole rest of the model_obj domain located below thet upper 
    #          part gets the HK values based on GLHMYPS 1.0 dataset. 
    #       2) This scenario is based on the assumption that more permeable layers are interlayed with less
    #          permeable ones - as mentiond in GLHYMPS 2.0 article, its values are usually 10 times higher than
    #          GLHYMPS 1.0 - that is why the SOILGRIDS thickness is repeated through the model_obj domain and assigned
    #          the GLHYMPS 1.0 or 2.0 values creating a sort of a zebra pattern.
    #       3) The last scenario takes into account the SOILGRIDS thickness and fills the upper aquifer part 
    #          with GLHYMPS 2.0 values. However, the rest of the aquifer is split into layers with either 
    #          GLHYMPS 1.0 or 2.0 based on an input list that defines the thickness and GLHYMPS type of each layer

    #   create a list where the glhymps 1 and glhymps 2 layers and respective layer indexes for each column will be stored
    #   only for inland part of the model domain.
    inland_lrs_lst = []

    #   loop through each model_obj domain column located in the inland part 
    #for r in range(len(soilgrids_thk_lst)):
    #   calculate the index of the coastline based on the x start and the offset of the topography
    cst_idx_geo = int(abs(model_obj.x_start - model_obj.cst_offset_plot) * 10)
    
    for r in range(cst_idx_geo):
        if len(soilgrids_thk_lst) > 0:
            try:
                #   calculate the number of layers that are part of the upper aquifer
                col_top_lay_n = int(round(soilgrids_thk_lst[r] / model_obj.del_lay))
            #   might be that the soilgrids list is shorter than the coastline, in that case take the last value
            except IndexError:
                col_top_lay_n = int(round(soilgrids_thk_lst[-1] / model_obj.del_lay))
        else:
            col_top_lay_n = 0
        #   fill the upper part of the aquifer in any case (same for both scenarios), the values are randomized
        #   based on the statistics of the GLHYMPS 2.0 list calculated in the previous step
        ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, r].tolist()) if x == 1]
        lrs_lst = []
        
        #lrs_lst.append(['glhymps_2', ibound_act_lay_idxs[:col_top_lay_n]])
        
        if abs(glh_2_val) > abs(glh_1_val):
            lrs_lst.append(['glhymps_2', ibound_act_lay_idxs[:col_top_lay_n]])
        else:            
            lrs_lst.append(['glhymps_1', ibound_act_lay_idxs[:col_top_lay_n]]) 
            
        for s in range(col_top_lay_n):
            #   try except - catch the IndexError in case the SOILGRIDS thickness is larger than the active thickness
            #   of the IBOUND array in the current column
            try:
                if not const_geo_hk_vals:
                    #model_obj.hk_arr[ibound_act_lay_idxs[s], 0, r] = round(abs(np.exp(np.random.normal(glh_2_mu, glh_2_std, size = 1))[0] * 3600 * 24), 4)
                    #model_obj.hk_arr[ibound_act_lay_idxs[s], 0, r] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                    hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                    if hk_cell_val > hk_top_mean / 10.:
                        model_obj.hk_arr[ibound_act_lay_idxs[s], 0, r] = hk_cell_val
                    else:
                        model_obj.hk_arr[ibound_act_lay_idxs[s], 0, r] = hk_top_mean / 10.
                else:
                    model_obj.hk_arr[ibound_act_lay_idxs[s], 0, r] = glh_2_val
            except IndexError:
                pass
    
        #   check if there are any ibound_act_lay_idxs indexes below the upper aquifer part
        if len(ibound_act_lay_idxs) > col_top_lay_n:
            #   split the active ibound cells in the column list into groups (geological layers) based on the % from the input list
            tot_cells = len(ibound_act_lay_idxs) - col_top_lay_n
            cell_cnt_lst = []
            geo_lay_thk_start = col_top_lay_n
                
            #   get the pct of thickness per model layer
            lay_thk_pct = round(1 / float(tot_cells) * 100., 0)
            
            #   decide the layer types in the rest of the model layers based on the majority of glhymps type in that layer
            lay_glh_lst = []

            # in very deep systems the pct can be equal to 0, in that case round to one more decimal
            if lay_thk_pct != 0.0:
                lay_thk_pct = math.ceil(1 / float(tot_cells) * 100.)
  
                #   create a list of 100 values (1 or 2 depending on the glhymps type)
                for i in range(len(inland_aqf_lrs)):
                    glh_pct = inland_aqf_lrs[i][0]
                    glh_type = inland_aqf_lrs[i][1]   

                    #   adapt the number of cells to be assigned
                    new_cell_n = int(round((tot_cells / 100) * glh_pct, 0))
                    
                    #print(glh_type, glh_pct, new_cell_n)
                    
                    for j in range(new_cell_n):
                        if glh_type == 'glhymps_1':
                            lay_glh_lst.append(1)
                        else:
                            lay_glh_lst.append(2)                        
                    
                    """                    
                    for j in range(glh_pct):
                        if glh_type == 'glhymps_1':
                            lay_glh_lst.append(1)
                        else:
                            lay_glh_lst.append(2)                    
            
                    #   split this list into even chunks based on the lay_thk_pct value
                    lay_glhymps_lst = [lay_glh_lst[i:i + int(lay_thk_pct)] for i in range(0, len(lay_glh_lst), int(lay_thk_pct))]         
                    """
                    
                lay_glhymps_lst = [lay_glh_lst[i:i + 1] for i in range(0, tot_cells, 1)]
                    
            else:
                
                #   create a list of 100 values (1 or 2 depending on the glhymps type)
                for i in range(len(inland_aqf_lrs)):
                    glh_pct = inland_aqf_lrs[i][0]
                    glh_type = inland_aqf_lrs[i][1]  
                    
                    #   adapt the number of cells to be assigned
                    new_cell_n = int(round((tot_cells / 100) * glh_pct, 0))
                    
                    #print(glh_type, glh_pct, new_cell_n)
                    
                    for j in range(new_cell_n):
                        if glh_type == 'glhymps_1':
                            lay_glh_lst.append(1)
                        else:
                            lay_glh_lst.append(2)                        
                
                lay_glhymps_lst = [lay_glh_lst[i:i + 1] for i in range(0, tot_cells, 1)]
                
            if lay_glhymps_lst != [[]] and [i for i in lay_glhymps_lst if i != []] != []:
                last_non_empty = [z[0] for z in lay_glhymps_lst if z != []][-1]
                
            else:
                glh_type_tofill = max([sublist for sublist in inland_aqf_lrs])[1]
                if glh_type_tofill == 'glhymps_2': 
                    last_non_empty = 2
                elif glh_type_tofill == 'glhymps_1': 
                    last_non_empty = 1            
            
            for x in range(len(lay_glhymps_lst)):
                if lay_glhymps_lst[x] == []:
                    lay_glhymps_lst[x] = [last_non_empty]
            
            #   make sure the list above has the right amount of elements (in case the % is 33.0 for example it appends one extra value to reach length of 100)
            if len(lay_glhymps_lst) > tot_cells:
                lay_glhymps_lst = lay_glhymps_lst[:-1]

            #   go through the list and depending on the majority of values in each sublist decide which GLHYMPS value will be filled in
            for l in range(len(lay_glhymps_lst)):
                lay = lay_glhymps_lst[l]
                glh_typ = max(set(lay), key = lay.count)
                #print glh_typ
                
                #   due to rounding of the  tot_cells above sometimes the ibound_act_lay is one index shorter
                try:
                    if not const_geo_hk_vals: 
                        if glh_typ == 1:
                            model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, r] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                        elif glh_typ == 2:
                            #model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, r] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                            hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                            if hk_cell_val > hk_top_mean / 10.:
                                model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, r] = hk_cell_val
                            else:
                                model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, r] = hk_top_mean / 10.
                    else:
                        if glh_typ == 1:
                            model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, r] = glh_1_val
                        elif glh_typ == 2:
                            model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, r] = glh_2_val         
                
                except IndexError:
                    #print('index error in column : ' + str(r) + 'len(ibound_act_lay_idxs) : ' + str(len(ibound_act_lay_idxs)) + ',  col_top_lay_n + l : ' + str(col_top_lay_n + l))
                    try:
                        if not const_geo_hk_vals: 
                            if glh_typ == 1:
                                model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l - 1], 0, r] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                            elif glh_typ == 2:
                                #model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l - 1], 0, r] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                if hk_cell_val > hk_top_mean / 10.:
                                    model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l - 1], 0, r] = hk_cell_val
                                else:
                                    model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l - 1], 0, r] = hk_top_mean / 10.
                        else:
                            if glh_typ == 1:
                                model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l - 1], 0, r] = glh_1_val
                            elif glh_typ == 2:
                                model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l - 1], 0, r] = glh_2_val                
                    except IndexError:
                        if not const_geo_hk_vals: 
                            if glh_typ == 1:
                                model_obj.hk_arr[ibound_act_lay_idxs[-1], 0, r] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                            elif glh_typ == 2:
                                #model_obj.hk_arr[ibound_act_lay_idxs[-1], 0, r] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                if hk_cell_val > hk_top_mean / 10.:
                                    model_obj.hk_arr[ibound_act_lay_idxs[-1], 0, r] = hk_cell_val
                                else:
                                    model_obj.hk_arr[ibound_act_lay_idxs[-1], 0, r] = hk_top_mean / 10.
                        else:
                            if glh_typ == 1:
                                model_obj.hk_arr[ibound_act_lay_idxs[-1], 0, r] = glh_1_val
                            elif glh_typ == 2:
                                model_obj.hk_arr[ibound_act_lay_idxs[-1], 0, r] = glh_2_val    

            glh_nr = 1
            glh_tx = 'glhymps_1'
            to_lst = []

            #   go through the list and depending on the majority of values in each sublist decide which GLHYMPS value will be filled in
            for l in range(tot_cells):
                if glh_typ == glh_nr:
                    to_lst.append(ibound_act_lay_idxs[col_top_lay_n + l])
                else:
                    lrs_lst.append([glh_tx, to_lst])
                    if glh_nr == 1:
                        glh_nr = 2
                        glh_tx = 'glhymps_2'
                    else:
                        glh_nr = 1
                        glh_tx = 'glhymps_1'                    
                    to_lst = []

        inland_lrs_lst.append([r, lrs_lst])
            
    #   define input into the next step for creating the offshore part of the model_obj domain
    tot_cells_cst = len([i for i, x in enumerate(model_obj.ibound_arr[:, 0, r].tolist()) if x == 1])
    
    #lay_num = len(inland_aqf_lrs) + 1
    lay_thk = []
    #   lay_thk_1 is going to be repeated for all the layers because they all have same thickness
    if len(soilgrids_thk_lst) > 0:
        try:
            lay_thk_1 = round((int(soilgrids_thk_lst[r] / model_obj.del_lay) / float(tot_cells_cst)) * 100., 2) 
        except IndexError:
            lay_thk_1 = round((int(soilgrids_thk_lst[-1] / model_obj.del_lay) / float(tot_cells_cst)) * 100., 2) 
    else:
        lay_thk_1 = 1
    lay_thk.append([lay_thk_1, 'glhymps_2'])
    for w in range(len(inland_aqf_lrs)):
        try:
            geo_lay_thk = (inland_aqf_lrs[w][0] * tot_cells) / 100
            lay_thk.append([round((geo_lay_thk / float(tot_cells_cst)) * 100., 2), inland_aqf_lrs[w][1]])
        # NameError: name 'tot_cells' is not defined
        except NameError: 
            lay_thk.append([0.0, inland_aqf_lrs[w][1]])
    lay_thk[-1][0] += round(100. - sum(i[0] for i in lay_thk), 2)

    #matplotlib.use('TkAgg')
    #show_arr = model_obj.hk_arr
    #lt.imshow(show_arr[:, 0, :])
    #plt.show()

    """       Once the inland part of the domain is filled, create the offshore part         """
    #   The offshore part of the model_obj domain is conceptualized as the unconsolidated sediment part of the contintental
    #   shelf. Since there are quite a lot of differences around the coastline and it is impossible to model_obj each individual
    #   case. Therefore the classification by Maria is used to classify the global continental shelf types based on the 
    #   sediment grain size deposited and the shape of the chronologically and stratigraphically different layers. There can
    #   be a present/absent clay aquitard layer in between these individual layers. This is all specified by various classes 
    #   that are specified below.
    #
    #       1) The first class to be specified is the sediment type based on Marias research. This will have effect on the
    #          upper composition of each layer. The following types are specified for variable sed_type 
    #           
    #          a) 'large' sediment size means that the top of each layer will have a coarder (*10) HK values than supposed
    #          b) 'medium' sediment size means that there is no change in permeabilities
    #          c) 'small' sediment size means that there is a clay (mud) layer on top of each layer
    #
    #          The other vairable specified in this class is the sed_num_cells which determines the number of cells (in %) for 
    #          each layer that gets the sediment type property specified above. 
    #          
    #       2) Number of different sediment layers off_lay_num (specified as GLHYMPS 1.0 or 2.0) can:
    #
    #          a) be specified based on the inland_aqf_scenario is use_inland_aqf_scenario - True/False
    #               - if True then just follow the num_layers from the inland part of the domain
    #               - if False then specify the off_lay_num manually 
    #          b) another variable here is the thickness of each layer (in % of the verical domain space) off_aqf_lrs which can
    #             either be:
    #               - equal to inland_aqf_lrs
    #               - specified manuall (but then the connection at the coastline will look a bit strange..)
    #
    #       3) The last class of variables to be specified concerns the sediment flux type, which can be
    #
    #          a) 'high' sediment flux means that the different layers are more extended in the horizontal direction, meaning that
    #             the axis_angle is larger 
    #          b) 'medium' sediment flux translates into an average axis_angle that is lower than in case of the previous type
    #          c) 'low' sediment flux means that the layers are shorter and shorter in the x direction the younger they are. No
    #             angle is necessary here as the layers do not extend and 'slide down' around the continental slope
    #   
    #          axis_angle is the angle between the coastline (verical x = 0) and the line that passes through the end of the 
    #          continental shelf point and the ATE point at the x = 0 line. The intersection of this line and each of the 
    #          horizontal layers lines defines the inflection point of each sediment layer.
        
    #   get the shelf index in and also the last layer of the upper aquifer part in that column
    shelf_edge_idx = model_obj.top_elev.index(cont_shelf_edge[1]) 
    shelf_edge_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edge_idx].tolist()) if x == 1]      

    #   First create different layers based on the layer number and the sediment influx, because that determines
    #   how the layers are positioned and thick
    
    #   fill the upper part of the aquifer in any case (same for both scenarios), the values are randomized
    #   based on the statistics of the GLHYMPS 2.0 list calculated in the previous step
    #ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, a].tolist()) if x == 1]
    #for s in xrange(col_top_lay_n):
        #   try except - catch the IndexError in case the SOILGRIDS thickness is larger than the active thickness
        #   of the IBOUND array in the current column
    #    try:
    #        model_obj.hk_arr[ibound_act_lay_idxs[s], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
    #    except IndexError:
    #        pass
    
    #   create a list where the glhymps 1 and glhymps 2 layers and respective layer indexes for each column will be stored
    #   only for offshore part of the model domain, the inland part will or wont have its own clay layers formed later
    offshore_lrs_lst = []
    
    #   first create profiles if the sed influx is low low sediment influx means that the layers are going
    #   to be split based on the ratio specified above, and will spread till the continental slope
    if sed_flux == 'low': 
        #   loop through each model_obj domain column located in the offshore
        #for a in range(cst_idx_geo, len(model_obj.top_elev)):
        for a in range(cst_idx_geo, model_obj.ibound_arr.shape[-1]):
            #   calculate the number of layers that are part of the upper aquifer, use the mean SOILGRIDS thickness value
            #col_top_lay_n = int(soilgrids_thk_mean / model_obj.del_lay) 
            if len(soilgrids_thk_lst) > 0:
                try:
                    col_top_lay_n = int(soilgrids_thk_lst[r] / model_obj.del_lay)
                except IndexError:
                    col_top_lay_n = int(soilgrids_thk_lst[-1] / model_obj.del_lay)             
            else:
                col_top_lay_n = 0
                
            ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, a].tolist()) if x == 1]

            #   create a list for to store the layer indexes for each of the sediment layers
            lrs_lst = []               
            
            #   first check if there is a edge of the continental shelf point found
            #       if yes, then create a line between the point and the continental slope point
            if cont_shelf_edge is not None:
                
                #   calculate the amount of layers at the continental shelf column
                tot_layers = len(inland_aqf_lrs)
                #   divide the area between the shelf edge and model_obj bottom into the tot_layers identical amount of cells

                #   find the limits for each of the geological layers 
                idx_1 = shelf_edge_ibound_act_lay_idxs[col_top_lay_n]                
                end_lay_idxs =  [idx_1]
                act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[idx_1, 0, :].tolist()) if x == 1]
                end_col_idxs = [act_cell_lst[-1]]                
                for b in range(tot_layers - 1):                    
                    lay_idx = idx_1 + int(round((inland_aqf_lrs[b][0] * (model_obj.ibound_arr.shape[0] - shelf_edge_ibound_act_lay_idxs[col_top_lay_n])) / 100, 0))              
                    end_lay_idxs.append(lay_idx)
                    act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[lay_idx - 1, 0, :].tolist()) if x == 1]
                    if len(act_cell_lst) > 0:
                        end_col_idxs.append(act_cell_lst[-1])
                        idx_1 = lay_idx
                #   append the last active layer and column as the end of the last layer
                end_lay_idxs.append(model_obj.ibound_arr.shape[0])
                end_col_idxs.append(int(model_obj.ncol) - 1)
 
                #   check if there are any ibound_act_lay_idxs indexes below the upper aquifer part
                #if len(ibound_act_lay_idxs) > col_top_lay_n:
                #   fill in the part till the shelf edge with the constant thickness, after that point fill in straight line
                if a < shelf_edge_idx:                    
                    #   the upper part of the offshore aquifer domain
                    ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, a].tolist()) if x == 1]
                    
                    """     add if glhymps_top is true or false     """                    
                    
                    lrs_lst.append(['glhymps_2', ibound_act_lay_idxs[:col_top_lay_n]]) 
                    for s in range(col_top_lay_n):
                        #   try except - catch the IndexError in case the SOILGRIDS thickness is larger than the active thickness
                        #   of the IBOUND array in the current column
                        try:
                            if not const_geo_hk_vals: 
                                #model_obj.hk_arr[ibound_act_lay_idxs[s], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                if hk_cell_val > hk_top_mean / 10.:
                                    model_obj.hk_arr[ibound_act_lay_idxs[s], 0, a] = hk_cell_val
                                else:
                                    model_obj.hk_arr[ibound_act_lay_idxs[s], 0, a] = hk_top_mean / 10.
                            else:
                                model_obj.hk_arr[ibound_act_lay_idxs[s], 0, a] = glh_2_val
                        except IndexError:
                            pass    
                        
                    #   split the active ibound cells in the column list into groups (geological layers) based on the % from the input list
                    tot_cells = len(ibound_act_lay_idxs) - col_top_lay_n
                    #cell_cnt_lst = []
                    #geo_lay_thk_start = col_top_lay_n
        
                    #   get the pct of thickness per model layer
                    lay_thk_pct = round(1 / float(tot_cells) * 100., 0)
                    
                    #   decide the layer types in the rest of the model layers based on the majority of glhymps type in that layer
                    lay_glh_lst = []


                    # in very deep systems the pct can be equal to 0, in that case round to one more decimal
                    if lay_thk_pct != 0.0:
                        lay_thk_pct = math.ceil(1 / float(tot_cells) * 100.)
          
                        #   create a list of 100 values (1 or 2 depending on the glhymps type)
                        for i in range(len(inland_aqf_lrs)):
                            glh_pct = inland_aqf_lrs[i][0]
                            glh_type = inland_aqf_lrs[i][1]   
        
                            #   adapt the number of cells to be assigned
                            new_cell_n = int(round((tot_cells / 100) * glh_pct, 0))
                            
                            #print(glh_type, glh_pct, new_cell_n)
                            
                            for j in range(new_cell_n):
                                if glh_type == 'glhymps_1':
                                    lay_glh_lst.append(1)
                                else:
                                    lay_glh_lst.append(2)                        
                            
                            """                    
                            for j in range(glh_pct):
                                if glh_type == 'glhymps_1':
                                    lay_glh_lst.append(1)
                                else:
                                    lay_glh_lst.append(2)                    
                    
                            #   split this list into even chunks based on the lay_thk_pct value
                            lay_glhymps_lst = [lay_glh_lst[i:i + int(lay_thk_pct)] for i in range(0, len(lay_glh_lst), int(lay_thk_pct))]         
                            """
                            
                        lay_glhymps_lst = [lay_glh_lst[i:i + 1] for i in range(0, tot_cells, 1)]
                            
                    else:
                        
                        #   create a list of 100 values (1 or 2 depending on the glhymps type)
                        for i in range(len(inland_aqf_lrs)):
                            glh_pct = inland_aqf_lrs[i][0]
                            glh_type = inland_aqf_lrs[i][1]  
                            
                            #   adapt the number of cells to be assigned
                            new_cell_n = int(round((tot_cells / 100) * glh_pct, 0))
                            
                            #print(glh_type, glh_pct, new_cell_n)
                            
                            for j in range(new_cell_n):
                                if glh_type == 'glhymps_1':
                                    lay_glh_lst.append(1)
                                else:
                                    lay_glh_lst.append(2)                        
                        
                        lay_glhymps_lst = [lay_glh_lst[i:i + 1] for i in range(0, tot_cells, 1)]
                        
                    if lay_glhymps_lst != [[]] and [i for i in lay_glhymps_lst if i != []] != []:
                        last_non_empty = [z[0] for z in lay_glhymps_lst if z != []][-1]
                        
                    else:
                        glh_type_tofill = max([sublist for sublist in inland_aqf_lrs])[1]
                        if glh_type_tofill == 'glhymps_2': 
                            last_non_empty = 2
                        elif glh_type_tofill == 'glhymps_1': 
                            last_non_empty = 1            
                    
                    for x in range(len(lay_glhymps_lst)):
                        if lay_glhymps_lst[x] == []:
                            lay_glhymps_lst[x] = [last_non_empty]
                    
                    #   make sure the list above has the right amount of elements (in case the % is 33.0 for example it appends one extra value to reach length of 100)
                    if len(lay_glhymps_lst) > tot_cells:
                        lay_glhymps_lst = lay_glhymps_lst[:-1]


                    """
                    #   create a list of 100 values (1 or 2 depending on the glhymps type)
                    for i in range(len(inland_aqf_lrs)):
                        glh_pct = inland_aqf_lrs[i][0]
                        glh_type = inland_aqf_lrs[i][1]                   
                        for j in range(glh_pct):
                            if glh_type == 'glhymps_1':
                                lay_glh_lst.append(1)
                            else:
                                lay_glh_lst.append(2)                    
                    
                    #   split this list into even chunks based on the lay_thk_pct value
                    lay_glhymps_lst = [lay_glh_lst[i:i + int(lay_thk_pct)] for i in range(0, len(lay_glh_lst), int(lay_thk_pct))]
                    """
                    
                    #   make sure the list above has the right amount of elements (in case the % is 33.0 for example it appends one extra value to reach length of 100)
                    if len(lay_glhymps_lst) > tot_cells:
                        diff_idx = len(lay_glhymps_lst) - tot_cells
                        rem_idxs = len(lay_glhymps_lst) / (diff_idx + 1)
                        idx_rem = int(rem_idxs)
                        for y in range(diff_idx):
                            try:
                                del lay_glhymps_lst[idx_rem]
                                idx_rem += int(rem_idxs)
                            except IndexError:
                                pass
                            
                    #   if some layers are missing (e.g. % is 6.0) then find out how many and at equidistant locations in the list repeat the glhymps sublist
                    elif len(lay_glhymps_lst) < tot_cells:
                        miss_lay_n = tot_cells - len(lay_glhymps_lst)
                        ins_idxs = int(round(len(lay_glhymps_lst) / (miss_lay_n + 1)))
                        idx_ins = ins_idxs
                        for x in range(miss_lay_n):
                            lay_glhymps_lst.insert(idx_ins, lay_glhymps_lst[idx_ins])
                            idx_ins += ins_idxs

                    #   loop through the list and assign the values to the lrs_lst
                    glh_nr = 1
                    glh_tx = 'glhymps_1'
                    to_lst = []

                    #   go through the list and depending on the majority of values in each sublist decide which GLHYMPS value will be filled in
                    for l in range(tot_cells):
                        lay = lay_glhymps_lst[l]
                        glh_typ = max(set(lay), key = lay.count)
                        #print glh_typ
                        if not const_geo_hk_vals: 
                            if glh_typ == 1:
                                model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                            elif glh_typ == 2:
                                #model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)                        
                                hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                if hk_cell_val > hk_top_mean / 10.:
                                    model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, a] = hk_cell_val
                                else:
                                    model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, a] = hk_top_mean / 10.
                        else:
                            if glh_typ == 1:
                                model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, a] = glh_1_val
                            elif glh_typ == 2:
                                model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, a] = glh_2_val  
                           
                        if glh_typ == glh_nr:
                            to_lst.append(ibound_act_lay_idxs[col_top_lay_n + l])
                        else:
                            lrs_lst.append([glh_tx, to_lst])
                            if glh_nr == 1:
                                glh_nr = 2
                                glh_tx = 'glhymps_2'
                            else:
                                glh_nr = 1
                                glh_tx = 'glhymps_1'                    
                            to_lst = []


                #   another case
                elif a < end_col_idxs[0]:
                    
                    #   find the upper part of the model column with glhymps 2 values
                    lay_idx_upp = ibound_act_lay_idxs.index(end_lay_idxs[0])
                    tot_cells = len(ibound_act_lay_idxs) - lay_idx_upp                  
                    
                    #   get the pct of thickness per model layer
                    lay_thk_pct = round(1 / float(tot_cells) * 100., 0)
                    #   decide the layer types in the rest of the model layers based on the majority of glhymps type in that layer
                    lay_glh_lst = []                    


                    # in very deep systems the pct can be equal to 0, in that case round to one more decimal
                    if lay_thk_pct != 0.0:
                        lay_thk_pct = math.ceil(1 / float(tot_cells) * 100.)
          
                        #   create a list of 100 values (1 or 2 depending on the glhymps type)
                        for i in range(len(inland_aqf_lrs)):
                            glh_pct = inland_aqf_lrs[i][0]
                            glh_type = inland_aqf_lrs[i][1]   
        
                            #   adapt the number of cells to be assigned
                            new_cell_n = int(round((tot_cells / 100) * glh_pct, 0))
                            
                            #print(glh_type, glh_pct, new_cell_n)
                            
                            for j in range(new_cell_n):
                                if glh_type == 'glhymps_1':
                                    lay_glh_lst.append(1)
                                else:
                                    lay_glh_lst.append(2)                        
                            
                            """                    
                            for j in range(glh_pct):
                                if glh_type == 'glhymps_1':
                                    lay_glh_lst.append(1)
                                else:
                                    lay_glh_lst.append(2)                    
                    
                            #   split this list into even chunks based on the lay_thk_pct value
                            lay_glhymps_lst = [lay_glh_lst[i:i + int(lay_thk_pct)] for i in range(0, len(lay_glh_lst), int(lay_thk_pct))]         
                            """
                            
                        lay_glhymps_lst = [lay_glh_lst[i:i + 1] for i in range(0, tot_cells, 1)]
                            
                    else:
                        
                        #   create a list of 100 values (1 or 2 depending on the glhymps type)
                        for i in range(len(inland_aqf_lrs)):
                            glh_pct = inland_aqf_lrs[i][0]
                            glh_type = inland_aqf_lrs[i][1]  
                            
                            #   adapt the number of cells to be assigned
                            new_cell_n = int(round((tot_cells / 100) * glh_pct, 0))
                            
                            #print(glh_type, glh_pct, new_cell_n)
                            
                            for j in range(new_cell_n):
                                if glh_type == 'glhymps_1':
                                    lay_glh_lst.append(1)
                                else:
                                    lay_glh_lst.append(2)                        
                        
                        lay_glhymps_lst = [lay_glh_lst[i:i + 1] for i in range(0, tot_cells, 1)]
                        
                    if lay_glhymps_lst != [[]] and [i for i in lay_glhymps_lst if i != []] != []:
                        last_non_empty = [z[0] for z in lay_glhymps_lst if z != []][-1]
                        
                    else:
                        glh_type_tofill = max([sublist for sublist in inland_aqf_lrs])[1]
                        if glh_type_tofill == 'glhymps_2': 
                            last_non_empty = 2
                        elif glh_type_tofill == 'glhymps_1': 
                            last_non_empty = 1            
                    
                    for x in range(len(lay_glhymps_lst)):
                        if lay_glhymps_lst[x] == []:
                            lay_glhymps_lst[x] = [last_non_empty]
                    
                    #   make sure the list above has the right amount of elements (in case the % is 33.0 for example it appends one extra value to reach length of 100)
                    if len(lay_glhymps_lst) > tot_cells:
                        lay_glhymps_lst = lay_glhymps_lst[:-1]

                    """
                    #   create a list of 100 values (1 or 2 depending on the glhymps type)
                    for i in range(0, len(inland_aqf_lrs)):
                        glh_pct = inland_aqf_lrs[i][0]
                        glh_type = inland_aqf_lrs[i][1]  
                        for j in range(glh_pct):
                            if glh_type == 'glhymps_1':
                                lay_glh_lst.append(1)
                            else:
                                lay_glh_lst.append(2)                    
                    """
                    
                    #   if some layers are missing (e.g. % is 6.0) then find out how many and at equidistant locations in the list repeat the glhymps sublist
                    if len(lay_glh_lst) < 100 and len(lay_glh_lst) > 0:
                        lay_glh_lst = list(itertools.chain.from_iterable((itertools.repeat(i, 2) for i in lay_glh_lst)))  
                        """
                        miss_lay_n = 100 - len(lay_glh_lst)
                        ins_idxs = len(lay_glh_lst) / (miss_lay_n + 1)
                        idx_ins = ins_idxs
                        for x in xrange(miss_lay_n):
                            lay_glh_lst.insert(idx_ins, lay_glh_lst[idx_ins])
                            idx_ins += ins_idxs
                        """

                    #   split this list into even chunks based on the lay_thk_pct value
                    #lay_glhymps_lst = [lay_glh_lst[i:i + int(lay_thk_pct)] for i in range(0, len(lay_glh_lst), int(lay_thk_pct))]
                    
                    #   make sure the list above has the right amount of elements (in case the % is 33.0 for example it appends one extra value to reach length of 100)
                    if len(lay_glhymps_lst) > tot_cells:
                        diff_idx = len(lay_glhymps_lst) - tot_cells
                        rem_idxs = len(lay_glhymps_lst) / (diff_idx + 1)
                        idx_rem = int(rem_idxs)
                        for y in range(diff_idx):
                            try:
                                del lay_glhymps_lst[idx_rem]
                                idx_rem += int(rem_idxs)
                            except IndexError:
                                pass
                            
                    #   if some layers are missing (e.g. % is 6.0) then find out how many and at equidistant locations in the list repeat the glhymps sublist
                    elif len(lay_glhymps_lst) < tot_cells:
                        miss_lay_n = tot_cells - len(lay_glhymps_lst)
                        ins_idxs = len(lay_glhymps_lst) / (miss_lay_n + 1)
                        idx_ins = ins_idxs
                        try:
                            for x in range(miss_lay_n):
                                lay_glhymps_lst.insert(idx_ins, lay_glhymps_lst[idx_ins])
                                idx_ins += ins_idxs
                        # in case the miss_lay_n is not an integer (TypeError: list indices must be integers or slices, not float)
                        except TypeError:
                            for x in range(miss_lay_n):
                                lay_glhymps_lst.insert(round(idx_ins), lay_glhymps_lst[round(idx_ins)])
                                idx_ins += ins_idxs
                                
                    #   create a list of final glhymps values and change the first layers to glhymps 2 if necessary
                    final_glh_lst = [max(set(sublist), key = sublist.count) for sublist in lay_glhymps_lst]
                    for x in range(lay_idx_upp):
                        final_glh_lst.insert(0, 2)
                    
                    #   go through the list and depending on the majority of values in each sublist decide which GLHYMPS value will be filled in
                    try:
                        for l in range(len(final_glh_lst)):
                            glh_typ = final_glh_lst[l]
                            if not const_geo_hk_vals: 
                                if glh_typ == 1:
                                    model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                                elif glh_typ == 2:
                                    #model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)                            
                                    hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                    if hk_cell_val > hk_top_mean / 10.:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_cell_val
                                    else:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_top_mean / 10.
                            else:
                                if glh_typ == 1:
                                    model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_1_val
                                elif glh_typ == 2:
                                    model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_2_val  
                    except IndexError:
                        for l in range(len(ibound_act_lay_idxs)):
                            glh_typ = final_glh_lst[l]
                            if not const_geo_hk_vals: 
                                if glh_typ == 1:
                                    model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                                elif glh_typ == 2:
                                    #model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)      
                                    hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                    if hk_cell_val > hk_top_mean / 10.:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_cell_val
                                    else:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_top_mean / 10.                                    
                            else:
                                if glh_typ == 1:
                                    model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_1_val
                                elif glh_typ == 2:
                                    model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_2_val  

                    #   loop through the list and assign the values to the lrs_lst
                    glh_nr = final_glh_lst[0]
                    if glh_nr == 1:
                        glh_tx = 'glhymps_1'
                    else:
                        glh_tx = 'glhymps_2'
                    to_lst = []
                    try:
                        for y in range(len(final_glh_lst)):
                            if final_glh_lst[y] == glh_nr:
                                to_lst.append(ibound_act_lay_idxs[y])
                            else:
                                lrs_lst.append([glh_tx, to_lst])
                                if glh_nr == 1:
                                    glh_nr = 2
                                    glh_tx = 'glhymps_2'
                                else:
                                    glh_nr = 1
                                    glh_tx = 'glhymps_1'                    
                                to_lst = []
                    except IndexError:
                        for y in range(len(ibound_act_lay_idxs)):
                            if final_glh_lst[y] == glh_nr:
                                to_lst.append(ibound_act_lay_idxs[y])
                            else:
                                lrs_lst.append([glh_tx, to_lst])
                                if glh_nr == 1:
                                    glh_nr = 2
                                    glh_tx = 'glhymps_2'
                                else:
                                    glh_nr = 1
                                    glh_tx = 'glhymps_1'                    
                                to_lst = []                        

                    #   get the upp layer limit of the top glhymps layer 
                    first_lay = 0
                    try:
                        for x in range(len(final_glh_lst)):
                            if final_glh_lst[x] == 1:
                                pass
                            else:
                                first_lay = ibound_act_lay_idxs[x - 1]
                                break    
                    except IndexError:
                        for x in range(len(ibound_act_lay_idxs)):
                            if final_glh_lst[x] == 1:
                                pass
                            else:
                                first_lay = ibound_act_lay_idxs[x - 1]
                                break                           
                        
                    first_col = a
                   
                #   also if the column is after the last end of a glhymps layer adjust the procedure
                elif a > end_col_idxs[-2]:                    
                    #   try except in case there are no active cells in the ibound column (happens in the end of the model IBOUND array)
                    try:
                        #   1) decide what layers can be active in the given model column
                        col_counter = 0
                        for b in range(1, len(end_col_idxs)):
                            if a > end_col_idxs[b]:
                                col_counter += 1                                        
                
                        #       get the glhymps type of the top model layers from the inland_aqf_lst
                        glh_typ_top = inland_aqf_lrs[col_counter][1]
                        if glh_typ_top == 'glhymps_2':
                            glh_typ_top_type = 2
                        else:
                            glh_typ_top_type = 1     
                            
                        try:
                            for l in range(len(final_glh_lst)):
                                glh_typ = final_glh_lst[l]
                                if not const_geo_hk_vals: 
                                    if glh_typ == 1:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                                    elif glh_typ == 2:
                                        #model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)  
                                        hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                        if hk_cell_val > hk_top_mean / 10.:
                                            model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_cell_val
                                        else:
                                            model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_top_mean / 10.                                           
                                else:
                                    if glh_typ == 1:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_1_val
                                    elif glh_typ == 2:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_2_val  
                        except IndexError:
                            for l in range(len(ibound_act_lay_idxs)):
                                glh_typ = final_glh_lst[l]
                                if not const_geo_hk_vals: 
                                    if glh_typ == 1:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                                    elif glh_typ == 2:
                                        #model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)    
                                        hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                        if hk_cell_val > hk_top_mean / 10.:
                                            model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_cell_val
                                        else:
                                            model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_top_mean / 10.                                           
                                else:
                                    if glh_typ == 1:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_1_val
                                    elif glh_typ == 2:
                                        model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_2_val  

                        #   loop through the list and assign the values to the lrs_lst
                        for y in range(len(ibound_act_lay_idxs)):
                            if glh_typ_top_type == 1:
                                lrs_lst.append(['glhymps_1', ibound_act_lay_idxs[:]])
                            else:
                                lrs_lst.append(['glhymps_2', ibound_act_lay_idxs[:]])
                            

                    #   this happens at the end of contintental slope part of the domain
                    except (ZeroDivisionError, IndexError):
                        #print('ok', a)
                        ibound_idx_grp = list(split(ibound_act_lay_idxs, 1))

                        for w in range(0, len(ibound_idx_grp)):
                            #   + 1 because the upper layer is already filled in by GLHYMPS 2.0
                            if not const_geo_hk_vals: 
                                if inland_aqf_lrs[w][1] == 'glhymps_2':
                                    #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                                    for t in range(len(ibound_idx_grp[w])):
                                        #model_obj.hk_arr[ibound_idx_grp[w][t], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                        hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                        if hk_cell_val > hk_top_mean / 10.:
                                            model_obj.hk_arr[ibound_idx_grp[w][t], 0, a] = hk_cell_val
                                        else:
                                            model_obj.hk_arr[ibound_idx_grp[w][t], 0, a] = hk_top_mean / 10.                                           
                                else:
                                    for t in range(len(ibound_idx_grp[w])):
                                        model_obj.hk_arr[ibound_idx_grp[w][t], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4) 
                            else:
                                if inland_aqf_lrs[w][1] == 'glhymps_2':
                                    #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                                    for t in range(len(ibound_idx_grp[w])):
                                        model_obj.hk_arr[ibound_idx_grp[w][t], 0, a] = glh_2_val
                                else:
                                    for t in range(len(ibound_idx_grp[w])):
                                        model_obj.hk_arr[ibound_idx_grp[w][t], 0, a] = glh_1_val

                else:
                    #   1) decide what layers can be active in the given model column
                    col_counter = 0
                    for b in range(1, len(end_col_idxs)):
                        if a > end_col_idxs[b]:
                            col_counter += 1                                        

                    #       get the glhymps type of the top model layers from the inland_aqf_lst
                    glh_typ_top = inland_aqf_lrs[col_counter][1]
                    if glh_typ_top == 'glhymps_2':
                        glh_typ_top_type = 2
                    else:
                        glh_typ_top_type = 1

                    #   at the beginning, check if we need to change the first_lay and first_col values
                    if a == end_col_idxs[col_counter] + 1:
                        first_lay = 0
                        #   make sure there are no remnants of the other glhymps type on top (usually one layer at the end of the stretch)
                        while final_glh_lst[0] != glh_typ_top_type:
                            final_glh_lst = final_glh_lst[1:]
                        for x in range(len(final_glh_lst)):
                            #print final_glh_lst[x]
                            if final_glh_lst[x] == glh_typ_top_type:
                                pass
                            else:
                                try:
                                    first_lay = ibound_act_lay_idxs[x + 1]
                                    break                    
                                except IndexError:
                                    pass
                        first_col = a

                    last_lay, last_col = end_lay_idxs[col_counter + 1], end_col_idxs[col_counter + 1]
                
                    #   calculate the cell step based on the distance between the end columns and end layers of the top layer
                    try:
                        cell_step = (last_col - first_col) / (last_lay - first_lay)          
                        #   based on the ascent decide how many cells will be added to the bottom layer of the top glhymps layer
                        add_lays = int(round((a - first_col) / float(cell_step), 0))
                    except (UnboundLocalError, ZeroDivisionError):
                        cell_step, add_lays, first_lay = 0, 0, 0

                    end_lay_new = first_lay + add_lays
                    
                    #   it can happen during the rounding up that a value is selected that doesnt exist in the actual list
                    try:
                        upp_lay_lim = ibound_act_lay_idxs.index(end_lay_new)
                    except ValueError:
                        upp_lay_lim = ibound_act_lay_idxs.index(min(ibound_act_lay_idxs, key = lambda x : abs(x - end_lay_new)))
                        
                    tot_cells = len(ibound_act_lay_idxs) - upp_lay_lim
                    
                    #   2) divide the column based on the % fractions from the inland_aqf_lst
                    #       get the pct of thickness per model layer
                    lay_thk_pct = round(1 / float(tot_cells) * 100., 0)
                    #       decide the layer types in the rest of the model layers based on the majority of glhymps type in that layer
                    lay_glh_lst = []                                



                    # in very deep systems the pct can be equal to 0, in that case round to one more decimal
                    if lay_thk_pct != 0.0:
                        lay_thk_pct = math.ceil(1 / float(tot_cells) * 100.)
          
                        #   create a list of 100 values (1 or 2 depending on the glhymps type)
                        for i in range(len(inland_aqf_lrs)):
                            glh_pct = inland_aqf_lrs[i][0]
                            glh_type = inland_aqf_lrs[i][1]   
        
                            #   adapt the number of cells to be assigned
                            new_cell_n = int(round((tot_cells / 100) * glh_pct, 0))
                            
                            #print(glh_type, glh_pct, new_cell_n)
                            
                            for j in range(new_cell_n):
                                if glh_type == 'glhymps_1':
                                    lay_glh_lst.append(1)
                                else:
                                    lay_glh_lst.append(2)                        
                            
                            """                    
                            for j in range(glh_pct):
                                if glh_type == 'glhymps_1':
                                    lay_glh_lst.append(1)
                                else:
                                    lay_glh_lst.append(2)                    
                    
                            #   split this list into even chunks based on the lay_thk_pct value
                            lay_glhymps_lst = [lay_glh_lst[i:i + int(lay_thk_pct)] for i in range(0, len(lay_glh_lst), int(lay_thk_pct))]         
                            """
                            
                        lay_glhymps_lst = [lay_glh_lst[i:i + 1] for i in range(0, tot_cells, 1)]
                            
                    else:
                        
                        #   create a list of 100 values (1 or 2 depending on the glhymps type)
                        for i in range(len(inland_aqf_lrs)):
                            glh_pct = inland_aqf_lrs[i][0]
                            glh_type = inland_aqf_lrs[i][1]  
                            
                            #   adapt the number of cells to be assigned
                            new_cell_n = int(round((tot_cells / 100) * glh_pct, 0))
                            
                            #print(glh_type, glh_pct, new_cell_n)
                            
                            for j in range(new_cell_n):
                                if glh_type == 'glhymps_1':
                                    lay_glh_lst.append(1)
                                else:
                                    lay_glh_lst.append(2)                        
                        
                        lay_glhymps_lst = [lay_glh_lst[i:i + 1] for i in range(0, tot_cells, 1)]
                        
                    if lay_glhymps_lst != [[]] and [i for i in lay_glhymps_lst if i != []] != []:
                        last_non_empty = [z[0] for z in lay_glhymps_lst if z != []][-1]
                        
                    else:
                        glh_type_tofill = max([sublist for sublist in inland_aqf_lrs])[1]
                        if glh_type_tofill == 'glhymps_2': 
                            last_non_empty = 2
                        elif glh_type_tofill == 'glhymps_1': 
                            last_non_empty = 1            
                    
                    for x in range(len(lay_glhymps_lst)):
                        if lay_glhymps_lst[x] == []:
                            lay_glhymps_lst[x] = [last_non_empty]
                    
                    #   make sure the list above has the right amount of elements (in case the % is 33.0 for example it appends one extra value to reach length of 100)
                    if len(lay_glhymps_lst) > tot_cells:
                        lay_glhymps_lst = lay_glhymps_lst[:-1]
                        
                        
                        
                    """
                    #   create a list of 100 values (1 or 2 depending on the glhymps type)
                    for i in range(col_counter + 1, len(inland_aqf_lrs)):
                        glh_pct = inland_aqf_lrs[i][0]
                        glh_type = inland_aqf_lrs[i][1]  
                        for j in range(glh_pct):
                            if glh_type == 'glhymps_1':
                                lay_glh_lst.append(1)
                            else:
                                lay_glh_lst.append(2)                    
                    """
                    
                    #   if some layers are missing (e.g. % is 6.0) then find out how many and at equidistant locations in the list repeat the glhymps sublist
                    if len(lay_glh_lst) < 100 and len(lay_glh_lst) > 0:
                        while len(lay_glh_lst) < 100:
                            lay_glh_lst = list(itertools.chain.from_iterable((itertools.repeat(i, 2) for i in lay_glh_lst)))  
                            #print(len(lay_glh_lst), i)
                            
                            """
                            miss_lay_n = 100 - len(lay_glh_lst)
                            ins_idxs = len(lay_glh_lst) / (miss_lay_n + 1)
                            idx_ins = ins_idxs
                            for x in xrange(miss_lay_n):
                                lay_glh_lst.insert(idx_ins, lay_glh_lst[idx_ins])
                                idx_ins += ins_idxs
                            """
                        
                    #   split this list into even chunks based on the lay_thk_pct value
                    #lay_glhymps_lst = [lay_glh_lst[i:i + int(lay_thk_pct)] for i in range(0, len(lay_glh_lst), int(lay_thk_pct))]
                    
                    #   make sure the list above has the right amount of elements (in case the % is 33.0 for example it appends one extra value to reach length of 100)
                    if len(lay_glhymps_lst) > tot_cells:
                        while len(lay_glhymps_lst) > tot_cells:
                            diff_idx = len(lay_glhymps_lst) - tot_cells
                            rem_idxs = len(lay_glhymps_lst) / (diff_idx + 1)
                            idx_rem = int(rem_idxs)
                            for y in range(diff_idx):
                                try:
                                    del lay_glhymps_lst[idx_rem]
                                    idx_rem += int(rem_idxs)
                                except IndexError:
                                    pass
                            
                    #   if some layers are missing (e.g. % is 6.0) then find out how many and at equidistant locations in the list repeat the glhymps sublist
                    elif len(lay_glhymps_lst) < tot_cells:
                        while len(lay_glhymps_lst) < tot_cells: 
                            miss_lay_n = tot_cells - len(lay_glhymps_lst)
                            ins_idxs = len(lay_glhymps_lst) / (miss_lay_n + 1)
                            idx_ins = int(ins_idxs)
                            for x in range(miss_lay_n):
                                lay_glhymps_lst.insert(idx_ins, lay_glhymps_lst[idx_ins])
                                idx_ins += int(ins_idxs)

                    #   create a list of final glhymps values and change the first layers to glhymps 2 if necessary
                    final_glh_lst = [max(set(sublist), key = sublist.count) for sublist in lay_glhymps_lst]

                    for x in range(upp_lay_lim):
                        final_glh_lst.insert(0, glh_typ_top_type)

                    for l in range(len(final_glh_lst)):
                        glh_typ = final_glh_lst[l]
                        #print glh_typ
                        if not const_geo_hk_vals: 
                            if glh_typ == 1:
                                model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                            elif glh_typ == 2:
                                #model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)    
                                hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                if hk_cell_val > hk_top_mean / 10.:
                                    model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_cell_val
                                else:
                                    model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = hk_top_mean / 10.          
                        else:
                            if glh_typ == 1:
                                model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_1_val
                            elif glh_typ == 2:
                                model_obj.hk_arr[ibound_act_lay_idxs[l], 0, a] = glh_2_val  

                    #   loop through the list and assign the values to the lrs_lst
                    glh_nr = final_glh_lst[0]
                    if glh_nr == 1:
                        glh_tx = 'glhymps_1'
                    else:
                        glh_tx = 'glhymps_2'
                    to_lst = []
                    for y in range(len(final_glh_lst)):
                        if final_glh_lst[y] == glh_nr:
                            to_lst.append(ibound_act_lay_idxs[y])
                        else:
                            lrs_lst.append([glh_tx, to_lst])
                            if glh_nr == 1:
                                glh_nr = 2
                                glh_tx = 'glhymps_2'
                            else:
                                glh_nr = 1
                                glh_tx = 'glhymps_1'                    
                            to_lst = []  
                            
                offshore_lrs_lst.append([a, lrs_lst])


    elif sed_flux == 'medium' or sed_flux == 'high':    

        #   first calculate the average slopes of the cont. shelf and cont. slope, as the elevation difference between the shelf edge and the
        #   elevation at the coastline
        avg_sl_shelf = round(100 * (model_obj.top_elev[cst_idx_geo] - model_obj.top_elev[shelf_edge_idx]) / ((shelf_edge_idx - cst_idx_geo) * 100.), 2)
        avg_sl_slope = round(100 * (model_obj.top_elev[shelf_edge_idx] - model_obj.top_elev[-1]) / ((model_obj.ncol - shelf_edge_idx) * 100.), 2)
            
        glhymps_2_top = False
        
        if sed_flux == 'medium':
            shelf_angle = 85.   #   this needs to be specified based on the sed_flux being either medium or high
        elif sed_flux == 'high':
            shelf_angle = 87.5
        
        """
        FIRST DEAL WITH THE AREA BETWEEN THE COASTLINE AND THE BEGINNING OF THE SHELF EDGES AREA
        """
        #   get the top points of each layer (not the top one..), different for each scenario
        cst_lay_pts = []
        
        #   calculate the thickness of the top glhymps_2 layer at the coastline, to be expanded into the offshore domain
        if len(soilgrids_thk_lst) > 0:
            col_top_lay_n = int(round(soilgrids_thk_lst[-1] / model_obj.del_lay))   
        else:
            col_top_lay_n = 0
        #   get the active layer indexes at the coastline
        cst_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, cst_offset_idx].tolist()) if x == 1]
        
        #   get the index of the shelf edge and also a list of active layers at that model column
        shelf_edge_idx = model_obj.top_elev.index(cont_shelf_edge[1]) 
        shelf_edge_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edge_idx].tolist()) if x == 1]   
        
        #   divide the area between the shelf edge and model_obj bottom into the tot_layers identical amount of cells based on the ratio 
        tot_layers = len(inland_aqf_lrs)
        tot_cells = len(cst_ibound_act_lay_idxs[col_top_lay_n:])
        #   append the last layer of the top glh layer, if possible - can happen that the total thickness is lower than the glh top
        #   layer thickness, in that case append an empty list
        try:
            cst_lay_pts.append(cst_ibound_act_lay_idxs[col_top_lay_n])
        except IndexError:
            cst_lay_pts.append([cst_ibound_act_lay_idxs[-1]])
        cell_cnt_lst = []
        geo_lay_thk_start = col_top_lay_n
        
        #   if there is no glhyps_2_top layer defined then still continue the more permeable sediment into the offshore domain, but only
        #   till the last inland model layer where this sediment layer occurs 
        if glhymps_2_top is False:
            try:
                last_glh_2_top_lay = cst_ibound_act_lay_idxs[col_top_lay_n]
            except IndexError:
                last_glh_2_top_lay = cst_ibound_act_lay_idxs[-1]
        
        #   loop through the individual glhymps 1 and 2 layers and find the corresponding bottom model layers
        for q in range(tot_layers):
            geo_lay_thk = int(round((inland_aqf_lrs[q][0] * tot_cells) / 100, 0))   #   thickness of glhymps layer expressed in model layers
            #   check if the calculated thickness is 0 then just go to the next layer
            if geo_lay_thk == 0:
                pass
            else:
                if q == len(inland_aqf_lrs) - 1:#   len(inland_aqf_lrs) - 1 because of python counting..
                    pass
                else:
                    #   if the soilgrids layer is too thick and there are not enough cells left to fill with all the other glhymps layers then skip them
                    try:
                        cst_lay_pts.append(cst_ibound_act_lay_idxs[geo_lay_thk_start + geo_lay_thk])
                        geo_lay_thk_start += geo_lay_thk    
                    except IndexError:
                        pass
                        
        #   calculate the angle between the aquifer bottom and the perpendicular line at the shelf edge
        aqf_bot_y_diff = round(model_obj.bot_elev[cst_idx_geo] - model_obj.bot_elev[shelf_edge_idx], 2)
        aqf_bot_x_diff = round((shelf_edge_idx - cst_idx_geo) * 100., 2)
        #   calculate the necessary angles for estimating the position of the past shelf edges
        shelf_edge_bot_angle = 90. - round(math.degrees(math.atan(aqf_bot_y_diff / aqf_bot_x_diff)), 2)
        cst_angle = round(180. - shelf_angle - shelf_edge_bot_angle, 2)   # the missing 3rd angle
        #   also calculate the necessary distances and ratios of the model part where the past shelf edges are located
        len_shelf = len(shelf_edge_ibound_act_lay_idxs) * 10.
        sin_ratio = len_shelf / math.sin(math.radians(cst_angle))
        len_bot = round(sin_ratio * math.sin(math.radians(shelf_angle)), 2)
        len_top = round(sin_ratio * math.sin(math.radians(shelf_edge_bot_angle)), 2)
        
        #   split the len_top into tot_layers parts equally and get the [layer, column] tuple for each point
        #       first calculate the perpendicular distance of the edge from the shelf edge
        dist_from_edge = round(math.sin(math.radians(shelf_angle)) * len_top, 2)
        dist_in_cols = int(round(dist_from_edge / 100.))
        
        #   check that the dist in cols between the shelf edge and coast is lower, if not set it
        if dist_in_cols > shelf_edge_idx - cst_idx_geo:
            dist_in_cols = shelf_edge_idx - cst_idx_geo
        
        #   for each glhymps layer calculate the column position of its shelf edge 
        shelf_edges_col_idx = []
        pct_start = inland_aqf_lrs[0][0]
        for i in range(len(inland_aqf_lrs) - 1):
            dist_col = (dist_in_cols * pct_start) / 100.
            shelf_edges_col_idx.append(shelf_edge_idx - int(round(dist_col, 0)))
            pct_start += inland_aqf_lrs[i + 1][0]
        
        #   this changes if the glhymps 2 layer is extended over the whole shelf domain, lay top designates the upper model layer from which we start
        #   to calculate the thickness of each glhymps layer 
        if glhymps_2_top is True:
            lay_top = shelf_edge_ibound_act_lay_idxs[col_top_lay_n] 
        else:
            if last_glh_2_top_lay not in shelf_edge_ibound_act_lay_idxs:
                lay_top = shelf_edge_ibound_act_lay_idxs[0]
            else:
                lay_top = shelf_edge_ibound_act_lay_idxs[shelf_edge_ibound_act_lay_idxs.index(last_glh_2_top_lay)]
        #   same for the bottom extent in terms of model layers
        lay_bot = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edge_idx - dist_in_cols].tolist()) if x == 1][-1]   
        
        #   next, calculate the shelf edges layer bottoms in terms of model layers
        shelf_edges_lay_idx = []
        pct_start = inland_aqf_lrs[0][0]
        for i in range(len(inland_aqf_lrs) - 1):
            dist_lay = ((lay_bot - lay_top) * pct_start) / 100.
            shelf_edges_lay_idx.append(lay_top + int(round(dist_lay, 0)))
            pct_start += inland_aqf_lrs[i + 1][0]        
        
        #   find the end of the layer based on the slope of the continental slope
        line_A_lst, line_B_lst, cst_point_lst = [], [], []
        line_thk_end = ((0.0, model_obj.bot_elev[cst_offset_idx]), (round(model_obj.x_end - model_obj.cst_offset_plot, 1), model_obj.bot_elev[-1]))   
        
        for i in range(len(shelf_edges_col_idx)):
            x_coord = round(model_obj.x_start - model_obj.cst_offset_plot + shelf_edges_col_idx[i] / 10., 1)# + 0.05
            y_coord = model_obj.top - shelf_edges_lay_idx[i] * 10.# + 5.
            #  try to calculate the y_coordinate at the coastline, in case the thickness is to shallow for all the layers to be stacked there
            #  assign the total bottom of the model as the top of these layers, they will then link with the shelf edges later on
            #y_coast = model_obj.top - (cst_lay_pts[i + 1]) * 10 
            try:
                y_coast = model_obj.top - (cst_lay_pts[i + 1]) * 10
            except IndexError:
                # TypeError: unsupported operand type(s) for -: 'float' and 'list' because clay_lay_pts = [[8]] (for example)
                try:
                    y_coast = model_obj.top - (cst_lay_pts[-1]) * 10 
                except TypeError:
                    y_coast = model_obj.top - (cst_lay_pts[-1][0]) * 10 
            #   calculate the y_coordinate of the line at the end of the model_obj domain
            y_coord_end = round(avg_sl_slope * (model_obj.x_end - model_obj.cst_offset_plot - x_coord) * 10., 1)
            line_A_lst.append((round(x_coord, 1), round(y_coord, 1)))
            line_B_lst.append((round(model_obj.x_end - model_obj.cst_offset_plot, 1), round(y_coord - y_coord_end, 1)))
            cst_point_lst.append((0.0, round(y_coast)))
        
        bot_point_lst = []
        for j in range(len(line_A_lst)):
            bot_pt = line_intersection((line_A_lst[j], line_B_lst[j]), line_thk_end)
            bot_point_lst.append([round(bot_pt[0], 1), round(bot_pt[1], 1)])
        
        
        #   if the shelf area is thinner than the end of the model domain it means that the bedrock line and the layer lines will never intersect offshore,
        #   instead they would intersect inland. Because of that, cross-section points will be calculated individually at the end of model domain
        shelf_active_layers = ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edges_col_idx[-1]].tolist()) if x == 1]       
        end_active_layers = ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, model_obj.ibound_arr.shape[-1] - 1].tolist()) if x == 1] 
        
        if len(end_active_layers) >= len(shelf_active_layers):
        
            #   check that the position of the intersection is in the offshore direction from the shelf break, if not then change the line_thk_end coordinates
            if bot_point_lst[-1][0] < line_A_lst[-1][0]:
                bot_point_lst = []
                
                x_end = round(model_obj.x_start + round((model_obj.ibound_arr.shape[-1] * (model_obj.del_col / 1000.)), 2), 1)
                y_top_end = model_obj.top_elev[-1]
                y_bot_end = model_obj.lay_elev[model_obj.ibound_arr.shape[0]]
                
                line_B_lst = []
                tot_lay_thk_cnt = 0
                tot_thk = abs(y_bot_end) - abs(y_top_end)
                for h in range(len(inland_aqf_lrs) - 1):
                    #   calculate the depth of the layer end
                    pct_incr = round((tot_thk * inland_aqf_lrs[h][0]) / 100., 1)
                    tot_lay_thk_cnt += pct_incr
                    line_B_lst.append((x_end, round(y_top_end - tot_lay_thk_cnt, 1)))
                    
                line_thk_end = ((x_end, y_top_end), (x_end, y_bot_end))
                for j in range(len(line_A_lst)):
                    bot_pt = line_intersection((line_A_lst[j], line_B_lst[j]), line_thk_end)
                    bot_point_lst.append([round(bot_pt[0], 1), round(bot_pt[1], 1)])            

        shelf_edge_0_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edges_col_idx[-1]].tolist()) if x == 1]   
        
        #   find the limits for each of the geological layers 
        if glhymps_2_top is True:
            idx_1 = shelf_edge_0_ibound_act_lay_idxs[col_top_lay_n]  
        else:
            idx_1 = shelf_edge_0_ibound_act_lay_idxs[0]
        
        end_lay_idxs =  [] #end_lay_idxs =  [idx_1]
        act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[idx_1, 0, :].tolist()) if x == 1]
        end_col_idxs = []#[act_cell_lst[-1]]                
        for b in range(tot_layers - 1):
            if glhymps_2_top is True:
                lay_idx = idx_1 + int(round((inland_aqf_lrs[b][0] * (model_obj.ibound_arr.shape[0] - shelf_edge_0_ibound_act_lay_idxs[col_top_lay_n])) / 100 , 0))          
            else:
                lay_idx = idx_1 + int(round((inland_aqf_lrs[b][0] * (model_obj.ibound_arr.shape[0] - shelf_edge_0_ibound_act_lay_idxs[0])) / 100 , 0))               
                            
            end_lay_idxs.append(lay_idx)
            try:
                act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[lay_idx - 1, 0, :].tolist()) if x == 1]
            except IndexError:
                act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[-1, 0, :].tolist()) if x == 1]
            
            #   if the list is empty then just paste the same value as the last one
            if act_cell_lst != []:
                end_col_idxs.append(act_cell_lst[-1])
                idx_1 = lay_idx
            else:
                end_col_idxs.append(end_col_idxs[-1])
                idx_1 = lay_idx                

            
        #   append the last active layer and column as the end of the last layer
        end_lay_idxs.append(model_obj.ibound_arr.shape[0])
        end_col_idxs.append(int(model_obj.ncol) - 1)
                
        
        #   choose the shelf_edge index, has to be larger than the coastal index
        #idx_shlf_lst = next(x[0] for x in enumerate(shelf_edges_col_idx[::-1]) if x[1] > cst_idx_geo)
        #first_shlf_idx = shelf_edges_col_idx[::-1][idx_shlf_lst]
        
        shelf_edges_col_idx = [i for i in shelf_edges_col_idx if i > cst_idx_geo]
        
        #   if the length of the shelf_edges_col_idx is empty it means that all the offshore is too shallow to identify a contintental shelf - so fill in the last index instead
        if len(shelf_edges_col_idx) == 0:
            shelf_edges_col_idx = [model_obj.ibound_arr.shape[-1] - 1]
        
        
        #   fill in the model domain between the coastline and the beginning of the shelf edges area
        for a in range(cst_idx_geo, shelf_edges_col_idx[-1]):    
            
            #   check if there are any ibound_act_lay_idxs indexes below the upper aquifer part
            #   fill in the part till the shelf edge with the constant thickness, after that point fill in straight line
            #   the upper part of the offshore aquifer domain

            ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, a].tolist()) if x == 1]
            lrs_lst = []  
            
            #   check if there are still cells from the upper glhymps 2 layer in the active cells offshore, if yes cut the list
            h = 0
            if last_glh_2_top_lay in ibound_act_lay_idxs:
                while last_glh_2_top_lay > ibound_act_lay_idxs[h]:
                    
                    hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                    if hk_cell_val > hk_top_mean / 10.:
                        model_obj.hk_arr[ibound_act_lay_idxs[h], 0, a] = hk_cell_val
                    else:
                        model_obj.hk_arr[ibound_act_lay_idxs[h], 0, a] = hk_top_mean / 10.
                    h += 1
                        
            ibound_act_lay_idxs = ibound_act_lay_idxs[h:]
            
            #   also check if the last glh 2 top layer index is actually higher than the last element of the active layer, then the full column is glh 2
            if last_glh_2_top_lay >= ibound_act_lay_idxs[-1]:
                for i in range(len(ibound_act_lay_idxs)):

                    hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                    if hk_cell_val > hk_top_mean / 10.:
                        model_obj.hk_arr[ibound_act_lay_idxs[i], 0, a] = hk_cell_val
                    else:
                        model_obj.hk_arr[ibound_act_lay_idxs[i], 0, a] = hk_top_mean / 10.
                    h += 1
                
                    #model_obj.hk_arr[ibound_act_lay_idxs[i], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)                
                ibound_act_lay_idxs = []
                continue
                    
            #   go shelf edge by edge and fill in the columns based on the line crossing
            pts_coord_start = cst_point_lst[:]
            ptc_coord_end = line_A_lst[:]
                    
            #   get the line of the column middle
            col_top_y = model_obj.botm[ibound_act_lay_idxs[0] - 1] # round(model_obj.top - (ibound_act_lay_idxs[0]) * 10 , 1)
            try:
                col_bot_y = model_obj.botm[ibound_act_lay_idxs[-1]] # round(model_obj.top - (ibound_act_lay_idxs[-1] + 1) * 10, 1) 
            except IndexError:
                col_bot_y = model_obj.botm[-1] # round(model_obj.top - (ibound_act_lay_idxs[-1] + 1) * 10, 1) 
            col_x =  round(model_obj.x_start - model_obj.cst_offset_plot + a / 10., 1)# + 0.05               
                    
            #   get the line crossing for each of the layer lines
            line_cross_lst, idx_cross_lst = [], []
            for x in range(len(pts_coord_start)):
                line_cross = line_intersection((pts_coord_start[x], ptc_coord_end[x]), ((col_x, col_top_y), (col_x, col_bot_y)))
                #print(round(line_cross[0], 1), model_obj.top - round(line_cross[1], 0))
                idx_cross_lst.append(int((model_obj.top - round(line_cross[1], 0)) / 10.))    
            #   reverse the list to start from the top layer 
            #idx_cross_lst = idx_cross_lst[::-1]    
            
            add_lay_lst = [i for i in idx_cross_lst if i < ibound_act_lay_idxs[0]]   #   add layers (empty) on top if necessary  - sometimes there is no intersection
            #   remove values from the list that are not indexes of active cells
            idx_cross_lst = [i for i in idx_cross_lst if i in ibound_act_lay_idxs]
            
            
            if len(idx_cross_lst) == 0:
                idx_cross_lst = [ibound_act_lay_idxs[x:x+len(pts_coord_start)] for x in range(0, len(ibound_act_lay_idxs), len(pts_coord_start))]
                idx_cross_lst = idx_cross_lst[0]
            
            
            #   fill in the HK array based on the list of indexes from previous step!
            grp_lay_idx = [ibound_act_lay_idxs[:ibound_act_lay_idxs.index(idx_cross_lst[0])]]        
            
            for y in range(1, len(idx_cross_lst)):
                grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[y - 1]) : ibound_act_lay_idxs.index(idx_cross_lst[y])])
            #   append the bottom layer till the end of the active model_obj domain
            grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[-1]):])
            
            #   add empty lists on top of the system if necessary
            for z in range(len(add_lay_lst)):
                grp_lay_idx.insert(0, [])
            
            #   make sure there are no duplicates or empty lists in the grp_lay_idx
            #grp_lay_idx = [x for x in grp_lay_idx if x != []]
            grp_lay_idx_sorted = []
            for sublist in grp_lay_idx:
                #   remove any empty list from the sorted one
                #if sublist not in grp_lay_idx_sorted or sublist != []:
                if sublist not in grp_lay_idx_sorted and len(sublist) != 0:
                    #print(sublist, len(sublist))
                    grp_lay_idx_sorted.append(sublist)
            

            #   check if the upper glhymps 2 layer (from soilgrids) is still located within the offshore model domain, if so
            #   then assign the right property to those cells and cut the respective sublist to be filled in by the other layers
            for w in range(len(grp_lay_idx_sorted)):
                #   + 1 because the upper layer is already filled in by GLHYMPS 2.0
                if inland_aqf_lrs[w][1] == 'glhymps_2':
                    lrs_lst.append(['glhymps_2', grp_lay_idx_sorted[w]])
                    #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                    for t in range(len(grp_lay_idx_sorted[w])):
                        hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                        if hk_cell_val > hk_top_mean / 10.:
                            model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, a] = hk_cell_val
                        else:
                            model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, a] = hk_top_mean / 10.                        
                        #model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                else:
                    lrs_lst.append(['glhymps_1', grp_lay_idx_sorted[w]])
                    for t in range(len(grp_lay_idx_sorted[w])):
                        model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)           
        
            """
            upp_lim = 0
            
            #   loop through all the groups and assign the right GLHYMPS values            
            for w in range(len(cell_cnt_lst)):
                #   + 1 because the upper layer is already filled in by GLHYMPS 2.0
                if inland_aqf_lrs[w][1] == 'glhymps_2':
                    lrs_lst.append(['glhymps_2', cell_cnt_lst[w]])
                    #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                    for t in range(len(cell_cnt_lst[w])):
                        model_obj.hk_arr[cell_cnt_lst[w][t], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                else:
                    lrs_lst.append(['glhymps_1', cell_cnt_lst[w]])
                    for t in range(len(cell_cnt_lst[w])):
                        model_obj.hk_arr[cell_cnt_lst[w][t], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)    
            """
            offshore_lrs_lst.append([a, lrs_lst])  
        
        """
        SECOND, DEAL WITH THE AREA WITH THE SHELF EDGES
        """
        
        edge_count = 1   
        
        shelf_edge_last_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edges_col_idx[0]].tolist()) if x == 1]   

        #   this changes if the glhymps 2 layer is extended over the whole shelf domain, lay top designates the upper model layer from which we start
        #   to calculate the thickness of each glhymps layer 
        if glhymps_2_top is True:
            lay_top = shelf_edge_last_ibound_act_lay_idxs[col_top_lay_n] 
        else:
            if last_glh_2_top_lay not in shelf_edge_last_ibound_act_lay_idxs:
                lay_top = shelf_edge_last_ibound_act_lay_idxs[0]
            else:
                lay_top = shelf_edge_last_ibound_act_lay_idxs[shelf_edge_last_ibound_act_lay_idxs.index(last_glh_2_top_lay)]

        #   find the limits for each of the geological layers 
        if lay_top <= last_glh_2_top_lay:
            idx_1 = last_glh_2_top_lay # shelf_edge_last_ibound_act_lay_idxs.index(last_glh_2_top_lay)
        else:
            idx_1 = shelf_edge_0_ibound_act_lay_idxs[0] 
            
        #end_lay_idxs =  []#[idx_1]
        act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[idx_1, 0, :].tolist()) if x == 1]
        end_col_idxs = []#[act_cell_lst[-1]]                
        for e in range(tot_layers - 1): 
            lay_idx = idx_1 + int(round((inland_aqf_lrs[e][0] * (model_obj.ibound_arr.shape[0] - shelf_edge_last_ibound_act_lay_idxs[0])) / 100, 0))
            #end_lay_idxs.append(lay_idx)
            try:
                act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[lay_idx - 1, 0, :].tolist()) if x == 1]
                end_col_idxs.append(act_cell_lst[-1])
            except IndexError:
                act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[-1, 0, :].tolist()) if x == 1]
                #end_col_idxs.append(int(model_obj.ncol) - 1)
                end_col_idxs.append(int(model_obj.ncol))
            idx_1 = lay_idx
        #   append the last active layer and column as the end of the last layer
        #.append(model_obj.nlay)
        #end_col_idxs.append(int(model_obj.ncol) - 1)
        end_col_idxs.append(int(model_obj.ncol))
        
        for b in range(shelf_edges_col_idx[-1], shelf_edge_idx + 1):#shelf_edges_col_idx[0]):
        #for b in range(first_shlf_idx, shelf_edge_idx + 1):    
            #   get the active layers in the column
            ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, b].tolist()) if x == 1]
            lrs_lst = []  
            
            #   check if there are still cells from the upper glhymps 2 layer in the active cells offshore, if yes cut the list
            h = 0
            if last_glh_2_top_lay in ibound_act_lay_idxs:
                while last_glh_2_top_lay > ibound_act_lay_idxs[h]:
                    hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                    if hk_cell_val > hk_top_mean / 10.:
                        model_obj.hk_arr[ibound_act_lay_idxs[h], 0, b] = hk_cell_val
                    else:
                        model_obj.hk_arr[ibound_act_lay_idxs[h], 0, b] = hk_top_mean / 10.                       
                    #model_obj.hk_arr[ibound_act_lay_idxs[h], 0, b] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                    h += 1
            ibound_act_lay_idxs = ibound_act_lay_idxs[h:]            
            
            #   go shelf edge by edge and fill in the columns based on the line crossing
            pts_coord_start = cst_point_lst[:-edge_count] + line_A_lst[-edge_count:]
            ptc_coord_end = line_A_lst[:-edge_count] + line_B_lst[-edge_count:]
                    
            #   get the line of the column middle
            col_top_y = round(model_obj.top - (ibound_act_lay_idxs[0] + 1) * 10 , 1)
            col_bot_y = round(model_obj.top - (ibound_act_lay_idxs[-1] + 1) * 10, 1) 
            col_x = round(model_obj.x_start - model_obj.cst_offset_plot + b / 10., 1)# + 0.05               
                    
            #   get the line crossing for each of the layer lines
            line_cross_lst, idx_cross_lst = [], []
            for x in range(len(pts_coord_start)):
                line_cross = line_intersection((pts_coord_start[x], ptc_coord_end[x]), ((col_x, col_top_y), (col_x, col_bot_y)))
                #print(round(line_cross[0], 1), model_obj.top - round(line_cross[1], 0))
                idx_cross_lst.append(int((model_obj.top - round(line_cross[1], 0)) / 10.))    
            #   reverse the list to start from the top layer 
            #idx_cross_lst = idx_cross_lst[::-1]    
        
            #   check that the indexes are not higher than the extend of the active cells
            for s in range(len(idx_cross_lst)):
                if idx_cross_lst[s] > ibound_act_lay_idxs[-1]:
                    idx_cross_lst[s] = ibound_act_lay_idxs[-1]
        
            #   fill in the HK array based on the list of indexes from previous step!
            try:
                grp_lay_idx = [ibound_act_lay_idxs[:ibound_act_lay_idxs.index(idx_cross_lst[0])]] 
            #   it can happen that the cross-section happens in the inactive zone above the ocean bottom, in that case fill in empty list
            except ValueError:
                grp_lay_idx = [[]] 
                
            for y in range(1, len(idx_cross_lst)):
                try:
                    grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[y - 1]) : ibound_act_lay_idxs.index(idx_cross_lst[y])])
                #   same as above, the first index will become the first active cell in the ibound_act_lay_idx
                except ValueError:
                    try:
                        grp_lay_idx.append(ibound_act_lay_idxs[: ibound_act_lay_idxs.index(idx_cross_lst[y])])
                    except ValueError:
                        grp_lay_idx.append(ibound_act_lay_idxs[: ibound_act_lay_idxs[0]])
            #   append the bottom layer till the end of the active model_obj domain
            grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[-1]):])
            
            #   make sure there are no duplicates or empty lists in the grp_lay_idx
            #grp_lay_idx = [x for x in grp_lay_idx if x != []]
            grp_lay_idx_sorted = []
            for sublist in grp_lay_idx:
                if sublist not in grp_lay_idx_sorted or sublist == []:
                    grp_lay_idx_sorted.append(sublist)
            
            #   check if the upper glhymps 2 layer (from soilgrids) is still located within the offshore model domain, if so
            #   then assign the right property to those cells and cut the respective sublist to be filled in by the other layers
            
            
            for w in range(len(grp_lay_idx_sorted)):
                #   + 1 because the upper layer is already filled in by GLHYMPS 2.0
                if inland_aqf_lrs[w][1] == 'glhymps_2':
                    lrs_lst.append(['glhymps_2', grp_lay_idx_sorted[w]])
                    #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                    for t in range(len(grp_lay_idx_sorted[w])):
                        #model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, b] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                        hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                        if hk_cell_val > hk_top_mean / 10.:
                            model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, b] = hk_cell_val
                        else:
                            model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, b] = hk_top_mean / 10.                                 
                else:
                    lrs_lst.append(['glhymps_1', grp_lay_idx_sorted[w]])
                    for t in range(len(grp_lay_idx_sorted[w])):
                        model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, b] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)           
        
            upp_lim = 0
            
            """
            
            for s in range(col_top_lay_n):
                try:
                    if ibound_act_lay_idxs[s] < last_top_glhymps_2_lay: 
                        model_obj.hk_arr[ibound_act_lay_idxs[s], 0, b] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                        upp_lim += 1
                        last_glhymps_2_top_col = a
                except IndexError:
                    pass  
            """
            
            #   check if the first sublist of the lrs_lst is glhymps_1, if yes check if there are indexes lower or equal to the upp_lim
            if lrs_lst[0][0] == 'glhymps_1':
                #   check if the upp_lim index is in the list of indexes
                try:
                    split_idx = lrs_lst[0][1].index(ibound_act_lay_idxs[upp_lim])                        
                    lrs_lst.insert(0, ['glhymps_2', lrs_lst[0][1][:split_idx]])
                    lrs_lst[1][1] = lrs_lst[1][1][split_idx:]
                except ValueError:
                    pass
                
            """
            if split_idx:
                lrs_lst.insert(0, ['glhymps_2', lrs_lst[0][1][:split_idx]])
                lrs_lst[1][1] = lrs_lst[1][1][split_idx:]
            """
            
            edge_count = 1
            #   check the counter
            for c in range(len(shelf_edges_col_idx) - 1):
                if b > shelf_edges_col_idx[c]:
                    edge_count += 1      
        
            offshore_lrs_lst.append([b, lrs_lst])     
        
        
        """
        IN THE LAST STEP, FILL IN THE REST OF THE OFFSHORE MODEL DOMAIN
        """
        
        #   Instead of just filling in between the lines and thus creating a disproporiantely large upper sediment zone we choose a different approach
        #   that takes into account the % specified for each of the glhymps layers. 
        #   First we need to find the ending points of each glhymps bottom layer (excpet the current last layer) based on the specified % thickness.
        #   To do that, we need to define the starting and ending point for the loop. First we start with the last column from the shelf edges zone and
        #   choose the first column from the offshore end of these glhymps layers as calculated earlier based on the average slope of the continental slope.
        
        #   define the starting and ending indexes
        st_idx = shelf_edges_col_idx[0]
        end_idx = end_col_idxs[0]
        
        #   create a counter for how many layers should be left out from the proportion calculation
        lrs_cnt = 1                     #   the counter
        lrs = tot_layers - lrs_cnt      #   total  number of layers to take into account
        
        #   create lists of new starting and ending points
        end_col_layer_pts = []
        offshore_aqf_thk_lst = []
        
        #   find the limits for each of the geological layers 
        #   create a list of columns designating the end of each glhymps layer in the offshore domain based on the average slope of the continental slope
        end_col_idx_lst = [int((abs(model_obj.x_start - model_obj.cst_offset_plot) + i[0]) * 10.) for i in bot_point_lst[::-1]]
        
        for t in range(len(end_col_idx_lst) - 1):
            #print(end_col_idx_lst[t])
            lrs = tot_layers - lrs_cnt 
            aqf_ratio_lst = []
            
            #   get the list of all active cells in the given column
            try:
                ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, end_col_idx_lst[t]].tolist()) if x == 1]    
                col_x =  round(model_obj.x_start - model_obj.cst_offset_plot + end_col_idx_lst[t] / 10., 1)# + 0.05   
            except IndexError:
                ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, -1].tolist()) if x == 1]    
                col_x =  round(model_obj.x_start - model_obj.cst_offset_plot + end_col_idx_lst[t] / 10., 1)# + 0.05                   
            
            #   next, change the % of the glhymps layer by distributing the % of the layer(s) not included over the rest of them
            active_lrs = inland_aqf_lrs[:lrs]
            pct_left = sum(i[0] for i in inland_aqf_lrs[lrs:])    
            tot_assigned, pct_ratio = 0, int(round(pct_left / len(active_lrs)))
            for j in range(len(active_lrs)):
                if j != len(active_lrs) - 1:
                    aqf_ratio_lst.append([active_lrs[j][0] + pct_ratio, active_lrs[j][1]])
                    tot_assigned += pct_ratio
                else:
                    aqf_ratio_lst.append([active_lrs[j][0] + (pct_left - tot_assigned), active_lrs[j][1]])
        
            #   find the bottom of each glhymps layer in model layers
            end_lay_idxs =  []     
            idx_1 = ibound_act_lay_idxs[0]    
            for e in range(lrs - 1): 
                lay_idx = idx_1 + int(round(aqf_ratio_lst[e][0] * len(ibound_act_lay_idxs) / 100, 0))
                end_lay_idxs.append(lay_idx)
                act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[lay_idx - 1, 0, :].tolist()) if x == 1]
                idx_1 = lay_idx
            #end_lay_idxs.append(ibound_act_lay_idxs[-1])
                
            #   get the corresponding depth of these points to calculate intersections later on
            end_lay_idxs_m = []
            for g in range(len(end_lay_idxs)):
                end_lay_idxs_m.append([col_x, model_obj.top - (end_lay_idxs[g] * 10.)])
            end_lay_idxs_m.append(bot_point_lst[len(end_lay_idxs)])
            
            end_col_layer_pts.append(end_lay_idxs_m)
            lrs_cnt += 1
            offshore_aqf_thk_lst.append(aqf_ratio_lst)
        
        #   insert the starting points from the past shelf edges part into the point list that we will work with further
        end_col_layer_pts.insert(0, line_A_lst)
        end_col_layer_pts.append([line_B_lst[0]])
        
        pts_coord_start = end_col_layer_pts[0][:]
        ptc_coord_end = end_col_layer_pts[1][:]    
            
        #edge_count = 1        
        
        #for z in range(shelf_edge_idx, end_col_idx_lst[0]):    
        for z in range(shelf_edge_idx + 1, model_obj.ncol):#)range(shelf_edges_col_idx[0], model_obj.ncol):    
            ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, z].tolist()) if x == 1]       
            lrs_lst = []  

            if ibound_act_lay_idxs != []:

                edge_count = 1
                #   check the counter
                for c in range(len(end_col_idx_lst) - 1):
                    if z > end_col_idx_lst[c]:
                        edge_count += 1
                
                #   check if there are still cells from the upper glhymps 2 layer in the active cells offshore, if yes cut the list
                h = 0
                if last_glh_2_top_lay in ibound_act_lay_idxs:
                    while last_glh_2_top_lay > ibound_act_lay_idxs[h]:
                        #model_obj.hk_arr[ibound_act_lay_idxs[h], 0, z] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                        hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                        if hk_cell_val > hk_top_mean / 10.:
                            model_obj.hk_arr[ibound_act_lay_idxs[h], 0, z] = hk_cell_val
                        else:
                            model_obj.hk_arr[ibound_act_lay_idxs[h], 0, z] = hk_top_mean / 10.            
                        h += 1
                ibound_act_lay_idxs = ibound_act_lay_idxs[h:]                  
                
                #   go shelf edge by edge and fill in the columns based on the line crossing
                #pts_coord_start = cst_point_lst[:-edge_count] + line_A_lst[-edge_count:]
                #ptc_coord_end = line_A_lst[:-edge_count] + line_B_lst[-edge_count:]
                        
                if edge_count <= 1:
                    pts_coord_start = end_col_layer_pts[0][:]
                    ptc_coord_end = end_col_layer_pts[1][:]    
                else:
                    pts_coord_start = end_col_layer_pts[edge_count - 1][:-1]
                    ptc_coord_end = end_col_layer_pts[edge_count][:]            
                
                #   get the line of the column middle
                col_top_y = round(model_obj.top - (ibound_act_lay_idxs[0] + 1) * 10 , 1)
                col_bot_y = round(model_obj.top - (ibound_act_lay_idxs[-1] + 1) * 10, 1) 
                col_x =  round(model_obj.x_start - model_obj.cst_offset_plot + z / 10., 1)# + 0.05               
                       
                try:
                    #   get the line crossing for each of the layer lines
                    line_cross_lst, idx_cross_lst = [], []
                    for x in range(len(pts_coord_start)):
                        line_cross = line_intersection((pts_coord_start[x], ptc_coord_end[x]), ((col_x, col_top_y), (col_x, col_bot_y)))
                        #print(round(line_cross[0], 1), model_obj.top - round(line_cross[1], 0))
                        idx_cross_lst.append(int((model_obj.top - round(line_cross[1], 0)) / 10.))    
                    #   reverse the list to start from the top layer 
                    #idx_cross_lst = idx_cross_lst[::-1]        
                    
                    if idx_cross_lst[-1] > ibound_act_lay_idxs[-1] and len(idx_cross_lst) > 1:
                        idx_cross_lst[-1] = ibound_act_lay_idxs[-1]
                    
                    #   check that the indexes are actually active cells, if not then replace with the first active cell 
                    for g in range(len(idx_cross_lst)):
                        if idx_cross_lst[g] not in ibound_act_lay_idxs and idx_cross_lst[g] < ibound_act_lay_idxs[0]:
                            idx_cross_lst[g] = ibound_act_lay_idxs[0]

                    #   check that the indexes are not higher than the extend of the active cells
                    for s in range(len(idx_cross_lst)):
                        if idx_cross_lst[s] > ibound_act_lay_idxs[-1]:
                            idx_cross_lst[s] = ibound_act_lay_idxs[-1]

                    
                    """
                    #   fill in the HK array based on the list of indexes from previous step!
                    grp_lay_idx = [ibound_act_lay_idxs[:ibound_act_lay_idxs.index(idx_cross_lst[0])]]        
                    for y in range(1, len(idx_cross_lst)):
                        grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[y - 1]) : ibound_act_lay_idxs.index(idx_cross_lst[y])])
                    #   append the bottom layer till the end of the active model_obj domain
                    grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[-1]):])
                    """

                    #   fill in the HK array based on the list of indexes from previous step!
                    try:
                        grp_lay_idx = [ibound_act_lay_idxs[:ibound_act_lay_idxs.index(idx_cross_lst[0])]] 
                    #   it can happen that the cross-section happens in the inactive zone above the ocean bottom, in that case fill in empty list
                    except ValueError:
                        grp_lay_idx = [[]] 
                        
                    for y in range(1, len(idx_cross_lst)):
                        try:
                            grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[y - 1]) : ibound_act_lay_idxs.index(idx_cross_lst[y])])
                        #   same as above, the first index will become the first active cell in the ibound_act_lay_idx
                        except ValueError:
                            try:
                                grp_lay_idx.append(ibound_act_lay_idxs[: ibound_act_lay_idxs.index(idx_cross_lst[y])])
                            except ValueError:
                                grp_lay_idx.append(ibound_act_lay_idxs[: ibound_act_lay_idxs[0]])
                    #   append the bottom layer till the end of the active model_obj domain
                    grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[-1]):])
                    
                    #   make sure there are no duplicates or empty lists in the grp_lay_idx
                    #grp_lay_idx = [x for x in grp_lay_idx if x != []]
                    grp_lay_idx_sorted = []
                    for sublist in grp_lay_idx:
                        if sublist not in grp_lay_idx_sorted or sublist == []:
                            grp_lay_idx_sorted.append(sublist)
                    
                    #   check if the upper glhymps 2 layer (from soilgrids) is still located within the offshore model domain, if so
                    #   then assign the right property to those cells and cut the respective sublist to be filled in by the other layers
                    
                    for w in range(len(grp_lay_idx_sorted)):
                        #   + 1 because the upper layer is already filled in by GLHYMPS 2.0
                        if inland_aqf_lrs[w][1] == 'glhymps_2':
                            lrs_lst.append(['glhymps_2', grp_lay_idx_sorted[w]])
                            #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                            for t in range(len(grp_lay_idx_sorted[w])):
                                #model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, z] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                
                                hk_cell_val = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                                if hk_cell_val > hk_top_mean / 10.:
                                    model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, z] = hk_cell_val
                                else:
                                    model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, z] = hk_top_mean / 10.                                            
                        else:
                            lrs_lst.append(['glhymps_1', grp_lay_idx_sorted[w]])
                            for t in range(len(grp_lay_idx_sorted[w])):
                                model_obj.hk_arr[grp_lay_idx_sorted[w][t], 0, z] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)     

                    #   check if the first sublist of the lrs_lst is glhymps_1, if yes check if there are indexes lower or equal to the upp_lim
                    if lrs_lst[0][0] == 'glhymps_1':
                        #   check if the upp_lim index is in the list of indexes
                        try:
                            split_idx = lrs_lst[0][1].index(ibound_act_lay_idxs[upp_lim])                        
                        except ValueError:
                            pass
                    try:
                        if split_idx:
                            lrs_lst.insert(0, ['glhymps_2', lrs_lst[0][1][:split_idx]])
                            lrs_lst[1][1] = lrs_lst[1][1][split_idx:]
                    #   error - NameError: name 'split_idx' is not defined
                    except NameError:
                        pass
                 
                except Exception:
                    grp_lay_idx = grp_lay_idx#[ibound_act_lay_idxs]
                    
                """
                for w in range(len(grp_lay_idx)):
                    #   + 1 because the upper layer is already filled in by GLHYMPS 2.0
                    if inland_aqf_lrs[w][1] == 'glhymps_2':
                        lrs_lst.append(['glhymps_2', grp_lay_idx[w]])
                        #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                        for t in range(len(grp_lay_idx[w])):
                            model_obj.hk_arr[grp_lay_idx[w][t], 0, z] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                    else:
                        lrs_lst.append(['glhymps_1', grp_lay_idx[w]])
                        for t in range(len(grp_lay_idx[w])):
                            model_obj.hk_arr[grp_lay_idx[w][t], 0, z] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)           
                """
                
                upp_lim = 0
                edge_count = 1
                #   check the counter
                for c in range(len(end_col_idx_lst) - 1):
                    if z > end_col_idx_lst[c]:
                        edge_count += 1
                
            offshore_lrs_lst.append([z, lrs_lst])     

    #backup_arr = model_obj.hk_arr    
    #model_obj.hk_arr = backup_arr

    #   The next step is to adjust the interfaces between the coarse/fine grained sediment layer (GlHYMPS 1 vs GLHYMPS 2).
    #   This is done for only the cases where the sed_type is either 'large' or 'small', it will adjust the top of either 
    #   the GLHYMPS 2 or 1 layers respectively. 
    
    #   do this for each layer - or randomly selected number of layers
    #lay_cnt = len([item for item in offshore_lrs_lst[0][1] if item[0] == 'glhymps_1'])
    
    #off_lay_start_range = [-2.5, 2.5]  
    #off_lay_thk_ratio = 60.
    #lay_pres_y1 = 10. 
    clay_kh_minmax = [0.01, 0.0001]
    #clay_cap_shelf = False       #   insert the clay capping layer on top of the continental shelf/slope, the same reworking
    #clay_cap_slope = False       #   parameter lay_pres_y1 will apply to these layers as well
    #clay_cap_thk = 20.   
    
    #n_lrs = 9
    #p_fact = 0.8
    #n_clay_lrs = 3

    """
    p_val = p_fact
    n_lrs = len(lay_lst_1[y])
    n_clay_lrs = sed_cells
    """

    #   define the layers that will get the clay property based on the p_fact value
    def choose_rand_clay_lrs(p_val, n_lrs, n_clay_lrs):
        #   first create a cummulative sum of weights 
        intervals = [int(round(i, 2) * 100) for i in np.linspace(p_fact, 1 - p_fact, n_lrs).tolist()]
        intervals_cumsum = []
        for i in range(1, len(intervals) + 1):
            intervals_cumsum.append(sum(intervals[: i]))
        #   from the create list of sum of weights select the desired number of clay layers
        def get_rand_layer(int_lst):
            rand_number = random.randint(0, intervals_cumsum[-1])
            rand_lay = 0
            #   check which layer that number corresponds to 
            for g in range(len(intervals_cumsum)):
                if rand_number >= intervals_cumsum[g]:
                    rand_lay = rand_lay + 1
                else:
                    break
            #print(rand_number, rand_lay, '____', intervals_cumsum)
            return rand_lay
                
        clay_lrs_idx = []
        j = 0
        while j < n_clay_lrs:
            r_lay = get_rand_layer(intervals_cumsum)
            if j == 0:
                clay_lrs_idx.append(r_lay)
                j += 1
            else:
                if r_lay in clay_lrs_idx:
                    pass
                else:
                    clay_lrs_idx.append(r_lay)         
                    j += 1
        return clay_lrs_idx
    

    #   create an array with the dimensions of the ibound_arr, set all values to 0, the clay cells will be marked as 1 in this 
    #   array and then it will be used to assign all the clay values at once
    sed_arr = model_obj.ibound_arr * 0    
    
    #   select a random start of the clay layer
    #off_lay_start = round(random.uniform(off_lay_start_range[0], off_lay_start_range[1]), 1)
    
    #   get the column index where the offshore sed_type adjustment starts
    off_sed_start_idx = int((abs(model_obj.x_start - model_obj.cst_offset_plot) + off_lay_start) * (1000. / model_obj.del_col))

    #   combine the two lists
    tot_clay_lst = inland_lrs_lst + offshore_lrs_lst

    #   find the starting index in the offshore_lrs_lst based on the distance from coast calculated above
    off_lays_idxs = [item[0] for item in tot_clay_lst]
    
    #print off_lays_idxs, off_sed_start_idx
    try:
        off_lay_idx_start = off_lays_idxs.index(off_sed_start_idx)
    #   sometimes it throws an error that the off_sed_start_idx is not in the off_lays_idxs list, probably because of islands.
    except ValueError:
        off_lay_idx_start = off_sed_start_idx

    end_col_glhymps_1_start = end_col_idxs[1:][::2]
    end_col_glhymps_1_end = end_col_idxs[2:][::2][:-1]
    end_col_glhymps = []
    for d in range(len(end_col_glhymps_1_end)):
        end_col_glhymps.append([end_col_glhymps_1_start[d], end_col_glhymps_1_end[d]])
    end_col_glhymps.append([model_obj.ncol, model_obj.ncol])
        
    #   counter for the different end layers 
    end_lay_cnt = 0

    #   define the GLHYMPS layer to which the clay will be inserted, do it based on a sample and average of 10 values
    """
    if abs(glh_2_val) > abs(glh_1_val):
        clay_lay = 'glhymps_1'
    else:
        clay_lay = 'glhymps_2'
    """
    
    lst_glh_1_vals, lst_glh_2_vals = [], []
    for w in range(10):
        lst_glh_1_vals.append(round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4))
        lst_glh_2_vals.append(round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4))

    if np.mean(lst_glh_2_vals) > np.mean(lst_glh_1_vals):
        clay_lay = 'glhymps_1'
    else:
        clay_lay = 'glhymps_2'

    if sed_flux == 'low':
        
        #   loop through the columns where the sed_type adjustment will be implemented        
        for i in range(off_lay_idx_start, len(tot_clay_lst)):
            #   check for the 'glhymps_1 columns', select those from the offshore_lrs_lst
            col_idx = tot_clay_lst[i][0]
            
            try:
                #   select only the parts where there is 'glhymps_1' sediment type
                lay_idxs = [item for item in tot_clay_lst[i][1] if item[0] == clay_lay]
                #   get the total amount of cells that should have the clay/silt properties, based on the % filled in the off_lay_thk_ratio                    
                lay_lst_1 = [item[1] for item in tot_clay_lst[i][1] if item[0] == clay_lay]
                #lay_lst = [item for sublist in lay_lst_1 for item in sublist]
                #Maybe in future substract the clay capping cells from the ones that are then changed in the column..
                for y in range(len(lay_lst_1)):
                    if lay_lst_1[y] != []:
                        sed_cells = int(round(len(lay_lst_1[y]) * (off_lay_thk_ratio / 100.), 0))
                        #   define the layers that will get the clay property based on the p_fact value
                        clay_idxs = choose_rand_clay_lrs(p_fact, len(lay_lst_1[y]), sed_cells)     
                        for x in range(len(clay_idxs)):
                            #print(lay_lst_1[y][x], lay_lst_1[y][clay_idxs[x]])
                            #sed_arr[lay_lst_1[y][x], 0, col_idx] = 1  
                            sed_arr[lay_lst_1[y][clay_idxs[x]], 0, col_idx] = 1  
            except IndexError:
                pass   
        
        
    if sed_flux == 'medium' or sed_flux == 'high':             
        
        #   loop through the columns where the sed_type adjustment will be implemented        
        for i in range(off_lay_idx_start, len(tot_clay_lst)):
            #   check for the 'glhymps_1 columns', select those from the offshore_lrs_lst
            col_idx = tot_clay_lst[i][0]
            
            try:
                #print(i, col_idx)
                #   select only the parts where there is 'glhymps_1' sediment type
                lay_idxs = [item for item in tot_clay_lst[i][1] if item[0] == clay_lay]
                #   get the total amount of cells that should have the clay/silt properties, based on the % filled in the off_lay_thk_ratio                    
                lay_lst_1 = [item[1] for item in tot_clay_lst[i][1] if item[0] == clay_lay]
                #lay_lst = [item for sublist in lay_lst_1 for item in sublist]
                #Maybe in future substract the clay capping cells from the ones that are then changed in the column..
                for y in range(len(lay_lst_1)):
                    if lay_lst_1[y] != []:
                        sed_cells = int(round(len(lay_lst_1[y]) * (off_lay_thk_ratio / 100.), 0))
                        #   define the layers that will get the clay property based on the p_fact value
                        clay_idxs = choose_rand_clay_lrs(p_fact, len(lay_lst_1[y]), sed_cells)     
                        for x in range(len(clay_idxs)):
                            #print(lay_lst_1[y][x], lay_lst_1[y][clay_idxs[x]])
                            #sed_arr[lay_lst_1[y][x], 0, col_idx] = 1  
                            sed_arr[lay_lst_1[y][clay_idxs[x]], 0, col_idx] = 1   
            except IndexError:
                pass     
    

    #   if the shelf clay cap is turned on then add the clay cells in the upper part of the model domain based on the thickness indicated
    #   first, calculate the amount of layers that will get the clay HK value based on the thickness specified in clay_cap_thk
    cap_cells_shelf = int(round(clay_cap_shelf_thk / 10., 0))
    cap_cells_slope = int(round(clay_cap_slope_thk / 10., 0))
    
    if clay_cap_shelf:
        #   do this for cells between the coastline and the shelf edge
        for x in range(cst_idx_geo, shelf_edge_idx):
            #   select the upper part of the model domain, number of cells determined above
            ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, x].tolist()) if x == 1][:cap_cells_shelf]
            for layer in ibound_act_lay_idxs:
                if model_obj.top_elev[x] < 0.0:
                    sed_arr[layer, 0, x] = 1  
                
    #   same for the clay cap in the continental slope area
    if clay_cap_slope:
        #   do this for cells between the shelf edge and the end of the model domain
        for y in range(shelf_edge_idx, model_obj.ncol):
            #   select the upper part of the model domain, number of cells determined above
            ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, y].tolist()) if x == 1][:cap_cells_slope]
            for layer in ibound_act_lay_idxs:
                sed_arr[layer, 0, y] = 1      
    
    #   select randomly an amount of cells based on the lay_pres_y1 preservation index, and assign them back to 0. 
    #   this simulates the reworking of the internal layers with time and creates openings in these layers..
    sed_cells_idxs = np.where(sed_arr == 1)
    #   create an empty list and loop through the array created above to assign the location of each clay cell to the list
    sed_cells_idxs_lst = []
    for f in range(sed_cells_idxs[0].shape[0]):
        sed_cells_idxs_lst.append([sed_cells_idxs[0][f], sed_cells_idxs[1][f], sed_cells_idxs[2][f]])                       
    #   calculate the number of cells to be selected based on the lay_pres_y1 preservation index
    rework_cells_cnt = int(round(len(sed_cells_idxs_lst) * lay_pres_y1 / 100., 0))
    rework_lst = random.sample(sed_cells_idxs_lst, rework_cells_cnt)            
    #   now go through the list of the reworked cells and assign those back to 0 in the clay_arr
    for g in range(len(rework_lst)):
        sed_arr[rework_lst[g][0], 0, rework_lst[g][2]] = 0            

    #   now loop through the clay_arr, and for each occurrence of 1 change the value to a random value in the clay range
    for i in range(sed_arr.shape[0]):
        for j in range(sed_arr.shape[-1]):
            if sed_arr[i, 0, j] == 1:
                if not const_geo_hk_vals: 
                    model_obj.hk_arr[i, 0, j] = random.uniform(clay_kh_minmax[0], clay_kh_minmax[1])    
                else:
                    model_obj.hk_arr[i, 0, j] = clay_val
    

     #   specify the position of different parts of the overall figure/movie                        
    fig = plt.figure(figsize = (20,10))
    ax1 = plt.subplot2grid((3, 2), (0, 0), colspan = 2)  #   overall figure
    ax2 = plt.subplot2grid((3, 2), (1, 0))  #   zoomed in coastal zone
    ax3 = plt.subplot2grid((3, 2), (1, 1))  #   zoomed in shelf edge area
    ax4 = plt.subplot2grid((3, 2), (2, 0))  #   color bar area 
    ax5 = plt.subplot2grid((3, 2), (2, 1))  #   text/legend area
    
    ax1.set_position([0.05, 0.5, 0.90, 0.45]) # [left, bottom, width, height]
    ax2.set_position([0.05, 0.175, 0.425, 0.25])
    ax3.set_position([0.525, 0.175, 0.425, 0.25])        
    ax4.set_position([0.05, 0.05, 0.3, 0.05])
    ax5.set_position([0.4, 0.025, 0.55, 0.1])
    
    #   create a plot if desired
    plot_hk_arr = model_obj.hk_arr
    plot_hk_arr[np.abs(plot_hk_arr) <= 0.000001] = np.nan


    #   check that there are no NaN cells in the active model domain, if yes fill with the average value of the whole Hk array
    #for i in range(plot_hk_arr.shape[0]):
    
    """
    for i in range(model_obj.ibound_arr.shape[0]):  
        
        act_cell_idx_lst = [k for k, x in enumerate(model_obj.ibound_arr[i, 0, :].tolist()) if x == 1]      
        
        if act_cell_idx_lst != []:
            nan_cells = [l for l, y in enumerate(plot_hk_arr[i, 0, act_cell_idx_lst[0] : act_cell_idx_lst[-1]]) if math.isnan(l)]
            print(nan_cells)
            
            if nan_cells != []:
                for j in nan_cells:
                    plot_hk_arr[i, 0, j] = round(np.nanmean(model_obj.hk_arr), 6)
                    model_obj.hk_arr[i, 0, j] = round(np.nanmean(model_obj.hk_arr), 6)
    """             
    """
    hk_nan = np.isnan(plot_hk_arr)
    ib_mask = np.ma.masked_equal(model_obj.ibound_arr, 1)
    
    np.where(np.logical_and(hk_nan, ib_mask.mask))[0]
    """
    
    for i in range(model_obj.ibound_arr.shape[0]): 
        for j in range(plot_hk_arr.shape[-1]):
            try:
                if model_obj.ibound_arr[i, 0, j] == 1:
                    if np.isnan(plot_hk_arr[i, 0, j]):
                        #print(i, 0, j, round(np.nanmean(model_obj.hk_arr), 6))    
                        plot_hk_arr[i, 0, j] = round(np.nanmean(model_obj.hk_arr), 6)
            #   IndexError: index 144 is out of bounds for axis 0 with size 144
            except IndexError:
                pass
    
        
    #plot_hk_arr = model_obj.hk_arr
    #plot_hk_arr[np.abs(plot_hk_arr) <= 0.] = np.nan

    #print np.nanmin(plot_hk_arr), np.nanmax(plot_hk_arr)

    #im = ax1.imshow(model_obj.ibound_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = 'Reds',
    #               extent = (model_obj.x_start - model_obj.cst_offset_plot, model_obj.x_end - model_obj.cst_offset_plot, model_obj.zbot, math.ceil((model_obj.top / 100.0) * 100.0)), vmin = 0., vmax = 1.)

    #   define the colormap
    cmap = cmx.viridis
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # define the bins and normalize
    unique_n = [0.0, 0.001, 0.1, 1.0, 5.0, 10., 20., 30.]# + list(np.arange(10.0, round(max(unique_nan) / 10) * 10, 10.))
    norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)            
    
    #x_lines_topelev = np.linspace(model_obj.x_start + (model_obj.del_col / (2 * 1000.)) - model_obj.cst_offset_plot, model_obj.x_end - model_obj.cst_offset_plot - (model_obj.del_col / (2 * 1000.)), model_obj.ncol)
    #y_lines = np.linspace(model_obj.top, model_obj.zbot, model_obj.nlay + 1)    
    
    x_end_new = round(model_obj.x_start - model_obj.cst_offset_plot, 1) + model_obj.ncol / 10.
    x_start_new = round(model_obj.x_start - model_obj.cst_offset_plot, 1)
    x_lines_topelev = np.linspace(x_start_new,  x_end_new, model_obj.ncol + 1)
    y_lines = np.linspace(model_obj.top, model_obj.zbot, model_obj.ibound_arr.shape[0] + 1)    
    
    #   reverse the list of layer ends at the coastline
    if sed_flux == 'medium' or sed_flux == 'high': 
        cst_lay_pts_plot = cst_lay_pts[::-1]
    else:
        cst_lay_pts_plot = []
    
    def plot_kh_arr(axis, x_start, x_end, y_start, y_end, title):
        im = axis.imshow(plot_hk_arr[:, 0, :], aspect = 'auto', interpolation = None, cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, model_obj.zbot, math.ceil((model_obj.top / 100.0) * 100.0)), vmin = 0., vmax = np.nanmax(plot_hk_arr))
        
        #lineA, = axis.plot(x_lines_topelev, model_obj.top_elev, c = 'black', linewidth = 3, label = 'GEBCO elevation')
        
        axis.set_xlim([x_start, x_end])    
        axis.set_ylim([y_end, y_start])
        
        axis.axhline(y = 0., linewidth = 1, color = 'k')
        axis.axvline(x = 0., linewidth = 1, color = 'k')
        axis.axvline(x = cont_shelf_edge[0], linewidth = 1, color = 'k')    
            
        """
        #   draw the shelf edge points
        for i in range(len(shelf_edges_col_idx)):
            x_coord = model_obj.x_start - model_obj.cst_offset_plot + shelf_edges_col_idx[i] / 10.# + 0.05
            y_coord = model_obj.top - shelf_edges_lay_idx[i] * 10.# + 5.
            #y_coord_coast = model_obj.top - (cst_lay_pts_plot[i] + 1) * 10.    
            #axis.plot(x_coord, y_coord, 'ro', markersize = 4)
            #axis.plot(0, y_coord_coast, 'ro', markersize = 4)
            
            axis.plot(line_A_lst[i][0], line_A_lst[i][1], 'ro', markersize = 4)
            axis.plot(cst_point_lst[i][0], cst_point_lst[i][1], 'ro', markersize = 4)       
            axis.plot([cst_point_lst[i][0], line_A_lst[i][0]], [cst_point_lst[i][1], line_A_lst[i][1]], color = 'red')
            axis.plot(bot_point_lst[i][0], bot_point_lst[i][1], 'bo', markersize = 4)    
            axis.plot([line_A_lst[i][0], bot_point_lst[i][0]], [line_A_lst[i][1], bot_point_lst[i][1]], color = 'blue')        
               
            #axis.plot([0., x_coord], [y_coord_coast, y_coord], color = 'red')
            #axis.plot(bot_point_lst[i][0], bot_point_lst[i][1], 'bo', markersize = 4)    
            #axis.plot([x_coord, bot_point_lst[i][0]], [y_coord, bot_point_lst[i][1]], color = 'blue')
        
        for x_line in end_col_idxs:
            axis.axvline(x = x_lines_topelev.tolist()[x_line], linestyle='--', linewidth = .5, color = 'k')            
        for y_line in end_lay_idxs:
            axis.axhline(y = y_lines.tolist()[y_line], linestyle='--', linewidth = .5, color = 'k')  
        """
        
        axis.set_title(title, fontsize = 22)
        axis.set_xlabel('distance from coast (km)', fontsize = 12)
        axis.set_ylabel('elevation (m asl.)', fontsize = 12)
        
        return im
    
    if clay_cap_shelf is True:
        clay_cap_shelf_str = 'yes'
        clay_cap_shelf_thk_str = str(clay_cap_shelf_thk)
    else:
        clay_cap_shelf_str = 'no'
        clay_cap_shelf_thk_str = 'None'
    if clay_cap_slope is True:
        clay_cap_slope_str = 'yes'
        clay_cap_slope_thk_str = str(clay_cap_slope_thk)
    else:
        clay_cap_slope_str = 'no'      
        clay_cap_slope_thk_str = 'None'
    
    txt = ax5.text(0.1, 0.25, 'Qs = %s' % (sed_flux) + '                    ' + 'pres_fact = %s' % str(round(lay_pres_y1, 1)) + '              '\
                   + 'stacking factor P = %s' % str(round(p_fact, 2)) + '\n'\
                   + 'clay cap slope/shelf (thk in m) = %s, %s (%s, %s)' %\
                   (clay_cap_shelf_str, clay_cap_slope_str, clay_cap_shelf_thk_str, clay_cap_slope_thk_str), fontsize = 14, fontweight = 'bold', clip_on = False)
    ax5.axis('off')
    
    y_zoom_max_cst = math.ceil(max(model_obj.top_elev[int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) : int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) + 100]))
    y_zoom_min_cst = math.floor(min(model_obj.bot_elev[int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) : int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) + 100]))
    try:
        y_zoom_max_shelf = math.ceil(max(model_obj.top_elev[shelf_edge_idx - 100 : shelf_edge_idx + 100]))
        y_zoom_min_shelf = math.floor(min(model_obj.bot_elev[shelf_edge_idx - 100 : shelf_edge_idx + 100]))
    #   ValueError: max() arg is an empty sequence
    #   can happen if the inland part of the model domain is too short
    except ValueError:
        y_zoom_max_shelf = math.ceil(max(model_obj.top_elev[0 : shelf_edge_idx + 100]))
        y_zoom_min_shelf = math.floor(min(model_obj.bot_elev[0 : shelf_edge_idx + 100]))        
        
    main_title = 'Horizontal hydraulic conductivity (Hk) distribution in the model domain of COSCAT %s, realization %s' %(str(model_obj.coscat_id), str(rand_seed_in))
    im1 = plot_kh_arr(ax1, x_start_new,  x_end_new, max(model_obj.top_elev), min(model_obj.bot_elev), main_title)
    im2 = plot_kh_arr(ax2, -5.0, 5.0, y_zoom_max_cst, y_zoom_min_cst, 'Hk, zoom in coastal zone')
    im3 = plot_kh_arr(ax3, cont_shelf_edge[0] - 10.0, cont_shelf_edge[0] + 10.0, max(0.0, y_zoom_max_shelf), y_zoom_min_shelf, 'Hk, zoom in shelf edge zone')
    
    #   plot the colorbar
    cbar = plt.colorbar(im1, cax = ax4, cmap = cmap, norm = norm, spacing = 'uniform', ticks = unique_n, boundaries = unique_n, orientation='horizontal')
    ax4.set_xticklabels([str(e) for e in unique_n])   
    ax4.set_title("Horizontal hydraulic conductivity (m/d)", fontsize = 14)#, y  = 1.025)
    ax4.tick_params(labelsize = 12)    
    
    rcParams['font.family'] = 'Garamond'
    rcParams['axes.facecolor'] = 'white'
    rcParams['savefig.facecolor'] = 'white'  
    
    
    #plt.show()
    
    #figname = '_coscat_%s' % (str(model_obj.coscat_id)) + '__%s' % (model_obj.cst_type) + '__GeoRandSc_%s' % str(rand_seed_in)
    
    plt.savefig(os.path.join(save_dir, figname))
    plt.savefig(os.path.join(summary_save_dir, figname))
    plt.close(fig)
    #del fig   
    
    plot_hk_arr = np.round(plot_hk_arr, 8)
    
    return plot_hk_arr
    
    """
    
    
    
    
    
    
    
    
    
    for g in xrange(lay_cnt):

        #   select a random start of the clay layer
        off_lay_start = round(random.uniform(off_lay_start_range[0], off_lay_start_range[1]), 1)
        
        #   get the column index where the offshore sed_type adjustment starts
        off_sed_start_idx = int((abs(model_obj.x_start - model_obj.cst_offset_plot) + off_lay_start) * (1000. / model_obj.del_col))
    
        #   find the starting index in the offshore_lrs_lst based on the distance from coast calculated above
        off_lays_idxs = [item[0] for item in offshore_lrs_lst]
        
        #print off_lays_idxs, off_sed_start_idx
        off_lay_idx_start = off_lays_idxs.index(off_sed_start_idx)
    
        end_col_glhymps_1_start = end_col_idxs[1:][::2]
        end_col_glhymps_1_end = end_col_idxs[2:][::2][:-1]
        end_col_glhymps = []
        for d in xrange(len(end_col_glhymps_1_end)):
            end_col_glhymps.append([end_col_glhymps_1_start[d], end_col_glhymps_1_end[d]])
        end_col_glhymps.append([model_obj.ncol, model_obj.ncol])
            
        #   counter for the different end layers 
        end_lay_cnt = 0
    
        if sed_flux == 'low':
    
            #   loop through the columns where the sed_type adjustment will be implemented        
            for i in range(off_lay_idx_start, len(offshore_lrs_lst)):
                #   check for the 'glhymps_1 columns', select those from the offshore_lrs_lst
                col_idx = offshore_lrs_lst[i][0]
                
                try:
                    #   select only the parts where there is 'glhymps_1' sediment type
                    lay_idxs = [item for item in offshore_lrs_lst[i][1] if item[0] == 'glhymps_1'][g]
                    #   get the total amount of cells that should have the clay/silt properties, based on the % filled in the off_lay_thk_ratio                    
                    lay_lst_1 = [item[1] for item in offshore_lrs_lst[i][1] if item[0] == 'glhymps_1']
                    lay_lst = [item for sublist in lay_lst_1 for item in sublist]
                    # Maybe in future substract the clay capping cells from the ones that are then changed in the column..
                    sed_cells = int(round(len(lay_lst) * (off_lay_thk_ratio / 100.), 0))
                    
                    #   find the indexes of the start and end of the layer at continental slope
                    start_lay_idx, end_lay_idx = end_col_glhymps[end_lay_cnt][0], end_col_glhymps[end_lay_cnt][1]
                    
                    #   in case the column is before the continental shelf edge then do not assign the clay to the upper layer 
                    #if col_idx <= shelf_edge_idx:
                        #   get the amount of layers to be adjusted, create a counter of total clay cells assigned in the column
                    
                    if col_idx <= start_lay_idx:            
                        tot_lrs = len(lay_idxs[1]) - 1
                        sed_cell_cnt = 0
                        
                        try:
                            sed_cells_per_layer = int(round(sed_cells / float(tot_lrs)))
                        #   this happens when there is no more glhymps_1 layer where the clay layer could be inserted
                        except ZeroDivisionError:
                            break
                        
                        all_sed_lrs_lst = []
                        for j in range(1, len(lay_idxs)):
                            #   check if it is the last layer (thats where any odd layer cell is added)
                            if j != len(lay_idxs) - 1:
                                try:
                                    #   select the cells from the upper glhymps_2 layer and from the current glhymps_1 layer
                                    first_idx = lay_idxs[1][0]
                                    #   based on the 
                                    lays = range(first_idx - (sed_cells_per_layer - 1) / 2, first_idx + (sed_cells_per_layer - 1) / 2 + 1)
                                    sed_cell_cnt += len(lays)
                                    all_sed_lrs_lst.append(lays)
                                except IndexError:
                                    print j
                                    pass
                            else:
                                sed_cells_per_layer = sed_cells - sed_cell_cnt
                                try:
                                    first_idx = lay_idxs[1][0]
                                    lays = range(first_idx - (sed_cells_per_layer - 1) / 2, first_idx + (sed_cells_per_layer - 1) / 2 + 1)
                                    #   add additional layer if necessary
                                    if sed_cells_per_layer - len(lays) != 0:
                                        lays.append(lays[-1] + 1)
                                    all_sed_lrs_lst.append(lays)  
                                except IndexError:
                                    print j
                                    pass
                    
                    else:
                        tot_lrs = len(lay_idxs)
                        sed_cell_cnt = 0  
                        
                        try:
                            sed_cells_per_layer = int(round(sed_cells / float(tot_lrs)))
                        #   this happens when there is no more glhymps_1 layer where the clay layer could be inserted
                        except ZeroDivisionError:
                            break
                        
                        all_sed_lrs_lst = []
                        for j in range(0, len(lay_idxs)):
                            #   check if it is the last layer (thats where any odd layer cell is added)
                            if j != len(lay_idxs) - 1:
                                try:
                                    #   select the cells from the upper glhymps_2 layer and from the current glhymps_1 layer
                                    first_idx = lay_idxs[1][0]
                                    #   based on the 
                                    lays = range(first_idx - (sed_cells_per_layer - 1) / 2, first_idx + (sed_cells_per_layer - 1) / 2 + 1)
                                    sed_cell_cnt += len(lays)
                                    all_sed_lrs_lst.append(lays)
                                except IndexError:
                                    print j
                                    pass
                            else:
                                try:
                                    sed_cells_per_layer = sed_cells - sed_cell_cnt
                                    first_idx = lay_idxs[1][0]
                                    lays = range(first_idx - (sed_cells_per_layer - 1) / 2, first_idx + (sed_cells_per_layer - 1) / 2 + 1)
                                    #   add additional layer if necessary
                                    if sed_cells_per_layer - len(lays) != 0:
                                        lays.append(lays[-1] + 1)
                                    all_sed_lrs_lst.append(lays)
                                except IndexError:
                                    print j
                                    pass
                                
                    if col_idx >= end_lay_idx:            
                        end_lay_cnt += 1
                    else:
                        pass
        
                    #   flatten the list and change the indexes in the clay_arr
                    lay_lst_to_arr = [item for sublist in all_sed_lrs_lst for item in sublist]
                        
                    for layer in lay_lst_to_arr:
                        #   first check if that cell is actually active in IBOUND
                        if model_obj.ibound_arr[layer, 0, col_idx] == 1:
                            sed_arr[layer, 0, col_idx] = 1   
                
                except IndexError:
                    pass

    #   select randomly an amount of cells based on the lay_pres_y1 preservation index, and assign them back to 0. 
    #   this simulates the reworking of the internal layers with time and creates openings in these layers..
    sed_cells_idxs = np.where(sed_arr == 1)
    #   create an empty list and loop through the array created above to assign the location of each clay cell to the list
    sed_cells_idxs_lst = []
    for f in xrange(sed_cells_idxs[0].shape[0]):
        sed_cells_idxs_lst.append([sed_cells_idxs[0][f], sed_cells_idxs[1][f], sed_cells_idxs[2][f]])                       
    #   calculate the number of cells to be selected based on the lay_pres_y1 preservation index
    rework_cells_cnt = int(round(len(sed_cells_idxs_lst) * lay_pres_y1 / 100., 0))
    rework_lst = random.sample(sed_cells_idxs_lst, rework_cells_cnt)            
    #   now go through the list of the reworked cells and assign those back to 0 in the clay_arr
    for g in xrange(len(rework_lst)):
        sed_arr[rework_lst[g][0], 0, rework_lst[g][2]] = 0            
    
    #   now loop through the clay_arr, and for each occurrence of 1 change the value to a random value in the clay range
    for i in xrange(sed_arr.shape[0]):
        for j in xrange(sed_arr.shape[-1]):
            if sed_arr[i, 0, j] == 1:
                model_obj.hk_arr[i, 0, j] = random.uniform(clay_kh_minmax[0], clay_kh_minmax[1])

















        
    if sed_flux == 'medium' or sed_flux == 'high': 
              
        #   first calculate the average slopes of the cont. shelf and cont. slope, as the elevation difference between the shelf edge and the
        #   elevation at the coastline
        #avg_sl_shelf = round(100 * (model_obj.top_elev[len(soilgrids_thk_lst) + 1] - model_obj.top_elev[shelf_edge_idx]) / ((shelf_edge_idx - len(soilgrids_thk_lst) + 1) * 100.), 2)
        avg_sl_slope = round(100 * (model_obj.top_elev[shelf_edge_idx] - model_obj.top_elev[-1]) / ((model_obj.ncol - shelf_edge_idx) * 100.), 2)
        
        math.degrees(math.atan((model_obj.top_elev[shelf_edge_idx] - model_obj.top_elev[-1]) / ((model_obj.ncol - shelf_edge_idx) * 100.)))
        
        #   then get the top points of each layer (not the top one..), different for each scenario
        cst_lay_pts = []
        
        #col_top_lay_n = int(soilgrids_thk_mean / model_obj.del_lay)
        col_top_lay_n = int(round(soilgrids_thk_lst[r] / model_obj.del_lay))     
        cst_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, cst_offset_idx].tolist()) if x == 1]
        
        shelf_edge_idx = model_obj.top_elev.index(cont_shelf_edge[1]) 
        shelf_edge_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edge_idx].tolist()) if x == 1]   
        
        #   divide the area between the shelf edge and model_obj bottom into the tot_layers identical amount of cells based on the ratio 
        tot_layers = len(inland_aqf_lrs)
        tot_cells = len(cst_ibound_act_lay_idxs[col_top_lay_n:])
        cst_lay_pts.append(cst_ibound_act_lay_idxs[col_top_lay_n])
        cell_cnt_lst = []
        geo_lay_thk_start = col_top_lay_n
        tot_cells = len(ibound_act_lay_idxs) - col_top_lay_n
        for q in range(tot_layers):
            geo_lay_thk = (inland_aqf_lrs[q][0] * tot_cells) / 100
            #   check if the calculated thickness is 0 then just go to the next layer
            if geo_lay_thk == 0:
                pass
            else:
                #   len(inland_aqf_lrs) - 1 because of python counting..
                if q == len(inland_aqf_lrs) - 1:
                    pass
                else:
                    cst_lay_pts.append(cst_ibound_act_lay_idxs[geo_lay_thk_start + geo_lay_thk])
                    geo_lay_thk_start += geo_lay_thk          
        
        #   get the shelf edge points for each layer, based on the angle between the perpendicular line at the newest shelf edge
        #   and the dip line (see pictures)
        
        #   calculate the angle between the aquifer bottom and the perpendicular line at the shelf edge
        aqf_bot_y_diff = round(model_obj.bot_elev[len(soilgrids_thk_lst) + 1] - model_obj.bot_elev[shelf_edge_idx], 2)
        aqf_bot_x_diff = round((shelf_edge_idx - len(soilgrids_thk_lst) + 1) * 100., 2)
        
        shelf_angle = 85.
        shelf_edge_bot_angle = 90. - round(math.degrees(math.atan(aqf_bot_y_diff / aqf_bot_x_diff)), 2)
        
        cst_angle = round(180. - shelf_angle - shelf_edge_bot_angle, 2)   # the missing 3rd angle
        
        len_shelf = len(shelf_edge_ibound_act_lay_idxs) * 10.
        sin_ratio = len_shelf / math.sin(math.radians(cst_angle))
        len_bot = round(sin_ratio * math.sin(math.radians(shelf_angle)), 2)
        len_top = round(sin_ratio * math.sin(math.radians(shelf_edge_bot_angle)), 2)
                   
        print len_bot, len_top
    
        #   split the len_top into tot_layers parts equally and get the [layer, column] tuple for each point
        #       first calculate the perpendicular distance of the edge from the shelf edge
        dist_from_edge = round(math.sin(math.radians(shelf_angle)) * len_top, 2)
        dist_in_cols = int(round(dist_from_edge / 100.))
 
        shelf_edges_col_idx = []
        pct_start = inland_aqf_lrs[0][0]
        for i in xrange(len(inland_aqf_lrs) - 1):
            dist_col = (dist_in_cols * pct_start) / 100.
            shelf_edges_col_idx.append(shelf_edge_idx - int(round(dist_col, 0)))
            pct_start += inland_aqf_lrs[i + 1][0]
        
        lay_top = shelf_edge_ibound_act_lay_idxs[col_top_lay_n] 
        lay_bot = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edge_idx - dist_in_cols].tolist()) if x == 1][-1]   
        
        shelf_edges_lay_idx = []
        pct_start = inland_aqf_lrs[0][0]
        for i in xrange(len(inland_aqf_lrs) - 1):
            dist_lay = ((lay_bot - lay_top) * pct_start) / 100.
            shelf_edges_lay_idx.append(lay_top + int(round(dist_lay, 0)))
            pct_start += inland_aqf_lrs[i + 1][0]        
        
        #   find the end of the layer based on the slope of the continental slope
        line_A_lst, line_B_lst, cst_point_lst = [], [], []
        line_thk_end = ((0.0, model_obj.bot_elev[cst_offset_idx]), (model_obj.x_end - model_obj.cst_offset_plot, model_obj.zbot))   
        
        for i in xrange(len(shelf_edges_col_idx)):
            x_coord = model_obj.x_start - model_obj.cst_offset_plot + shelf_edges_col_idx[i] / 10.# + 0.05
            y_coord = model_obj.top - shelf_edges_lay_idx[i] * 10.# + 5.
     
             #  calculate the y_coordinate at the coastline
            y_coast = model_obj.top - (cst_lay_pts[i + 1]) * 10 
          
            #   calculate the y_coordinate of the line at the end of the model_obj domain
            y_coord_end = round(avg_sl_slope * (model_obj.x_end - model_obj.cst_offset_plot - x_coord) * 10., 1)
    
            line_A_lst.append((round(x_coord, 1), round(y_coord, 1)))
            line_B_lst.append((model_obj.x_end - model_obj.cst_offset_plot, round(y_coord - y_coord_end, 1)))
            cst_point_lst.append((0.0, round(y_coast)))
    
        bot_point_lst = []
        for j in xrange(len(line_A_lst)):
            bot_pt = line_intersection((line_A_lst[j], line_B_lst[j]), line_thk_end)
            bot_point_lst.append([round(bot_pt[0], 1), round(bot_pt[1], 1)])


        for a in range(len(soilgrids_thk_lst), shelf_edges_col_idx[-1]):
              
            tot_layers = len(inland_aqf_lrs)
            shelf_edge_0_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edges_col_idx[-1]].tolist()) if x == 1]   

            #   find the limits for each of the geological layers 
            idx_1 = shelf_edge_0_ibound_act_lay_idxs[col_top_lay_n]                
            end_lay_idxs =  [idx_1]
            act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[idx_1, 0, :].tolist()) if x == 1]
            end_col_idxs = [act_cell_lst[-1]]                
            for b in xrange(tot_layers - 1):                    
                lay_idx = idx_1 + (inland_aqf_lrs[b][0] * (model_obj.nlay - shelf_edge_0_ibound_act_lay_idxs[col_top_lay_n])) / 100              
                end_lay_idxs.append(lay_idx)
                act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[lay_idx - 1, 0, :].tolist()) if x == 1]
                end_col_idxs.append(act_cell_lst[-1])
                idx_1 = lay_idx
            #   append the last active layer and column as the end of the last layer
            end_lay_idxs.append(model_obj.nlay)
            end_col_idxs.append(int(model_obj.ncol) - 1)
                
            #   check if there are any ibound_act_lay_idxs indexes below the upper aquifer part
            #   fill in the part till the shelf edge with the constant thickness, after that point fill in straight line
            #   the upper part of the offshore aquifer domain
            ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, a].tolist()) if x == 1]
            for s in xrange(col_top_lay_n):
                #   try except - catch the IndexError in case the SOILGRIDS thickness is larger than the active thickness
                #   of the IBOUND array in the current column
                try:
                    model_obj.hk_arr[ibound_act_lay_idxs[s], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                except IndexError:
                    pass    

            #   bottom part of the aquifer domain
            #   loop through all the groups and assign the right GLHYMPS values       
            cell_cnt_lst = []
            geo_lay_thk_start = col_top_lay_n
            tot_cells = len(ibound_act_lay_idxs) - col_top_lay_n
            for q in range(tot_layers):
                geo_lay_thk = (inland_aqf_lrs[q][0] * tot_cells) / 100
                #   check if the calculated thickness is 0 then just go to the next layer
                if geo_lay_thk == 0:
                    pass
                else:
                    #   len(inland_aqf_lrs) - 1 because of python counting..
                    if q == len(inland_aqf_lrs) - 1:
                        #   if it is the last layer make sure that all the remaining cells in the column are select and no empty
                        #   cells are created in the process
                        cell_cnt_lst.append(ibound_act_lay_idxs[geo_lay_thk_start :])
                    else:
                        cell_cnt_lst.append(ibound_act_lay_idxs[geo_lay_thk_start : geo_lay_thk_start + geo_lay_thk])
                        geo_lay_thk_start += geo_lay_thk                        
                
                
            #   loop through all the groups and assign the right GLHYMPS values            
            for w in range(len(cell_cnt_lst)):
                #   + 1 because the upper layer is already filled in by GLHYMPS 2.0
                if inland_aqf_lrs[w][1] == 'glhymps_2':
                    #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                    for t in xrange(len(cell_cnt_lst[w])):
                        model_obj.hk_arr[cell_cnt_lst[w][t], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                else:
                    for t in xrange(len(cell_cnt_lst[w])):
                        model_obj.hk_arr[cell_cnt_lst[w][t], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)    
                
        edge_count = 1        
        for b in range(shelf_edges_col_idx[-1], shelf_edges_col_idx[0]):
              
            tot_layers = len(inland_aqf_lrs)
            shelf_edge_last_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, shelf_edges_col_idx[0]].tolist()) if x == 1]   

            #   find the limits for each of the geological layers 
            idx_1 = shelf_edge_0_ibound_act_lay_idxs[col_top_lay_n]                
            end_lay_idxs =  [idx_1]
            act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[idx_1, 0, :].tolist()) if x == 1]
            end_col_idxs = [act_cell_lst[-1]]                
            for e in xrange(tot_layers - 1):                    
                lay_idx = idx_1 + (inland_aqf_lrs[e][0] * (model_obj.nlay - shelf_edge_last_ibound_act_lay_idxs[col_top_lay_n])) / 100              
                end_lay_idxs.append(lay_idx)
                act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[lay_idx - 1, 0, :].tolist()) if x == 1]
                end_col_idxs.append(act_cell_lst[-1])
                idx_1 = lay_idx
            #   append the last active layer and column as the end of the last layer
            end_lay_idxs.append(model_obj.nlay)
            end_col_idxs.append(int(model_obj.ncol) - 1)
                
            #   get the active layers in the column
            col_top_lay_n = int(soilgrids_thk_lst[r] / model_obj.del_lay)
            ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, b].tolist()) if x == 1]
    
            #   go shelf edge by edge and fill in the columns based on the line crossing
            pts_coord_start = line_A_lst[:edge_count] + cst_point_lst[edge_count:]
            ptc_coord_end = line_B_lst[:edge_count] + line_A_lst[edge_count:]
                    
            #   get the line of the column middle
            col_top_y = round(model_obj.top - (ibound_act_lay_idxs[col_top_lay_n] + 1) * 10 , 1)
            col_bot_y = round(model_obj.top - (ibound_act_lay_idxs[-1] + 1) * 10, 1) 
            col_x =  round(model_obj.x_start - model_obj.cst_offset_plot + b / 10., 1)# + 0.05               
                    
            #   get the line crossing for each of the layer lines
            line_cross_lst, idx_cross_lst = [], []
            for x in xrange(len(pts_coord_start)):
                line_cross = line_intersection((pts_coord_start[x], ptc_coord_end[x]), ((col_x, col_top_y), (col_x, col_bot_y)))
                #print round(line_cross[0], 1), model_obj.top - round(line_cross[1], 0)
                idx_cross_lst.append(int((model_obj.top - round(line_cross[1], 0)) / 10.))    
            #   reverse the list to start from the top layer 
            #idx_cross_lst = idx_cross_lst[::-1]    

            #   fill in the HK array based on the list of indexes from previous step!
            grp_lay_idx = [ibound_act_lay_idxs[:ibound_act_lay_idxs.index(idx_cross_lst[0])]]        
            for y in xrange(1, len(idx_cross_lst)):
                grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[y - 1]) : ibound_act_lay_idxs.index(idx_cross_lst[y])])
            #   append the bottom layer till the end of the active model_obj domain
            grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[-1]):])
             
            for w in range(len(grp_lay_idx)):
                #   + 1 because the upper layer is already filled in by GLHYMPS 2.0
                if inland_aqf_lrs[w][1] == 'glhymps_2':
                    #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                    for t in xrange(len(grp_lay_idx[w])):
                        model_obj.hk_arr[grp_lay_idx[w][t], 0, b] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                else:
                    for t in xrange(len(grp_lay_idx[w])):
                        model_obj.hk_arr[grp_lay_idx[w][t], 0, b] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)           
                    
            #   fill in the upper part with the upper layer
            for s in xrange(col_top_lay_n):
                try:
                    model_obj.hk_arr[ibound_act_lay_idxs[s], 0, b] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                except IndexError:
                    pass   
    
            #   check the counter
            for c in range(len(shelf_edges_col_idx) - 1):
                if b > shelf_edges_col_idx[c]:
                    edge_count =+ 1                         

        for z in range(shelf_edges_col_idx[0], model_obj.ncol):
            #   get the active layers in the column
            col_top_lay_n = int(soilgrids_thk_lst[r] / model_obj.del_lay)
            ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, z].tolist()) if x == 1]        
    
            #   sometimes there are no active layers at the far edge of the model_obj domain
            if ibound_act_lay_idxs == []:
                pass
            else:
                #   go shelf edge by edge and fill in the columns based on the line crossing
                pts_coord_start = line_A_lst
                ptc_coord_end = bot_point_lst
                        
                #   get the line of the column middle
                try:
                    col_top_y = round(model_obj.top - (ibound_act_lay_idxs[col_top_lay_n] + 1) * 10 , 1)
                except IndexError:
                    col_top_y = round(model_obj.top - (ibound_act_lay_idxs[0] + 1) * 10 , 1)
                    
                col_bot_y = round(model_obj.top - (ibound_act_lay_idxs[-1] + 1) * 10, 1) 
                col_x =  round(model_obj.x_start - model_obj.cst_offset_plot + z / 10., 1)# + 0.05    
        
                #   get the line crossing for each of the layer lines
                line_cross_lst, idx_cross_lst = [], []
                for x in xrange(len(pts_coord_start)):
                    try:
                        line_cross = line_intersection((pts_coord_start[x], ptc_coord_end[x]), ((col_x, col_top_y), (col_x, col_bot_y)))
                        print round(line_cross[0], 1), model_obj.top - round(line_cross[1], 0)
                        idx_cross_lst.append(int((model_obj.top - round(line_cross[1], 0)) / 10.) + 1)    
                    except Exception:
                        #   set the idx_cross_lst to nlay so it always fits the next conditions
                        idx_cross_lst = [model_obj.nlay]
                        break
                #   reverse the list to start from the top layer 
                #idx_cross_lst = idx_cross_lst[::-1] 
        
                #   if the last element of the intersection layer indexes is bigger than the last index of the active column it means
                #   that it is the end of the model_obj domain and there is only one layer present
                if idx_cross_lst[0] >= ibound_act_lay_idxs[-1]:
                    grp_lay_idx = [ibound_act_lay_idxs]
        
                else:
                    #   fill in the HK array based on the list of indexes from previous step!
                    grp_lay_idx = [ibound_act_lay_idxs[:ibound_act_lay_idxs.index(idx_cross_lst[0])]]        
                    for y in xrange(1, len(idx_cross_lst)):
                        try:
                            grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(idx_cross_lst[y - 1]) : ibound_act_lay_idxs.index(idx_cross_lst[y])])
                        #   when index is out of the active layer indeces in the model_obj column
                        except ValueError:
                            pass
                    #   append the bottom layer till the end of the active model_obj domain
                    grp_lay_idx.append(ibound_act_lay_idxs[ibound_act_lay_idxs.index(grp_lay_idx[-1][-1]):])
                        
                for w in range(len(grp_lay_idx)):
                    #   + 1 because the upper layer is already filled in by GLHYMPS 2.0
                    if inland_aqf_lrs[w][1] == 'glhymps_2':
                        #   for every even layer fill in with GLHMYPS 2.0 (even is 0, 2, 4..)
                        for t in xrange(len(grp_lay_idx[w])):
                            model_obj.hk_arr[grp_lay_idx[w][t], 0, z] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                    else:
                        for t in xrange(len(grp_lay_idx[w])):
                            model_obj.hk_arr[grp_lay_idx[w][t], 0, z] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)  
    
                for s in xrange(col_top_lay_n):
                    #   try except - catch the IndexError in case the SOILGRIDS thickness is larger than the active thickness
                    #   of the IBOUND array in the current column
                    try:
                        if ibound_act_lay_idxs[s] <= shelf_edge_ibound_act_lay_idxs[col_top_lay_n]:    
                            print("upper")
                            model_obj.hk_arr[ibound_act_lay_idxs[s], 0, z] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                        else:          
                            pass
                            #model_obj.hk_arr[ibound_act_lay_idxs[s], 0, a] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                    except IndexError:
                        pass
                        
    

    # Concentration and Flow Plot
    #fig = plt.figure(figsize = (20, 20))
    #ax1 = plt.subplot2grid((2, 1), (0, 0))  #   overall concentration profiles
    #ax2 = plt.subplot2grid((2, 1), (1, 0))  #   color bar area 
    
    #ax1.set_position([0.05, 0.15, 0.9, 0.8])
    #ax2.set_position([0.3, 0.025, 0.4, 0.05])        
    
    
    #   specify the position of different parts of the overall figure/movie                        
    fig = plt.figure(figsize = (20,10))
    ax1 = plt.subplot2grid((3, 2), (0, 0), colspan = 2)  #   overall figure
    ax2 = plt.subplot2grid((3, 2), (1, 0))  #   zoomed in coastal zone
    ax3 = plt.subplot2grid((3, 2), (1, 1))  #   zoomed in shelf edge area
    ax4 = plt.subplot2grid((3, 2), (2, 0))  #   color bar area 
    ax5 = plt.subplot2grid((3, 2), (2, 1))  #   text/legend area
    
    ax1.set_position([0.05, 0.475, 0.90, 0.45]) # [left, bottom, width, height]
    ax2.set_position([0.05, 0.15, 0.425, 0.25])
    ax3.set_position([0.525, 0.15, 0.425, 0.25])        
    ax4.set_position([0.15, 0.025, 0.3, 0.05])
    ax5.set_position([0.55, 0.025, 0.3, 0.05])
    
    #   create a plot if desired
    plot_hk_arr = model_obj.hk_arr
    plot_hk_arr[np.abs(plot_hk_arr) == 0.] = np.nan
    #   find all the unique values in the KH array, remove the ones that are marked as no-value
    unique = np.unique(model_obj.hk_arr)
    unique_nan = [value for value in unique if not math.isnan(value)]
    #   define the colormap
    cmap = cmx.viridis
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # define the bins and normalize
    
    #unique_n = [0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 2.5, 5.0] + list(np.arange(10.0, round(max(unique_nan) / 10) * 10, 10.))
    unique_n = [0.0, 0.001, 0.1, 1.0, 5.0, 10., 20., 30.]# + list(np.arange(10.0, round(max(unique_nan) / 10) * 10, 10.))
    norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)            
    
    x_lines_topelev = np.linspace(model_obj.x_start + (model_obj.del_col / (2 * 1000.)) - model_obj.cst_offset_plot, model_obj.x_end - model_obj.cst_offset_plot - (model_obj.del_col / (2 * 1000.)), model_obj.ncol)
    y_lines = np.linspace(model_obj.top, model_obj.zbot, model_obj.nlay + 1)    
    
    #   reverse the list of layer ends at the coastline
    if sed_flux == 'medium' or sed_flux == 'high': 
        cst_lay_pts_plot = cst_lay_pts[::-1]
    else:
        cst_lay_pts_plot = []
    
    def plot_kh_arr(axis, x_start, x_end, y_start, y_end, title):
        im = axis.imshow(plot_hk_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (model_obj.x_start - model_obj.cst_offset_plot, model_obj.x_end - model_obj.cst_offset_plot, model_obj.zbot, math.ceil((model_obj.top / 100.0) * 100.0)), vmin = 0., vmax = np.amax(plot_hk_arr))
        
        lineA, = axis.plot(x_lines_topelev, model_obj.top_elev, c = 'black', linewidth = 3, label = 'GEBCO elevation')
        
        axis.set_xlim([x_start, x_end])    
        axis.set_ylim([y_end, y_start])
        
        axis.axhline(y = 0., linewidth = 1, color = 'k')
        axis.axvline(x = 0., linewidth = 1, color = 'k')
        axis.axvline(x = cont_shelf_edge[0], linewidth = 1, color = 'k')    
            

        #for x_line in end_col_idxs:
        #    axis.axvline(x = x_lines_topelev.tolist()[x_line], linestyle='--', linewidth = .5, color = 'k')            
        #for y_line in end_lay_idxs:
        #    axis.axhline(y = y_lines.tolist()[y_line], linestyle='--', linewidth = .5, color = 'k')  
            
        axis.set_title(title, fontsize = 22)
        axis.set_xlabel('distance from coast (km)', fontsize = 12)
        axis.set_ylabel('elevation (m asl.)', fontsize = 12)
        
        return im
    
    
    txt = ax5.text(0.1, 0.25, 'Qs = %s' % (sed_flux) + '      ' + 'S/M = %s' % (sed_type) + '      ' + 'Geo_sc = %s' % str(inland_aqf_scenario),\
                   fontsize = 20, fontweight = 'bold', clip_on = False)
    ax5.axis('off')
    
    y_zoom_max_cst = math.ceil(max(model_obj.top_elev[int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) : int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) + 100]))
    y_zoom_min_cst = math.floor(min(model_obj.bot_elev[int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) : int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) + 100]))
    y_zoom_max_shelf = math.ceil(max(model_obj.top_elev[shelf_edge_idx - 100 : shelf_edge_idx + 100]))
    y_zoom_min_shelf = math.floor(min(model_obj.bot_elev[shelf_edge_idx - 100 : shelf_edge_idx + 100]))
    
    main_title = 'Horizontal hydraulic conductivity (Hk) distribution in the model_obj domain'
    im1 = plot_kh_arr(ax1, model_obj.x_start - model_obj.cst_offset_plot, model_obj.x_end - model_obj.cst_offset_plot, max(model_obj.top_elev), min(model_obj.bot_elev), main_title)
    im2 = plot_kh_arr(ax2, -5.0, 5.0, y_zoom_max_cst, y_zoom_min_cst, 'Hk, zoom in coastal zone')
    im3 = plot_kh_arr(ax3, cont_shelf_edge[0] - 10.0, cont_shelf_edge[0] + 10.0, max(0.0, y_zoom_max_shelf), y_zoom_min_shelf, 'Hk, zoom in shelf edge zone')
    
    #   plot the colorbar
    cbar = plt.colorbar(im1, cax = ax4, cmap = cmap, norm = norm, spacing = 'uniform', ticks = unique_n, boundaries = unique_n, orientation='horizontal')
    ax4.set_xticklabels([str(e) for e in unique_n])   
    ax4.set_title("Horizontal hydraulic conductivity (m/d)", fontsize = 14)#, y  = 1.025)
    ax4.tick_params(labelsize = 12)    
    
    rcParams['font.family'] = 'Garamond'
    rcParams['axes.facecolor'] = 'white'
    rcParams['savefig.facecolor'] = 'white'  
    
    
    plt.show()
    
    figname = '_Qs_%s' % (sed_flux) + '__SM_%s' % (sed_type) + '__GeoSc_%s' % str(inland_aqf_scenario)
    
    plt.savefig(os.path.join(save_dir, figname))
    plt.close()
    del fig

    return model_obj.hk_arr
    """







