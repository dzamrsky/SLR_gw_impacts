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
    y_lines = np.linspace(model_obj.top, model_obj.zbot, model_obj.nlay + 1)     
    
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
                
    except ValueError:
        print('Not possible to find shelf break based on slope analysis - the shelf decline is too steep probably and matches the slope of the continental slope')
        #   instead, return the value where the bathymetry reaches the 125.m depth 
        shelf_val  = min(topo_lst, key=lambda x:abs(x-(-125.)))     
        shelf_idx = topo_lst.index(shelf_val)        
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
if the upper function fails then try to find either the position of the second landmass on the right side of the model domain.
if that is not found then just give the last column as the shelf break, if there is no point that is lower than -120m asl.

"""

def find_shelf_break_fresh(model_obj, topo_lst, avg = True):
    
    idx_neg = 0
    while topo_lst[idx_neg] < 0.:
        idx_neg += 1

    topo_lst = topo_lst[idx_neg:]

    cst_val  = next((x for x in topo_lst if x < 0.0), None)
    cst_idx = topo_lst.index(cst_val) - 1
    
    sb_point = None
    
    #   1) try to find position of the other land mass, first find the minimum point of the topo list
    idx_min = topo_lst.index(min(topo_lst))
    for i in range(idx_min, len(topo_lst)):
        if topo_lst[i] >= 0.0:
            sb_point = [(i - cst_idx) / 10., topo_lst[i]]
            break

    #   2) if that didnt work then try to find a point below -120. m asl
    if sb_point is None:
        for i in range(idx_min, len(topo_lst)):
            if topo_lst[i] <= -120.0:
                sb_point = [(i - cst_idx) / 10., topo_lst[i]]
                break        

    #   if it still didnt find the shelf break point just give the lowes point then
    if sb_point is None:
        sb_point = [(idx_min - cst_idx) / 10., topo_lst[idx_min]]

    return sb_point


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
save_dir = topo_sc_dir
summary_save_dir = geo_summary_csv_dir
const_geo_hk_vals
clay_cap_shelf
clay_cap_slope
figname
trim_middle = True


"""
 
def create_geology_profile(model_obj, rand_seed_in, inland_aqf_lrs, p_fact, off_lay_thk_ratio, sed_type, lay_pres_y1, clay_cap_shelf_thk, clay_cap_slope_thk,\
                           off_lay_start, sed_flux, save_dir, summary_save_dir, const_geo_hk_vals, clay_cap_shelf, clay_cap_slope, trim_middle, figname):           
           
    random.seed(rand_seed_in)             
           
    """                    First prepare the datasets and lists                     """           
    cont_shelf_edge = find_shelf_break_fresh(model_obj, model_obj.top_elev, True) #   try to find the shelf break              
              
    #   find the coastal index - that will determine the limit between the inland and offshore domain
    #cst_offset_idx_clip = (model_obj.cst_idx - model_obj.idx_start) * 5 + 1  
     
    cst_offset_idx_clip =  int(abs(round(model_obj.x_start, 1) * 10) - 1)
    #cst_offset_val = model_obj.top_elev[cst_offset_idx_clip]
       
    cst_offset_val  = next((x for x in model_obj.top_elev if x < 0.0), None)
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
        if model_obj.hk_vals_bot[x] < 0. or model_obj.hk_vals_bot[x] > 100.:
           model_obj.hk_vals_bot[x] = np.mean([i for i in model_obj.hk_vals_bot if i > 0 and i < 100.])
    
    #   calculate the statistics for both top and bottom aquifer layers
    hk_top_mean = np.nanmean(model_obj.hk_vals_top, axis=0)
    hk_top_std = np.nanstd(model_obj.hk_vals_top, axis=0)
    hk_bot_mean = np.nanmean(model_obj.hk_vals_bot, axis=0)
    hk_bot_std = np.nanstd(model_obj.hk_vals_bot, axis=0)    

    glh_1_val = round(hk_bot_mean, 4)
    glh_2_val = round(hk_top_mean, 4)
    clay_val = 0.0001
    
    #   create the HK_ARR array that will be filled in the next steps
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
        try:
            #   calculate the number of layers that are part of the upper aquifer
            col_top_lay_n = int(round(soilgrids_thk_lst[r] / model_obj.del_lay))
        #   might be that the soilgrids list is shorter than the coastline, in that case take the last value
        except IndexError:
            col_top_lay_n = int(round(soilgrids_thk_lst[-1] / model_obj.del_lay))
        #   fill the upper part of the aquifer in any case (same for both scenarios), the values are randomized
        #   based on the statistics of the GLHYMPS 2.0 list calculated in the previous step
        ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, r].tolist()) if x == 1]
        lrs_lst = []
        
        lrs_lst.append(['glhymps_2', ibound_act_lay_idxs[:col_top_lay_n]]) 
        for s in range(col_top_lay_n):
            #   try except - catch the IndexError in case the SOILGRIDS thickness is larger than the active thickness
            #   of the IBOUND array in the current column
            try:
                if not const_geo_hk_vals:
                    #model_obj.hk_arr[ibound_act_lay_idxs[s], 0, r] = round(abs(np.exp(np.random.normal(glh_2_mu, glh_2_std, size = 1))[0] * 3600 * 24), 4)
                    model_obj.hk_arr[ibound_act_lay_idxs[s], 0, r] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
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
            
            #   decide the layer types in the rest of the model layers based on the majority of glhymps type in that layer
            lay_glh_lst = []

            #   get the pct of thickness per model layer
            lay_thk_pct = round(1 / float(tot_cells) * 100., 0)

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
    
    
            if lay_glhymps_lst != [[]]:
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


            """
            if tot_cells <= 100:
                #   get the pct of thickness per model layer
                lay_thk_pct = round(1 / float(tot_cells) * 100., 0)
                
                #   decide the layer types in the rest of the model layers based on the majority of glhymps type in that layer
                lay_glh_lst = []
      
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
                
            
            #   make sure the list above has the right amount of elements (in case the % is 33.0 for example it appends one extra value to reach length of 100)
            if len(lay_glhymps_lst) > tot_cells:
                lay_glhymps_lst = lay_glhymps_lst[:-1]

            else:
                glh_type_tofill = max([sublist for sublist in inland_aqf_lrs])[1]
                if glh_type_tofill == 'glhymps_2': 
                    last_non_empty = 2
                elif glh_type_tofill == 'glhymps_1': 
                    last_non_empty = 1            
            
            for x in range(len(lay_glhymps_lst)):
                if lay_glhymps_lst[x] == []:
                    lay_glhymps_lst[x] = [last_non_empty]
                    
            """

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
                            model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, r] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
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
                                model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l - 1], 0, r] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
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
                                model_obj.hk_arr[ibound_act_lay_idxs[-1], 0, r] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
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
    try:
        lay_thk_1 = round((int(soilgrids_thk_lst[r] / model_obj.del_lay) / float(tot_cells_cst)) * 100., 2) 
    except IndexError:
        lay_thk_1 = round((int(soilgrids_thk_lst[-1] / model_obj.del_lay) / float(tot_cells_cst)) * 100., 2) 
    lay_thk.append([lay_thk_1, 'glhymps_2'])
    for w in range(len(inland_aqf_lrs)):
        geo_lay_thk = (inland_aqf_lrs[w][0] * tot_cells) / 100
        lay_thk.append([round((geo_lay_thk / float(tot_cells_cst)) * 100., 2), inland_aqf_lrs[w][1]])
    lay_thk[-1][0] += round(100. - sum(i[0] for i in lay_thk), 2)


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
    #shelf_edge_idx = model_obj.top_elev.index(cont_shelf_edge[1]) 
    shelf_edge_idx = int(abs(model_obj.x_start) + cont_shelf_edge[0]) * 10
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
    
    glhymps_2_top = False
    
    #   find the lowest elevation in the offshore domain, this will be the axis of symetry for the geological layers, also will change
    #   the bottom elevation of the model domain, i.e. it will be also the deepest point of the offshore sediments. There is no need
    #   to try to find the continental shelf edge as these models usually do not have one and only consist of a very shallow sea. 

    #   first calculate the average slopes of the cont. shelf and cont. slope, as the elevation difference between the shelf edge and the
    #   elevation at the coastline
    avg_sl_shelf = round(100 * (model_obj.top_elev[cst_idx_geo] - model_obj.top_elev[shelf_edge_idx]) / ((shelf_edge_idx - cst_idx_geo) * 100.), 2)
    avg_sl_slope = round(100 * (model_obj.top_elev[shelf_edge_idx] - model_obj.top_elev[-1]) / ((model_obj.ncol - shelf_edge_idx) * 100.), 2)
        
    #   get the top points of each layer (not the top one..), different for each scenario
    cst_lay_pts = []
    
    #   calculate the thickness of the top glhymps_2 layer at the coastline, to be expanded into the offshore domain
    col_top_lay_n = int(round(soilgrids_thk_lst[-1] / model_obj.del_lay))     
    #   get the active layer indexes at the coastline
    cst_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, cst_offset_idx].tolist()) if x == 1]
        
    #   no matter the sediment flux quantity, the shape of the geological layers will be the same for this model concept
    #   fill the are between the coast and the min_elev_idx with the glhymps layers based on the % fractions..
    #   look for the lowest elevation between the coastal point and the end of model domain
    min_elev_val = min(model_obj.top_elev[cst_idx_geo:])
    min_elev_idx = model_obj.top_elev[cst_idx_geo:].index(min_elev_val) + cst_idx_geo
    min_elev_ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, min_elev_idx].tolist()) if x == 1]

    #   calculate the amount of layers at the continental shelf column
    tot_layers = len(inland_aqf_lrs)
    #   divide the area between the shelf edge and model_obj bottom into the tot_layers identical amount of cells

    #   find the limits for each of the geological layers 
    idx_1 = min_elev_ibound_act_lay_idxs[col_top_lay_n]                
    end_lay_idxs =  [idx_1]
    act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[idx_1, 0, :].tolist()) if x == 1]
    end_col_idxs = [min_elev_idx - 1]                
    for b in range(tot_layers - 1):                    
        lay_idx = idx_1 + int(round((inland_aqf_lrs[b][0] * (min_elev_ibound_act_lay_idxs[-1] - min_elev_ibound_act_lay_idxs[col_top_lay_n])) / 100, 0))              
        end_lay_idxs.append(lay_idx)
        act_cell_lst = [i for i, x in enumerate(model_obj.ibound_arr[lay_idx - 1, 0, :].tolist()) if x == 1]
        end_col_idxs.append(min_elev_idx - 1)
        idx_1 = lay_idx
    #   append the last active layer and column as the end of the last layer
    end_lay_idxs.append(min_elev_ibound_act_lay_idxs[-1])
    end_col_idxs.append(min_elev_idx - 1)
    
    """ maybe this should be changed in the other script as well - model.nlay to min_elev_ibound_act_lay_idxs[-1]  """
    
    for a in range(cst_idx_geo, min_elev_idx):        
        
        #   the upper part of the offshore aquifer domain
        ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, a].tolist()) if x == 1]
        lrs_lst = []  
        """     add if glhymps_top is true or false     """                    
        
        lrs_lst.append(['glhymps_2', ibound_act_lay_idxs[:col_top_lay_n]]) 
        for s in range(col_top_lay_n):
            #   try except - catch the IndexError in case the SOILGRIDS thickness is larger than the active thickness
            #   of the IBOUND array in the current column
            try:
                if not const_geo_hk_vals: 
                    model_obj.hk_arr[ibound_act_lay_idxs[s], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                else:
                    model_obj.hk_arr[ibound_act_lay_idxs[s], 0, a] = glh_2_val
            except IndexError:
                pass    
            
        #   split the active ibound cells in the column list into groups (geological layers) based on the % from the input list
        tot_cells = len(ibound_act_lay_idxs) - col_top_lay_n        
    
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


        if lay_glhymps_lst != [[]]:
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
                    model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, a] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)                        
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

        offshore_lrs_lst.append([a, lrs_lst])



    #   now comes the part where the IBOUND, model_obj.zbot and model_obj.bot_elev lists are going to be changed and adapted
    #   the slope will be the same as in the other part of the model. Do this in sort of a mirror way
    col_offset = 0
    for b in range(min_elev_idx, model_obj.ibound_arr.shape[-1]):
        
        #   get the active layers in the column, but limit the bottom the the column on the left of the lowest elevation index
        ibound_act_lay_idxs = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, b].tolist()) if x == 1]
        lrs_lst = []  
        
        if ibound_act_lay_idxs != []:
            
            #   now get the extent of the column in the left side of the model domain and use that as a bottom limit
            bot_lim_lay_idx = [i for i, x in enumerate(model_obj.ibound_arr[:, 0, min_elev_idx - col_offset].tolist()) if x == 1][-1]
            
            if bot_lim_lay_idx in ibound_act_lay_idxs:
                ibound_act_lay_idxs = ibound_act_lay_idxs[:ibound_act_lay_idxs.index(bot_lim_lay_idx) + 1]
    
            #   set the rest of the IBOUND array to 0
            model_obj.ibound_arr[ibound_act_lay_idxs[-1]:, :, b] = 0
    
            lrs_lst.append(['glhymps_2', ibound_act_lay_idxs[:col_top_lay_n]]) 
            for s in range(col_top_lay_n):
                #   try except - catch the IndexError in case the SOILGRIDS thickness is larger than the active thickness
                #   of the IBOUND array in the current column
                try:
                    if not const_geo_hk_vals: 
                        model_obj.hk_arr[ibound_act_lay_idxs[s], 0, b] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)
                    else:
                        model_obj.hk_arr[ibound_act_lay_idxs[s], 0, b] = glh_2_val
                except IndexError:
                    pass    
                
            #   split the active ibound cells in the column list into groups (geological layers) based on the % from the input list
            tot_cells = len(ibound_act_lay_idxs) - col_top_lay_n        
        
            if tot_cells != 0:
        
                #   get the pct of thickness per model layer
                lay_thk_pct = round(1 / float(tot_cells) * 100., 0)
                
                #   decide the layer types in the rest of the model layers based on the majority of glhymps type in that layer
                lay_glh_lst = []
          
                #   create a list of 100 values (1 or 2 depending on the glhymps type)
                for i in range(len(inland_aqf_lrs)):
                    glh_pct = inland_aqf_lrs[i][0]
                    glh_type = inland_aqf_lrs[i][1]                   
                    for j in range(glh_pct):
                        if glh_type == 'glhymps_1':
                            lay_glh_lst.append(1)
                        else:
                            lay_glh_lst.append(2)                    
                
                #   sometimes the lay_thk_pct is close to 0%, if that is the case set it to 1
                if lay_thk_pct <= 0.0:
                    lay_thk_pct = 1.0
                
                #   split this list into even chunks based on the lay_thk_pct value
                lay_glhymps_lst = [lay_glh_lst[i:i + int(lay_thk_pct)] for i in range(0, len(lay_glh_lst), int(lay_thk_pct))]
                
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
                            model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, b] = round(abs(np.random.normal(hk_bot_mean, hk_bot_std)), 4)
                        elif glh_typ == 2:
                            model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, b] = round(abs(np.random.normal(hk_top_mean, hk_top_std)), 4)                        
                    else:
                        if glh_typ == 1:
                            model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, b] = glh_1_val
                        elif glh_typ == 2:
                            model_obj.hk_arr[ibound_act_lay_idxs[col_top_lay_n + l], 0, b] = glh_2_val  
                       
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
                
        col_offset += 1

        offshore_lrs_lst.append([b, lrs_lst])

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
    off_lay_idx_start = off_lays_idxs.index(off_sed_start_idx)

    end_col_glhymps_1_start = end_col_idxs[1:][::2]
    end_col_glhymps_1_end = end_col_idxs[2:][::2][:-1]
    end_col_glhymps = []
    for d in range(len(end_col_glhymps_1_end)):
        end_col_glhymps.append([end_col_glhymps_1_start[d], end_col_glhymps_1_end[d]])
    end_col_glhymps.append([model_obj.ncol, model_obj.ncol])
        
    #   counter for the different end layers 
    end_lay_cnt = 0

    #   loop through the columns where the sed_type adjustment will be implemented        
    for i in range(off_lay_idx_start, len(tot_clay_lst)):
        #   check for the 'glhymps_1 columns', select those from the offshore_lrs_lst
        col_idx = tot_clay_lst[i][0]
        
        try:
            #   select only the parts where there is 'glhymps_1' sediment type
            lay_idxs = [item for item in tot_clay_lst[i][1] if item[0] == 'glhymps_1']
            #   get the total amount of cells that should have the clay/silt properties, based on the % filled in the off_lay_thk_ratio                    
            lay_lst_1 = [item[1] for item in tot_clay_lst[i][1] if item[0] == 'glhymps_1']
            #lay_lst = [item for sublist in lay_lst_1 for item in sublist]
            """ Maybe in future substract the clay capping cells from the ones that are then changed in the column.."""
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
                    
    if trim_middle:
        model_obj.hk_arr = model_obj.hk_arr[:, :, : min_elev_idx]

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
    for i in range(plot_hk_arr.shape[0]):
        for j in range(plot_hk_arr.shape[-1]):
            if model_obj.ibound_arr[i, 0, j] == 1:
                if np.isnan(plot_hk_arr[i, 0, j]):
                    #print(i, 0, j, round(np.nanmean(model_obj.hk_arr), 6))    
                    plot_hk_arr[i, 0, j] = round(np.nanmean(model_obj.hk_arr), 6)

    #   trim the final array, check for empty columns and layers at left, right, top, bottom
    left_lim, top_lim = 0, 0
    
    for i in range(plot_hk_arr.shape[0]):
        act_cells_lst = [i for i, x in enumerate(model_obj.ibound_arr[i, 0, :].tolist()) if x == 1][:cap_cells_shelf]
        if len(act_cells_lst) == 0:
            top_lim += 1
        else:
            break
    bot_lim = int(top_lim)
    for i in range(top_lim, plot_hk_arr.shape[0]):
        act_cells_lst = [i for i, x in enumerate(model_obj.ibound_arr[i, 0, :].tolist()) if x == 1][:cap_cells_shelf]
        if len(act_cells_lst) != 0 and len(act_cells_lst) != []:
            bot_lim += 1
        else:
            break    
    
    for j in range(plot_hk_arr.shape[-1]):
        act_cells_lst = [j for j, x in enumerate(model_obj.ibound_arr[:, 0, j].tolist()) if x == 1][:cap_cells_shelf]
        if len(act_cells_lst) == 0:
            left_lim += 1
        else:
            break
    right_lim = int(left_lim)
    for j in range(top_lim, plot_hk_arr.shape[-1]):
        act_cells_lst = [j for j, x in enumerate(model_obj.ibound_arr[:, 0, j].tolist()) if x == 1][:cap_cells_shelf]
        if len(act_cells_lst) != 0:
            right_lim += 1
        else:
            break       
    
    #   check the right hand side, can happen that there is a weird looking full column
    for g in range(min_elev_idx + 1, model_obj.ibound_arr.shape[-1] - 1):
        try:
            last_lay_left = [j for j, x in enumerate(model_obj.ibound_arr[:, 0, g - 1].tolist()) if x == 1][-1]
            last_lay_right = [j for j, x in enumerate(model_obj.ibound_arr[:, 0, g].tolist()) if x == 1][-1]
            if last_lay_right > last_lay_left or last_lay_right == []:
                break
        #   happens when the whole column is inactive
        except IndexError:
            break
    
    right_lim = g
    
    plot_hk_arr = plot_hk_arr[top_lim : bot_lim, :, left_lim : right_lim]
    
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
    
    #x_end_new = round(model_obj.x_start - model_obj.cst_offset_plot, 1) + model_obj.ncol / 10.
    
    x_end_new = round(model_obj.x_start - model_obj.cst_offset_plot, 1) + plot_hk_arr.shape[-1] / 10.
    x_start_new = round(model_obj.x_start - model_obj.cst_offset_plot, 1)
    x_lines_topelev = np.linspace(x_start_new,  x_end_new, model_obj.ncol + 1)
    y_lines = np.linspace(model_obj.top, model_obj.zbot, model_obj.nlay + 1)    
    
    print('x_end_new = ' + str(x_end_new))
    
    #   reverse the list of layer ends at the coastline
    cst_lay_pts_plot = []
    
    def plot_kh_arr(axis, x_start, x_end, y_start, y_end, title):
        #im = axis.imshow(plot_hk_arr[:, 0, :], aspect = 'auto', interpolation = None, cmap = cmap, norm = norm,
        #               extent = (x_start_new, x_end, model_obj.zbot, math.ceil((model_obj.top / 100.0) * 100.0)), vmin = 0., vmax = np.nanmax(plot_hk_arr))

        im = axis.imshow(plot_hk_arr[:, 0, :], aspect = 'auto', interpolation = None, cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end, 10 * ((model_obj.top / 100.0) * 10.0 - bot_lim), math.ceil((model_obj.top / 100.0) * 100.0)),\
                       vmin = 0., vmax = np.nanmax(plot_hk_arr))
        
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
    
    txt = ax5.text(0.1, 0.25, 'Qs = %s' % (sed_flux) + '                    ' + 'pres_fact = %s' % str(round(lay_pres_y1, 1)) + '                   '\
                   + 'stacking factor P = %s' % str(round(p_fact, 2)) + '\n' + 'mud_pct = %s' % str(round(off_lay_thk_ratio, 1)) + '           ' + 'clay cap slope/shelf (thk in m) = %s, %s (%s, %s)' %\
                   (clay_cap_shelf_str, clay_cap_slope_str, clay_cap_shelf_thk_str, clay_cap_slope_thk_str), fontsize = 14, fontweight = 'bold', clip_on = False)
    ax5.axis('off')
    
    y_zoom_max_cst = math.ceil(max(model_obj.top_elev[int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) : int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) + 100]))
    y_zoom_min_cst = math.floor(min(model_obj.bot_elev[int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) : int(abs(model_obj.x_start - model_obj.cst_offset_plot + 5.0) * 10) + 100]))
    y_zoom_max_shelf = math.ceil(max(model_obj.top_elev[shelf_edge_idx - 100 : shelf_edge_idx + 100]))
    y_zoom_min_shelf = math.floor(min(model_obj.bot_elev[shelf_edge_idx - 100 : shelf_edge_idx + 100]))
    
    main_title = 'Horizontal hydraulic conductivity (Hk) distribution in the model domain of COSCAT %s, coastal type %s, realization %s' %(str(model_obj.coscat_id), model_obj.cst_type, str(rand_seed_in))
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
    

    
    #   trim the output array 
    final_array = plot_hk_arr[:, :, :]       
    
    return final_array, top_lim, bot_lim, left_lim, right_lim, min_elev_idx
    
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







