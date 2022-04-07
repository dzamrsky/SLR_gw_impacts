# -*- coding: utf-8 -*-
"""
    Toolbox for creating 2D cross-sectional SEAWAT models, using the Flopy library
"""

import math
import numpy as np
import flopy
import flopy.utils.binaryfile as bf
import os
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors
import matplotlib.patches as mpatches
import matplotlib
from time import time
import random
import matplotlib.animation as animation
from matplotlib import rcParams
#from matplotlib.ticker import AutoMinorLocator
import itertools
import matplotlib as mpl
import matplotlib.cm as cmx
#from collections import OrderedDict

#   first create necessary dictionaries (with K values for soil and GLIM classes)
k_soil_dict = dict()
for soil_id in np.arange(1, 12, 1):
    k_soil_dict[soil_id] = dict()
    k_soil_dict[soil_id]['soil name'] = ['']
    k_soil_dict[soil_id]['k_val'] = []
    k_soil_dict[soil_id]['plot_cl'] = ['']

#   first define dictionary of hydraulic conductivity values
k_soil_dict[1] = {'soil_name': 'Clay', 'k_val': 0.01123, 'plot_cl': 'brown'}
k_soil_dict[2] = {'soil_name': 'Silty clay', 'k_val': 0.0864, 'plot_cl': 'chocolate'}
k_soil_dict[3] = {'soil_name': 'Silty clay loam', 'k_val': 0.1469, 'plot_cl': 'crimson'}
k_soil_dict[4] = {'soil_name': 'Clay loam', 'k_val': 0.216, 'plot_cl': 'sienna'}
k_soil_dict[5] = {'soil_name': 'Silt loam', 'k_val': 0.622, 'plot_cl': 'tan'}
k_soil_dict[6] = {'soil_name': 'Sandy clay', 'k_val': 0.19, 'plot_cl': 'tbd'}
k_soil_dict[7] = {'soil_name': 'Sandy clay loam', 'k_val': 0.5443, 'plot_cl': 'orangered'}
k_soil_dict[8] = {'soil_name': 'Sandy loam', 'k_val': 2.946, 'plot_cl': 'orange'}
k_soil_dict[9] = {'soil_name': 'Sand', 'k_val': 15.2064, 'plot_cl': 'yellow'}
k_soil_dict[10] = {'soil_name': 'Loamy sand', 'k_val': 13.5043, 'plot_cl': 'coral'}
k_soil_dict[11] = {'soil_name': 'Loam', 'k_val': 0.6048, 'plot_cl': 'wheat'}
# and you can save them like:
#np.save(fullpath.npy, my_soil)

"""
k_glim_dict = dict()
glim_id_lst = np.arange(-1, 16, 1)
glim_id_lst = np.delete(glim_id_lst, [1])

for glim_id in glim_id_lst:
    k_glim_dict[glim_id] = dict()
    k_glim_dict[glim_id]['glim name'] = ['']
    k_glim_dict[glim_id]['k_val'] = []
    k_soil_dict[soil_id]['plot_cl'] = ['']

#   first define dictionary of hydraulic conductivity values
k_glim_dict[-1] = {'glim_name': 'NoData', 'k_val': -1.}
k_glim_dict[1] = {'glim_name': 'su - Unconsolidated sediments', 'k_val': 8.459E-02}
k_glim_dict[2] = {'glim_name': 'ss - Siliciclastic sedimentary rocks', 'k_val': 5.337E-04}
k_glim_dict[3] = {'glim_name': 'py - Pyroclastics', 'k_val': }
k_glim_dict[4] = {'glim_name': 'sm - Mixed sedimentary rocks', 'k_val': }
k_glim_dict[5] = {'glim_name': 'sc - Carbonate sedimentary rocks', 'k_val': }
k_glim_dict[6] = {'glim_name': 'ev - Evaporites', 'k_val': }
k_glim_dict[7] = {'glim_name': 'va - Acid volcanic rocks', 'k_val': }
k_glim_dict[8] = {'glim_name': 'vi - Intermediate volcanic rocks', 'k_val': }
k_glim_dict[9] = {'glim_name': 'vb - Basic volcanic rocks', 'k_val': }
k_glim_dict[10] = {'glim_name': 'pa - Acid plutonic rocks', 'k_val': }
k_glim_dict[11] = {'glim_name': 'pi - Intermediate plutonic rocks', 'k_val': }
k_glim_dict[12] = {'glim_name': 'pb - Basic plutonic rocks', 'k_val': }
k_glim_dict[13] = {'glim_name': 'mt - Metamorphics', 'k_val': }
k_glim_dict[14] = {'glim_name': 'wb - Water bodies', 'k_val': }
k_glim_dict[15] = {'glim_name': 'ig - Ice and glaciers', 'k_val': }
k_glim_dict[16] = {'glim_name': 'nd - No data', 'k_val': }
# and you can save them like:
#np.save(fullpath.npy, my_soil)
"""

#   Define a model class including all necessary methods associated with creating
#   2D cross-sectional models using SEAWAT packages via flopy
class cs_model(object):

    #   initialize the cross-section object and assign necessary variables (topography + bathymetry, Pelletier 2016 thickness etc.)
    def __init__(self, id_cs, topo_vals, nasa_vals, nasa_x, nasa_y, cst_thick_est, cst_plain_end, soil_thk_vals, soil_type_vals,\
                 wtd_depth_vals, offshore_sed_vals, pcr_rch_vals, watergap_rch_vals, p_min_et_vals, k_soil_vals, drn_rate_vals,\
                 glhymps_top_lay_vals, glhymps_bot_lay_vals, topo_res, cst_bound = False):
        self.id_cs = id_cs
        self.topo_vals = topo_vals
        self.nasa_vals = nasa_vals
        self.anchor_x  = nasa_x
        self.anchor_y  = nasa_y
        self.cst_thick = cst_thick_est
        self.cst_plain_end = cst_plain_end
        self.cst_bound = cst_bound      #   indicates if the seaward boundary should be placed at coastline
        self.soil_thk = soil_thk_vals
        self.soil_type = soil_type_vals
        self.wtd_vals = wtd_depth_vals
        self.off_sed = offshore_sed_vals
        self.pcr_rch = pcr_rch_vals
        self.watergap_rch = watergap_rch_vals
        self.p_min_et = p_min_et_vals
        self.k_soil = k_soil_vals
        self.drn_dens = drn_rate_vals
        self.glhymps_top_lay = glhymps_top_lay_vals
        self.glhymps_bot_lay = glhymps_bot_lay_vals
        self.topo_res = topo_res
        self.mid_idx = int((len(self.topo_vals) - 1) / 2)

    #   method that processes the input data - extends the original lists
    def pre_process_data(self, del_col):
        self.vals_topo_pts_x = np.linspace(-200., 200., 801)         #   keep constant for now
        self.vals_topo_2_pts_x, self.vals_topo_2_pts_y = [], []
        self.vals_soil_thk_2_pts_x, self.vals_soil_thk_2_pts_y = [], []
        self.vals_wtd_depth_2_pts_x, self.vals_wtd_depth_2_pts_y = [], []

        n = 5                   #   number of points to be inserted between the topographic points (= nr of columns)
        delcol = del_col / 1000.
        #   loop through the original list
        for i in range(1, len(self.topo_vals)):
            #   calculate the gradient between the two topographic points
            gradient_topo = ( self.topo_vals[i - 1] - self.topo_vals[i] ) / float(n)
            gradient_soil = ( self.soil_thk[i - 1] - self.soil_thk[i] ) / float(n)
            gradient_wtd = ( self.wtd_vals[i - 1] - self.wtd_vals[i] ) / float(n)
            for j in range(n):
                #   calculate the x and y value for each model column (middle of the cell)
                self.vals_topo_2_pts_x.append(round(self.vals_topo_pts_x[i - 1] + j * delcol + delcol / 2., 2))
                self.vals_soil_thk_2_pts_x.append(round(self.vals_topo_pts_x[i - 1] + j * delcol + delcol / 2., 2))
                self.vals_wtd_depth_2_pts_x.append(round(self.vals_topo_pts_x[i - 1] + j * delcol + delcol / 2., 2))
                self.vals_topo_2_pts_y.append(round(self.topo_vals[i - 1] - j * gradient_topo - gradient_topo / 2., 1))
                self.vals_soil_thk_2_pts_y.append(round(self.soil_thk[i - 1] - j * gradient_soil - gradient_soil / 2., 1))
                self.vals_wtd_depth_2_pts_y.append(round(self.wtd_depth[i - 1] - j * gradient_wtd - gradient_wtd / 2., 1))

    #   Method to get the top and bottom elevation lists for the whole model domain, raw data from dbase
    #   the difference here is that the position of the coastline is a parameter for the function
    #   cst_idx_shift = index that indicates how far from the middle of the cross-section topo_lst
    #                   we look for the position of the coastline
    def get_top_bot_lst_find_cst(self, cst_idx_shift, fos_point, sea_depth_val, wtd_cst = False):
        #   create output lists
        out_top_lst, out_bot_lst, out_wtd_lst = [], [], []

        #   first check if there is any land in the zone around the middle of the topo list (it should be around 400th index)
        cst_topo_lst = self.topo_vals[cst_idx_shift : self.mid_idx + (self.mid_idx - cst_idx_shift)]
        
        #   check if there are any islands or other landmass in the coastal stretch - the length of uninterupted below sea level elevations has to be more than 1km and deeper than 5m
        import itertools
        vals = self.topo_vals[cst_idx_shift : self.mid_idx + (self.mid_idx - cst_idx_shift)]
        splitted = [list(g) for i, g in itertools.groupby(vals,lambda x: x<0)]
        #   check that the first list is positive numbers
        try:
            if len(splitted) <= 2 and splitted[0][0] > 0:
                #   if all values are negative then try to find the parts with values above sea level
                if all(i < sea_depth_val for i in cst_topo_lst):
                    #   if it is true then find the land part in the left part of the self domain (assuming the right one is sea/ocean anyway)
                    pos_vals = [x for x in self.topo_vals[:self.mid_idx] if x > 0] 
                    pos_vals_idx = [self.topo_vals.index(pos_val) for pos_val in pos_vals]
                    self.topo_vals[pos_vals_idx[0] : pos_vals_idx[-1] + 1]
                    #   find the coastline index, take the last land cell index
                    cst_offset_idx = pos_vals_idx[-1] + 1
                    cst_offset_val = self.topo_vals[cst_offset_idx]
                    self.cst_offset = (cst_offset_idx - self.mid_idx) * self.topo_res / 1000.
                    self.cst_idx = cst_offset_idx 
        
                else:
                    try:
                        #   first check where is the bathymetry location which is lower than that sea level
                        cst_offset_val_deep  = next((x for x in self.topo_vals[cst_idx_shift:] if x < -32.), None)
                        cst_offset_idx_deep = self.topo_vals[cst_idx_shift:].index(cst_offset_val_deep) - 1
                        #   second, go from there in the opposite direction (landward) and look for first cell above current sea level (0)
                        cst_offset_val  = next((x for x in list(reversed(self.topo_vals[:cst_offset_idx_deep + cst_idx_shift])) if x >= sea_depth_val), None)
                        cst_offset_idx = self.topo_vals[cst_idx_shift:].index(cst_offset_val) 
                        #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                        self.cst_offset = ((cst_idx_shift + cst_offset_idx) - self.mid_idx) * self.topo_res / 1000.
                        self.cst_idx = cst_idx_shift + cst_offset_idx   
                    except ValueError:
                    
                        cst_offset_val = next((x for x in self.topo_vals[cst_idx_shift:] if x < sea_depth_val), None)
                        cst_offset_idx = self.topo_vals[cst_idx_shift:].index(cst_offset_val) - 1
                        self.cst_offset = ((cst_idx_shift + cst_offset_idx) - self.mid_idx) * self.topo_res / 1000.
                        self.cst_idx = cst_idx_shift + cst_offset_idx 

            if len(splitted) > 2 and splitted[0][0] > 0 and splitted[1][0] < 0:
                if all(i < sea_depth_val for i in cst_topo_lst):
                    #   if it is true then find the land part in the left part of the self domain (assuming the right one is sea/ocean anyway)
                    pos_vals = [x for x in self.topo_vals[:self.mid_idx] if x > 0] 
                    pos_vals_idx = [self.topo_vals.index(pos_val) for pos_val in pos_vals]
                    self.topo_vals[pos_vals_idx[0] : pos_vals_idx[-1] + 1]
                    #   find the coastline index, take the last land cell index
                    cst_offset_idx = pos_vals_idx[-1] + 1
                    cst_offset_val = self.topo_vals[cst_offset_idx]
                    self.cst_offset = (cst_offset_idx - self.mid_idx) * self.topo_res / 1000.
                    self.cst_idx = cst_offset_idx 
        
                else:
                    try:
                        #   first check where is the bathymetry location which is lower than that sea level
                        cst_offset_val_deep  = next((x for x in self.topo_vals[cst_idx_shift : cst_idx_shift + len(splitted[0]) + len(splitted[1])] if x <= min(splitted[1])), None)
                        cst_offset_idx_deep = self.topo_vals[cst_idx_shift : cst_idx_shift + len(splitted[0]) + len(splitted[1])].index(cst_offset_val_deep) - 1
                        #   second, go from there in the opposite direction (landward) and look for first cell above current sea level (0)
                        cst_offset_val  = next((x for x in list(reversed(self.topo_vals[:cst_offset_idx_deep + cst_idx_shift])) if x >= sea_depth_val), None)
                        cst_offset_idx = self.topo_vals[cst_idx_shift:].index(cst_offset_val) 
                        #cst_offset = ((cst_look_up_idx + cst_offset_idx) - 400) * 0.5 
                        self.cst_offset = ((cst_idx_shift + cst_offset_idx) - self.mid_idx) * self.topo_res / 1000.
                        self.cst_idx = cst_idx_shift + cst_offset_idx   
                    except ValueError:
                    
                        cst_offset_val = next((x for x in self.topo_vals[cst_idx_shift:] if x < sea_depth_val), None)
                        cst_offset_idx = self.topo_vals[cst_idx_shift:].index(cst_offset_val) - 1
                        self.cst_offset = ((cst_idx_shift + cst_offset_idx) - self.mid_idx) * self.topo_res / 1000.
                        self.cst_idx = cst_idx_shift + cst_offset_idx 
                        
        except IndexError:
            cst_offset_val = next((x for x in self.topo_vals[self.mid_idx:] if x < sea_depth_val), None)
            cst_offset_idx = self.topo_vals.index(cst_offset_val) - 1
            self.cst_offset = (cst_offset_idx - self.mid_idx) * self.topo_res / 1000.
            self.cst_idx = cst_offset_idx 

        if wtd_cst == True:
            for g in range(self.mid_idx, len(self.wtd_vals)):
                if math.isnan(self.wtd_vals[g]):
                    print(g)
                    self.cst_offset = ((cst_idx_shift + cst_offset_idx) - self.mid_idx) * self.topo_res / 1000.
                    self.cst_idx = g
                    break
        
        if not hasattr(self, 'cst_idx'):
            self.cst_offset = 0.
            self.cst_idx = self.mid_idx
        
        #   get index of coastal plain end/start 2 because the distance between cross-section points is 0.5km
        idx_start = self.cst_idx - int((1000 / self.topo_res) * self.cst_plain_end)
        anchor_x_idx = self.cst_idx - int((1000 / self.topo_res) * self.anchor_x)         #   find the x index of the anchor point
        sed_depth = self.anchor_y 
        
        #   append the depth indicated by the nasa dataset to the bottom list
        for a in range(idx_start, anchor_x_idx):
            #   in some cases (usually the 'min_min' scenarios the elevation can reach negative values on both side)
            #   if that happens then cut the self from the left (inland) side..
            if self.topo_vals[a] < 0:
                idx_start += 1
                pass
            else:
                out_bot_lst.append(self.topo_vals[a] - self.nasa_vals[a])
        #   calculate the increment per 0.5km (distance between the cross-section points) for the bedrock slope
        depth_incr_ate = abs(self.anchor_y - (- self.cst_thick)) / (self.cst_idx - anchor_x_idx)

        #   now complete the bottom list with values between the coast and the FOS point (if there is one)
        if fos_point is not None:

            #   form the anchor point till the coastline follow the gradient between these two points
            for b in range(anchor_x_idx, self.cst_idx):
                sed_depth = sed_depth - depth_incr_ate #   add the depth_incr to the coastal thickness estimation
                out_bot_lst.append(sed_depth)

            fos_x_from_cst = fos_point[0]
            #   fos_x_idx = self.cst_idx + int(fos_x_from_cst / 0.5)
            fos_x_idx = self.mid_idx + int(fos_x_from_cst / (self.topo_res / 1000.)) # the fos point is estimated before the coastline is shifted, thats why 400
            fos_y = fos_point[1]
            #   calculate the increment between the coastline and the fos point
            depth_incr_fos = abs(- self.cst_thick - fos_y) / (fos_x_idx - self.cst_idx)
            sed_depth = - self.cst_thick
            #   loop through the points located between coast and fos
            for c in range(self.cst_idx, fos_x_idx):
                sed_depth = sed_depth - depth_incr_fos #   add the depth_incr_fos to the coastal thickness estimation
                if sed_depth > self.topo_vals[c]:
                    idx_end = c #- 2 #   -2 to get the right index
                    break
                else:
                    out_bot_lst.append(sed_depth)
                    idx_end = fos_x_idx
                    continue
            
            print ('idx_start, idx_end :' + str(idx_start) + ', ' + str(idx_end))
            out_top_lst = self.topo_vals[idx_start : idx_end]   #   list with top elevations is only a clip from the overall topography
            out_wtd_lst = self.wtd_vals[idx_start : idx_end]
            
            out_pcr_rch_lst = self.pcr_rch[idx_start : idx_end]
            out_watergap_rch_lst = self.watergap_rch[idx_start : idx_end]
            out_p_min_et_lst = self.p_min_et[idx_start : idx_end]
            out_k_soil = self.k_soil[idx_start : idx_end]
            out_drn_dens = self.drn_dens[idx_start : idx_end]
            out_glhymps_top_lay = self.glhymps_top_lay[idx_start : idx_end]
            out_glhymps_bot_lay = self.glhymps_bot_lay[idx_start : idx_end]
            out_soil_thk = self.soil_thk[idx_start : idx_end]

        else:
            print('No FOS point detected..')
            for a in range(anchor_x_idx, len(self.topo_vals) - 1):
                sed_depth = sed_depth - depth_incr_ate              #   add the depth_incr to the coastal thickness estimation
                #print ('sed_depth, topo_vals[a] :' + str(sed_depth) + ', ' + str(self.topo_vals[a]))
                if sed_depth > self.topo_vals[a]:
                    idx_end = a #- 2 #   -2 to get the right index
                    break
                else:
                    out_bot_lst.append(sed_depth)
                    idx_end = len(self.topo_vals) - 1
                    continue
            
            print ('idx_start, idx_end :' + str(idx_start) + ', ' + str(idx_end))
            out_top_lst = self.topo_vals[idx_start : idx_end]   #   list with top elevations is only a clip from the overall topography
            out_wtd_lst = self.wtd_vals[idx_start : idx_end]
            out_pcr_rch_lst = self.pcr_rch[idx_start : idx_end]
            out_watergap_rch_lst = self.watergap_rch[idx_start : idx_end]
            out_p_min_et_lst = self.p_min_et[idx_start : idx_end]
            out_k_soil = self.k_soil[idx_start : idx_end]
            out_drn_dens = self.drn_dens[idx_start : idx_end]
            out_glhymps_top_lay = self.glhymps_top_lay[idx_start : idx_end]
            out_glhymps_bot_lay = self.glhymps_bot_lay[idx_start : idx_end]
            out_soil_thk = self.soil_thk[idx_start : idx_end]
            
        
        #   assign the final lists
        self.top_elev = out_top_lst
        self.bot_elev = out_bot_lst
        self.idx_start = idx_start
        self.idx_end = idx_end
        self.wtd_elev = out_wtd_lst
        self.rch_pcr = out_pcr_rch_lst
        self.rch_watergap = out_watergap_rch_lst
        self.rch_p_min_et = out_p_min_et_lst
        self.hk_soil = out_k_soil
        self.drn_density = out_drn_dens
        self.glhymps_top_lay = out_glhymps_top_lay
        self.glhymps_bot_lay = out_glhymps_bot_lay
        self.soil_thk = out_soil_thk
        
        
    #   Method that gives a list of top and bottom elevation for each column depending on the del_col
    #   Added an option for smoothing out the topography depending on the ratio between the del_col and the
    #   distance between the cross-section points. THE TOP AND BOT LISTS MUST BE CREATED BEFORE RUNNING THIS METHOD!
    def get_top_bot_col_lst(self, del_col, cs_points_dist, cst_idx, smooth = True):
        #   run the get_top_bot_lst method above first
        #   cs_model.get_top_bot_lst_find_cst(self, cst_idx)
        #print ('Len top_elev: ' + str(len(self.top_elev)))
        #print ('Len bot_elev: ' + str(len(self.bot_elev)))
        out_top_lst, out_bot_lst, out_wtd_lst = [], [], []               #   define the output lists
        #out_glhymps_top_lay_lst, out_glhymps_bot_lay_lst = [], []
        ratio = int(cs_points_dist / del_col)            #   find out the ratio between the cross-section points distance and the del_col
        if smooth:                                      #   if the smoothin option is on assign to each column a specific elevation
            for x in range(int(round(ratio/2., 0))):   #   append the first topo values X times depending on the ratio
                out_top_lst.append(self.top_elev[0])
                out_bot_lst.append(self.bot_elev[0])
                out_wtd_lst.append(self.wtd_elev[0])
            for a in range(1, len(self.top_elev)):     #   now go on with the smoothing
            #for a in xrange(0, len(self.top_elev)):     #   now go on with the smoothing
                for b in range(ratio):
                    out_top_lst.append(round(self.top_elev[a-1] - b * (float(self.top_elev[a-1] - self.top_elev[a]) / (ratio)), 1))
                    out_bot_lst.append(round(self.bot_elev[a-1] - b * (float(self.bot_elev[a-1] - self.bot_elev[a]) / (ratio)), 1))
                    out_wtd_lst.append(round(self.wtd_elev[a-1] - b * (float(self.wtd_elev[a-1] - self.wtd_elev[a]) / (ratio)), 1))
            for x in range(int(round(ratio/2., 0))):
                out_top_lst.append(self.top_elev[a])
                out_bot_lst.append(self.bot_elev[a])
                out_wtd_lst.append(self.wtd_elev[a])
        #   if the smoothing is turned off give the coarse topography back loop through the lists and append the values to the output lists
        else:
            for a in range(len(self.top_elev)):
                for b in range(ratio):
                    out_top_lst.append(self.top_elev[a])
            for a in range(len(self.bot_elev)):
                for b in range(ratio):
                    out_bot_lst.append(self.bot_elev[a])
            for a in range(len(self.wtd_elev)):
                for b in range(ratio):
                    out_wtd_lst.append(self.wtd_elev[a])
        #   overwrite the top and bot elevation lists
        self.top_elev = out_top_lst
        self.bot_elev = out_bot_lst
        self.wtd_elev = out_wtd_lst
        #   return the final output lists
        #return out_top_lst, out_bot_lst

    #   do the same for other input lists
    def get_top_bot_col_lst_top_sys_geo(self, del_col, cs_points_dist, cst_idx, smooth = True):

        out_rch_pcr_lst, out_rch_watergap_lst, out_p_min_et_lst = [], [], []              
        out_hk_soil_lst, out_drn_density_lst = [], []         
        out_glhymps_top_lay_lst, out_glhymps_bot_lay_lst, out_soil_thk_lst = [], [], []            
        
        #ratio = int(cs_points_dist / del_col)            #   find out the ratio between the cross-section points distance and the del_col
        ratio = int(self.topo_res / del_col)            #   find out the ratio between the cross-section points distance and the del_col
        if smooth:                                      #   if the smoothin option is on assign to each column a specific elevation
            for x in range(int(round(ratio/2., 0))):   #   append the first topo values X times depending on the ratio
                out_rch_pcr_lst.append(self.rch_pcr[0])
                out_rch_watergap_lst.append(self.rch_watergap[0])
                out_p_min_et_lst.append(self.p_min_et[0])   
                out_hk_soil_lst.append(self.hk_soil[0])
                out_drn_density_lst.append(self.drn_density[0])  
                out_glhymps_top_lay_lst.append(self.glhymps_top_lay[0])
                out_glhymps_bot_lay_lst.append(self.glhymps_bot_lay[0])  
                out_soil_thk_lst.append(self.soil_thk[0])
                
            for a in range(1, len(self.rch_pcr)):     #   now go on with the smoothing
            #for a in xrange(0, len(self.top_elev)):     #   now go on with the smoothing
                for b in range(ratio):
                    out_rch_pcr_lst.append(self.rch_pcr[a-1] - b * (float(self.rch_pcr[a-1] - self.rch_pcr[a]) / (ratio)))
                    out_rch_watergap_lst.append(self.rch_watergap[a-1] - b * (float(self.rch_watergap[a-1] - self.rch_watergap[a]) / (ratio)))
                    out_p_min_et_lst.append(self.p_min_et[a-1] - b * (float(self.p_min_et[a-1] - self.p_min_et[a]) / (ratio)))
                    out_hk_soil_lst.append(self.hk_soil[a-1] - b * (float(self.hk_soil[a-1] - self.hk_soil[a]) / (ratio)))
                    out_drn_density_lst.append(self.drn_density[a-1] - b * (float(self.drn_density[a-1] - self.drn_density[a]) / (ratio)))
                    out_glhymps_top_lay_lst.append(self.glhymps_top_lay[a-1] - b * (float(self.glhymps_top_lay[a-1] - self.glhymps_top_lay[a]) / (ratio)))
                    out_glhymps_bot_lay_lst.append(self.glhymps_bot_lay[a-1] - b * (float(self.glhymps_bot_lay[a-1] - self.glhymps_bot_lay[a]) / (ratio)))                       
                    out_soil_thk_lst.append(self.soil_thk[a-1] - b * (float(self.soil_thk[a-1] - self.soil_thk[a]) / (ratio)))
                    
            for x in range(int(round(ratio/2., 0))):
                out_rch_pcr_lst.append(self.rch_pcr[a])
                out_rch_watergap_lst.append(self.rch_watergap[a])
                out_p_min_et_lst.append(self.p_min_et[a])   
                out_hk_soil_lst.append(self.hk_soil[a])
                out_drn_density_lst.append(self.drn_density[a])     
                out_glhymps_top_lay_lst.append(self.glhymps_top_lay[a])
                out_glhymps_bot_lay_lst.append(self.glhymps_bot_lay[a])   
                out_soil_thk_lst.append(self.soil_thk[a])
                
        #   if the smoothing is turned off give the coarse topography back loop through the lists and append the values to the output lists
        else:
            for a in range(len(self.rch_pcr)):
                for b in range(ratio):
                    out_rch_pcr_lst.append(self.rch_pcr[a])  
            for a in range(len(self.rch_watergap)):
                for b in range(ratio):
                    out_rch_watergap_lst.append(self.rch_watergap[a])
            for a in range(len(self.p_min_et)):
                for b in range(ratio):
                    out_p_min_et_lst.append(self.p_min_et[a])
            for a in range(len(self.hk_soil)):
                for b in range(ratio):
                    out_hk_soil_lst.append(self.hk_soil[a])
            for a in range(len(self.drn_density)):
                for b in range(ratio):
                    out_drn_density_lst.append(self.drn_density[a])                
            for a in range(len(self.glhymps_top_lay)):
                for b in range(ratio):
                    out_glhymps_top_lay_lst.append(self.glhymps_top_lay[a])                       
            for a in range(len(self.glhymps_bot_lay)):
                for b in range(ratio):
                    out_glhymps_bot_lay_lst.append(self.glhymps_bot_lay[a])   
            for a in range(len(self.soil_thk)):
                for b in range(ratio):
                    out_soil_thk_lst.append(self.soil_thk[a])   

        #   overwrite the top and bot elevation lists
        self.rch_pcr = out_rch_pcr_lst
        self.rch_watergap = out_rch_watergap_lst
        self.rch_p_min_et = out_p_min_et_lst
        self.hk_soil = out_hk_soil_lst
        self.drn_density = out_drn_density_lst
        self.glhymps_top_lay = out_glhymps_top_lay_lst
        self.glhymps_bot_lay = out_glhymps_bot_lay_lst
        self.soil_thk = out_soil_thk_lst
        #   return the final output lists
        #return out_top_lst, out_bot_lst

    #   Method that creates a list with bottom elevations of each cell in the model domain the script is made to work only
    #   for a constant del_lay value but can be updated for a list of del_lay values in case of grid refinement
    def modflow_bot_elev_list(self, del_lay, start_from_sea_level = True):
        self.delv = del_lay
        out_lay_elev_lst = []               #   create the output list

        #   start from the sea level
        if start_from_sea_level:
            #   first get the highest and lowest elevation from both lists
            top_elev = round(int(math.ceil(max(self.top_elev) / 10.) * 10 ), -1)
            bot_elev = int(math.ceil(min(self.bot_elev)))
            out_lay_elev_lst.append(top_elev)   #   append the top elevation overall to the list
            #   calculate the amount of layers necessary above and below sea level
            num_lay_above = int(math.ceil(max(self.top_elev) / del_lay))
            num_lay_below = abs(int(math.ceil(min(self.bot_elev) / del_lay)))
            num_lay = num_lay_above + num_lay_below
            #   create a loop to fill the layer elevation list
            for a in range(1, num_lay):
                lay_bot = out_lay_elev_lst[0] - a * del_lay
                out_lay_elev_lst.append(lay_bot)
            ##  add the bottom of the model domain
            if out_lay_elev_lst[-1] > bot_elev:
                out_lay_elev_lst.append(out_lay_elev_lst[0] - (a + 1) * del_lay)
            print('Total number of layers :     ' + str(len(out_lay_elev_lst) - 1))
            print('Bottom of the model is :     ' + str(out_lay_elev_lst[-1]))
            print('Layer distribution is uniform!')
            #   assign the final list
            self.lay_elev = out_lay_elev_lst

        else:
            #   first get the highest and lowest elevation from both lists
            top_elev = max(self.top_elev)
            bot_elev = min(self.bot_elev)
            out_lay_elev_lst.append(top_elev)   #   append the top elevation overall to the list
            #   get the total nr. of layers necessary based on the overall highest and lowest elevation
            num_lay = int(math.ceil(float(top_elev - bot_elev) / del_lay))
            #   create a loop to fill the layer elevation list
            for a in range(1, num_lay):
                lay_bot = out_lay_elev_lst[0] - a * del_lay
                out_lay_elev_lst.append(lay_bot)
            ##  add the bottom of the model domain
            if out_lay_elev_lst[-1] > bot_elev:
                out_lay_elev_lst.append(out_lay_elev_lst[0] - (a + 1) * del_lay)
            print('Total number of layers :     ' + str(len(out_lay_elev_lst) - 1))
            print('Bottom of the model is :     ' + str(out_lay_elev_lst[-1]))
            print('Layer distribution is uniform!')
            #   assign the final list
            self.lay_elev = out_lay_elev_lst
            #   return the final output list
            #return out_lay_elev_lst

    #   Method that creates a list with column discretizations (width of the columns)
    #   cs_points_dist = distance between the cross-section points in METERS!
    def modflow_col_list_constant(self, cs_points_dist, del_col):
        #   transform the indexes into distances
        #   this wasnt working properly and the coastline in the graphs was in the wrong positions
        #   therefore, first find the index of the first cell below the sea evel and adapt the
        #   x_start and x_end values.
        #self.x_coast_idx = self.top_elev.index([i for i in self.top_elev if i < 0][0])
        self.x_coast_idx = self.cst_idx
        #self.x_start = -1 * (self.x_coast_idx * del_col) / 1000. #   convert to km
        self.x_start = self.idx_start / (1000 / self.topo_res) - 200.
        #self.x_end = (self.idx_end * cs_points_dist) / (1000 / self.topo_res)  - 200.
        self.x_end = (self.idx_end) / (1000 / self.topo_res)  - 200.
        #self.x_end = (len(self.top_elev) * del_col + self.x_start * 1000) / 1000.
        #self.out_col_width_lst = int(math.ceil((self.x_end - self.x_start) * 1000 / del_col)) * [del_col] #   create the output list
        self.out_col_width_lst = len(self.top_elev) * [del_col] #   create the output list
        self.col_width = self.out_col_width_lst
        print(self.x_start, self.x_end)

        """
        self.x_coast_idx = self.top_elev.index([i for i in self.top_elev if i < 0][0])
        x_start = (self.x_coast_idx * cs_points_dist) / 1000. - 200. #   convert to km
        x_end = (self.idx_end * cs_points_dist) / 1000. - 200.
        out_col_width_lst = int(math.ceil((x_end - x_start) * 1000 / del_col)) * [del_col] #   create the output list
        self.x_start = x_start
        self.x_end = x_end
        self.col_width = out_col_width_lst
        print self.x_start, self.x_end
        """
        #   return the list
        #return out_col_width_lst, [x_start, x_end]

    #   Function that creates the IBOUND array for the MODFLOW model, the IBOUND values are:
    #       ibound < 0 - constant head cell (here = -1)
    #       ibound = 0 - inactive cell
    #       ibound > 0 - variable head (here = 1)
    #   The variable head cells are always only between the topography and sed. thickness estimation, the function
    #   loops through all the columns, and check the position of the middle of each cell and assigns IBOUND value
    def modflow_create_IBOUND_arr(self): #in_top_elev_lst, in_bot_elev_lst, in_col_width_lst, in_lay_bound):
        #   create the initial IBOUND array with the right dimensions
        out_ibound_arr = np.zeros((len(self.lay_elev) - 1, 1, len(self.col_width)), dtype=np.int32)
        dist_from_x0 = 0.0                  #   define starting point in X  direction
        #  loop through all the columns and assign the right IBOUND value to each cell in the column
        for c in range(len(self.col_width)):
            #  calculate the middle x coordinate of the column
            dist_from_x0 += float(self.col_width[c] / 1000.0)
            #   for the given column get the top and bottom elevation
            top_elev_col = self.top_elev[c]
            bot_elev_col = self.bot_elev[c]
            #   loop layer by layer and find its middle y value
            for d in range(len(self.lay_elev) - 1):
                #   find the middle y value
                mid_y = self.lay_elev[d] - ((self.lay_elev[d] - self.lay_elev[d + 1]) / 2.)
                #   check if the middle y value of the layer is between the top and bottom elevation for given column
                if mid_y < top_elev_col and mid_y > bot_elev_col:
                    out_ibound_arr[d, 0, c] = 1
        self.ibound_arr = out_ibound_arr
        #   return the final ibound array
        #return out_ibound_arr

    def modflow_create_IBOUND_arr_FOS(self, fos_x, fos_y): #in_top_elev_lst, in_bot_elev_lst, in_col_width_lst, in_lay_bound):
        #   create the initial IBOUND array with the right dimensions
        out_ibound_arr = np.zeros((len(self.lay_elev) - 1, 1, len(self.col_width)), dtype=np.int32)
        dist_from_x0 = 0.0                  #   define starting point in X  direction
        #  loop through all the columns and assign the right IBOUND value to each cell in the column
        for c in range(len(self.col_width)):
            #  calculate the middle x coordinate of the column
            dist_from_x0 += float(self.col_width[c] / 1000.0)
            #   for the given column get the top and bottom elevation
            top_elev_col = self.top_elev[c]
            bot_elev_col = self.bot_elev[c]
            #   loop layer by layer and find its middle y value
            for d in range(len(self.lay_elev) - 1):
                #   find the middle y value
                mid_y = self.lay_elev[d] - ((self.lay_elev[d] - self.lay_elev[d + 1]) / 2.)
                #   check if the middle y value of the layer is between the top and bottom elevation for given column
                if mid_y < top_elev_col and mid_y > bot_elev_col:
                    out_ibound_arr[d, 0, c] = 1
        self.ibound_arr = out_ibound_arr

        #   return the final ibound array
        #return out_ibound_arr


    #   Function that finds the topographical divide in the ibound array and defines a new IBOUND and dimensions
    #   the topo_diff is the elevation difference in meters
    def find_topo_divide(self, topo_diff, bot_limit, del_col, inland_cutoff_elev, fos_point, ibound_act_cells_limit, offshore_dist_limit):
        """
        
        topo_diff = 10.        
        bot_limit = -10000.
        del_col = del_col
        inland_cutoff_elev = model.top_elev[model.idx_start]
        fos_point = fos[1]
        ibound_act_cells_limit = 2
        offshore_dist_limit = 200000. 
        
        #   go through the top elevation list and find the divide
        for a in xrange(len(self.top_elev) - 1):
            if self.top_elev[a] > self.top_elev[a + 1] or self.top_elev[a] - self.top_elev[a + 1] > -topo_diff:
                pass
            else:
                break
        idx_fresh_bound = a + 1
        #   if there is no fresh water boundary found, set the starting index to 0 (initial position)
        if idx_fresh_bound == len(self.top_elev) - 1:
            idx_fresh_bound = 0
        divide_val = self.top_elev[idx_fresh_bound]
        #   also if the water divide is in the ocean make the col_start index 0
        if divide_val < 0.:
            col_idx = 0
        #   find the index of the divide value in the column list and the layer list
        else:
            col_idx = self.top_elev.index(divide_val)

        for c in xrange(1, len(self.lay_elev)):
            #   check where the new top value falls in the layer elevation list and trim it
            if self.lay_elev[c - 1] > divide_val > self.lay_elev[c]:
                lay_idx = c - 1
            else:
                lay_idx = 0
        """
    
        #   first check if there is a topographical divide or trim the values at a certain elevation                
        #cst_offset_val_2  = next((x for x in self.top_elev if x < 0.0), None)
        #cst_offset_idx_2 = self.top_elev.index(cst_offset_val_2) # - 1

        cst_offset_idx_2 = (self.cst_idx + 1 - self.idx_start) * int(self.topo_res / del_col)   
        #cst_offset_idx_2 = (self.cst_idx + 1 - self.idx_start) * 5                       
        #cst_offset_val_2 = self.top_elev[cst_offset_idx_2]    

        #if cst_offset_idx_2 < (self.cst_idx - self.idx_start) * 5:
        if cst_offset_idx_2 < (self.cst_idx - self.idx_start) *  int(self.topo_res / del_col):

            for a in range((self.cst_idx - self.idx_start) * 5, cst_offset_idx_2, -1):
                #print a, self.top_elev[a], self.top_elev[a - 1]
                if (self.top_elev[a] >= self.top_elev[a - 1] or self.top_elev[a - 1] - self.top_elev[a] < topo_diff):# and self.top_elev[a] < inland_cutoff_elev:
                    #print a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1]
                    idx_fresh_bound = 0                
                    pass
                else:
                    #print a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1]
                    idx_fresh_bound = a                            
                    break
            
        else:
            
            for a in range(cst_offset_idx_2, 0, -1):
                #print a
                if (self.top_elev[a - 1] >= self.top_elev[a] or self.top_elev[a] - self.top_elev[a - 1] < topo_diff):# and self.top_elev[a] < inland_cutoff_elev:
                    #print a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1]
                    idx_fresh_bound = 0                
                    pass
                else:
                    #print a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1]
                    idx_fresh_bound = a                            
                    break  
            
        #divide_val = self.top_elev[idx_fresh_bound]
        
        #   if there is no fresh water boundary found, set the starting index to 0 (initial position)
    
        if idx_fresh_bound == len(self.top_elev) - 1:
            idx_fresh_bound = 0
            
        divide_val = self.top_elev[idx_fresh_bound]
        
        #   also if the water divide is in the ocean make the col_start index 0
        if divide_val < 0.:
            col_idx = 0
        #   find the index of the divide value in the column list and the layer list
        else:
            col_idx = self.top_elev.index(divide_val)
        
        for c in range(1, len(self.lay_elev)):
            #   check where the new top value falls in the layer elevation list and trim it
            if self.lay_elev[c] > divide_val:
                self.lay_idx = 0
                pass
            else:
                self.lay_idx = c - 1
                break     
    
        #   adjust the lay_idx based on the highest elevation in the trimmed list
        max_elev = max(self.top_elev[col_idx:])
        for i in range(len(self.lay_elev)):
            if self.lay_elev[i] >= max_elev:
                self.lay_idx = i
            else:
                pass
                            
        #except UnboundLocalError:
        #    col_idx = 0
        #    self.lay_idx = 0
        #    print('fresh bound doesnt exist')
    
        #   go through the IBOUND array between the coast and the new inland boundary and check for columns with 
        #   a number of inactive cells lower than the limit value for considering them part of the self domain (2)
        for b in range(cst_offset_idx_2, col_idx, -1):
            #   get the column from the IBOUND array
            ibound_col = self.ibound_arr[:, 0, b].tolist()
            #   count the number of active cells in the column and check that there are more than the limit
            ibound_act_cells = [i for i in ibound_col if i != 0]
            #print ibound_act_cells
            if len(ibound_act_cells) > ibound_act_cells_limit:
                pass
            else:
                col_idx = b
                #print b
                break       

        #   now check the offshore part of the self domain, if there is no FOS point limit the distance from coast 
        #   to the input limit value. This is mainly to crop the huge offshore parts of the Delta cross-sections and 
        #   some others.
        if fos_point is None:
            #   calculate the limit column index offshore
            col_idx_off_limit = cst_offset_idx_2 + int(offshore_dist_limit / del_col)
            if col_idx_off_limit > len(self.bot_elev):
                col_idx_off_limit = len(self.bot_elev) - 1
            """
            for c in range(cst_offset_idx_2, len(self.top_elev), 1):
                #   get the column from the IBOUND array
                ibound_col = self.ibound_arr[:, 0, c].tolist()                
                #   count the number of active cells in the column and check that there are more than the limit
                ibound_act_cells = [i for i in ibound_col if i != 0]
                if len(ibound_act_cells) > ibound_act_cells_limit and c <= col_idx_off_limit:
                    pass
                else:
                    col_idx_off_limit = c
                    print c
                    break        
            """
            print (len(self.bot_elev), col_idx_off_limit)
            divide_bot_val = self.bot_elev[col_idx_off_limit - 1]
            for d in range(1, len(self.lay_elev)):
                #   check where the new top value falls in the layer elevation list and trim it
                if self.lay_elev[d] > divide_bot_val:
                    self.lay_end_idx = d
                    pass
                else:
                    self.lay_end_idx = d + 1
                    break      
    
            nlay_end_idx = self.lay_end_idx

        else:
            #   set the end index for columns to be where the top elevation is below -1000m
            try:
                col_idx_off_limit = self.top_elev.index([i for i in self.top_elev if i < bot_limit][0])
                #self.ncol = col_idx_off_limit
                #print '1a', col_idx_off_limit, self.ncol
            except IndexError:
                col_idx_off_limit = len(self.top_elev)
                print ('1b', col_idx_off_limit)
                
            try:
                nlay_end_val = min(self.bot_elev, key=lambda x:abs(x - ([i for i in self.top_elev if i < bot_limit][0])))
                nlay_end_idx = self.bot_elev.index(nlay_end_val)
                #nlay_end_val = divide_bot_val
                #nlay_end_idx = self.lay_end_idx # self.bot_elev.index(nlay_end_val)                
                self.nlay = nlay_end_idx
                self.bot_elev = self.bot_elev[: nlay_end_idx]
                print ('2a', nlay_end_val, nlay_end_idx, self.nlay, len(self.bot_elev))
            except IndexError:
                nlay_end_val = self.bot_elev[-1]
                nlay_end_idx = self.bot_elev.index(nlay_end_val)
                print ('2b', nlay_end_val, nlay_end_idx)   
            
            
        #   depending if the costal boundary is selected find the proper end index if not then just go on
        if not self.cst_bound:
            #   trim the ibound to new values
            #out_ibound_arr = self.ibound_arr[lay_idx : nlay_end_idx, :, col_idx : ncol_end_idx]
            print (col_idx, col_idx_off_limit)
            out_ibound_arr = self.ibound_arr[self.lay_idx : nlay_end_idx, :, col_idx : col_idx_off_limit]
            self.ibound_arr = out_ibound_arr
            self.lay_idx = self.lay_idx
            self.col_idx = col_idx
            self.end_idx = col_idx_off_limit
            
            self.top_elev = self.top_elev[self.col_idx : self.end_idx]            
            self.bot_elev = self.bot_elev[self.col_idx : self.end_idx]              
            #   self.idx_start and end are for trimming the other input data listes (e.g. recharge)            
            self.idx_start = self.idx_start + (self.col_idx / 5)
            self.idx_end = self.end_idx
            self.wtd_elev = self.wtd_elev[self.col_idx : self.end_idx]  
            self.ncol = self.ibound_arr.shape[-1]
            self.cst_idx = cst_offset_idx_2 - self.col_idx		
            self.x_start = 0 - (self.cst_idx * del_col / 1000.)
            #   update the x_end
            x_end2 = self.x_start + self.ncol * del_col / 1000.  # (self.end_idx * del_col + self.x_start * 1000) / 1000.
            #print self.end_idx, del_col, self.x_start, (self.end_idx * del_col + self.x_start * 1000) / 1000.
            #print self.x_end, x_end2
            #out_col_width_lst = int(math.ceil((x_end2 - self.x_start) * 1000 / del_col)) * [del_col] #   create the output list
            self.x_end = x_end2
            #self.col_width = out_col_width_lst
            self.col_width = len(self.top_elev) * [del_col] #out_col_width_lst
            print (self.x_start, self.x_end)             
            
            self.rch_pcr = self.rch_pcr[self.col_idx : self.end_idx]  
            self.rch_watergap = self.rch_watergap[self.col_idx : self.end_idx]  
            self.rch_p_min_et = self.rch_p_min_et[self.col_idx : self.end_idx]  
            self.hk_soil = self.hk_soil[self.col_idx : self.end_idx]  
            self.drn_density = self.drn_density[self.col_idx : self.end_idx]  
            self.glhymps_top_lay = self.glhymps_top_lay[self.col_idx : self.end_idx]  
            self.glhymps_bot_lay = self.glhymps_bot_lay[self.col_idx : self.end_idx]  
            self.soil_thk = self.soil_thk[self.col_idx : self.end_idx]  
            self.soil_type = self.soil_type[self.col_idx : self.end_idx]  
            



    def find_topo_divide_SRM_v2(self, topo_diff, bot_limit, del_col, del_lay, fos_point, max_elev, ibound_act_thk_limit_meters, offshore_dist_limit_min, offshore_dist_limit_max, inland_dist_limit_min, inland_dist_limit_max):
        """       
        topo_diff = 10.        
        bot_limit = -2000.
        del_col = del_col

        fos_point = fos[1]
        ibound_act_thk_limit_meters = 20
        del_lay = 10
        offshore_dist_limit_min = 50000. 
        offshore_dist_limit_max = 100000. 
        inland_dist_limit_min = 10000.
        inland_dist_limit_max = 50000.
        max_elev = 100.
        
        """
        
        """ FIRST DEAL WITH THE INLAND PART OF THE self DOMAIN """
        #   1) find the coast position from the topography list, the self.cst_idx is at this point the absolute coastline index 
        #      using the whole span of the 200km inland and 200km offshore cross-section. In this step we calculate the coastline
        #      index as column index relative to the start of the inland self boundary (col = 0)
        cst_idx_inl_bnd = int(self.cst_idx + 1 - self.idx_start) * int(self.topo_res / del_col)  
        inl_topo_lst = self.top_elev[: cst_idx_inl_bnd]
            
        #   2) check if the inland part of the self domain is longer than 10km, if not then do not trim the inland self domain.
        #      Also, if it stretches further than 50km then only take into account the interval between 10-50km from the coastline
        #      for finding the inland boundary - water divide. 
        inl_len_m = cst_idx_inl_bnd * self.topo_res
        if inl_len_m > inland_dist_limit_min:
            cst_bnd = cst_idx_inl_bnd - int(inland_dist_limit_min / self.topo_res)
        else:
            cst_bnd = cst_idx_inl_bnd
        if inl_len_m > inland_dist_limit_max:
            inl_bnd = cst_idx_inl_bnd - int(inland_dist_limit_max / self.topo_res)
        else:
            inl_bnd = 0
        inl_topo_lst = inl_topo_lst[inl_bnd : cst_bnd]

        #   3) Loop through the inland elevation, direction from coastline towards inland boundary. 
        for a in range(cst_bnd, inl_bnd, -1):
            #print a
            if (self.top_elev[a - 1] >= self.top_elev[a] or self.top_elev[a] - self.top_elev[a - 1] < topo_diff):# and self.top_elev[a] < inland_cutoff_elev:
                #print(a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1])
                idx_fresh_bound = inl_bnd                
                pass
            else:
                #print(a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1])
                idx_fresh_bound = a                            
                break          
        
        #   4) Check if there is an elevation above max_elev in the new inland cutout part, if so cut from the first cell with that elevation
        for b in range(cst_bnd, idx_fresh_bound, -1):
            if self.top_elev[b] <= max_elev:
                #print(b, self.top_elev[b])
                pass
            else:
                idx_fresh_bound = cst_bnd - b       
        #       define max elevation cell in the self domain, extend also into the first 10km inland just in case there is a higher elevated cell there
        max_elev_self = max(self.top_elev[idx_fresh_bound : cst_idx_inl_bnd])
        #       if there is no fresh water boundary found, set the starting index to 0 (initial position)
        if idx_fresh_bound == len(self.top_elev) - 1:
            idx_fresh_bound = 0

        
        #   5) Trim the layer elevation list so it covers only the new active self domain
        for c in range(1, len(self.lay_elev)):
        #       check where the new top value falls in the layer elevation list and trim it
            if self.lay_elev[c] > max_elev_self:
                self.lay_idx = 0
                pass
            else:
                self.lay_idx = c - 1
                break     
        #       adjust the lay_idx based on the highest elevation in the trimmed list
        for i in range(len(self.lay_elev)):
            if self.lay_elev[i] >= max_elev_self:
                self.lay_idx = i
            else:
                pass        
    
        #   6) Go through the IBOUND array between the coast and the new inland boundary and check for columns with 
        #      a number of inactive cells lower than the limit value for considering them part of the self domain (2)
        ibound_act_cells_limit = int(ibound_act_thk_limit_meters / del_lay)
        for b in range(cst_idx_inl_bnd, idx_fresh_bound, -1):
            #   get the column from the IBOUND array
            ibound_col = self.ibound_arr[:, 0, b].tolist()
            #   count the number of active cells in the column and check that there are more than the limit
            ibound_act_cells = [i for i in ibound_col if i != 0]
            #print(ibound_act_cells)
            if len(ibound_act_cells) >= ibound_act_cells_limit:
                pass
            else:
                idx_fresh_bound = b
                #print(b)
                break
        col_idx = idx_fresh_bound 

        """ NEXT, DEAL WITH THE OFFSHORE PART OF THE self DOMAIN """
        #   7) if there is no FOS point limit the distance from coast to the input limit value. This is mainly to crop the huge offshore parts 
        #      of the Delta cross-sections and other such cases.

        #   calculate the limit column index offshore
        col_idx_off_limit = cst_idx_inl_bnd + int(offshore_dist_limit_max / del_col)
        if col_idx_off_limit > len(self.bot_elev):
            col_idx_off_limit = len(self.bot_elev) - 1
        print (len(self.bot_elev), col_idx_off_limit)
        divide_bot_val = self.bot_elev[col_idx_off_limit - 1]
        for d in range(1, len(self.lay_elev)):
            #   check where the new top value falls in the layer elevation list and trim it
            if self.lay_elev[d] > divide_bot_val:
                self.lay_end_idx = d
                pass
            else:
                self.lay_end_idx = d + 1
                break      
        nlay_end_idx = self.lay_end_idx

        """
        if fos_point is None:
            #   calculate the limit column index offshore
            col_idx_off_limit = cst_idx_inl_bnd + int(offshore_dist_limit_max / del_col)
            if col_idx_off_limit > len(self.bot_elev):
                col_idx_off_limit = len(self.bot_elev) - 1
            print (len(self.bot_elev), col_idx_off_limit)
            divide_bot_val = self.bot_elev[col_idx_off_limit - 1]
            for d in range(1, len(self.lay_elev)):
                #   check where the new top value falls in the layer elevation list and trim it
                if self.lay_elev[d] > divide_bot_val:
                    self.lay_end_idx = d
                    pass
                else:
                    self.lay_end_idx = d + 1
                    break      
            nlay_end_idx = self.lay_end_idx

        else:
            #   set the end index for columns to be where the top elevation is below -1000m
            try:
                col_idx_off_limit = self.top_elev.index([i for i in self.top_elev if i < bot_limit][0])
            except IndexError:
                col_idx_off_limit = len(self.top_elev)
                print ('1b', col_idx_off_limit)
        
            try:
                nlay_end_val = min(self.bot_elev, key=lambda x:abs(x - ([i for i in self.top_elev if i < bot_limit][0])))
                nlay_end_idx = self.bot_elev.index(nlay_end_val)     
                self.nlay = nlay_end_idx
                self.bot_elev = self.bot_elev[: nlay_end_idx]
                print ('2a', nlay_end_val, nlay_end_idx, self.nlay, len(self.bot_elev))
            except IndexError:
                nlay_end_val = self.bot_elev[-1]
                nlay_end_idx = -1 #self.bot_elev.index(nlay_end_val)
                print ('2b', nlay_end_val, nlay_end_idx)   
        """

        #   depending if the costal boundary is selected find the proper end index if not then just go on
        if not self.cst_bound:
            #   trim the ibound to new values
            #out_ibound_arr = self.ibound_arr[lay_idx : nlay_end_idx, :, col_idx : ncol_end_idx]
            print (col_idx, col_idx_off_limit)
            out_ibound_arr = self.ibound_arr[self.lay_idx : -1, :, col_idx : col_idx_off_limit]
            self.ibound_arr = out_ibound_arr
            self.lay_idx = self.lay_idx
            self.col_idx = col_idx
            self.end_idx = col_idx_off_limit                
                
            self.top_elev = self.top_elev[self.col_idx : self.end_idx]            
            self.bot_elev = self.bot_elev[self.col_idx : self.end_idx]              
            #   self.idx_start and end are for trimming the other input data listes (e.g. recharge)            
            self.idx_start = self.col_idx #self.idx_start + (self.col_idx / 5)
            self.idx_end = self.end_idx
            self.wtd_elev = self.wtd_elev[self.col_idx : self.end_idx]  
            self.ncol = self.ibound_arr.shape[-1]
            self.cst_idx = cst_idx_inl_bnd - self.col_idx		
            self.x_start = 0 - (self.cst_idx * del_col / 1000.)
            #   update the x_end
            x_end2 = self.x_start + self.ncol * del_col / 1000.  # (self.end_idx * del_col + self.x_start * 1000) / 1000.
            #print self.end_idx, del_col, self.x_start, (self.end_idx * del_col + self.x_start * 1000) / 1000.
            #print self.x_end, x_end2
            #out_col_width_lst = int(math.ceil((x_end2 - self.x_start) * 1000 / del_col)) * [del_col] #   create the output list
            self.x_end = x_end2
            #self.col_width = out_col_width_lst
            self.col_width = len(self.top_elev) * [del_col] #out_col_width_lst
            print (self.x_start, self.x_end)             
            
            self.rch_pcr = self.rch_pcr[self.col_idx : self.end_idx]  
            self.rch_watergap = self.rch_watergap[self.col_idx : self.end_idx]  
            self.rch_p_min_et = self.rch_p_min_et[self.col_idx : self.end_idx]  
            self.hk_soil = self.hk_soil[self.col_idx : self.end_idx]  
            self.drn_density = self.drn_density[self.col_idx : self.end_idx]  
            self.glhymps_top_lay = self.glhymps_top_lay[self.col_idx : self.end_idx]  
            self.glhymps_bot_lay = self.glhymps_bot_lay[self.col_idx : self.end_idx]  
            self.soil_thk = self.soil_thk[self.col_idx : self.end_idx]  
            self.soil_type = self.soil_type[self.col_idx : self.end_idx]     




    def find_topo_divide_SRM(self, topo_diff, bot_limit, del_col, del_lay, fos_point, ibound_act_thk_limit_meters, offshore_dist_limit_min, offshore_dist_limit_max, inland_dist_limit_min, inland_dist_limit_max):
        """       
        topo_diff = 10.        
        bot_limit = -2000.
        del_col = del_col
        inland_cutoff_elev = model.top_elev[model.idx_start]
        fos_point = fos[1]
        ibound_act_thk_limit_meters = 20
        del_lay = 10
        offshore_dist_limit_min = 50000. 
        offshore_dist_limit_max = 100000. 
        inland_dist_limit_min = 10000.
        inland_dist_limit_max = 50000.
        """

        #   trim the topo elevation list to only include the inland elevations from 50km to 10km dist from coast
        #inland_topo_lst = model.top_elev[max(model.idx_start, model.cst_idx - int(inland_dist_limit_max / del_col)) : model.cst_idx - int(inland_dist_limit_min / del_col)]



        #   first check if there is a topographical divide or trim the values at a certain elevation                
        cst_offset_idx_2 = (self.cst_idx + 1 - self.idx_start) * int(self.topo_res / del_col)   

        if cst_offset_idx_2 < (self.cst_idx - self.idx_start) *  int(self.topo_res / del_col):

            for a in range((self.cst_idx - self.idx_start) * 5, cst_offset_idx_2, -1):
                #print a, self.top_elev[a], self.top_elev[a - 1]
                if (self.top_elev[a] >= self.top_elev[a - 1] or self.top_elev[a - 1] - self.top_elev[a] < topo_diff):# and self.top_elev[a] < inland_cutoff_elev:
                    #print a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1]
                    idx_fresh_bound = 0                
                    pass
                else:
                    #print a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1]
                    idx_fresh_bound = a                            
                    break
            
        else:
            
            for a in range(cst_offset_idx_2, 0, -1):
                #print a
                if (self.top_elev[a - 1] >= self.top_elev[a] or self.top_elev[a] - self.top_elev[a - 1] < topo_diff):# and self.top_elev[a] < inland_cutoff_elev:
                    #print a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1]
                    idx_fresh_bound = 0                
                    pass
                else:
                    #print a, self.top_elev[a], '----', a-1, self.top_elev[a - 1], '----', self.top_elev[a] - self.top_elev[a - 1]
                    idx_fresh_bound = a                            
                    break  
            
        #divide_val = self.top_elev[idx_fresh_bound]
        
        #   if there is no fresh water boundary found, set the starting index to 0 (initial position)
    
        if idx_fresh_bound == len(self.top_elev) - 1:
            idx_fresh_bound = 0
            
        divide_val = self.top_elev[idx_fresh_bound]
        
        #   also if the water divide is in the ocean make the col_start index 0
        if divide_val < 0.:
            col_idx = 0
        #   find the index of the divide value in the column list and the layer list
        else:
            col_idx = self.top_elev.index(divide_val)
        
        for c in range(1, len(self.lay_elev)):
            #   check where the new top value falls in the layer elevation list and trim it
            if self.lay_elev[c] > divide_val:
                self.lay_idx = 0
                pass
            else:
                self.lay_idx = c - 1
                break     
    
        #   adjust the lay_idx based on the highest elevation in the trimmed list
        max_elev = max(self.top_elev[col_idx:])
        for i in range(len(self.lay_elev)):
            if self.lay_elev[i] >= max_elev:
                self.lay_idx = i
            else:
                pass
                            
        #except UnboundLocalError:
        #    col_idx = 0
        #    self.lay_idx = 0
        #    print('fresh bound doesnt exist')
    
        #   go through the IBOUND array between the coast and the new inland boundary and check for columns with 
        #   a number of inactive cells lower than the limit value for considering them part of the self domain (2)
        ibound_act_cells_limit = int(ibound_act_thk_limit_meters / del_lay)

        """   NOT SURE ? """
        cst_offset_idx_2 = max(cst_offset_idx_2 - 100, 100)
        
        for b in range(cst_offset_idx_2, col_idx, -1):
            #   get the column from the IBOUND array
            ibound_col = self.ibound_arr[:, 0, b].tolist()
            #   count the number of active cells in the column and check that there are more than the limit
            ibound_act_cells = [i for i in ibound_col if i != 0]
            #print ibound_act_cells
            if len(ibound_act_cells) >= ibound_act_cells_limit:
                pass
            else:
                col_idx = b
                #print b
                break       

        #   now check the offshore part of the self domain, if there is no FOS point limit the distance from coast 
        #   to the input limit value. This is mainly to crop the huge offshore parts of the Delta cross-sections and 
        #   some others.
        if fos_point is None:
            #   calculate the limit column index offshore
            col_idx_off_limit = cst_offset_idx_2 + int(offshore_dist_limit_max / del_col)
            if col_idx_off_limit > len(self.bot_elev):
                col_idx_off_limit = len(self.bot_elev) - 1
            """
            for c in range(cst_offset_idx_2, len(self.top_elev), 1):
                #   get the column from the IBOUND array
                ibound_col = self.ibound_arr[:, 0, c].tolist()                
                #   count the number of active cells in the column and check that there are more than the limit
                ibound_act_cells = [i for i in ibound_col if i != 0]
                if len(ibound_act_cells) > ibound_act_cells_limit and c <= col_idx_off_limit:
                    pass
                else:
                    col_idx_off_limit = c
                    print c
                    break        
            """
            print (len(self.bot_elev), col_idx_off_limit)
            divide_bot_val = self.bot_elev[col_idx_off_limit - 1]
            for d in range(1, len(self.lay_elev)):
                #   check where the new top value falls in the layer elevation list and trim it
                if self.lay_elev[d] > divide_bot_val:
                    self.lay_end_idx = d
                    pass
                else:
                    self.lay_end_idx = d + 1
                    break      
    
            nlay_end_idx = self.lay_end_idx

        else:
            #   set the end index for columns to be where the top elevation is below -1000m
            try:
                col_idx_off_limit = self.top_elev.index([i for i in self.top_elev if i < bot_limit][0])
                #self.ncol = col_idx_off_limit
                #print '1a', col_idx_off_limit, self.ncol
            except IndexError:
                col_idx_off_limit = len(self.top_elev)
                print ('1b', col_idx_off_limit)
                
            try:
                nlay_end_val = min(self.bot_elev, key=lambda x:abs(x - ([i for i in self.top_elev if i < bot_limit][0])))
                nlay_end_idx = self.bot_elev.index(nlay_end_val)
                #nlay_end_val = divide_bot_val
                #nlay_end_idx = self.lay_end_idx # self.bot_elev.index(nlay_end_val)                
                self.nlay = nlay_end_idx
                self.bot_elev = self.bot_elev[: nlay_end_idx]
                print ('2a', nlay_end_val, nlay_end_idx, self.nlay, len(self.bot_elev))
            except IndexError:
                nlay_end_val = self.bot_elev[-1]
                nlay_end_idx = self.bot_elev.index(nlay_end_val)
                print ('2b', nlay_end_val, nlay_end_idx)   
            
        #   depending if the costal boundary is selected find the proper end index if not then just go on
        if not self.cst_bound:
            #   trim the ibound to new values
            #out_ibound_arr = self.ibound_arr[lay_idx : nlay_end_idx, :, col_idx : ncol_end_idx]
            print (col_idx, col_idx_off_limit)
            out_ibound_arr = self.ibound_arr[self.lay_idx : nlay_end_idx, :, col_idx : col_idx_off_limit]
            self.ibound_arr = out_ibound_arr
            self.lay_idx = self.lay_idx
            self.col_idx = col_idx
            self.end_idx = col_idx_off_limit
            
            self.top_elev = self.top_elev[self.col_idx : self.end_idx]            
            self.bot_elev = self.bot_elev[self.col_idx : self.end_idx]              
            #   self.idx_start and end are for trimming the other input data listes (e.g. recharge)            
            self.idx_start = self.idx_start + (self.col_idx / 5)
            self.idx_end = self.end_idx
            self.wtd_elev = self.wtd_elev[self.col_idx : self.end_idx]  
            self.ncol = self.ibound_arr.shape[-1]
            self.cst_idx = cst_offset_idx_2 - self.col_idx		
            self.x_start = 0 - (self.cst_idx * del_col / 1000.)
            #   update the x_end
            x_end2 = self.x_start + self.ncol * del_col / 1000.  # (self.end_idx * del_col + self.x_start * 1000) / 1000.
            #print self.end_idx, del_col, self.x_start, (self.end_idx * del_col + self.x_start * 1000) / 1000.
            #print self.x_end, x_end2
            #out_col_width_lst = int(math.ceil((x_end2 - self.x_start) * 1000 / del_col)) * [del_col] #   create the output list
            self.x_end = x_end2
            #self.col_width = out_col_width_lst
            self.col_width = len(self.top_elev) * [del_col] #out_col_width_lst
            print (self.x_start, self.x_end)             
            
            self.rch_pcr = self.rch_pcr[self.col_idx : self.end_idx]  
            self.rch_watergap = self.rch_watergap[self.col_idx : self.end_idx]  
            self.rch_p_min_et = self.rch_p_min_et[self.col_idx : self.end_idx]  
            self.hk_soil = self.hk_soil[self.col_idx : self.end_idx]  
            self.drn_density = self.drn_density[self.col_idx : self.end_idx]  
            self.glhymps_top_lay = self.glhymps_top_lay[self.col_idx : self.end_idx]  
            self.glhymps_bot_lay = self.glhymps_bot_lay[self.col_idx : self.end_idx]  
            self.soil_thk = self.soil_thk[self.col_idx : self.end_idx]  
            self.soil_type = self.soil_type[self.col_idx : self.end_idx]  


    """
    #   check that there are no disconnected areas in the IBOUND array -
    def check_ibound_arr(self):
        self.new_start_idx = 0
        for j in xrange(len(self.cst_idx):
            #   check if there are any active cells in the column, if yes then move to the next one, if not then remember that column number as a starting column
            active_cells = [i for i in self.ibound_arr[:, 0, i].tolist() if i == 1]
            #   check the 1st column with ibound_val = 1 (active cell)
            if active_cells != []:
                pass
            else:
                self.new_start_idx = j + 1
    """            
    
    #   Method that trims the IBOUND array based on the specified maximum depth of the cross-section
    def dis_input(self, nrow, delc, del_col, del_lay, nper, perlen, nstp, laycbd, max_depth, cs_points_dist):
        #   assign all the necessary information
        self.nlay = len(self.lay_elev[self.lay_idx:]) - 1                   #   number of layers in model domain / -1 to remove the top elevation in the lay_elev_lst
        self.nrow = nrow                                                    #   number of rows in model domain (2D so = 1)
        #self.ncol = len(self.col_width[self.col_idx : self.end_idx])        #   number of column in model domain
        self.nper = nper                                                    #   number of model stress periods
        #self.delr = self.col_width[self.col_idx :]# self.end_idx]             #   an array of spacings along a row, a bit confusing - check the MODFLOW manual in case..
        #self.delr = self.col_width # self.end_idx]        
        self.delr = len(self.top_elev) * [del_col]        
        self.delc = delc                                                    #   an array of spacings along a column, using the list created earlier (2D so = 1)
        
        #self.top = round(np.nanmax(self.top_elev) / 10., 0) * 10
        self.top = self.lay_elev[0]
        
        #try:
        #    self.top = math.ceil((self.top / 10.0) * 10.0)
        #except AttributeError:
        #    self.top = round(np.nanmax(self.top_elev) / 10., 0) * 10

        #y_max_new = self.top
        self.bottom = round(np.nanmin(self.bot_elev) / 10., 0) * 10
        
        try:
            lay_top_idx = self.lay_elev.index(self.top)
        except ValueError:
            self.lay_elev = np.linspace(self.top, self.lay_elev[0], int((self.top - self.lay_elev[0]) / self.del_lay)).tolist() + self.lay_elev
            lay_top_idx = 0
        
        try:
            lay_bot_idx = self.lay_elev.index(self.bottom)
            self.lay_elev = self.lay_elev[lay_top_idx : lay_bot_idx + 1]
        #   happens when writing the package for the second time in the script..
        except ValueError:
            self.lay_elev = self.lay_elev[:]
            
        #self.top = self.lay_elev[self.lay_idx]                              #   top elevation of the model domain - check if it should be like this!
        #self.botm = self.lay_elev[self.lay_idx + 1:]                        #   an array of the bottom elevation for each model cell
        self.top = self.lay_elev[0]                              #   top elevation of the model domain - check if it should be like this!
        self.botm = self.lay_elev[1:]                        #   an array of the bottom elevation for each model cell
        
        self.perlen = perlen                                                #   an array of the stress period lengths
        self.nstp = nstp                                                    #   number of time steps in each stress period
        self.laycbd = laycbd                                                #   0 indicates no confining bed... -> this needs to be checked!!
        self.zbot = self.top - self.botm[-1]                                #   total depth of the model
        self.del_lay = del_lay
        self.del_col = del_col
        self.ncol = self.ibound_arr.shape[-1]
        
        #   set the end index for columns to be where the top elevation is below the max_depth
        try:
            self.idx_end = self.bot_elev.index([i for i in self.bot_elev if i < max_depth][0])
            self.ncol = self.idx_end
            self.delr = self.delr[: self.idx_end]
        except IndexError:
            self.idx_end = self.ncol
        try:
            self.nlay_end_val = min(self.lay_elev, key = lambda x : abs(x - ([i for i in self.lay_elev if i < max_depth][0])))
            self.nlay_end_idx = self.lay_elev.index(self.nlay_end_val)
            self.nlay = self.nlay_end_idx
            self.botm = self.botm[: self.nlay_end_idx]
        except IndexError:
            self.nlay_end_val = self.botm[-1]
            self.nlay_end_idx = self.lay_elev.index(self.nlay_end_val)
        
        #   set the end index for columns to be where the top elevation is below the max_depth
        try:
            self.end_idx = self.bot_elev.index([i for i in self.bot_elev if i < max_depth][0])
            #self.ncol = self.end_idx
            #self.delr = self.delr[: self.end_idx]
        except IndexError:
            self.end_idx = self.ncol - self.col_idx
        try:
            self.nlay_end_val = min(self.lay_elev, key = lambda x : abs(x - ([i for i in self.lay_elev if i < max_depth][0])))
            self.nlay_end_idx = self.lay_elev.index(self.nlay_end_val)
            self.nlay = self.nlay_end_idx
            self.botm = self.botm[: self.nlay_end_idx]
        except IndexError:
            self.nlay_end_val = self.botm[-1]
            self.nlay_end_idx = self.lay_elev.index(self.nlay_end_val)

        #   now trim the IBOUND array accordingly and change all necessary variables
        self.zbot = self.botm[-1]
        #self.ibound_arr = self.ibound_arr[: nlay_end_idx, :, 0 : self.end_idx]
        #self.col_width = cs_model.modflow_col_list_constant(self, cs_points_dist, del_col)   #   500. is the cross-section point distance
        self.x0, self.x1 = 0, self.ibound_arr.shape[-1] * del_col

    #  Define the SCONC and STRT arrays (starting concentration, starting head)
    def bas_input(self, sconc_arr, strt_arr, sea_level, fresh_salt = None):
        if sconc_arr is None and strt_arr is None:# or fresh_salt is None:
            print('no sconc_arr and strt_arr')
            if fresh_salt == 'salt':
                #self.strt_arr = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.float64)
                #self.sconc_arr = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.float64) * 35.
                self.sconc_arr = self.ibound_arr * 35.
                self.strt_arr = self.ibound_arr * 1.
            elif fresh_salt == 'fresh':
                #self.strt_arr = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.float64)
                #self.sconc_arr = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.float64) * 0.
                self.sconc_arr = self.ibound_arr * 0.
                self.strt_arr = self.ibound_arr * 1.
            else:
                #self.strt_arr = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.float64)
                #self.sconc_arr = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.float64) * 35.                          
                self.sconc_arr = self.ibound_arr * 0.
                self.strt_arr = self.ibound_arr * sea_level
                #   find the indexes of columns where the elevation is lower than sea level, assign the concentration  of 35.0 to these cells in the sconc_array
                for a in range(self.ibound_arr.shape[2]):
                    if self.top_elev[a] < 0:
                        self.sconc_arr[:, 0, a] = 35.0
        else:
            #self.strt_arr = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.float64)
            #self.sconc_arr = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.float64) * 35.                 
            self.sconc_arr = sconc_arr
            self.strt_arr = strt_arr
        """
        elif fresh_salt == 'salt':
            if sconc_arr is None:
                self.sconc_arr = self.ibound_arr * 35.
                self.strt_arr = self.ibound_arr * 1.
            else:
                self.sconc_arr = sconc_arr
                self.strt_arr = strt_arr
        elif fresh_salt == 'fresh':
            if sconc_arr is None:
                self.sconc_arr = self.ibound_arr * 0.
                self.strt_arr = self.ibound_arr * 1.
            else:
                self.sconc_arr = sconc_arr
                self.strt_arr = strt_arr
        """



    #   Define the input for the LPF package - homogenous geology
    def lpf_input_simple_geo(self, Hk_val, Vk_val, laytyp_val):
        self.hk_arr = Hk_val * self.ibound_arr
        self.vk_arr = Vk_val * self.ibound_arr
        self.laytyp = laytyp_val


    #   Define the input for the LPF package - heterogeneous geology
    #   rand_seed   : value that ensures the reproducability of the randomization process
    #   Hk_val      : horizontal conductivity values
    #   Vk_val      : vertical conductivity values
    #   laytyp_val  : MODFLOW layer type
    #   sand_pct    : sand percentage in the sedimental body
    #   clay_pct    : clay percentage in the sedimental body
    #   n_aqt_in    : number of aquitard layers inland
    #   aqt_in_x0   : the inland limit of the extent of aquitard layers (in -km from coastline)
    #   n_aqt_off   : number of aquitard layers offshore
    #   aqt_off_x1  : the offshore limit of the extent of aquitard layers (in km from coastline)
    def lpf_input_complex_geo(self, rand_seed, Hk_val, Vk_val, laytyp_val, sand_pct, clay_pct, n_aqt_in, aqt_in_x1, n_aqt_off, aqt_off_x0,
                              top_soil = True, top_offshore = True, rand_aqt_inland = True, rand_aqt_offshore = True):

        def split(a, n):
            k, m = divmod(len(a), n)
            return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

        #   set the random.seed
        random.seed(a = rand_seed)

        #   create an array of hydraulic conductivity with all values = aquifer Kh
        self.hk_arr = Hk_val * self.ibound_arr
        self.vk_arr = Vk_val * self.ibound_arr
        self.laytyp = laytyp_val

        #Hk_val = 10.
        #Vk_val = 1.
        #laytyp_val = 0

        #   find out how thick is the soil layer and what type of soil it is
        self.soil_thk = self.soil_thk[self.idx_start : self.idx_end]
        self.soil_type = self.soil_type[self.idx_start : self.idx_end]

        self.glhymps_top_lay = self.glhymps_top_lay[self.idx_start : self.idx_end]
        self.glhymps_bot_lay = self.glhymps_bot_lay[self.idx_start : self.idx_end]




        soil_type_lst = [i for i in self.soil_type if i > -1]
        soil_thk_lst = [j for j in self.soil_thk if j > -1]

        #   repeat the elements of each list n times depending on the col width
        soil_thk_lst = np.repeat(soil_thk_lst, int(500. / self.del_col)).tolist()
        soil_type_lst = np.repeat(soil_type_lst, int(500. / self.del_col)).tolist()

        start_idx = int(abs(self.x_start * 1000. / self.del_col))

        #   create the top soil layer based on soilgrids
        if top_soil:
            #   loop through all the self columns
            for a in range(len(soil_type_lst)):
                #   select the column from the ibound array
                col_a = self.ibound_arr[:, :, a].tolist()
                indices = [i for i, x in enumerate(col_a) if x == [1]]
                #   find how many layers will have the soil K value based on the thickness
                try:
                    idx_count = int(round(soil_thk_lst[a] / (100 * self.del_lay))) # * 100 because the thickness is in cm
                except IndexError:  #   in case the list of soil types is longer than soil thicknesses, take the last thickness
                    idx_count = int(round(soil_thk_lst[-1] / (100 * self.del_lay)))
                try:
                    #   for the number of indices change the K values in the arrays
                    for b in range(idx_count):
                        self.hk_arr[indices[b], :, a] = k_soil_dict[soil_type_lst[a]]['k_val']
                        self.vk_arr[indices[b], :, a] = k_soil_dict[soil_type_lst[a]]['k_val'] / 10.  # this has to be changed later on!!!
                #   in case no cell is active in the given column of the ibound array
                except IndexError:
                    continue

        #   create a clay (mud) layer on top of the continental shelf
        if top_offshore:
            #   loop through all the self columns
            for c in range(a, self.ncol):
                #   select the column from the ibound array
                col_c = self.ibound_arr[:, :, c].tolist()
                indices = [i for i, x in enumerate(col_c) if x == [1]]
                #   take the last thickness
                idx_count = int(round(soil_thk_lst[-1] / (100 * self.del_lay)))
                try:
                    #   for the number of indices change the K values in the arrays
                    if self.top_elev[c] > -120.:
                        for b in range(idx_count):
                            self.hk_arr[indices[b], :, c] = k_soil_dict[1]['k_val']
                            self.vk_arr[indices[b], :, c] = k_soil_dict[1]['k_val'] / 10.  # this has to be changed later on!!!
                #   in case no cell is active in the given column of the ibound array
                except IndexError:
                    continue

        #   create a random amount of layers in the inland part of the self domain
        if rand_aqt_inland:
            aqt_in_x0 = self.x_start + 5.1
            #   define the domain where the aquitards occure - between the soil layer and the bottom
            #   in the list below the range of layers will be stored - based on the sand/clay ratio
            aqt_bot_top_lst = []
            #for i in range(0, self.x_coast_idx):
            for i in range(0, start_idx):
                #   find the number of layers with the aquifer Hk - this will make sure aquitard is in the aquifer and not the soil layer
                col_i = self.hk_arr[:, :, i].tolist()
                indices = [a for a, x in enumerate(col_i) if x == [Hk_val]]
                #   append to the bot_top_lst the slices of the indices list, depending on total number of aquitard layers
                inner_lst = list(split(indices, n_aqt_in))
                """
                inner_lst = [i]
                idx_top_aqt = 0
                for j in xrange(n_aqt_in):
                    if j != n_aqt_in - 1:
                        inner_lst.append(indices[idx_top_aqt : len(indices) / n_aqt_in + idx_top_aqt + 1])
                        idx_top_aqt += len(indices) / n_aqt_in + 1
                    else:
                        inner_lst.append(indices[idx_top_aqt :])
                """
                aqt_bot_top_lst.append(inner_lst)
            #   loop through the number of aquitard layers inland
            for j in range(n_aqt_in):
                #   get random extent in the horizontal direction
                aqt_extent_1 = round(random.uniform(aqt_in_x0, self.x_start/2), 1)
                aqt_extent_2 = round(random.uniform(self.x_start/4, aqt_in_x1), 1)
                #   find the column indexes of the start and end of the aquitard layer
                aqt_start = min(aqt_extent_1, aqt_extent_2)
                aqt_end = max(aqt_extent_1, aqt_extent_2)
                aqt_col_start = int(abs((self.x_start - aqt_start) * 1000. / self.del_col))
                aqt_col_end = int(abs((self.x_start - aqt_end) * 1000. / self.del_col))

                random_layer = False
                linear_layer = True

                #   randomly distributed layer in vertical direction
                if random_layer:
                    #   loop through the aquitard layer column extent and change the hydraylic conductivity to clay in corresponding cells
                    for g in range(aqt_col_start, aqt_col_end):
                        aqt_lay_indices = aqt_bot_top_lst[g][j]
                        #   split the indices based on the sand/clay ratio - to make sure the full thickness is possible. If it is 50/50
                        #   then the upper limit of the aquitard layer should be in the upper half of the indices so full thickness is possible
                        clay_ratio = 100 / clay_pct
                        aqt_lay_start = random.choice(aqt_lay_indices[: (-1 * (len(aqt_lay_indices) / clay_ratio))])
                        aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                        #   change the values in the array
                        self.hk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val']
                        self.vk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val'] / 10.  # this has to be changed later on!!!
                #   linearly distributed layer in vertical direction
                elif linear_layer:
                    #   select the column in the middle
                    aqt_lay_indices = aqt_bot_top_lst[len(aqt_bot_top_lst) / 2][j]

                    if aqt_lay_indices == []:
                        for g in range(len(aqt_bot_top_lst) / 2, len(aqt_bot_top_lst)):
                            not_empty = 0
                            for h in range(n_aqt_in):
                                if aqt_bot_top_lst[g][h] != []:
                                    not_empty += 1
                                else:
                                    pass
                            if not_empty == n_aqt_in:
                                aqt_lay_indices = aqt_bot_top_lst[g][j]
                                break
                    try:
                        #   split the indices based on the sand/clay ratio - to make sure the full thickness is possible. If it is 50/50
                        #   then the upper limit of the aquitard layer should be in the upper half of the indices so full thickness is possible
                        clay_ratio = 100 / clay_pct
                        aqt_lay_start = random.choice(aqt_lay_indices[: (-1 * (len(aqt_lay_indices) / clay_ratio))])
                        aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])

                    except IndexError:
                        aqt_lay_start = aqt_lay_indices[0]
                        aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                        #   express the position of the first layer as distance from the upper limit
                        start_ratio = int(round(aqt_lay_indices.index(aqt_lay_start) / float(len(aqt_lay_indices)) * 100))
                        #   loop through the aquitard layer column extent and change the hydraylic conductivity to clay in corresponding cells
                        for g in range(aqt_col_start, aqt_col_end):
                            try:
                                aqt_lay_indices = aqt_bot_top_lst[g][j]
                                #   find the starting index based on the ratio
                                try:
                                    start_idx_g = int(round(start_ratio / float(len(aqt_lay_indices))))
                                #   if the list is empty the lenght is also empty
                                except ZeroDivisionError:
                                    continue
                                aqt_lay_start = aqt_lay_indices[start_idx_g]
                                aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                                #   change the values in the array
                                self.hk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val']
                                self.vk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val'] / 10.  # this has to be changed later on!!!
                            #   at the edges of the model sometimes the domain is so thin it cant be split into desire number of aquitard layers
                            except IndexError:
                                continue


        #   create a random amount of layers in the offshore part of the self domain
        if rand_aqt_offshore:
            aqt_off_x1 = self.x_end - 15.1
            #   define the domain where the aquitards occure - between the soil layer and the bottom
            #   in the list below the range of layers will be stored - based on the sand/clay ratio
            aqt_bot_top_lst = []
            #for i in range(self.x_coast_idx, self.ncol):
            for i in range(start_idx, self.ncol):
                #   find the number of layers with the aquifer Hk - this will make sure aquitard is in the aquifer and not the soil layer
                col_i = self.hk_arr[:, :, i].tolist()
                indices = [a for a, x in enumerate(col_i) if x == [Hk_val]]
                #   append to the bot_top_lst the slices of the indices list, depending on total number of aquitard layers
                inner_lst = list(split(indices, n_aqt_off))
                """
                inner_lst = [i]
                idx_top_aqt = 0
                for j in xrange(n_aqt_off):
                    if j != n_aqt_off - 1:
                        inner_lst.append(indices[idx_top_aqt : len(indices) / n_aqt_off + idx_top_aqt + 1])
                        idx_top_aqt += len(indices) / n_aqt_off + 1
                    else:
                        inner_lst.append(indices[idx_top_aqt :])
                """
                aqt_bot_top_lst.append(inner_lst)

            #   loop through the number of aquitard layers inland
            for x in range(n_aqt_off):
                #   get random extent in the horizontal direction
                aqt_extent_1 = round(random.uniform(aqt_off_x0, 5.0), 1)
                aqt_extent_2 = round(random.uniform(aqt_off_x1, self.x_end), 1)
                #   find the column indexes of the start and end of the aquitard layer
                aqt_start = min(aqt_extent_1, aqt_extent_2)
                aqt_end = max(aqt_extent_1, aqt_extent_2)

                #aqt_col_start = self.x_coast_idx + int(abs((aqt_start) * 1000. / self.del_col))
                #aqt_col_end = self.x_coast_idx + int(abs((aqt_end) * 1000. / self.del_col))

                aqt_col_start = start_idx + int(abs((aqt_start) * 1000. / self.del_col))
                aqt_col_end = start_idx + int(abs((aqt_end) * 1000. / self.del_col))

                random_layer = False
                linear_layer = True

                #   randomly distributed layer in vertical direction
                if random_layer:
                    #   loop through the aquitard layer column extent and change the hydraylic conductivity to clay in corresponding cells
                    for g in range(aqt_col_start, aqt_col_end):
                        #aqt_lay_indices = aqt_bot_top_lst[g - self.x_coast_idx][x]
                        print (len(aqt_bot_top_lst) / 2, len(aqt_bot_top_lst), x)
                        aqt_lay_indices = aqt_bot_top_lst[g - start_idx][x]
                        #   split the indices based on the sand/clay ratio - to make sure the full thickness is possible. If it is 50/50
                        #   then the upper limit of the aquitard layer should be in the upper half of the indices so full thickness is possible
                        clay_ratio = 100 / clay_pct
                        try:
                            aqt_lay_start = random.choice(aqt_lay_indices[: (-1 * (len(aqt_lay_indices) / clay_ratio))])
                            aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                            #   change the values in the array
                            self.hk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val']
                            self.vk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val'] / 10.  # this has to be changed later on!!!
                        #   at the edges of the model sometimes the domain is so thin it cant be split into desire number of aquitard layers
                        except IndexError:
                            continue
                #   linearly distributed layer in vertical direction
                elif linear_layer:
                    #   select the column in the middle
                    aqt_lay_indices = aqt_bot_top_lst[len(aqt_bot_top_lst) / 2][x]
                    #   split the indices based on the sand/clay ratio - to make sure the full thickness is possible. If it is 50/50
                    #   then the upper limit of the aquitard layer should be in the upper half of the indices so full thickness is possible
                    clay_ratio = 100 / clay_pct
                    aqt_lay_start = random.choice(aqt_lay_indices[: (-1 * (len(aqt_lay_indices) / clay_ratio))])
                    aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                    #   express the position of the first layer as distance from the upper limit
                    start_ratio = int(round(aqt_lay_indices.index(aqt_lay_start) / float(len(aqt_lay_indices)) * 100))
                    #   loop through the aquitard layer column extent and change the hydraylic conductivity to clay in corresponding cells
                    for g in range(aqt_col_start, aqt_col_end):
                        try:
                            #aqt_lay_indices = aqt_bot_top_lst[g - self.x_coast_idx][x]
                            aqt_lay_indices = aqt_bot_top_lst[g - start_idx][x]
                            #   find the starting index based on the ratio
                            try:
                                start_idx_g = int(round(start_ratio / float(len(aqt_lay_indices))))
                            #   if the list is empty the lenght is also empty
                            except ZeroDivisionError:
                                continue
                            aqt_lay_start = aqt_lay_indices[start_idx_g]
                            aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                            #   change the values in the array
                            self.hk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val']
                            self.vk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val'] / 10.  # this has to be changed later on!!!
                        #   at the edges of the model sometimes the domain is so thin it cant be split into desire number of aquitard layers
                        except IndexError:
                            continue






    #   Define the input for the LPF package - heterogeneous geology
    #   rand_seed   : value that ensures the reproducability of the randomization process
    #   Hk_val      : horizontal conductivity values
    #   Vk_val      : vertical conductivity values
    #   laytyp_val  : MODFLOW layer type
    #   sand_pct    : sand percentage in the sedimental body
    #   clay_pct    : clay percentage in the sedimental body
    #   n_aqt_in    : number of aquitard layers inland
    #   aqt_in_x0   : the inland limit of the extent of aquitard layers (in -km from coastline)
    #   n_aqt_off   : number of aquitard layers offshore
    #   aqt_off_x1  : the offshore limit of the extent of aquitard layers (in km from coastline)
    def lpf_input_complex_geo_GLHYMPS_rand(self, rand_seed, laytyp_val, sand_pct, clay_pct, n_aqt_in, aqt_in_x1, n_aqt_off, aqt_off_x0,
                              soil_thk = 10.0, aq_vals_smooth_hor = False, top_soil = False, top_offshore = True, rand_aqt_inland = True, rand_aqt_offshore = True):

        def split(a, n):
            k, m = divmod(len(a), n)
            return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

        #self.soil_thk = self.soil_thk[self.idx_start : self.idx_end]
        soil_thk_lst = [p for p in self.soil_thk if p > -1]
        #   repeat the elements of each list n times depending on the col width
        soil_thk_lst = np.repeat(soil_thk_lst, int(500. / self.del_col)).tolist()

        #   set the random.seed
        random.seed(a = rand_seed)

        self.hk_arr = np.zeros([self.nlay, 1, self.ncol])


        #self.soil_thk = self.soil_thk[self.idx_start : len(self.top_elev)]  
        #self.soil_type = self.soil_type[self.idx_start : len(self.top_elev)]  
        
        self.soil_thk = self.soil_thk[self.idx_start : self.cst_idx]  
        self.soil_type = self.soil_type[self.idx_start : self.cst_idx]          

        #   get rid of the negative values - why are they even there?
        self.hk_vals_bot = self.glhymps_bot_lay
        self.hk_vals_top = self.glhymps_top_lay

        #   recalculate the values based on the GLHYMPS 2 formula
        self.hk_vals_top = [i * 3600 *24 for i in self.glhymps_top_lay]

        #   replace the negative values
        for x in range(len(self.hk_vals_top)):
            if self.hk_vals_top[x] < 0.:
               self.hk_vals_top[x] = np.mean([i for i in self.hk_vals_top if i > 0])
            if self.hk_vals_bot[x] < 0.:
               self.hk_vals_bot[x] = np.mean([i for i in self.hk_vals_bot if i > 0])

        #   create an array of hydraulic conductivity with all values = aquifer Kh
        self.laytyp = laytyp_val

        #   find out how thick is the soil layer and what type of soil it is
        self.top_lay_thk = self.soil_thk[self.idx_start : self.idx_end]
        #top_lay_thk_lst = [j for j in self.top_lay_thk if j > -1]

        #start_idx = int(abs(self.x_start * 1000. / self.del_col))
        start_idx = (self.cst_idx - self.idx_start) * 5 

        #   loop through the raw list and assign each value as many times as the ratio 500./col_width
        col_ratio = int(500. / self.del_col)

        hk_top_mean = np.nanmean(self.hk_vals_top, axis=0)
        hk_top_std = np.nanstd(self.hk_vals_top, axis=0)

        hk_bot_mean = np.nanmean(self.hk_vals_bot, axis=0)
        hk_bot_std = np.nanstd(self.hk_vals_bot, axis=0)

        """
        #   get the top layer thickness based on Soilgrids
        col_idx = 0
        for r in xrange(len(self.soil_thk)):
            print r, col_idx
            if col_idx <= start_idx:
                top_lay_thk = round(self.soil_thk[r] / 100., 0)    # /100. because the raw values are in centimeters
                col_top_lay_n = int(top_lay_thk / self.del_lay)

                ibound_col_idxs = [i for i, e in enumerate(self.ibound_arr[:,:,col_idx]) if e != 0]

                neg_vals = True
                while neg_vals:
                    try:
                        #   added abs to only get positive values
                        rand_hk_top_vals = abs(np.random.normal(hk_top_mean, hk_top_std, len(ibound_col_idxs)))
                        rand_hk_bot_vals = abs(np.random.normal(hk_bot_mean, hk_bot_std, len(ibound_col_idxs)))
                        if sum(n < 0 for n in rand_hk_top_vals) == 0 and sum(n < 0 for n in rand_hk_bot_vals) == 0:
                            neg_vals = False
                        else:
                            pass
                    except ValueError:
                        rand_hk_top_vals = np.random.normal(hk_top_mean, 0.001, len(ibound_col_idxs))
                        rand_hk_bot_vals = np.random.normal(hk_bot_mean, 0.001, len(ibound_col_idxs))
                        if sum(n < 0 for n in rand_hk_top_vals) == 0 and sum(n < 0 for n in rand_hk_bot_vals) == 0:
                            neg_vals = False
                        else:
                            pass
                    #   check that there are no negative values in both randomized samples, if yes then randomize again
                    #   loop through the rows in the column that are active
                    for g in xrange(len(ibound_col_idxs)):
                        if g < col_top_lay_n:
                            if rand_hk_top_vals[g] > 0:
                                self.hk_arr[ibound_col_idxs[g], 0, col_idx] = rand_hk_top_vals[g]
                                g =+ 1
                            else:
                                print rand_hk_top_vals[g], r, g, col_idx
                        else:
                            if rand_hk_bot_vals[g] > 0:
                                self.hk_arr[ibound_col_idxs[g], 0, col_idx] = rand_hk_bot_vals[g]
                                g =+ 1
                            else:
                                print rand_hk_bot_vals[g], r, g, col_idx
                    col_idx += 1
            else:
                #for k in xrange(col_ratio):
                ibound_col_idxs = [i for i, e in enumerate(self.ibound_arr[:,:,col_idx]) if e != 0]
                rand_hk_bot_vals = np.random.normal(hk_bot_mean, hk_bot_std, len(ibound_col_idxs))

                neg_vals = True
                while neg_vals:
                    try:

                        rand_hk_bot_vals = abs(np.random.normal(hk_bot_mean, hk_bot_std, len(ibound_col_idxs)))
                        #print rand_hk_bot_vals
                        if sum(n < 0 for n in rand_hk_bot_vals) == 0:
                            neg_vals = False
                        else:
                            pass
                    except ValueError:
                        rand_hk_bot_vals = np.random.normal(hk_bot_mean, 0.001, len(ibound_col_idxs))
                        if sum(n < 0 for n in rand_hk_bot_vals) == 0:
                            neg_vals = False
                        else:
                            pass
                    #   loop through the rows in the column that are active
                    for g in xrange(len(ibound_col_idxs)):
                            if rand_hk_bot_vals[g] > 0:
                                self.hk_arr[ibound_col_idxs[g], 0, col_idx] = rand_hk_bot_vals[g]
                                g =+ 1
                            else:
                                print rand_hk_bot_vals[g], r, g, col_idx
                    col_idx += 1


        """
        #   get the top layer thickness based on Soilgrids
        col_idx = 0
        for r in range(len(self.soil_thk)):
        #for r in xrange(self.ibound_arr.shape[-1] / col_ratio):
            if col_idx <= start_idx:
                top_lay_thk = round(self.soil_thk[r] / 100., 0)    # /100. because the raw values are in centimeters
                col_top_lay_n = int(top_lay_thk / self.del_lay)
                for k in range(col_ratio):
                    ibound_col_idxs = [i for i, e in enumerate(self.ibound_arr[:,:,col_idx]) if e != 0]

                    neg_vals = True
                    while neg_vals:
                        try:
                            #   added abs to only get positive values
                            rand_hk_top_vals = abs(np.random.normal(hk_top_mean, hk_top_std, len(ibound_col_idxs)))
                            rand_hk_bot_vals = abs(np.random.normal(hk_bot_mean, hk_bot_std, len(ibound_col_idxs)))
                            if sum(n < 0 for n in rand_hk_top_vals) == 0 and sum(n < 0 for n in rand_hk_bot_vals) == 0:
                                neg_vals = False
                            else:
                                pass
                        except ValueError:
                            rand_hk_top_vals = np.random.normal(hk_top_mean, 0.001, len(ibound_col_idxs))
                            rand_hk_bot_vals = np.random.normal(hk_bot_mean, 0.001, len(ibound_col_idxs))
                            if sum(n < 0 for n in rand_hk_top_vals) == 0 and sum(n < 0 for n in rand_hk_bot_vals) == 0:
                                neg_vals = False
                            else:
                                pass
                    #   check that there are no negative values in both randomized samples, if yes then randomize again
                    #   loop through the rows in the column that are active
                    for g in range(len(ibound_col_idxs)):
                        if g < col_top_lay_n:
                            if rand_hk_top_vals[g] > 0:
                                self.hk_arr[ibound_col_idxs[g], 0, col_idx] = rand_hk_top_vals[g]
                                g =+ 1
                            else:
                                print (rand_hk_top_vals[g], r, g, col_idx)
                        else:
                            if rand_hk_bot_vals[g] > 0:
                                self.hk_arr[ibound_col_idxs[g], 0, col_idx] = rand_hk_bot_vals[g]
                                g =+ 1
                            else:
                                print (rand_hk_bot_vals[g], r, g, col_idx)
                    col_idx += 1
            else:
                pass
            
        for q in range(col_idx, self.ibound_arr.shape[2]):     
            
            ibound_col_idxs = [i for i, e in enumerate(self.ibound_arr[:,:,q]) if e != 0]
            rand_hk_bot_vals = np.random.normal(hk_bot_mean, hk_bot_std, len(ibound_col_idxs))

            neg_vals = True
            while neg_vals:
                try:

                    rand_hk_bot_vals = abs(np.random.normal(hk_bot_mean, hk_bot_std, len(ibound_col_idxs)))
                    #print rand_hk_bot_vals
                    if sum(n < 0 for n in rand_hk_bot_vals) == 0:
                        neg_vals = False
                    else:
                        pass
                except ValueError:
                    rand_hk_bot_vals = np.random.normal(hk_bot_mean, 0.001, len(ibound_col_idxs))
                    if sum(n < 0 for n in rand_hk_bot_vals) == 0:
                        neg_vals = False
                    else:
                        pass
            #   loop through the rows in the column that are active
            for g in range(len(ibound_col_idxs)):
                    if rand_hk_bot_vals[g] > 0:
                        self.hk_arr[ibound_col_idxs[g], 0, q] = rand_hk_bot_vals[g]
                        g =+ 1
                    else:
                        print (rand_hk_bot_vals[g], r, g, q)
          
        
            
            
            
        """      
                for k in xrange(col_ratio):
                    ibound_col_idxs = [i for i, e in enumerate(self.ibound_arr[:,:,col_idx]) if e != 0]
                    rand_hk_bot_vals = np.random.normal(hk_bot_mean, hk_bot_std, len(ibound_col_idxs))

                    neg_vals = True
                    while neg_vals:
                        try:

                            rand_hk_bot_vals = abs(np.random.normal(hk_bot_mean, hk_bot_std, len(ibound_col_idxs)))
                            #print rand_hk_bot_vals
                            if sum(n < 0 for n in rand_hk_bot_vals) == 0:
                                neg_vals = False
                            else:
                                pass
                        except ValueError:
                            rand_hk_bot_vals = np.random.normal(hk_bot_mean, 0.001, len(ibound_col_idxs))
                            if sum(n < 0 for n in rand_hk_bot_vals) == 0:
                                neg_vals = False
                            else:
                                pass
                    #   loop through the rows in the column that are active
                    for g in xrange(len(ibound_col_idxs)):
                            if rand_hk_bot_vals[g] > 0:
                                self.hk_arr[ibound_col_idxs[g], 0, col_idx] = rand_hk_bot_vals[g]
                                g =+ 1
                            else:
                                print rand_hk_bot_vals[g], r, g, col_idx
                    col_idx += 1
        #self.hk_arr[self.hk_arr == 0.0] = np.nan
        #plt.imshow(self.hk_arr[:,0,:])
        
        
        #   now loop through the rest of the initial list and assign Hk values depending on the type
        #       1) if it is False then just assign the raw values, no averaging or smoothening
        if aq_vals_smooth_hor is False:
            #   first fill the Hk array with aquifer values from GLHYMPS and GLHYMPS 2.0
            k_vals_top, k_vals_bot = [], []
            for i in xrange(len(self.hk_vals_top)):
                for j in xrange(col_ratio):
                    k_vals_top.append(round(self.hk_vals_top[i], 2))
                    k_vals_bot.append(round(self.hk_vals_top[i], 2))

        #       2) if it is True then smoothen out the values in the columns between the cross-sectional points
        elif aq_vals_smooth_hor is True:
            #   first fill the Hk array with aquifer values from GLHYMPS and GLHYMPS 2.0
            k_vals_top, k_vals_bot = [], []
            for j in xrange(3):
                k_vals_top.append(round(self.hk_vals_top[0], 2))
                k_vals_bot.append(round(self.hk_vals_top[0], 2))

            for i in xrange(len(self.hk_vals_top) - 1):
                for j in xrange(col_ratio - 1):
                    #   calculate the smoothen value for given column
                    k_val_top = self.hk_vals_top[i] - j * ((self.hk_vals_top[i] - self.hk_vals_top[i + 1]) / (col_ratio - 1))
                    k_val_bot = self.hk_vals_top[i] - j * ((self.hk_vals_bot[i] - self.hk_vals_bot[i + 1]) / (col_ratio - 1))
                    #   append to the final lists
                    k_vals_top.append(round(k_val_top, 2))
                    k_vals_bot.append(round(k_val_bot, 2))
                #   append the cross-section points values as well
                k_vals_top.append(round(self.hk_vals_top[i + 1], 2))
                k_vals_bot.append(round(self.hk_vals_bot[i + 1], 2))
            #   at the end append the final two values
            for j in xrange(2):
                k_vals_top.append(round(self.hk_vals_top[-1], 2))
                k_vals_bot.append(round(self.hk_vals_top[-1], 2))


        #   next, based on the preferences fill in the Hk values per column, either smoothened out between the end of the more permeable
        #   upper layer and the bottom of the self or just fill in the top and bottom Hk values based on the thickness
        if aq_vals_smooth_ver is False:

            #   get the top layer thickness based on Soilgrids
            col_idx = 0
            for j in xrange(len(self.soil_thk)):
                top_lay_thk = round(self.soil_thk[j] / 100., 0)    # /100. because the raw values are in centimeters
                col_top_lay_n = int(top_lay_thk / self.del_lay)
                for k in xrange(col_ratio):
                    ibound_col_idxs = [i for i, e in enumerate(self.ibound_arr[:,:,col_idx]) if e != 0]
                    #   loop through the rows in the column that are active
                    for g in xrange(len(ibound_col_idxs)):
                        if g < col_top_lay_n:
                            self.hk_arr[ibound_col_idxs[g], 0, col_idx] = k_vals_top[col_idx]
                            g =+ 1
                        else:
                            self.hk_arr[ibound_col_idxs[g] : ibound_col_idxs[-1], 0, col_idx] = k_vals_bot[col_idx]
                    col_idx += 1

            plt.imshow(self.hk_arr[:,0,:])

            #   loop through each column and check the number of layers in there
            for i in xrange(len(k_vals_top)):
                ibound_col_idxs = [i for i, e in enumerate(self.ibound_arr[:,:,i]) if e != 0]
        """

        self.vk_arr = self.hk_arr / 10.

        #   create the top soil layer - this is different from the previous function because here we use the SOILGRIDS thickness
        #   for the highly conductive part of the unconsolidates aquifer itself.. So in this case the top soil layer would only
        #   be few top meters at the surface representing the unsaturated zone (sort of), this is by default turned off
        if top_soil:
            #   loop through all the self columns
            for a in range(len(soil_type_lst)):
                #   select the column from the ibound array
                col_a = self.ibound_arr[:, :, a].tolist()
                indices = [i for i, x in enumerate(col_a) if x == [1]]
                #   find how many layers will have the soil K value based on the thickness
                try:
                    idx_count = int(round(soil_thk_lst[a] / (100 * self.del_lay))) # * 100 because the thickness is in cm
                except IndexError:  #   in case the list of soil types is longer than soil thicknesses, take the last thickness
                    idx_count = int(round(soil_thk_lst[-1] / (100 * self.del_lay)))
                try:
                    #   for the number of indices change the K values in the arrays
                    for b in range(idx_count):
                        self.hk_arr[indices[b], :, a] = k_soil_dict[soil_type_lst[a]]['k_val']
                        self.vk_arr[indices[b], :, a] = k_soil_dict[soil_type_lst[a]]['k_val'] / 10.  # this has to be changed later on!!!
                #   in case no cell is active in the given column of the ibound array
                except IndexError:
                    continue

        #   create a clay (mud) layer on top of the continental shelf
        if top_offshore:
            #   loop through all the self columns
            for c in range(self.cst_idx, self.ncol):
                #   select the column from the ibound array
                col_c = self.ibound_arr[:, :, c].tolist()
                indices = [i for i, x in enumerate(col_c) if x == [1]]
                #   take the last thickness
                idx_count = int(round(soil_thk_lst[-1] / (100 * self.del_lay)))
                try:
                    #   for the number of indices change the K values in the arrays
                    if self.top_elev[c] > -120.:
                        for b in range(idx_count):
                            self.hk_arr[indices[b], :, c] = k_soil_dict[1]['k_val']
                            self.vk_arr[indices[b], :, c] = k_soil_dict[1]['k_val'] / 10.  # this has to be changed later on!!!
                #   in case no cell is active in the given column of the ibound array
                except IndexError:
                    continue

        #   create a random amount of layers in the inland part of the self domain
        if rand_aqt_inland:
            aqt_in_x0 = self.x_start + 5.1
            #   define the domain where the aquitards occure - between the soil layer and the bottom
            #   in the list below the range of layers will be stored - based on the sand/clay ratio
            aqt_bot_top_lst = []
            #for i in range(0, self.x_coast_idx):
            for i in range(0, start_idx):

                if not top_soil:

                   #   find the number of layers with the aquifer Hk - this will make sure aquitard is in the aquifer and not the soil layer
                    col_i = self.hk_arr[:, :, i].tolist()
                    indices = [a for a, x in enumerate(col_i) if x != [0.]]#[2:]
                    #   append to the bot_top_lst the slices of the indices list, depending on total number of aquitard layers
                    inner_lst = list(split(indices, n_aqt_in))
                    """
                    inner_lst = [i]
                    idx_top_aqt = 0
                    for j in xrange(n_aqt_in):
                        if j != n_aqt_in - 1:
                            inner_lst.append(indices[idx_top_aqt : len(indices) / n_aqt_in + idx_top_aqt + 1])
                            idx_top_aqt += len(indices) / n_aqt_in + 1
                        else:
                            inner_lst.append(indices[idx_top_aqt :])
                    """
                    aqt_bot_top_lst.append(inner_lst)


                else:
                   #   find the number of layers with the aquifer Hk - this will make sure aquitard is in the aquifer and not the soil layer
                    col_i = self.hk_arr[:, :, i].tolist()
                    indices = [a for a, x in enumerate(col_i) if x == [Hk_val]]
                    #   append to the bot_top_lst the slices of the indices list, depending on total number of aquitard layers
                    inner_lst = list(split(indices, n_aqt_in))
                    """
                    inner_lst = [i]
                    idx_top_aqt = 0
                    for j in xrange(n_aqt_in):
                        if j != n_aqt_in - 1:
                            inner_lst.append(indices[idx_top_aqt : len(indices) / n_aqt_in + idx_top_aqt + 1])
                            idx_top_aqt += len(indices) / n_aqt_in + 1
                        else:
                            inner_lst.append(indices[idx_top_aqt :])
                    """
                    aqt_bot_top_lst.append(inner_lst)

            #   loop through the number of aquitard layers inland
            for j in range(n_aqt_in):
                #   get random extent in the horizontal direction
                aqt_extent_1 = round(random.uniform(aqt_in_x0, self.x_start/2), 1)
                aqt_extent_2 = round(random.uniform(self.x_start/4, aqt_in_x1), 1)
                #   find the column indexes of the start and end of the aquitard layer
                aqt_start = min(aqt_extent_1, aqt_extent_2)
                aqt_end = max(aqt_extent_1, aqt_extent_2)
                aqt_col_start = int(abs((self.x_start - aqt_start) * 1000. / self.del_col))
                aqt_col_end = int(abs((self.x_start - aqt_end) * 1000. / self.del_col))

                random_layer = False
                linear_layer = True

                #   randomly distributed layer in vertical direction
                if random_layer:
                    #   loop through the aquitard layer column extent and change the hydraylic conductivity to clay in corresponding cells
                    for g in range(aqt_col_start, aqt_col_end):
                        aqt_lay_indices = aqt_bot_top_lst[g][j]
                        #   split the indices based on the sand/clay ratio - to make sure the full thickness is possible. If it is 50/50
                        #   then the upper limit of the aquitard layer should be in the upper half of the indices so full thickness is possible
                        clay_ratio = 100 / clay_pct
                        aqt_lay_start = random.choice(aqt_lay_indices[: (-1 * (len(aqt_lay_indices) / clay_ratio))])
                        aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                        #   change the values in the array
                        self.hk_arr[aqt_lay_start : aqt_lay_end, :, g] = 1e-05
                        self.vk_arr[aqt_lay_start : aqt_lay_end, :, g] = 1e-03  # this has to be changed later on!!!

                #   linearly distributed layer in vertical direction
                elif linear_layer:
                    #   select the column in the middle
                    aqt_lay_indices = aqt_bot_top_lst[len(aqt_bot_top_lst) / 2][j]

                    if aqt_lay_indices == []:
                        for g in range(len(aqt_bot_top_lst) / 2, len(aqt_bot_top_lst)):
                            not_empty = 0
                            for h in range(n_aqt_in):
                                if aqt_bot_top_lst[g][h] != []:
                                    not_empty += 1
                                else:
                                    pass
                            if not_empty == n_aqt_in:
                                aqt_lay_indices = aqt_bot_top_lst[g][j]
                                break

                    clay_ratio = clay_pct / 100.
                    try:
                        #   split the indices based on the sand/clay ratio - to make sure the full thickness is possible. If it is 50/50
                        #   then the upper limit of the aquitard layer should be in the upper half of the indices so full thickness is possible
                        aqt_lay_start = random.choice(aqt_lay_indices[: int((-1 * (len(aqt_lay_indices) / (100 * clay_ratio))))])
                        aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                    except IndexError:
                        aqt_lay_start = aqt_lay_indices[0]
                        nr_aqt_layers = int(round(len(aqt_lay_indices) * clay_ratio, 0))
                        aqt_lay_end = aqt_lay_indices[nr_aqt_layers - 1] + (aqt_lay_start - aqt_lay_indices[0])

                    #   express the position of the first layer as distance from the upper limit
                    start_ratio = round(aqt_lay_indices.index(aqt_lay_start) / float(len(aqt_lay_indices)))# * 100))
                    #   loop through the aquitard layer column extent and change the hydraylic conductivity to clay in corresponding cells
                    for g in range(aqt_col_start, aqt_col_end):
                        try:
                            aqt_lay_indices = aqt_bot_top_lst[g][j]
                            #   find the starting index based on the ratio
                            try:
                                start_idx_g = int(round(start_ratio / float(len(aqt_lay_indices))))
                            #   if the list is empty the lenght is also empty
                            except ZeroDivisionError:
                                continue
                            aqt_lay_start = aqt_lay_indices[start_idx_g]
                            aqt_lay_end = aqt_lay_indices[int(round(len(aqt_lay_indices) * clay_ratio, 0))] + (aqt_lay_start - aqt_lay_indices[0])

                            #aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) * clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                            #   change the values in the array
                            self.hk_arr[aqt_lay_start : aqt_lay_end, :, g] = 1e-05
                            self.vk_arr[aqt_lay_start : aqt_lay_end, :, g] = 1e-03  # this has to be changed later on!!!
                        #   at the edges of the self sometimes the domain is so thin it cant be split into desire number of aquitard layers
                        except IndexError:
                            continue

        #   create a random amount of layers in the offshore part of the self domain
        if rand_aqt_offshore:
            aqt_off_x1 = self.x_end - 15.1
            #   define the domain where the aquitards occure - between the soil layer and the bottom
            #   in the list below the range of layers will be stored - based on the sand/clay ratio
            aqt_bot_top_lst = []
            #for i in range(self.x_coast_idx, self.ncol):
            for i in range(start_idx, self.ncol):
                #   find the number of layers with the aquifer Hk - this will make sure aquitard is in the aquifer and not the soil layer
                col_i = self.hk_arr[:, :, i].tolist()
                indices = [a for a, x in enumerate(col_i) if x !=[0.0]]
                #   append to the bot_top_lst the slices of the indices list, depending on total number of aquitard layers
                inner_lst = list(split(indices, n_aqt_off))
                """
                inner_lst = [i]
                idx_top_aqt = 0
                for j in xrange(n_aqt_off):
                    if j != n_aqt_off - 1:
                        inner_lst.append(indices[idx_top_aqt : len(indices) / n_aqt_off + idx_top_aqt + 1])
                        idx_top_aqt += len(indices) / n_aqt_off + 1
                    else:
                        inner_lst.append(indices[idx_top_aqt :])
                """
                aqt_bot_top_lst.append(inner_lst)

            #   loop through the number of aquitard layers inland
            for x in range(n_aqt_off):
                #   get random extent in the horizontal direction
                aqt_extent_1 = round(random.uniform(aqt_off_x0, 5.0), 1)
                aqt_extent_2 = round(random.uniform(aqt_off_x1, self.x_end), 1)
                #   find the column indexes of the start and end of the aquitard layer
                aqt_start = min(aqt_extent_1, aqt_extent_2)
                aqt_end = max(aqt_extent_1, aqt_extent_2)

                #aqt_col_start = self.x_coast_idx + int(abs((aqt_start) * 1000. / self.del_col))
                #aqt_col_end = self.x_coast_idx + int(abs((aqt_end) * 1000. / self.del_col))

                aqt_col_start = start_idx + int(abs((aqt_start) * 1000. / self.del_col))
                aqt_col_end = start_idx + int(abs((aqt_end) * 1000. / self.del_col))

                random_layer = False
                linear_layer = True

                #   randomly distributed layer in vertical direction
                if random_layer:
                    #   loop through the aquitard layer column extent and change the hydraylic conductivity to clay in corresponding cells
                    for g in range(aqt_col_start, aqt_col_end):
                        #aqt_lay_indices = aqt_bot_top_lst[g - self.x_coast_idx][x]
                        aqt_lay_indices = aqt_bot_top_lst[g - start_idx][x]
                        #   split the indices based on the sand/clay ratio - to make sure the full thickness is possible. If it is 50/50
                        #   then the upper limit of the aquitard layer should be in the upper half of the indices so full thickness is possible
                        clay_ratio = 100 / clay_pct
                        try:
                            aqt_lay_start = random.choice(aqt_lay_indices[: (-1 * (len(aqt_lay_indices) / clay_ratio))])
                            aqt_lay_end = aqt_lay_indices[len(aqt_lay_indices) / clay_ratio] + (aqt_lay_start - aqt_lay_indices[0])
                            #   change the values in the array
                            self.hk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val']
                            self.vk_arr[aqt_lay_start : aqt_lay_end, :, g] = k_soil_dict[1]['k_val'] / 10.  # this has to be changed later on!!!
                        #   at the edges of the self sometimes the domain is so thin it cant be split into desire number of aquitard layers
                        except IndexError:
                            continue
                #   linearly distributed layer in vertical direction
                elif linear_layer:
                    #   select the column in the middle
                    aqt_lay_indices = aqt_bot_top_lst[len(aqt_bot_top_lst) / 2][x]
                    #   split the indices based on the sand/clay ratio - to make sure the full thickness is possible. If it is 50/50
                    #   then the upper limit of the aquitard layer should be in the upper half of the indices so full thickness is possible

                    clay_ratio = clay_pct / 100.
                    try:
                        #   split the indices based on the sand/clay ratio - to make sure the full thickness is possible. If it is 50/50
                        #   then the upper limit of the aquitard layer should be in the upper half of the indices so full thickness is possible
                        aqt_lay_start = random.choice(aqt_lay_indices[: int((-1 * (len(aqt_lay_indices) / (100 * clay_ratio))))])
                        aqt_lay_end = aqt_lay_indices[int(len(aqt_lay_indices) / clay_ratio)] + (aqt_lay_start - aqt_lay_indices[0])
                    except IndexError:
                        aqt_lay_start = aqt_lay_indices[0]
                        nr_aqt_layers = int(round(len(aqt_lay_indices) * clay_ratio, 0))
                        aqt_lay_end = aqt_lay_indices[nr_aqt_layers - 1] + (aqt_lay_start - aqt_lay_indices[0])

                    #   express the position of the first layer as distance from the upper limit
                    start_ratio = int(round(aqt_lay_indices.index(aqt_lay_start) / float(len(aqt_lay_indices)) * 100))
                    #   loop through the aquitard layer column extent and change the hydraylic conductivity to clay in corresponding cells
                    for g in range(aqt_col_start, aqt_col_end):
                        try:
                            #aqt_lay_indices = aqt_bot_top_lst[g - self.x_coast_idx][x]
                            aqt_lay_indices = aqt_bot_top_lst[g - start_idx][x]
                            #   find the starting index based on the ratio
                            try:
                                start_idx_g = int(round(start_ratio / float(len(aqt_lay_indices))))
                            #   if the list is empty the lenght is also empty
                            except ZeroDivisionError:
                                continue
                            aqt_lay_start = aqt_lay_indices[start_idx_g]
                            aqt_lay_end = aqt_lay_indices[int(round(len(aqt_lay_indices) * clay_ratio, 0))] + (aqt_lay_start - aqt_lay_indices[0])
                            #   change the values in the array
                            self.hk_arr[aqt_lay_start : aqt_lay_end, :, g] = 1e-05
                            self.vk_arr[aqt_lay_start : aqt_lay_end, :, g] = 1e-03  # this has to be changed later on!!!
                        #   at the edges of the self sometimes the domain is so thin it cant be split into desire number of aquitard layers
                        except IndexError:
                            continue



    #   Define the input for the GHB package
    def ghb_input(self, sea_lvl, cond_fact, clay_cond_val, cond_limit_val, clay_incr = True, offshore_last_col_ghb = False, inland_ghb = False,\
                  offshore_last_col_ghb_fresh = False, offshore_ghb = False, wtd_inl_lvl = True):
        #   because laycon = 0 we need to specify the vcont and tran (optional)
        #tran = np.zeros(((self.nlay - 1), 1, self.ncol), dtype = np.float)
        tran = np.zeros(((self.ibound_arr.shape[0] - 1), 1, self.ibound_arr.shape[-1]), dtype = np.float)
        #tran = np.dot(self.hk_arr, self.delv)
        tran = np.dot(self.hk_arr, self.del_lay)
        tran = np.round(tran, 2)
        tran[tran == 0.] = 0.01
        #   define and initiate the output arrays and lists
        self.ghb_arr = self.ibound_arr * 1
        self.ghb_input_lst = []
        if clay_incr:
            self.cond_val = (tran * cond_fact) / self.delc
            self.cond_val_no_incr = (tran * cond_fact) / self.delc
            #   for all cells that are lower than 0.01 (which is the conductance of clay cell * 100 - to be on the safe side) assign the clay_cond_val
            self.cond_val[self.cond_val < cond_limit_val] = clay_cond_val
        else:
            self.cond_val = (tran * cond_fact) / self.delc
        self.ssmdata = []
        #  create the SSM dictionary where the ssm input will be written to
        itype = flopy.mt3d.Mt3dSsm.itype_dict()

        if offshore_ghb:
            #   assign the sea water boundary at the last column cells, only in active cells
            #for a in range(1, self.ncol):
            #for a in range(1, self.ibound_arr.shape[-1]):
            cst_col = len(self.top_elev)
            for f in reversed(self.top_elev):
                if f < sea_lvl:
                    cst_col -= 1
                else:
                    break            
            for a in range(cst_col, len(self.top_elev)):
                #   get the index of the first non-zero ibound cell
                try:
                    top_cell_lay_idx = self.ibound_arr[:, 0, a].tolist().index(1)
                    topo_val = self.top_elev[a]
                    if topo_val < sea_lvl:
                        #self.ghb_input_lst.append([top_cell_lay_idx, 0, a, self.top_elev[a], self.cond_val[top_cell_lay_idx, 0, a]])   #   sea level is 0.0m (present condition)
                        self.ghb_input_lst.append([top_cell_lay_idx, 0, a, sea_lvl, self.cond_val[top_cell_lay_idx, 0, a]])   #   sea level is 0.0m (present condition)
                        self.ghb_arr[top_cell_lay_idx, 0, a] = -1
                        self.ssmdata.append([top_cell_lay_idx, 0, a, 35.0, itype['GHB']])
                except ValueError:
                    pass

        if inland_ghb:
            #   add a GHB cell in the cells on the fresh-water boundary, only in active cells
            #for b in range(self.nlay):
            for b in range(self.ibound_arr.shape[0]):
                if self.ibound_arr[b, 0, 0] == 1:
                    if wtd_inl_lvl:
                        self.ghb_input_lst.append([b, 0, 0, max(0, self.wtd_elev[0]), self.cond_val[b,0,0]])
                    else:
                        self.ghb_input_lst.append([b, 0, 0, self.top_elev[0], self.cond_val[b,0,0]])
                    self.ghb_arr[b, 0, 0] = -1
                    self.ssmdata.append([b, 0, 0, 0.0, itype['GHB']])

        if offshore_last_col_ghb:
            #   add GHB cells in the last offshore model column 
            #for c in range(self.nlay):
            for c in range(self.ibound_arr.shape[0]):
                #if self.ibound_arr[c, 0, self.ncol - 1] == 1:
                #    self.ghb_input_lst.append([c, 0, self.ncol - 1, sea_lvl, self.cond_val[c][0][self.ncol - 1]])
                #    self.ghb_arr[c, 0, self.ncol - 1] = -1
                #    self.ssmdata.append([c, 0, self.ncol - 1, 35.0, itype['GHB']])      
                if self.ibound_arr[c, 0, self.ibound_arr.shape[-1] - 1] == 1:
                    self.ghb_input_lst.append([c, 0, self.ibound_arr.shape[-1] - 1, sea_lvl, self.cond_val[c][0][self.ibound_arr.shape[-1] - 1]])
                    #self.ghb_input_lst.append([c, 0, self.ibound_arr.shape[-1] - 1, self.top_elev[-1], self.cond_val[c][0][self.ibound_arr.shape[-1] - 1]])
                    self.ghb_arr[c, 0, self.ibound_arr.shape[-1] - 1] = -1
                    self.ssmdata.append([c, 0, self.ibound_arr.shape[-1] - 1, 35.0, itype['GHB']])      
                    
        if offshore_last_col_ghb_fresh:
            #   add GHB cells in the last offshore model column 
            #for c in range(self.nlay):
            for c in range(self.ibound_arr.shape[0]):
                #if self.ibound_arr[c, 0, self.ncol - 1] == 1 and self.top_elev[-1] >= 0.0:
                #    self.ghb_input_lst.append([c, 0, self.ncol - 1, self.top_elev[-1], self.cond_val[c][0][self.ncol - 1]])
                #    self.ghb_arr[c, 0, self.ncol - 1] = -1
                #    self.ssmdata.append([c, 0, self.ncol - 1, 0.0, itype['GHB']])              
                if self.ibound_arr[c, 0, self.ibound_arr.shape[-1] - 1] == 1 and self.top_elev[-1] >= 0.0:
                    self.ghb_input_lst.append([c, 0, self.ibound_arr.shape[-1] - 1, self.top_elev[-1], self.cond_val[c][0][self.ibound_arr.shape[-1] - 1]])
                    self.ghb_arr[c, 0, self.ibound_arr.shape[-1] - 1] = -1
                    self.ssmdata.append([c, 0, self.ibound_arr.shape[-1] - 1, 0.0, itype['GHB']])                         
        #   write the final output dictionary, inlcude each stress period
        self.ghb_arr_in = {}
        for d in range(len(self.perlen)):
            self.ghb_arr_in[d] = self.ghb_input_lst

    #   Define the input for the CHD package
    def chd_input(self, sea_lvl, offshore_last_col_chd = False):
        #   define and initiate the output arrays and lists
        self.chd_arr = self.ibound_arr * 1
        self.chd_input_lst = []
        if not self.ssmdata:
            self.ssmdata = []
        #  create the SSM dictionary where the ssm input will be written to
        itype = flopy.mt3d.Mt3dSsm.itype_dict()

        #   assign the sea water boundary at the last column cells, only in active cells
        for a in range(1, self.ncol):
            #   get the index of the first non-zero ibound cell
            try:
                top_cell_lay_idx = self.ibound_arr[:, 0, a].tolist().index(1)
                topo_val = self.top_elev[a]
                if topo_val < sea_lvl:
                    self.chd_input_lst.append([top_cell_lay_idx, 0, a, sea_lvl, sea_lvl])   #   sea level is 0.0m (present condition)
                    self.chd_arr[top_cell_lay_idx, 0, a] = -1
                    self.ssmdata.append([top_cell_lay_idx, 0, a, 35.0, itype['CHD']])
            except ValueError:
                pass

        if offshore_last_col_chd:
            #   add GHB cells in the last offshore model column 
            #for c in range(self.nlay):
            for c in range(self.ibound_arr.shape[0]):
                if self.ibound_arr[c, 0, self.ncol - 1] == 1:
                    self.chd_input_lst.append([c, 0, self.ncol - 1, sea_lvl, sea_lvl])
                    self.chd_arr[c, 0, self.ncol - 1] = -1
                    self.ssmdata.append([c, 0, self.ncol - 1, 35.0, itype['CHD']])            

        #   write the final output dictionary, inlcude each stress period
        self.chd_arr_in = {}
        for d in range(len(self.perlen)):
            self.chd_arr_in[d] = self.chd_input_lst    


    #   Define the input for the RCH package
    def rch_input(self, rch_rate, rch_type, sea_level):

        self.nrchop = 3  #  3 - rchrge applied to the highest active cell in each column
        self.rch_arr = np.array([[0.0] * 1 * self.ncol], dtype = np.float32)
        self.irch_arr = np.zeros((1, 1, self.ncol), dtype=np.float)

        if rch_rate is not None:
            self.rch_rate = rch_rate
            #   loop through the top_elev_lst and assign the precipitation value to cells above sea level
            for a in range(len(self.top_elev)):
                if self.top_elev[a] >= 0:
                    self.rch_arr[0][a] = rch_rate
                    self.irch_arr[0][0][a] = rch_rate

        else:
            self.rch_rate = None
            if rch_type == 'pcr':
                self.rch_lst = self.pcr_rch
            elif rch_type == 'watergap':
                self.rch_lst = self.watergap_rch
            elif rch_type == 'p_min_et':
                self.rch_lst = self.p_min_et
                
            #   repeat each list value multiple times to match the dimensions of the model discretization
            self.rch_lst = list(itertools.chain.from_iterable(itertools.repeat(x, int(500./self.del_col)) for x in self.rch_lst))
            self.rch_lst = self.rch_lst[:len(self.top_elev)]            
            
            #   loop through the top_elev_lst and assign the precipitation value to cells above sea level
            for a in range(len(self.top_elev)):
                if self.top_elev[a] >= sea_level:
                    #print a, self.rch_lst[a], self.rch_lst[a] * 0.001
                    self.rch_arr[0][a] = self.rch_lst[a] * 0.001
                    self.irch_arr[0][0][a] = self.rch_lst[a] * 0.001

    #   Define the input for the RCH package
    def rch_input_from_lst(self, rch_lst, sea_level):

        self.nrchop = 3  #  3 - rchrge applied to the highest active cell in each column
        #self.rch_arr = np.array([[0.0] * 1 * self.ncol], dtype = np.float32)
        #self.irch_arr = np.zeros((1, 1, self.ncol), dtype=np.float)
        self.rch_arr = np.array([[0.0] * 1 * self.ibound_arr.shape[-1]], dtype = np.float32)
        self.irch_arr = np.zeros((1, 1, self.ibound_arr.shape[-1]), dtype=np.float)
        
        #   loop through the top_elev_lst and assign the precipitation value to cells above sea level
        for a in range(self.ibound_arr.shape[-1]):
            #print(a)
            if self.top_elev[a] >= sea_level:
                #print a, self.rch_lst[a], self.rch_lst[a] * 0.001
                self.rch_arr[0][a] = rch_lst[a]# * 0.001
                self.irch_arr[0][0][a] = rch_lst[a]# * 0.001

                self.rch_extent = a

    #   Define the input for the RCH package
    def rch_input_from_lst_vk_limit(self, rch_lst, sea_level):

        self.nrchop = 3  #  3 - rchrge applied to the highest active cell in each column
        #self.rch_arr = np.array([[0.0] * 1 * self.ncol], dtype = np.float32)
        #self.irch_arr = np.zeros((1, 1, self.ncol), dtype=np.float)
        self.rch_arr = np.array([[0.0] * 1 * self.ibound_arr.shape[-1]], dtype = np.float32)
        self.irch_arr = np.zeros((1, 1, self.ibound_arr.shape[-1]), dtype=np.float)
        
        #   loop through the top_elev_lst and assign the precipitation value to cells above sea level
        #for a in range(self.ibound_arr.shape[-1]):
        for a in range(len(self.top_elev)):
            ibound_act_lay_idxs = [k for k, x in enumerate(self.ibound_arr[:, 0, a].tolist()) if x == 1]    
            #print(a)
            if self.top_elev[a] >= sea_level:
                #print a, self.rch_lst[a], self.rch_lst[a] * 0.001
                #if a < 50:
                #    print(self.vk_arr[ibound_act_lay_idxs[0], 0, a], rch_lst[a], min(self.vk_arr[ibound_act_lay_idxs[0], 0, a], rch_lst[a]))
                try:
                    self.rch_arr[0][a] = min(self.vk_arr[ibound_act_lay_idxs[0], 0, a], rch_lst[a])# * 0.001
                    self.irch_arr[0][0][a] = min(self.vk_arr[ibound_act_lay_idxs[0], 0, a], rch_lst[a])# * 0.001
                    self.rch_extent = a
                except IndexError:
                    self.rch_arr[0][a] = rch_lst[a]# * 0.001
                    self.irch_arr[0][0][a] =rch_lst[a]# * 0.001
                    self.rch_extent = a                    
                
    #   Define the input for the DRN package
    def drn_input(self, sea_level, dpth_bel_surf):
        #self.drn_factor = self.rch_rate * 10             #   why 10? is it ok?

        #  drainage is assigned only to cells that receive recharge - cells with elev above sea level
        drn_input_lst = []
        #for i in range(self.ncol):
        #for i in range(self.ibound_arr.shape[-1]):
        for i in range(len(self.top_elev)):
            ibound_col_lst = self.ibound_arr[:, 0, i].tolist()
            #   check the 1st column with ibound_val = 1 (active cell)
            try:
                drn_lay = ibound_col_lst.index(1)
                ##  now check if the elevation is below sea level, if so assign the cell to ghb list
                if self.top_elev[i] >= sea_level:
                    conc_cell = ((self.hk_arr[drn_lay, 0, i] / 10.) * self.del_col) / self.del_lay
                    drn_input_lst.append([drn_lay, 0, i, self.top_elev[i] - dpth_bel_surf, max(1., conc_cell)])
                    #drn_input_lst.append([drn_lay, 0, i, self.top_elev[i] - dpth_bel_surf, self.cond_val[drn_lay][0][i]])
                else:
                    pass
            except ValueError:
                pass

        #   write the final output dictionary, inlcude each stress period
        self.drn_arr_in = {}
        for c in range(len(self.perlen)):
            self.drn_arr_in[c] = drn_input_lst

    #   Define the input for the DRN package
    def drn_input_wtd(self, sea_level):
        #if self.rch_rate is not None:
        #    self.drn_factor = self.rch_rate # * 10             #   why 10? is it ok?
        #else:
        #    self.drn_factor = self.rch_lst

        #  drainage is assigned only to cells that receive recharge - cells with elev above sea level
        drn_input_lst = []
        #for i in range(self.ncol):
        #for i in range(self.ibound_arr.shape[-1]):
        for i in range(len(self.top_elev)):
            ibound_col_lst = self.ibound_arr[:, 0, i].tolist()
            #   check the 1st column with ibound_val = 1 (active cell)
            try:
                drn_lay = ibound_col_lst.index(1)
                ##  now check if the elevation is below sea level, if so assign the cell to ghb list
                if self.top_elev[i] >= sea_level:
                    if self.wtd_elev[i] > 0:
                        drn_input_lst.append([drn_lay, 0, i, self.top_elev[i] - self.wtd_elev[i], self.cond_val[drn_lay][0][i]])
                    else:
                        drn_input_lst.append([drn_lay, 0, i, self.top_elev[i], self.cond_val[drn_lay][0][i]])
                else:
                    pass
            except ValueError:
                #drn_input_lst.append([drn_lay, 0, i, self.top_elev[i] - 0.5, self.cond_val[drn_lay][0][i]])
                pass

        #   write the final output dictionary, inlcude each stress period
        self.drn_arr_in = {}
        for c in range(len(self.perlen)):
            self.drn_arr_in[c] = drn_input_lst


    #   Define the input for the RIV package
    def riv_input(self, m_factor = '05'):
        #   select the right m_factor conductance set
        if m_factor == '01':
            self.riv_cond = self.riv_cond_01
        elif m_factor == '05':
            self.riv_cond = self.riv_cond_05
        elif m_factor == '10':
            self.riv_cond = self.riv_cond_10

        #   save the cell by cell budget, 0 is No, anything else is yes
        riv_input_lst = []

        for i in range(len(self.riv_cond)):
            if math.isnan(self.riv_cond[i]):
                pass
            else:
                ibound_col_lst = self.ibound_arr[:, 0, i].tolist()
                try:
                    riv_lay = ibound_col_lst.index(1)

                    riv_input_lst.append([riv_lay, 0, i, self.riv_head_elev[i], self.riv_cond[i], self.riv_bot_elev[i]])
                    print (i, self.riv_cond[i])
                except ValueError:   # ValueError: 1 is not in list
                    pass

        #   write the final output dictionary, inlcude each stress period
        self.riv_arr_in = {}
        for c in range(len(self.perlen)):
            self.riv_arr_in[c] = riv_input_lst

    #   Define the input for the WEL package
    def wel_input(self, well_def, str_per):
        #   the well definition has to be in the following format [[[layer, row, column, pumping rate], [...],...]]
        #   the total amount of lists in the list has to correspond to the number of stress periods that have any pumping activity
        #   define the dictionary with all the wells in the system, this is the output
        self.lrcq = {}
        for i in range(len(str_per)):
            self.lrcq[str_per[i]] = well_def[i]

    #   Define the input for the PCG package
    def pcg_input(self, hclose_val, rclose_val):
        self.hclose_val = hclose_val
        self.rclose_val = rclose_val

    #   Define the input for the OC package
    def oc_input(self, th_time_step):

        self.ihedfm = 1          # a code for the format in which heads will be printed.
        self.iddnfm = 0          # a code for the format in which drawdowns will be printed.
        self.extension = ['oc','hds','ddn','cbc']
        self.unitnumber = [14, 30, 0, 50]
        #   create the dictionary that defines how to write the output file
        self.spd = {(0, 0): ['save head', 'save budget']}
        for t in range(0, self.nper):
            per = t #+ 1
            #   xrange allows to iterate through the list with specified step size - 25
            #   to save space on disk, every 10th timestep is saved
            for g in range(0, self.nstp[t] + 1, th_time_step):
                self.spd[(per, int(g))] = ['save head', 'save budget']
                self.spd[(per, int(g) + 1)] = []

            self.spd[(per, int(g) + 1)] = ['save head', 'save budget']
            self.spd[(per, int(g) - 1)] = ['save head', 'save budget']

    #   Define the input for the OC package
    def oc_input_start_end(self):

        self.ihedfm = 1          # a code for the format in which heads will be printed.
        self.iddnfm = 0          # a code for the format in which drawdowns will be printed.
        self.extension = ['oc','hds','ddn','cbc']
        self.unitnumber = [14, 30, 0, 50]
        #   create the dictionary that defines how to write the output file
        self.spd = {(0, 0): ['save head', 'save budget']}
        self.spd[(0, 0)] = ['save head', 'save budget'] 
        self.spd[(0, self.nstp[0])] = ['save head', 'save budget']  
        self.spd[(0, self.nstp[0])] = ['save head', 'save budget']        

    #   Define the input for the BTN package
    def btn_input(self, prsity_val, dt0_val, th_time_step, ifmtcn_val, nprs_val, chkmas_val, nprobs_val, nprmas_val):
        #   the nprs parameter defines how many transport steps are going to be exported to the UCN file
        timprs_sp1 = np.linspace(1., self.perlen[0], th_time_step[0], endpoint = True)
        if self.nper > 1:
            timprs_sp2 = np.linspace(self.perlen[0], self.perlen[0] + self.perlen[1], th_time_step[1], endpoint = True)
            self.timprs = np.concatenate((timprs_sp1, timprs_sp2[1:]), axis = 0)
        else:
            self.timprs = timprs_sp1
        self.dt0 = dt0_val
        self.porosity = prsity_val
        self.ifmtcn = ifmtcn_val
        self.nprs = nprs_val
        self.chkmas = chkmas_val
        self.nprobs = nprobs_val
        self.nprmas = nprmas_val

    #   Define the input for the ADV package
    def adv_input(self, mixelm_val):
        self.mixelm = mixelm_val

    #   Define the input for the DSP package
    def dsp_input(self, al_val, trpt_val, trpv_val, dmcoef_val):
        self.al = al_val
        self.trpt = trpt_val
        self.trpv = trpv_val
        self.dmcoef = dmcoef_val

    #   Define the input for the GCG package
    def gcg_input(self, iter1_val, mxiter_val, isolve_val, cclose_val):
        self.iter1 = iter1_val
        self.mxiter = mxiter_val
        self.isolve = isolve_val
        self.cclose = cclose_val

    #   Define the input for the VDF package
    def vdf_input(self, iwtable_val, densemin_val, densemax_val, denseref_val, denseslp_val, firstdt_val):
        self.iwtable = iwtable_val
        self.densemin = densemin_val
        self.densemax = densemax_val
        self.denseref = denseref_val
        self.denseslp = denseslp_val
        self.firstdt = firstdt_val

    #   Define the input for the SSM package
    def ssm_input(self):
        #   define the two SSM input lists/dictionary
        self.ssm_rch_in = np.array([[0.0] * 1 * self.ncol], dtype = np.float32)
        self.ssmdata_dict = {0: self.ssmdata,\
                             1: self.ssmdata}

    #   write all the MODFLOW/MT3D/SEAWAT packages
    def write_all_input(self, modelname, mf_exe_dir, mt3d_exe_dir, swat_exe_dir, folder_name):
        self.modelname = modelname
        self.foldername = folder_name
        #   1) create the models via flopy
        #self.mf = flopy.modflow.Modflow(modelname, exe_name = mf_exe_dir, model_ws = folder_name)
        #self.mt = flopy.mt3d.Mt3dms(modelname, 'nam_mt3dms', self.mf, model_ws = folder_name, exe_name = mt3d_exe_dir)
        self.mswt = flopy.seawat.Seawat(modelname, 'nam_swt',  model_ws = folder_name, exe_name = swat_exe_dir)
        #self.mswt = flopy.seawat.Seawat(modelname, 'nam_swt', model_ws = folder_name, exe_name = swat_exe_dir)

        """MF"""

        self.nlay = self.ibound_arr.shape[0]
        self.botm = self.botm[:self.nlay]
        self.hk_arr = self.hk_arr[:self.nlay, :, :]
        self.vk_arr = self.vk_arr[:self.nlay, :, :]

        #   2) write the DIS (discretization) package
        self.dis = flopy.modflow.ModflowDis(self.mswt, self.nlay, self.nrow, self.ncol, self.nper, self.delr, self.delc,\
                                            self.laycbd, self.top, self.botm, self.perlen, self.nstp)
        #   3) write the BAS package
        self.bas = flopy.modflow.ModflowBas(self.mswt, self.ibound_arr, self.strt_arr)
        #   4) write the LPF package
        self.lpf = flopy.modflow.ModflowLpf(self.mswt, laytyp = self.laytyp, hk = self.hk_arr, vka = self.vk_arr,  ipakcb = 1)
        #   5) write the GHB package
        self.ghb = flopy.modflow.ModflowGhb(self.mswt, ipakcb = 1, stress_period_data = self.ghb_arr_in)  #   ipakcb - write output in cbc file
        #   6) write the RCH package
        self.rch = flopy.modflow.ModflowRch(self.mswt, nrchop = self.nrchop, ipakcb = 1, rech = self.rch_arr, irch = self.irch_arr)
        #   7) write the DRN package
        if len(self.drn_arr_in[0]) != 0:
            self.drn = flopy.modflow.ModflowDrn(self.mswt, ipakcb = 1, stress_period_data = self.drn_arr_in)
        #   8) write the WEL package
        #self.wel = flopy.modflow.ModflowWel(self.mf, stress_period_data = self.lrcq)
        #   9) write the PCG package
        self.pcg = flopy.modflow.ModflowPcg(self.mswt, hclose = self.hclose_val, rclose = self.rclose_val)
        #   10) write the OC package
        self.oc = flopy.modflow.ModflowOc(self.mswt, stress_period_data = self.spd, compact = True)
        #   11) write the RIV package
        #self.riv = flopy.modflow.ModflowRiv(self.mswt, ipakcb = 1, stress_period_data = self.riv_arr_in)
        #   12) write the CHD package
        #self.chd = flopy.modflow.ModflowChd(self.mswt, stress_period_data = self.chd_arr_in)  #   ipakcb - write output in cbc file
        
        """MT"""

        #   11) write the BTN package
        self.btn = flopy.mt3d.Mt3dBtn(self.mswt, nprs = self.nprs, timprs = self.timprs, prsity = self.porosity, sconc = self.sconc_arr, ifmtcn = self.ifmtcn,
                             chkmas = self.chkmas, nprobs = self.nprobs, nprmas = self.nprmas, dt0 = self.dt0)
        #   12) write the ADV package
        self.adv = flopy.mt3d.Mt3dAdv(self.mswt, mixelm = self.mixelm, mxpart = 2000000)
        #   13) write the ADV package
        self.dsp = flopy.mt3d.Mt3dDsp(self.mswt, al = self.al, trpt = self.trpt, trpv = self.trpv, dmcoef = self.dmcoef)
        #   14) write the ADV package
        self.gcg = flopy.mt3d.Mt3dGcg(self.mswt, iter1 = self.iter1, mxiter = self.mxiter, isolve = self.isolve, cclose = self.cclose)
        #   15) write the ADV package
        self.vdf = flopy.seawat.SeawatVdf(self.mswt, iwtable = self.iwtable, densemin = self.densemin, densemax = self.densemax,\
                                          denseref = self.denseref, denseslp = self.denseslp, firstdt = self.firstdt)
        #   16) write the SSM package
        self.ssm = flopy.mt3d.Mt3dSsm(self.mswt, crch = self.ssm_rch_in, stress_period_data = self.ssmdata_dict)
        #   17) write all the input files
        #self.mf.write_input()
        #self.mt.write_input()
        self.mswt.write_input()


    #   write all the MODFLOW/MT3D/SEAWAT packages
    def write_all_input_no_lpf(self, modelname, mf_exe_dir, mt3d_exe_dir, swat_exe_dir, folder_name):
        self.modelname = modelname
        self.foldername = folder_name
        #   1) create the models via flopy
        self.mf = flopy.modflow.Modflow(modelname, exe_name = mf_exe_dir, model_ws = folder_name)
        self.mt = flopy.mt3d.Mt3dms(modelname, 'nam_mt3dms', self.mf, model_ws = folder_name, exe_name = mt3d_exe_dir)
        self.mswt = flopy.seawat.Seawat(modelname, 'nam_swt', self.mf, self.mt, model_ws = folder_name, exe_name = swat_exe_dir)
        #self.mswt = flopy.seawat.Seawat(modelname, 'nam_swt', model_ws = folder_name, exe_name = swat_exe_dir)

        """MF"""

        #   2) write the DIS (discretization) package
        self.dis = flopy.modflow.ModflowDis(self.mf, self.nlay, self.nrow, self.ncol, self.nper,self.delr, self.delc,\
                                            self.laycbd, self.top, self.botm, self.perlen, self.nstp)
        #   3) write the BAS package
        self.bas = flopy.modflow.ModflowBas(self.mf, self.ibound_arr, self.strt_arr)
        #   4) write the LPF package
        #self.lpf = flopy.modflow.ModflowLpf(self.mf, laytyp = self.laytyp, hk = self.hk_arr, vka = self.vk_arr)
        #   5) write the GHB package
        self.ghb = flopy.modflow.ModflowGhb(self.mf, ipakcb = 1, stress_period_data = self.ghb_arr_in)  #   ipakcb - write output in cbc file
        #   6) write the RCH package
        self.rch = flopy.modflow.ModflowRch(self.mf, nrchop = self.nrchop, ipakcb = 1, rech = self.rch_arr, irch = self.irch_arr)
        #   7) write the DRN package
        self.drn = flopy.modflow.ModflowDrn(self.mf, ipakcb = 1, stress_period_data = self.drn_arr_in)
        #   8) write the WEL package
        #self.wel = flopy.modflow.ModflowWel(self.mf, stress_period_data = self.lrcq)
        #   9) write the PCG package
        self.pcg = flopy.modflow.ModflowPcg(self.mf, hclose = self.hclose_val, rclose = self.rclose_val)
        #   10) write the OC package
        self.oc = flopy.modflow.ModflowOc(self.mf, self.ihedfm, self.iddnfm, cboufm='(20i5)', stress_period_data = self.spd)
        #   11) write the RIV package
        self.riv = flopy.modflow.ModflowRiv(self.mf, ipakcb = 1, stress_period_data = self.riv_arr_in)

        """MT"""

        #   11) write the BTN package
        self.btn = flopy.mt3d.Mt3dBtn(self.mt, nprs = self.nprs, timprs = self.timprs, prsity = self.porosity, sconc = self.sconc_arr, ifmtcn = self.ifmtcn,
                             chkmas = self.chkmas, nprobs = self.nprobs, nprmas = self.nprmas, dt0 = self.dt0)
        #   12) write the ADV package
        self.adv = flopy.mt3d.Mt3dAdv(self.mt, mixelm = self.mixelm, mxpart = 2000000)
        #   13) write the ADV package
        self.dsp = flopy.mt3d.Mt3dDsp(self.mt, al = self.al, trpt = self.trpt, trpv = self.trpv, dmcoef = self.dmcoef)
        #   14) write the ADV package
        self.gcg = flopy.mt3d.Mt3dGcg(self.mt, iter1 = self.iter1, mxiter = self.mxiter, isolve = self.isolve, cclose = self.cclose)
        #   15) write the ADV package
        self.vdf = flopy.seawat.SeawatVdf(self.mt, iwtable = self.iwtable, densemin = self.densemin, densemax = self.densemax,\
                                          denseref = self.denseref, denseslp = self.denseslp, firstdt = self.firstdt)
        #   16) write the SSM package
        self.ssm = flopy.mt3d.Mt3dSsm(self.mt, crch = self.ssm_rch_in, stress_period_data = self.ssmdata_dict)
        #   17) write all the input files
        self.mf.write_input()
        self.mt.write_input()
        self.mswt.write_input()


    #   Run the model and measure the run time
    def run_model(self):
        t0 = time()
        #   run the model
        v = self.mswt.run_model(silent = False, report = True)
        for idx in range(-3, 0):
            print(v[1][idx])
        #   stop measuring time and calculate total run time
        t1 = time()
        self.run_time = t1 - t0

    #   Read model output; heads, concentrations and cell budget flow
    def read_output(self):
        #self.ml_results = flopy.modflow.Modflow.load(self.modelname + ".nam", model_ws = self.foldername, verbose = False, check = False, exe_name = "mfnwt")
        self.hdsobj = bf.HeadFile(os.path.join(self.foldername, self.modelname + '.hds'))#, model = self.ml_results)#, precision = 'double')
        self.h = self.hdsobj.get_alldata()
        self.ucnobj = bf.UcnFile(os.path.join(self.foldername, 'MT3D001.UCN'))#, precision = 'double')
        self.time_steps = self.ucnobj.get_times()
        self.cbbobj = bf.CellBudgetFile(os.path.join(self.foldername, self.modelname + '.cbc'))#, precision = 'double')
        self.times_heads = self.cbbobj.get_times()

    #   Check that the heads in GHB cells are equal to the head elevation specified for them
    def check_ghb_head(self, time_step, tolerance):
        succ_cnt, fail_cnt, max_diff = 0, 0, 0
        self.ghb_head_diff_lst = []
        #   loop through the GHB cells and check if the condition is met, calculate how many cells are ok and how many not
        for ghb_cell in self.ghb_input_lst:
            #   get infor for the given cell and compare it to the Heads output array
            ghb_cell_lay, ghb_cell_row, ghb_cell_col = ghb_cell[0], ghb_cell[1], ghb_cell[2]
            ghb_cell_head_elev = ghb_cell[3]
            cell_head_val = self.h[time_step][ghb_cell_lay, ghb_cell_row, ghb_cell_col]
            #   calculate the head difference
            head_diff = ghb_cell_head_elev - cell_head_val # abs(cell_head_val) - abs(ghb_cell_head_elev)
            self.ghb_head_diff_lst.append([ghb_cell_lay, ghb_cell_row, ghb_cell_col, head_diff])
            #   check if the criteria is met
            #if head_diff < tolerance:
            if abs(head_diff) < tolerance:
                succ_cnt += 1
            else:
                fail_cnt += 1
            #if head_diff > max_diff:
            if abs(head_diff) > max_diff:
                max_diff = head_diff
            else:
                pass
        #   finally based on the sizes of both the success and failure lists decide if the check is succesful
        if fail_cnt == 0:
            self.ghb_head_check = True
            self.ghb_head_check_max_diff = max_diff
        else:
            self.ghb_head_check = False
            self.ghb_head_check_max_diff = max_diff
        self.ghb_succ, self.ghb_fail = succ_cnt, fail_cnt

    #   function that creates and saves a plot of total mass in the model domain and of changes of the amount of
    #   fresh/brackish and seawater concentration cells in the whole model domain.
    def model_output_fresh_volume_graphs(self):
        #   create empty lists
        fresh_cnt, brackish_cnt, salt_cnt, tot_mass = [], [], [], []
        fresh_comparison, brackish_comparison, salt_comparison = [], [], []
        #   get all the time steps from the UCN file
        ts_yrs = [round(i / 365.25, 0) for i in self.time_steps]
        #   get the count of fresh/brackish and salt concentration cells at the first time step
        c = self.ucnobj.get_data(totim = self.time_steps[0])
        c[c > 35.0] = -1        #   all concentrations above 35.0 to no value
        fresh_cnt_ref = ((c > -1.) & (c <= 0.05)).sum()
        brackish_cnt_ref = ((c > 0.05) & (c <= 15.0)).sum()
        salt_cnt_ref = ((c > 15.0) & (c <= 35.0)).sum()
        fresh_comparison.append(0)
        brackish_comparison.append(0)
        salt_comparison.append(0)
        fresh_cnt.append(fresh_cnt_ref)
        brackish_cnt.append(brackish_cnt_ref)
        salt_cnt.append(salt_cnt_ref)
        #   calculate the total mass in the system
        c[c == -1] = 0.0        #   all negative values to 0.0 so it doesnt count negative concentrations..
        total_mass = (c * (self.del_lay * self.del_col)).sum()
        tot_mass.append(total_mass)

        #   loop through all the time steps
        for i in range(1, len(self.time_steps)):
            #   read in the concentration profile
            c = self.ucnobj.get_data(totim = self.time_steps[i])
            c[c > 35.0] = -1        #   all concentrations above 35.0 to no value
            #   get the number of cells for the given time step
            fresh = ((c > -1.) & (c <= 0.05)).sum()
            brackish = ((c > 0.05) & (c <= 15.0)).sum()
            salt = ((c > 15.0) & (c <= 35.0)).sum()
            #   calculate the total mass in the system
            c[c == -1] = 0.0        #   all negative values to 0.0 so it doesnt count negative concentrations..
            total_mass = (c * (self.del_lay * self.del_col)).sum()
            #   append all the output to the lists
            fresh_cnt.append(fresh)
            brackish_cnt.append(brackish)
            salt_cnt.append(salt)
            tot_mass.append(total_mass)
            fresh_comparison.append(fresh_cnt_ref - fresh)
            brackish_comparison.append(brackish_cnt_ref - brackish)
            salt_comparison.append(salt_cnt_ref - salt)

        #   create the output plots
        f = plt.figure(figsize=(20, 10))
        #   first the plot with fresh/brackish/salt cells count
        ax1 = f.add_subplot(3, 1, 1)
        f.add_subplot()
        ax1.set_title('Change in number of fresh, brackish and sea water concentration cells in the model domain')
        ax1.plot(ts_yrs, fresh_cnt, color = 'b')
        ax1.plot(ts_yrs, brackish_cnt, color = 'g')
        ax1.plot(ts_yrs, salt_cnt, color = 'r')
        #   also add the graph with differences in cell count
        ax2 = f.add_subplot(3, 1, 2)
        f.add_subplot()
        ax2.set_title('Changes in cell count compared to 1st time step (1st ts - nth ts)')
        ax2.plot(ts_yrs, fresh_comparison, color = 'b')
        ax2.plot(ts_yrs, brackish_comparison, color = 'g')
        ax2.plot(ts_yrs, salt_comparison, color = 'r')

        #   also add the total mass graph in the system
        ax3 = f.add_subplot(3, 1, 3)
        f.add_subplot()
        ax3.set_title('Total mass in the model domain')
        ax3.plot(ts_yrs, tot_mass, color = 'black')
        #   save the output
        figname = 'model_output_graphs.png'
        plt.savefig(os.path.join(self.foldername, figname))
        plt.close()
        del f

    #   Plotting of 2D concentration profiles at given time
    def plot_conc_profiles(self, time_yrs_lst):
        """
        #   transform the pro
        plot_profile_days = []
        plot_profile_days = [x * 365.25 for x in time_yrs_lst]
        #   define the colorbar
        cmap = plt.cm.jet
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
        # define the bins and normalize
        bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        #   create the mask for plotting the aquitard layers
        #   more info - http://stackoverflow.com/questions/24539296/outline-a-region-in-a-graph
        mapimg = (self.hk_arr == 0.01123)
        ver_seg = np.where(mapimg[:,1:] != mapimg[:,:-1])
        hor_seg = np.where(mapimg[1:,:] != mapimg[:-1,:])
        l = []
        for p in zip(*hor_seg):
            l.append((p[1], p[0]+1))
            l.append((p[1]+1, p[0]+1))
            l.append((np.nan,np.nan))
        for p in zip(*ver_seg):
            l.append((p[1]+1, p[0]))
            l.append((p[1]+1, p[0]+1))
            l.append((np.nan, np.nan))
        segments = np.array(l)
        segments[:,0] = self.x_start + (self.x_end - self.x_start) * segments[:,0] / mapimg.shape[1]
        segments[:,1] = min(self.top_elev) + (max(self.top_elev) - min(self.top_elev)) * segments[:,1] / mapimg.shape[0]
        #   loop through the list created above and plot results for each time specified
        for z in xrange(len(plot_profile_days)):
            #   first find the index of the corresponding time step in the time list
            ts_conc_indx = min(range(len(self.time_steps)), key=lambda i: abs(self.time_steps[i] - plot_profile_days[z]))
            #   read the concentration for the time step
            concentration = self.ucnobj.get_data(totim = self.time_steps[ts_conc_indx])
            concentration[np.abs(concentration) >= 999.] = np.nan
            #   get the exact time since the start of the model
            time_title = round(self.time_steps[ts_conc_indx] / 365.25, 0)
            print(ts_conc_indx), str(int(time_title))
            # Concentration and Flow Plot
            fig1 = plt.figure(figsize=(20, 5))
            ax1 = fig1.add_subplot(1, 1, 1)
            ax1.axhline(y = 0., linewidth = 2, color = 'k')
            ax1.axvline(x = 0., linewidth = 2, color = 'k')
            print self.x_start, self.x_end, self.zbot, math.ceil((self.top / 100.0) * 100.0)
            im1 = ax1.imshow(concentration[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                           extent = (self.x_start, self.x_end, self.zbot, math.ceil((self.top / 100.0) * 100.0)), vmin = 0, vmax = 35)
            ax1.plot(segments[:,0], segments[:,1], color = (1,0,0,.5), linewidth = 3)
            y, x, z = self.dis.get_node_coordinates()
            x = np.linspace(self.x_start, self.x_end, self.end_idx)
            X, Z = np.meshgrid(x, z[:, 0, 0])
            ax1.set_title('Concentration at ' + str(int(time_title)) + ' years')
            cbaxes = fig1.add_axes([0.04, 0.1, 0.03, 0.8])
            cbar = plt.colorbar(im1, cax = cbaxes, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)#, format = '%1f')
            cbar.set_label('Concentration of salt g/l', rotation = 270, fontsize = 9, labelpad= 15)
            cbar.ax.tick_params(labelsize = 8)
            figname = '_conc_' + str(int(time_title)) + '_yrs.png'
            plt.savefig(os.path.join(self.foldername, figname))
            plt.close()
            del fig1
        """
        #   calculate the offset of the coastline for plotting - this happens because the first offset is calculated from the
        #   raw topo data (from database with 0.5km spacing). When creating the top_elev interpolation the coast can be shifted by a bit
        cst_offset_val  = next((x for x in self.top_elev if x < 0.0), None)
        cst_offset_idx = self.top_elev.index(cst_offset_val) - 1
        self.cst_offset_plot = round(self.x_start + cst_offset_idx / (1000. /self.del_col), 2) #  (1000. / model.del_col) - to get the distance in km
        #   transform the pro
        plot_profile_days = []
        plot_profile_days = [x * 365.25 for x in time_yrs_lst]
        #   define the colorbar
        cmap = plt.cm.jet
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
        # define the bins and normalize
        bounds = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 35.0]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        #   create the mask for plotting the aquitard layers
        #   more info - http://stackoverflow.com/questions/24539296/outline-a-region-in-a-graph
        mapimg = (self.hk_arr == 0.01123)

        plt.figure()
        plt.imshow(mapimg[:,0,:])

        ver_seg = np.where(mapimg[:,0,1:] != mapimg[:,0,:-1])   #   this is empty for some reason!
        hor_seg = np.where(mapimg[1:,0,:] != mapimg[:-1,0,:])
        l = []
        for p in zip(*hor_seg):
            l.append((p[1], p[0]+1))
            l.append((p[1]+1, p[0]+1))
            l.append((np.nan,np.nan))
        for p in zip(*ver_seg):
            l.append((p[1]+1, p[0]))
            l.append((p[1]+1, p[0]+1))
            l.append((np.nan, np.nan))
        segments = np.array(l)
        #segments[:,0] = self.x_start - self.cst_offset_plot + (self.x_end - self.x_start - self.cst_offset_plot) * segments[:,0] / (mapimg.shape[-1])
        segments[:,0] = self.x_start + (self.x_end - self.x_start) * segments[:,0] / (mapimg.shape[-1])        
        segments[:,1] = min(self.top_elev) + (max(self.top_elev) - min(self.top_elev) + self.del_lay) * -1 * (segments[:,1] / mapimg.shape[0]) + (max(self.top_elev) - min(self.top_elev) + self.del_lay)  #


        #   loop through the list created above and plot results for each time specified
        for z in range(len(plot_profile_days)):
            #   first find the index of the corresponding time step in the time list
            ts_conc_indx = min(range(len(self.time_steps)), key=lambda i: abs(self.time_steps[i] - plot_profile_days[z]))
            #   read the concentration for the time step
            concentration = self.ucnobj.get_data(totim = self.time_steps[ts_conc_indx])
            concentration[np.abs(concentration) >= 999.] = np.nan
            #   get the exact time since the start of the self
            time_title = round(self.time_steps[ts_conc_indx] / 365.25, 0)
            print(ts_conc_indx), str(int(time_title))
            # Concentration and Flow Plot
            fig1 = plt.figure(figsize=(20, 5))
            ax1 = fig1.add_subplot(1, 1, 1)
            ax1.axhline(y = 0., linewidth = 2, color = 'k')
            ax1.axvline(x = 0., linewidth = 2, color = 'k')
            print (self.x_start, self.x_end, self.zbot, math.ceil((self.top / 100.0) * 100.0))
            im1 = ax1.imshow(concentration[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                           extent = (self.x_start, self.x_end, self.zbot, math.ceil((self.top / 100.0) * 100.0)), vmin = 0, vmax = 35)
                           #extent = (self.x_start - self.cst_offset_plot, self.x_end - self.cst_offset_plot, self.zbot, math.ceil((self.top / 100.0) * 100.0)), vmin = 0, vmax = 35)
            ax1.plot(segments[:,0], segments[:,1], color = 'white', linewidth = .5)
            y, x, z = self.dis.get_node_coordinates()
            x = np.linspace(self.x_start, self.x_end, self.end_idx)
            X, Z = np.meshgrid(x, z[:, 0, 0])
            ax1.set_title('Concentration at ' + str(int(time_title)) + ' years')
            cbaxes = fig1.add_axes([0.04, 0.1, 0.03, 0.8])
            cbar = plt.colorbar(im1, cax = cbaxes, cmap = cmap, norm = norm, spacing = 'uniform', ticks = bounds, boundaries = bounds)#, format = '%1f')
            cbar.set_label('Concentration of salt g/l', rotation = 270, fontsize = 9, labelpad= 15)
            cbar.ax.tick_params(labelsize = 8)
            figname = '_conc_' + str(int(time_title)) + '_yrs.png'
            plt.savefig(os.path.join(self.foldername, figname), dpi = 300, facecolor = 'w', edgecolor = 'w',
                        orientation = 'portrait', papertype = None, format = None,
                        transparent = False, bbox_inches = 'tight', pad_inches = 0.1, frameon = None)
            plt.close(fig1)
            #del fig1

    #   Plotting of 2D concentration profiles at given time
    def plot_flow_profiles(self, time_yrs_lst):
        #   transform the pro
        plot_profile_days = []
        plot_profile_days = [x * 365.25 for x in time_yrs_lst]
        #   loop through the list created above and plot results for each time specified
        for z in range(len(plot_profile_days)):
            #   first find the index of the corresponding time step in the time list
            ts_heads_indx = min(range(len(self.times_heads)), key=lambda i: abs(self.times_heads[i] - plot_profile_days[z]))
            #   read the heads for the time step
            heads = self.hdsobj.get_data(totim = self.times_heads[ts_heads_indx])
            heads[np.abs(heads) >= 999.] = np.nan
            #   get the exact time since the start of the model
            time_title = round(self.times_heads[ts_heads_indx] / 365.25, 0)
            print(ts_heads_indx), str(int(time_title))
            # Flow for the time step, prepare for flow arrows
            qx1 = self.cbbobj.get_data(text = 'flow right face', totim = self.times_heads[ts_heads_indx])[0]
            qz1 = self.cbbobj.get_data(text = 'flow lower face', totim = self.times_heads[ts_heads_indx])[0]
            qx1_avg = np.empty(qx1.shape, dtype = qx1.dtype)
            qz1_avg = np.empty(qz1.shape, dtype = qz1.dtype)
            qx1_avg[:, :, :] = 0.5 * (qx1[:, :, 0 : self.ncol] + qx1[:, :, : self.ncol])
            qx1_avg[:, :, 0] = 0.5 * qx1[:, :, 0]
            qz1_avg[:, :, :] = 0.5 * (qz1[0 : self.nlay, :, :] + qz1[: self.nlay, :, :])
            qz1_avg[0, :, :] = 0.5 * qz1[0, :, :]
            # Concentration and Flow Plot
            fig1 = plt.figure(figsize=(20, 5))
            ax1 = fig1.add_subplot(1, 1, 1)
            print (self.x_start, self.x_end, self.zbot, math.ceil((self.top / 100.0) * 100.0))
            im1 = ax1.imshow(heads[:, 0, :], aspect = 'auto', interpolation = 'nearest',
                           extent = (self.x_start, self.x_end, self.zbot, math.ceil((self.top / 100.0) * 100.0)), vmin = np.amin(heads), vmax = np.amax(heads))
                           #extent = (self.x_start - self.cst_offset_plot, self.x_end - self.cst_offset_plot, self.zbot, math.ceil((self.top / 100.0) * 100.0)), vmin = np.amin(heads), vmax = np.amax(heads))
            y, x, z = self.dis.get_node_coordinates()
            x = np.linspace(self.x_start, self.x_end, self.end_idx)
            X, Z = np.meshgrid(x, z[:, 0, 0])
            iskip = 3
            ax1.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip],
                       qx1_avg[::iskip, 0, ::iskip], -qz1_avg[::iskip, 0, ::iskip],
                       color = 'black', scale = 9, headwidth = 3, headlength = 2,
                       headaxislength=2, width = 0.0015)
            ax1.axhline(y = 0., linewidth = 2, color = 'k')
            ax1.axvline(x = 0., linewidth = 2, color = 'k')
            ax1.set_title('BlaBla Heads and flow pattern at ' + str(int(time_title)) + ' years')
            cbaxes = fig1.add_axes([0.04, 0.1, 0.03, 0.8])
            plt.colorbar(im1, cax = cbaxes)
            figname = '_heads_flow_' + str(int(time_title)) + '_yrs.png'
            plt.savefig(os.path.join(self.foldername, figname))
            plt.close(fig1)
            #del fig1
            
            
    #   plot the input parameter values and the HK array
    #   zoom = [start_x, end_x]     - zooms into the specified X domain            
    def plot_bound_topsys_input(self, out_dir, title, sea_lvl, zoom = False):
    
        #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
        fig = plt.figure(figsize = (20, 24))
        #   first specify the position of all the boundary and top system plots, the bottom one is the HK geology plot
        ax1 = plt.subplot2grid((7, 2), (0, 0))  #    
        ax2 = plt.subplot2grid((7, 2), (1, 0), sharex = ax1)  #   
        ax3 = plt.subplot2grid((7, 2), (2, 0), sharex = ax1)  #   
        ax4 = plt.subplot2grid((7, 2), (3, 0), sharex = ax1)  #  
        ax5 = plt.subplot2grid((7, 2), (4, 0), sharex = ax1)  #  
        ax6 = plt.subplot2grid((7, 2), (5, 0), sharex = ax1)  #   
        ax7 = plt.subplot2grid((7, 2), (6, 0), sharex = ax1)  #   
        #   this is the part with the colorbars and legends
        ax1b = plt.subplot2grid((7, 2), (0, 1))  #   
        ax2b = plt.subplot2grid((7, 2), (1, 1))  #   
        ax3b = plt.subplot2grid((7, 2), (2, 1))  #   
        ax4b = plt.subplot2grid((7, 2), (3, 1))  #  
        ax5b = plt.subplot2grid((7, 2), (4, 1))  # 
        ax6b = plt.subplot2grid((7, 2), (5, 1))  #  
        ax7b = plt.subplot2grid((7, 2), (6, 1))  # 
        #   set the exact positions for each axis
        ax1.set_position([0.05, 0.825, 0.825, 0.1]) # [left, bottom, width, height]
        ax2.set_position([0.05, 0.725, 0.825, 0.1])
        ax3.set_position([0.05, 0.625, 0.825, 0.1])        
        ax4.set_position([0.05, 0.525, 0.825, 0.1]) # [left, bottom, width, height]
        ax5.set_position([0.05, 0.425, 0.825, 0.1])
        ax6.set_position([0.05, 0.325, 0.825, 0.1])       
        ax7.set_position([0.05, 0.025, 0.825, 0.3]) # [left, bottom, width, height]
        ax1b.set_position([0.885, 0.825, 0.07, 0.1]) # [left, bottom, width, height]
        ax2b.set_position([0.885, 0.725, 0.07, 0.1])
        ax3b.set_position([0.885, 0.625, 0.07, 0.1])       
        ax4b.set_position([0.885, 0.525, 0.07, 0.1]) # [left, bottom, width, height]
        ax5b.set_position([0.885, 0.425, 0.07, 0.1])
        ax6b.set_position([0.885, 0.325, 0.07, 0.1])       
        ax7b.set_position([0.915, 0.025, 0.025, 0.3]) # [left, bottom, width, height]
        
        #   define the extent of the figures
        x_start_new = self.x_start   #   -(self.cst_idx - self.idx_start) / 2.
        x_end_new = self.x_end
        y_max_new = math.ceil((self.top / 100.0) * 100.0)
        y_min_new = min(self.bot_elev)    
        
        #   based on the zoom being True or False decide what to plot
        if zoom:
            if x_start_new > zoom[0]: 
                x_start_new_zoom = x_start_new
                y_max_new_zoom = math.ceil((self.top / 100.0) * 100.0)
            else:
                x_start_new_zoom = zoom[0]
                y_max_new_zoom = math.ceil((self.top_elev[int(abs(x_start_new_zoom - zoom[0]) * 10)] / 100.0) * 100.0)
            if x_end_new < zoom[1]:                
                x_end_new_zoom = x_end_new
                y_min_new_zoom = min(self.bot_elev)
            else:
                x_end_new_zoom = zoom[1]
                y_min_new_zoom = math.floor((self.bot_elev[int((self.cst_idx - self.idx_start) * 5 + (zoom[1] * (1000. / self.del_col)))] / 100.) * 100.)       
                
            x_start_plt, x_end_plt = x_start_new_zoom, x_end_new_zoom
            y_min_plt, y_max_plt = y_min_new_zoom, y_max_new_zoom 
            fig_out_dir = os.path.join(out_dir, '_bound_topsys_input_zoom.png') 
              
        else:
            x_start_plt, x_end_plt = x_start_new, x_end_new
            y_min_plt, y_max_plt = y_min_new, y_max_new         
            fig_out_dir = os.path.join(out_dir, '_bound_topsys_input.png')
    
        #   start with the RIV input in ax1, self.riv_arr_in = [riv_lay, 0, i, self.riv_head_elev[i], self.riv_cond[i], self.riv_bot_elev[i]]
        riv_x, riv_cond, riv_depth = [], [], []
        for riv_cell in self.riv_arr_in[0]:
            #   the RIV cells will be visualized as points and the distance along the x axis is set to be as middle of each river cell
            riv_x.append(round(riv_cell[2] * 0.1 + x_start_plt + 0.05, 2))
            riv_cond.append(round(riv_cell[4], 1))
            riv_head = round(riv_cell[3], 1)
            riv_bot = round(riv_cell[5], 1)
            print (riv_head, riv_bot)
            riv_depth.append(riv_head - riv_bot)
    
        ax1.plot(riv_x, riv_cond, 'bo', markersize = 5, markeredgewidth = 0.0)
        #ax1.set_ylabel('Conductance', color = 'b')
        yticks = ax1.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax1.tick_params('y', colors = 'k')
        ax1.set_xlim([x_start_plt, x_end_plt])
    
        #   next plot the RIV depth   
        ax2.plot(riv_x, riv_depth, 'ro', markersize = 5, markeredgewidth = 0.0)
        yticks = ax2.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)    
        #ax2.set_ylabel('Riv. depth', color = 'r')
        ax2.tick_params('y', colors = 'k')
    
        #   next plot the RCH rates
        rch_x, rch_rate = [], []
        for i in range(self.rch_arr.shape[-1]):
            rch_x.append(round(i * 0.1 + 0.05 + x_start_plt, 2))
            rch_rate.append(self.rch_arr[0, i] * 1000)
    
        ax3.plot(rch_x, rch_rate, 'go-', linewidth = 1, markersize = 5, markeredgewidth = 0.0)
        yticks = ax3.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        #ax3.set_ylabel('Rch rate (mm/d)', color = 'k')
        ax3.tick_params('y', colors = 'k')
    
        #   next the GHB head elevation - initial, the concentration of the GHB cells and the GHB conductivity 
        ghb_x, ghb_head, ghb_conc, ghb_cond = [], [], [], []    
        for j in range(len(self.ghb_input_lst)):
            ghb_x.append(round(self.ghb_input_lst[j][2] * 0.1 + 0.05 + x_start_plt, 2))
            ghb_head.append(round(self.ghb_input_lst[j][3], 2))
            ghb_cond.append(round(self.ghb_input_lst[j][4], 2))
            ghb_conc.append(self.ssmdata[j][3])
            
        ax4.plot(ghb_x, ghb_head, 'ko-', linewidth = 1, markersize = 5, markeredgewidth = 0.0)
        yticks = ax4.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        #ax4.set_ylabel('GHB head', color = 'k')
        ax4.tick_params('y', colors = 'k')
    
        ax5.plot(ghb_x, ghb_cond, 'ko-', linewidth = 1, markersize = 5, markeredgewidth = 0.0)
        yticks = ax5.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        #ax5.set_ylabel('GHB cond', color = 'k')
        ax5.tick_params('y', colors = 'k')
        
        ax6.plot(ghb_x, ghb_conc, 'ko-', linewidth = 1, markersize = 5, markeredgewidth = 0.0)
        yticks = ax6.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        #ax6.set_ylabel('GHB conc', color = 'k')
        ax6.tick_params('y', colors = 'k')
    
        #   create a plot if desired
        plot_hk_arr = self.hk_arr
        plot_hk_arr[np.abs(plot_hk_arr) == 0.] = np.nan
        #   create empty lists and fill them for each value
        #color_patches = [], [], []
        #   create a list with Hk values for the discretiazion
        unique_n = [0., 1E-06, 1E-03, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,\
                    18.0, 20.0, 22.0, 24.0, 26.0, 30.0]
        unique_n_label = ['0.', '1E-06', '1E-03', '0.1', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '4.0', '5.0', '6.0', '7.0', '8.0',\
                          '9.0', '10.0', '12.0', '14.0', '16.0', '18.0', '20.0', '22.0', '24.0', '26.0', '30.0']
    
        #   define the colormap
        cmap = plt.cm.viridis
        #cmaplist = [cmap(i) for i in range(cmap.N)]
        #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
        norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)
        
        im1 = ax7.imshow(plot_hk_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_plt, x_end_plt, min(self.bot_elev), math.ceil((self.top / 100.0) * 100.0)), vmin = 0.0, vmax = 30.0)
        y, x, z = self.dis.get_node_coordinates()
        #x = np.linspace(self.x_start, self.x_end, self.end_idx + 1)
        #x_lines_topelev = np.linspace(self.x_start + (self.del_col / (2 * 1000.)), self.x_end - (self.del_col / (2 * 1000.)), self.ncol)
        #lineA, = ax1.plot(x_lines_topelev, self.top_elev, c = 'black', linewidth = 3, label = 'GEBCO elevation')
        ax7.set_xlim([x_start_plt, x_end_plt])
        #ax1.set_ylim([min(self.bot_elev), max(self.top_elev)])
        ax7.axhline(y = 0., linewidth = 1, color = 'k')
        ax7.axvline(x = 0., linewidth = 1, color = 'k')    
    
        ax1.axhline(y = sea_lvl, linewidth = 3, color = 'red', zorder = 2)    
        ax2.axhline(y = sea_lvl, linewidth = 3, color = 'red', zorder = 2)        
    
        #   plot all the texts
        ax1b.text(0.1, 0.5, 'RIV input data \n Conductance (m2/d)')
        ax1b.axis('off')    
        ax2b.text(0.1, 0.5, 'RIV input data \n River depth (m)')
        ax2b.axis('off')
        ax3b.text(0.1, 0.5, 'RCH input data \n Recharge rate (mm/d)')
        ax3b.axis('off') 
        ax4b.text(0.1, 0.5, 'GHB input data \n Head elevation (m asl.)')
        ax4b.axis('off')
        ax5b.text(0.1, 0.5, 'GHB input data \n Concentration (g/l)')
        ax5b.axis('off')
        ax6b.text(0.1, 0.5, 'GHB input data \n Conductance (m2/d)')
        ax6b.axis('off')
        
        #   colorbar for the HK array
        cbar = plt.colorbar(im1, cax = ax7b, cmap = cmap, norm = norm, spacing = 'uniform', ticks = unique_n, boundaries = unique_n)
        cbar.ax.set_yticklabels(unique_n_label)
        cbar.ax.tick_params(labelsize = 12)        
        ax_cb = cbar.ax
        ax_cb.text(-.75, 0.85, 'Horizontal hydraulic conductivity (m/d)', fontsize = 14, rotation = 90)
    
        #   set the font 
        rcParams['font.family'] = 'Garamond'
        rcParams['axes.facecolor'] = 'white'
        rcParams['savefig.facecolor'] = 'white'    
    
        #   set the gridlines and constant lines in the plot
        x_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_plt / 5.) * 5., math.ceil(x_end_plt), 5.0), nbins = None)    
        x_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_plt / 1.) * 1., math.ceil(x_end_plt), 1.0), nbins = None)    
        ax7.xaxis.set_major_locator(x_major_locator)
        ax7.xaxis.set_minor_locator(x_minor_locator)
     
        y_major_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(y_min_plt / 500.) * 500., math.ceil(y_max_plt), 500.0), nbins = None)    
        y_minor_locator = matplotlib.ticker.FixedLocator(np.arange(math.floor(y_min_plt / 100.) * 100. , math.ceil(y_max_plt), 100.0), nbins = None)    
        ax7.yaxis.set_major_locator(y_major_locator)
        ax7.yaxis.set_minor_locator(y_minor_locator)
         
        ax1.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
        ax1.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)
        ax2.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
        ax2.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
        ax3.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
        ax3.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
        ax4.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
        ax4.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
        ax5.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
        ax5.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
        ax6.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
        ax6.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)  
        ax7.grid(which = 'major', color = 'grey', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = 1)    
        ax7.grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = .5, alpha = 0.5, zorder = 1)     
    
        ax1.set_title(title, fontsize = 24, y = 1.15)
    
        if os.path.exists(fig_out_dir):
            print('Figure   ' + fig_out_dir + '   already exists.')
        else: 
            #   save the figure
            plt.savefig(fig_out_dir, dpi = 300)
            plt.close(fig)
            #del fig
            

    #   Plotting of 2D concentration profiles at given time
    def plot_profile_yrs(self, dest_dir, time_step_yrs, in_dict, in_cbc_dict, sea_lvl, type_profile = 'conc', zoom = [-10., 10.]):#(self, time_yrs_lst):
    
        #   first check if the output folder exists, if not create it
        out_dir = os.path.join(dest_dir, '_' + type_profile + '_pngs')
        if os.path.exists(out_dir):
            pass
        else:
            os.makedirs(out_dir)
            
        #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
        fig = plt.figure(figsize = (20, 20))
        ax1 = plt.subplot2grid((2, 2), (0, 1))  #   overall concentration profiles
        ax2 = plt.subplot2grid((2, 2), (1, 1))  #   zoomed in concentration profile
        ax3 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)  #   color bar area 
        ax1.set_position([0.135, 0.55, 0.85, 0.4]) # [left, bottom, width, height]
        ax2.set_position([0.135, 0.1, 0.85, 0.4])
        ax3.set_position([0.04, 0.2, 0.025, 0.6])        
        
        #   define the extent of the figures
        x_start_new = self.x_start #-(self.cst_idx - self.idx_start) / 2.
        x_end_new = self.x_end
        y_max_new = math.ceil((self.top / 100.0) * 100.0)
        y_min_new = min(self.bot_elev)    
        
        if zoom:
            if x_start_new > zoom[0]: 
                x_start_new_zoom = x_start_new
                y_max_new_zoom = math.ceil((self.top / 100.0) * 100.0)
            else:
                x_start_new_zoom = zoom[0]
                y_max_new_zoom = math.ceil((self.top_elev[int(abs(x_start_new_zoom - zoom[0]) * 10)] / 100.0) * 100.0)
        
            if x_end_new < zoom[1]:                
                x_end_new_zoom = x_end_new
                y_min_new_zoom = min(self.bot_elev)
            else:
                x_end_new_zoom = zoom[1]
                y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / self.del_col)) + (zoom[1] * (1000. / self.del_col)))] * 100)
                #y_min_new_zoom = math.floor((self.bot_elev[int((self.cst_idx - self.idx_start) * 5 + (zoom[1] * (1000. / self.del_col)))] / 100.) * 100.)       
    
        #   add the location of the aquitard layers    
        mapimg = (self.hk_arr == 1e-05)
           
        #   read in the concentration profile from the dictionary, find the closest time step to the input tume_step_yrs
        real_time_yrs = min(in_dict.keys(), key = lambda k: abs(k - time_step_yrs))
        plot_arr = in_dict[min(in_dict.keys(), key = lambda k: abs(k - time_step_yrs))][type_profile]
        
        if type_profile == 'conc':
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
            fig_out_dir = os.path.join(out_dir, '_conc_' + str(real_time_yrs) + '_yrs.png')    
        
        elif type_profile == 'head': 
            vmin = np.nanmin(plot_arr)
            vmax = np.nanmax(plot_arr)
            title = 'Heads profile and flow pattern at ' + str(real_time_yrs) + ' years since simulation start'
    
            key_max = max(in_dict.keys(), key=(lambda k: np.nanmax(in_dict[k]['head'])))
            key_min = min(in_dict.keys(), key=(lambda k: np.nanmin(in_dict[k]['head'])))
            
            vmax = math.ceil(np.nanmax(in_dict[key_max]['head']))
            vmin = math.floor(np.nanmin(in_dict[key_min]['head']))
            
            cmap = plt.cm.plasma
            cmaplist = [cmap(i) for i in range(cmap.N)]
            cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
            #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)        
            
            bounds = np.linspace(vmin, vmax, 20)
            bounds = [round(x, 1) for x in bounds]
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
            cbar_title = 'Head elevation (m asl.)'
            
            #   get the flow arrays from the dictionary
            qx1 = in_cbc_dict[min(in_dict.keys(), key = lambda k: abs(k - time_step_yrs))]['q_right']
            qz1 = in_cbc_dict[min(in_dict.keys(), key = lambda k: abs(k - time_step_yrs))]['q_bottom']
            
            qx1[np.abs(qx1) == 0.] = np.nan        
            qz1[np.abs(qz1) == 0.] = np.nan        
            
            qx1_avg = np.empty(qx1.shape, dtype = qx1.dtype)
            qz1_avg = np.empty(qz1.shape, dtype = qz1.dtype)
            qx1_avg[:, :, :] = 0.5 * (qx1[:, :, 0 : self.ncol] + qx1[:, :, : self.ncol])
            qx1_avg[:, :, 0] = 0.5 * qx1[:, :, 0]
            qz1_avg[:, :, :] = 0.5 * (qz1[0 : self.nlay, :, :] + qz1[: self.nlay, :, :])
            qz1_avg[0, :, :] = 0.5 * qz1[0, :, :]
            
            y, x, z = self.dis.get_node_coordinates()
            x = np.linspace(x_start_new, x_end_new, plot_arr.shape[-1])
            X, Z = np.meshgrid(x, z[:, 0, 0])
            iskip = 3
    
            #   create output name for the individual figure
            fig_out_dir = os.path.join(out_dir, '_head_' + str(real_time_yrs) + '_yrs.png')           
            
        im1 = ax1.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax1.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
        
        ax2.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax2.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
        
        if type_profile == 'head':
            
            ax1.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip], qx1_avg[::iskip, 0, ::iskip], -qz1_avg[::iskip, 0, ::iskip],
               color = 'white', scale = 9, headwidth = 2, headlength = 2, headaxislength = 2, width = 0.0015)
            ax2.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip], qx1_avg[::iskip, 0, ::iskip], -qz1_avg[::iskip, 0, ::iskip],
               color = 'white', scale = 9, headwidth = 3, headlength = 2, headaxislength = 2, width = 0.0015)    
        
        ax1.set_xlim([x_start_new, x_end_new])
        ax1.set_ylim([y_min_new, y_max_new])
           
        ax2.set_xlim([x_start_new_zoom, x_end_new_zoom])
        ax2.set_ylim([y_min_new_zoom, y_max_new_zoom])
        
        #   set the gridlines and constant lines in the plot
        #   add constant lines with elevation = 0m asl. and coastline 
        ax1.axhline(y = 0., linewidth = 2, color = 'k', zorder = 2)
        ax1.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)
        ax2.axhline(y = 0., linewidth = 2, color = 'k', zorder = 2)
        ax2.axvline(x = 0., linewidth = 2, color = 'k', zorder = 2)    
        
        ax1.axhline(y = sea_lvl, linewidth = 3, color = 'red', zorder = 2)    
        ax2.axhline(y = sea_lvl, linewidth = 3, color = 'red', zorder = 2)    
    
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
    
        x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 5.) * 5., math.ceil(x_end_new_zoom), 5.0), nbins = None)    
        x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 1.) * 1., math.ceil(x_end_new_zoom), 1.0), nbins = None)    
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
            plt.close(fig)
            #del fig
            
            


    def plot_Kh_arr(self):
        #   calculate the offset of the coastline for plotting - this happens because the first offset is calculated from the
        #   raw topo data (from database with 0.5km spacing). When creating the top_elev interpolation the coast can be shifted by a bit
        cst_offset_val  = next((x for x in self.top_elev if x < 0.0), None)
        cst_offset_idx = self.top_elev.index(cst_offset_val) - 1
        self.cst_offset_plot = round(self.x_start + cst_offset_idx / (1000. /self.del_col), 2) #  (1000. / model.del_col) - to get the distance in km
        #   create a plot if desired
        plot_hk_arr = self.hk_arr
        plot_hk_arr[np.abs(plot_hk_arr) == 0.] = np.nan
        #   find all the unique values in the KH array, remove the ones that are marked as no-value
        unique = np.unique(self.hk_arr)
        unique_n = [value for value in unique if not math.isnan(value)]
        #   create empty lists and fill them for each value
        color_lst, color_name_lst, color_patches = [], [], []
        for i in range(len(unique_n)):
            for key in k_soil_dict:
                if k_soil_dict[key]['k_val'] == unique_n[i]:
                    color_lst.append(k_soil_dict[key]['plot_cl'])
                    color_name_lst.append(k_soil_dict[key]['soil_name'])
                    color_patches.append(mpatches.Patch(color = k_soil_dict[key]['plot_cl'],\
                    label = k_soil_dict[key]['soil_name'] + ' (Kh = ' + str(k_soil_dict[key]['k_val']) + ' m/d) '))
        #   add aquifer color
        color_lst.append('green')
        color_name_lst.append('Sandy aquifer')
        color_patches.append(mpatches.Patch(color = 'green', label = 'Sandy aquifer ' + '(Kh = ' + str(unique_n[-1]) + ' m/d)'))
        unique_n.append(20.)
        #   define the colormap
        cmap, norm = from_levels_and_colors(unique_n, color_lst)

        # Concentration and Flow Plot
        fig1 = plt.figure(figsize=(20, 5))
        ax1 = fig1.add_subplot(1, 1, 1)
        im1 = ax1.imshow(plot_hk_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (self.x_start, self.x_end - self.cst_offset_plot, self.zbot, math.ceil((self.top / 100.0) * 100.0)), vmin = 0., vmax = np.amax(plot_hk_arr))
                       #extent = (self.x_start - self.cst_offset_plot, self.x_end - self.cst_offset_plot, self.zbot, math.ceil((self.top / 100.0) * 100.0)), vmin = 0., vmax = np.amax(plot_hk_arr))
        y, x, z = self.dis.get_node_coordinates()
        #x = np.linspace(self.x_start, self.x_end, self.end_idx + 1)
        x_lines_topelev = np.linspace(self.x_start + (self.del_col / (2 * 1000.)), self.x_end - (self.del_col / (2 * 1000.)), self.ncol)
        lineA, = ax1.plot(x_lines_topelev, self.top_elev, c = 'black', linewidth = 3, label = 'GEBCO elevation')
        ax1.set_xlim([self.x_start, self.x_end])
        ax1.set_ylim([min(self.top_elev), max(self.top_elev)])
        ax1.axhline(y = 0., linewidth = 1, color = 'k')
        ax1.axvline(x = 0., linewidth = 1, color = 'k')
        ax1.set_title('Horizontal hydraulic conductivity distribution in the model domain')
        figname = '_kh.png'
        ax1.legend(handles = [lineA] + color_patches, loc='best')#, bbox_to_anchor=(0.85, 0.9))
        plt.savefig(os.path.join(self.foldername, figname))
        plt.close(fig1)
        #del fig1



    def plot_Kh_arr_GLHMYPS_rand(self, top_soil = False):
        #   calculate the offset of the coastline for plotting - this happens because the first offset is calculated from the
        #   raw topo data (from database with 0.5km spacing). When creating the top_elev interpolation the coast can be shifted by a bit
        cst_offset_val  = next((x for x in self.top_elev if x < 0.0), None)
        cst_offset_idx = self.top_elev.index(cst_offset_val) - 1
        self.cst_offset_plot = round(self.x_start + cst_offset_idx / (1000. /self.del_col), 2) #  (1000. / self.del_col) - to get the distance in km
        #   create a plot if desired
        plot_hk_arr = self.hk_arr
        plot_hk_arr[np.abs(plot_hk_arr) == 0.] = np.nan
        #   find all the unique values in the KH array, remove the ones that are marked as no-value
        unique = np.unique(self.hk_arr)
        unique_nan = [value for value in unique if not math.isnan(value)]
        #   define the colormap
        cmap = cmx.viridis
        cmaplist = [cmap(i) for i in range(cmap.N)]
        # define the bins and normalize
        unique_n = [0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 2.5, 5.0] + list(np.arange(10.0, round(max(unique_nan) / 10) * 10, 10.))
        norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)            

        # Concentration and Flow Plot
        fig = plt.figure(figsize = (20, 20))
        ax1 = plt.subplot2grid((2, 1), (0, 0))  #   overall concentration profiles
        ax2 = plt.subplot2grid((2, 1), (1, 0))  #   color bar area 

        ax1.set_position([0.05, 0.15, 0.9, 0.8])
        ax2.set_position([0.3, 0.025, 0.4, 0.05])        
            
        im1 = ax1.imshow(plot_hk_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (self.x_start, self.x_end, self.zbot, math.ceil((self.top / 100.0) * 100.0)), vmin = 0., vmax = np.amax(plot_hk_arr))
                       #extent = (self.x_start - self.cst_offset_plot, self.x_end - self.cst_offset_plot, self.zbot, math.ceil((self.top / 100.0) * 100.0)), vmin = 0., vmax = np.amax(plot_hk_arr))
        y, x, z = self.dis.get_node_coordinates()
        #x = np.linspace(self.x_start, self.x_end, self.end_idx + 1)
        x_lines_topelev = np.linspace(self.x_start + (self.del_col / (2 * 1000.)), self.x_end - (self.del_col / (2 * 1000.)), self.ncol)
        lineA, = ax1.plot(x_lines_topelev, self.top_elev, c = 'black', linewidth = 3, label = 'GEBCO elevation')
        ax1.set_xlim([self.x_start, self.x_end])
        #ax1.set_ylim([min(self.top_elev), max(self.top_elev)])
        ax1.set_ylim([min(self.bot_elev), max(self.top_elev)])
        ax1.axhline(y = 0., linewidth = 1, color = 'k')
        ax1.axvline(x = 0., linewidth = 1, color = 'k')
        ax1.set_title('Horizontal hydraulic conductivity distribution in the self domain', fontsize = 26)
        ax1.set_xlabel('distance from coast (km)', fontsize = 14)
        ax1.set_ylabel('elevation (m asl.)', fontsize = 14)

        figname = '_kh.png'

        #   plot the colorbar
        cbar = plt.colorbar(im1, cax = ax2, cmap = cmap, norm = norm, spacing = 'uniform', ticks = unique_n, boundaries = unique_n, orientation='horizontal')
        ax2.set_xticklabels([str(e) for e in unique_n])   
        ax2.set_title("Horizontal hydraulic conductivity (m/d)", fontsize = 14)#, y  = 1.025)
        ax2.tick_params(labelsize = 12)    

        rcParams['font.family'] = 'Garamond'
        rcParams['axes.facecolor'] = 'white'
        rcParams['savefig.facecolor'] = 'white'  

        plt.savefig(os.path.join(self.foldername, figname))
        plt.close(fig)
        #del fig


    #   function that compares model output for two different time steps and returns the min, max, avg difference
    #   in both heads and concentrations between these two outputs
    #def compare_ts_profiles(self, ts_1, ts_2):
    #    conc_ts_1 =
    #    conc_ts_2 =
    #    head_ts_1 =
    #   head_ts_2 =

    def out_video_graphs_conc(self, dict_in):

        #   create matplotlib objects necessary to export the video
        fig = plt.figure(figsize = (20, 5))
        ax = fig.add_subplot(111)

        #   define the template for showing time
        time_template = 'Time (years): = %.1f'

        # ims is a list of lists, each row is a list of artists to draw in the
        # current frame; here we are just animating one artist, the image, in each frame
        ims = []
        fresh_cnt, brackish_cnt, salt_cnt, tot_mass = [], [], [], []
        fresh_comparison, brackish_comparison, salt_comparison = [], [], []

        conc_0 = dict_in[0]['conc']
        conc_0[np.abs(conc_0) > 35.] = np.nan
        fresh_0 = ((conc_0 > -1.) & (conc_0 <= 0.05)).sum()
        brackish_0 = ((conc_0 > 0.05) & (conc_0 <= 15.0)).sum()
        salt_0 = ((conc_0 > 15.0) & (conc_0 <= 35.0)).sum()
        #   calculate the total mass in the system
        total_mass_0 = (conc_0 * (self.del_lay * self.del_col)).sum()

        fresh_comparison.append(0)
        brackish_comparison.append(0)
        salt_comparison.append(0)
        fresh_cnt.append(fresh_0)
        brackish_cnt.append(brackish_0)
        salt_cnt.append(salt_0)
        tot_mass.append(total_mass_0)

        for key in range(1, len(dict_in.keys())):

            time_key = dict_in.keys()[key]
            #print(key), time_key

            #time_key = key
            conc_key = dict_in[time_key]['conc']
            conc_key[np.abs(conc_key) > 35.5] = np.nan  #   35.5 because sometimes the concentrations are > 35.0

            fresh = ((conc_key > -1.) & (conc_key <= 0.05)).sum()
            brackish = ((conc_key > 0.05) & (conc_key <= 15.0)).sum()
            salt = ((conc_key > 15.0) & (conc_key <= 35.0)).sum()
            fresh_comparison.append(0)
            brackish_comparison.append(0)
            salt_comparison.append(0)
            fresh_cnt.append(fresh_0 - fresh)
            brackish_cnt.append(brackish_0 - brackish)
            salt_cnt.append(salt_0 - salt)

            #   calculate the total mass in the system
            total_mass = (conc_key * (self.del_lay * self.del_col)).sum()
            tot_mass.append(total_mass)

            #print conc_key

            #   draw the concentration profile for given time step
            im = ax.imshow(conc_key[:, 0, :], aspect='auto', interpolation='none',
                               extent = (self.x_start, self.x_end, self.botm[-1], self.top), vmin = 0, vmax = 35, animated = True)
            ax.axhline(y = 0, xmin = self.x_start, xmax = self.x_end, lw = 0.5, color = 'grey')
            ax.axvline(x = 0, ymin = self.botm[-1], ymax = self.top, lw = 0.5, color = 'grey')
            #   assign the annotation for given time step
            t = ax.annotate(time_template % (time_key), (self.x_start + 3000, self.botm[-1] + 100))
            #   append to the final list of artists
            ims.append([im, t])

        #   create the animation itself, interval is miliseconds between each frame,
        #   and repeat_delay is delay between the loop starts again
        ani = animation.ArtistAnimation(fig, ims, interval = 150, blit = False, repeat_delay = 1000)

        #   save the final video
        fig.tight_layout()
        ani.save(os.path.join(self.foldername, 'conc_through_time.mp4'))
        del fig

        #   create the output plots
        f = plt.figure(figsize=(20, 10))
        #   first the plot with fresh/brackish/salt cells count
        ax1 = f.add_subplot(3, 1, 1)
        f.add_subplot()
        ax1.set_title('Change in number of fresh, brackish and sea water concentration cells in the model domain')
        ax1.plot(dict_in.keys(), fresh_cnt, color = 'b')
        ax1.plot(dict_in.keys(), brackish_cnt, color = 'g')
        ax1.plot(dict_in.keys(), salt_cnt, color = 'r')
        #   also add the graph with differences in cell count
        ax2 = f.add_subplot(3, 1, 2)
        f.add_subplot()
        ax2.set_title('Changes in cell count compared to 1st time step (1st ts - nth ts)')
        ax2.plot(dict_in.keys(), fresh_comparison, color = 'b')
        ax2.plot(dict_in.keys(), brackish_comparison, color = 'g')
        ax2.plot(dict_in.keys(), salt_comparison, color = 'r')

        #   also add the total mass graph in the system
        ax3 = f.add_subplot(3, 1, 3)
        f.add_subplot()
        ax3.set_title('Total mass in the model domain')
        ax3.plot(dict_in.keys(), tot_mass, color = 'black')
        #   save the output
        figname = 'model_output_graphs.png'
        plt.savefig(os.path.join(self.foldername, figname))
        plt.close(f)
        #del f


    #   Plotting of 2D concentration profiles at given time
    def plot_conc_profile_yrs(self, dest_dir, x_topo, plot_arr, sea_lvl, real_time_yrs, cst_shift, zoom = [-10., 10.], conc_diff = False):#(self, time_yrs_lst):
    
        #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
        fig = plt.figure(figsize = (20, 12))
        ax1 = plt.subplot2grid((2, 2), (0, 1))  #   overall concentration profiles
        ax2 = plt.subplot2grid((2, 2), (1, 1))  #   zoomed in concentration profile
        ax3 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)  #   color bar area 
        ax1.set_position([0.135, 0.55, 0.85, 0.4]) # [left, bottom, width, height]
        ax2.set_position([0.135, 0.1, 0.85, 0.4])
        ax3.set_position([0.04, 0.2, 0.025, 0.6])        
        
        #   define the extent of the figures
        x_start_new = self.x_start #-(self.cst_idx - self.idx_start) / 2.
        x_end_new = self.x_end
        #y_max_new = math.ceil((self.top / 10.0) * 10.0)
        y_max_new = self.top
        y_min_new = min(self.bot_elev)    
        
        if zoom:
            if x_start_new > zoom[0]: 
                x_start_new_zoom = x_start_new
                y_max_new_zoom = math.ceil((self.top / 10.0) * 10.0)
            else:
                x_start_new_zoom = zoom[0]
                y_max_new_zoom = math.ceil(self.top_elev[int(abs(x_start_new_zoom - zoom[0]) * 10)] / 10.0) * 10.0
        
            if x_end_new < zoom[1]:                
                x_end_new_zoom = x_end_new
                y_min_new_zoom = min(self.bot_elev)
            else:
                x_end_new_zoom = zoom[1]
                #y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / max(100., self.del_col_res))) + (zoom[1] * (1000. / max(100., self.del_col_res))))] / 100.) * 100.
                #y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / self.del_col_res)) + (zoom[1] * (1000. / self.del_col_res)))] / 100.) * 100.
                try:
                    y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / self.del_col_res)) + (zoom[1] * (1000. / self.del_col_res)))] / 10.) * 10.
                except IndexError:
                    y_min_new_zoom = math.floor(self.bot_elev[-1] / 10.) * 10.
                
                #y_min_new_zoom = math.floor((self.bot_elev[int((self.cst_idx - self.idx_start) * 5 + (zoom[1] * (1000. / self.del_col)))] / 100.) * 100.)       
    
        #   add the location of the aquitard layers    
        mapimg = (self.hk_arr == 1e-04)
           
        #   read in the concentration profile from the dictionary, find the closest time step to the input tume_step_yrs
        if not conc_diff:
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
            fig_out_dir = os.path.join(dest_dir, '_conc_' + str(real_time_yrs) + '_yrs.png')   

        else:
            vmin = np.nanmin(plot_arr)
            vmax = np.nanmax(plot_arr)
            title = 'Difference in concentration profile between ' + str(int(real_time_yrs - 100.)) + ' and '+ str(real_time_yrs) + ' years since simulation start'
            cmap = plt.cm.jet
            cmaplist = [cmap(i) for i in range(cmap.N)]
            cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
            # define the bins and normalize
            bounds = np.linspace(-35., -5., 4).tolist() + [-2.5, -1.0, -0.5, -0.1, -0.05, 0.0, 0.05, 0.1, 0.5, 1.0, 2.5] + np.linspace(5.0, 35., 4).tolist()
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
            cbar_title = 'Difference in concentration of \n salt (g/l)'
            plot_arr[np.abs(plot_arr) >= 999.] = np.nan            
            #   create output name for the individual figure
            fig_out_dir = os.path.join(dest_dir, '_conc_diff_' + str(int(real_time_yrs - 100.)) + '_to_' + str(real_time_yrs) + '_yrs.png')    
        
        im1 = ax1.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax1.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
        
        ax2.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax2.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
        
        ax1.set_xlim([x_start_new, x_end_new])
        ax1.set_ylim([y_min_new, y_max_new])
           
        ax2.set_xlim([x_start_new_zoom, x_end_new_zoom])
        ax2.set_ylim([y_min_new_zoom, y_max_new_zoom])
        
        #   set the gridlines and constant lines in the plot
        #   add constant lines with elevation = 0m asl. and coastline 
        ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
 
        #ax1.axvline(x = 0. + cst_shift, linewidth = 2, color = 'k', zorder = 2)
        #ax2.axvline(x = 0. + cst_shift, linewidth = 2, color = 'k', zorder = 2)    
       
        ax1.axhline(y = sea_lvl, linewidth = 2, color = 'k', zorder = 2)    
        ax2.axhline(y = sea_lvl, linewidth = 2, color = 'k', zorder = 2)    
    
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
    
        x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 5.) * 5., math.ceil(x_end_new_zoom), 5.0), nbins = None)    
        x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 1.) * 1., math.ceil(x_end_new_zoom), 1.0), nbins = None)    
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

        #ax1.plot(x_topo, self.top_elev, c = 'lime', linewidth = 1.0)
        #ax2.plot(x_topo, self.top_elev, c = 'lime', linewidth = 1.0)
        
        #   set the font 
        rcParams['font.family'] = 'Garamond'
        rcParams['axes.facecolor'] = 'white'
        rcParams['savefig.facecolor'] = 'white'    
        
        if os.path.exists(fig_out_dir):
            print('Figure   ' + fig_out_dir + '   already exists.')
        else: 
            #   save the figure
            plt.savefig(fig_out_dir, dpi = 300)
            plt.close(fig)
            #del fig


    #   Plotting of 2D concentration profiles at given time
    def plot_conc_flow_vectors_profile_yrs(self, dest_dir, plot_arr, sea_lvl, real_time_yrs, cst_shift, qx1, qz1, zoom = [-10., 10.], conc_diff = False):#(self, time_yrs_lst):
    
        #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
        fig = plt.figure(figsize = (20, 12))
        ax1 = plt.subplot2grid((2, 2), (0, 1))  #   overall concentration profiles
        ax2 = plt.subplot2grid((2, 2), (1, 1))  #   zoomed in concentration profile
        ax3 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)  #   color bar area 
        ax1.set_position([0.135, 0.55, 0.85, 0.4]) # [left, bottom, width, height]
        ax2.set_position([0.135, 0.1, 0.85, 0.4])
        ax3.set_position([0.04, 0.2, 0.025, 0.6])        
        
        #   define the extent of the figures
        x_start_new = self.x_start #-(self.cst_idx - self.idx_start) / 2.
        x_end_new = self.x_end
        #y_max_new = math.ceil((self.top / 10.0) * 10.0)
        y_max_new = self.top
        y_min_new = min(self.bot_elev)    
        
        if zoom:
            if x_start_new > zoom[0]: 
                x_start_new_zoom = x_start_new
                y_max_new_zoom = math.ceil((self.top / 10.0) * 10.0)
            else:
                x_start_new_zoom = zoom[0]
                y_max_new_zoom = math.ceil(self.top_elev[int(abs(x_start_new_zoom - zoom[0]) * 10)] / 10.0) * 10.0
        
            if x_end_new < zoom[1]:                
                x_end_new_zoom = x_end_new
                y_min_new_zoom = min(self.bot_elev)
            else:
                x_end_new_zoom = zoom[1]
                #y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / max(100., self.del_col_res))) + (zoom[1] * (1000. / max(100., self.del_col_res))))] / 100.) * 100.
                #y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / self.del_col_res)) + (zoom[1] * (1000. / self.del_col_res)))] / 100.) * 100.
                #y_min_new_zoom = math.floor((self.bot_elev[int((self.cst_idx - self.idx_start) * 5 + (zoom[1] * (1000. / self.del_col)))] / 100.) * 100.)       
                try:
                    y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / self.del_col_res)) + (zoom[1] * (1000. / self.del_col_res)))] / 10.) * 10.
                except IndexError:
                    y_min_new_zoom = math.floor(self.bot_elev[-1] / 10.) * 10.
        
        #   add the location of the aquitard layers    
        mapimg = (self.hk_arr == 1e-04)
           
        #   read in the concentration profile from the dictionary, find the closest time step to the input tume_step_yrs
        if not conc_diff:
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
            fig_out_dir = os.path.join(dest_dir, '_conc_' + str(real_time_yrs) + '_yrs.png')   
        
        else:
            vmin = np.nanmin(plot_arr)
            vmax = np.nanmax(plot_arr)
            title = 'Difference in concentration profile between ' + str(int(real_time_yrs - 100.)) + ' and '+ str(real_time_yrs) + ' years since simulation start'
            cmap = plt.cm.jet
            cmaplist = [cmap(i) for i in range(cmap.N)]
            cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
            # define the bins and normalize
            bounds = np.linspace(-35., -5., 4).tolist() + [-2.5, -1.0, -0.5, -0.1, -0.05, 0.0, 0.05, 0.1, 0.5, 1.0, 2.5] + np.linspace(5.0, 35., 4).tolist()
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
            cbar_title = 'Difference in concentration of \n salt (g/l)'
            plot_arr[np.abs(plot_arr) >= 999.] = np.nan            
            #   create output name for the individual figure
            fig_out_dir = os.path.join(dest_dir, '_conc_diff_' + str(int(real_time_yrs - 100.)) + '_to_' + str(real_time_yrs) + '_yrs.png')    
        
        if qx1 is not None and qz1 is not None:        
        
            qx1[np.abs(qx1) == 0.] = np.nan        
            qz1[np.abs(qz1) == 0.] = np.nan        
            
            qx1_avg = np.empty(qx1.shape, dtype = qx1.dtype)
            qz1_avg = np.empty(qz1.shape, dtype = qz1.dtype)
            qx1_avg[:, :, :] = 0.5 * (qx1[:, :, 0 : self.ncol] + qx1[:, :, : self.ncol])
            qx1_avg[:, :, 0] = 0.5 * qx1[:, :, 0]
            qz1_avg[:, :, :] = 0.5 * (qz1[0 : self.nlay, :, :] + qz1[: self.nlay, :, :])
            qz1_avg[0, :, :] = 0.5 * qz1[0, :, :]
            
            y, x, z = self.dis.get_node_coordinates()
            z = np.linspace(y_max_new, y_min_new, plot_arr.shape[0])
            x = np.linspace(x_start_new, x_end_new, plot_arr.shape[-1])
            X, Z = np.meshgrid(x, z)
            #X, Z = np.meshgrid(x, z[:, 0, 0])
            
            #   define the distances for each flowline in x and z directions
            x_q_dst = 1000
            y_q_dst = 100
            iskip_x = int(x_q_dst / self.del_col_res)
            iskip_z = int(y_q_dst / self.del_lay_res)
        
        im1 = ax1.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax1.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
        
        ax2.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax2.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
        
        if qx1 is not None and qz1 is not None:              
            ax1.quiver(X[::iskip_z, ::iskip_x], Z[::iskip_z, ::iskip_x], qx1_avg[::iskip_z, 0, ::iskip_x], -qz1_avg[::iskip_z, 0, ::iskip_x],
               color = 'white', scale = None, headwidth = 3, headlength = 3, headaxislength = 2, width = 0.0015)
            ax2.quiver(X[::iskip_z, ::iskip_x], Z[::iskip_z, ::iskip_x], qx1_avg[::iskip_z, 0, ::iskip_x], -qz1_avg[::iskip_z, 0, ::iskip_x],
               color = 'white', scale = None, headwidth = 3, headlength = 3, headaxislength = 2, width = 0.0015)    
        
        ax1.set_xlim([x_start_new, x_end_new])
        ax1.set_ylim([y_min_new, y_max_new])
           
        ax2.set_xlim([x_start_new_zoom, x_end_new_zoom])
        ax2.set_ylim([y_min_new_zoom, y_max_new_zoom])
        
        #   set the gridlines and constant lines in the plot
        #   add constant lines with elevation = 0m asl. and coastline 
        ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
         
        ax1.axvline(x = 0. + cst_shift, linewidth = 2, color = 'k', zorder = 2)
        ax2.axvline(x = 0. + cst_shift, linewidth = 2, color = 'k', zorder = 2)    
           
        ax1.axhline(y = sea_lvl, linewidth = 2, color = 'k', zorder = 2)    
        ax2.axhline(y = sea_lvl, linewidth = 2, color = 'k', zorder = 2)    
        
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
        
        x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 5.) * 5., math.ceil(x_end_new_zoom), 5.0), nbins = None)    
        x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 1.) * 1., math.ceil(x_end_new_zoom), 1.0), nbins = None)    
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
            plt.close(fig)
            #del fig



    #   Plotting of 2D concentration profiles at given time
    def plot_head_profile_yrs(self, dest_dir, x_topo, plot_arr, sea_lvl, real_time_yrs, cst_shift, qx1, qz1, zoom = [-10., 10.], head_diff = False):#(self, time_yrs_lst):
    
        #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
        fig = plt.figure(figsize = (20, 12))
        ax1 = plt.subplot2grid((2, 2), (0, 1))  #   overall concentration profiles
        ax2 = plt.subplot2grid((2, 2), (1, 1))  #   zoomed in concentration profile
        ax3 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)  #   color bar area 
        ax1.set_position([0.135, 0.55, 0.85, 0.4]) # [left, bottom, width, height]
        ax2.set_position([0.135, 0.1, 0.85, 0.4])
        ax3.set_position([0.04, 0.2, 0.025, 0.6])        
        
        #   define the extent of the figures
        x_start_new = self.x_start #-(self.cst_idx - self.idx_start) / 2.
        x_end_new = self.x_end
        #y_max_new = math.ceil((self.top / 10.0) * 10.0)
        y_max_new = self.top
        y_min_new = min(self.bot_elev)    
        
        if zoom:
            if x_start_new > zoom[0]: 
                x_start_new_zoom = x_start_new
                y_max_new_zoom = math.ceil((self.top / 10.0) * 10.0)
            else:
                x_start_new_zoom = zoom[0]
                y_max_new_zoom = math.ceil(self.top_elev[int(abs(x_start_new_zoom - zoom[0]) * 10)] / 10.0) * 10.0
        
            if x_end_new < zoom[1]:                
                x_end_new_zoom = x_end_new
                y_min_new_zoom = min(self.bot_elev)
            else:
                x_end_new_zoom = zoom[1]
                #y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / max(100., self.del_col_res))) + (zoom[1] * (1000. / max(100., self.del_col_res))))] / 100.) * 100.
                try:
                    y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / self.del_col_res)) + (zoom[1] * (1000. / self.del_col_res)))] / 10.) * 10.
                except IndexError:
                    y_min_new_zoom = math.floor(self.bot_elev[-1] / 10.) * 10.
                #y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / self.del_col)) + (zoom[1] * (1000. / self.del_col)))] / 100.) * 100.
                #y_min_new_zoom = math.floor((self.bot_elev[int((self.cst_idx - self.idx_start) * 5 + (zoom[1] * (1000. / self.del_col)))] / 100.) * 100.)       
    
        #   add the location of the aquitard layers    
        mapimg = (self.hk_arr == 1e-04)
           
        if not head_diff:
            #   read in the concentration profile from the dictionary, find the closest time step to the input tume_step_yrs
            title = 'Heads profile and flow pattern at ' + str(real_time_yrs) + ' years since simulation start'
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
            #   create output name for the individual figure
            fig_out_dir = os.path.join(dest_dir, '_head_' + str(real_time_yrs) + '_yrs.png')        

        else:
            title = 'Difference in heads between ' + str(int(real_time_yrs - 100.)) + ' and '+ str(real_time_yrs) + ' years since simulation start'
            plot_arr[np.abs(plot_arr) >= 999.] = np.nan
            vmax = 10.
            vmin = -10.
            cmap = plt.cm.plasma
            cmaplist = [cmap(i) for i in range(cmap.N)]
            cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) 
            #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)        
            bounds = np.linspace(vmin, -1.0, 5).tolist() + [-0.5, -0.1, -0.05, 0.0, 0.05, 0.1, 0.5] + np.linspace(1.0, vmax, 5).tolist()
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
            cbar_title = 'Difference in head elevation (m asl.)'        
            #   create output name for the individual figure
            fig_out_dir = os.path.join(dest_dir, '_head_diff_' + str(int(real_time_yrs - 100.)) + '_to_'+ str(real_time_yrs) + '_yrs.png')  
        
        if qx1 is not None and qz1 is not None:        
        
            qx1[np.abs(qx1) == 0.] = np.nan        
            qz1[np.abs(qz1) == 0.] = np.nan        
            
            qx1_avg = np.empty(qx1.shape, dtype = qx1.dtype)
            qz1_avg = np.empty(qz1.shape, dtype = qz1.dtype)
            qx1_avg[:, :, :] = 0.5 * (qx1[:, :, 0 : self.ncol] + qx1[:, :, : self.ncol])
            qx1_avg[:, :, 0] = 0.5 * qx1[:, :, 0]
            qz1_avg[:, :, :] = 0.5 * (qz1[0 : self.nlay, :, :] + qz1[: self.nlay, :, :])
            qz1_avg[0, :, :] = 0.5 * qz1[0, :, :]
            
            y, x, z = self.dis.get_node_coordinates()
            x = np.linspace(x_start_new, x_end_new, plot_arr.shape[-1])
            X, Z = np.meshgrid(x, z[:, 0, 0])
            iskip = 3

   
            
        im1 = ax1.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax1.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
        
        ax2.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax2.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
            
        if qx1 is not None and qz1 is not None:              
            ax1.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip], qx1_avg[::iskip, 0, ::iskip], -qz1_avg[::iskip, 0, ::iskip],
               color = 'white', scale = 9, headwidth = 2, headlength = 2, headaxislength = 2, width = 0.0015)
            ax2.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip], qx1_avg[::iskip, 0, ::iskip], -qz1_avg[::iskip, 0, ::iskip],
               color = 'white', scale = 9, headwidth = 3, headlength = 2, headaxislength = 2, width = 0.0015)    
        
        ax1.set_xlim([x_start_new, x_end_new])
        ax1.set_ylim([y_min_new, y_max_new])
           
        ax2.set_xlim([x_start_new_zoom, x_end_new_zoom])
        ax2.set_ylim([y_min_new_zoom, y_max_new_zoom])
        
        #   set the gridlines and constant lines in the plot
        #   add constant lines with elevation = 0m asl. and coastline 
        ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
 
        #ax1.axvline(x = 0. + cst_shift, linewidth = 2, color = 'k', zorder = 2)
        #ax2.axvline(x = 0. + cst_shift, linewidth = 2, color = 'k', zorder = 2)    
       
        ax1.axhline(y = sea_lvl, linewidth = 2, color = 'k', zorder = 2)    
        ax2.axhline(y = sea_lvl, linewidth = 2, color = 'k', zorder = 2)    
    
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
    
        x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 5.) * 5., math.ceil(x_end_new_zoom), 5.0), nbins = None)    
        x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 1.) * 1., math.ceil(x_end_new_zoom), 1.0), nbins = None)    
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

        #ax1.plot(x_topo, self.top_elev, c = 'lime', linewidth = 1.0)
        #ax2.plot(x_topo, self.top_elev, c = 'lime', linewidth = 1.0)
        
        #   set the font 
        rcParams['font.family'] = 'Garamond'
        rcParams['axes.facecolor'] = 'white'
        rcParams['savefig.facecolor'] = 'white'    
        
        if os.path.exists(fig_out_dir):
            print('Figure   ' + fig_out_dir + '   already exists.')
        else: 
            #   save the figure
            plt.savefig(fig_out_dir, dpi = 300)
            plt.close(fig)
            #del fig



    #   Plotting of 2D concentration profiles at given time
    def plot_HK_arr_res(self, dest_dir, x_topo, plot_arr, sea_lvl, real_time_yrs, cst_shift, zoom = [-10., 10.]):#(self, time_yrs_lst):
    
        #   define the layout of the figure - two plots below each other and a colorbar on the side with concentration legend
        fig = plt.figure(figsize = (20, 12))
        ax1 = plt.subplot2grid((2, 2), (0, 1))  #   overall concentration profiles
        ax2 = plt.subplot2grid((2, 2), (1, 1))  #   zoomed in concentration profile
        ax3 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)  #   color bar area 
        ax1.set_position([0.135, 0.55, 0.85, 0.4]) # [left, bottom, width, height]
        ax2.set_position([0.135, 0.1, 0.85, 0.4])
        ax3.set_position([0.04, 0.2, 0.025, 0.6])        
        
        #   define the extent of the figures
        x_start_new = self.x_start #-(self.cst_idx - self.idx_start) / 2.
        x_end_new = self.x_end
        #y_max_new = math.ceil((self.top / 10.0) * 10.0)
        y_max_new = self.top
        y_min_new = min(self.bot_elev)    
        
        if zoom:
            if x_start_new > zoom[0]: 
                x_start_new_zoom = x_start_new
                y_max_new_zoom = math.ceil((self.top / 100.0) * 100.0)
            else:
                x_start_new_zoom = zoom[0]
                y_max_new_zoom = math.ceil(self.top_elev[int(abs(x_start_new_zoom - zoom[0]) * 10)] / 100.0) * 100.0
        
            if x_end_new < zoom[1]:                
                x_end_new_zoom = x_end_new
                y_min_new_zoom = min(self.bot_elev)
            else:
                x_end_new_zoom = zoom[1]
                try:
                    y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / self.del_col_res)) + (zoom[1] * (1000. / self.del_col_res)))] / 100.) * 100.
                except IndexError:
                    y_min_new_zoom = -200.
                #y_min_new_zoom = math.floor(self.bot_elev[int(abs(self.x_start * (1000. / self.del_col)) + (zoom[1] * (1000. / self.del_col)))] / 100.) * 100.
                #y_min_new_zoom = math.floor((self.bot_elev[int((self.cst_idx - self.idx_start) * 5 + (zoom[1] * (1000. / self.del_col)))] / 100.) * 100.)       
    
        #   add the location of the aquitard layers    
        mapimg = (self.hk_arr == 1e-04)
           
        #   read in the concentration profile from the dictionary, find the closest time step to the input tume_step_yrs
        title = 'HK (horizontal hydraulic conductivity), del_col = ' + str(self.del_col_res) + ' del_lay = ' + str(self.del_lay_res)
        plot_arr[np.abs(plot_arr) >= 999.] = np.nan
        vmax = 30.0
        vmin = 0.0
        #   define the colormap
        cmap = cmx.viridis
        cmaplist = [cmap(i) for i in range(cmap.N)]
        # define the bins and normalize
        unique_n = [0.0, 0.001, 0.1, 1.0, 5.0, 10., 20., 30.]# + list(np.arange(10.0, round(max(unique_nan) / 10) * 10, 10.))
        norm = matplotlib.colors.BoundaryNorm(unique_n, cmap.N)                          
        #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)        
        bounds = unique_n
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        cbar_title = 'Head elevation (m asl.)'
        #   create output name for the individual figure
        fig_out_dir = os.path.join(dest_dir, '_HK_resolution_' + str(self.del_col_res) + '_' + str(self.del_lay_res) + '.png')        
            
        im1 = ax1.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax1.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)
        
        ax2.imshow(plot_arr[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = cmap, norm = norm,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), vmin = vmin, vmax = vmax)
        ax2.imshow(mapimg[:, 0, :], aspect = 'auto', interpolation = 'nearest', cmap = plt.cm.gray,
                       extent = (x_start_new, x_end_new, y_min_new, y_max_new), alpha = 0.25)

        ax1.plot(x_topo, self.top_elev, c = 'red', linewidth = 2.0)
        ax2.plot(x_topo, self.top_elev, c = 'red', linewidth = 2.0)
            
        ax1.set_xlim([x_start_new, x_end_new])
        ax1.set_ylim([y_min_new, y_max_new])
           
        ax2.set_xlim([x_start_new_zoom, x_end_new_zoom])
        ax2.set_ylim([y_min_new_zoom, y_max_new_zoom])

        #   set the gridlines and constant lines in the plot
        #   add constant lines with elevation = 0m asl. and coastline 
        ax1.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax1.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axhline(y = 0., linewidth = 1, color = 'k', zorder = 2)
        ax2.axvline(x = 0., linewidth = 1, color = 'k', zorder = 2)    
 
        ax1.axvline(x = 0. + cst_shift, linewidth = 2, color = 'k', zorder = 2)
        ax2.axvline(x = 0. + cst_shift, linewidth = 2, color = 'k', zorder = 2)    
       
        ax1.axhline(y = sea_lvl, linewidth = 2, color = 'k', zorder = 2)    
        ax2.axhline(y = sea_lvl, linewidth = 2, color = 'k', zorder = 2)    
    
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
    
        x_major_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 5.) * 5., math.ceil(x_end_new_zoom), 5.0), nbins = None)    
        x_minor_locator_zoom = matplotlib.ticker.FixedLocator(np.arange(math.floor(x_start_new_zoom / 1.) * 1., math.ceil(x_end_new_zoom), 1.0), nbins = None)    
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
            plt.close(fig)
            #del fig



"""
    #   create concentration videos from the concentration file
    def out_conc_video(self):
        #   get the index of the coastline and the extent of the domain in the x direction
        def get_coast_lay_col(delcol, numcol, lay_top_elev, ibound):
            #   find the first negative layer bottom elevation, get the index of the last positive value (- 1)
            first_neg_idx = len([i for i in self.top_elev if i > 0.0]) - 1
            #   select the layer from the IBOUND array
            sea_level_lay = self.ibound_arr[first_neg_idx, 0, :].tolist()
            #   get the index of the last active cell in the layer - that is the index of the coastline
            cst_idx = len(sea_level_lay) - 2 - sea_level_lay[::-1].index(1)
            #   get the x0 and x1 values
            x0 = self.x_start
            x1 = self.x_end
            #   give the output
            return cst_idx, x0, x1



        #   look for all the UCN files in the rootdir
        for subdir, dirs, files in os.walk(rootdir):
            for file in files:
                if file.endswith(".UCN"):
                    print os.path.join(subdir, file)
                    model_name = os.path.basename(subdir)

                    #   get the simulation number
                    #sim_num = subdir.split((os.sep))[-1]
                    #print sim_num,'...'

                    model_ws = subdir
                    ml = flopy.modflow.Modflow.load(model_name + ".nam", model_ws=model_ws, verbose=False, check=False, exe_name="mfnwt")
                    dis = flopy.modflow.ModflowDis.load(model_ws + '\\' + model_name + ".dis", ml)
                    delr = dis.delr[0]
                    nlay = dis.nlay
                    ncol = dis.ncol
                    top = dis.top[0]
                    botm_lst = dis.botm[:, 0, 0].tolist()
                    botm_lst = [top] + botm_lst
                    bottom = botm_lst[-1]
                    delv = (top - bottom)  / nlay

                    #   read in the .bas file
                    bas_file = flopy.modflow.ModflowBas.load(os.path.join(subdir, model_name + '.bas'), ml)
                    ibound_arr = bas_file.ibound[:,:,:]
                    print os.path.join(subdir, model_name + '.bas')

                    #   split the model name by '_'
                    model_name_split = model_name.split('_')
                    model_thickness = int(model_name_split[-1])
                    #   if there is 'offshore' in the model split, find the x0 and x1
                    if 'offshore' in model_name_split:
                        get_x = get_coast_lay_col(delr, ncol, botm_lst, ibound_arr)
                        cst_idx, x0, x1 = get_x[0], get_x[1], get_x[2]
                    else:
                        cst_idx = ncol
                        x0 = -delr * ncol
                        x1 = 0.

                    #   load the concnetration file into an array c
                    ucnobj = bf.UcnFile(os.path.join(subdir, 'MT3D001.UCN'),model=ml)
                    time_steps = ucnobj.get_times()

                    #   create matplotlib objects necessary to export the video
                    fig = plt.figure(figsize = (20, 5))
                    ax = fig.add_subplot(111)

                    #   define data generator, in this case the concentration profile in time
                    def data_gen(i):
                        c = ucnobj.get_data(totim = time_steps[i])
                        c[np.abs(c) > 36.] = np.nan
                        return c

                    #   assign the dimensions for plotting
                    # change to distance from coast!
                    y, x, z = dis.get_node_coordinates()
                    X, Z = np.meshgrid(x, z[:, 0, 0])

                    #   define the template for showing time
                    time_template = 'Time (years): = %.1f'

                    # ims is a list of lists, each row is a list of artists to draw in the
                    # current frame; here we are just animating one artist, the image, in each frame
                    ims = []
                    for j in range(0, len(time_steps)):
                        #   draw the concentration profile for given time step
                        im = ax.imshow(data_gen(j)[:,0,:], aspect='auto', interpolation='none',
                                           extent=(x0, x1, bottom, top), vmin=0, vmax=35, animated = True)
                        ax.axhline(y = 0, xmin = x0, xmax = x1, lw = 0.5, color = 'grey')
                        ax.axvline(x = 0, ymin = bottom, ymax = top, lw = 0.5, color = 'grey')
                        #   assign the annotation for given time step
                        t = ax.annotate(time_template % (round(time_steps[j] / 365.25, 0)), (x0 + 3000, bottom + 100))
                        #   append to the final list of artists
                        ims.append([im, t])

                    #   create the animation itself, interval is miliseconds between each frame,
                    #   and repeat_delay is delay between the loop starts again
                    ani = animation.ArtistAnimation(fig, ims, interval = 150, blit = False, repeat_delay = 1000)

                    #   save the final video
                    fig.tight_layout()
                    ani.save(os.path.join(rootdir, model_name + 'conc_through_time.mp4'))
                    del fig
"""


























