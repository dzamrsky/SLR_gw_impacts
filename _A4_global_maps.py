# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 11:03:51 2020

@author: daniel
"""


"""    ---------------------------------------------------------------    """
import pandas as pd
import numpy as np
import os

#   define paths
data_out_dir = r'g:\Water_Nexus\_A4_paper\_data_output'
#data_out_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models_OUT_files'

main_dir = r'g:\Water_Nexus\_A4_paper\_data_output\RCP_'
#main_dir = r'/projects/0/qt16165/_dzamrsky/_A4_SRM_models/_SLR_models_OUT_files/RCP_'

#   ------------------------------------------------------------------------------------
#   Create csv files showing mean and stdev (across all 3 DEMs) for each year 
#   ------------------------------------------------------------------------------------

out_cols = ['fid', 'coscat', 'srm', 'cs_srm_id', 'inl_ext_km',\
            'IFGW_2000_mu', 'IFGW_2000_std', 'IFGW_2050_mu', 'IFGW_2050_std',\
            'IFGW_2100_mu', 'IFGW_2100_std', 'IFGW_2200_mu', 'IFGW_2200_std',\
            'IFGW_2300_mu', 'IFGW_2300_std', 'IFGW_2400_mu', 'IFGW_2400_std',\
            'IFGW_2500_mu', 'IFGW_2500_std']

for rcp in ['26', '45', '85']:

    if not os.path.exists(main_dir + rcp + '_mu_std.csv'):
        df = pd.DataFrame(columns = out_cols)
        df.set_index('fid')
        df.to_csv(main_dir + rcp + '_mu_std_v2.csv' , encoding = 'utf-8', index = False)
    
    rcp_in_dir = main_dir + rcp + '_results.csv'
    rcp_in = pd.read_csv(rcp_in_dir)
    fid = 1
    
    for index, row in rcp_in.iterrows():
        
        cs_id = row['coscat']
        srm_id = row['srm']
        cs_srm_id = row['cs_srm_id']
        inl_ext = row['inl_ext_km']
        
        IFGW_2000_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2000'], row['IFGW_co_2000'], row['IFGW_ge_2000']] if x != -1])), 1)
        IFGW_2000_std = round(np.std(np.array([x for x in [row['IFGW_me_2000'], row['IFGW_co_2000'], row['IFGW_ge_2000']] if x != -1])), 1)    
        IFGW_2050_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2050'], row['IFGW_co_2050'], row['IFGW_ge_2050']] if x != -1])), 1)
        IFGW_2050_std = round(np.std(np.array([x for x in [row['IFGW_me_2050'], row['IFGW_co_2050'], row['IFGW_ge_2050']] if x != -1])), 1)    
        IFGW_2100_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2100'], row['IFGW_co_2100'], row['IFGW_ge_2100']] if x != -1])), 1)
        IFGW_2100_std = round(np.std(np.array([x for x in [row['IFGW_me_2100'], row['IFGW_co_2100'], row['IFGW_ge_2100']] if x != -1])), 1)    
        IFGW_2200_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2200'], row['IFGW_co_2200'], row['IFGW_ge_2200']] if x != -1])), 1)
        IFGW_2200_std = round(np.std(np.array([x for x in [row['IFGW_me_2200'], row['IFGW_co_2200'], row['IFGW_ge_2200']] if x != -1])), 1)    
        IFGW_2300_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2300'], row['IFGW_co_2300'], row['IFGW_ge_2300']] if x != -1])), 1)
        IFGW_2300_std = round(np.std(np.array([x for x in [row['IFGW_me_2300'], row['IFGW_co_2300'], row['IFGW_ge_2300']] if x != -1])), 1)    
        IFGW_2400_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2400'], row['IFGW_co_2400'], row['IFGW_ge_2400']] if x != -1])), 1)
        IFGW_2400_std = round(np.std(np.array([x for x in [row['IFGW_me_2400'], row['IFGW_co_2400'], row['IFGW_ge_2400']] if x != -1])), 1)        
        IFGW_2500_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2500'], row['IFGW_co_2500'], row['IFGW_ge_2500']] if x != -1])), 1)
        IFGW_2500_std = round(np.std(np.array([x for x in [row['IFGW_me_2500'], row['IFGW_co_2500'], row['IFGW_ge_2500']] if x != -1])), 1)        

        #IFGW_2000_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2000'], row['IFGW_co_2000'], row['IFGW_ge_2000']] if x != -1])), 1)
        #IFGW_2000_std = round(np.std(np.array([x for x in [row['IFGW_me_2000'], row['IFGW_co_2000'], row['IFGW_ge_2000']] if x != -1])), 1)    
        #IFGW_2050_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2050'], row['IFGW_co_2050'], row['IFGW_ge_2050']] if x != -1])), 1)
        #IFGW_2050_std = round(np.std(np.array([x for x in [row['IFGW_me_2050'], row['IFGW_co_2050'], row['IFGW_ge_2050']] if x != -1])), 1)    
        #IFGW_2100_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2100'], row['IFGW_co_2100'], row['IFGW_ge_2100']] if x != -1])), 1)
        #IFGW_2100_std = round(np.std(np.array([x for x in [row['IFGW_me_2100'], row['IFGW_co_2100'], row['IFGW_ge_2100']] if x != -1])), 1)    
        #IFGW_2200_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2150'], row['IFGW_co_2150'], row['IFGW_ge_2150']] if x != -1])), 1)
        #IFGW_2200_std = round(np.std(np.array([x for x in [row['IFGW_me_2150'], row['IFGW_co_2150'], row['IFGW_ge_2150']] if x != -1])), 1)    
        #IFGW_2300_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2200'], row['IFGW_co_2200'], row['IFGW_ge_2200']] if x != -1])), 1)
        #IFGW_2300_std = round(np.std(np.array([x for x in [row['IFGW_me_2200'], row['IFGW_co_2200'], row['IFGW_ge_2200']] if x != -1])), 1)    
        #IFGW_2400_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2250'], row['IFGW_co_2250'], row['IFGW_ge_2250']] if x != -1])), 1)
        #IFGW_2400_std = round(np.std(np.array([x for x in [row['IFGW_me_2250'], row['IFGW_co_2250'], row['IFGW_ge_2250']] if x != -1])), 1)        
        #IFGW_2500_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2300'], row['IFGW_co_2300'], row['IFGW_ge_2300']] if x != -1])), 1)
        #IFGW_2500_std = round(np.std(np.array([x for x in [row['IFGW_me_2300'], row['IFGW_co_2300'], row['IFGW_ge_2300']] if x != -1])), 1)        
    
        f = open(main_dir + rcp + '_mu_std.csv', 'a')
        f.write(str(fid) + ',' + str(cs_id) + ',' + str(srm_id) + ',' + cs_srm_id + ',' + str(inl_ext) + ',' + str(IFGW_2000_mu) + ',' + str(IFGW_2000_std)\
                 + ',' + str(IFGW_2050_mu) + ',' + str(IFGW_2050_std) + ',' + str(IFGW_2100_mu) + ',' + str(IFGW_2100_std)\
                 + ',' + str(IFGW_2200_mu) + ',' + str(IFGW_2200_std) + ',' + str(IFGW_2300_mu) + ',' + str(IFGW_2300_std) + ',' + str(IFGW_2400_mu) + ',' + str(IFGW_2400_std)\
                 + ',' + str(IFGW_2500_mu) + ',' + str(IFGW_2500_std))
        f.write("\n")
        f.close()                  
    
        fid += 1
    
#   ------------------------------------------------------------------------------------
#   Also create csv files showing difference between future predictions and current state
#   ------------------------------------------------------------------------------------
        
out_cols = ['fid', 'coscat', 'srm', 'cs_srm_id', 'inl_ext_km',\
            'IFGW_2000_DIFF_mu', 'IFGW_2000_DIFF_std', 'IFGW_2050_DIFF_mu', 'IFGW_2050_DIFF_std',\
            'IFGW_2100_DIFF_mu', 'IFGW_2100_DIFF_std', 'IFGW_2200_DIFF_mu', 'IFGW_2200_DIFF_std',\
            'IFGW_2300_DIFF_mu', 'IFGW_2300_DIFF_std', 'IFGW_2400_DIFF_mu', 'IFGW_2400_DIFF_std',\
            'IFGW_2500_DIFF_mu', 'IFGW_2500_DIFF_std']    
    
for rcp in ['26', '45', '85']:

    if not os.path.exists(main_dir + rcp + '_mu_std_DIFF_2000.csv'):
        df = pd.DataFrame(columns = out_cols)
        df.set_index('fid')
        df.to_csv(main_dir + rcp + '_mu_std_DIFF_2000.csv' , encoding = 'utf-8', index = False)
    
    rcp_in_dir = main_dir + rcp + '_results.csv'
    rcp_in = pd.read_csv(rcp_in_dir)
    fid = 1
    
    for index, row in rcp_in.iterrows():
        
        cs_id = row['coscat']
        srm_id = row['srm']
        cs_srm_id = row['cs_srm_id']
        inl_ext = row['inl_ext_km']
        
        IFGW_2000_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2000'], row['IFGW_co_2000'], row['IFGW_ge_2000']] if x != -1])), 1)
        IFGW_2000_std = round(np.std(np.array([x for x in [row['IFGW_me_2000'], row['IFGW_co_2000'], row['IFGW_ge_2000']] if x != -1])), 1)
        
        IFGW_2050_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2050'], row['IFGW_co_2050'], row['IFGW_ge_2050']] if x != -1])) - IFGW_2000_mu, 1)
        IFGW_2050_std = round(np.std(np.array([x for x in [row['IFGW_me_2050'], row['IFGW_co_2050'], row['IFGW_ge_2050']] if x != -1])), 1)    
        IFGW_2100_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2100'], row['IFGW_co_2100'], row['IFGW_ge_2100']] if x != -1])) - IFGW_2000_mu, 1)
        IFGW_2100_std = round(np.std(np.array([x for x in [row['IFGW_me_2100'], row['IFGW_co_2100'], row['IFGW_ge_2100']] if x != -1])), 1)    
        IFGW_2200_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2200'], row['IFGW_co_2200'], row['IFGW_ge_2200']] if x != -1])) - IFGW_2000_mu, 1)
        IFGW_2200_std = round(np.std(np.array([x for x in [row['IFGW_me_2200'], row['IFGW_co_2200'], row['IFGW_ge_2200']] if x != -1])), 1)    
        IFGW_2300_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2300'], row['IFGW_co_2300'], row['IFGW_ge_2300']] if x != -1])) - IFGW_2000_mu, 1)
        IFGW_2300_std = round(np.std(np.array([x for x in [row['IFGW_me_2300'], row['IFGW_co_2300'], row['IFGW_ge_2300']] if x != -1])), 1)    
        IFGW_2400_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2400'], row['IFGW_co_2400'], row['IFGW_ge_2400']] if x != -1])) - IFGW_2000_mu, 1)
        IFGW_2400_std = round(np.std(np.array([x for x in [row['IFGW_me_2400'], row['IFGW_co_2400'], row['IFGW_ge_2400']] if x != -1])), 1)        
        IFGW_2500_mu = round(np.mean(np.array([x for x in [row['IFGW_me_2500'], row['IFGW_co_2500'], row['IFGW_ge_2500']] if x != -1])) - IFGW_2000_mu, 1)
        IFGW_2500_std = round(np.std(np.array([x for x in [row['IFGW_me_2500'], row['IFGW_co_2500'], row['IFGW_ge_2500']] if x != -1])), 1)        
    
        f = open(main_dir + rcp + '_mu_std_DIFF_2000.csv', 'a')
        f.write(str(fid) + ',' + str(cs_id) + ',' + str(srm_id) + ',' + cs_srm_id + ',' + str(inl_ext) + ',' + str(0) + ',' + str(0)\
                 + ',' + str(IFGW_2050_mu) + ',' + str(IFGW_2050_std) + ',' + str(IFGW_2100_mu) + ',' + str(IFGW_2100_std)\
                 + ',' + str(IFGW_2200_mu) + ',' + str(IFGW_2200_std) + ',' + str(IFGW_2300_mu) + ',' + str(IFGW_2300_std) + ',' + str(IFGW_2400_mu) + ',' + str(IFGW_2400_std)\
                 + ',' + str(IFGW_2500_mu) + ',' + str(IFGW_2500_std))
        f.write("\n")
        f.close()                  
    
        fid += 1
    
"""    ---------------------------------------------------------------    """

import matplotlib.pyplot as plt
import math
import os
import matplotlib

import os
os.environ["PROJ_LIB"] = r'c:\Miniconda3\pkgs\proj4-5.2.0-h6538335_1006\Library\share'

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
from matplotlib.patches import PathPatch
from matplotlib.colors import Normalize
from matplotlib import colors
from matplotlib import colorbar
import matplotlib.cm as cmx
import matplotlib.patches as mpatches
import rasterio

#   ------------------------------------------------------------------------------------
#   Create final plots (worldmaps) for each RCP and year  
#   ------------------------------------------------------------------------------------

srm_ifgw_diff_dir = os.path.join(data_out_dir, '_SRM_IFGW_diff_2000.shp') #r'g:\Water_Nexus\_A4_paper\_data_output\_SRM_IFGW_diff_2000.shp'
basemap_raster_bw = os.path.join(data_out_dir, 'background_bw.tif') #r'g:\Water_Nexus\_A4_paper\_data_output\background_bw.tif'
basemap_raster_col = os.path.join(data_out_dir, 'background_col.tif')#r'g:\_ORIGINAL_DATA\natural_earth\NE2_LR_LC_SR\NE2_LR_LC_SR.tif'


#from mpl_toolkits.basemap import Basemap
#import matplotlib.pyplot as plt
import numpy as np
#from itertools import chain
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import pandas as pd

ts_name_lst = [['RCP_26_I_1', 'RCP_26_2050_diff_ifgw'], ['RCP_26_I_2', 'RCP_26_2100_diff_ifgw'], ['RCP_26_I_3', 'RCP_26_2200_diff_ifgw'],\
               ['RCP_26_I_4', 'RCP_26_2300_diff_ifgw'], ['RCP_26_I_5', 'RCP_26_2400_diff_ifgw'], ['RCP_26_I_6', 'RCP_26_2500_diff_ifgw'],\
               ['RCP_45_I_1', 'RCP_45_2050_diff_ifgw'], ['RCP_45_I_2', 'RCP_45_2100_diff_ifgw'], ['RCP_45_I_3', 'RCP_45_2200_diff_ifgw'],\
               ['RCP_45_I_4', 'RCP_45_2300_diff_ifgw'], ['RCP_45_I_5', 'RCP_45_2400_diff_ifgw'], ['RCP_45_I_6', 'RCP_45_2500_diff_ifgw'],\
               ['RCP_85_I_1', 'RCP_85_2050_diff_ifgw'], ['RCP_85_I_2', 'RCP_85_2100_diff_ifgw'], ['RCP_85_I_3', 'RCP_85_2200_diff_ifgw'],\
               ['RCP_85_I_4', 'RCP_85_2300_diff_ifgw'], ['RCP_85_I_5', 'RCP_85_2400_diff_ifgw'], ['RCP_85_I_6', 'RCP_85_2500_diff_ifgw']]

ts_pct_lst = []

"""
for i in range(len(ts_name_lst)):
    print(ts_name_lst[i][1])
    ts_cnt = [0, 0, 0, 0, 0]
    
    for shapedict,shape in zip(m.SRM_IFGW_diff_info, m.SRM_IFGW_diff):
        ifgw = shapedict[ts_name_lst[i][0]]
        xx, yy = zip(*shape)
        # show part of track where storm > Cat 4 as thick red.
        if ifgw is not None:
            if ifgw > -5.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'forestgreen')
                ts_cnt[0] += 1
            elif ifgw < -5.0 and ifgw > -10.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'gold')
                ts_cnt[1] += 1
            elif ifgw < -10.0 and ifgw > -25.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'hotpink')
                ts_cnt[2] += 1
            elif ifgw < -25.0 and ifgw > -50.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'red')
                ts_cnt[3] += 1
            elif ifgw < -50.0 and ifgw > -100.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'darkred')
                ts_cnt[4] += 1
        else:
            m.plot(xx, yy, linewidth = .1, color = 'k')
            
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    ts_pct_lst.append(ts_pct)
"""

for i in range(len(ts_name_lst)):

    #   plot the basemap and the results
    fig = plt.figure(figsize=(8 , 4), edgecolor='w')
    ax = fig.add_subplot(111)
    
    """
    # A good LCC projection for USA plots
    m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 75,\
                llcrnrlon = -140, urcrnrlon = 175, lat_ts = 20, resolution = 'i')
    
    m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
    m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
    """
    
    # A good LCC projection for USA plots
    m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = 0, urcrnrlat = 60,\
                llcrnrlon = -25, urcrnrlon = 95, lat_ts = 20, resolution = 'i')
    m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
    m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)

    m.drawcountries(linewidth = 0.1, linestyle = 'solid', color = 'grey') 
    m.drawcoastlines(linewidth = 0.1, linestyle = 'solid', color = 'grey')
    m.shadedrelief(scale = 0.2, alpha = 0.25)
    #m.etopo(scale = 0.2, alpha = 0.5)
    #   read in the shapefile
    m.readshapefile(os.path.join(data_out_dir, '_SRM_IFGW_diff_2000'), 'SRM_IFGW_diff')
    
    print(ts_name_lst[i][1])
    ts_cnt = [0, 0, 0, 0, 0]
    
    for shapedict,shape in zip(m.SRM_IFGW_diff_info, m.SRM_IFGW_diff):
        ifgw = shapedict[ts_name_lst[i][0]]
        xx, yy = zip(*shape)
        # show part of track where storm > Cat 4 as thick red.
        if ifgw is not None:
            if ifgw > -5.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'forestgreen')
                ts_cnt[0] += 1
            elif ifgw < -5.0 and ifgw > -10.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'gold')
                ts_cnt[1] += 1
            elif ifgw < -10.0 and ifgw > -25.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'hotpink')
                ts_cnt[2] += 1
            elif ifgw < -25.0 and ifgw > -50.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'red')
                ts_cnt[3] += 1
            elif ifgw < -50.0 and ifgw > -100.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'darkred')
                ts_cnt[4] += 1
        else:
            m.plot(xx, yy, linewidth = .1, color = 'k')
            
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    ts_pct_lst.append(ts_pct)
            
    plt.savefig(os.path.join(data_out_dir, ts_name_lst[i][1] + '_zoom_300dpi.png'), dpi = 300, facecolor = 'w', edgecolor = 'w',
        orientation = 'portrait', papertype = None, format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()




"""  -------------------------------------------------------------------------
            FIGURE WITH 3 WORLD MAPS FOR THE MAIN ARTICLE
--------------------------------------------------------------------------  """

ts_cnt = [0, 0, 0, 0, 0]

fig = plt.figure(figsize=(8 , 8), edgecolor='w') 
                                                                                                                                                                                             
ax1_1 = plt.subplot2grid((4, 1), (0, 0), fig = fig)             
ax1_2 = plt.subplot2grid((4, 1), (1, 0), fig = fig)             
ax1_3 = plt.subplot2grid((4, 1), (2, 0), fig = fig)             
ax1_4 = plt.subplot2grid((4, 1), (3, 0), fig = fig)     
        
ax1_1.set_position([0.025, 0.7, 0.95, 0.295]) # [left, bottom, width, height]                     
ax1_2.set_position([0.025, 0.4, 0.95, 0.295])
ax1_3.set_position([0.025, 0.1, 0.95, 0.295])  
ax1_4.set_position([0.25, 0.025, 0.5, 0.05])  

m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 75,\
            llcrnrlon = -140, urcrnrlon = 175, lat_ts = 20, resolution = 'i', ax = ax1_1)
m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawcountries(linewidth = 0.1, linestyle = 'solid', color = 'grey') 
m.drawcoastlines(linewidth = 0.1, linestyle = 'solid', color = 'grey')
m.shadedrelief(scale = 0.2, alpha = 0.25)
#m.etopo(scale = 0.2, alpha = 0.5)
#   read in the shapefile
m.readshapefile(os.path.join(data_out_dir, '_SRM_IFGW_diff_2000'), 'SRM_IFGW_diff')
for shapedict,shape in zip(m.SRM_IFGW_diff_info, m.SRM_IFGW_diff):
    ifgw = shapedict[ts_name_lst[1][0]]
    xx, yy = zip(*shape)
    # show part of track where storm > Cat 4 as thick red.
    if ifgw is not None:
        if ifgw > -5.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'forestgreen')
            ts_cnt[0] += 1
        elif ifgw < -5.0 and ifgw > -10.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'gold')
            ts_cnt[1] += 1
        elif ifgw < -10.0 and ifgw > -25.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'hotpink')
            ts_cnt[2] += 1
        elif ifgw < -25.0 and ifgw > -50.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'red')
            ts_cnt[3] += 1
        elif ifgw < -50.0 and ifgw > -100.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'darkred')
            ts_cnt[4] += 1
    else:
        m.plot(xx, yy, linewidth = .1, color = 'k')

# A good LCC projection for USA plots
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 75,\
            llcrnrlon = -140, urcrnrlon = 175, lat_ts = 20, resolution = 'i', ax = ax1_2)
m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawcountries(linewidth = 0.1, linestyle = 'solid', color = 'grey') 
m.drawcoastlines(linewidth = 0.1, linestyle = 'solid', color = 'grey')
m.shadedrelief(scale = 0.2, alpha = 0.25)
#m.etopo(scale = 0.2, alpha = 0.5)
#   read in the shapefile
m.readshapefile(os.path.join(data_out_dir, '_SRM_IFGW_diff_2000'), 'SRM_IFGW_diff')
for shapedict,shape in zip(m.SRM_IFGW_diff_info, m.SRM_IFGW_diff):
    ifgw = shapedict[ts_name_lst[7][0]]
    xx, yy = zip(*shape)
    # show part of track where storm > Cat 4 as thick red.
    if ifgw is not None:
        if ifgw > -5.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'forestgreen')
            ts_cnt[0] += 1
        elif ifgw < -5.0 and ifgw > -10.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'gold')
            ts_cnt[1] += 1
        elif ifgw < -10.0 and ifgw > -25.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'hotpink')
            ts_cnt[2] += 1
        elif ifgw < -25.0 and ifgw > -50.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'red')
            ts_cnt[3] += 1
        elif ifgw < -50.0 and ifgw > -100.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'darkred')
            ts_cnt[4] += 1
    else:
        m.plot(xx, yy, linewidth = .1, color = 'k')

# A good LCC projection for USA plots
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 75,\
            llcrnrlon = -140, urcrnrlon = 175, lat_ts = 20, resolution = 'i', ax = ax1_3)
m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawcountries(linewidth = 0.1, linestyle = 'solid', color = 'grey') 
m.drawcoastlines(linewidth = 0.1, linestyle = 'solid', color = 'grey')
m.shadedrelief(scale = 0.2, alpha = 0.25)
#m.etopo(scale = 0.2, alpha = 0.5)
#   read in the shapefile
m.readshapefile(os.path.join(data_out_dir, '_SRM_IFGW_diff_2000'), 'SRM_IFGW_diff')
for shapedict,shape in zip(m.SRM_IFGW_diff_info, m.SRM_IFGW_diff):
    ifgw = shapedict[ts_name_lst[13][0]]
    xx, yy = zip(*shape)
    # show part of track where storm > Cat 4 as thick red.
    if ifgw is not None:
        if ifgw > -5.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'forestgreen', label = '< -5%')
            ts_cnt[0] += 1
        elif ifgw < -5.0 and ifgw > -10.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'gold', label = '-5% < AND < -10%')
            ts_cnt[1] += 1
        elif ifgw < -10.0 and ifgw > -25.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'hotpink', label = '-10% < AND < -25%')
            ts_cnt[2] += 1
        elif ifgw < -25.0 and ifgw > -50.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'red', label = '-25% < AND < -50%')
            ts_cnt[3] += 1
        elif ifgw < -50.0 and ifgw > -100.0:
            m.plot(xx, yy, linewidth = 1.25, color = 'darkred', label = '> -50%')
            ts_cnt[4] += 1
    else:
        m.plot(xx, yy, linewidth = .1, color = 'k')

h, l = ax1_3.get_legend_handles_labels() 
handle_list, label_list = [], []
for h, l in zip(h, l):
    if l not in label_list:
        handle_list.append(h)
        label_list.append(l)

ax1_4.legend(handle_list, label_list, fontsize = 7, frameon = False, ncol = 5, loc = 'lower center', title = "IFGV difference compared to 2000 [%]", title_fontsize = 8)
ax1_4.axis('off')

canvas = FigureCanvas(fig)
canvas.print_figure(os.path.join(data_out_dir, '_A4_figure_v2.png'), dpi=300)


"""  -------------------------------------------------------------------------
    HISTOGRAM WITH AVERAGE IFGW DIFFERENCE TO 2000, ALL DEMS COMBINED
--------------------------------------------------------------------------  """

fig = plt.figure(figsize=(8 , 3), edgecolor='w') 
                                                                                                                                                                                             
ax1_1 = plt.subplot2grid((3, 3), (0, 0), fig = fig)             
ax1_2 = plt.subplot2grid((3, 3), (0, 1), fig = fig)             
ax1_3 = plt.subplot2grid((3, 3), (0, 2), fig = fig)             
ax2_1 = plt.subplot2grid((3, 3), (1, 0), fig = fig)             
ax2_2 = plt.subplot2grid((3, 3), (1, 1), fig = fig)             
ax2_3 = plt.subplot2grid((3, 3), (1, 2), fig = fig)     
ax3 = plt.subplot2grid((3, 3), (2, 0), colspan = 3, fig = fig)     

ax1_1.set_position([0.1, 0.925, 0.266, 0.05])                      
ax1_2.set_position([0.366, 0.925, 0.266, 0.05])
ax1_3.set_position([0.632, 0.925, 0.266, 0.05])  
ax2_1.set_position([0.1, 0.25, 0.266, 0.65])  # [left, bottom, width, height]
ax2_2.set_position([0.366, 0.25, 0.266, 0.65])
ax2_3.set_position([0.632, 0.25, 0.266, 0.65])
ax3.set_position([0.1, 0., 0.8, 0.055])

txt = ax1_1.text(0.4, 0.5, 'RCP 26', fontsize = 7, fontweight = 'bold', va = 'center')
txt.set_clip_on(False)
ax1_1.axis('off')   
txt = ax1_2.text(0.4, 0.5, 'RCP 45', fontsize = 7, fontweight = 'bold', va = 'center')
txt.set_clip_on(False)
ax1_2.axis('off')   
txt = ax1_3.text(0.4, 0.5, 'RCP 85', fontsize = 7, fontweight = 'bold', va = 'center')
txt.set_clip_on(False)
ax1_3.axis('off')   

rcp_26 = ts_pct_lst[:6]
# Data
r = [0, 1, 2, 3, 4, 5]
raw_data = {'greenBars': [round(i[0], 4) for i in rcp_26], 'goldBars': [round(i[1], 4) for i in rcp_26], 'pinkBars': [round(i[2], 4) for i in rcp_26],\
            'redBars': [round(i[3], 4) for i in rcp_26], 'darkredBars': [round(i[4], 4) for i in rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()

# plot
barWidth = 0.85
names = ('2050','2100','2200','2300','2400', '2500')
# Create green Bars
ax2_1.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_1.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_1.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_1.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_1.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)
# Custom x axis
ax2_1.set_xticks(r)
ax2_1.set_xticklabels(names, fontsize = 7)
ax2_1.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_1.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax2_1.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax2_1.set_xlabel("Year", fontsize = 7)
ax2_1.set_ylim(1, 0)

rcp_45 = ts_pct_lst[6:12]
# Data
r = [0, 1, 2, 3, 4, 5]
raw_data = {'greenBars': [round(i[0], 4) for i in rcp_45], 'goldBars': [round(i[1], 4) for i in rcp_45], 'pinkBars': [round(i[2], 4) for i in rcp_45],\
            'redBars': [round(i[3], 4) for i in rcp_45], 'darkredBars': [round(i[4], 4) for i in rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# plot
barWidth = 0.85
# Create green Bars
ax2_2.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)
# Custom x axis
ax2_2.set_xticks(r)
ax2_2.set_xticklabels(names, fontsize = 7)
ax2_2.get_yaxis().set_ticks([])
ax2_2.set_xlabel("Year", fontsize = 7)
ax2_2.set_ylim(1, 0)

rcp_85 = ts_pct_lst[12:]
# Data
r = [0, 1, 2, 3, 4, 5]
raw_data = {'greenBars': [round(i[0], 4) for i in rcp_85], 'goldBars': [round(i[1], 4) for i in rcp_85], 'pinkBars': [round(i[2], 4) for i in rcp_85],\
            'redBars': [round(i[3], 4) for i in rcp_85], 'darkredBars': [round(i[4], 4) for i in rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# plot
barWidth = 0.85
# Create green Bars
ax2_3.bar(r, greenBars, color = 'forestgreen', width = barWidth, label = '< -5%')
ax2_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth, label = '-5% < AND < -10%')
ax2_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth, label = '-10% < AND < -25%')
ax2_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth, label = '-25% < AND < -50%')
ax2_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth, label = '> -50%')
# Custom x axis
ax2_3.set_xticks(r)
ax2_3.set_xticklabels(names, fontsize = 7)
ax2_3.get_yaxis().set_ticks([])
ax2_3.set_ylim(1, 0)
ax2_3.set_xlabel("Year", fontsize = 7)
h, l = ax2_3.get_legend_handles_labels() 
ax3.legend(h, l, fontsize = 7, frameon = False, ncol = 5, loc = 'lower center', title = "IFGW difference compared to 2000 [%]", title_fontsize = 8)
ax3.axis('off')

canvas = FigureCanvas(fig)
canvas.print_figure(os.path.join(data_out_dir, 'hist_IFGW_diff_2000.png'), dpi=300)



"""  -------------------------------------------------------------------------
****** VERSION ONLY TILL 2300    HISTOGRAM WITH AVERAGE IFGW DIFFERENCE TO 2000, ALL DEMS COMBINED
--------------------------------------------------------------------------  """

fig = plt.figure(figsize=(8 , 3), edgecolor='w') 
                                                                                                                                                                                             
ax1_1 = plt.subplot2grid((3, 3), (0, 0), fig = fig)             
ax1_2 = plt.subplot2grid((3, 3), (0, 1), fig = fig)             
ax1_3 = plt.subplot2grid((3, 3), (0, 2), fig = fig)             
ax2_1 = plt.subplot2grid((3, 3), (1, 0), fig = fig)             
ax2_2 = plt.subplot2grid((3, 3), (1, 1), fig = fig)             
ax2_3 = plt.subplot2grid((3, 3), (1, 2), fig = fig)     
ax3 = plt.subplot2grid((3, 3), (2, 0), colspan = 3, fig = fig)     

ax1_1.set_position([0.1, 0.925, 0.266, 0.05])                      
ax1_2.set_position([0.366, 0.925, 0.266, 0.05])
ax1_3.set_position([0.632, 0.925, 0.266, 0.05])  
ax2_1.set_position([0.1, 0.25, 0.266, 0.65])  # [left, bottom, width, height]
ax2_2.set_position([0.366, 0.25, 0.266, 0.65])
ax2_3.set_position([0.632, 0.25, 0.266, 0.65])
ax3.set_position([0.1, 0., 0.8, 0.055])

txt = ax1_1.text(0.4, 0.5, 'RCP 26', fontsize = 7, fontweight = 'bold', va = 'center')
txt.set_clip_on(False)
ax1_1.axis('off')   
txt = ax1_2.text(0.4, 0.5, 'RCP 45', fontsize = 7, fontweight = 'bold', va = 'center')
txt.set_clip_on(False)
ax1_2.axis('off')   
txt = ax1_3.text(0.4, 0.5, 'RCP 85', fontsize = 7, fontweight = 'bold', va = 'center')
txt.set_clip_on(False)
ax1_3.axis('off')   

rcp_26 = ts_pct_lst[:4]
# Data
r = [0, 1, 2, 3]
raw_data = {'greenBars': [round(i[0], 4) for i in rcp_26], 'goldBars': [round(i[1], 4) for i in rcp_26], 'pinkBars': [round(i[2], 4) for i in rcp_26],\
            'redBars': [round(i[3], 4) for i in rcp_26], 'darkredBars': [round(i[4], 4) for i in rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()

# plot
barWidth = 0.85
names = ('2050','2100','2200','2300')
# Create green Bars
ax2_1.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_1.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_1.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_1.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_1.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)
# Custom x axis
ax2_1.set_xticks(r)
ax2_1.set_xticklabels(names, fontsize = 7)
ax2_1.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_1.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax2_1.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax2_1.set_xlabel("Year", fontsize = 7)
ax2_1.set_ylim(1, 0)

rcp_45 = ts_pct_lst[6:10]
# Data
r = [0, 1, 2, 3]
raw_data = {'greenBars': [round(i[0], 4) for i in rcp_45], 'goldBars': [round(i[1], 4) for i in rcp_45], 'pinkBars': [round(i[2], 4) for i in rcp_45],\
            'redBars': [round(i[3], 4) for i in rcp_45], 'darkredBars': [round(i[4], 4) for i in rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# plot
barWidth = 0.85
# Create green Bars
ax2_2.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)
# Custom x axis
ax2_2.set_xticks(r)
ax2_2.set_xticklabels(names, fontsize = 7)
ax2_2.get_yaxis().set_ticks([])
ax2_2.set_xlabel("Year", fontsize = 7)
ax2_2.set_ylim(1, 0)

rcp_85 = ts_pct_lst[12:16]
# Data
r = [0, 1, 2, 3]
raw_data = {'greenBars': [round(i[0], 4) for i in rcp_85], 'goldBars': [round(i[1], 4) for i in rcp_85], 'pinkBars': [round(i[2], 4) for i in rcp_85],\
            'redBars': [round(i[3], 4) for i in rcp_85], 'darkredBars': [round(i[4], 4) for i in rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# plot
barWidth = 0.85
# Create green Bars
ax2_3.bar(r, greenBars, color = 'forestgreen', width = barWidth, label = '< -5%')
ax2_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth, label = '-5% < AND < -10%')
ax2_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth, label = '-10% < AND < -25%')
ax2_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth, label = '-25% < AND < -50%')
ax2_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth, label = '> -50%')
# Custom x axis
ax2_3.set_xticks(r)
ax2_3.set_xticklabels(names, fontsize = 7)
ax2_3.get_yaxis().set_ticks([])
ax2_3.set_ylim(1, 0)
ax2_3.set_xlabel("Year", fontsize = 7)
h, l = ax2_3.get_legend_handles_labels() 
ax3.legend(h, l, fontsize = 7, frameon = False, ncol = 5, loc = 'lower center', title = "IFGW difference compared to 2000 [%]", title_fontsize = 8)
ax3.axis('off')

canvas = FigureCanvas(fig)
canvas.print_figure(os.path.join(data_out_dir, 'hist_IFGW_diff_2000_v2.png'), dpi=300)




"""  -------------------------------------------------------------------------
            Create the plots comparing DEM results for each RCP
--------------------------------------------------------------------------  """

import pandas as pd

#   read in the results and land use csv
lu_pop_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\_land_use_population_GDP.csv')
rcp_26_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_26_results.csv')
rcp_45_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_45_results.csv')
rcp_85_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_85_results.csv')

merit_rcp_26, coastal_rcp_26, gebco_rcp_26 = [], [], []
merit_rcp_45, coastal_rcp_45, gebco_rcp_45 = [], [], []
merit_rcp_85, coastal_rcp_85, gebco_rcp_85 = [], [], []

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300', 'IFGW_me_2400', 'IFGW_me_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the merit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_26['IFGW_me_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    merit_rcp_26.append(ts_pct)        
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_co_2050', 'IFGW_co_2100', 'IFGW_co_2200', 'IFGW_co_2300', 'IFGW_co_2400', 'IFGW_co_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the corit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_co = row_rcp_26['IFGW_co_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_co
            if val_2000_co is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    coastal_rcp_26.append(ts_pct)         

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_ge_2050', 'IFGW_ge_2100', 'IFGW_ge_2200', 'IFGW_ge_2300', 'IFGW_ge_2400', 'IFGW_ge_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the gerit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_ge = row_rcp_26['IFGW_ge_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_ge
            if val_2000_ge is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    gebco_rcp_26.append(ts_pct)      
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300', 'IFGW_me_2400', 'IFGW_me_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the merit list
            row_rcp_45 = rcp_45_csv.loc[rcp_45_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_45['IFGW_me_2000'].values[0]
            ifgw = row_rcp_45[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    merit_rcp_45.append(ts_pct)        
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_co_2050', 'IFGW_co_2100', 'IFGW_co_2200', 'IFGW_co_2300', 'IFGW_co_2400', 'IFGW_co_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the corit list
            row_rcp_45 = rcp_45_csv.loc[rcp_45_csv['cs_srm_id'] == row['id_srm']]
            val_2000_co = row_rcp_45['IFGW_co_2000'].values[0]
            ifgw = row_rcp_45[ts].values[0] - val_2000_co
            if val_2000_co is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    coastal_rcp_45.append(ts_pct)         

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_ge_2050', 'IFGW_ge_2100', 'IFGW_ge_2200', 'IFGW_ge_2300', 'IFGW_ge_2400', 'IFGW_ge_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the gerit list
            row_rcp_45 = rcp_45_csv.loc[rcp_45_csv['cs_srm_id'] == row['id_srm']]
            val_2000_ge = row_rcp_45['IFGW_ge_2000'].values[0]
            ifgw = row_rcp_45[ts].values[0] - val_2000_ge
            if val_2000_ge is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    gebco_rcp_45.append(ts_pct)              

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300', 'IFGW_me_2400', 'IFGW_me_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the merit list
            row_rcp_85 = rcp_85_csv.loc[rcp_85_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_85['IFGW_me_2000'].values[0]
            ifgw = row_rcp_85[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    merit_rcp_85.append(ts_pct)        
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_co_2050', 'IFGW_co_2100', 'IFGW_co_2200', 'IFGW_co_2300', 'IFGW_co_2400', 'IFGW_co_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the corit list
            row_rcp_85 = rcp_85_csv.loc[rcp_85_csv['cs_srm_id'] == row['id_srm']]
            val_2000_co = row_rcp_85['IFGW_co_2000'].values[0]
            ifgw = row_rcp_85[ts].values[0] - val_2000_co
            if val_2000_co is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    coastal_rcp_85.append(ts_pct)         

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_ge_2050', 'IFGW_ge_2100', 'IFGW_ge_2200', 'IFGW_ge_2300', 'IFGW_ge_2400', 'IFGW_ge_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the gerit list
            row_rcp_85 = rcp_85_csv.loc[rcp_85_csv['cs_srm_id'] == row['id_srm']]
            val_2000_ge = row_rcp_85['IFGW_ge_2000'].values[0]
            ifgw = row_rcp_85[ts].values[0] - val_2000_ge
            if val_2000_ge is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    gebco_rcp_85.append(ts_pct)      


fig = plt.figure(figsize=(8 , 8), edgecolor='w') 
#   column titles
ax1_1 = plt.subplot2grid((6, 5), (0, 1), fig = fig)             
ax1_2 = plt.subplot2grid((6, 5), (0, 2), fig = fig)             
ax1_3 = plt.subplot2grid((6, 5), (0, 3), fig = fig)
#   first row of plots
ax2_1 = plt.subplot2grid((6, 5), (1, 4), fig = fig)             
ax2_2 = plt.subplot2grid((6, 5), (1, 1), fig = fig)             
ax2_3 = plt.subplot2grid((6, 5), (1, 2), fig = fig)     
ax2_4 = plt.subplot2grid((6, 5), (1, 3), fig = fig) 
#   second row of plots
ax3_1 = plt.subplot2grid((6, 5), (2, 4), fig = fig)             
ax3_2 = plt.subplot2grid((6, 5), (2, 1), fig = fig)             
ax3_3 = plt.subplot2grid((6, 5), (2, 2), fig = fig)     
ax3_4 = plt.subplot2grid((6, 5), (2, 3), fig = fig) 
#   third row of plots
ax4_1 = plt.subplot2grid((6, 5), (3, 4), fig = fig)             
ax4_2 = plt.subplot2grid((6, 5), (3, 1), fig = fig)             
ax4_3 = plt.subplot2grid((6, 5), (3, 2), fig = fig)     
ax4_4 = plt.subplot2grid((6, 5), (3, 3), fig = fig) 

#   x axis title
#ax5_1 = plt.subplot2grid((6, 5), (4, 1), colspan = 3, fig = fig) 
#   legend
ax6_1 = plt.subplot2grid((6, 5), (5, 1), colspan = 3, fig = fig)     

#   define the dimensions of each plot 
ax1_1.set_position([0.075, 0.95, 0.275, 0.05])  # [left, bottom, width, height]                     
ax1_2.set_position([0.358, 0.95, 0.275, 0.05])
ax1_3.set_position([0.641, 0.95, 0.275, 0.05])  

ax2_1.set_position([0.925, 0.666, 0.05, 0.275])  # [left, bottom, width, height]                     
ax2_2.set_position([0.075, 0.666, 0.275, 0.275])
ax2_3.set_position([0.358, 0.666, 0.275, 0.275])  
ax2_4.set_position([0.641, 0.666, 0.275, 0.275])

ax3_1.set_position([0.925, 0.383, 0.05, 0.275])  # [left, bottom, width, height]                     
ax3_2.set_position([0.075, 0.383, 0.275, 0.275])
ax3_3.set_position([0.358, 0.383, 0.275, 0.275])  
ax3_4.set_position([0.641, 0.383, 0.275, 0.275])

ax4_1.set_position([0.925, 0.1, 0.05, 0.275])  # [left, bottom, width, height]                     
ax4_2.set_position([0.075, 0.1, 0.275, 0.275])
ax4_3.set_position([0.358, 0.1, 0.275, 0.275])  
ax4_4.set_position([0.641, 0.1, 0.275, 0.275])

#ax5_1.set_position([0.075, 0.055, 0.875, 0.025])  
ax6_1.set_position([0.075, 0.01, 0.875, 0.05])

#   define all the text plots
txt = ax1_1.text(0.5, 0.5, 'Coastal DEM', fontsize = 9, fontweight = 'bold', va = 'center', ha = 'center')
txt.set_clip_on(False)
ax1_1.axis('off')   
txt = ax1_2.text(0.5, 0.5, 'GEBCO', fontsize = 9, fontweight = 'bold', va = 'center', ha = 'center')
txt.set_clip_on(False)
ax1_2.axis('off')   
txt = ax1_3.text(0.5, 0.5, 'MERIT DEM', fontsize = 9, fontweight = 'bold', va = 'center', ha = 'center')
txt.set_clip_on(False)
ax1_3.axis('off')  

txt = ax2_1.text(0.5, 0.5, 'RCP 26', fontsize = 9, fontweight = 'bold', rotation=270, va = 'center')
txt.set_clip_on(False)
ax2_1.axis('off')   
txt = ax3_1.text(0.5, 0.5, 'RCP 45', fontsize = 9, fontweight = 'bold',rotation=270, va = 'center')
txt.set_clip_on(False)
ax3_1.axis('off')   
txt = ax4_1.text(0.5, 0.5, 'RCP 85', fontsize = 9, fontweight = 'bold', rotation=270, va = 'center')
txt.set_clip_on(False)
ax4_1.axis('off')   

barWidth = 0.85
r = [0, 1, 2, 3, 4, 5]
# Data
raw_data = {'greenBars': [round(i[0], 4) for i in coastal_rcp_26], 'goldBars': [round(i[1], 4) for i in coastal_rcp_26], 'pinkBars': [round(i[2], 4) for i in coastal_rcp_26],\
            'redBars': [round(i[3], 4) for i in coastal_rcp_26], 'darkredBars': [round(i[4], 4) for i in coastal_rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax2_2.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in gebco_rcp_26], 'goldBars': [round(i[1], 4) for i in gebco_rcp_26], 'pinkBars': [round(i[2], 4) for i in gebco_rcp_26],\
            'redBars': [round(i[3], 4) for i in gebco_rcp_26], 'darkredBars': [round(i[4], 4) for i in gebco_rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax2_3.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in merit_rcp_26], 'goldBars': [round(i[1], 4) for i in merit_rcp_26], 'pinkBars': [round(i[2], 4) for i in merit_rcp_26],\
            'redBars': [round(i[3], 4) for i in merit_rcp_26], 'darkredBars': [round(i[4], 4) for i in merit_rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax2_4.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_4.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_4.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_4.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_4.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in coastal_rcp_45], 'goldBars': [round(i[1], 4) for i in coastal_rcp_45], 'pinkBars': [round(i[2], 4) for i in coastal_rcp_45],\
            'redBars': [round(i[3], 4) for i in coastal_rcp_45], 'darkredBars': [round(i[4], 4) for i in coastal_rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax3_2.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax3_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax3_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax3_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax3_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in gebco_rcp_45], 'goldBars': [round(i[1], 4) for i in gebco_rcp_45], 'pinkBars': [round(i[2], 4) for i in gebco_rcp_45],\
            'redBars': [round(i[3], 4) for i in gebco_rcp_45], 'darkredBars': [round(i[4], 4) for i in gebco_rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax3_3.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax3_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax3_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax3_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax3_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in merit_rcp_45], 'goldBars': [round(i[1], 4) for i in merit_rcp_45], 'pinkBars': [round(i[2], 4) for i in merit_rcp_45],\
            'redBars': [round(i[3], 4) for i in merit_rcp_45], 'darkredBars': [round(i[4], 4) for i in merit_rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax3_4.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax3_4.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax3_4.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax3_4.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax3_4.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in coastal_rcp_85], 'goldBars': [round(i[1], 4) for i in coastal_rcp_85], 'pinkBars': [round(i[2], 4) for i in coastal_rcp_85],\
            'redBars': [round(i[3], 4) for i in coastal_rcp_85], 'darkredBars': [round(i[4], 4) for i in coastal_rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax4_2.bar(r, greenBars, color = 'forestgreen', width = barWidth, label = '< -5%')
ax4_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth, label = '-5% < AND < -10%')
ax4_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth, label = '-10% < AND < -25%')
ax4_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth, label = '-25% < AND < -50%')
ax4_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth, label = '> -50%')

raw_data = {'greenBars': [round(i[0], 4) for i in gebco_rcp_85], 'goldBars': [round(i[1], 4) for i in gebco_rcp_85], 'pinkBars': [round(i[2], 4) for i in gebco_rcp_85],\
            'redBars': [round(i[3], 4) for i in gebco_rcp_85], 'darkredBars': [round(i[4], 4) for i in gebco_rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax4_3.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax4_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax4_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax4_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax4_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in merit_rcp_85], 'goldBars': [round(i[1], 4) for i in merit_rcp_85], 'pinkBars': [round(i[2], 4) for i in merit_rcp_85],\
            'redBars': [round(i[3], 4) for i in merit_rcp_85], 'darkredBars': [round(i[4], 4) for i in merit_rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax4_4.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax4_4.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax4_4.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax4_4.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax4_4.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

# Custom x axis
names = ('2050','2100','2200','2300','2400', '2500')
ax2_2.set_xticks(r)
ax2_2.set_xticklabels([])
ax2_2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_2.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax2_2.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax2_2.set_ylim(1, 0)

ax2_3.set_xticks(r)
ax2_3.set_xticklabels([])
ax2_3.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_3.set_yticklabels([])
ax2_3.set_ylim(1, 0)

ax2_4.set_xticks(r)
ax2_4.set_xticklabels([])
ax2_4.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_4.set_yticklabels([])
ax2_4.set_ylim(1, 0)

# Custom x axis
ax3_2.set_xticks(r)
ax3_2.set_xticklabels([])
ax3_2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax3_2.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax3_2.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax3_2.set_ylim(1, 0)

ax3_3.set_xticks(r)
ax3_3.set_xticklabels([])
ax3_3.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax3_3.set_yticklabels([])
ax3_3.set_ylim(1, 0)

ax3_4.set_xticks(r)
ax3_4.set_xticklabels([])
ax3_4.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax3_4.set_yticklabels([])
ax3_4.set_ylim(1, 0)

# Custom x axis
ax4_2.set_xticks(r)
ax4_2.set_xticklabels(names, fontsize = 7)
ax4_2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax4_2.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax4_2.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax4_2.set_ylim(1, 0)
ax4_2.set_xlabel("Year", fontsize = 7)

ax4_3.set_xticks(r)
ax4_3.set_xticklabels(names, fontsize = 7)
ax4_3.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax4_3.set_yticklabels([])
ax4_3.set_ylim(1, 0)
ax4_3.set_xlabel("Year", fontsize = 7)

ax4_4.set_xticks(r)
ax4_4.set_xticklabels(names, fontsize = 7)
ax4_4.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax4_4.set_yticklabels([])
ax4_4.set_ylim(1, 0)
ax4_4.set_xlabel("Year", fontsize = 7)

h, l = ax4_2.get_legend_handles_labels() 
ax6_1.legend(h, l, fontsize = 7, frameon = False, ncol = 5, loc = 'lower center', title = "IFGW difference compared to 2000 [%]", title_fontsize = 8)
ax6_1.axis('off')

canvas = FigureCanvas(fig)
canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_urban_diff_2000.png'), dpi=300)
#canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_agriculture_diff_2000.png'), dpi=300)
#canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_nature_diff_2000.png'), dpi=300)
#canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_ALL_diff_2000.png'), dpi=300)






"""  -------------------------------------------------------------------------
**********NEWEST VERSIOn till 2300            Create the plots comparing DEM results for each RCP
--------------------------------------------------------------------------  """

import pandas as pd

#   read in the results and land use csv
lu_pop_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\_land_use_population_GDP.csv')
rcp_26_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_26_results.csv')
rcp_45_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_45_results.csv')
rcp_85_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_85_results.csv')

merit_rcp_26, coastal_rcp_26, gebco_rcp_26 = [], [], []
merit_rcp_45, coastal_rcp_45, gebco_rcp_45 = [], [], []
merit_rcp_85, coastal_rcp_85, gebco_rcp_85 = [], [], []

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        #if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        if row['nature'] > 0.:
            #   open the rcp result list and then append to the merit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_26['IFGW_me_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    merit_rcp_26.append(ts_pct)        
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_co_2050', 'IFGW_co_2100', 'IFGW_co_2200', 'IFGW_co_2300']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        #if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        if row['nature'] > 0.:
            #   open the rcp result list and then append to the corit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_co = row_rcp_26['IFGW_co_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_co
            if val_2000_co is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    coastal_rcp_26.append(ts_pct)         

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_ge_2050', 'IFGW_ge_2100', 'IFGW_ge_2200', 'IFGW_ge_2300']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        #if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        if row['nature'] > 0.:
            #   open the rcp result list and then append to the gerit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_ge = row_rcp_26['IFGW_ge_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_ge
            if val_2000_ge is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    gebco_rcp_26.append(ts_pct)      
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        #if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        if row['nature'] > 0.:
            #   open the rcp result list and then append to the merit list
            row_rcp_45 = rcp_45_csv.loc[rcp_45_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_45['IFGW_me_2000'].values[0]
            ifgw = row_rcp_45[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    merit_rcp_45.append(ts_pct)        
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_co_2050', 'IFGW_co_2100', 'IFGW_co_2200', 'IFGW_co_2300']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        #if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        if row['nature'] > 0.:
            #   open the rcp result list and then append to the corit list
            row_rcp_45 = rcp_45_csv.loc[rcp_45_csv['cs_srm_id'] == row['id_srm']]
            val_2000_co = row_rcp_45['IFGW_co_2000'].values[0]
            ifgw = row_rcp_45[ts].values[0] - val_2000_co
            if val_2000_co is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    coastal_rcp_45.append(ts_pct)         

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_ge_2050', 'IFGW_ge_2100', 'IFGW_ge_2200', 'IFGW_ge_2300']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        #if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        if row['nature'] > 0.:
            #   open the rcp result list and then append to the gerit list
            row_rcp_45 = rcp_45_csv.loc[rcp_45_csv['cs_srm_id'] == row['id_srm']]
            val_2000_ge = row_rcp_45['IFGW_ge_2000'].values[0]
            ifgw = row_rcp_45[ts].values[0] - val_2000_ge
            if val_2000_ge is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    gebco_rcp_45.append(ts_pct)              

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        #if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        if row['nature'] > 0.:
            #   open the rcp result list and then append to the merit list
            row_rcp_85 = rcp_85_csv.loc[rcp_85_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_85['IFGW_me_2000'].values[0]
            ifgw = row_rcp_85[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    merit_rcp_85.append(ts_pct)        
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_co_2050', 'IFGW_co_2100', 'IFGW_co_2200', 'IFGW_co_2300']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        #if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        if row['nature'] > 0.:
            #   open the rcp result list and then append to the corit list
            row_rcp_85 = rcp_85_csv.loc[rcp_85_csv['cs_srm_id'] == row['id_srm']]
            val_2000_co = row_rcp_85['IFGW_co_2000'].values[0]
            ifgw = row_rcp_85[ts].values[0] - val_2000_co
            if val_2000_co is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    coastal_rcp_85.append(ts_pct)         

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_ge_2050', 'IFGW_ge_2100', 'IFGW_ge_2200', 'IFGW_ge_2300']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        #if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        if row['nature'] > 0.:
            #   open the rcp result list and then append to the gerit list
            row_rcp_85 = rcp_85_csv.loc[rcp_85_csv['cs_srm_id'] == row['id_srm']]
            val_2000_ge = row_rcp_85['IFGW_ge_2000'].values[0]
            ifgw = row_rcp_85[ts].values[0] - val_2000_ge
            if val_2000_ge is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    gebco_rcp_85.append(ts_pct)      


fig = plt.figure(figsize=(8 , 8), edgecolor='w') 
#   column titles
ax1_1 = plt.subplot2grid((6, 5), (0, 1), fig = fig)             
ax1_2 = plt.subplot2grid((6, 5), (0, 2), fig = fig)             
ax1_3 = plt.subplot2grid((6, 5), (0, 3), fig = fig)
#   first row of plots
ax2_1 = plt.subplot2grid((6, 5), (1, 4), fig = fig)             
ax2_2 = plt.subplot2grid((6, 5), (1, 1), fig = fig)             
ax2_3 = plt.subplot2grid((6, 5), (1, 2), fig = fig)     
ax2_4 = plt.subplot2grid((6, 5), (1, 3), fig = fig) 
#   second row of plots
ax3_1 = plt.subplot2grid((6, 5), (2, 4), fig = fig)             
ax3_2 = plt.subplot2grid((6, 5), (2, 1), fig = fig)             
ax3_3 = plt.subplot2grid((6, 5), (2, 2), fig = fig)     
ax3_4 = plt.subplot2grid((6, 5), (2, 3), fig = fig) 
#   third row of plots
ax4_1 = plt.subplot2grid((6, 5), (3, 4), fig = fig)             
ax4_2 = plt.subplot2grid((6, 5), (3, 1), fig = fig)             
ax4_3 = plt.subplot2grid((6, 5), (3, 2), fig = fig)     
ax4_4 = plt.subplot2grid((6, 5), (3, 3), fig = fig) 

#   x axis title
#ax5_1 = plt.subplot2grid((6, 5), (4, 1), colspan = 3, fig = fig) 
#   legend
ax6_1 = plt.subplot2grid((6, 5), (5, 1), colspan = 3, fig = fig)     

#   define the dimensions of each plot 
ax1_1.set_position([0.075, 0.95, 0.275, 0.05])  # [left, bottom, width, height]                     
ax1_2.set_position([0.358, 0.95, 0.275, 0.05])
ax1_3.set_position([0.641, 0.95, 0.275, 0.05])  

ax2_1.set_position([0.925, 0.666, 0.05, 0.275])  # [left, bottom, width, height]                     
ax2_2.set_position([0.075, 0.666, 0.275, 0.275])
ax2_3.set_position([0.358, 0.666, 0.275, 0.275])  
ax2_4.set_position([0.641, 0.666, 0.275, 0.275])

ax3_1.set_position([0.925, 0.383, 0.05, 0.275])  # [left, bottom, width, height]                     
ax3_2.set_position([0.075, 0.383, 0.275, 0.275])
ax3_3.set_position([0.358, 0.383, 0.275, 0.275])  
ax3_4.set_position([0.641, 0.383, 0.275, 0.275])

ax4_1.set_position([0.925, 0.1, 0.05, 0.275])  # [left, bottom, width, height]                     
ax4_2.set_position([0.075, 0.1, 0.275, 0.275])
ax4_3.set_position([0.358, 0.1, 0.275, 0.275])  
ax4_4.set_position([0.641, 0.1, 0.275, 0.275])

#ax5_1.set_position([0.075, 0.055, 0.875, 0.025])  
ax6_1.set_position([0.075, 0.01, 0.875, 0.05])

#   define all the text plots
txt = ax1_1.text(0.5, 0.5, 'Coastal DEM', fontsize = 9, fontweight = 'bold', va = 'center', ha = 'center')
txt.set_clip_on(False)
ax1_1.axis('off')   
txt = ax1_2.text(0.5, 0.5, 'GEBCO', fontsize = 9, fontweight = 'bold', va = 'center', ha = 'center')
txt.set_clip_on(False)
ax1_2.axis('off')   
txt = ax1_3.text(0.5, 0.5, 'MERIT DEM', fontsize = 9, fontweight = 'bold', va = 'center', ha = 'center')
txt.set_clip_on(False)
ax1_3.axis('off')  

txt = ax2_1.text(0.5, 0.5, 'RCP 26', fontsize = 9, fontweight = 'bold', rotation=270, va = 'center')
txt.set_clip_on(False)
ax2_1.axis('off')   
txt = ax3_1.text(0.5, 0.5, 'RCP 45', fontsize = 9, fontweight = 'bold',rotation=270, va = 'center')
txt.set_clip_on(False)
ax3_1.axis('off')   
txt = ax4_1.text(0.5, 0.5, 'RCP 85', fontsize = 9, fontweight = 'bold', rotation=270, va = 'center')
txt.set_clip_on(False)
ax4_1.axis('off')   

barWidth = 0.85
r = [0, 1, 2, 3]
# Data
raw_data = {'greenBars': [round(i[0], 4) for i in coastal_rcp_26], 'goldBars': [round(i[1], 4) for i in coastal_rcp_26], 'pinkBars': [round(i[2], 4) for i in coastal_rcp_26],\
            'redBars': [round(i[3], 4) for i in coastal_rcp_26], 'darkredBars': [round(i[4], 4) for i in coastal_rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax2_2.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in gebco_rcp_26], 'goldBars': [round(i[1], 4) for i in gebco_rcp_26], 'pinkBars': [round(i[2], 4) for i in gebco_rcp_26],\
            'redBars': [round(i[3], 4) for i in gebco_rcp_26], 'darkredBars': [round(i[4], 4) for i in gebco_rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax2_3.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in merit_rcp_26], 'goldBars': [round(i[1], 4) for i in merit_rcp_26], 'pinkBars': [round(i[2], 4) for i in merit_rcp_26],\
            'redBars': [round(i[3], 4) for i in merit_rcp_26], 'darkredBars': [round(i[4], 4) for i in merit_rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax2_4.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_4.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_4.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_4.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_4.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in coastal_rcp_45], 'goldBars': [round(i[1], 4) for i in coastal_rcp_45], 'pinkBars': [round(i[2], 4) for i in coastal_rcp_45],\
            'redBars': [round(i[3], 4) for i in coastal_rcp_45], 'darkredBars': [round(i[4], 4) for i in coastal_rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax3_2.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax3_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax3_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax3_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax3_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in gebco_rcp_45], 'goldBars': [round(i[1], 4) for i in gebco_rcp_45], 'pinkBars': [round(i[2], 4) for i in gebco_rcp_45],\
            'redBars': [round(i[3], 4) for i in gebco_rcp_45], 'darkredBars': [round(i[4], 4) for i in gebco_rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax3_3.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax3_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax3_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax3_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax3_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in merit_rcp_45], 'goldBars': [round(i[1], 4) for i in merit_rcp_45], 'pinkBars': [round(i[2], 4) for i in merit_rcp_45],\
            'redBars': [round(i[3], 4) for i in merit_rcp_45], 'darkredBars': [round(i[4], 4) for i in merit_rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax3_4.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax3_4.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax3_4.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax3_4.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax3_4.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in coastal_rcp_85], 'goldBars': [round(i[1], 4) for i in coastal_rcp_85], 'pinkBars': [round(i[2], 4) for i in coastal_rcp_85],\
            'redBars': [round(i[3], 4) for i in coastal_rcp_85], 'darkredBars': [round(i[4], 4) for i in coastal_rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax4_2.bar(r, greenBars, color = 'forestgreen', width = barWidth, label = '< -5%')
ax4_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth, label = '-5% < AND < -10%')
ax4_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth, label = '-10% < AND < -25%')
ax4_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth, label = '-25% < AND < -50%')
ax4_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth, label = '> -50%')

raw_data = {'greenBars': [round(i[0], 4) for i in gebco_rcp_85], 'goldBars': [round(i[1], 4) for i in gebco_rcp_85], 'pinkBars': [round(i[2], 4) for i in gebco_rcp_85],\
            'redBars': [round(i[3], 4) for i in gebco_rcp_85], 'darkredBars': [round(i[4], 4) for i in gebco_rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax4_3.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax4_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax4_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax4_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax4_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in merit_rcp_85], 'goldBars': [round(i[1], 4) for i in merit_rcp_85], 'pinkBars': [round(i[2], 4) for i in merit_rcp_85],\
            'redBars': [round(i[3], 4) for i in merit_rcp_85], 'darkredBars': [round(i[4], 4) for i in merit_rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax4_4.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax4_4.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax4_4.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax4_4.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax4_4.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

# Custom x axis
names = ('2050','2100','2200','2300')
ax2_2.set_xticks(r)
ax2_2.set_xticklabels([])
ax2_2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_2.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax2_2.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax2_2.set_ylim(1, 0)

ax2_3.set_xticks(r)
ax2_3.set_xticklabels([])
ax2_3.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_3.set_yticklabels([])
ax2_3.set_ylim(1, 0)

ax2_4.set_xticks(r)
ax2_4.set_xticklabels([])
ax2_4.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_4.set_yticklabels([])
ax2_4.set_ylim(1, 0)

# Custom x axis
ax3_2.set_xticks(r)
ax3_2.set_xticklabels([])
ax3_2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax3_2.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax3_2.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax3_2.set_ylim(1, 0)

ax3_3.set_xticks(r)
ax3_3.set_xticklabels([])
ax3_3.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax3_3.set_yticklabels([])
ax3_3.set_ylim(1, 0)

ax3_4.set_xticks(r)
ax3_4.set_xticklabels([])
ax3_4.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax3_4.set_yticklabels([])
ax3_4.set_ylim(1, 0)

# Custom x axis
ax4_2.set_xticks(r)
ax4_2.set_xticklabels(names, fontsize = 7)
ax4_2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax4_2.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax4_2.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax4_2.set_ylim(1, 0)
ax4_2.set_xlabel("Year", fontsize = 7)

ax4_3.set_xticks(r)
ax4_3.set_xticklabels(names, fontsize = 7)
ax4_3.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax4_3.set_yticklabels([])
ax4_3.set_ylim(1, 0)
ax4_3.set_xlabel("Year", fontsize = 7)

ax4_4.set_xticks(r)
ax4_4.set_xticklabels(names, fontsize = 7)
ax4_4.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax4_4.set_yticklabels([])
ax4_4.set_ylim(1, 0)
ax4_4.set_xlabel("Year", fontsize = 7)

h, l = ax4_2.get_legend_handles_labels() 
ax6_1.legend(h, l, fontsize = 7, frameon = False, ncol = 5, loc = 'lower center', title = "IFGW difference compared to 2000 [%]", title_fontsize = 8)
ax6_1.axis('off')

canvas = FigureCanvas(fig)
#canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_urban_diff_2000.png'), dpi=300)
#canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_agriculture_diff_2000.png'), dpi=300)
#canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_nature_diff_2000.png'), dpi=300)
canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_ALL_diff_2000_v2.png'), dpi=300)








        

"""  -------------------------------------------------------------------------
            Create the plots comparing DEM results for each RCP
--------------------------------------------------------------------------  """

import pandas as pd

#   read in the results and land use csv
lu_pop_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\_land_use_population_GDP.csv')
rcp_26_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_26_results.csv')
rcp_45_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_45_results.csv')
rcp_85_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_85_results.csv')

merit_rcp_26, coastal_rcp_26, gebco_rcp_26 = [], [], []
merit_rcp_45, coastal_rcp_45, gebco_rcp_45 = [], [], []
merit_rcp_85, coastal_rcp_85, gebco_rcp_85 = [], [], []

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300', 'IFGW_me_2400', 'IFGW_me_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
            #   open the rcp result list and then append to the merit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_26['IFGW_me_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        








import fiona
import rasterio
from rasterio.mask import mask
#from shapely.geometry import mapping
import numpy as np
import pandas as pd
import shapely.wkt
import shapely.geometry

from shapely.geometry import shape, mapping
from shapely.ops import unary_union
from shapely.ops import cascaded_union
from shapely.geometry import Polygon, Point

#import geopandas as gpd
#shapefile = gpd.read_file(r'g:\Water_Nexus\_A4_paper\_data_output\coast_buffer_10km.shp')


shapefile = fiona.open(r'g:\Water_Nexus\_A4_paper\_data_output\coast_buffer_10km.shp')
#shapefile = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\coast_buffer_10km.csv')
full_dataset = rasterio.open(r'g:\_ORIGINAL_DATA\Land Cover GLC-SHARE\glc_shv10_DOM.Tif')
population_dataset = rasterio.open(r'g:\_ORIGINAL_DATA\_GPW_population_2020\gpw_v4_population_count_rev11_2020_2pt5_min.asc')
gdp_dataset = rasterio.open(r'g:\_ORIGINAL_DATA\GDP_2018\GDP_2015_5arcmin_v2.tiff')


ts_name_lst = [['RCP_26_I_1', 'RCP_26_2050_diff_ifgw'], ['RCP_26_I_2', 'RCP_26_2100_diff_ifgw'], ['RCP_26_I_3', 'RCP_26_2200_diff_ifgw'],\
               ['RCP_26_I_4', 'RCP_26_2300_diff_ifgw'], ['RCP_26_I_5', 'RCP_26_2400_diff_ifgw'], ['RCP_26_I_6', 'RCP_26_2500_diff_ifgw'],\
               ['RCP_45_I_1', 'RCP_45_2050_diff_ifgw'], ['RCP_45_I_2', 'RCP_45_2100_diff_ifgw'], ['RCP_45_I_3', 'RCP_45_2200_diff_ifgw'],\
               ['RCP_45_I_4', 'RCP_45_2300_diff_ifgw'], ['RCP_45_I_5', 'RCP_45_2400_diff_ifgw'], ['RCP_45_I_6', 'RCP_45_2500_diff_ifgw'],\
               ['RCP_85_I_1', 'RCP_85_2050_diff_ifgw'], ['RCP_85_I_2', 'RCP_85_2100_diff_ifgw'], ['RCP_85_I_3', 'RCP_85_2200_diff_ifgw'],\
               ['RCP_85_I_4', 'RCP_85_2300_diff_ifgw'], ['RCP_85_I_5', 'RCP_85_2400_diff_ifgw'], ['RCP_85_I_6', 'RCP_85_2500_diff_ifgw']]


out_csv_cols = ['RCP', 'TIME', 'POP_TOTAL_1', 'POP_TOTAL_2', 'POP_TOTAL_3', 'POP_TOTAL_4', 'POP_TOTAL_5', 'GDP_TOTAL_1',\
                'GDP_TOTAL_2', 'GDP_TOTAL_3', 'GDP_TOTAL_4', 'GDP_TOTAL_5']
out_csv_data = []

for ts_name in ts_name_lst:
    print(ts_name[0])
    lst_1, lst_2, lst_3, lst_4, lst_5 = [], [], [], [], []
    for feature in shapefile:
        #print(feature['properties']['id_srm'])
        if feature['properties'][ts_name[0]] >= -5.:
            lst_1.append(feature['properties']['id_srm'])
        elif feature['properties'][ts_name[0]] < -5. and feature['properties'][ts_name[0]] >= -10.:
           lst_2.append(feature['properties']['id_srm'])      
        elif feature['properties'][ts_name[0]] < -10. and feature['properties'][ts_name[0]] >= -25.:
            lst_3.append(feature['properties']['id_srm'])   
        elif feature['properties'][ts_name[0]] < -25. and feature['properties'][ts_name[0]] >= -50.:
            lst_4.append(feature['properties']['id_srm'])    
        elif feature['properties'][ts_name[0]] >= -100.:
            lst_5.append(feature['properties']['id_srm'])   

    #   create the merged geometry
    polys_1 = []
    for srm in lst_1:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_1.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_1.append(Polygon(feature['geometry']['coordinates'][0]))

    #   create the merged geometry
    polys_2 = []
    for srm in lst_2:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_2.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_2.append(Polygon(feature['geometry']['coordinates'][0]))    

    #   create the merged geometry
    polys_3 = []
    for srm in lst_3:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_3.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_3.append(Polygon(feature['geometry']['coordinates'][0]))
                     
    #   create the merged geometry
    polys_4 = []
    for srm in lst_4:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_4.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_4.append(Polygon(feature['geometry']['coordinates'][0]))
                     
    #   create the merged geometry
    polys_5 = []
    for srm in lst_5:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_5.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_5.append(Polygon(feature['geometry']['coordinates'][0]))

    #   merge the geometry 
    merged_1 = unary_union(polys_1)
    out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_1)], crop=True)
    out_image[out_image < 0] = 0
    total_population_1 = np.sum(out_image)
    out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_1)], crop=True)
    out_image[out_image < 0] = 0
    total_gdp_1 = np.nansum(out_image)
    
    merged_2 = unary_union(polys_2)
    out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_2)], crop=True)
    out_image[out_image < 0] = 0
    total_population_2 = np.sum(out_image)
    out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_2)], crop=True)
    out_image[out_image < 0] = 0
    total_gdp_2 = np.nansum(out_image)
    
    merged_3 = unary_union(polys_3)
    out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_3)], crop=True)
    out_image[out_image < 0] = 0
    total_population_3 = np.sum(out_image)
    out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_3)], crop=True)
    out_image[out_image < 0] = 0
    total_gdp_3 = np.nansum(out_image)
    
    if len(polys_4) > 0:
        merged_4 = unary_union(polys_4)
        out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_4)], crop=True)
        out_image[out_image < 0] = 0
        total_population_4 = np.sum(out_image)
        out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_4)], crop=True)
        out_image[out_image < 0] = 0
        total_gdp_4 = np.nansum(out_image)
    else:
        total_population_4 = 0
        total_gdp_4 = 0

    if len(polys_5) > 0:    
        merged_5 = unary_union(polys_5)
        out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_5)], crop=True)
        out_image[out_image < 0] = 0
        total_population_5 = np.sum(out_image)
        out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_5)], crop=True)
        out_image[out_image < 0] = 0
        total_gdp_5 = np.nansum(out_image)
    else:
        total_population_5 = 0
        total_gdp_5 = 0

    rcp = ts_name[0].split('_I')[0]
    time = int(ts_name[1].split('_')[2])

    out_csv_data.append([rcp, time, total_population_1, total_population_2, total_population_3, total_population_4, total_population_5,\
                         total_gdp_1, total_gdp_2, total_gdp_3, total_gdp_4, total_gdp_5])


    
my_df = pd.DataFrame(out_csv_data)
my_df.to_csv(r'g:\Water_Nexus\_A4_paper\_data_output\POPULATION_GDP_affected.csv', index=False, header=out_csv_cols)





#   ------------------------------------------------------------------------------------
#   Do the same but for the positions of SWW - Salt water wedge  
#   ------------------------------------------------------------------------------------


#   Salinity treshold - 17.5 ppm

out_cols = ['fid', 'coscat', 'srm', 'cs_srm_id', 'inl_ext_km',\
            'SWW_2000_mu', 'SWW_2000_std', 'SWW_2050_mu', 'SWW_2050_std',\
            'SWW_2100_mu', 'SWW_2100_std', 'SWW_2200_mu', 'SWW_2200_std',\
            'SWW_2300_mu', 'SWW_2300_std', 'SWW_2400_mu', 'SWW_2400_std',\
            'SWW_2500_mu', 'SWW_2500_std']

for rcp in ['26', '45', '85']:

    if not os.path.exists(main_dir + rcp + '_SWW_mu_std.csv'):
        df = pd.DataFrame(columns = out_cols)
        df.set_index('fid')
        df.to_csv(main_dir + rcp + '_SWW_mu_std.csv' , encoding = 'utf-8', index = False)
    
    rcp_in_dir = main_dir + rcp + '_results.csv'
    rcp_in = pd.read_csv(rcp_in_dir)
    fid = 1
    
    for index, row in rcp_in.iterrows():
        
        cs_id = row['coscat']
        srm_id = row['srm']
        cs_srm_id = row['cs_srm_id']
        inl_ext = row['inl_ext_km']
        
        SWW_2000_mu = round(np.mean(np.array([x for x in [row['SWW_me_2000'], row['SWW_co_2000'], row['SWW_ge_2000']] if x != -1])), 1)
        SWW_2000_std = round(np.std(np.array([x for x in [row['SWW_me_2000'], row['SWW_co_2000'], row['SWW_ge_2000']] if x != -1])), 1)    
        SWW_2050_mu = round(np.mean(np.array([x for x in [row['SWW_me_2050'], row['SWW_co_2050'], row['SWW_ge_2050']] if x != -1])), 1)
        SWW_2050_std = round(np.std(np.array([x for x in [row['SWW_me_2050'], row['SWW_co_2050'], row['SWW_ge_2050']] if x != -1])), 1)    
        SWW_2100_mu = round(np.mean(np.array([x for x in [row['SWW_me_2100'], row['SWW_co_2100'], row['SWW_ge_2100']] if x != -1])), 1)
        SWW_2100_std = round(np.std(np.array([x for x in [row['SWW_me_2100'], row['SWW_co_2100'], row['SWW_ge_2100']] if x != -1])), 1)    
        SWW_2200_mu = round(np.mean(np.array([x for x in [row['SWW_me_2200'], row['SWW_co_2200'], row['SWW_ge_2200']] if x != -1])), 1)
        SWW_2200_std = round(np.std(np.array([x for x in [row['SWW_me_2200'], row['SWW_co_2200'], row['SWW_ge_2200']] if x != -1])), 1)    
        SWW_2300_mu = round(np.mean(np.array([x for x in [row['SWW_me_2300'], row['SWW_co_2300'], row['SWW_ge_2300']] if x != -1])), 1)
        SWW_2300_std = round(np.std(np.array([x for x in [row['SWW_me_2300'], row['SWW_co_2300'], row['SWW_ge_2300']] if x != -1])), 1)    
        SWW_2400_mu = round(np.mean(np.array([x for x in [row['SWW_me_2400'], row['SWW_co_2400'], row['SWW_ge_2400']] if x != -1])), 1)
        SWW_2400_std = round(np.std(np.array([x for x in [row['SWW_me_2400'], row['SWW_co_2400'], row['SWW_ge_2400']] if x != -1])), 1)        
        SWW_2500_mu = round(np.mean(np.array([x for x in [row['SWW_me_2500'], row['SWW_co_2500'], row['SWW_ge_2500']] if x != -1])), 1)
        SWW_2500_std = round(np.std(np.array([x for x in [row['SWW_me_2500'], row['SWW_co_2500'], row['SWW_ge_2500']] if x != -1])), 1)        
    
        f = open(main_dir + rcp + '_SWW_mu_std.csv', 'a')
        f.write(str(fid) + ',' + str(cs_id) + ',' + str(srm_id) + ',' + cs_srm_id + ',' + str(inl_ext) + ',' + str(SWW_2000_mu) + ',' + str(SWW_2000_std)\
                 + ',' + str(SWW_2050_mu) + ',' + str(SWW_2050_std) + ',' + str(SWW_2100_mu) + ',' + str(SWW_2100_std)\
                 + ',' + str(SWW_2200_mu) + ',' + str(SWW_2200_std) + ',' + str(SWW_2300_mu) + ',' + str(SWW_2300_std) + ',' + str(SWW_2400_mu) + ',' + str(SWW_2400_std)\
                 + ',' + str(SWW_2500_mu) + ',' + str(SWW_2500_std))
        f.write("\n")
        f.close()                  
    
        fid += 1
    
#   ------------------------------------------------------------------------------------
#   Also create csv files showing difference between future predictions and current state
#   ------------------------------------------------------------------------------------
        
out_cols = ['fid', 'coscat', 'srm', 'cs_srm_id', 'inl_ext_km',\
            'SWW_2000_DIFF_mu', 'SWW_2000_DIFF_std', 'SWW_2050_DIFF_mu', 'SWW_2050_DIFF_std',\
            'SWW_2100_DIFF_mu', 'SWW_2100_DIFF_std', 'SWW_2200_DIFF_mu', 'SWW_2200_DIFF_std',\
            'SWW_2300_DIFF_mu', 'SWW_2300_DIFF_std', 'SWW_2400_DIFF_mu', 'SWW_2400_DIFF_std',\
            'SWW_2500_DIFF_mu', 'SWW_2500_DIFF_std']    
    
for rcp in ['26', '45', '85']:

    if not os.path.exists(main_dir + rcp + '_SWW_mu_std_DIFF_2000.csv'):
        df = pd.DataFrame(columns = out_cols)
        df.set_index('fid')
        df.to_csv(main_dir + rcp + '_SWW_mu_std_DIFF_2000.csv' , encoding = 'utf-8', index = False)
    
    rcp_in_dir = main_dir + rcp + '_results.csv'
    rcp_in = pd.read_csv(rcp_in_dir)
    fid = 1
    
    for index, row in rcp_in.iterrows():
        
        cs_id = row['coscat']
        srm_id = row['srm']
        cs_srm_id = row['cs_srm_id']
        inl_ext = row['inl_ext_km']
        
        SWW_2000_mu = round(np.mean(np.array([x for x in [row['SWW_me_2000'], row['SWW_co_2000'], row['SWW_ge_2000']] if x != -1])), 1)
        SWW_2000_std = round(np.std(np.array([x for x in [row['SWW_me_2000'], row['SWW_co_2000'], row['SWW_ge_2000']] if x != -1])), 1)
    
        SWW_2050_mu = round(np.mean(np.array([x for x in [row['SWW_me_2050'], row['SWW_co_2050'], row['SWW_ge_2050']] if x != -1])) - SWW_2000_mu, 1)
        SWW_2050_std = round(np.std(np.array([x for x in [row['SWW_me_2050'], row['SWW_co_2050'], row['SWW_ge_2050']] if x != -1])), 1)    
        SWW_2100_mu = round(np.mean(np.array([x for x in [row['SWW_me_2100'], row['SWW_co_2100'], row['SWW_ge_2100']] if x != -1])) - SWW_2000_mu, 1)
        SWW_2100_std = round(np.std(np.array([x for x in [row['SWW_me_2100'], row['SWW_co_2100'], row['SWW_ge_2100']] if x != -1])), 1)    
        SWW_2200_mu = round(np.mean(np.array([x for x in [row['SWW_me_2200'], row['SWW_co_2200'], row['SWW_ge_2200']] if x != -1])) - SWW_2000_mu, 1)
        SWW_2200_std = round(np.std(np.array([x for x in [row['SWW_me_2200'], row['SWW_co_2200'], row['SWW_ge_2200']] if x != -1])), 1)    
        SWW_2300_mu = round(np.mean(np.array([x for x in [row['SWW_me_2300'], row['SWW_co_2300'], row['SWW_ge_2300']] if x != -1])) - SWW_2000_mu, 1)
        SWW_2300_std = round(np.std(np.array([x for x in [row['SWW_me_2300'], row['SWW_co_2300'], row['SWW_ge_2300']] if x != -1])), 1)    
        SWW_2400_mu = round(np.mean(np.array([x for x in [row['SWW_me_2400'], row['SWW_co_2400'], row['SWW_ge_2400']] if x != -1])) - SWW_2000_mu, 1)
        SWW_2400_std = round(np.std(np.array([x for x in [row['SWW_me_2400'], row['SWW_co_2400'], row['SWW_ge_2400']] if x != -1])), 1)        
        SWW_2500_mu = round(np.mean(np.array([x for x in [row['SWW_me_2500'], row['SWW_co_2500'], row['SWW_ge_2500']] if x != -1])) - SWW_2000_mu, 1)
        SWW_2500_std = round(np.std(np.array([x for x in [row['SWW_me_2500'], row['SWW_co_2500'], row['SWW_ge_2500']] if x != -1])), 1)        
    
        f = open(main_dir + rcp + '_SWW_mu_std_DIFF_2000.csv', 'a')
        f.write(str(fid) + ',' + str(cs_id) + ',' + str(srm_id) + ',' + cs_srm_id + ',' + str(inl_ext) + ',' + str(0) + ',' + str(0)\
                 + ',' + str(SWW_2050_mu) + ',' + str(SWW_2050_std) + ',' + str(SWW_2100_mu) + ',' + str(SWW_2100_std)\
                 + ',' + str(SWW_2200_mu) + ',' + str(SWW_2200_std) + ',' + str(SWW_2300_mu) + ',' + str(SWW_2300_std) + ',' + str(SWW_2400_mu) + ',' + str(SWW_2400_std)\
                 + ',' + str(SWW_2500_mu) + ',' + str(SWW_2500_std))
        f.write("\n")
        f.close()                  
    
        fid += 1


"""    ---------------------------------------------------------------    """

import matplotlib.pyplot as plt
import math
import os
import matplotlib

import os
os.environ["PROJ_LIB"] = r'c:\Miniconda3\pkgs\proj4-5.2.0-h6538335_1006\Library\share'

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
from matplotlib.patches import PathPatch
from matplotlib.colors import Normalize
from matplotlib import colors
from matplotlib import colorbar
import matplotlib.cm as cmx
import matplotlib.patches as mpatches
import rasterio

#   ------------------------------------------------------------------------------------
#   Create final plots (worldmaps) for each RCP and year  
#   ------------------------------------------------------------------------------------

srm_ifgw_diff_dir = os.path.join(data_out_dir, '_SRM_SWW_diff_2000.shp') #r'g:\Water_Nexus\_A4_paper\_data_output\_SRM_IFGW_diff_2000.shp'
basemap_raster_bw = os.path.join(data_out_dir, 'background_bw.tif') #r'g:\Water_Nexus\_A4_paper\_data_output\background_bw.tif'
basemap_raster_col = os.path.join(data_out_dir, 'background_col.tif')#r'g:\_ORIGINAL_DATA\natural_earth\NE2_LR_LC_SR\NE2_LR_LC_SR.tif'


#from mpl_toolkits.basemap import Basemap
#import matplotlib.pyplot as plt
import numpy as np
#from itertools import chain
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import pandas as pd

ts_name_lst = [['RCP_26_SWW', 'RCP_26_2050_diff_sww'], ['RCP_26_S_1', 'RCP_26_2100_diff_sww'], ['RCP_26_S_2', 'RCP_26_2200_diff_sww'],\
               ['RCP_26_S_3', 'RCP_26_2300_diff_sww'], ['RCP_26_S_4', 'RCP_26_2400_diff_sww'], ['RCP_26_S_5', 'RCP_26_2500_diff_sww'],\
               ['RCP_45_SWW', 'RCP_45_2050_diff_sww'], ['RCP_45_S_1', 'RCP_45_2100_diff_sww'], ['RCP_45_S_2', 'RCP_45_2200_diff_sww'],\
               ['RCP_45_S_3', 'RCP_45_2300_diff_sww'], ['RCP_45_S_4', 'RCP_45_2400_diff_sww'], ['RCP_45_S_5', 'RCP_45_2500_diff_sww'],\
               ['RCP_85_SWW', 'RCP_85_2050_diff_sww'], ['RCP_85_S_1', 'RCP_85_2100_diff_sww'], ['RCP_85_S_2', 'RCP_85_2200_diff_sww'],\
               ['RCP_85_S_3', 'RCP_85_2300_diff_sww'], ['RCP_85_S_4', 'RCP_85_2400_diff_sww'], ['RCP_85_S_5', 'RCP_85_2500_diff_sww']]

ts_pct_lst = []

for i in range(len(ts_name_lst)):

    #   plot the basemap and the results
    fig = plt.figure(figsize=(8 , 4), edgecolor='w')
    ax = fig.add_subplot(111)
    
    """
    # A good LCC projection for USA plots
    m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 75,\
                llcrnrlon = -140, urcrnrlon = 175, lat_ts = 20, resolution = 'i')
    
    m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
    m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
    """
    
    # A good LCC projection for USA plots
    m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = 0, urcrnrlat = 60,\
                llcrnrlon = -25, urcrnrlon = 95, lat_ts = 20, resolution = 'i')
    m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
    m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)

    m.drawcountries(linewidth = 0.1, linestyle = 'solid', color = 'grey') 
    m.drawcoastlines(linewidth = 0.1, linestyle = 'solid', color = 'grey')
    m.shadedrelief(scale = 0.2, alpha = 0.25)
    #m.etopo(scale = 0.2, alpha = 0.5)
    #   read in the shapefile
    m.readshapefile(os.path.join(data_out_dir, '_SRM_SWW_diff_2000'), 'SRM_SWW_diff')
    
    print(ts_name_lst[i][1])
    ts_cnt = [0, 0, 0, 0, 0]
    
    for shapedict,shape in zip(m.SRM_IFGW_diff_info, m.SRM_IFGW_diff):
        try:
            sww = float(shapedict[ts_name_lst[i][0]])
            xx, yy = zip(*shape)
            # show part of track where storm > Cat 4 as thick red.
            if ifgw is not None:
                if sww <= 1.0:
                    m.plot(xx, yy, linewidth = 1.25, color = 'forestgreen')
                    ts_cnt[0] += 1
                elif sww > 1.0 and sww <= 2.5:
                    m.plot(xx, yy, linewidth = 1.25, color = 'gold')
                    ts_cnt[1] += 1
                elif sww > 2.5 and sww <= 5.0:
                    m.plot(xx, yy, linewidth = 1.25, color = 'hotpink')
                    ts_cnt[2] += 1
                elif sww > 5.0 and sww <= 10.0:
                    m.plot(xx, yy, linewidth = 1.25, color = 'red')
                    ts_cnt[3] += 1
                elif sww > 10.0:
                    m.plot(xx, yy, linewidth = 1.25, color = 'darkred')
                    ts_cnt[4] += 1
            else:
                m.plot(xx, yy, linewidth = .1, color = 'k')
        except ValueError:
            pass
            
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    ts_pct_lst.append(ts_pct)
            
    plt.savefig(os.path.join(data_out_dir, ts_name_lst[i][1] + '_zoom_300dpi.png'), dpi = 300, facecolor = 'w', edgecolor = 'w',
        orientation = 'portrait', papertype = None, format = None,
        transparent = False, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()




"""  -------------------------------------------------------------------------
            FIGURE WITH 3 WORLD MAPS FOR THE MAIN ARTICLE
--------------------------------------------------------------------------  """

fig = plt.figure(figsize=(8 , 8), edgecolor='w') 
                                                                                                                                                                                             
ax1_1 = plt.subplot2grid((4, 1), (0, 0), fig = fig)             
ax1_2 = plt.subplot2grid((4, 1), (1, 0), fig = fig)             
ax1_3 = plt.subplot2grid((4, 1), (2, 0), fig = fig)             
ax1_4 = plt.subplot2grid((4, 1), (3, 0), fig = fig)     
        
ax1_1.set_position([0.025, 0.7, 0.95, 0.295]) # [left, bottom, width, height]                     
ax1_2.set_position([0.025, 0.4, 0.95, 0.295])
ax1_3.set_position([0.025, 0.1, 0.95, 0.295])  
ax1_4.set_position([0.25, 0.025, 0.5, 0.05])  

m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 75,\
            llcrnrlon = -140, urcrnrlon = 175, lat_ts = 20, resolution = 'i', ax = ax1_1)
m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawcountries(linewidth = 0.1, linestyle = 'solid', color = 'grey') 
m.drawcoastlines(linewidth = 0.1, linestyle = 'solid', color = 'grey')
m.shadedrelief(scale = 0.2, alpha = 0.25)
#m.etopo(scale = 0.2, alpha = 0.5)
#   read in the shapefile
m.readshapefile(os.path.join(data_out_dir, '_SRM_SWW_diff_2000'), 'SRM_SWW_diff')
for shapedict,shape in zip(m.SRM_SWW_diff_info, m.SRM_SWW_diff):
    # show part of track where storm > Cat 4 as thick red.
    try:
        sww = float(shapedict[ts_name_lst[1][0]])
        xx, yy = zip(*shape)
        # show part of track where storm > Cat 4 as thick red.
        if sww is not None:
            if sww <= 1.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'forestgreen')
                ts_cnt[0] += 1
            elif sww > 1.0 and sww <= 2.5:
                m.plot(xx, yy, linewidth = 1.25, color = 'gold')
                ts_cnt[1] += 1
            elif sww > 2.5 and sww <= 5.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'hotpink')
                ts_cnt[2] += 1
            elif sww > 5.0 and sww <= 10.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'red')
                ts_cnt[3] += 1
            elif sww > 10.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'darkred')
                ts_cnt[4] += 1
        else:
            m.plot(xx, yy, linewidth = .1, color = 'k')
    except ValueError:
        pass

# A good LCC projection for USA plots
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 75,\
            llcrnrlon = -140, urcrnrlon = 175, lat_ts = 20, resolution = 'i', ax = ax1_2)
m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawcountries(linewidth = 0.1, linestyle = 'solid', color = 'grey') 
m.drawcoastlines(linewidth = 0.1, linestyle = 'solid', color = 'grey')
m.shadedrelief(scale = 0.2, alpha = 0.25)
#m.etopo(scale = 0.2, alpha = 0.5)
#   read in the shapefile
m.readshapefile(os.path.join(data_out_dir, '_SRM_SWW_diff_2000'), 'SRM_SWW_diff')
for shapedict,shape in zip(m.SRM_SWW_diff_info, m.SRM_SWW_diff):
    try:
        sww = float(shapedict[ts_name_lst[7][0]])
        xx, yy = zip(*shape)
        # show part of track where storm > Cat 4 as thick red.
        if sww is not None:
            if sww <= 1.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'forestgreen', label = '< 1')
                ts_cnt[0] += 1
            elif sww > 1.0 and sww <= 2.5:
                m.plot(xx, yy, linewidth = 1.25, color = 'gold', label = '1 < AND < 2.5')
                ts_cnt[1] += 1
            elif sww > 2.5 and sww <= 5.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'hotpink', label = '2.5 < AND < 5')
                ts_cnt[2] += 1
            elif sww > 5.0 and sww <= 10.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'red', label = '5 < AND < 10')
                ts_cnt[3] += 1
            elif sww > 10.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'darkred', label = '> 10')
                ts_cnt[4] += 1
        else:
            m.plot(xx, yy, linewidth = .1, color = 'k')
    except ValueError:
        pass

# A good LCC projection for USA plots
m = Basemap(projection = 'cyl', ellps = 'wgs84', llcrnrlat = -60, urcrnrlat = 75,\
            llcrnrlon = -140, urcrnrlon = 175, lat_ts = 20, resolution = 'i', ax = ax1_3)
m.drawparallels(np.arange(-90.,91.,30.), labels = [True, False, False, False], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawmeridians(np.arange(-180.,181.,30.), labels= [False, False, False, True], linewidth = 0.25, color = 'grey', dashes = [2,2], fontsize = 7)
m.drawcountries(linewidth = 0.1, linestyle = 'solid', color = 'grey') 
m.drawcoastlines(linewidth = 0.1, linestyle = 'solid', color = 'grey')
m.shadedrelief(scale = 0.2, alpha = 0.25)
#m.etopo(scale = 0.2, alpha = 0.5)
#   read in the shapefile
m.readshapefile(os.path.join(data_out_dir, '_SRM_SWW_diff_2000'), 'SRM_SWW_diff')
for shapedict,shape in zip(m.SRM_SWW_diff_info, m.SRM_SWW_diff):
    try:
        sww = float(shapedict[ts_name_lst[13][0]])
        xx, yy = zip(*shape)
        # show part of track where storm > Cat 4 as thick red.
        if sww is not None:
            if sww <= 1.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'forestgreen')
                ts_cnt[0] += 1
            elif sww > 1.0 and sww <= 2.5:
                m.plot(xx, yy, linewidth = 1.25, color = 'gold')
                ts_cnt[1] += 1
            elif sww > 2.5 and sww <= 5.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'hotpink')
                ts_cnt[2] += 1
            elif sww > 5.0 and sww <= 10.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'red')
                ts_cnt[3] += 1
            elif sww > 10.0:
                m.plot(xx, yy, linewidth = 1.25, color = 'darkred')
                ts_cnt[4] += 1
        else:
            m.plot(xx, yy, linewidth = .1, color = 'k')
    except ValueError:
        pass

h, l = ax1_2.get_legend_handles_labels() 
handle_list, label_list = [], []
for h, l in zip(h, l):
    if l not in label_list:
        handle_list.append(h)
        label_list.append(l)

ax1_4.legend(handle_list, label_list, fontsize = 7, frameon = False, ncol = 5, loc = 'lower center', title = "SWW difference compared to 2000 [km inland]", title_fontsize = 8)
ax1_4.axis('off')

canvas = FigureCanvas(fig)
canvas.print_figure(os.path.join(data_out_dir, '_SWW_A4_figure.png'), dpi=300)


"""  -------------------------------------------------------------------------
    HISTOGRAM WITH AVERAGE IFGW DIFFERENCE TO 2000, ALL DEMS COMBINED
--------------------------------------------------------------------------  """

fig = plt.figure(figsize=(8 , 3), edgecolor='w') 
                                                                                                                                                                                             
ax1_1 = plt.subplot2grid((3, 3), (0, 0), fig = fig)             
ax1_2 = plt.subplot2grid((3, 3), (0, 1), fig = fig)             
ax1_3 = plt.subplot2grid((3, 3), (0, 2), fig = fig)             
ax2_1 = plt.subplot2grid((3, 3), (1, 0), fig = fig)             
ax2_2 = plt.subplot2grid((3, 3), (1, 1), fig = fig)             
ax2_3 = plt.subplot2grid((3, 3), (1, 2), fig = fig)     
ax3 = plt.subplot2grid((3, 3), (2, 0), colspan = 3, fig = fig)     

ax1_1.set_position([0.1, 0.925, 0.266, 0.05])                      
ax1_2.set_position([0.366, 0.925, 0.266, 0.05])
ax1_3.set_position([0.632, 0.925, 0.266, 0.05])  
ax2_1.set_position([0.1, 0.25, 0.266, 0.65])  # [left, bottom, width, height]
ax2_2.set_position([0.366, 0.25, 0.266, 0.65])
ax2_3.set_position([0.632, 0.25, 0.266, 0.65])
ax3.set_position([0.1, 0., 0.8, 0.055])

txt = ax1_1.text(0.4, 0.5, 'RCP 26', fontsize = 7, fontweight = 'bold', va = 'center')
txt.set_clip_on(False)
ax1_1.axis('off')   
txt = ax1_2.text(0.4, 0.5, 'RCP 45', fontsize = 7, fontweight = 'bold', va = 'center')
txt.set_clip_on(False)
ax1_2.axis('off')   
txt = ax1_3.text(0.4, 0.5, 'RCP 85', fontsize = 7, fontweight = 'bold', va = 'center')
txt.set_clip_on(False)
ax1_3.axis('off')   

rcp_26 = ts_pct_lst[:6]
# Data
r = [0, 1, 2, 3, 4, 5]
raw_data = {'greenBars': [round(i[0], 4) for i in rcp_26], 'goldBars': [round(i[1], 4) for i in rcp_26], 'pinkBars': [round(i[2], 4) for i in rcp_26],\
            'redBars': [round(i[3], 4) for i in rcp_26], 'darkredBars': [round(i[4], 4) for i in rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()

# plot
barWidth = 0.85
names = ('2050','2100','2200','2300','2400', '2500')
# Create green Bars
ax2_1.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_1.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_1.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_1.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_1.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)
# Custom x axis
ax2_1.set_xticks(r)
ax2_1.set_xticklabels(names, fontsize = 7)
ax2_1.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_1.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax2_1.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax2_1.set_xlabel("Year", fontsize = 7)
ax2_1.set_ylim(1, 0)

rcp_45 = ts_pct_lst[6:12]
# Data
r = [0, 1, 2, 3, 4, 5]
raw_data = {'greenBars': [round(i[0], 4) for i in rcp_45], 'goldBars': [round(i[1], 4) for i in rcp_45], 'pinkBars': [round(i[2], 4) for i in rcp_45],\
            'redBars': [round(i[3], 4) for i in rcp_45], 'darkredBars': [round(i[4], 4) for i in rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# plot
barWidth = 0.85
# Create green Bars
ax2_2.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)
# Custom x axis
ax2_2.set_xticks(r)
ax2_2.set_xticklabels(names, fontsize = 7)
ax2_2.get_yaxis().set_ticks([])
ax2_2.set_xlabel("Year", fontsize = 7)
ax2_2.set_ylim(1, 0)

rcp_85 = ts_pct_lst[12:]
# Data
r = [0, 1, 2, 3, 4, 5]
raw_data = {'greenBars': [round(i[0], 4) for i in rcp_85], 'goldBars': [round(i[1], 4) for i in rcp_85], 'pinkBars': [round(i[2], 4) for i in rcp_85],\
            'redBars': [round(i[3], 4) for i in rcp_85], 'darkredBars': [round(i[4], 4) for i in rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# plot
barWidth = 0.85
# Create green Bars
ax2_3.bar(r, greenBars, color = 'forestgreen', width = barWidth, label = '< -5%')
ax2_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth, label = '-5% < AND < -10%')
ax2_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth, label = '-10% < AND < -25%')
ax2_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth, label = '-25% < AND < -50%')
ax2_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth, label = '> -50%')
# Custom x axis
ax2_3.set_xticks(r)
ax2_3.set_xticklabels(names, fontsize = 7)
ax2_3.get_yaxis().set_ticks([])
ax2_3.set_ylim(1, 0)
ax2_3.set_xlabel("Year", fontsize = 7)
h, l = ax2_3.get_legend_handles_labels() 
ax3.legend(h, l, fontsize = 7, frameon = False, ncol = 5, loc = 'lower center', title = "IFGW difference compared to 2000 [%]", title_fontsize = 8)
ax3.axis('off')

canvas = FigureCanvas(fig)
canvas.print_figure(os.path.join(data_out_dir, 'hist_IFGW_diff_2000.png'), dpi=300)



"""  -------------------------------------------------------------------------
            Create the plots comparing DEM results for each RCP
--------------------------------------------------------------------------  """

import pandas as pd

#   read in the results and land use csv
lu_pop_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\_land_use_population_GDP.csv')
rcp_26_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_26_results.csv')
rcp_45_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_45_results.csv')
rcp_85_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_85_results.csv')

merit_rcp_26, coastal_rcp_26, gebco_rcp_26 = [], [], []
merit_rcp_45, coastal_rcp_45, gebco_rcp_45 = [], [], []
merit_rcp_85, coastal_rcp_85, gebco_rcp_85 = [], [], []

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300', 'IFGW_me_2400', 'IFGW_me_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the merit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_26['IFGW_me_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    merit_rcp_26.append(ts_pct)        
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_co_2050', 'IFGW_co_2100', 'IFGW_co_2200', 'IFGW_co_2300', 'IFGW_co_2400', 'IFGW_co_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the corit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_co = row_rcp_26['IFGW_co_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_co
            if val_2000_co is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    coastal_rcp_26.append(ts_pct)         

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_ge_2050', 'IFGW_ge_2100', 'IFGW_ge_2200', 'IFGW_ge_2300', 'IFGW_ge_2400', 'IFGW_ge_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the gerit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_ge = row_rcp_26['IFGW_ge_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_ge
            if val_2000_ge is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    gebco_rcp_26.append(ts_pct)      
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300', 'IFGW_me_2400', 'IFGW_me_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the merit list
            row_rcp_45 = rcp_45_csv.loc[rcp_45_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_45['IFGW_me_2000'].values[0]
            ifgw = row_rcp_45[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    merit_rcp_45.append(ts_pct)        
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_co_2050', 'IFGW_co_2100', 'IFGW_co_2200', 'IFGW_co_2300', 'IFGW_co_2400', 'IFGW_co_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the corit list
            row_rcp_45 = rcp_45_csv.loc[rcp_45_csv['cs_srm_id'] == row['id_srm']]
            val_2000_co = row_rcp_45['IFGW_co_2000'].values[0]
            ifgw = row_rcp_45[ts].values[0] - val_2000_co
            if val_2000_co is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    coastal_rcp_45.append(ts_pct)         

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_ge_2050', 'IFGW_ge_2100', 'IFGW_ge_2200', 'IFGW_ge_2300', 'IFGW_ge_2400', 'IFGW_ge_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the gerit list
            row_rcp_45 = rcp_45_csv.loc[rcp_45_csv['cs_srm_id'] == row['id_srm']]
            val_2000_ge = row_rcp_45['IFGW_ge_2000'].values[0]
            ifgw = row_rcp_45[ts].values[0] - val_2000_ge
            if val_2000_ge is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    gebco_rcp_45.append(ts_pct)              

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300', 'IFGW_me_2400', 'IFGW_me_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the merit list
            row_rcp_85 = rcp_85_csv.loc[rcp_85_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_85['IFGW_me_2000'].values[0]
            ifgw = row_rcp_85[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    merit_rcp_85.append(ts_pct)        
        
#   loop through each row of the lu_pop_csv
for ts in ['IFGW_co_2050', 'IFGW_co_2100', 'IFGW_co_2200', 'IFGW_co_2300', 'IFGW_co_2400', 'IFGW_co_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the corit list
            row_rcp_85 = rcp_85_csv.loc[rcp_85_csv['cs_srm_id'] == row['id_srm']]
            val_2000_co = row_rcp_85['IFGW_co_2000'].values[0]
            ifgw = row_rcp_85[ts].values[0] - val_2000_co
            if val_2000_co is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    coastal_rcp_85.append(ts_pct)         

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_ge_2050', 'IFGW_ge_2100', 'IFGW_ge_2200', 'IFGW_ge_2300', 'IFGW_ge_2400', 'IFGW_ge_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
        #if row['agriculture'] > 10.:
        #if row['nature'] > 90.:
        #if row['nature'] > 0.:
            #   open the rcp result list and then append to the gerit list
            row_rcp_85 = rcp_85_csv.loc[rcp_85_csv['cs_srm_id'] == row['id_srm']]
            val_2000_ge = row_rcp_85['IFGW_ge_2000'].values[0]
            ifgw = row_rcp_85[ts].values[0] - val_2000_ge
            if val_2000_ge is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        
        
    ts_pct = [i/sum(ts_cnt) for i in ts_cnt]
    gebco_rcp_85.append(ts_pct)      


fig = plt.figure(figsize=(8 , 8), edgecolor='w') 
#   column titles
ax1_1 = plt.subplot2grid((6, 5), (0, 1), fig = fig)             
ax1_2 = plt.subplot2grid((6, 5), (0, 2), fig = fig)             
ax1_3 = plt.subplot2grid((6, 5), (0, 3), fig = fig)
#   first row of plots
ax2_1 = plt.subplot2grid((6, 5), (1, 4), fig = fig)             
ax2_2 = plt.subplot2grid((6, 5), (1, 1), fig = fig)             
ax2_3 = plt.subplot2grid((6, 5), (1, 2), fig = fig)     
ax2_4 = plt.subplot2grid((6, 5), (1, 3), fig = fig) 
#   second row of plots
ax3_1 = plt.subplot2grid((6, 5), (2, 4), fig = fig)             
ax3_2 = plt.subplot2grid((6, 5), (2, 1), fig = fig)             
ax3_3 = plt.subplot2grid((6, 5), (2, 2), fig = fig)     
ax3_4 = plt.subplot2grid((6, 5), (2, 3), fig = fig) 
#   third row of plots
ax4_1 = plt.subplot2grid((6, 5), (3, 4), fig = fig)             
ax4_2 = plt.subplot2grid((6, 5), (3, 1), fig = fig)             
ax4_3 = plt.subplot2grid((6, 5), (3, 2), fig = fig)     
ax4_4 = plt.subplot2grid((6, 5), (3, 3), fig = fig) 

#   x axis title
#ax5_1 = plt.subplot2grid((6, 5), (4, 1), colspan = 3, fig = fig) 
#   legend
ax6_1 = plt.subplot2grid((6, 5), (5, 1), colspan = 3, fig = fig)     

#   define the dimensions of each plot 
ax1_1.set_position([0.075, 0.95, 0.275, 0.05])  # [left, bottom, width, height]                     
ax1_2.set_position([0.358, 0.95, 0.275, 0.05])
ax1_3.set_position([0.641, 0.95, 0.275, 0.05])  

ax2_1.set_position([0.925, 0.666, 0.05, 0.275])  # [left, bottom, width, height]                     
ax2_2.set_position([0.075, 0.666, 0.275, 0.275])
ax2_3.set_position([0.358, 0.666, 0.275, 0.275])  
ax2_4.set_position([0.641, 0.666, 0.275, 0.275])

ax3_1.set_position([0.925, 0.383, 0.05, 0.275])  # [left, bottom, width, height]                     
ax3_2.set_position([0.075, 0.383, 0.275, 0.275])
ax3_3.set_position([0.358, 0.383, 0.275, 0.275])  
ax3_4.set_position([0.641, 0.383, 0.275, 0.275])

ax4_1.set_position([0.925, 0.1, 0.05, 0.275])  # [left, bottom, width, height]                     
ax4_2.set_position([0.075, 0.1, 0.275, 0.275])
ax4_3.set_position([0.358, 0.1, 0.275, 0.275])  
ax4_4.set_position([0.641, 0.1, 0.275, 0.275])

#ax5_1.set_position([0.075, 0.055, 0.875, 0.025])  
ax6_1.set_position([0.075, 0.01, 0.875, 0.05])

#   define all the text plots
txt = ax1_1.text(0.5, 0.5, 'Coastal DEM', fontsize = 9, fontweight = 'bold', va = 'center', ha = 'center')
txt.set_clip_on(False)
ax1_1.axis('off')   
txt = ax1_2.text(0.5, 0.5, 'GEBCO', fontsize = 9, fontweight = 'bold', va = 'center', ha = 'center')
txt.set_clip_on(False)
ax1_2.axis('off')   
txt = ax1_3.text(0.5, 0.5, 'MERIT DEM', fontsize = 9, fontweight = 'bold', va = 'center', ha = 'center')
txt.set_clip_on(False)
ax1_3.axis('off')  

txt = ax2_1.text(0.5, 0.5, 'RCP 26', fontsize = 9, fontweight = 'bold', rotation=270, va = 'center')
txt.set_clip_on(False)
ax2_1.axis('off')   
txt = ax3_1.text(0.5, 0.5, 'RCP 45', fontsize = 9, fontweight = 'bold',rotation=270, va = 'center')
txt.set_clip_on(False)
ax3_1.axis('off')   
txt = ax4_1.text(0.5, 0.5, 'RCP 85', fontsize = 9, fontweight = 'bold', rotation=270, va = 'center')
txt.set_clip_on(False)
ax4_1.axis('off')   

barWidth = 0.85
r = [0, 1, 2, 3, 4, 5]
# Data
raw_data = {'greenBars': [round(i[0], 4) for i in coastal_rcp_26], 'goldBars': [round(i[1], 4) for i in coastal_rcp_26], 'pinkBars': [round(i[2], 4) for i in coastal_rcp_26],\
            'redBars': [round(i[3], 4) for i in coastal_rcp_26], 'darkredBars': [round(i[4], 4) for i in coastal_rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax2_2.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in gebco_rcp_26], 'goldBars': [round(i[1], 4) for i in gebco_rcp_26], 'pinkBars': [round(i[2], 4) for i in gebco_rcp_26],\
            'redBars': [round(i[3], 4) for i in gebco_rcp_26], 'darkredBars': [round(i[4], 4) for i in gebco_rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax2_3.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in merit_rcp_26], 'goldBars': [round(i[1], 4) for i in merit_rcp_26], 'pinkBars': [round(i[2], 4) for i in merit_rcp_26],\
            'redBars': [round(i[3], 4) for i in merit_rcp_26], 'darkredBars': [round(i[4], 4) for i in merit_rcp_26]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax2_4.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax2_4.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax2_4.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax2_4.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax2_4.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in coastal_rcp_45], 'goldBars': [round(i[1], 4) for i in coastal_rcp_45], 'pinkBars': [round(i[2], 4) for i in coastal_rcp_45],\
            'redBars': [round(i[3], 4) for i in coastal_rcp_45], 'darkredBars': [round(i[4], 4) for i in coastal_rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax3_2.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax3_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax3_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax3_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax3_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in gebco_rcp_45], 'goldBars': [round(i[1], 4) for i in gebco_rcp_45], 'pinkBars': [round(i[2], 4) for i in gebco_rcp_45],\
            'redBars': [round(i[3], 4) for i in gebco_rcp_45], 'darkredBars': [round(i[4], 4) for i in gebco_rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax3_3.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax3_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax3_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax3_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax3_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in merit_rcp_45], 'goldBars': [round(i[1], 4) for i in merit_rcp_45], 'pinkBars': [round(i[2], 4) for i in merit_rcp_45],\
            'redBars': [round(i[3], 4) for i in merit_rcp_45], 'darkredBars': [round(i[4], 4) for i in merit_rcp_45]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax3_4.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax3_4.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax3_4.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax3_4.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax3_4.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in coastal_rcp_85], 'goldBars': [round(i[1], 4) for i in coastal_rcp_85], 'pinkBars': [round(i[2], 4) for i in coastal_rcp_85],\
            'redBars': [round(i[3], 4) for i in coastal_rcp_85], 'darkredBars': [round(i[4], 4) for i in coastal_rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax4_2.bar(r, greenBars, color = 'forestgreen', width = barWidth, label = '< -5%')
ax4_2.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth, label = '-5% < AND < -10%')
ax4_2.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth, label = '-10% < AND < -25%')
ax4_2.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth, label = '-25% < AND < -50%')
ax4_2.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth, label = '> -50%')

raw_data = {'greenBars': [round(i[0], 4) for i in gebco_rcp_85], 'goldBars': [round(i[1], 4) for i in gebco_rcp_85], 'pinkBars': [round(i[2], 4) for i in gebco_rcp_85],\
            'redBars': [round(i[3], 4) for i in gebco_rcp_85], 'darkredBars': [round(i[4], 4) for i in gebco_rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax4_3.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax4_3.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax4_3.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax4_3.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax4_3.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

raw_data = {'greenBars': [round(i[0], 4) for i in merit_rcp_85], 'goldBars': [round(i[1], 4) for i in merit_rcp_85], 'pinkBars': [round(i[2], 4) for i in merit_rcp_85],\
            'redBars': [round(i[3], 4) for i in merit_rcp_85], 'darkredBars': [round(i[4], 4) for i in merit_rcp_85]}
df = pd.DataFrame(raw_data)
greenBars = df['greenBars'].values.tolist()
goldBars = df['goldBars'].values.tolist()
pinkBars = df['pinkBars'].values.tolist()
redBars = df['redBars'].values.tolist()
darkredBars = df['darkredBars'].values.tolist()
# Create Bars
ax4_4.bar(r, greenBars, color = 'forestgreen', width = barWidth)
ax4_4.bar(r, goldBars, bottom = greenBars, color = 'gold', width = barWidth)
ax4_4.bar(r, pinkBars, bottom = [i + j for i, j in zip(greenBars, goldBars)], color = 'hotpink', width = barWidth)
ax4_4.bar(r, redBars, bottom = [i + j + k for i, j, k in zip(greenBars, goldBars, pinkBars)], color = 'red', width = barWidth)
ax4_4.bar(r, darkredBars, bottom = [i + j + k + l for i, j, k, l in zip(greenBars, goldBars, pinkBars, redBars)], color = 'darkred', width = barWidth)

# Custom x axis
names = ('2050','2100','2200','2300','2400', '2500')
ax2_2.set_xticks(r)
ax2_2.set_xticklabels([])
ax2_2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_2.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax2_2.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax2_2.set_ylim(1, 0)

ax2_3.set_xticks(r)
ax2_3.set_xticklabels([])
ax2_3.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_3.set_yticklabels([])
ax2_3.set_ylim(1, 0)

ax2_4.set_xticks(r)
ax2_4.set_xticklabels([])
ax2_4.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax2_4.set_yticklabels([])
ax2_4.set_ylim(1, 0)

# Custom x axis
ax3_2.set_xticks(r)
ax3_2.set_xticklabels([])
ax3_2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax3_2.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax3_2.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax3_2.set_ylim(1, 0)

ax3_3.set_xticks(r)
ax3_3.set_xticklabels([])
ax3_3.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax3_3.set_yticklabels([])
ax3_3.set_ylim(1, 0)

ax3_4.set_xticks(r)
ax3_4.set_xticklabels([])
ax3_4.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax3_4.set_yticklabels([])
ax3_4.set_ylim(1, 0)

# Custom x axis
ax4_2.set_xticks(r)
ax4_2.set_xticklabels(names, fontsize = 7)
ax4_2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax4_2.set_yticklabels(['20', '40', '60', '80', '100'], fontsize = 7)
ax4_2.set_ylabel("Affected SRMs [% of all SRMs]", fontsize = 7)
ax4_2.set_ylim(1, 0)
ax4_2.set_xlabel("Year", fontsize = 7)

ax4_3.set_xticks(r)
ax4_3.set_xticklabels(names, fontsize = 7)
ax4_3.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax4_3.set_yticklabels([])
ax4_3.set_ylim(1, 0)
ax4_3.set_xlabel("Year", fontsize = 7)

ax4_4.set_xticks(r)
ax4_4.set_xticklabels(names, fontsize = 7)
ax4_4.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax4_4.set_yticklabels([])
ax4_4.set_ylim(1, 0)
ax4_4.set_xlabel("Year", fontsize = 7)

h, l = ax4_2.get_legend_handles_labels() 
ax6_1.legend(h, l, fontsize = 7, frameon = False, ncol = 5, loc = 'lower center', title = "IFGW difference compared to 2000 [%]", title_fontsize = 8)
ax6_1.axis('off')

canvas = FigureCanvas(fig)
canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_urban_diff_2000.png'), dpi=300)
#canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_agriculture_diff_2000.png'), dpi=300)
#canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_nature_diff_2000.png'), dpi=300)
#canvas.print_figure(os.path.join(data_out_dir, 'hist_DEM_ALL_diff_2000.png'), dpi=300)



        

"""  -------------------------------------------------------------------------
            Create the plots comparing DEM results for each RCP
--------------------------------------------------------------------------  """

import pandas as pd

#   read in the results and land use csv
lu_pop_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\_land_use_population_GDP.csv')
rcp_26_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_26_results.csv')
rcp_45_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_45_results.csv')
rcp_85_csv = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\RCP_85_results.csv')

merit_rcp_26, coastal_rcp_26, gebco_rcp_26 = [], [], []
merit_rcp_45, coastal_rcp_45, gebco_rcp_45 = [], [], []
merit_rcp_85, coastal_rcp_85, gebco_rcp_85 = [], [], []

#   loop through each row of the lu_pop_csv
for ts in ['IFGW_me_2050', 'IFGW_me_2100', 'IFGW_me_2200', 'IFGW_me_2300', 'IFGW_me_2400', 'IFGW_me_2500']:
    ts_cnt = [0, 0, 0, 0, 0]
    for index, row in lu_pop_csv.iterrows():
        if row['urban'] > 0. and row['pop_total'] > 100000:
            #   open the rcp result list and then append to the merit list
            row_rcp_26 = rcp_26_csv.loc[rcp_26_csv['cs_srm_id'] == row['id_srm']]
            val_2000_me = row_rcp_26['IFGW_me_2000'].values[0]
            ifgw = row_rcp_26[ts].values[0] - val_2000_me
            if val_2000_me is not None:
                if ifgw > -5.0:
                    ts_cnt[0] += 1
                elif ifgw < -5.0 and ifgw > -10.0:
                    ts_cnt[1] += 1
                elif ifgw < -10.0 and ifgw > -25.0:
                    ts_cnt[2] += 1
                elif ifgw < -25.0 and ifgw > -50.0:
                    ts_cnt[3] += 1
                elif ifgw < -50.0 and ifgw > -100.0:
                    ts_cnt[4] += 1        








import fiona
import rasterio
from rasterio.mask import mask
#from shapely.geometry import mapping
import numpy as np
import pandas as pd
import shapely.wkt
import shapely.geometry

from shapely.geometry import shape, mapping
from shapely.ops import unary_union
from shapely.ops import cascaded_union
from shapely.geometry import Polygon, Point

#import geopandas as gpd
#shapefile = gpd.read_file(r'g:\Water_Nexus\_A4_paper\_data_output\coast_buffer_10km.shp')


shapefile = fiona.open(r'g:\Water_Nexus\_A4_paper\_data_output\coast_buffer_10km.shp')
#shapefile = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\coast_buffer_10km.csv')
full_dataset = rasterio.open(r'g:\_ORIGINAL_DATA\Land Cover GLC-SHARE\glc_shv10_DOM.Tif')
population_dataset = rasterio.open(r'g:\_ORIGINAL_DATA\_GPW_population_2020\gpw_v4_population_count_rev11_2020_2pt5_min.asc')
gdp_dataset = rasterio.open(r'g:\_ORIGINAL_DATA\GDP_2018\GDP_2015_5arcmin_v2.tiff')


ts_name_lst = [['RCP_26_I_1', 'RCP_26_2050_diff_ifgw'], ['RCP_26_I_2', 'RCP_26_2100_diff_ifgw'], ['RCP_26_I_3', 'RCP_26_2200_diff_ifgw'],\
               ['RCP_26_I_4', 'RCP_26_2300_diff_ifgw'], ['RCP_26_I_5', 'RCP_26_2400_diff_ifgw'], ['RCP_26_I_6', 'RCP_26_2500_diff_ifgw'],\
               ['RCP_45_I_1', 'RCP_45_2050_diff_ifgw'], ['RCP_45_I_2', 'RCP_45_2100_diff_ifgw'], ['RCP_45_I_3', 'RCP_45_2200_diff_ifgw'],\
               ['RCP_45_I_4', 'RCP_45_2300_diff_ifgw'], ['RCP_45_I_5', 'RCP_45_2400_diff_ifgw'], ['RCP_45_I_6', 'RCP_45_2500_diff_ifgw'],\
               ['RCP_85_I_1', 'RCP_85_2050_diff_ifgw'], ['RCP_85_I_2', 'RCP_85_2100_diff_ifgw'], ['RCP_85_I_3', 'RCP_85_2200_diff_ifgw'],\
               ['RCP_85_I_4', 'RCP_85_2300_diff_ifgw'], ['RCP_85_I_5', 'RCP_85_2400_diff_ifgw'], ['RCP_85_I_6', 'RCP_85_2500_diff_ifgw']]


out_csv_cols = ['RCP', 'TIME', 'POP_TOTAL_1', 'POP_TOTAL_2', 'POP_TOTAL_3', 'POP_TOTAL_4', 'POP_TOTAL_5', 'GDP_TOTAL_1',\
                'GDP_TOTAL_2', 'GDP_TOTAL_3', 'GDP_TOTAL_4', 'GDP_TOTAL_5']
out_csv_data = []

for ts_name in ts_name_lst:
    print(ts_name[0])
    lst_1, lst_2, lst_3, lst_4, lst_5 = [], [], [], [], []
    for feature in shapefile:
        #print(feature['properties']['id_srm'])
        if feature['properties'][ts_name[0]] >= -5.:
            lst_1.append(feature['properties']['id_srm'])
        elif feature['properties'][ts_name[0]] < -5. and feature['properties'][ts_name[0]] >= -10.:
           lst_2.append(feature['properties']['id_srm'])      
        elif feature['properties'][ts_name[0]] < -10. and feature['properties'][ts_name[0]] >= -25.:
            lst_3.append(feature['properties']['id_srm'])   
        elif feature['properties'][ts_name[0]] < -25. and feature['properties'][ts_name[0]] >= -50.:
            lst_4.append(feature['properties']['id_srm'])    
        elif feature['properties'][ts_name[0]] >= -100.:
            lst_5.append(feature['properties']['id_srm'])   

    #   create the merged geometry
    polys_1 = []
    for srm in lst_1:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_1.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_1.append(Polygon(feature['geometry']['coordinates'][0]))

    #   create the merged geometry
    polys_2 = []
    for srm in lst_2:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_2.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_2.append(Polygon(feature['geometry']['coordinates'][0]))    

    #   create the merged geometry
    polys_3 = []
    for srm in lst_3:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_3.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_3.append(Polygon(feature['geometry']['coordinates'][0]))
                     
    #   create the merged geometry
    polys_4 = []
    for srm in lst_4:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_4.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_4.append(Polygon(feature['geometry']['coordinates'][0]))
                     
    #   create the merged geometry
    polys_5 = []
    for srm in lst_5:
        for feature in shapefile:
             if feature['properties']['id_srm'] == srm:
                 if len(feature['geometry']['coordinates']) > 1:
                     for i in range(len(feature['geometry']['coordinates'])):
                         try:
                             polys_5.append(Polygon(feature['geometry']['coordinates'][i][0]))
                         except TypeError:
                             pass
                 else:
                     polys_5.append(Polygon(feature['geometry']['coordinates'][0]))

    #   merge the geometry 
    merged_1 = unary_union(polys_1)
    out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_1)], crop=True)
    out_image[out_image < 0] = 0
    total_population_1 = np.sum(out_image)
    out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_1)], crop=True)
    out_image[out_image < 0] = 0
    total_gdp_1 = np.nansum(out_image)
    
    merged_2 = unary_union(polys_2)
    out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_2)], crop=True)
    out_image[out_image < 0] = 0
    total_population_2 = np.sum(out_image)
    out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_2)], crop=True)
    out_image[out_image < 0] = 0
    total_gdp_2 = np.nansum(out_image)
    
    merged_3 = unary_union(polys_3)
    out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_3)], crop=True)
    out_image[out_image < 0] = 0
    total_population_3 = np.sum(out_image)
    out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_3)], crop=True)
    out_image[out_image < 0] = 0
    total_gdp_3 = np.nansum(out_image)
    
    if len(polys_4) > 0:
        merged_4 = unary_union(polys_4)
        out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_4)], crop=True)
        out_image[out_image < 0] = 0
        total_population_4 = np.sum(out_image)
        out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_4)], crop=True)
        out_image[out_image < 0] = 0
        total_gdp_4 = np.nansum(out_image)
    else:
        total_population_4 = 0
        total_gdp_4 = 0

    if len(polys_5) > 0:    
        merged_5 = unary_union(polys_5)
        out_image, out_transform = mask(population_dataset, [shapely.geometry.mapping(merged_5)], crop=True)
        out_image[out_image < 0] = 0
        total_population_5 = np.sum(out_image)
        out_image, out_transform = mask(gdp_dataset, [shapely.geometry.mapping(merged_5)], crop=True)
        out_image[out_image < 0] = 0
        total_gdp_5 = np.nansum(out_image)
    else:
        total_population_5 = 0
        total_gdp_5 = 0

    rcp = ts_name[0].split('_I')[0]
    time = int(ts_name[1].split('_')[2])

    out_csv_data.append([rcp, time, total_population_1, total_population_2, total_population_3, total_population_4, total_population_5,\
                         total_gdp_1, total_gdp_2, total_gdp_3, total_gdp_4, total_gdp_5])


    
my_df = pd.DataFrame(out_csv_data)
my_df.to_csv(r'g:\Water_Nexus\_A4_paper\_data_output\POPULATION_GDP_affected.csv', index=False, header=out_csv_cols)
















"""
# Define a polygon feature geometry with one attribute
schema = {
    'geometry': 'Polygon',
    'properties': {'id': 'int'},
}

# Write a new Shapefile
with fiona.open(r'g:\Water_Nexus\_A4_paper\_data_output\my_shp2.shp', 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here
    c.write({
        'geometry': mapping(merged_1),
        'properties': {'id': 123},
    })    


import xarray as xr

in_data = xr.open_dataset(r'g:\_ORIGINAL_DATA\GDP_2018\GDP_PPP_1990_2015_5arcmin_v2.nc')

gdp = in_data.sel(time = 2015)['GDP_PPP'].values


#   First transform the .mat grid into tiff 
gw_rch_dir_out = r'g:\_ORIGINAL_DATA\GDP_2018\GDP_2015_5arcmin_v2.tiff'

xmin, ymin, xmax, ymax = -180., -90., 180., 90.
nrows, ncols = np.shape(gdp)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(gw_rch_dir_out, gdp.shape[1], gdp.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(gdp)
output_raster.FlushCache()
output_raster = None
"""