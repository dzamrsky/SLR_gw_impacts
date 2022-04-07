# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:06:02 2017

@author: daniel
"""

#   import the libraries
import sys
sys.path.append(r'g:\Water_Nexus\_A3\_scripts_py3')
#sys.path.append(r'home/dzamrsky/_A4/_scripts/_fc_scripts')
import ws_datatools_py3 as ws              # database tools
import DTBSE_sed_thick_functions_py3 as db_sed_thick_func  # old sediment thickness estimation script with all functions
import DTBSE_functions_py3 as db_func
import numpy as np
import os
import pandas as pd

#   specify database connection info
dtb_name = "'cs_db_v1'"
dtb_host = "'localhost'"
dtb_user = "'postgres'"
dtb_pass = "'postgres'"

cs_points_raw_tb_name = "ne_10m_coastline_points_5km_sed_thick_est"
cs_nodes_raw_tb_name = "ne_coastline_nodes_test"
cs_points_tb_name = "cs_points"

#   there will be a list (with check boxes) 
read_data_lst = ['topo', 'glim', 'cs_soil_thick', 'cs_soil_type', 'cs_wtd_depth', 'cs_offshore_sed']


csv_dir = r'g:\Water_Nexus\_A4_GUM\_GIS_data'
#csv_dir = r'/projects/0/qt16165/_dzamrsky/_A4_input_data/_csv_input_files'

tb_gebco_topo = os.path.join(csv_dir, 'cs_gebco_2014_avg.csv')
tb_gebco_topo_glc = os.path.join(csv_dir, 'cs_glc.csv')
tb_gebco_topo_nasa_thk = os.path.join(csv_dir, 'cs_nasa_thick_avg.csv')
tb_gebco_topo_pcr_thk = os.path.join(csv_dir, 'cs_pcrglob_thick.csv')
tb_gebco_topo_sed_thk_est = os.path.join(csv_dir, 'cs_sed_thick_est.csv')
tb_gebco_topo_soil_thk = os.path.join(csv_dir, 'cs_soil_thick.csv')
tb_gebco_topo_soil_type = os.path.join(csv_dir, 'cs_soil_type.csv')
tb_gebco_topo_wtd_depth = os.path.join(csv_dir, 'cs_wtd_depth.csv')
tb_gebco_topo_sed_offshore = os.path.join(csv_dir, 'cs_offshore_sed.csv')
tb_gebco_topo_k_soil = os.path.join(csv_dir, 'cs_k_soil.csv')
tb_gebco_topo_drn_rate = os.path.join(csv_dir, 'cs_drn_rate.csv')
tb_gebco_topo_glhymps_1 = os.path.join(csv_dir, 'cs_glhymps_layer_1.csv')
tb_gebco_topo_glhymps_2 = os.path.join(csv_dir, 'cs_glhymps_layer_2.csv')
tb_glim = os.path.join(csv_dir, 'cs_glim_litho_v1.csv')

"""      lon = x coordinate, lat = y coordinate      """

#   define class whose only input are lon, lat coordinates - as manual input now
#   but later as coordinates of a point that user clicks on the map in the web service
class point(object):
    
    #   initialize the object, only two attributes are lat, lon
    def __init__(self, lon, lat, n_pts, max_dst, id_cs = None):
        #   first option is that no cross-section ID is specified, instead click on the global map input coordinates are assigned
        if id_cs is None:
            self.lon = lon        
            self.lat = lat
            self.cs_id = None
        else:
            self.cs_id = id_cs
        #   the rest is the same for both cases
        self.n_pts = n_pts
        self.max_dst = max_dst
        #   get a database connection cursor
        self.db_connect = ws.dbase_tools(dtb_name, dtb_user, dtb_host, dtb_pass)

    #   method that searches in a radius around the point and finds the closest
    #   coastal point in the database, print out the ID of the point, its coordinates
    #   and the distance from the clicked point (in meters). 
    def closest_cs_point(self):
        if self.cs_id is None:
            print('None')
        #   build the SQL command that will first select all the points within a radius of the given coordinates and then get the closest one out of that list
            sql_get_cs_point = "SELECT gid, ST_X(geom), ST_Y(geom), ST_AsText(geom) FROM %s WHERE ST_DWithin(%s.geom, (SELECT ST_MakePoint(%s, %s)) , 1)\
                                ORDER BY %s.geom <-> (SELECT ST_MakePoint(%s, %s)) LIMIT 1 ;" % (cs_points_raw_tb_name, cs_points_raw_tb_name,\
                                self.lon, self.lat, cs_points_raw_tb_name, self.lon, self.lat)  
        else:
            print(str(self.cs_id))
            sql_get_cs_point = "SELECT gid, ST_X(geom), ST_Y(geom), ST_AsText(geom) FROM %s WHERE gid = %s" % (cs_points_raw_tb_name, self.cs_id)
        #   execute the querry and try to get the output
        self.db_connect.cur.execute(sql_get_cs_point)
        try:
            cs_point = self.db_connect.cur.fetchall()[0]
            #   assign the output of the sql
            self.cs_id = cs_point[0]
            self.cs_lon = cs_point[1]
            self.cs_lat = cs_point[2]
            self.cs_geom = cs_point[3]
        except:
            print('No coastal point in 1 degree radius around the point or with the given ID number. Please select new point.')
            
        
    #   read in the all necessary data for the cross-section from the main database
    #   [1:] to skip the id_cs that also gets loaded by the get_data function
    #   adjust_to_ocean - if true then the data is adjusted so that the ocean is always on the right hand side of the cross-section
    def get_input_data(self, adjust_to_ocean = True):
        #print self.cs_id
        if self.cs_id is not None:
            """
            self.dta_topo    = self.db_connect.get_data(self.cs_id, 'cs_gebco_2014_avg')[1:]
            self.dta_glim    = self.db_connect.get_data(self.cs_id, 'cs_glim_litho')[1:]
            self.dta_glc     = self.db_connect.get_data(self.cs_id, 'cs_glc')[1:]
            self.dta_thk_pel = self.db_connect.get_data(self.cs_id, 'cs_nasa_thick_avg')[1:]   #   thickness as estimated by Pelletier et. al (2016)
            self.dta_thk_pcr = self.db_connect.get_data(self.cs_id, 'cs_pcrglob_thick')[1:]    #   thickness as estimated by the PCR model input
            self.dta_thk_est = self.db_connect.get_data(self.cs_id, 'cs_sed_thick_est')[1:]    #   thickness as estimated by my own estimation method
 
            self.dta_soil_thk = self.db_connect.get_data(self.cs_id, 'cs_soil_thick')[1:]      #    soil thickness from soilgrids estimation
            self.dta_soil_type = self.db_connect.get_data(self.cs_id, 'cs_soil_type')[1:]      #    soil type derived from soilgrids
            self.dta_wtd_depth = self.db_connect.get_data(self.cs_id, 'cs_wtd_depth')[1:]      #    water table depth from Earth2Observe
            self.dta_offshore_sed = self.db_connect.get_data(self.cs_id, 'cs_offshore_sed')[1:]#    offshore sediment type from the 

            self.dta_pcr_rch = self.db_connect.get_data(self.cs_id, 'cs_rch_pcr')[1:]          #    GW recharge from PCRGLOBWB
            self.dta_watergap_rch = self.db_connect.get_data(self.cs_id, 'cs_rch_watergap')[1:]
            self.dta_p_min_et = self.db_connect.get_data(self.cs_id, 'cs_p_min_et')[1:]
            self.dta_k_soil = self.db_connect.get_data(self.cs_id, 'cs_k_soil')[1:]
            self.dta_drn_rate = self.db_connect.get_data(self.cs_id, 'cs_drn_rate')[1:]
            
            self.dta_riv_cond_01 = self.db_connect.get_data(self.cs_id, 'cs_riv_cond_0_1')[1:]            
            self.dta_riv_cond_05 = self.db_connect.get_data(self.cs_id, 'cs_riv_cond_0_5')[1:]                     
            self.dta_riv_cond_10 = self.db_connect.get_data(self.cs_id, 'cs_riv_cond_1_0')[1:]                     
            self.dta_riv_width = self.db_connect.get_data(self.cs_id, 'cs_gaia_riv_width')[1:]            
            self.dta_riv_bot_elev = self.db_connect.get_data(self.cs_id, 'cs_riv_bot_surf_elev')[1:]                     
            self.dta_riv_head_elev = self.db_connect.get_data(self.cs_id, 'cs_riv_head_surf_elev')[1:]   
            self.dta_glhymps_top_lay = self.db_connect.get_data(self.cs_id, 'cs_glhymps_layer_1')[1:]                     
            self.dta_glhymps_bot_lay = self.db_connect.get_data(self.cs_id, 'cs_glhymps_layer_2')[1:]   
            """
            
            self.dta_topo = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo)[1:] 
            self.dta_glim = self.db_connect.get_data_csv( self.cs_id, tb_glim)[1:] 
            self.dta_glc = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_glc)[1:] 
            self.dta_thk_pel = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_nasa_thk)[1:]    #   thickness as estimated by Pelletier et. al (2016)
            self.dta_thk_pcr = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_pcr_thk)[1:]    #   thickness as estimated by the PCR model input
            self.dta_thk_est = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_sed_thk_est)[1:]     #   thickness as estimated by my own estimation method
            self.dta_soil_thk = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_soil_thk)[1:]       #    soil thickness from soilgrids estimation
            self.dta_soil_type = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_soil_type)[1:]      #    soil type derived from soilgrids
            self.dta_wtd_depth = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_wtd_depth)[1:]      #    water table depth from Earth2Observe
            self.dta_offshore_sed = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_sed_offshore)[1:] #    offshore sediment type from the 
            self.dta_k_soil = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_k_soil)[1:] 
            self.dta_drn_rate = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_drn_rate)[1:] 
            self.dta_glhymps_top_lay = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_glhymps_1)[1:]                      
            self.dta_glhymps_bot_lay = self.db_connect.get_data_csv(self.cs_id, tb_gebco_topo_glhymps_2)[1:]    
            
            
            if adjust_to_ocean:
                #   check the position of the land, use the function from the DTBSE_sed_thick_functions
                #slope_est = db_sed_thick_func.est_slope(self.cs_id, self.dta_topo, self.dta_glim, self.db_connect.cur, 400)
                slope_est = db_sed_thick_func.est_slope(self.cs_id, self.dta_topo, self.dta_glim, 400)
                self.land_pos = slope_est[-1]
                #   transform to position of the ocean
                if self.land_pos == 'right':
                    #   reverse the list of values so it seems like ocean is on the right hand side
                    self.dta_topo = list(reversed(self.dta_topo))
                    self.dta_thk_pel = list(reversed(self.dta_thk_pel))
                    self.dta_thk_pcr = list(reversed(self.dta_thk_pcr))
                    self.dta_glim = list(reversed(self.dta_glim))
                    self.dta_glc = list(reversed(self.dta_glc))
                    self.dta_soil_thk = list(reversed(self.dta_soil_thk))
                    self.dta_soil_type = list(reversed(self.dta_soil_type))
                    self.dta_wtd_depth = list(reversed(self.dta_wtd_depth))
                    self.dta_offshore_sed = list(reversed(self.dta_offshore_sed)) 
                    #self.dta_pcr_rch = list(reversed(self.dta_pcr_rch)) 
                    #self.dta_watergap_rch = list(reversed(self.dta_watergap_rch)) 
                    #self.dta_p_min_et = list(reversed(self.dta_p_min_et)) 
                    self.dta_k_soil = list(reversed(self.dta_k_soil)) 
                    self.dta_drn_rate = list(reversed(self.dta_drn_rate)) 
                    #self.dta_riv_cond_01 = list(reversed(self.dta_riv_cond_01))
                    #self.dta_riv_cond_05 = list(reversed(self.dta_riv_cond_05))                  
                    #self.dta_riv_cond_10 = list(reversed(self.dta_riv_cond_10))                   
                    #self.dta_riv_width = list(reversed(self.dta_riv_width))       
                    #self.dta_riv_bot_elev = list(reversed(self.dta_riv_bot_elev))          
                    #self.dta_riv_head_elev = list(reversed(self.dta_riv_head_elev))  
                    self.dta_glhymps_top_lay = list(reversed(self.dta_glhymps_top_lay))                 
                    self.dta_glhymps_bot_lay = list(reversed(self.dta_glhymps_bot_lay))     
                    
                else:
                    pass     
            try:
                #   get the info from the sediment thickness estimation dbase table
                self.cs_plain_width = float(self.dta_thk_est[2])
                self.anchor_dist_to_cst = 200.0 - float(self.dta_thk_est[3])
                self.anchor_depth = float(self.dta_thk_est[4])
                self.sed_thick_est_avg = float(self.dta_thk_est[-2])
                self.sed_thick_est_stdev = float(self.dta_thk_est[-1]) 
            except TypeError:       #   in case the sed_thick_est is empty - no estimation was successful so no final value
                self.sed_thick_est_avg = None
                print('No sediment thickness estimation for cross-section ID : ' + str(self.cs_id))

    #   get the thickness (below the surface level) and assign it to the database for the cross-section
    #   to save time first check if the information for the given cross-section is already present in the database
    def get_crs_sec_pts(self):
        if self.cs_id is not None and self.sed_thick_est_avg is not None:
            #   check that the cross-section point for the given coastal point are not in the dbase already
            sql_check = "SELECT * FROM %s WHERE id_cs = %s" % (cs_points_tb_name, self.cs_id)
            self.db_connect.cur.execute(sql_check)
            check = self.db_connect.cur.fetchall()
            #   if the check is not empty than it is already in the database
            if check:
                print('The thickness is already estimated for coastal profile ID = ' + str(self.cs_id) + ', loading from database..')
                print('The thickness is:    ' + str(self.sed_thick_est_avg) + 'm')
            #   otherwise find the thickness estimation along the profile and write it to the database
            else:
                print('Calculating the thickness estimation for coastal profile ID = ' + str(self.cs_id) + ', loading from database..')
                #   1) find the closest node to the coastal point
                closest_node_wgs = self.db_connect.find_closest_node(self.cs_geom, cs_points_raw_tb_name, cs_nodes_raw_tb_name)
                self.closest_node_lon_wgs = closest_node_wgs[0][2]
                self.closest_node_lat_wgs = closest_node_wgs[0][3]
                #   2) transform all the coordinates into UTM and calculate the distances 
                self.cs_point_utm = self.db_connect.wgs_to_utm(self.cs_lon, self.cs_lat)
                self.closest_node_utm = self.db_connect.wgs_to_utm(self.closest_node_lon_wgs, self.closest_node_lat_wgs)  
                self.utm_zone_num = self.cs_point_utm[0][2]
                self.utm_zone_let = self.cs_point_utm[0][3]
                #print('UTM zone : ' + str(self.utm_zone_num) + self.utm_zone_let)
                #   3) check that zone_num is between 1 and 60, case for areas that are split between the zones 1 and 60 (e.g. certain islands in Pacific..)
                if self.utm_zone_num > 60:
                    self.utm_zone_num = self.utm_zone_num - 60
                #    4) set the secondary ID to 1, this id number starts at 1 for each cross-section
                cs_point_id = 1
                #   5) loop through all the cross-section points and find the estimated sediment thickness
                for a in range(self.n_pts):
                    b_length = self.max_dst - a * (self.max_dst / self.n_pts)       #   get the distance of the cross-section point
                    #print b_length
                #   6) only estimate the thickness for points that are within the coastal plain
                    if b_length <= (self.cs_plain_width * 1000):                    #   * 1000 to convert to meters
                #   7) check if the current cross-section point is either the anchor point or the end of the coastal plain
                        if b_length == (self.anchor_dist_to_cst * 1000):
                            nasa_point_bool = "'1'"
                        else:
                            nasa_point_bool = "'0'"
                        if b_length == (self.anchor_dist_to_cst * 1000):
                            cst_plain_end_point_bool = "'1'"
                        else:
                            cst_plain_end_point_bool = "'0'"
                #   8) calculate the UTM coordinates of the cross-section points on both sides of the cross-section
                        c_coords_utm = db_func.coord_from_triangle(self.cs_point_utm[0][0], self.cs_point_utm[0][1], self.closest_node_utm[0][0], self.closest_node_utm[0][1], b_length)
                #   9) transform the coordinates back to wgs84, position of the equator is taken into account
                        c_coords_wgs = db_func.equator_position_angles(c_coords_utm[1], c_coords_utm[0], c_coords_utm[3], c_coords_utm[2], self.utm_zone_num, self.utm_zone_let)
                #   10) calculate and assign all values based on the continent position
                        if self.land_pos == 'right':
                            c_coord_wgs_x_land, c_coord_wgs_y_land = c_coords_wgs[0][1], c_coords_wgs[0][0]
                        else:
                            c_coord_wgs_x_land, c_coord_wgs_y_land = c_coords_wgs[1][1], c_coords_wgs[1][0]        
                #    11) get the estimated thickness for the given point, based on the distance from the coast and the thickness indicated by the NASA dataset
                        pt_nasa_thick = self.dta_thk_pel[a]
                        #print self.cs_id, self.dta_thk_pel[a-20: a], pt_nasa_thick
                #    12) for the points located between the anchor point and the end of the coastal plain take the Pelletier thickness (even if = 50m). The STDEV value is 0.0
                        if b_length >= (self.anchor_dist_to_cst * 1000):
                            sed_thick_avg = pt_nasa_thick
                            last_nasa_thick_rel_sea_level = self.dta_topo[a] - pt_nasa_thick   #   bottom of sediment relative to sea level
                            last_nasa_dist = a * 0.5
                            sed_thick_stdev = 0.0                   
                #   13) if the nasa thickness is higher than 50m, calculate the thickness (difference between elevation and
                #       the estimation thickness line) for the given point. The estimation line is given by two points, the
                #       last_nasa_thick value that represents the last value lowet than 50m, and the estimated thickness at  the coastline.
                        else:
                            try:
                                x1 = last_nasa_dist
                                y1 = last_nasa_thick_rel_sea_level                      #   Y coordinate (relative to sea level) of the last nasa_thick value
                            except:
                                x1 = self.anchor_dist_to_cst
                                y1 = self.anchor_depth
                            #x1 = last_nasa_dist
                            #y1 = last_nasa_thick_rel_sea_level
                            #print x1, y1, self.anchor_dist_to_cst, self.anchor_depth
                            x2 = 200.                                               #   X coordinate of the coastline (constant at 200.m)
                            y2 = self.dta_topo[400] - self.sed_thick_est_avg        #   Y coordinate at the coastline, estimated
                            x3 = 200. - b_length / 1000.
                            y3 = y1 + np.round((abs(x1 - x3) / abs(x1 - x2)) * (y2 - y1))
                            sed_thick_avg = self.dta_topo[a] - y3
                            try:
                                sed_thick_stdev = np.round(self.sed_thick_est_stdev - ((self.sed_thick_est_stdev / self.anchor_dist_to_cst) * (200. - x3)))
                            except ZeroDivisionError:
                                continue
                #   14) create the string for the geom of the point
                        point_geom = "ST_SETSRID(ST_MAKEPOINT (%s, %s), 4326)" % (c_coord_wgs_x_land, c_coord_wgs_y_land)
                #   15) insert the values to the database for each cross-section point
                        sql_write_to_cs_points_land = "INSERT INTO cs_points (id_cs, id_point, sed_thick_avg, sed_thick_stdev, x_coord_wgs84,\
                                                       y_coord_wgs84, geom, nasa_point, cst_plain_end) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)"\
                                                       % (self.cs_id, cs_point_id, sed_thick_avg, sed_thick_stdev, c_coord_wgs_x_land, c_coord_wgs_y_land,\
                                                          point_geom, nasa_point_bool, cst_plain_end_point_bool)
                        self.db_connect.cur.execute(sql_write_to_cs_points_land)
                        self.db_connect.conn.commit()
                        print('The thickness is:    ' + str(sed_thick_avg) + 'm')
                #   16) increase the value + 1
                        cs_point_id = cs_point_id + 1
































    