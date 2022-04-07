# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 11:58:49 2017

@author: daniel
"""

import psycopg2
import math
import utm
from collections import Counter
import numpy as np
import gdal
import pandas as pd
import os

#import sys
#sys.path.append(r'g:\Water_Nexus\_A3\_scripts_py3')
#import DTBSE_sed_thick_functions_py3  # old sediment thickness estimation script with all functions



#   Define a class that enables access to a database and is able to perform
#   a set of functions in the dbase (create or update tables, get data etc..)
class dbase_tools(object):
    
    #   Initiate the object 
    def __init__(self, db_name, db_user, db_host, db_pass):
        self.db_name = db_name
        self.db_user = db_user
        self.db_host = db_host
        self.db_pass = db_pass
    
        #   automatically connect to the database, define the connection string
        conn_string = str("dbname=%s user=%s host=%s password=%s") % (self.db_name, self.db_user, self.db_host, self.db_pass)
        print('Not connecting to ', conn_string)
        #   try to connect to the database
        """
        try:
            self.conn = psycopg2.connect(conn_string)
            print("Successfully connected to database : " + str(self.db_name))
            #   set the cursor
            self.cur = self.conn.cursor()
            True
        except:
            print("An error occurred - unable to connect to database : " + str(self.db_name))   
            False
        """
        
    #   Method that selects data from given table for the id_cs
    def get_data(self, id_cs, db_table):
        #   define the SQL command, execute it and fetch the results
        sql_select_id = "SELECT * FROM %s WHERE id_cs = %s" % (db_table, id_cs)
        self.cur.execute(sql_select_id)
        return self.cur.fetchall()[0]

    #   Method that selects data from given table for the id_cs
    def get_data_csv(self, cs_id, tb_dir):
        #   open the csv file
        tb = pd.read_csv(tb_dir)
        select_id = tb.loc[tb['id_cs'] == cs_id].values[0].tolist()
        del tb
        return select_id
        
    #   Method that selects data from given table for the id_cs
    def get_data_condition(self, columns, condition, db_table):
        #   define the SQL command, execute it and fetch the results
        sql_select_id = "SELECT %s FROM %s WHERE %s" % (columns, db_table, condition)
        self.cur.execute(sql_select_id)
        return self.cur.fetchall()

    #   Method that selects data from given table for the id_cs
    def get_ids_coscat(self, columns, db_table):
        #   define the SQL command, execute it and fetch the results
        sql_select_id = "SELECT DISTINCT %s FROM %s" % (columns, db_table)
        self.cur.execute(sql_select_id)
        return self.cur.fetchall()

    #   Function that reads in the extent of the coastal plain, and the location of the Anchor point.
    def read_cst_plain_anchor_from_db(self, tb_name, cs_id):
        #   create the sql command
        sql_select_cst_plain_anchor = "SELECT cst_plain_width, nasa_point_dist, nasa_point_depth FROM %s WHERE id_cs = %s" % (tb_name, cs_id)
        self.cur.execute(sql_select_cst_plain_anchor)
        read_vals = self.cur.fetchall()[0]
        cst_plain, anchor_point_x, anchor_point_y = read_vals[0], read_vals[1], read_vals[2]
        return cst_plain, anchor_point_x, anchor_point_y        
    
    #   Method to create a list of columns to be inserted into the newly created table.
    #   Number of columns depends on the distance between individual cross-section points
    #   (for each a column is created) and the total distance from the coast in both directions.
    def create_col_string(self, dist, col_cnt, col_type):
        #   create an empty string first, then create a column for each distance 
        #   negative values represent the offshore cross-section points etc.
        col_str = ""
        for i in range(col_cnt + 1):
            dist_to_cs = i - (col_cnt/2)
            if dist_to_cs < 0:
                col_str += 'dist_minus_' + str((abs(dist_to_cs) * dist)) + 'm ' + col_type + ','
            elif dist_to_cs == 0:
                col_str += 'dist_' + str(dist_to_cs * dist) + 'm ' + col_type + ','
            else:
                col_str += 'dist_plus_' + str(dist_to_cs * dist) + 'm ' + col_type + ','
        #  remove the last comma in the string to get correct SQL command
        col_str = col_str[:-1]
        return col_str       

    #   Method that creates a table in the database. Intended for creation of values 
    #   extracted from input data at cross-section points location.
    def create_table(self, tb_name, id_name, tb_cols):
        #   define the SQL command, execute it and commit the changes to the database
        sql_create_table = "CREATE TABLE IF NOT EXISTS %s (%s serial PRIMARY KEY, %s);" % (tb_name, id_name, tb_cols)
        self.cur.execute(sql_create_table)
        self.conn.commit()

    #   Method to insert values into existing table, creates a new row
    def insert_values(self, db_table, db_columns, ins_vals):
        #    create the SQL string to insert the values
        sql_insert = "INSERT INTO %s (%s) VALUES (%s)" % (db_table, db_columns, ins_vals)
        self.cur.execute(sql_insert)
        self.conn.commit()

    #   Method that updates existing records (row) in a database table. 
    def sql_update(self, db_table, db_columns, ins_vals, where_cond):
        #    create the SQL string to update the values
        sql_update = "UPDATE %s SET %s = %s WHERE %s" % (db_table, db_columns, ins_vals, where_cond)
        self.cur.execute(sql_update)
        self.conn.commit()

    #   Find all coastal points in the appropriate database table
    def find_coastal_points(self, cs_pts_table):
        #   create the SQL string
        sql_str_select_cs_pts = "SELECT id_cs, ST_AsText((ST_Dump(%s.geom)).geom), ST_X(geom), ST_Y(geom) \
                                 AS the_POINT_geom FROM %s ORDER BY gid" % (cs_pts_table, cs_pts_table)
        self.cur.execute(sql_str_select_cs_pts)
        return self.cur.fetchall()
        
    #   Find half the angle between the coastal point and closest node if these two overlay
    def find_angle(self, ax, ay, cs_pt_x, cs_pt_y, cx, cy):        
        #   create arrays with the point coordinates
        a = np.array([ax, ay])
        b = np.array([cs_pt_x, cs_pt_y])
        c = np.array([cx, cy])
        #   substract the arrays
        ba = a - b
        bc = c - b
        #   compute the final angle
        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        #   output is the angle in degrees
        return np.degrees(angle)        
        
    #   Find the closest nodes to a coastal point if the coastal point lies on top of the closest node
    def find_help_nodes(self, tb_cstline, tb_nodes, closest_node_geom):
        #   define the SQL string to find the ID of the coastline line feature that intersects the node
        sql_str_select_nodes = "SELECT gid, geom FROM %s WHERE ST_DWithin(ST_GeomFromText('%s'), %s.geom, 1) ORDER BY \
                                ST_Distance(ST_GeomFromText('%s'), %s.geom) LIMIT 1"\
                                % (tb_cstline, closest_node_geom, tb_cstline, closest_node_geom, tb_cstline)
        cur3 = self.cur
        cur3.execute(sql_str_select_nodes)
        cstline = cur3.fetchall()
        #   find all the nodes that lie on the same costline line feature
        sql_str_select_all_nodes = "SELECT gid, ST_AsText(geom) FROM %s WHERE ST_DWithin(%s.geom, (SELECT geom FROM %s WHERE gid = %s), 0.001)"\
                                    % (tb_nodes, tb_nodes, tb_cstline, cstline[0][0])
        cur3.execute(sql_str_select_all_nodes)     
        all_nodes_id = cur3.fetchall()
        return all_nodes_id

    #   Find the closest node to a coastal point
    def find_closest_node(self, cs_pt_geom, tb_cstline, tb_nodes):
        #   define the SQL string
        sql_str_select_node = "SELECT gid, ST_AsText(geom), ST_X(geom), ST_Y(geom) FROM ne_coastline_nodes_test\
                                ORDER BY ne_coastline_nodes_test.geom <-> ST_GeomFromText('%s') LIMIT 1;" % (cs_pt_geom) 
        cur2 = self.cur
        cur2.execute(sql_str_select_node)
        node = cur2.fetchall()[0]
        #   check if the node position is identical to the position of the coastal point
        if cs_pt_geom == node[1]:
            #   get the X and Y coordinates 
            sql_get_coords = "SELECT ST_X(ST_GeomFromText('%s')), ST_Y(ST_GeomFromText('%s'))" % (cs_pt_geom, cs_pt_geom)
            cur2.execute(sql_get_coords)
            cs_pt_coords = cur2.fetchall()[0]
            #   get the two closest nodes
            other_nodes = self.find_help_nodes(tb_cstline, tb_nodes, node[1])
            sql_get_coords_2 = "SELECT ST_X(ST_GeomFromText('%s')), ST_Y(ST_GeomFromText('%s')),\
                                ST_X(ST_GeomFromText('%s')), ST_Y(ST_GeomFromText('%s'))"\
                                % (other_nodes[1][1], other_nodes[1][1], other_nodes[-2][1], other_nodes[-2][1])
            cur2.execute(sql_get_coords_2)
            helper_pt_coords = cur2.fetchall()
            #   and find the angle 
            angle = self.find_angle(helper_pt_coords[0][0], helper_pt_coords[0][1], cs_pt_coords[0], cs_pt_coords[1], helper_pt_coords[0][2], helper_pt_coords[0][3])
            #   give the closest node (2nd in the list, not necessarily closest but it doesnt matter here since we also give the angle)            
            closest_node = other_nodes[1]
            cs_angle = angle
        else:
            closest_node = node
            cs_angle = 90
        #   give the closest node coordinates and the angle 
        return closest_node, cs_angle
        
    #   get the geometry (and translate it to coordinates) from the raw points table based on the ID 
    def get_point_geom(self, cs_points_raw_tb_name):
        #   get the geometry of the cs_point, define the SQL
        sql_get_pt_geom = "SELECT ST_AsText((ST_Dump(%s.geom)).geom), ST_X(geom), ST_Y(geom) AS the_POINT_geom FROM %s \
                              WHERE id_cs = %s" % (cs_points_raw_tb_name, cs_points_raw_tb_name, self.cs_id)        
        cur = self.cur
        cur.execute(sql_get_pt_geom)
        pt_geom_info = cur.fetchall()[0]
        self.pt_lat, self.pt_lon = pt_geom_info[1], pt_geom_info[2]     
        
    #   transform the coordinates of a point to UTM coordinates, get the UTM zone as well
    def wgs_to_utm(self, wgs_lat, wgs_lon):
        point_utm = utm.from_latlon(wgs_lon, wgs_lat)
        utm_zone = str(point_utm[2]) + point_utm[3]
        return point_utm, utm_zone

    #   find out on which side is the ocean (from the elevation data)
    def ocean_side(self, topo_vals, cst_idx):
        #      a) deal with the continental shelf, first look for the "real" coastal point with elev = 0 or as
        #         close to 0 as possible (staying in positive values) which is at the same time closest to the middle of the list of values
        cs_point_gebco = topo_vals[cst_idx]
        if cs_point_gebco < 0:
            ##  loop through the neighbouring points and check if they are equal to or larger than 0
            ##  also check in which direction from the mid-list value the coastal point actually is,
            ##  assign that to cont_direction (continental) - either left or right (east or west)
            for c in range(len(topo_vals) / 2):
                next_val_left = topo_vals[cst_idx - c]
                next_val_right = topo_vals[cst_idx + c]
                if next_val_left >= 0:
                    coastal_point_indx = cst_idx - c
                    cont_direction = 'left'
                    break
                elif next_val_right >= 0:
                    coastal_point_indx = cst_idx + c
                    cont_direction = 'right'
                    break
                else:
                    pass
        ##  if the coastal point itself is = 0 assign it
        elif cs_point_gebco == 0:
            coastal_point_indx = cst_idx
            ##  check which direction is the continental part - 5 pixels from the coastal point
            next_val_left = topo_vals[coastal_point_indx - 25 : coastal_point_indx]
            next_val_right = topo_vals[coastal_point_indx : coastal_point_indx + 25]
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
            for z in range(len(topo_vals) / 2):
                next_val_left = topo_vals[cst_idx - z]
                next_val_right = topo_vals[cst_idx + z]
                if next_val_left < 0:
                    coastal_point_indx = cst_idx - z
                    cont_direction = 'right'
                    break
                elif next_val_right < 0:
                    coastal_point_indx = cst_idx + z
                    cont_direction = 'left'
                    break
                else:
                    pass
        self.cont_position = cont_direction
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
