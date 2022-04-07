#-------------------------------------------------------------------------------
# http://gis.stackexchange.com/questions/75407/how-to-create-linestrings-with-a-definite-angle-and-length-that-are-fixed-to-a
# http://gis.stackexchange.com/questions/59169/how-to-draw-perpendicular-lines-in-qgis/59196#59196
# Export postgreSQL table to csv file:
#   ->  COPY cs_gebco2014 TO 'g:\Water_Nexus\Databases\coastal_dbase\cs_gebco2014_test.csv ' DELIMITER ',' CSV HEADER;
# http://www.gis.usu.edu/~chrisg/python/2009/
# http://www.movable-type.co.uk/scripts/latlong.html
#
#
#
#-------------------------------------------------------------------------------

"""
DELETE FROM public.cs_gebco_2014;
DELETE FROM public.cs_glc;
DELETE FROM public.cs_hwsd;
DELETE FROM public.cs_seabed_litho;
DELETE FROM public.cs_seabed_thick;
DELETE FROM public.cs_aq_thick;
DELETE FROM public.cs_aq_kval;
DELETE FROM public.cs_gebco_2014_045_deg;
DELETE FROM public.cs_gebco_2014_135_deg;
DELETE FROM public.cs_gebco_2014_avg;
DELETE FROM public.cs;

"""

""" ************************************************************************ """
"""                     IMPORT ALL NECESSARY LIBRARIES                       """
""" ************************************************************************ """

import psycopg2
#from psycopg2 import IntegrityError
#import sys
#import os
#import ogr
import gdal
import math
import utm
#import csv
from time import time
#import matplotlib.pyplot as plt
from collections import Counter

""" ************************************************************************ """
"""                     DEFINE ALL INPUT FILES AND DIRS                      """
""" ************************************************************************ """

##  format is: dir, table name in postgreSQL database
input_data = [#[r'g:\_ORIGINAL_DATA\GEBCO 2014\GEBCO_2014_2D.nc', 'cs_gebco_2014', 'avg', 'int']]#,
             [r'g:\_ORIGINAL_DATA\PCR_GLOBWB\pcr_globwb_thickness.tif', 'cs_aq_thick', 'avg', 'real']]#,
             #[r'g:\_ORIGINAL_DATA\_GEOLOGY_DATA\_NASA_Global_soil_regolith_sediment_thickness\data\average_soil_and_sedimentary-deposit_thickness.tif', 'nasa_thick_avg', 'avg', 'real']]


#   define the table in the database that contains the coastline points shapefile
#   this is not the cs table where the info is written to!
coastline_points_raw_table = '_ne_points_deltas'#'ne_points_5km'#'ne_points_5km_test_raster_edge'##'ne_coatline_points_5km_test_avg'#'ne_coatline_points_5km_test_avg'

#  determine how many points to create (25 means point per 1km)
#  and the distance of the farthest point from point_5km
n_points = 799
max_dist = 399500

#   define the number of points and distance from the cross-section point
n_points_avg = 5
avg_dist = 2500


""" ************************************************************************ """
"""                     DEFINE ALL FUNCTIONS NECESSARY                       """
""" ************************************************************************ """

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
    px = int((point_coord_x - gt[0]) / gt[1]) #x pixel
    py = int((point_coord_y - gt[3]) / gt[5]) #y pixel
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

#x = 180.75
#y = -16.9780

#read_raster_val(rb, gt, 43200, 21600, x, y)

""" ************************************************************************ """
"""            WRITE THE INPUT DATA INTO POSTGRESQL DATABASE                 """
""" ************************************************************************ """

#del dtbase_con, dtbase_cur, dtbase_cur2, dtbase_cur3, dtbase_cur4

##   connect to the database
db_name = "'cs_db_deltas'"
db_user = "'postgres'"
#db_user = "'super_postgres'"
db_host = "'localhost'"
db_pass = "'postgres'"
#db_pass = "'Sup3rm4n'"
##   connect to coastal database
dtbase_con = connect_to_dtbase(db_name, db_user, db_host, db_pass)
#dtbase_cur = dtbase_con[1]

##  start the loop for importing data into database, raster by raster
for raster_file in input_data:

    #raster_file = input_data[0]
    print(raster_file[0])

    ##  first open the raster file via GDAL
    raster_in = gdal.Open(raster_file[0])

    ##  assign the geotransform, band and noval
    gt = raster_in.GetGeoTransform()
    rb = raster_in.GetRasterBand(1)
    noval = rb.GetNoDataValue()
    rx = raster_in.RasterXSize
    ry = raster_in.RasterYSize

    ##  define the database cursor - necessary to move around the dbase tables
    dtbase_cur = dtbase_con[1]

    ##   first run query to find all the 5km points
    sql_string_select_5km_points = "SELECT id_cs, ST_AsText((ST_Dump(%s.geom)).geom),\
                                    ST_X(geom), ST_Y(geom) AS the_POINT_geom FROM %s ORDER BY gid"\
                                    % (coastline_points_raw_table, coastline_points_raw_table)

    ##  run the sql and get all the coastline points from the raw table
    dtbase_cur.execute(sql_string_select_5km_points)
    points_5km = dtbase_cur.fetchall()

    ##  create table in database that will contain the extracted values from raster
    ##  dist_columns runs the function to get a list of columns for the created table
    dist_columns = create_column_string(2 * n_points, 500, raster_file[3])
    sql_create_table(dtbase_con[0], dtbase_cur, raster_file[1], 'id_cs', dist_columns)

    ##  delete the cursor
    del dtbase_cur

    ##  create an empty list for the csv file output
    to_csv = [[],[],[],[],[]]
    to_csv_2 = [[],[],[],[],[]]
    to_csv_3 = [[],[],[],[],[]]
    to_csv_4 = [[],[],[],[],[]]

    ##  start measuring the time for writing the raster info to database
    t0 =  time()

    #pbar = ProgressBar()

    ##   run loop for each point found in the points_5km table
    for point_5km in points_5km[5000:]:

        #point_5km = points_5km[0]

        ##  get id and coordinates
        point_id = point_5km[0]
        point_geom = point_5km[1]
        point_x, point_y = point_5km[2], point_5km[3]

        #   check if the point already exists in the database
        dtbase_cur0 = dtbase_con[1]
        sql_check_if_exists = "SELECT id_cs FROM %s WHERE id_cs = %s" % (raster_file[1], point_id)
        dtbase_cur0 = dtbase_con[1]
        dtbase_cur0.execute(sql_check_if_exists)
        exists = dtbase_cur0.fetchall()
    
        if exists:
            print ('Data for cs_id = ' + str(point_id) + ' already exist, moving to next..')
            continue
        else:
            print ('Extracting data for cs_id = ' + str(point_id))

        ##  create the sql string to find closest node to point_5km
        sql_string_select_points = "SELECT gid, ST_AsText(geom), ST_X(geom), ST_Y(geom) FROM ne_10km_coastline_nodes\
                                    ORDER BY ne_10km_coastline_nodes.geom <->\
                                    ST_GeomFromText('%s') LIMIT 1;" % (point_geom)

        ##  open new database cursor
        dtbase_cur2 = dtbase_con[1]

        ##  run the sql and get results
        dtbase_cur2.execute(sql_string_select_points)
        closest_node = dtbase_cur2.fetchall()

        ##  assign the results of the query to variables below
        closest_node_id = closest_node[0][0]
        closest_node_geom = closest_node[0][1]
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
        dist_to_node = math.sqrt((point_5km_utm[0] - closest_node_utm[0]) ** 2 + (point_5km_utm[1] - closest_node_utm[1]) ** 2)

        ##  extract value on the coastline point
        pixel_val_cs_point = read_raster_val(rb, gt, rx, ry, point_x, point_y)

        ##  average values and different angle cross-section only for GEBCO
        if raster_file[1] == 'cs_gebco_2014':

            ##  do the same for average points around this coastal point
            sum_avg = pixel_val_cs_point
            cnt_avg = 1

            ##  calculate the average value for each cross-section point
            for c in range(n_points_avg):
                avg_point_dist = avg_dist - c * (max_dist/n_points)
                avg_point_coord = avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], point_5km_utm[0], point_5km_utm[1], avg_point_dist)
                ##  transform to wgs coordinates
                avg_point_coord_wgs = equator_position(avg_point_coord[1], avg_point_coord[0], avg_point_coord[3], avg_point_coord[2], zone_num, zone_let)

                ##  get the pixel value for the cross-section avg point
                avg_pixel_val_1 = read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs[0][1], avg_point_coord_wgs[0][0])
                sum_avg += avg_pixel_val_1

                avg_pixel_val_2 = read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs[1][1], avg_point_coord_wgs[1][0])
                sum_avg += avg_pixel_val_2

                cnt_avg += 2

            ##  write the average value into the dbase table
            sql_insert(dtbase_con[0], dtbase_cur2, 'cs_gebco_2014_avg', 'id_cs', str(point_id))
            sql_update(dtbase_con[0], dtbase_cur2, 'cs_gebco_2014_avg', 'dist_0m', int(sum_avg / cnt_avg), 'id_cs = ' + str(point_id))

            ##  write the id_cs point into the database and assign the point_5km value to the dist_0m column
            sql_insert(dtbase_con[0], dtbase_cur2, 'cs_gebco_2014_045_deg', 'id_cs', str(point_id))
            sql_insert(dtbase_con[0], dtbase_cur2, 'cs_gebco_2014_135_deg', 'id_cs', str(point_id))

            sql_update(dtbase_con[0], dtbase_cur2, 'cs_gebco_2014_045_deg', 'dist_0m', pixel_val_cs_point, 'id_cs = ' + str(point_id))
            sql_update(dtbase_con[0], dtbase_cur2, 'cs_gebco_2014_135_deg', 'dist_0m', pixel_val_cs_point, 'id_cs = ' + str(point_id))

        else:
            pass

        ##  insert the cs id to the table of the raster info
        sql_insert(dtbase_con[0], dtbase_cur2, raster_file[1], 'id_cs', str(point_id))

        dtbase_cur3 = dtbase_con[1]

        ##  check if the point location falls on a novalue pixel of the raster
        ##  if yes then get the average/max freq value from neighbouring pixels
        if pixel_val_cs_point == noval:
            ##  get the raster coordinates
            px_cs_point = int((point_x - gt[0]) / gt[1]) #x pixel
            py_cs_point = int((point_y - gt[3]) / gt[5]) #y pixel
            ##  find the nearest value for the point
            near_val_cs_point = find_closest_neighbour(px_cs_point, py_cs_point, rb, noval, rx, ry, raster_file[2])
            ##  insert the value at coastline point to the database (dist_0m)
            sql_update(dtbase_con[0], dtbase_cur3, raster_file[1], 'dist_0m', near_val_cs_point, 'id_cs = ' + str(point_id))
        ##  if not insert to database the extracted value
        else:
            sql_update(dtbase_con[0], dtbase_cur3, raster_file[1], 'dist_0m', pixel_val_cs_point, 'id_cs = ' + str(point_id))

        dtbase_cur4 = dtbase_con[1]

        ##  write the point into the database tables - cs and cs_gebco2014
        ##  check if the point is already written to cs table
        sql_string_check_cs_id = "SELECT id_cs FROM cs WHERE id_cs = %s" % (str(point_id))
        dtbase_cur4.execute(sql_string_check_cs_id)
        id_check = dtbase_cur4.fetchall()

        if id_check:
            #print ('cs id : ' + str(point_id) + ' exists already..')
            pass
        else:
            sql_insert(dtbase_con[0], dtbase_cur4, 'cs', 'id, x_coord, y_coord', str(point_id) + ',' + str(point_x) + ',' + str(point_y))

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

                dtbase_cur5 = dtbase_con[1]

                ##  write the values to the database table
                #   create string so it matches the column name in format dist_plus/minus_XXXXm
                col_name_land = 'dist_plus_' + str(int(b_length)) + 'm'
                col_name_sea = 'dist_minus_' + str(int(b_length)) + 'm'
                #   write the values to the database
                sql_update(dtbase_con[0], dtbase_cur5, raster_file[1], col_name_land, pixel_val_land, 'id_cs = ' + str(point_id))
                sql_update(dtbase_con[0], dtbase_cur5, raster_file[1], col_name_sea, pixel_val_sea, 'id_cs = ' + str(point_id))

                ##  do this only in case of GEBCO
                if raster_file[1] == 'cs_gebco_2014':

                    ##  get the values for different-angle cross-sections
                    c_coords_utm_angle_045 = coord_from_triangle_alfa(point_5km_utm[0], point_5km_utm[1], closest_node_utm[0], closest_node_utm[1], b_length, 45)
                    c_coords_utm_angle_135 = coord_from_triangle_alfa(point_5km_utm[0], point_5km_utm[1], closest_node_utm[0], closest_node_utm[1], b_length, 135)

                    c_coords_wgs_angle_045 = equator_position_angles(c_coords_utm_angle_045[0], c_coords_utm_angle_045[1], c_coords_utm_angle_045[2], c_coords_utm_angle_045[3], zone_num, zone_let)
                    c_coords_wgs_angle_135 = equator_position_angles(c_coords_utm_angle_135[0], c_coords_utm_angle_135[1], c_coords_utm_angle_135[2], c_coords_utm_angle_135[3], zone_num, zone_let)

                    ##  check closest node coords vs coastal point coords.
                    c_coord_wgs_x_land_angle_045, c_coord_wgs_y_land_angle_045, c_coord_wgs_x_sea_angle_045, c_coord_wgs_y_sea_angle_045 = \
                    c_coords_wgs_angle_045[0][1], c_coords_wgs_angle_045[0][0], c_coords_wgs_angle_045[1][1], c_coords_wgs_angle_045[1][0]

                    c_coord_wgs_x_land_angle_135, c_coord_wgs_y_land_angle_135, c_coord_wgs_x_sea_angle_135, c_coord_wgs_y_sea_angle_135 = \
                    c_coords_wgs_angle_135[0][1], c_coords_wgs_angle_135[0][0], c_coords_wgs_angle_135[1][1], c_coords_wgs_angle_135[1][0]

                    pixel_val_land_angle_045 = read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_land_angle_045, c_coord_wgs_y_land_angle_045)
                    pixel_val_sea_angle_045 = read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_sea_angle_045, c_coord_wgs_y_sea_angle_045)

                    pixel_val_land_angle_135 = read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_land_angle_135, c_coord_wgs_y_land_angle_135)
                    pixel_val_sea_angle_135 = read_raster_val(rb, gt, rx, ry, c_coord_wgs_x_sea_angle_135, c_coord_wgs_y_sea_angle_135)

                    sql_update(dtbase_con[0], dtbase_cur5, 'cs_gebco_2014_045_deg', col_name_land, pixel_val_land_angle_045, 'id_cs = ' + str(point_id))
                    sql_update(dtbase_con[0], dtbase_cur5, 'cs_gebco_2014_045_deg', col_name_sea, pixel_val_sea_angle_045, 'id_cs = ' + str(point_id))
                    sql_update(dtbase_con[0], dtbase_cur5, 'cs_gebco_2014_135_deg', col_name_land, pixel_val_land_angle_135, 'id_cs = ' + str(point_id))
                    sql_update(dtbase_con[0], dtbase_cur5, 'cs_gebco_2014_135_deg', col_name_sea, pixel_val_sea_angle_135, 'id_cs = ' + str(point_id))

                    """
                    ##  append the info from above to the to_csv list
                    to_csv_3[0].append(point_id)
                    to_csv_3[1].append(b_length)
                    to_csv_3[2].append(c_coord_wgs_x_land_angle_045)
                    to_csv_3[3].append(c_coord_wgs_y_land_angle_045)
                    to_csv_3[4].append(pixel_val_land_angle_045)
                    ##  also for the sea
                    to_csv_3[0].append(point_id)
                    to_csv_3[1].append(b_length * -1)
                    to_csv_3[2].append(c_coord_wgs_x_sea_angle_045)
                    to_csv_3[3].append(c_coord_wgs_y_sea_angle_045)
                    to_csv_3[4].append(pixel_val_sea_angle_045)

                    ##   write results to .csv file
                    with open(r'g:\Water_Nexus\Databases\coastal_dbase\tests\cs_points_test_045_deg.csv','wb') as h:
                        out_3 = csv.writer(h, delimiter=',',quoting=csv.QUOTE_ALL)
                        out_3.writerow(["point_ID","dist_to_cs", "X", "Y", "elev"])
                        out_3.writerows(zip(*to_csv_3))

                    ##  append the info from above to the to_csv list
                    to_csv_4[0].append(point_id)
                    to_csv_4[1].append(b_length)
                    to_csv_4[2].append(c_coord_wgs_x_land_angle_135)
                    to_csv_4[3].append(c_coord_wgs_y_land_angle_135)
                    to_csv_4[4].append(pixel_val_land_angle_135)
                    ##  also for the sea
                    to_csv_4[0].append(point_id)
                    to_csv_4[1].append(b_length * -1)
                    to_csv_4[2].append(c_coord_wgs_x_sea_angle_135)
                    to_csv_4[3].append(c_coord_wgs_y_sea_angle_135)
                    to_csv_4[4].append(pixel_val_sea_angle_135)

                    ##   write results to .csv file
                    with open(r'g:\Water_Nexus\Databases\coastal_dbase\tests\cs_points_test_135_deg.csv','wb') as h:
                        out_4 = csv.writer(h, delimiter=',',quoting=csv.QUOTE_ALL)
                        out_4.writerow(["point_ID","dist_to_cs", "X", "Y", "elev"])
                        out_4.writerows(zip(*to_csv_4))
                    """

                    ##  define sum and count of avg_values for final average value
                    sum_avg_land, cnt_avg_land = 0, 0
                    sum_avg_sea, cnt_avg_sea = 0, 0

                    ##  calculate the average value for each cross-section point
                    for b in range(n_points_avg):
                        avg_point_dist = avg_dist - b * (max_dist/n_points)

                        avg_point_coord_land = avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[1], c_coords_utm[0], avg_point_dist)
                        avg_point_coord_sea = avg_cs_point_coord(point_5km_utm[0], point_5km_utm[1], c_coords_utm[3], c_coords_utm[2], avg_point_dist)

                        ##  transform to wgs coordinates
                        #print ('avg_point_coord_land :  '), avg_point_coord_land
                        #print ('avg_point_coord_sea :  '), avg_point_coord_sea
                        avg_point_coord_wgs_land = equator_position_angles(avg_point_coord_land[1], avg_point_coord_land[0], avg_point_coord_land[3], avg_point_coord_land[2], zone_num, zone_let)
                        avg_point_coord_wgs_sea = equator_position_angles(avg_point_coord_sea[1], avg_point_coord_sea[0], avg_point_coord_sea[3], avg_point_coord_sea[2], zone_num, zone_let)


                        ##  get the pixel value for the cross-section avg point for the land pixel
                        avg_pixel_val_land_1 = read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_land[0][1], avg_point_coord_wgs_land[0][0])
                        sum_avg_land += avg_pixel_val_land_1

                        avg_pixel_val_land_2 = read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_land[1][1], avg_point_coord_wgs_land[1][0])
                        sum_avg_land += avg_pixel_val_land_2

                        cnt_avg_land += 2

                        ##  get the pixel value for the cross-section avg point for the sea pixel
                        avg_pixel_val_sea_1 = read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_sea[0][1], avg_point_coord_wgs_sea[0][0])
                        sum_avg_sea += avg_pixel_val_sea_1

                        avg_pixel_val_sea_2 = read_raster_val(rb, gt, rx, ry, avg_point_coord_wgs_sea[1][1], avg_point_coord_wgs_sea[1][0])
                        sum_avg_sea += avg_pixel_val_sea_2

                        cnt_avg_sea += 2

                        """
                        ##  append the info from above to the to_csv list
                        to_csv_2[0].append(point_id)
                        to_csv_2[1].append(b_length)
                        to_csv_2[2].append(avg_point_coord_wgs_land[0][1])
                        to_csv_2[3].append(avg_point_coord_wgs_land[0][0])
                        to_csv_2[4].append(avg_pixel_val_land_1)
                        ##  also for the sea
                        to_csv_2[0].append(point_id)
                        to_csv_2[1].append(b_length * -1)
                        to_csv_2[2].append(avg_point_coord_wgs_sea[0][1])
                        to_csv_2[3].append(avg_point_coord_wgs_sea[0][0])
                        to_csv_2[4].append(avg_pixel_val_sea_1)
                        """

                    ##  update the values in the avg database table
                    sql_update(dtbase_con[0], dtbase_cur4, 'cs_gebco_2014_avg', col_name_land, int(sum_avg_land / cnt_avg_land), 'id_cs = ' + str(point_id))
                    sql_update(dtbase_con[0], dtbase_cur4, 'cs_gebco_2014_avg', col_name_sea, int(sum_avg_sea / cnt_avg_sea), 'id_cs = ' + str(point_id))
                    #print ('cs_gebco_2014_avg', col_name_land, int(sum_avg_land / cnt_avg_land), 'id_cs = ' + str(point_id))
                    #print ('cs_gebco_2014_avg', col_name_sea, int(sum_avg_sea / cnt_avg_sea), 'id_cs = ' + str(point_id))

                    """
                    ##   write results to .csv file
                    with open(r'g:\Water_Nexus\Databases\coastal_dbase\tests\cs_points_test_NL_2.csv','wb') as h:
                        out_2 = csv.writer(h, delimiter=',',quoting=csv.QUOTE_ALL)
                        out_2.writerow(["point_ID","dist_to_cs", "X", "Y", "elev"])
                        out_2.writerows(zip(*to_csv_2))
                    """

                else:
                    pass

                """
                ##  append the info from above to the to_csv list
                to_csv[0].append(point_id)
                to_csv[1].append(b_length)
                to_csv[2].append(c_coord_wgs_x_land)
                to_csv[3].append(c_coord_wgs_y_land)
                to_csv[4].append(pixel_val_land)
                ##  also for the sea
                to_csv[0].append(point_id)
                to_csv[1].append(b_length * -1)
                to_csv[2].append(c_coord_wgs_x_sea)
                to_csv[3].append(c_coord_wgs_y_sea)
                to_csv[4].append(pixel_val_sea)
                """

            ##  ZeroDivisionError happens when point is located on top of the node (closest node)
            ##  this is the case in small islands where equidistant points are created..
            except ZeroDivisionError:
                continue
                print('Division by zero for id ' + str(point_id + 1))

        #del dtbase_cur2, dtbase_cur3, dtbase_cur4, dtbase_cur5

    t1 = time()

    print('Runtime : %f' % (t1 - t0) + '        Number of coastal points : %i' % len(points_5km))

    #pbar.finish()

"""
##   write results to .csv file
with open(r'g:\Water_Nexus\Databases\coastal_dbase\tests\cs_points_test_NL.csv','wb') as f:
    out = csv.writer(f, delimiter=',',quoting=csv.QUOTE_ALL)
    out.writerow(["point_ID","dist_to_cs", "X", "Y", "elev"])
    out.writerows(zip(*to_csv))

"""











