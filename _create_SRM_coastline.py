# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 10:54:04 2020

@author: daniel

"""


""" ---------------------------------------------------------------------------
          The script below joins the srm_id to the cut coastline polyline
    ----------------------------------------------------------------------- """

line_shp_in = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\split_lines5.shp'
pts_shp_in = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\cs_reg_id_VALID_only.shp'
line_shp_out = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\SRM_coastline_world.shp'

#line_shp_in = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\test\line.shp'
#pts_shp_in = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\test\points.shp'
#line_shp_out = r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\test\line_points.shp'

# this is the testing for a small area
line_shp = gpd.read_file(line_shp_in)  # loads road network
pt_shp = gpd.read_file(pts_shp_in)  # loads points

import fiona
lines = fiona.open(line_shp_in)
lst_lines_all = []
with fiona.open(line_shp_in) as lines:
   for line in lines:
       # print("vertices: ", line['geometry']['coordinates'])
        all_coords = line['geometry']['coordinates']
        line_x, line_y = [], []
        for coords in all_coords:
            line_x.append(coords[0])
            line_y.append(coords[1])
        lst_lines_all.append([line_x, line_y])
            
pts_lst_x, pts_lst_y, pts_lst = [], [], []

for index, row in pt_shp.iterrows():
    pt_nodes = row['geometry'].bounds
    pts_lst_x.append(round(pt_nodes[0], 6))
    pts_lst_y.append(round(pt_nodes[1], 6))
    pts_lst.append(row['cst_reg__1'])

line_lst = []

for a in range(len(lst_lines_all)):
    x_lst = [round(x, 6) for x in lst_lines_all[a][0]]
    y_lst = [round(x, 6) for x in lst_lines_all[a][1]]

    x_match = list(set(x_lst) & set(pts_lst_x))
    y_match = list(set(y_lst) & set(pts_lst_y))    
    
    #print(a, x_match, y_match)

    if len(x_match) > 1:
        srm_id_1 = pts_lst[pts_lst_x.index(x_match[0])]
        srm_id_2 = pts_lst[pts_lst_x.index(x_match[1])]
        
        #print(srm_id_1, srm_id_2)
        
        if srm_id_1 == srm_id_2:
            line_lst.append(srm_id_1)
        else:
            line_lst.append(-1)
    else:
        line_lst.append(-1)

line_shp['id_srm'] = 0
counter = 0
if len(line_lst) == line_shp.id_srm.count():
    for index, row in line_shp.iterrows():
        line_shp.loc[index, 'id_srm'] = line_lst[counter]
        counter += 1

# save the GeoDataFrame
line_shp.to_file(driver = 'ESRI Shapefile', filename = line_shp_out)




""" ---------------------------------------------------------------------------
          Extract raster values to each polygon feature (buffer around SRM)
    ----------------------------------------------------------------------- """

import fiona
import rasterio
from rasterio.mask import mask
#from shapely.geometry import mapping
import numpy as np
import pandas as pd
import shapely.wkt
import shapely.geometry
#import geopandas as gpd
#shapefile = gpd.read_file(r'g:\Water_Nexus\_A4_paper\_data_output\coast_buffer_10km.shp')

shapefile = pd.read_csv(r'g:\Water_Nexus\_A4_paper\_data_output\coast_buffer_10km.csv')
full_dataset = rasterio.open(r'g:\_ORIGINAL_DATA\Land Cover GLC-SHARE\glc_shv10_DOM.Tif')
population_dataset = rasterio.open(r'g:\_ORIGINAL_DATA\_GPW_population_2020\gpw_v4_population_count_rev11_2020_2pt5_min.asc')

gdp_dataset = rasterio.open(r'g:\_ORIGINAL_DATA\GDP_2018\GDP_PPP_2015_5arcmin.tif')

lst_out = []
# let's grab a single shapely geometry to check
for index, row in shapefile.iterrows():
    
    id_srm = row['id_srm']
    #geometry = row['geometry']
    geometry_wkt = row['WKT']
    g1 = shapely.wkt.loads(geometry_wkt)
    g2 = shapely.geometry.mapping(g1)
    
    out_image, out_transform = mask(full_dataset, [g2], crop=True)
    out_image.shape
    
    # transform to GeoJSON format
    #feature = [mapping(geometry)] # can also do this using polygon.__geo_interface__
    #out_image, out_transform = mask(full_dataset, feature, crop=True)
    #out_image.shape

    total_nonzero = np.count_nonzero(out_image)
    urban = round(100 * np.count_nonzero(out_image == 1) / total_nonzero, 1)
    agriculture = round(100 * np.count_nonzero(out_image == 2) / total_nonzero, 1)
    nature = round(100 * np.count_nonzero(out_image > 2) / total_nonzero, 1)

    print(id_srm, urban, agriculture, nature)

    out_image, out_transform = mask(population_dataset, [g2], crop=True)
    out_image.shape
    out_image[out_image < 0] = 0

    total_nonzero = np.sum(out_image)
    print(id_srm, round(total_nonzero, 0))
    
    
    out_image, out_transform = mask(gdp_dataset, [g2], crop=True)
    out_image.shape
    out_image[out_image == 0] = np.nan        
    #print(id_srm, round(np.nanmean(out_image), 0))
    gdp_ppp_avg = np.nansum(out_image)
    
    lst_out.append([id_srm, urban, agriculture, nature, round(total_nonzero, 0), round(gdp_ppp_avg, 0), round(round(gdp_ppp_avg, 0) / (round(total_nonzero, 0) * 365.25), 1)])

my_df = pd.DataFrame(lst_out)
my_df.to_csv(r'g:\Water_Nexus\_A4_paper\_data_output\_land_use_population_GDP.csv', index = False, header = ['id_srm', 'urban', 'agriculture', 'nature', 'pop_total', 'gdp_2015','gdp_2015_$/d_per_person'])



"""

import xarray as xr
gdp_dataset = xr.open_dataset(r'g:\_ORIGINAL_DATA\GDP_2018\GDP_per_capita_PPP_1990_2015_v2.nc')

gdp_dataset = xr.open_dataset(r'g:\_ORIGINAL_DATA\GDP_2018\GDP_PPP_30arcsec_v3.nc')


gdp_arr = gdp_dataset.sel(time = 2015.0)['GDP_per_capita_PPP'].values

import scipy.io
import numpy as np
import gdal
from osgeo import osr

xmin, ymin, xmax, ymax = -180., -90., 180., 90. 
nrows, ncols = np.shape(gdp_arr)
xres = (xmax - xmin)/float(ncols)
yres = (ymax - ymin)/float(nrows)
geotransform = (xmin, xres, 0, ymax, 0, -yres)  
output_raster = gdal.GetDriverByName('GTiff').Create(r'g:\_ORIGINAL_DATA\GDP_2018\GDP_PPP_2015_5arcminc.tif', gdp_arr.shape[1], gdp_arr.shape[0], 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)      
output_raster.SetProjection(srs.ExportToWkt())
output_raster.GetRasterBand(1).WriteArray(gdp_arr)
output_raster.FlushCache()
output_raster = None

"""

"""
import geopandas as gpd
from shapely.ops import split, snap

line_shp = gpd.read_file(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\_SRM_coastline_AF.shp')  # loads road network
pt_shp = gpd.read_file(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\_split_coastline_cs_pts_AF.shp')  # loads points


def split_line_by_nearest_points(gdf_line, gdf_points, tolerance):


    # union all geometries
    line = gdf_line.geometry.unary_union
    coords = gdf_points.geometry.unary_union

    # snap and split coords on line
    # returns GeometryCollection
    split_line = split(line, snap(coords, line, tolerance))
    #split_line = split(line, coords)

    # transform Geometry Collection to GeoDataFrame
    segments = [feature for feature in split_line]

    gdf_segments = gpd.GeoDataFrame(
        list(range(len(segments))), geometry=segments)
    gdf_segments.columns = ['index', 'geometry']

    return gdf_segments


line_shp = gpd.read_file(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\line_test_2.shp')  # loads road network
pt_shp = gpd.read_file(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\points_test_2.shp')  # loads points

gdf_split = split_line_by_nearest_points(line_shp, pt_shp, 0.001)
gdf_split_2 = split_line_by_nearest_points(gdf_split, pt_shp, 0.001)
gdf_split_3 = split_line_by_nearest_points(gdf_split_2, pt_shp, 0.001)

gdf_split_3.to_file(driver = 'ESRI Shapefile', filename= r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\split_test_2.shp')


import geopandas as gpd
import fiona
from shapely.ops import split
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points
from shapely.ops import snap
#   global shapefiles
line_shp = gpd.read_file(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\line_test.shp')  # loads road network
pt_shp = gpd.read_file(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\points_test.shp')  # loads points

pts_lst_x, pts_lst_y, pts_lst = [], [], []

for index, row in pt_shp.iterrows():
    pt_nodes = row['geometry'].bounds
    pts_lst_x.append(pt_nodes[0])
    pts_lst_y.append(pt_nodes[1])
    pts_lst.append(row['cst_reg__1'])

with fiona.open(r'g:\Water_Nexus\_A4_GUM\_GIS_data\_coastal_points\line_test.shp') as lines:
   for line in lines:
        print("vertices: ", line['geometry']['coordinates'])
        geom_orig = line['geometry']['coordinates']
        geom_rounded = []
        for tup in geom_orig:
            geom_rounded.append((round(tup[0], 4),round(tup[1], 4)))
        line_shp = LineString(geom_rounded)
        pts = []
        #   get all the points that intersect the line
        for i in range(len(pts_lst)):
            #if pts_lst[i] != '-1':
            point = Point(round(pts_lst_x[i], 4), round(pts_lst_y[i], 4))
            point = Point(pts_lst_x[i], pts_lst_y[i])
            
            
            if line_shp.distance(point) < 1e-04:
                print(line_shp.distance(point), pts_lst[i])
                print(line_shp.interpolate(line_shp.project(point)).intersects(line_shp))
                pts.append(list(point.coords)[0])        
        
        print(line_shp.interpolate(line_shp.project(point)).intersects(line_shp))
        
        print(line_shp.interpolate(line_shp.project(Point(pts_lst_x[i], pts_lst_y[i]))).wkt)
        
        
        line_shp = LineString(round(i, 4) for i in line['geometry']['coordinates'])
        pts = []
        #   get all the points that intersect the line
        for i in range(len(pts_lst)):
            if pts_lst[i] != '-1':
                point = Point(round(pts_lst_x[i], 4), round(pts_lst_y[i], 4))
                if line_shp.distance(point) < 1e-08:
                    print(line_shp.distance(point), pts_lst[i])
                    print(point.intersects(line_shp))
                    pts.append(list(point.coords)[0])

        # Length along line that is closest to the point
        print(line_shp.project(Point(pts[0])))
        # Now combine with interpolated point on line
        np = line_shp.interpolate(line_shp.project(Point(pts[0])))
        print(np)  # POINT (5 7)        
                
        
        print(pts)        
        
        print(Point(pts[0]).distance(line_shp))
        np = nearest_points(line_shp, Point(pts[0]))[0]
        print(np[1])
        result = snap(line_shp, Point(pts[0]), 0.0000005)
        np[0].intersection(line_shp)
        
        #   loop through that list and split the line at each point
        result = split(line_shp, Point(pts[0]))
        print(result.wkt)
        
        all_coords = line['geometry']['coordinates']
        line_x, line_y = [], []
        for coords in all_coords:
            line_x.append(coords[0])
            line_y.append(coords[1])
        lst_lines_all.append([line_x, line_y])

pt = Point((1, 1))
line = LineString([(0,0), (2,2)])
result = split(line, pt)
result.wkt

"""