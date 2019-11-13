import numpy as np
import os
# import glob
import matplotlib.pyplot as plt
from pyresample import kd_tree, geometry
import xarray as xr
import utm
import time
import argparse  # For command line interface
import math
import netCDF4 as nc4
'''
This is a Resampler for one year of swath data.
Take data and annotation files as arguments
Store data array from data file
Parse annotation file for number of longitude and latitude lines, starting longitude and latitude for area extent, and spacing
Use number of longitude and latitude lines to reshape the the data array
Create utm arrays from the latlons
Find minimum data value and make a new data array that fills values lower than minimum with np.nan
Use parsed starting longitude and latitude corners to define area extent of original area
'''


#
# Parses data from anonotation file
# This function assumes every annotation contains information in the same format
#
def get_variables(annotation):
    variables = []
    with open(annotation, "r") as f:
        lines = f.readlines()
        idx = 0
        for line in lines:
            if line.startswith(
                    "GRD Latitude Lines"):  # Start of data should start at index 62 but checking just in case
                startofdata = idx
                break
            idx += 1

        # Assuming all data is on successive lines and there are 14 data points
        for r in range(14):
            ln = lines[startofdata + r].split()

            # The value will be at the 5th index for the first 6 lines and at the 6th index for the last 8
            if r < 6:
                variables.append(ln[5])
            else:
                variables.append(ln[6])

    # Convert the list of "strings" into floats
    variables = list(map(float, variables))
    return variables


def resample(grid, variables, commonDomainExtent,resolution):
    # Preliminary set-up is swath specific (does not use common domain extent)

    # lat/lon_lines required to be int not float
    lat_lines = int(variables[0])
    lon_lines = int(variables[1])

    lat_start = variables[2]
    lon_start = variables[3]
    lat_space = variables[4]
    lon_space = variables[5]

    lats = np.linspace(lat_start + 0.5 * lat_space, (lat_start - 0.5 * lat_space) + (lat_space * (lat_lines + 1)),
                       lat_lines)
    lons = np.linspace(lon_start + 0.5 * lon_space, (lon_start - 0.5 * lon_space) + (lon_space * (lon_lines + 1)),
                       lon_lines)

    # print("        The grid shape is ("+str(lat_lines)+","+str(lon_lines)+")")
    # print("        The grid shape is (" + str(np.shape(grid)[0]) + "," + str(np.shape(grid)[1]) + ")")

    # Original Area definition in swath geometry:
    Lons,Lats=np.meshgrid(lons,lats)
    Lons=np.reshape(Lons,(np.size(Lons),1))
    Lats=np.reshape(Lats,(np.size(Lats),1))
    grid=np.reshape(grid,(np.size(grid),1))

    # Remove nans so averaging is ubiquitous
    non_nans = np.nonzero(~np.isnan(grid))
    Lons=Lons[non_nans]
    Lats=Lats[non_nans]
    grid=grid[non_nans]

    area_original=geometry.SwathDefinition(lons=Lons,lats=Lats)

    # print("        Grid size and orginal area size must be the same size")
    # print("        Grid size is " + str(grid.size))
    # print('        The original area size is  ' + str(area_original.size))

    min_lon = commonDomainExtent[0]
    max_lon = commonDomainExtent[2]
    min_lat = commonDomainExtent[1]
    max_lat = commonDomainExtent[3]

    print('            Preparing geometry for the resample')
    (x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers, n_xcells, n_ycells,zone) = create_utms(resolution, min_lat,
                                                                                                   min_lon, max_lat,
                                                                                                   max_lon)

    print("            Grid will be reprojected to UTM zone " + str(zone)+" N")

    # Define new area for the resampled data
    area_id_new = 'WGS84'
    description_new = 'UTM ' + str(zone) + "N"
    proj_id_new = 'WGS84'
    proj_string_new = 'EPSG:326' + str(zone)
    width_new = n_xcells
    height_new = n_ycells
    area_extent_new = (x_utm_corners[0], y_utm_corners[-1], x_utm_centers[-1], y_utm_corners[0])  # may need to revisit
    area_new = geometry.AreaDefinition(area_id_new, description_new, proj_id_new, proj_string_new, width_new,
                                       height_new, area_extent_new)

    # print('        The new AreaDefinitions shape is ' + str(area_new.shape))


    print("            Resampling via pyresample...")
    wf = lambda r: 1
    start_time = time.time()

    radiusOfInfluence=math.sqrt(2 * (resolution / 2) ** 2)
    nNeighborsMax = int((2*radiusOfInfluence)**2/(np.ceil(radiusOfInfluence/3))**2)

    result = kd_tree.resample_custom(area_original, grid, area_new,
                                     radius_of_influence=radiusOfInfluence, fill_value=np.nan, neighbours=nNeighborsMax,
                                     weight_funcs=wf)

    print("            Result calculation time: --- %s seconds ---" % (time.time() - start_time))
    print("            Result shape: " +str(result.shape))


    return (result, area_new,str(zone), x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers)


#
# Gets the max and min values in x and y directions for utm
# as well as the change in both x and y.
# Returns the information as a tuple
#
def getCommonZoneNumber(ll_lat, ll_lon, ur_lat, ur_lon):
    test1=utm.from_latlon(ll_lat, ll_lon)[2]
    test2=utm.from_latlon(ur_lat, ur_lon)[2]
    zone=int(np.mean([test1,test2]))
    return(zone)



def get_utm_range(ll_lat, ll_lon, ur_lat, ur_lon):
    #get a common zone for the center point
    zone=getCommonZoneNumber(ll_lat, ll_lon, ur_lat, ur_lon)

    x_low = utm.from_latlon(ll_lat, ll_lon,force_zone_number=zone)[0]
    x_high = utm.from_latlon(ur_lat, ur_lon,force_zone_number=zone)[0]
    dx = x_high - x_low

    y_low = utm.from_latlon(ll_lat, ll_lon,force_zone_number=zone)[1]
    y_high = utm.from_latlon(ur_lat, ur_lon,force_zone_number=zone)[1]
    dy = y_high - y_low

    return (x_low, x_high, dx, y_low, y_high, dy,zone)


#
# Increases an area in all directions by a set factor given the area's corners
# Returns the new corners as a tuple
#
def scale_dimensions(factor, x_low, x_high, y_low, y_high, dx, dy):
    x_lower = x_low - factor * dx
    x_higher = x_high + factor * dx
    y_lower = y_low - factor * dy
    y_higher = y_high + factor * dy

    return (x_lower, x_higher, y_lower, y_higher)


#
# Creates utm coordinate arrays to use for saving with the xarray Dataset
#
def create_utms(resolution, ll_lat, ll_lon, ur_lat, ur_lon):
    (x_low, x_high, dx, y_low, y_high, dy,zone) = get_utm_range(ll_lat, ll_lon, ur_lat, ur_lon)
    (x_lower, x_higher, y_lower, y_higher) = scale_dimensions(.2, x_low, x_high, y_low, y_high, dx, dy)

    # Simplifying dimensions
    real_x_low = round(x_lower, -3)
    real_y_low = round(y_lower, -3)

    Dx = x_higher - real_x_low
    Dy = y_higher - real_y_low

    n_xcells = (math.ceil(Dx / resolution))
    n_ycells = (math.ceil(Dy / resolution))

    real_x_high = real_x_low + n_xcells * resolution
    real_y_high = real_y_low + n_ycells * resolution

    # print(dx)
    # print(dy)
    # print(Dx)
    # print(Dy)
    # print(n_xcells)
    # print(n_ycells)

    n_xcorners = n_xcells + 1
    n_ycorners = n_ycells + 1

    x_utm_corners = np.linspace(real_x_low, real_x_high, n_xcorners)
    y_utm_corners = np.linspace(real_y_low, real_y_high, n_ycorners)

    x_utm_centers = np.linspace(real_x_low + .5 * resolution, \
                                real_x_high - .5 * resolution, n_xcells)
    y_utm_centers = np.linspace(real_y_low + .5 * resolution, \
                                real_y_high - .5 * resolution, n_ycells)
    # x_utms = np.linspace(real_x_low, real_x_high, 76)
    # y_utms = np.linspace(real_y_low, real_y_high, 258)

    # print(x_utm_corners.shape, y_utm_corners.shape)
    # print("utm arrays created")
    return (x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers, n_xcells, n_ycells,zone)


def zoneToSpatialRef(zone):
    spatial_ref='PROJCS["WGS 84 / UTM zone '+zone+'N"'
    spatial_ref+='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],'
    spatial_ref+='PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],'
    spatial_ref+='UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],'
    central_meridian=6*int(zone)-183
    spatial_ref+='PARAMETER["central_meridian",'+str(central_meridian)+'],'
    spatial_ref+='PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],'
    spatial_ref+='AUTHORITY["EPSG","326'+zone+'"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'
    return(spatial_ref)

#
# Saves resampled swath data in netCDF format using xarray
# Requires the resampled data, utm arrays, and the new resampled area to save the swath as an xarray Dataset
# The resolution and full path of the data file are used for
#
def save_resample(dataFolder,indexID,swathID,resampled_grid, utm_x, utm_y, new_area,utm_zone, resolution):
    output_folder = os.path.join(dataFolder,'Resampled_'+str(resolution)+'m',indexID)
    year = '20' + swathID.split('_')[2][:2]
    output_file = indexID+'_'+year+'_'+ '{:03d}'.format(resolution)+'.nc'
    if output_file in os.listdir(output_folder):
        os.remove(os.path.join(output_folder,output_file))

    new_grid_lon, new_grid_lat = new_area.get_lonlats()

    swath = xr.Dataset( \
        {'elevation': (['y', 'x'], resampled_grid),'projection':''}, \
        coords={'y': utm_y, \
                'x': utm_x, \
                'longitude': (('y', 'x'), new_grid_lon), \
                'latitude': (('y', 'x'), new_grid_lat)})

    swath['elevation'].attrs['units'] = 'meters'
    swath['elevation'].attrs['_FillValue'] = np.NaN

    swath['projection'].astype(str)

    swath['x'].attrs['units'] = 'meters'
    swath['x'].attrs['long_name'] = 'projection_grid_x_centers'
    swath['x'].attrs['standard_name'] = 'projection_x_coordinate'
    swath['x'].attrs['axis'] = 'X'
    swath['x'].attrs['valid_range'] = (np.min(utm_x), np.max(utm_x))

    swath['y'].attrs['units'] = 'meters'
    swath['y'].attrs['long_name'] = 'projection_grid_y_centers'
    swath['y'].attrs['standard_name'] = 'projection_x_coordinate'
    swath['y'].attrs['axis'] = 'Y'
    swath['y'].attrs['valid_range'] = (np.min(utm_y), np.max(utm_y))

    swath['longitude'].attrs['units'] = 'degrees_east'
    swath['longitude'].attrs['long_name'] = 'longitude'
    swath['longitude'].attrs['standard_name'] = 'longitude'
    swath['longitude'].attrs['valid_range'] = (-180, 180)

    swath['latitude'].attrs['units'] = 'degrees_north'
    swath['latitude'].attrs['long_name'] = 'latitude'
    swath['latitude'].attrs['standard_name'] = 'latitude'
    swath['latitude'].attrs['valid_range'] = (-90, 90)

    swath['projection'].attrs['grid_boundary_top_projected_y'] = np.max(utm_y)
    swath['projection'].attrs['grid_boundary_bottom_projected_y'] = np.min(utm_y)
    swath['projection'].attrs['grid_boundary_right_projected_x'] = np.max(utm_x)
    swath['projection'].attrs['grid_boundary_left_projected_x'] = np.min(utm_y)
    swath['projection'].attrs['parent_grid_cell_row_subset_start'] = 0
    swath['projection'].attrs['parent_grid_cell_row_subset_end'] = np.shape(resampled_grid)[0]
    swath['projection'].attrs['parent_grid_cell_column_subset_start'] = 0
    swath['projection'].attrs['parent_grid_cell_column_subset_end'] = np.shape(resampled_grid)[1]
    swath['projection'].attrs['spatial_ref'] = zoneToSpatialRef(utm_zone)
    swath['projection'].attrs['proj4text'] = "+proj=utm +zone=" + utm_zone + " +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    swath['projection'].attrs['srid'] = "urn:ogc:def:crs:EPSG::326" + utm_zone
    swath['projection'].attrs['GeoTransform'] = str(np.min(utm_x)) + ' ' + str(resolution) + ' 0 ' + str(np.max(utm_y)) + ' 0 ' + str(resolution)
    swath['projection'].attrs['grid_mapping_name'] = 'transverse_mercator'
    swath['projection'].attrs['latitude_of_projection_origin'] = 0
    swath['projection'].attrs['longitude_of_projection_origin'] = 6 * int(utm_zone) - 183
    swath['projection'].attrs['scaling_factor'] = 0.9996
    swath['projection'].attrs['false_easting'] = 500000
    swath['projection'].attrs['false_northing'] = 0
    swath['projection'].attrs['semimajor_radius'] = 6378273.0
    swath['projection'].attrs['semiminor_radius'] = 6356889.449
    swath['projection'].attrs['units'] = 'meters'

    swath.attrs['Conventions'] = 'CF-1.7'
    swath.attrs['title'] = 'NASA/JPL GLISTIN-A Elevation Data of Greenland Glaciers and Coastline'
    swath.attrs['references'] = 'Fenty et al 2019. GLISTIN-A Elevation. Nature: Scientific Data'
    swath.attrs['program'] = ''
    swath.attrs['cdr_variable'] = ''
    swath.attrs['software_version_id'] = ''
    swath.attrs['Metadata_Link'] = ''
    swath.attrs['product_version'] = ''
    swath.attrs['spatial_resolution'] = str(resolution) + 'm'
    swath.attrs['standard_name_vocabulary'] = ''
    swath.attrs['id'] = swathID
    swath.attrs['naming_authority'] = ''
    swath.attrs['license'] = 'No constraints on data access or use'
    swath.attrs['summary'] = 'This data set provides elevation measurements above the WGS84 geoid.'
    swath.attrs['keywords'] = 'EARTH SCIENCE > CRYOSPHERE > GLACIER, Geographic Region > Arctic'
    swath.attrs['keywords_vocabulary'] = 'NASA Global Change Master Directory (GCMD) Keywords, Version 7.0.0'
    swath.attrs['cdm_data_type'] = 'Grid'
    swath.attrs['project'] = 'NASA/JPL EVS-2 Oceans Melting Greenland'
    swath.attrs['creator_url'] = ''
    swath.attrs['creator_email'] = ''
    swath.attrs['institution'] = 'JPL > Jet Propulsion Laboratory'
    swath.attrs['processing_level'] = 'Level 2'
    swath.attrs['geospatial_lat_min'] = np.min(np.min(new_grid_lat))
    swath.attrs['geospatial_lat_max'] = np.max(np.max(new_grid_lat))
    swath.attrs['geospatial_lat_units'] = "degrees_north"
    swath.attrs['geospatial_lon_min'] = np.min(np.min(new_grid_lon))
    swath.attrs['geospatial_lon_max'] = np.max(np.max(new_grid_lon))
    swath.attrs['geospatial_lon_units'] = "degrees_east"
    swath.attrs['time_coverage_duration'] = ''
    swath.attrs['time_coverage_resolution'] = ''
    swath.attrs['contributor_name'] = 'Ian Fenty, Matthew Gonzalgo, Forrest Graham, Joshua Willis, Michael Wood'
    swath.attrs['contributor_role'] = 'author,author,author,author,PrincipalInvestigator'
    swath.attrs['acknowledgement'] = 'This project was supported by ...'
    swath.attrs['source'] = ''
    swath.attrs['history'] = ''
    swath.attrs['platform'] = ''
    swath.attrs['sensor'] = 'Glacier and Land Ice Surface Topography Interferometer'
    swath.attrs['date_created'] = ''
    swath.attrs['time_coverage_start'] = ''
    swath.attrs['time_coverage_end'] = ''
    swath.attrs['comment'] = ''

    # print(swath)

    swath.to_netcdf(os.path.join(output_folder,output_file))




def resampler(dataFolder,indexID,swathID,commonDomainExtent,resolution):
    # parser = argparse.ArgumentParser()

    #
    # Format for adding an argument:
    # parser.add_argument("-x", "--fullname", action="store", help="comment", default=default value, dest="variable_name",  type=datatype, required=T/F)
    #
    # parser.add_argument("-d", "--data", action="store", help="Complete path to data file", dest="data", type=str,
    #                     required=True)
    # parser.add_argument("-t", "--test", action="store",
    #                     help="Varible to scale down data and decrease runtime for testing", default=1, dest="test",
    #                     type=int, required=False)
    # parser.add_argument("-r", "--res", action="store", help="Resolution", default=10, dest="res", type=int,
    #                     required=False)
    # parser.add_argument("-s", "--save", action="store", help="Path to save", default="/Users/gonzalgo/Desktop/",
    #                     dest="save", type=str, required=False)
    # parser.add_argument("-x", "--suffix", action="store", help="filename suffix", default="", dest="suf", type=str,
    #                     required=False)

    # args = parser.parse_args()

    year = '20' + swathID.split('_')[2][:2]
    print('    Resampling swath ' + indexID+ '_' + year)

    print("        Reading metadata")
    metadataPath = os.path.join(dataFolder,'Raw',year,'Metadata',swathID+'_metadata.txt')
    vars = get_variables(metadataPath)

    lat_lines = int(vars[0])
    lon_lines = int(vars[1])
    print('            Original data has ' + str(lat_lines) + ' lat lines and ' + str(lon_lines) + ' lon lines')

    print('        Reading in binary data from file...')
    dataPath = os.path.join(dataFolder, 'Raw', year, 'Data', swathID + '.hgt.grd')
    g = np.fromfile(dataPath, dtype='<f4')
    grid = np.reshape(g, (lat_lines, lon_lines))
    # print('        The grid min should be -10000.0, check: '+str(grid.min()))
    print('            Setting nan values to np.nan')
    grid_nan = np.where(grid > grid.min(), grid, np.nan)
    # grid_nan=np.zeros((lat_lines,lon_lines)) #for testing

    print('        Resampling at resolution of '+str(resolution)+'m')
    if indexID not in os.listdir(os.path.join(dataFolder,'Resampled_'+str(resolution)+'m')):
        os.mkdir(os.path.join(dataFolder,'Resampled_'+str(resolution)+'m',indexID))

    (resampled_grid, new_area,utm_zone, x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers) = \
        resample(grid_nan, vars, commonDomainExtent,resolution)

    save_resample(dataFolder,indexID,swathID,resampled_grid, x_utm_centers, y_utm_centers, new_area,utm_zone, resolution)

    # print("    Resampling complete for "+name)