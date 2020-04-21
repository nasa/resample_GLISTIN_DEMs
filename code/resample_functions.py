# Copyright 2020, by the California Institute of Technology.
# ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
# Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
# This software may be subject to U.S. export control laws.
# By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations.
# User has the responsibility to obtain export licenses, or other export authority as may be required before exporting
# such information to foreign countries or providing access to foreign persons.

# Codes last tested 05 April 2020 by MW and IF

import xarray as xr
import os
import numpy as np
from pyresample import kd_tree, geometry
import utm
import time
import math
from pyproj import Proj, transform
import swath_references as ref

########################################################################################################################
# These are functions used throughout all main steps

def create_directory_structure(dataFolder,resolution,fileIndices,projection):
    if 'EPSG' in projection:
        baseDirectory = 'Resampled_'+str(resolution)+'m_'+projection.split(':')[1]
    else:
        baseDirectory = 'Resampled_'+str(resolution)+'m'

    if baseDirectory not in os.listdir(dataFolder):
        os.mkdir(os.path.join(dataFolder,baseDirectory))

    for fileIndex in fileIndices:
        if 'OMG_Ice_GLISTIN-A_L3_'+'{:02d}'.format(fileIndex) not in os.listdir(os.path.join(dataFolder,baseDirectory)):
            os.mkdir(os.path.join(dataFolder,baseDirectory,'OMG_Ice_GLISTIN-A_L3_'+'{:02d}'.format(fileIndex)))


def read_metadata_dictionary(dataFolder,fileID):
    swathID = ref.fileNameToSwathID(fileID)
    # print('        Original data file: '+swathID)
    metadata_dictionary={}
    with open(os.path.join(dataFolder,'Raw',fileID.split('_')[-2][:4],'Metadata',swathID+'_metadata.txt'), "r") as f:
        lines = f.readlines()
        idx = 0
        for line in lines:
            if line.startswith("GRD Latitude Lines"):  # Start of data should start at index 62 but checking just in case
                startofdata = idx
                break
            idx += 1
        # Assuming all data is on successive lines and there are 14 data points
        for r in range(14):
            ln = lines[startofdata + r].split()
            # The value will be at the 5th index for the first 6 lines and at the 6th index for the last 8
            if r < 2:
                metadata_dictionary[' '.join(ln[:3])]=int(ln[5])
            elif r>=2 and r<6:
                metadata_dictionary[' '.join(ln[:3])] = float(ln[5])
            else:
                metadata_dictionary[' '.join(ln[:4])] = float(ln[6])

    return(metadata_dictionary)

def reproject_point(point, inputCRS, outputCRS):
    x,y = transform(inputCRS, outputCRS, point[1], point[0])
    return ([x,y])



########################################################################################################################
# step 1: for a given swath, find an extent which encompasses all available DEMs

def find_common_index_extent(dataFolder,fileIndex,printStatus,useMetadata=False):
    if useMetadata:
        if printStatus:
            print('    Step 1: Finding a common extent for all DEMs with index '+str(fileIndex))

        min_lon = 360
        max_lon = -360
        min_lat = 90
        max_lat = -90
        addFileData=True
        for year in [2016,2017,2018,2019]:
            if year==2016:
                if fileIndex in ref.fileIndicesMissingIn2016():
                    addFileData=False
            else:
                addFileData=True
            if addFileData:
                fileID = ref.indexAndYearToFileID(fileIndex, year)
                metadata_dictionary = read_metadata_dictionary(dataFolder,fileID)

                min_swath_lon = metadata_dictionary['GRD Starting Longitude']
                max_swath_lon = metadata_dictionary['GRD Starting Longitude'] + metadata_dictionary['GRD Longitude Samples'] * metadata_dictionary['GRD Longitude Spacing']
                min_swath_lat = metadata_dictionary['GRD Starting Latitude'] + metadata_dictionary['GRD Latitude Lines'] * metadata_dictionary['GRD Latitude Spacing']
                max_swath_lat = metadata_dictionary['GRD Starting Latitude']
                min_lon = np.min([min_lon, min_swath_lon])
                max_lon = np.max([max_lon, max_swath_lon])
                min_lat = np.min([min_lat, min_swath_lat])
                max_lat = np.max([max_lat, max_swath_lat])

        if printStatus:
            print('        Longitude extents -> Min: '+'{:.06f}'.format(min_lon)+'  Max: '+'{:.06f}'.format(max_lon))
            print('        Latitude extents  -> Min:  ' + '{:.06f}'.format(min_lat) + '  Max:  ' + '{:.06f}'.format(max_lat))
    else:
        if printStatus:
            print('    Step 1: Finding a common extent for all DEMs with index ' + str(fileIndex))

        saved_extent = ref.indexToCommonExtent(fileIndex)
        min_lon = saved_extent[0]
        max_lon = saved_extent[1]
        min_lat = saved_extent[2]
        max_lat = saved_extent[3]

        if printStatus:
            print('        Longitude extents -> Min: ' + '{:.06f}'.format(min_lon) + '  Max: ' + '{:.06f}'.format(
                max_lon))
            print('        Latitude extents  -> Min:  ' + '{:.06f}'.format(min_lat) + '  Max:  ' + '{:.06f}'.format(
                max_lat))

    return(min_lon,max_lon,min_lat,max_lat)


########################################################################################################################
# step 2: read in the swath and create the input geometry

def read_swath_and_create_geometry(dataFolder,fileIndex,year,printStatus):
    if printStatus:
        print('    Step 2: Reading in the binary grid and creating the original geometry')
        print('        Reading in binary data from file')

    fileID = ref.indexAndYearToFileID(fileIndex, year)
    metadata_dictionary = read_metadata_dictionary(dataFolder, fileID)

    swathID = ref.fileNameToSwathID(fileID)
    dataPath = os.path.join(dataFolder, 'Raw', str(year), 'Data', swathID + '.hgt.grd')
    g = np.fromfile(dataPath, dtype='<f4')
    grid = np.reshape(g, (metadata_dictionary['GRD Latitude Lines'], metadata_dictionary['GRD Longitude Samples']))
    if printStatus:
        print('        Preparing original geometry of the swath from metadata')
    grid = np.where(grid > grid.min(), grid, np.nan)

    min_swath_lon = metadata_dictionary['GRD Starting Longitude']
    max_swath_lon = metadata_dictionary['GRD Starting Longitude'] + metadata_dictionary['GRD Longitude Samples'] * \
                    metadata_dictionary['GRD Longitude Spacing']
    min_swath_lat = metadata_dictionary['GRD Starting Latitude'] + metadata_dictionary['GRD Latitude Lines'] * \
                    metadata_dictionary['GRD Latitude Spacing']
    max_swath_lat = metadata_dictionary['GRD Starting Latitude']

    lats = np.linspace(min_swath_lat, max_swath_lat, metadata_dictionary['GRD Latitude Lines'])
    lons = np.linspace(min_swath_lon, max_swath_lon, metadata_dictionary['GRD Longitude Samples'])

    grid = np.flipud(grid)

    if printStatus:
        # print("            The grid shape is (" + str(len(lats)) + "," + str(len(lons)) + ")")
        print("            The grid shape is (" + str(np.shape(grid)[0]) + "," + str(np.shape(grid)[1]) + ")")

    # Original Area definition in swath geometry:
    Lons, Lats = np.meshgrid(lons, lats)
    Lons = np.reshape(Lons, (np.size(Lons), ))
    Lats = np.reshape(Lats, (np.size(Lats), ))
    grid = np.reshape(grid, (np.size(grid), ))

    # Remove nans so averaging is ubiquitous
    non_nans = np.invert(np.isnan(grid))
    if printStatus:
        print('            Removed '+str(np.sum(np.isnan(grid)))+' nan points out of '+str(np.size(grid))+' grid points')
    Lons = Lons[non_nans]
    Lats = Lats[non_nans]
    grid = grid[non_nans]

    area_original = geometry.SwathDefinition(lons=Lons, lats=Lats)
    return(area_original,grid)



########################################################################################################################
# step 3: create the output geometry for the swath

def create_output_geometry(common_index_extent,resolution,projection,printStatus):
    if printStatus:
        print('    Step 3: Creating the destination geometry')

    if 'EPSG' not in projection:
        zone = utm.from_latlon(np.mean([common_index_extent[2],common_index_extent[3]]),
                               np.mean([common_index_extent[0],common_index_extent[1]]))[2]
        projection = 'EPSG:326'+str(zone)
        if printStatus:
            print('        Destination geometry: '+projection)
    else:
        if printStatus:
            print('        Destination geometry: ' + projection)

    #######

    ll_corner = reproject_point([common_index_extent[0], common_index_extent[2]], 4326, int(projection.split(':')[1]))
    lr_corner = reproject_point([common_index_extent[1], common_index_extent[2]], 4326, int(projection.split(':')[1]))
    ur_corner = reproject_point([common_index_extent[1], common_index_extent[3]], 4326, int(projection.split(':')[1]))
    ul_corner = reproject_point([common_index_extent[0], common_index_extent[3]], 4326, int(projection.split(':')[1]))

    buffer_dist = 10*resolution
    left_x = np.min([ll_corner[0],ul_corner[0],ur_corner[0],lr_corner[0]])-buffer_dist
    right_x = np.max([lr_corner[0],ur_corner[0],ll_corner[0],ul_corner[0]])+buffer_dist
    bottom_y = np.min([lr_corner[1],ll_corner[1],ur_corner[1],ul_corner[1]])-buffer_dist
    top_y = np.max([ul_corner[1],ur_corner[1],ll_corner[1],lr_corner[1]])+buffer_dist

    left_x = np.round(left_x, -3)
    right_x = np.round(right_x, -3)
    bottom_y = np.round(bottom_y, -3)
    top_y = np.round(top_y, -3)

    x = np.arange(left_x,right_x,resolution)
    y = np.arange(bottom_y,top_y,resolution)
    xstep=resolution
    ystep=resolution

    area_id_new = 'WGS84'
    description_new = projection
    proj_id_new = 'WGS84'
    proj_string_new = projection
    width_new = len(x)
    height_new = len(y)
    y = np.flip(y)

    min_x = np.min(x) - xstep / 2
    max_x = np.max(x) + xstep / 2
    min_y = np.min(y) - ystep / 2
    max_y = np.max(y) + ystep / 2

    area_extent_new = (min_x, min_y, max_x, max_y)
    area_new = geometry.AreaDefinition(area_id_new, description_new, proj_id_new, proj_string_new, width_new,
                                       height_new, area_extent_new)

    return(area_new,x,y,projection)


########################################################################################################################
# step 4: resample the DEM at the requested resolution and reproject into the desired coordinate system

def calculate_resampled_grid(area_original,grid,area_new,resolution,printStatus):
    if printStatus:
        print("    Step 4: Resampling via pyresample...")

    wf = lambda r: 1.0
    start_time = time.time()

    radiusOfInfluence=math.sqrt(2) * (resolution / 2.0)
    nNeighborsMax = int((np.pi*(radiusOfInfluence+1)**2)/9)

    output =  kd_tree.resample_custom(area_original, grid, area_new, radius_of_influence=radiusOfInfluence,
                                      fill_value=np.nan, neighbours=nNeighborsMax, weight_funcs=wf, with_uncert=True)
    if printStatus:
        print("        Resampling calculation time: %s seconds" % (time.time() - start_time))


    result = output[0]
    stddev = output[1]
    count = output[2]

    #if printStatus:

    lon, lat = area_new.get_lonlats()

    #change the result to be in steps of 20 cm - this is the vertical resolution of the data
    result *=0.5
    result = np.round(result,1)
    result *=2

    return(result, stddev, count,lon,lat)


########################################################################################################################
# step 6: save the swath as an nc file which can be read into QGIS

def save_resample(dataFolder,year,fileIndex,resolution,projection,x,y,output_grid, stdev_grid, count_grid,lon_grid,lat_grid):

    if 'EPSG:326' in projection:
        baseDirectory = 'Resampled_' + str(resolution) + 'm'
    else:
        baseDirectory = 'Resampled_' + str(resolution) + 'm_' + projection.split(':')[1]

    output_folder = os.path.join(dataFolder,baseDirectory,'OMG_Ice_GLISTIN-A_L3_'+ '{:02d}'.format(fileIndex))
    output_file = ref.indexAndYearToFileID(fileIndex, year)+'.nc'
    if output_file in os.listdir(output_folder):
        os.remove(os.path.join(output_folder,output_file))


    swath = xr.Dataset( \
        {'elevation': (['y', 'x'], output_grid),
        'standard_deviation': (['y', 'x'], stdev_grid),
        'count': (['y', 'x'], count_grid),
         'projection':'',
         'longitude': (('y', 'x'), lon_grid),
         'latitude': (('y', 'x'), lat_grid)}, \
        coords={'y': y, \
                'x': x })

    swath['projection'].attrs['EPSG'] = projection

    swath.to_netcdf(os.path.join(output_folder,output_file))