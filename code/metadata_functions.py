# Copyright 2020, by the California Institute of Technology.
# ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
# Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
# This software may be subject to U.S. export control laws.
# By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations.
# User has the responsibility to obtain export licenses, or other export authority as may be required before exporting
# such information to foreign countries or providing access to foreign persons.

# Codes last tested 05 April 2020 by MW and IF

import uuid
import urllib3
import os
import xarray as xr
import numpy as np
import netCDF4 as nc4
import imp
import datetime
import swath_references as ref
from bs4 import BeautifulSoup

########################################################################################################################
# these are some tools to help generate the metadata

########################################################################################################################
#step 1: read in the variables from the regridded file

def read_regridded_swath(filePath):
    data = xr.open_dataset(filePath)

    variables = []
    variableNames = []
    coordinates = []
    coordinateNames = []

    for data_var in data.data_vars.keys():
        if np.size(np.array(data[data_var])) > 1:
            variables.append(np.array(data[data_var]))
            variableNames.append(data_var)

    for coord in data.coords.keys():
        coordinates.append(np.array(data[coord]))
        coordinateNames.append(coord)


    projection = data['projection'].attrs['EPSG']

    return(variables,variableNames,coordinates,coordinateNames,projection)


########################################################################################################################
#step 2: generate a new swath with the variable and coordinates

def generate_new_dataset(variables,variableNames,coordinates,coordinateNames,projection,addGeoid):
    data_vars = {}
    coords = {}

    for vv in range(len(variables)):
        if variableNames[vv] == 'elevation':
            elevation = variables[vv]
            elevation[np.isnan(elevation)] = nc4.default_fillvals['f8']
            data_vars['elevation'] = (['y', 'x'], elevation)

    for vv in range(len(variables)):
        if variableNames[vv] == 'count':
            count = variables[vv]
            data_vars['elevation_count'] = (['y', 'x'], count)
        if variableNames[vv] == 'standard_deviation':
            standard_deviation = variables[vv]
            standard_deviation[np.isnan(standard_deviation)] = nc4.default_fillvals['f8']
            standard_deviation[elevation == nc4.default_fillvals['f8']] = nc4.default_fillvals['f8']
            data_vars['elevation_standardDeviation'] = (['y', 'x'], standard_deviation)
        if variableNames[vv] == 'longitude':
            longitude = variables[vv]
            data_vars['longitude'] = (['y', 'x'], longitude)
        if variableNames[vv] == 'latitude':
            latitude = variables[vv]
            data_vars['latitude'] = (['y', 'x'], latitude)
        if variableNames[vv] == 'geoid':
            if addGeoid:
                geoid = variables[vv]
                data_vars['geoid'] = (['y', 'x'], geoid)

    if addGeoid:
        quality_flag = np.zeros_like(elevation)
        indices1=np.logical_and(elevation-geoid>-5,standard_deviation>20)
        quality_flag[indices1] = 1
        indices2 = np.logical_and(elevation - geoid < -5, standard_deviation < 20)
        quality_flag[indices2] = 2
        indices3 = np.logical_and(elevation - geoid < -5, standard_deviation > 20)
        quality_flag[indices3] = 3
        quality_flag[count==0] = 4
    else:
        quality_flag = np.zeros_like(elevation)
        quality_flag[standard_deviation>20]=1
        quality_flag[count==0]=2
    quality_flag = quality_flag.astype(np.int)
    data_vars['elevation_qualityFlag'] = (['y', 'x'], quality_flag)

    data_vars['projection'] = chr(0)


    for cc in range(len(coordinates)):
        if coordinateNames[cc] in ['x']:
            x = coordinates[cc]
            coords[coordinateNames[cc]] = x
        if coordinateNames[cc] in ['y']:
            y = coordinates[cc]
            coords[coordinateNames[cc]] = y

    dataset = xr.Dataset(data_vars=data_vars, coords=coords)

    return(dataset)


########################################################################################################################
#step 3: add in the global metadata values

def read_timespan_from_metadata(dataFolder,fileID):
    swathID=ref.fileNameToSwathID(fileID)
    year='20'+swathID.split('_')[2][:2]
    metadataFile=os.path.join(dataFolder,'Raw',year,'Metadata',swathID+'_metadata.txt')
    f=open(metadataFile)
    lines=f.read()
    f.close()
    lines=lines.split('\n')
    for line in lines:
        if 'Start Time of Acquisition' in line:
            date=line.split()[-3]
            yr=date.split('-')[2]
            mo = date.split('-')[1]
            if mo=='Mar':
                mo='03'
            if mo=='Apr':
                mo='04'
            dy = date.split('-')[0]
            dy = '{:02d}'.format(int(dy))
            time=line.split()[-2]
            minTime=yr+'-'+mo+'-'+dy+'T'+time+'Z'
        if 'Stop Time of Acquisition' in line:
            date = line.split()[-3]
            yr = date.split('-')[2]
            mo = date.split('-')[1]
            if mo == 'Mar':
                mo = '03'
            if mo == 'Apr':
                mo = '04'
            dy = date.split('-')[0]
            dy = '{:02d}'.format(int(dy))
            time = line.split()[-2]
            maxTime = yr + '-' + mo + '-' + dy + 'T' + time + 'Z'

    ymd = minTime.split('T')[0]
    hms = minTime.split('T')[1][:-1]
    startTime = datetime.datetime(int(ymd.split('-')[0]), int(ymd.split('-')[1]), int(ymd.split('-')[2]),
                               int(hms.split(':')[0]), int(hms.split(':')[1]), int(hms.split(':')[2]))

    ymd = maxTime.split('T')[0]
    hms = maxTime.split('T')[1][:-1]
    endTime = datetime.datetime(int(ymd.split('-')[0]), int(ymd.split('-')[1]), int(ymd.split('-')[2]),
                                  int(hms.split(':')[0]), int(hms.split(':')[1]), int(hms.split(':')[2]))

    duration = endTime - startTime
    duration_in_s = duration.total_seconds()
    days = divmod(duration_in_s, 86400)
    hours = divmod(days[1], 3600)
    minutes = divmod(hours[1], 60)
    seconds = divmod(minutes[1], 1)
    durationString = 'P0Y0M0DT' + str(int(hours[0])) + 'H' + str(int(minutes[0])) + 'M' + str(int(seconds[0])) + 'S'

    return(minTime,maxTime,durationString)

def main_attribute_dictionary(dataFolder,regridded_filepath,resolution,elevation,longitude,latitude,projection):

    minLon = float(np.min(np.min(longitude)))
    maxLon = float(np.max(np.max(longitude)))
    lonRes = np.mean(np.mean(np.diff(longitude)))

    minLat = float(np.min(np.min(latitude)))
    maxLat = float(np.max(np.max(latitude)))
    latRes = np.abs(np.mean(np.mean(np.diff(longitude, axis=0))))

    minElev = np.min(np.min(elevation[elevation < nc4.default_fillvals['f8']]))
    maxElev = np.max(np.max(elevation[elevation < nc4.default_fillvals['f8']]))

    regridded_file = os.path.split(regridded_filepath)[-1]
    fileID = regridded_file[:-3]

    minTime, maxTime, timeDuration = read_timespan_from_metadata(dataFolder, fileID)

    creationDate = os.stat(regridded_filepath).st_birthtime
    creationDate = datetime.datetime.utcfromtimestamp(creationDate)
    modifiedDate = datetime.datetime.utcnow()

    attributeDictionary = {'title': 'OMG GLISTIN-A Elevation Data of Greenland Glaciers and Coastline',
                           'summary': 'Elevation measurements of the Greenland coastline measured by the Glacier and Land Ice Surface ' \
                                      'Topography Interferometer 2016-2019. These fields have been bin-averaged at '+str(resolution)+'m and ' \
                                      'reformatted to meet the specifications of the JPL PO.DAAC',
                           'keywords': 'Ice Sheet Elevation, Ice Sheet Mass Balance, Ice Sheet Measurements, Ice Sheet Thickness, Ice Sheet Topography',
                           'keywords_vocabulary': 'NASA Global Change Master Directory (GCMD) Science Keywords',
                           'Conventions': 'CF-1.7, ACDD-1.3',
                           'id': 'OMG_Ice_GLISTIN-A_L3',
                           'uuid': str(uuid.uuid4()),
                           'naming_authority': 'gov.nasa.jpl',
                           'cdm_data_type': 'Grid',
                           'history': 'GLISTIN-A elevation grid created by bin-averaging swath data from input product '+ref.fileNameToSwathID(fileID)+'.grd at a '+str(resolution)+' meter posting, and outputting to netCDF file format.',
                           'source': 'Elevation data from the GLISTIN-A radar interferometer.',
                           'platform': 'Gulfstream III',
                           'platform_vocabulary': 'GCMD platform keywords',
                           'instrument': 'Glacier and Ice Surface Topography Interferometer',
                           'instrument_vocabulary': 'GCMD instrument keywords',
                           'processing_level': 'L3',
                           'comment': 'Grid was generated using resample_pipeline from github/nasa/glistin.' \
                                      ' These open-source scripts can be used to generate similar grids at different resolutions.',
                           'standard_name_vocabulary': 'NetCDF Climate and Forecast (CF) Metadata Convention',
                           'acknowledgement': 'This research was carried out by the Jet Propulsion Laboratory, managed by the California Institute of Technology under a contract with the National Aeronautics and Space Administration.',
                           'license': 'Public Domain',
                           'metadata_link': 'http://podaac.jpl.nasa.gov/ws/metadata/dataset/?format=iso&shortName=' + fileID,
                           'product_version': '1.0',
                           'references': 'DOI:10.5067/OMGEV-GLNA1',
                           'creator_name': 'OMG Science Team',
                           'creator_email': 'omg-science@jpl.nasa.gov',
                           'creator_url': 'http://dx.doi.org/10.5067/OMGEV-GLNA1',
                           'creator_type': 'group',
                           'creator_institution': 'NASA Jet Propulsion Laboratory (JPL)',
                           'institution': 'NASA Jet Propulsion Laboratory (JPL)',
                           'project': 'Oceans Melting Greenland (OMG)',
                           'program': 'NASA Earth Venture Suborbital-2 (EVS-2)',
                           'contributor_name': 'Ian Fenty, Michael Wood, Brian Bachman, Matthew Gonzalgo, Forrest Graham, Yunling Lou, Delwyn Moller, Ronald Muellerschoen, Joshua Willis, Yang Zheng',
                           'contributor_role': 'principal investigator,author,author,author,author,author,author,author,principal investigator,author',
                           'publisher_name': 'Physical Oceanography Distributed Active Archive Center (PO.DAAC)',
                           'publisher_email': 'podaac@podaac.jpl.nasa.gov',
                           'publisher_url': 'http://dx.doi.org/10.5067/OMGEV-GLNA1',
                           'publisher_type': 'group',
                           'publisher_institution': 'NASA Jet Propulsion Laboratory (JPL)',
                           'geospatial_bounds': 'POLYGON((' + str(minLon) + ' ' + str(minLat) + ',' + str(
                               maxLon) + ' ' + str(minLat) + ',' + str(maxLon) + ' ' + str(maxLat) + ',' + str(
                               minLon) + ' ' + str(maxLat) + '))',
                           'geospatial_bounds_crs': 'WGS1984',
                           'geospatial_bounds_vertical_crs': 'ESPG:' + projection.split(':')[-1],
                           'geospatial_lat_min': minLat,
                           'geospatial_lat_max': maxLat,
                           'geospatial_lat_units': 'degrees_north',
                           'geospatial_lat_resolution': latRes,
                           'geospatial_lon_min': minLon,
                           'geospatial_lon_max': maxLon,
                           'geospatial_lon_units': 'degrees_east',
                           'geospatial_lon_resolution': lonRes,
                           'geospatial_vertical_min': minElev,
                           'geospatial_vertical_max': maxElev,
                           'geospatial_vertical_units': 'meters',
                           'geospatial_vertical_positive': 'up',
                           'geospatial_vertical_resolution': '0.2',
                           'time_coverage_start': minTime,
                           'time_coverage_end': maxTime,
                           'time_coverage_duration': timeDuration,
                           'time_coverage_resolution': '',
                           'date_created': 'T'.join(str(creationDate).split()) + 'Z',
                           'date_modified': '',
                           'date_issued': 'PO.DAAC',
                           'date_metadata_modified': 'T'.join(str(modifiedDate).split()) + 'Z'}
    return (attributeDictionary)

def write_global_attributes(dataset,dataFolder,regridded_filepath,resolution,variables,variableNames,coordinates,coordinateNames,projection):

    elevation = variables[variableNames.index('elevation')]
    longitude = variables[variableNames.index('longitude')]
    latitude = variables[variableNames.index('latitude')]
    # longitude = coordinates[coordinateNames.index('longitude')]
    # latitude = coordinates[coordinateNames.index('latitude')]

    attributes=['title','summary','keywords','keywords_vocabulary','Conventions','id','uuid',
                'naming_authority','cdm_data_type','history','source','platform','platform_vocabulary',
                'instrument','instrument_vocabulary','processing_level','comment','standard_name_vocabulary',
                'acknowledgement','license','product_version','references',
                'creator_name','creator_email','creator_url','creator_type','creator_institution',
                'institution','project','program','contributor_name','contributor_role',
                'publisher_name','publisher_email','publisher_url','publisher_type',
                'publisher_institution','geospatial_bounds','geospatial_bounds_crs','geospatial_bounds_vertical_crs',
                'geospatial_lat_min','geospatial_lat_max','geospatial_lat_units','geospatial_lat_resolution',
                'geospatial_lon_min','geospatial_lon_max','geospatial_lon_units','geospatial_lon_resolution',
                'geospatial_vertical_min','geospatial_vertical_max','geospatial_vertical_resolution','geospatial_vertical_units','geospatial_vertical_positive',
                'time_coverage_start','time_coverage_end','time_coverage_duration',
                'date_created']

    globalAttrDict = main_attribute_dictionary(dataFolder,regridded_filepath,resolution,elevation,longitude,latitude,projection)
    for attribute in attributes:
        if attribute in globalAttrDict.keys():
            dataset.attrs[attribute] = globalAttrDict[attribute]
    return(dataset)


########################################################################################################################
#step 4: add in the global metadata values

def write_variables_attributes(dataset,variables,variableNames,coordinates,coordinateNames,projection,addGeoid):

    elevation = variables[variableNames.index('elevation')]
    standard_deviation = variables[variableNames.index('standard_deviation')]
    count = variables[variableNames.index('count')]
    longitude = variables[variableNames.index('longitude')]
    latitude = variables[variableNames.index('latitude')]
    # longitude = coordinates[coordinateNames.index('longitude')]
    # latitude = coordinates[coordinateNames.index('latitude')]

    EPSG = projection.split(':')[1]

    dataset['elevation'].attrs['long_name'] = 'Surface elevation relative to the WGS1984 ellipsoid'
    dataset['elevation'].attrs['standard_name'] = 'surface_altitude'
    dataset['elevation'].attrs['units'] = 'meters'
    dataset['elevation'].attrs['coverage_content_type'] = 'physicalMeasurement'
    dataset['elevation'].attrs['valid_min'] = np.min(elevation[elevation!=nc4.default_fillvals['f8']])
    dataset['elevation'].attrs['valid_max'] = np.max(elevation[elevation != nc4.default_fillvals['f8']])
    dataset['elevation'].attrs['coordinates'] = 'latitude longitude'
    dataset['elevation'].attrs['scale_factor'] = 1
    dataset['elevation'].attrs['add_offset'] = 0
    dataset['elevation'].attrs['_FillValue'] = nc4.default_fillvals['f8']
    dataset['elevation'].attrs['grid_mapping'] = 'projection'
    dataset['elevation'].attrs['datum'] = '+ellps=urn:ogc:def:crs:EPSG::' + EPSG
    dataset['elevation'].attrs['comment'] = 'Elevation is referenced to the WGS1984 geoid. To obtain height above mean sea level, use correction provided in the geoid variable.'

    dataset['elevation_standardDeviation'].attrs['long_name'] = 'Standard deviation of elevation measurements within averaging bin'
    dataset['elevation_standardDeviation'].attrs['standard_name'] = 'standard_deviation'
    dataset['elevation_standardDeviation'].attrs['units'] = 'meters'
    dataset['elevation_standardDeviation'].attrs['coverage_content_type'] = 'auxiliaryInformation'
    dataset['elevation_standardDeviation'].attrs['valid_min'] = np.min(standard_deviation[standard_deviation!=nc4.default_fillvals['f8']])
    dataset['elevation_standardDeviation'].attrs['valid_max'] = np.max(standard_deviation[standard_deviation != nc4.default_fillvals['f8']])
    dataset['elevation_standardDeviation'].attrs['coordinates'] = 'latitude longitude'
    dataset['elevation_standardDeviation'].attrs['scale_factor'] = 1
    dataset['elevation_standardDeviation'].attrs['add_offset'] = 0
    dataset['elevation_standardDeviation'].attrs['_FillValue'] = nc4.default_fillvals['f8']
    dataset['elevation_standardDeviation'].attrs['grid_mapping'] = 'projection'
    dataset['elevation_standardDeviation'].attrs['datum'] = '+ellps=urn:ogc:def:crs:EPSG::' + EPSG
    dataset['elevation_standardDeviation'].attrs['comment'] = ''

    dataset['elevation_count'].attrs['long_name'] = 'Number of elevation measurements within averaging bin'
    dataset['elevation_count'].attrs['standard_name'] = 'count'
    dataset['elevation_count'].attrs['units'] = 'meters'
    dataset['elevation_count'].attrs['coverage_content_type'] = 'auxiliaryInformation'
    dataset['elevation_count'].attrs['valid_min'] = np.min(count)
    dataset['elevation_count'].attrs['valid_max'] = np.max(count)
    dataset['elevation_count'].attrs['coordinates'] = 'latitude longitude'
    dataset['elevation_count'].attrs['scale_factor'] = 1
    dataset['elevation_count'].attrs['add_offset'] = 0
    dataset['elevation_count'].attrs['_FillValue'] = nc4.default_fillvals['f8']
    dataset['elevation_count'].attrs['grid_mapping'] = 'projection'
    dataset['elevation_count'].attrs['datum'] = '+ellps=urn:ogc:def:crs:EPSG::' + EPSG
    dataset['elevation_count'].attrs['comment'] = ''

    dataset['elevation_qualityFlag'].attrs['long_name'] = 'Quality of average elevation from three criteria'
    dataset['elevation_qualityFlag'].attrs['standard_name'] = 'quality_flag'
    dataset['elevation_qualityFlag'].attrs['units'] = 'meters'
    dataset['elevation_qualityFlag'].attrs['coverage_content_type'] = 'qualityInformation'
    if addGeoid:
        dataset['elevation_qualityFlag'].attrs['flag_values'] = [0,1,2,3,4]
        dataset['elevation_qualityFlag'].attrs['flag_meanings'] = 'elevation_standardDeviation<20;elevation-geoid>-5 ' \
                                                                  'elevation_standardDeviation>20;elevation-geoid>-5 ' \
                                                                  'elevation_standardDeviation<20;elevation-geoid<-5 ' \
                                                                  'elevation_standardDeviation>20;elevation-geoid<-5 ' \
                                                                  'elevation_count=0'
        dataset['elevation_qualityFlag'].attrs['_FillValue'] = 4
    else:
        dataset['elevation_qualityFlag'].attrs['flag_values'] = [0,1,2]
        dataset['elevation_qualityFlag'].attrs['flag_meanings'] = 'elevation_standardDeviation<20 elevation_standardDeviation<20 elevation_count=0'
        dataset['elevation_qualityFlag'].attrs['_FillValue'] = 2
    dataset['elevation_qualityFlag'].attrs['coordinates'] = 'latitude longitude'
    dataset['elevation_qualityFlag'].attrs['grid_mapping'] = 'projection'
    dataset['elevation_qualityFlag'].attrs['datum'] = '+ellps=urn:ogc:def:crs:EPSG::' + EPSG
    dataset['elevation_qualityFlag'].attrs['comment'] = 'Recommended to use elevation_qualityFlag=0 for robust results.'

    dataset['longitude'].attrs['long_name'] = 'longitude'
    dataset['longitude'].attrs['standard_name'] = 'longitude'
    dataset['longitude'].attrs['units'] = 'degrees_east'
    dataset['longitude'].attrs['axis'] = 'X'
    dataset['longitude'].attrs['coverage_content_type'] = 'coordinate'
    dataset['longitude'].attrs['valid_min'] = np.min(np.min(longitude))
    dataset['longitude'].attrs['valid_max'] = np.max(np.max(longitude))
    dataset['longitude'].attrs['_FillValue'] = nc4.default_fillvals['f8']

    dataset['latitude'].attrs['long_name'] = 'latitude'
    dataset['latitude'].attrs['standard_name'] = 'latitude'
    dataset['latitude'].attrs['units'] = 'degrees_north'
    dataset['latitude'].attrs['axis'] = 'Y'
    dataset['latitude'].attrs['coverage_content_type'] = 'coordinate'
    dataset['latitude'].attrs['valid_min'] = np.min(np.min(latitude))
    dataset['latitude'].attrs['valid_max'] = np.max(np.max(latitude))
    dataset['latitude'].attrs['_FillValue'] = nc4.default_fillvals['f8']


    return (dataset)


########################################################################################################################
#step 5: add in the projection attributes

def write_projection_attributes(dataset,coordinates,coordinateNames,projection,resolution):
    x = coordinates[coordinateNames.index('x')]
    y = coordinates[coordinateNames.index('y')]

    EPSG = projection.split(':')[1]

    url = "http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(int(EPSG))
    http = urllib3.PoolManager()
    response = http.request('GET', url)
    soup = BeautifulSoup(response.data,features="html.parser")
    wkt = soup.get_text()
    # spatial_ref = wkt.replace(" ", "")
    spatial_ref = " ".join(wkt.split())

    dataset['projection'].attrs['grid_boundary_top_projected_y'] = np.max(y)
    dataset['projection'].attrs['grid_boundary_bottom_projected_y'] = np.min(y)
    dataset['projection'].attrs['grid_boundary_right_projected_x'] = np.max(x)
    dataset['projection'].attrs['grid_boundary_left_projected_x'] = np.min(x)
    dataset['projection'].attrs['parent_grid_cell_row_subset_start'] = int(0.0)
    dataset['projection'].attrs['parent_grid_cell_row_subset_end'] = int(float(len(y)))
    dataset['projection'].attrs['parent_grid_cell_column_subset_start'] = int(0.0)
    dataset['projection'].attrs['parent_grid_cell_column_subset_end'] = int(float(len(x)))
    dataset['projection'].attrs['spatial_ref']= spatial_ref
    dataset['projection'].attrs['proj4text'] = '+proj=utm +zone='+EPSG[-2:]+' +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
    dataset['projection'].attrs['srid'] ='urn:ogc:def:crs:EPSG::' + EPSG
    dataset['projection'].attrs['GeoTransform'] = str(np.min(x))+' '+str(resolution)+' 0 '+str(np.max(y))+' 0 '+str(-resolution)
    # if EPSG=='3413':
    #     dataset['projection'].attrs['grid_mapping_name'] = "polar_stereographic"
    # if EPSG[:3]=='326':
    #     dataset['projection'].attrs['grid_mapping_name'] = "transverse_mercator"
    dataset['projection'].attrs['latitude_of_projection_origin'] = int(0.0)
    dataset['projection'].attrs['longitude_of_projection_origin'] = 6*int(EPSG[-2:])-183
    dataset['projection'].attrs['scaling_factor'] = 0.9996
    dataset['projection'].attrs['false_easting'] = 500000.0
    dataset['projection'].attrs['false_northing'] = 0.0
    dataset['projection'].attrs['semimajor_radius'] = 6378137
    dataset['projection'].attrs['semiminor_radius'] = 6356752.3142
    dataset['projection'].attrs['units'] = "meters"
    return(dataset)


########################################################################################################################
#step 6: add in the coordinate attributes

def write_coordinates_attributes(dataset, coordinates,coordinateNames):

    x = coordinates[coordinateNames.index('x')]
    y = coordinates[coordinateNames.index('y')]

    dataset['x'].attrs['long_name'] = 'Cartesian x-coordinate'
    dataset['x'].attrs['standard_name'] = 'projection_x_coordinate'
    dataset['x'].attrs['units'] = 'meters'
    dataset['x'].attrs['axis'] = 'X'
    dataset['x'].attrs['coverage_content_type'] = 'coordinate'
    dataset['x'].attrs['valid_min'] = np.min(x)
    dataset['x'].attrs['valid_max'] = np.max(x)
    dataset['x'].attrs['_FillValue'] = nc4.default_fillvals['f8']
    dataset['x'].attrs['comment'] = 'Projected horizontal coordinates of the grid'

    dataset['y'].attrs['long_name'] = 'Cartesian y-coordinate'
    dataset['y'].attrs['standard_name'] = 'projection_y_coordinate'
    dataset['y'].attrs['units'] = 'meters'
    dataset['y'].attrs['axis'] = 'Y'
    dataset['y'].attrs['coverage_content_type'] = 'coordinate'
    dataset['y'].attrs['valid_min'] = np.min(y)
    dataset['y'].attrs['valid_max'] = np.max(y)
    dataset['y'].attrs['_FillValue'] = nc4.default_fillvals['f8']
    dataset['y'].attrs['comment'] = 'Projected vertical coordinates of the grid'

    return(dataset)


########################################################################################################################
#step 7: add in the geoid attributes if geoid is added

def write_geoid_attributes(dataset,variables,variableNames,projection):

    geoid = variables[variableNames.index('geoid')]
    EPSG = projection.split(':')[1]

    dataset['geoid'].attrs['long_name'] = 'GOCO05C Geoid - WGS84 Ellipsoid difference'
    dataset['geoid'].attrs['standard_name'] = 'geoid_height_above_reference_ellipsoid'
    dataset['geoid'].attrs['units'] = 'meters'
    dataset['geoid'].attrs['coverage_content_type'] = 'physicalMeasurement'
    dataset['geoid'].attrs['valid_min'] = np.min(geoid)
    dataset['geoid'].attrs['valid_max'] = np.max(geoid)
    dataset['geoid'].attrs['coordinates'] = 'latitude longitude'
    dataset['geoid'].attrs['scale_factor'] = 1
    dataset['geoid'].attrs['add_offset'] = 0
    dataset['geoid'].attrs['_FillValue'] = nc4.default_fillvals['f8']
    dataset['geoid'].attrs['grid_mapping'] = 'projection'
    dataset['geoid'].attrs['datum'] = '+ellps=urn:ogc:def:crs:EPSG::' + EPSG
    dataset['geoid'].attrs['comment'] = 'This field can be subtracted from the elevation field to obtain elevation referenced to mean sea level. Geoid correction is from Fecher et al 2017.'

    return(dataset)


########################################################################################################################

def add_DEM_metadata(dataFolder, year, fileIndex, resolution, projection, addGeoid, printStatus=True):
    if printStatus:
        print("    Step 6: Adding metadata to file...")

    if 'EPSG' in projection:
        baseDirectory = 'Resampled_' + str(resolution) + 'm_' + projection.split(':')[1]
    else:
        baseDirectory = 'Resampled_' + str(resolution) + 'm'

    regridded_folder = os.path.join(dataFolder, baseDirectory,'OMG_Ice_GLISTIN-A_L3_' + '{:02d}'.format(fileIndex))
    regridded_file = ref.indexAndYearToFileID(fileIndex, year)+'.nc'
    regridded_filepath = os.path.join(regridded_folder,regridded_file)

    #step 1: read in the variables from the regridded file
    variables,variableNames,coordinates,coordinateNames,projection = read_regridded_swath(regridded_filepath)

    #step 2: generate a new swath with the variable and coordinates
    dataset = generate_new_dataset(variables,variableNames,coordinates,coordinateNames,projection,addGeoid)

    #step 3: add in the global metadata values
    dataset = write_global_attributes(dataset,dataFolder,regridded_filepath,resolution,variables,variableNames,coordinates,coordinateNames,projection)

    #step 4: add in the global metadata values
    dataset = write_variables_attributes(dataset,variables,variableNames,coordinates,coordinateNames,projection,addGeoid)

    #step 5: add in the projection information
    dataset = write_projection_attributes(dataset, coordinates,coordinateNames, projection, resolution)

    #step 6: add in the coordinate attributes
    dataset = write_coordinates_attributes(dataset, coordinates, coordinateNames)

    #step 7: add in the geoid attributes if geoid is added
    if addGeoid:
        dataset = write_geoid_attributes(dataset, variables, variableNames, projection)

    os.remove(regridded_filepath)  # remove the old file so that the new one can be saved in its place
    dataset.to_netcdf(regridded_filepath)
