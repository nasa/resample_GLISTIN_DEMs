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
import netCDF4 as nc4
import numpy as np
from scipy.interpolate import interp2d
import swath_references as ref


def read_GOCO05c_geoid(dataFolder):
    outputFolder = dataFolder
    geoidFileExists = False
    if 'Geoid' in os.listdir(outputFolder):
        if 'GOCO05C_720.nc' in os.listdir(os.path.join(outputFolder, 'Geoid')):
            geoidFileExists = True

    if geoidFileExists:
        mcFile=os.path.join(dataFolder,'Geoid','GOCO05C_720.nc')
        data = nc4.Dataset(mcFile)
        x=data['x'][:]
        y=data['y'][:]
        geoid=data['geoid'][:,:]
        geoidGeometry = [x,y,geoid]
    else:
        raise ValueError('Geoid file not found - check that file has been downloaded correctly using instructions on github.com/NASA/glistin')
    return(geoidGeometry)

def read_lat_long(dataset):
    lon = np.array(dataset['longitude'])
    lat = np.array(dataset['latitude'])
    epsg = dataset['projection'].attrs['EPSG']
    return(lon,lat,epsg)

def create_geoid_subset(geoidGeometry,lon,lat):
    geoidLon = geoidGeometry[0]
    geoidLat = geoidGeometry[1]
    gocoGeoid = geoidGeometry[2]

    minX = np.min(np.min(lon))
    maxX = np.max(np.max(lon))

    minY = np.min(np.min(lat))
    maxY = np.max(np.max(lat))

    xIndices = np.logical_and(geoidLon>minX-0.05,geoidLon<maxX+0.05)
    yIndices = np.logical_and(geoidLat > minY - 0.05, geoidLat < maxY + 0.05)

    gXsubset = geoidLon[xIndices]
    gYsubset = geoidLat[yIndices]
    geoidSubset = gocoGeoid[yIndices,:]
    geoidSubset = geoidSubset[:, xIndices]

    gXsubset = gXsubset.astype(np.float64)
    gYsubset = gYsubset.astype(np.float64)
    geoidSubset = geoidSubset.astype(np.float64)

    return(gXsubset, gYsubset, geoidSubset)

def interpolate_geoid_on_swath(gridX, gridY,set_int):
    outGrid = np.zeros_like(gridX)
    for i in range(np.shape(gridX)[0]):
        for j in range(np.shape(gridX)[1]):
            outGrid[i,j] = set_int(gridX[i,j],gridY[i,j])
    return(outGrid)

def add_geoid_layer(dataFolder, year, fileIndex, resolution,projection, geoidGeometry, printStatus=True):
    if printStatus:
        print("    Step 5: Adding geoid layer to file...")

    if 'EPSG' in projection:
        baseDirectory = 'Resampled_' + str(resolution) + 'm_' + projection.split(':')[1]
    else:
        baseDirectory = 'Resampled_' + str(resolution) + 'm'

    regridded_folder = os.path.join(dataFolder, baseDirectory, 'OMG_Ice_GLISTIN-A_L3_' + '{:02d}'.format(fileIndex))
    regridded_file = ref.indexAndYearToFileID(fileIndex, year) + '.nc'
    regridded_filepath = os.path.join(regridded_folder, regridded_file)

    print('        Reading the data geometry')
    dataset = xr.open_dataset(regridded_filepath)
    lon, lat, espg = read_lat_long(dataset)

    print('        Finding a geoid subset around the swath')
    # step 3: find geoid subset around the swath extent
    gXsubset, gYsubset, geoidSubset = create_geoid_subset(geoidGeometry, lon, lat)

    print('        Making an interpolation grid for the geoid')
    # step 4: make an interpolation function for the geoid subset
    set_int = interp2d(gXsubset, gYsubset, geoidSubset)

    print('        Interpolating the geoid onto the grid points')
    # step 5: for each point in the swath, find the value in the geoid grid
    geoid = interpolate_geoid_on_swath(lon, lat, set_int)

    print('        Geoid added to the dataset')
    dataset['geoid'] = (['y', 'x'], geoid)
    os.remove(regridded_filepath) #remove the old file so that the new one can be saved in its place
    dataset.to_netcdf(regridded_filepath)
