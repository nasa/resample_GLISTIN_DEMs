# Copyright 2020, by the California Institute of Technology.
# ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
# Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
# This software may be subject to U.S. export control laws.
# By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations.
# User has the responsibility to obtain export licenses, or other export authority as may be required before exporting
# such information to foreign countries or providing access to foreign persons.

# Codes last tested 05 April 2020 by MW and IF

import resample_functions as rf
import numpy as np
import swath_references as ref
import argparse
from pathlib import Path

def resample_DEM(dataFolder,year,swathIndex,resolution,projection,printStatus=True):
    if printStatus:
        print('\nResampling DEM with index '+str(swathIndex)+' in year '+str(year)+' at a resolution of '+str(resolution)+'m in projection '+projection)
    
    #step 1: for a given swath, find an extent which encompasses all available DEMs 2016-2019
    common_index_extent = rf.find_common_index_extent(dataFolder,swathIndex,printStatus)

    #step 2: read in the swath and create the input geometry
    swath_area_definition, masked_dem_points = rf.read_swath_and_create_geometry(dataFolder,swathIndex,year,printStatus)

    #step 3: create the output geometry for the swath
    output_area_defintion, x, y, projection = rf.create_output_geometry(common_index_extent,resolution,projection,printStatus)

    #step 4: resample the DEM at the requested resolution and reproject into the desired coordinate system
    output_grid, stdev_grid, count_grid, lon_grid, lat_grid = rf.calculate_resampled_grid(swath_area_definition, masked_dem_points, output_area_defintion, resolution, printStatus)

    #step 5: save the swath as an nc file which can be read into QGIS
    rf.save_resample(dataFolder,year,swathIndex,resolution,projection,x,y,output_grid, stdev_grid, count_grid,lon_grid,lat_grid)


def resample_GLISTIN_DEMs(dataFolder,years,swathIndices,resolution,projection, addGeoid=False):
    print('\nRunning resample_GLISTIN_DEMs with the following parameters:')
    print('    dataFolder (-d)   : '+ str(dataFolder))
    print('    resolution (-r)   : '+ str(resolution)+' m (Default: 50 m)')
    if len(swathIndices) == 81:
        print('    swathIndices (-i) : 1-81 (Default: -1=all swaths)')
    else:
        print('    swathIndices (-i) : '+ str(swathIndices)[1:-1]+' (Default: 1-81)')

    print('    years (-y)        : '+ str(years)[1:-1]+' (Default: 2016,2017,2018,2019)')
    print('    projection (-p)   : '+ projection+' (Default: \'UTM\')')
    print('    addGeoid (-g)     : '+ str(addGeoid)+' (Default: 0=False)\n')

    if dataFolder.exists():
        print('Searching for downloaded swath data from ' + str(dataFolder))
    else:
        print('Specified directory with downloaded swath data does not exist: ' + str(dataFolder))
        print('... aborting program')
        exit()

    rf.create_directory_structure(dataFolder, resolution, swathIndices, projection)

    if addGeoid:
        import geoid_functions as gf
        geoidGeometry = gf.read_GOCO05c_geoid(dataFolder)

    counter=0
    runProcess=True
    for swathIndex in swathIndices:
        for year in years:
            if year==2016:
                if swathIndex in ref.swathIndicesMissingIn2016():
                    runProcess=False
            else:
                runProcess=True
            if runProcess:
                counter+=1
                resample_DEM(dataFolder, year, swathIndex, resolution, projection)

                if addGeoid:
                    if len(geoidGeometry)>0:
                        gf.add_geoid_layer(dataFolder, year, swathIndex, resolution, projection, geoidGeometry)

                import metadata_functions as mf
                mf.add_DEM_metadata(dataFolder, year, swathIndex, resolution, projection, addGeoid)

    print('\nResample Summary:')
    print('    Resampled '+str(counter)+' DEM(s) at a resolution of '+str(resolution)+' m in the '+projection+' projection')
    print('    Files saved in ' + str(dataFolder))
    print('    Requested years: ' + str(years)[1:-1])
    print('    Requested indicies: ' + str(swathIndices)[1:-1])



if __name__ == '__main__':
    parser = argparse.ArgumentParser()


    parser.add_argument("-d", "--dataFolder", action="store", help="(Required) Directory where resampled data will be stored.", dest="dataFolder", type=str, required=True)

    parser.add_argument("-r", "--resolution", action="store",
                        help="(Optional) Resolution for the resampling. Default resolution is 50m.",
                        default=50,
                        dest="resolution", type=int, required=False)

    parser.add_argument("-i", "--swathIndices", action="store",
                        help="(Optional) List of swath indices to resample. See data availability map on github.com/NASA/glistin to choose indices of interest. "
                             "A value of -1 will choose all files. Default value is -1. ", default=-1,
                        dest="swathIndices", type=int, nargs='+', required=False)

    parser.add_argument("-y", "--years", action="store", help="(Optional) List of years to resample. A value of -1 will choose all years (2016-2019 as available). "
                                                              "Default value is -1.", default=-1, dest="years", type=int, nargs='+',
                        required=False)

    parser.add_argument("-p", "--projection", action="store",
                        help="(Optional) The projection of the output DEM. Input with an EPSG reference code as EPSG:XXXX or 'UTM'. "
                             "'UTM' is the default value and will choose the UTM zone that corresponds to the "
                             "center lat/long of the grid. This data spans UTM zones 19N to 27N.", default='UTM',
                        dest="projection", type=str,
                        required=False)

    parser.add_argument("-g", "--geoid", action="store",
                        help="(Optional) Choose 1 if you would like to add a geoid correction to the file. Must be downloaded manually prior to running this script using the "
                             "instructions provided on github.com/NASA/glistin). Default value is 0 (i.e. do not include a geoid file).", default=0,
                        dest="addGeoid", type=int,
                        required=False)


    args = parser.parse_args()

    dataFolder = args.dataFolder
    dataFolder = Path(dataFolder)
    resolution = args.resolution
    swathIndices = args.swathIndices
    years = args.years
    projection = args.projection
    addGeoid = args.addGeoid

    if swathIndices==-1 or -1 in swathIndices:
        swathIndices=np.arange(1,82)

    if years==-1 or -1 in years:
        years=[2016,2017,2018,2019]

    if addGeoid==1:
        addGeoid=True
    else:
        addGeoid=False

    resample_GLISTIN_DEMs(dataFolder,years,swathIndices,resolution,projection,addGeoid=addGeoid)
