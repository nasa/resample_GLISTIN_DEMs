
import os
import numpy as np
import argparse
import resampler as rs
import swath_references as sr



def createSwathListDictionary(indexList,yearSet):
    swathDictionary={}
    allSwathsDictionary={}
    for index in indexList:
        swathSet=[]
        allSwathsSet=[]
        for year in range(2016,2020):
            fileID = 'OMG_GLISTIN_Elevation_' + '{:02d}'.format(index) + '_' + str(year)
            if fileID not in sr.fileIDsMissingIn2016():
                swathID = sr.fileNameToSwathID(fileID, 50)
                allSwathsSet.append(swathID)
                if year in yearSet:
                    swathSet.append(swathID)
        swathDictionary[index]=swathSet
        allSwathsDictionary[index]=allSwathsSet
    return(swathDictionary,allSwathsDictionary)


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


def commonDomainForSwathsInList(dataFolder,swathList):
    min_lon = 360
    max_lon = -360
    min_lat = 90
    max_lat = -90
    for swathID in swathList:
        year = '20'+swathID.split('_')[2][:2]
        metadataFile = os.path.join(dataFolder,'Raw',year,'Metadata',swathID+'_metadata.txt')
        variables = get_variables(metadataFile)
        min_swath_lon = np.min([variables[7], variables[9], variables[11], variables[13]])
        max_swath_lon = np.max([variables[7], variables[9], variables[11], variables[13]])
        min_swath_lat = np.min([variables[6], variables[8], variables[10], variables[12]])
        max_swath_lat = np.max([variables[6], variables[8], variables[10], variables[12]])
        min_lon=np.min([min_lon,min_swath_lon])
        max_lon = np.max([max_lon, max_swath_lon])
        min_lat = np.min([min_lat, min_swath_lat])
        max_lat = np.max([max_lat, max_swath_lat])
    return([min_lon,min_lat,max_lon,max_lat])

def resamplePipeline(dataFolder,indexList=np.arange(81)+1,resolution=50,yearSet=[2016,2017,2018,2019]):
    #step 0: generate list of swaths to resample for each index in indexList
    swathDictionary,allSwathsDictionary = createSwathListDictionary(indexList,yearSet)

    if 'Resampled_'+str(resolution)+'m' not in os.listdir(dataFolder):
        os.mkdir(os.path.join(dataFolder,'Resampled_'+str(resolution)+'m'))

    for index in indexList:
        indexID = 'OMG_GLISTIN_Elevation_' + '{:02d}'.format(index)
        swathList = swathDictionary[index]
        print('Resampling '+str(len(swathList))+' swaths for '+indexID)

        #step 1: get a common domain for each swath
        commonDomainExtent = commonDomainForSwathsInList(dataFolder, allSwathsDictionary[index])
        print('    The common domain for all '+str(len(allSwathsDictionary[index]))+' swaths is:')
        print('        Lat range: min = ' + '{:0.6f}'.format(commonDomainExtent[1]) + ', max = ' + '{:0.6f}'.format(commonDomainExtent[3]))
        print('        Lon range: min = ' + '{:0.6f}'.format(commonDomainExtent[0]) + ', max = ' + '{:0.6f}'.format(commonDomainExtent[2]))

        #step 2: resample the swaths at the prescribed resolution within the common swath
        for swathID in swathList:
            rs.resampler(dataFolder, indexID, swathID, commonDomainExtent, resolution)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    #
    # Format for adding an argument:
    # parser.add_argument("-x", "--fullname", action="store", help="comment", default=default value, dest="variable_name",  type=datatype, required=T/F)
    #
    parser.add_argument("-d", "--dataFolder", action="store", help="Directory where data will be downloaded.", dest="dataFolder", type=str, required=True)

    parser.add_argument("-r", "--resolution", action="store",
                        help="Resolutions for the resampling - the default resolution is 50m",
                        default=50,
                        dest="resolution", type=int, required=False)

    parser.add_argument("-i", "--indexList", action="store",
                        help="List of file indices to download - a value of -1 (default) will choose all files. See data availability map to choose indices of interest.", default=-1,
                        dest="indexList", type=int, nargs='+', required=False)

    parser.add_argument("-y", "--yearSet", action="store", help="Years of file to download, a value of -1 (default) will choose all years.", default=-1, dest="yearSet", type=int, nargs='+',
                        required=False)


    args = parser.parse_args()

    dataFolder = args.dataFolder
    resolution = args.resolution
    indexList = args.indexList
    yearSet = args.yearSet

    if indexList==-1 or -1 in indexList:
        indexList=np.arange(1,82)

    if yearSet==-1 or -1 in yearSet:
        yearSet=[2016,2017,2018,2019]

    resamplePipeline(dataFolder,indexList=indexList,resolution=resolution,yearSet=yearSet)