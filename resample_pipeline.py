
import os
import numpy as np
import resampler as rs



def createSwathListDictionary(dataFolder,indexList):
    indexFile=dataFolder+'/GLISTIN File List - By Index.csv'
    f=open(indexFile)
    lines=f.read()
    f.close()
    lines=lines.split('\n')
    lines.pop(0)
    indexToSwathList={}
    for line in lines:
        if len(line)>0:
            line=line.split(',')
            index=int(line[0])
            if index in indexList:
                swathSet=[]
                for ll in line:
                    if len(ll)>2:
                        swathSet.append(ll)
                indexToSwathList[index]=swathSet
    return(indexToSwathList)


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

def resamplePipeline(dataFolder,indexList=np.arange(81)+1,resolutionSet=[50]):
    #step 0: generate list of swaths to resample for each index in indexList
    indexToSwathList = createSwathListDictionary(dataFolder,indexList)

    for resolution in resolutionSet:
        if 'Resampled_'+str(resolution)+'m' not in os.listdir(dataFolder):
            os.mkdir(os.path.join(dataFolder,'Resampled_'+str(resolution)+'m'))

    for index in indexList:
        indexID = 'OMG_GLISTIN_Elevation_' + '{:02d}'.format(index)
        print('Resampling swaths for '+indexID)

        #step 1: get a common domain for each swath
        swathList = indexToSwathList[index]
        commonDomainExtent = commonDomainForSwathsInList(dataFolder, swathList)
        print('    The common domain for all '+str(len(swathList))+' swaths is :')
        print('        Lat range: min = ' + '{:0.6f}'.format(commonDomainExtent[1]) + ', max = ' + '{:0.6f}'.format(commonDomainExtent[3]))
        print('        Lon range: min = ' + '{:0.6f}'.format(commonDomainExtent[0]) + ', max = ' + '{:0.6f}'.format(commonDomainExtent[2]))

        #step 2: resample the swaths at the prescribed resolution within the common swath
        for swathID in swathList:
            rs.resampler(dataFolder, indexID, swathID, commonDomainExtent, resolutionSet)


dataFolder='/Volumes/mhwood/UCI/Research/Data Repository/Greenland/DEM/GLISTIN/Strips'
resamplePipeline(dataFolder,indexList=np.arange(1,5))