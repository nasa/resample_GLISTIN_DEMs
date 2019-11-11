
import os
import numpy as np
import requests
import argparse
import swath_references as sr



def createDirectoryStructure(dataFolder,yearSet):

    if 'Raw' not in os.listdir(dataFolder):
        os.mkdir(os.path.join(dataFolder,'Raw'))

    for year in yearSet:
        if str(year) not in os.listdir(os.path.join(dataFolder,'Raw')):
            os.mkdir(os.path.join(dataFolder,'Raw',str(year)))
        if 'Data' not in os.listdir(os.path.join(dataFolder,'Raw',str(year))):
            os.mkdir(os.path.join(dataFolder,'Raw',str(year),'Data'))
        if 'Metadata' not in os.listdir(os.path.join(dataFolder,'Raw',str(year))):
            os.mkdir(os.path.join(dataFolder,'Raw',str(year),'Metadata'))

def generateDownloadLinks(indexList,yearSet):
    dataDownloadLinks=[]
    metadataDownloadLinks=[]
    years=[]
    swathIDs=[]
    for index in indexList:
        for year in yearSet:
            fileID = 'OMG_GLISTIN_Elevation_'+ '{:02d}'.format(index)+'_'+str(year)
            if fileID not in sr.fileIDsMissingIn2016():
                years.append(year)
                resolution=50 #resolution is set for the resampling later, it is not used here
                swathID = sr.fileNameToSwathID(fileID,resolution)
                swathIDs.append(swathID)
                dataLink = sr.swathIDtoDataLink(swathID)
                metadatalink = sr.swathIDtoMetadataLink(swathID)
                dataDownloadLinks.append(dataLink)
                metadataDownloadLinks.append(metadatalink)
            else:
                print('Note: '+fileID[:-5]+' is not available in 2016')
    return(dataDownloadLinks,metadataDownloadLinks,years,swathIDs)

def downloadDataGrid(dataFolder,year,swathID,downloadLink):
    outputFile = os.path.join(dataFolder,'Raw',str(year),'Data',swathID+'.hgt.grd')
    with requests.get(downloadLink, stream=True) as r:
        r.raise_for_status()
        with open(outputFile, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

def downloadMetadataFile(dataFolder,year,swathID,downloadLink):
    outputFile = os.path.join(dataFolder, 'Raw', str(year), 'Metadata', swathID + '_metadata.txt')
    r = requests.get(downloadLink)
    txt = r.text
    f=open(outputFile,'w')
    f.write(txt)
    f.close()

def downloadDataSet(dataFolder,indexList,yearSet,downloadData,downloadMetadata):
    createDirectoryStructure(dataFolder,yearSet)

    dataLinks, metadataLinks, years, swathIDs = generateDownloadLinks(indexList,yearSet)

    for dd in range(len(dataLinks)):
        print('Working on swath '+swathIDs[dd])
        if downloadData:
            print('    Downloading the data grid')
            downloadDataGrid(dataFolder,years[dd],swathIDs[dd],dataLinks[dd])
        if downloadMetadata:
            print('    Downloading the metadata')
            downloadMetadataFile(dataFolder,years[dd],swathIDs[dd],metadataLinks[dd])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    #
    # Format for adding an argument:
    # parser.add_argument("-x", "--fullname", action="store", help="comment", default=default value, dest="variable_name",  type=datatype, required=T/F)
    #
    parser.add_argument("-d", "--dataFolder", action="store", help="Directory where data will be downloaded.", dest="dataFolder", type=str, required=True)

    parser.add_argument("-i", "--indexList", action="store",
                        help="List of file indices to download - a value of -1 (default) will choose all files. See data availability map to choose indices of interest.", default=-1,
                        dest="indexList", type=int, nargs='+', required=False)

    parser.add_argument("-y", "--yearSet", action="store", help="Years of file to download, a value of -1 (default) will choose all years.", default=-1, dest="yearSet", type=int, nargs='+',
                        required=False)

    parser.add_argument("-a", "--data", action="store", help="Choose 1 if you would like to download data, otherwise choose 0.", default=1,
                        dest="downloadData", type=int, required=False)

    parser.add_argument("-m", "--metadata", action="store", help="Choose 1 if you would like to download metadata, otherwise choose 0.", default=1,
                        dest="downloadMetadata", type=int,
                        required=False)

    args = parser.parse_args()

    dataFolder = args.dataFolder
    indexList = args.indexList
    yearSet = args.yearSet
    downloadData = args.downloadData
    downloadMetadata = args.downloadMetadata

    if indexList==-1 or -1 in indexList:
        indexList=np.arange(1,82)

    if yearSet==-1 or -1 in yearSet:
        yearSet=[2016,2017,2018,2019]

    if downloadData==1:
        downloadData=True
    else:
        downloadData=False

    if downloadMetadata == 1:
        downloadMetadata = True
    else:
        downloadMetadata = False

    downloadDataSet(dataFolder,indexList,yearSet,downloadData,downloadMetadata)

