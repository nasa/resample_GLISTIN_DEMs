# Copyright 2020, by the California Institute of Technology.
# ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
# Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
# This software may be subject to U.S. export control laws.
# By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations.
# User has the responsibility to obtain export licenses, or other export authority as may be required before exporting
# such information to foreign countries or providing access to foreign persons.

# Codes last tested 05 April 2020 by MW and IF

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import argparse


def save_geoid_as_nc(dataFolder,geoidFile):

    f=open(os.path.join(dataFolder,'Geoid',geoidFile))
    lines=f.read()
    f.close()

    lines=lines.split('\n')

    firstLineFound = False
    for ll in range(len(lines)):
        line=lines[ll]
        if not firstLineFound:
            if 'end_of_head' in line:
                startLine=ll
                firstLineFound=True

    points=[]
    for ll in range(startLine+1,len(lines)):
        line=lines[ll].split()
        if len(line)>2:
            x=float(line[0])
            y=float(line[1])
            z=float(line[2])
            points.append([x,y,z])

    points=np.array(points)
    minX=np.min(points[:,0])
    maxX=np.max(points[:,0])
    minY=np.min(points[:,1])
    maxY=np.max(points[:,1])

    step = 0.1
    x=np.arange(minX,maxX,step)
    y=np.arange(minY,maxY+step,step)

    X,Y=np.meshgrid(x,y)
    Z=np.zeros_like(X)

    nNans=0
    for xi in range(len(x)):
        xLine = np.flipud(points[np.round(points[:,0],2)==np.round(x[xi],2),:])
        if len(xLine)>0:
            Z[:,xi] = xLine[:,2]
        else:
            nNans+=1

    x=x-360.0

    outputFile=os.path.join(dataFolder,'Geoid','GOCO05C_720.nc')
    rootgrp=nc4.Dataset(outputFile, "w", format="NETCDF4")
    rootgrp.createDimension("x", np.shape(Z)[1])
    rootgrp.createDimension("y", np.shape(Z)[0])
    xvar= rootgrp.createVariable("x", "f4", ("x",))
    yvar= rootgrp.createVariable("y", "f4", ("y",))
    geo = rootgrp.createVariable("geoid", "f4", ("y", "x"))
    xvar[:]=x
    yvar[:]=y
    geo[:, :] = Z
    rootgrp.close()




folder='/Users/mhwood/Documents/Research/Data Repository/Greenland/Geoid/'
gocoFile='GOCO05c_b6f3824cd4fa51cbc2e8eb18cca7b18f56b2ff267b4f551574404192eb222708.gdf'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()


    parser.add_argument("-d", "--dataFolder", action="store", help="Directory where data will be downloaded.", dest="dataFolder", type=str, required=True)

    parser.add_argument("-g", "--geoidFile", action="store", help="File path of GOCO05c geoid.", dest="geoidFile", type=str, required=True)


    args = parser.parse_args()

    dataFolder = args.dataFolder
    geoidFile = args.geoidFile

    save_geoid_as_nc(dataFolder,geoidFile)
