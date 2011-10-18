#!/usr/bin/env python
"""
Created 2011-03-28, by David Segersson
"""

import os, sys, re
from os import path
from optparse import OptionParser

import numpy
from osgeo import gdal
from osgeo.gdalconst import *

from pyAirviro_dev.geo import raster

import pdb

#Docstrings for the option parser
usage = "usage: %prog <raster-file> [-b band]"
version="%prog 1.0"


def read(filename,bandIndex=None,dataType=numpy.float):
    """
    Read raster using GDAL drivers
    @param filename: raster filename
    @param bandIndex: index of band to read, None gives all bands
    @param dataType: convert raster to type, default=numpy.float
    """
    #ArgGIS grid: AIG
    #GeoTIFF: GTiff
    #GIF
    #BMP

    # register all of the drivers
    gdal.AllRegister()
    ds = gdal.Open(filename, GA_ReadOnly)
    if ds is None:
        print 'Could not open ' + filename
        sys.exit(1)
    
    ncols= ds.RasterXSize
    nrows = ds.RasterYSize
    nbands = ds.RasterCount
    
    #Info used for georeferencing
    geoTransform = ds.GetGeoTransform()
    xul=geoTransform[0] #top left x
    cellsizeX=geoTransform[1] #w-e pixel resolution
    rot1=geoTransform[2] #rotation, 0 if image is "north up"
    yul=geoTransform[3] #top left y
    rot2=geoTransform[4] #rotation, 0 if image is "north up"
    cellsizeY=geoTransform[5] #n-s pixel resolution
    
    
    xll=xul
    yll=yul-nrows*cellsizeY

    if rot1!=0 or rot2!=0:
        print 'Rotated rasters are not supported by pyAirviro.geo.raster'
        sys.exit(1)
        
    if abs(cellsizeX)!=abs(cellsizeY):
        print 'Non-homogenous cellsizes are not supported by pyAirviro.geo.raster'
        sys.exit(1)
    
    #Creates a list with band indices to be read
    rasters=[]
    if bandIndex is None:
        bandRange=range(1,nbands+1)
    else:
        if isinstance(bandIndex,list):
            bandRange=bandIndex
        else:
            bandRange=[bandIndex]

    for bi in bandRange:
        band = ds.GetRasterBand(bi)
        nodata=band.GetNoDataValue()
        rast=raster.raster(Xll=xll,Yll=yll,Ncols=ncols,Nrows=nrows,Cellsize=abs(cellsizeX),Nodata=nodata)


        xBlockSize,yBlockSize=band.GetBlockSize()
        data =None
        for i in range(0, nrows, yBlockSize):
            
            if i + yBlockSize < nrows:
                numRows = yBlockSize 
            else:
                numRows=nrows-i  #last block is not complete in y direction
            colBlock=None
            for j in range(0, ncols, xBlockSize):
                if j + xBlockSize < ncols:
                    numCols = xBlockSize
                else:
                    numCols=ncols-j #last block is not complete in x direction
                dataBlock = band.ReadAsArray(j, i, numCols, numRows)
                #Read whole array
                dataBlock = dataBlock.astype(dataType)

                # Process block here if large file ,before
                # reading the next block
                
                if colBlock is None:
                    colBlock=dataBlock
                else:
                    colBlock=numpy.concatenate((colBlock,dataBlock),1)
            if data is None:
                data=colBlock
            else:
                data=numpy.concatenate((data,colBlock),0)

        rast.data=data
        rasters.append(rast)
    return rasters

    
def write(rasters,filename,driverName="GTiff",dataType=GDT_Float32,proj=None):
    """
    writes list of rasters as bands in a single GDAL multi layer raster
    @param rasters: list of rasters (instances of pyAirviro.geo.raster)
    @param driverName: format specifier for GDAL raster, default='GTiff'
    @param dataType: GDAL data type, default=GDT_Float32
    @param proj: GDAL projection definition, default=None
    """
    
    #Checks that all rasters are of equal extent and dimensions
    propList=[{'ncols':rast.ncols,'nrows':rast.nrows,'xll':rast.xll,'yll':rast.yll,'cellsize':rast.cellsize} for rast in rasters]

    for ind,rastProps in enumerate(propList):
        for key in propList[0].keys():
            if rastProps[key]!=propList[0][key]:
                print "Error: %s is differs between first raster and raster %i" %(key,ind)
                sys.exit(1)
                   
    ncols=propList[0]['ncols']
    nrows=propList[0]['nrows']
    cellsize=propList[0]['cellsize']
    xll=propList[0]['xll']
    yll=propList[0]['yll']
    nbands=len(rasters)
    
    driver = gdal.GetDriverByName(driverName)
    outDataset = driver.Create(filename, ncols,nrows, nbands, dataType)
    if outDataset is None:
        print "Error: could not create output raster"
        sys.exit(1)
    
    geotransform = [xll, cellsize, 0, yll+nrows*cellsize, 0, cellsize]
    outDataset.SetGeoTransform(geotransform)
    
    for rast in rasters:    
        outBand = outDataset.GetRasterBand(1)
        outBand.WriteArray(rast.data, 0, 0)
        outBand.SetNoDataValue(rast.nodata)
        outBand.FlushCache()

    if proj !=None:
        outDataset.SetProjection(proj)
    
    
def main():
    parser=OptionParser(usage= usage, version=version)
        
    parser.add_option("-b", "--iband",
                      action="store", dest="iband", default=None,
                      help="Bands to read, specify list e.g. '0 1 2', one index used for each specified raster, default: read all")
    
    parser.add_option("-o", "--dstfile",
                      action="store", dest="dstfile", default=None,
                      help="Destination file")

    parser.add_option("-s","--summary",
                      action="store_true",dest="summarize",default=False,
                      help="Summarize rasters and write report to stdout")

    parser.add_option("-e","--evaluate",
                      action="store",dest="expression",default=None,
                      help="Expression using variables r1.0 -rn.n to represent the rasters 0-n in the argument list with band 0-n from")

    parser.add_option("-f","--format",
                      action="store",dest="format",default="GTiff",
                      help="Output raster format identifier, e.g. 'GTiff','AIG','GIF', default=%default")
                      

    (options, args) = parser.parse_args()


    if len(args)<1:
        parser.error("Wrong number of arguments")
        sys.exit(1)
        
    srcfiles=[]        
    for arg in args:
        srcfiles.append(path.abspath(args[0]))

    if options.dstfile is not None:
        dstfile=path.abspath(options.dstfile)
        
    if options.iband is not None:
        iband=[int(i) for i in options.iband.split()]
    else:
        iband=None

    if iband is not None and len(iband)!=len(srcfiles):
        print "Error: Number of rasters given as arguments not equal to number of band indices specified"
        sys.exit(1)

    rastDict={}
    for i in range(len(srcfiles)):
        if iband is not None:
            bandList=read(srcfiles[i],bandIndex=iband[i])
            rastDict["r"+str(i)+"."+str(iband[i])]=bandList[0]
        else:
            bandList=read(srcfiles[i],bandIndex=None)
            for j,band in enumerate(bandList):
                rastDict["r"+str(i)+"."+str(j)]=bandList[0]
                
    if options.summarize:
        print "Raster"+24*' '+'Band Sum'
        for key,rast in enumerate(bandList):
            srcfile=srcfiles[int(key[1])]
            bandIndex=int(key.split(".")[1])
            print "%-30s%-5i%f" %(path.basename(srcfile),bandIndex,rast.sum())

    if options.expression is not None:
        exp=options.expression
        pattern="(r[0-9]*?\.[0-9]*)"
        expMod=re.sub(pattern,"rastDict[\"\g<1>\"]",exp)

        rast=eval(expMod)
    
        write([rast],dstfile,driverName=options.format)


    #write(rasters,filename,format="GTiff",dataType=GDT_Float32,proj=None)    


if __name__=="__main__":
    main()



# def main():

#     # register all of the drivers
#     gdal.AllRegister()

#     # open the image
#     ds = gdal.Open('aster.img', GA_ReadOnly)
#     if ds is None:
#         print 'Could not open raster'
#         sys.exit(1)

#     # get image size
#     rows = ds.RasterYSize
#     cols = ds.RasterXSize
#     bands = ds.RasterCount
#     # get georeference info
#     transform = ds.GetGeoTransform()
#     xOrigin = transform[0]
#     yOrigin = transform[3]
#     pixelWidth = transform[1]
#     pixelHeight = transform[5]
#     # create a list to store band data in
#     OS Python week 4: Reading raster data [21]
#     bandList = []
#     # read in bands and store all the data in bandList
#     for i in range(bands):
#         band = ds.GetRasterBand(i+1)
#         data = band.ReadAsArray(0, 0, cols, rows)
#         bandList.append(data)
#         # loop through the coordinates
#         for i in range(3):
#             # get x,y
#             x = xValues[i]
#             y = yValues[i]
#             # compute pixel offset
#             xOffset = int((x - xOrigin) / pixelWidth)
#             yOffset = int((y - yOrigin) / pixelHeight)
#             # create a string to print out
#             s = str(x) + ' ' + str(y) + ' ' + str(xOffset) + ' ' + str(yOffset) + ' .
#             # loop through the bands and get the pixel value
#             for j in range(bands):
#             data = bandList[j]
#             value = data[yOffset, xOffset] # math matrix notation order
#             s = s + str(value) + ' .
#             # print out the data string
#             OS Python week 4: Reading raster data [22]
#             print s
#             # figure out how long the script took to run
#             endTime = time.time()
#             print 'The script took ' + str(endTime - startTime) + ' seconds'


# if __name__=__main__:
#     main()
