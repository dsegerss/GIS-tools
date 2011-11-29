#!/usr/bin/env python
# -*- coding: latin-1 -*-
doc="""
********************************************************************************
Name: reclassRaster.py
Created: 25 Oct 2011
Author: David Segersson

Description
--------------------------------------------------------------------------------
Reclass rasters without reading the complete raster into memory
Makes it possible to reclassify very large rasters

"""
#Standard modules
from os import path
import sys,gc
from optparse import OptionParser
import pdb

import numpy as np
#pyAirviro-modules
from pyAirviro_dev.other import logger2, dataTable
from  pyAirviro_dev.other.utilities import ProgressBar

try:
    from osgeo import gdal
    from osgeo.gdalconst import *
    __gdal_loaded__=True
except:
    __gdal_loaded__=False
    
#Docstrings for the option parser
usage = "usage: %prog [options] "
version="%prog 1.0"


dataTypes={"Float32":GDT_Float32,
           "Int16":GDT_Int16}

def reclassBlock(block,classDict):
    classes=np.unique1d(block)
    outBlock=block[:,:]
    for c in classes:
        try:
            val=classDict[c]
        except KeyError:
            raise IOError("No value specified for reclass of "+
                         " category %i" %c)
        outBlock=np.where(block==c,val,outBlock)
    return outBlock


def resampleBlock(block,cellFactor,method):
    orgNrows=block.shape[0]
    orgNcols=block.shape[1]

    if orgNcols%cellFactor!=0 or orgNrows%cellFactor!=0:
        raise ValueError("Raster dimensions have to be dividable with the cellFactor for resampling")
    
    newNcols=orgNcols/cellFactor
    newNrows=orgNrows/cellFactor
    newBlock=np.zeros((newNrows,newNcols))
    
    if method =="sum":
        for row in range(newNrows):
            for col in range(newNcols):
                cells=block[row*cellFactor:(row+1)*cellFactor,
                                col*cellFactor:(col+1)*cellFactor]
                newBlock[row,col]=cells.sum()
    elif method=="mean":
        for row in range(newNrows):
            #print row/float(rast.nrows)
            for col in range(newNcols):
                cells=block[row*cellFactor:(row+1)*cellFactor,
                                col*cellFactor:(col+1)*cellFactor]                
                newBlock[row,col]=cells.mean()
    elif method=="majority":
        for row in range(newNrows):
            for col in range(newNcols):
                cells=block[row*cellFactor:(row+1)*cellFactor,
                                col*cellFactor:(col+1)*cellFactor]                
                maxCountVal=-9e99
                maxCount=0
                for val in cells.unique():
                    count=(cells==val).sum()
                    if count>maxCount:
                        maxCount=count
                        maxCountVal=val
                newBlock[row,col]=maxCountVal
    else:
        raise IOError("Method %s does not exist" %method)
    return newBlock


def main():
    #-----------Setting up and unsing option parser-----------------------
    parser=OptionParser(usage= usage, version=version)

    parser.add_option("-d","--doc",
                      action="store_true", dest="doc",
                      help="Prints more detailed documentation and exit")
                      
    parser.add_option("-l", "--loglevel",
                      action="store",dest="loglevel",default=2,
                      help="Sets the loglevel (0-3 where 3=full logging)")
    
    parser.add_option("-o","--output",
                      action="store",dest="outfileName",default=None,
                      help="Output file")

    parser.add_option("-i","--input",
                      action="store",dest="infileName",
                      help="Input raster")

    parser.add_option("--reclassify",
                      action="store",dest="classTable",
                      help="Tab-separated table with code "+
                      "and z0 value for each landuse class")

    parser.add_option("--resample",
                      action="store",dest="cellFactor",
                      help="Resample grid to lower resolution by dividing cellsize with an integer factor")

    parser.add_option("--resamplingMethod",
                      action="store",dest="resamplingMethod",
                      help="Choose between 'mean', 'sum' or 'majority', default is %default",
                      default="mean")
                      
    parser.add_option("--bandIndex",
                      action="store",dest="bandIndex",
                      help="Band index to read from",
                      default=1)

    parser.add_option("--dataType",
                      action="store",dest="dataType",
                      help="Output raster data type",
                      default="Float32")

        
    (options, args) = parser.parse_args()
    #------------Setting up logging capabilities -----------
    rootLogger=logger2.RootLogger(int(options.loglevel))
    logger=rootLogger.getLogger(sys.argv[0])
    # ------------------------------------------------------

    if options.doc:
        print doc
        sys.exit()
    
    if len(args) > 0:
        parser.error("Incorrect number of arguments")

    if options.infileName is not None:
        inFilePath=path.abspath(options.infileName)
        if not path.exists(inFilePath):
            logger.error("Input raster does not exist")
            sys.exit(1)

    if options.outfileName is not None:
        outFilePath=path.abspath(options.outfileName)
        
    if options.classTable is not None:
        classTablePath=path.abspath(options.classTable)
        classTable=dataTable.DataTable(desc=[{"id":"code","type":int},
                                             {"id":"z0","type":float}],
                                       keys=["code"])
        classTable.read(classTablePath)
        logger.debug("Successfully read landuse class table")

        classDict={}
        for row in classTable.data:
            classDict[row[classTable.colIndex("code")]]=row[classTable.colIndex("z0")]

    if not __gdal_loaded__:
        raise OSError("Function readGDAL needs GDAL with python bindings")        

    # register all of the drivers
    gdal.AllRegister()
    ds = gdal.Open(inFilePath, GA_ReadOnly)
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

    proj = ds.GetProjection()
    
    xll=xul
    yll=yul+nrows*cellsizeY #cellsizeY should be a negative value

    if rot1!=0 or rot2!=0:
        print 'Rotated rasters are not supported by pyAirviro.geo.raster'
        sys.exit(1)
        
    if abs(cellsizeX)!=abs(cellsizeY):
        print 'Non-homogenous cellsizes are not supported by pyAirviro.geo.raster'
        sys.exit(1)

    bandIndex=int(options.bandIndex)

    band = ds.GetRasterBand(bandIndex)
    nodata=band.GetNoDataValue()

    if nodata is None:
        nodata=-9999


    #Processing of data is made for blocks of the following size
    procXBlockSize=ncols #Blocks should cover all columns for simpler processing

    
    if options.cellFactor is not None:
        cellFactor=int(options.cellFactor)
    else:
        cellFactor=1
        
    procYBlockSize=cellFactor 
    newNcols=ncols/cellFactor
    newNrows=nrows/cellFactor
    newCellsizeX=cellsizeX*cellFactor
    newCellsizeY=cellsizeY*cellFactor
    try:
        dataType=dataTypes[options.dataType]
    except KeyError:
        logger.error("Unknown datatype choose between: %s" %",".join(dataTypes.keys()))
        sys.exit(1)

    #Creates a raster dataset with 1 band
    driver=ds.GetDriver()
    outDataset = driver.Create(outFilePath, newNcols,newNrows, 1, dataType)
    if outDataset is None:
        print "Error: could not create output raster"
        sys.exit(1)

    geotransform = [xll, newCellsizeX, 0, yul, 0, newCellsizeY]
    outDataset.SetGeoTransform(geotransform)
    outDataset.SetProjection(proj)
    outBand = outDataset.GetRasterBand(1)
    outBand.SetNoDataValue(nodata) #Set nodata-value

    rowsOffset=0
    pg=ProgressBar(nrows,sys.stdout)
    for i in range(0, nrows, procYBlockSize):
        pg.update(i)
        data = band.ReadAsArray(0, i, procXBlockSize, procYBlockSize)
     
        if options.classTable is not None:
            try:
                data=reclassBlock(data,classDict)
            except IOError as (errno,errMsg):
                logger.error(errMsg)
                sys.exit(1)

        if options.cellFactor is not None:
            try:
                data=resampleBlock(data[:,:],cellFactor,options.resamplingMethod)
            except ValueError as (errno,errMsg):
                 logger.error(errMsg)
                 sys.exit(1)
            except IOError as (errno,errMsg):
                 logger.error(errMsg)
                 sys.exit(1)

        outBand.WriteArray(data,0,rowsOffset) #Write block to raster
        rowsOffset+=procYBlockSize/cellFactor #Update offset
        outBand.FlushCache() #Write data to disk
        
    pg.finished()

if __name__=="__main__":
    main()
#     try:
#         main()
#     except:
#         import pdb, sys
#         e, m, tb = sys.exc_info()
#         pdb.post_mortem(tb)
