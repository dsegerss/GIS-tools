#!/usr/bin/env python
# -*- coding: latin-1 -*-
doc="""
********************************************************************************
Name: reclassRaster.py
Created: 25 Oct 2011
Author: David Segersson

Description
--------------------------------------------------------------------------------
Process rasters without reading the complete raster into memory (good for large rasters)
- Resample to a coarser resolution (using an integer factor of the original resolution)
- Reclass using a reclassification table in form of a tab-separated textfile
- Print a summary of the raster statistics
- Convert into Shape-format (fishnet), a filter can be given to produce a smaller 'sparse' raster
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
    from osgeo import ogr
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

    
def block2vector(block,layer,xll,yll,cellsizeX,cellsizeY,nodata,filter=None):
    nrows,ncols=block.shape
    #start loop at first row, first col
    yul=yll-nrows*cellsizeY #cellsizeY is negative
    
    for row in np.arange(nrows):
        for col in np.arange(ncols):
            if filter is not None:
                if block[row,col]<=filter:
                    continue
            polygon = ogr.Geometry(ogr.wkbPolygon)
            ring = ogr.Geometry(ogr.wkbLinearRing)
            
            #Keep in mind that cellsizeY is assumed to be negative
            ring.AddPoint(xll+col*cellsizeX,yul+(row+1)*cellsizeY) #lower left corner
            ring.AddPoint(xll+col*cellsizeX,yul+(row+2)*cellsizeY) #upper left corner
            ring.AddPoint(xll+col*cellsizeX+cellsizeX,yul+(row+2)*cellsizeY) #upper right corner
            ring.AddPoint(xll+col*cellsizeX+cellsizeX,yul+(row+1)*cellsizeY) #lower right corner
            ring.AddPoint(xll+col*cellsizeX,yul+(row+1)*cellsizeY) #close ring
            #ring.CloseRings()
            polygon.AddGeometry(ring)
            featureDefn = layer.GetLayerDefn()
            feature = ogr.Feature(featureDefn)
            feature.SetGeometry(polygon)
            feature.SetField('value', block[row,col])
            layer.CreateFeature(feature)
            polygon.Destroy()
            feature.Destroy()


def reclassBlock(block,classDict):
    """Reclassify values in the data block accoring to a mappings given in a dictionary
    @param block: numpy array with data
    @param classDict: dictionary with oldvalue:newvalue pairs
    """
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

def updateGridSummary(block,gridSummary,nodata):
    """
    Update the grid summary with the data from a numpy array
    @param block: numpy array with data
    @param gridSummary: dictionary with entries for sum, mean, number of negative values, number of nodata values
    @param nodata: nodata value
    """
    gridSummary["nnodata"]+=(block==nodata).sum()    
    block=np.where(block==nodata,0,block)
    gridSummary["sum"]+=block.sum()
    gridSummary["nnegative"]+=(block<0).sum()
    return gridSummary

def printGridSummary(gridSummary):
    """Print summary to stdout
    @param gridSummary: dictionary with entries for sum, mean, number of negative values, number of nodata values, nodata value
    """
    g=gridSummary
    g["mean"]=g["sum"]/float(g["nrows"]*g["ncols"]-g["nnodata"])
    print 40*"_"
    print "(xmin,xmax): (%f,%f)" %(g["xll"],g["xll"]+g["ncols"]*g["cellsizeX"])
    print "(ymin,ymax): (%f,%f)" %(g["yll"],g["yll"]+g["nrows"]*g["cellsizeY"])
    print "(cellsizeX,cellsizeY): (%f,%f)" %(g["cellsizeX"],g["cellsizeY"])
    print "(ncols,nrows): (%i,%i)" %(g["ncols"],g["nrows"])
    print "nodata value: %f" %g["nodatavalue"]
    print 40*"-"+"\nStatistics:"
    print "Sum: %f" %g["sum"]
    print "Mean: %f" %g["mean"]
    print "Number of nodata values: %i" %g["nnodata"]
    print "Number of negative values: %i" %g["nnegative"]
    print 40*"_"


def resampleBlock(block,cellFactor,method,nodata):
    """
    Resample the numpy array to a coarser grid
    @param block: a numpy array
    @param cellFactor: an integer factor to divide the cell side by
    @param method: method to use when resampling, sum, mean or majority can be used
    @param nodata: nodata value
    """
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

    parser.add_option("--summarize",
                      action="store_true",dest="summarize",
                      help="Print a summary of the input grid properties")
                      
    parser.add_option("--bandIndex",
                      action="store",dest="bandIndex",
                      help="Band index to read from",
                      default=1)

    parser.add_option("--dataType",
                      action="store",dest="dataType",
                      help="Output raster data type",
                      default="Float32")

    parser.add_option("--toShape",
                      action="store_true",dest="toShape",
                      help="Path to output shape file",
                      default=None)

    parser.add_option("--filter",
                      action="store",dest="filter",
                      help="Filter out data equal or below limit in shape output",
                      default=None)
        
    (options, args) = parser.parse_args()
    
    #------------Setting up logging capabilities -----------
    rootLogger=logger2.RootLogger(int(options.loglevel))
    logger=rootLogger.getLogger(sys.argv[0])

    #------------Process and validate options---------------
    if options.doc:
        print doc
        sys.exit()
    
    if len(args) > 0:
        parser.error("Incorrect number of arguments")

    #validate infile path
    if options.infileName is not None:
        inFilePath=path.abspath(options.infileName)
        if not path.exists(inFilePath):
            logger.error("Input raster does not exist")
            sys.exit(1)
    else:
        parser.error("No input data specified")

    #validate outfile path
    if options.outfileName is not None:
        outFilePath=path.abspath(options.outfileName)
        if options.toShape and ".shp" not in outFilePath:
            parser.error("Output shape has to to be specified with .shp extension")
            
    else:
        outFilePath=None

    #Validate filter option and convert filter to numeric value if present
    if not options.toShape and options.filter is not None:
        parser.error("Filter option only allowed together with shape output")
    elif options.toShape:
        if options.filter is not None:
            filter=float(options.filter)
        else:
            filter=None

    #read and process reclass table file
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

    #Assure that gdal is present
    if not __gdal_loaded__:
        raise OSError("Function readGDAL needs GDAL with python bindings")        

    # register all of the raster drivers
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

    #Calculate lower left corner
    xll=xul
    yll=yul+nrows*cellsizeY #cellsizeY should be a negative value

    #Rotated rasters not handled...yet
    if rot1!=0 or rot2!=0:
        print 'Rotated rasters are not supported by pyAirviro.geo.raster'
        sys.exit(1)
        
    if abs(cellsizeX)!=abs(cellsizeY):
        print 'Non-homogenous cellsizes are not supported by pyAirviro.geo.raster'
        sys.exit(1)

    bandIndex=int(options.bandIndex)

    band = ds.GetRasterBand(bandIndex)
    nodata=band.GetNoDataValue()

    #If no nodata value is present in raster, set to -9999 for completeness
    if nodata is None:
        nodata=-9999

    #Processing of data is made for blocks of the following size
    #Important - blocks are set to cover all columns for simpler processing
    #This might not always be optimal for read/write speed
    procXBlockSize=ncols 

    #process option for resampling
    if options.cellFactor is not None:
        cellFactor=int(options.cellFactor)
    else:
        cellFactor=1
        
    #Set output raster dimensions and cellsize
    procYBlockSize=cellFactor 
    newNcols=ncols/cellFactor
    newNrows=nrows/cellFactor
    newCellsizeX=cellsizeX*cellFactor
    newCellsizeY=cellsizeY*cellFactor #is a negative value as before

    #process option for dataType
    try:
        dataType=dataTypes[options.dataType]
    except KeyError:
        logger.error("Unknown datatype choose between: %s" %",".join(dataTypes.keys()))
        sys.exit(1)

    #Create and configure output raster data source
    if outFilePath is not None and not options.toShape:
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

    #Create and inititialize output vector data source
    if options.toShape:
        shapeDriver = ogr.GetDriverByName('ESRI Shapefile')
        if path.exists(outFilePath):
            shapeDriver.DeleteDataSource(outFilePath)
        shapeFile = shapeDriver.CreateDataSource(outFilePath)
        if shapeFile is None:
            logger.error("Could not open output shapefile %s" %outFilePath)
            sys.exit(1)    
        layer=shapeFile.CreateLayer(outFilePath,geom_type=ogr.wkbPolygon)
        fieldDefn = ogr.FieldDefn('value', ogr.OFTReal)
        layer.CreateField(fieldDefn)

    #inititialize input grid summary
    inputGridSummary={"sum":0,                 
                      "mean":0,
                      "nnodata":0,
                      "nnegative":0,
                      "xll":xll,
                      "yll":yll,
                      "ncols":ncols,
                      "nrows":nrows,
                      "cellsizeX":cellsizeX,
                      "cellsizeY":cellsizeY,
                      "nodatavalue":nodata}
    
    outputGridSummary={"sum":0,
                       "mean":0,
                       "nnodata":0,
                       "nnegative":0,
                       "xll":xll,
                       "yll":yll,
                       "ncols":newNcols,
                       "nrows":newNrows,
                       "cellsizeX":newCellsizeX,
                       "cellsizeY":newCellsizeY,
                       "nodatavalue":nodata}
        

    #Loop over block of raster (at least one row in each block)
    rowsOffset=0
    pg=ProgressBar(nrows,sys.stdout)
    for i in range(0, nrows, procYBlockSize):
        pg.update(i)
        data = band.ReadAsArray(0, i, procXBlockSize, procYBlockSize)

        if options.summarize:
            gridSummary=updateGridSummary(data,inputGridSummary,nodata)
     
        if options.classTable is not None:
            try:
                data=reclassBlock(data,classDict)
            except IOError as (errno,errMsg):
                logger.error(errMsg)
                sys.exit(1)

        if options.cellFactor is not None:
            try:
                data=resampleBlock(data[:,:],cellFactor,options.resamplingMethod,nodata)
            except ValueError as (errno,errMsg):
                 logger.error(errMsg)
                 sys.exit(1)
            except IOError as (errno,errMsg):
                 logger.error(errMsg)
                 sys.exit(1)

        if outFilePath is not None:
            if options.toShape:
                blockYll=yul+i*newCellsizeY #newCellsizeY is negative
                blockXll=xll
                block2vector(data,layer,blockXll,blockYll,newCellsizeX,newCellsizeY,nodata,filter)
            else:
                outBand.WriteArray(data,0,rowsOffset) #Write block to raster
                outBand.FlushCache() #Write data to disk

        if options.summarize:
            gridSummary=updateGridSummary(data,outputGridSummary,nodata)

        rowsOffset+=procYBlockSize/cellFactor #Update offset
        

    if options.toShape:
        shapeFile.Destroy()
        
    pg.finished()
    
    if options.summarize:
        print "\nInput raster summary"
        printGridSummary(inputGridSummary)
        print "\nOutput raster summary"
        printGridSummary(outputGridSummary)
       
    

if __name__=="__main__":
    main()
#       try:
#          main()
#      except:
#          import pdb, sys
#          e, m, tb = sys.exc_info()
#          pdb.post_mortem(tb)
