#!/usr/bin/env python
# -*- coding: latin-1 -*-
doc="""
********************************************************************************
Name: processRaster.py
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
import re
from optparse import OptionParser
import pdb

import numpy as np

#pyAirviro-modules
from pyAirviro.other import logger, datatable
from  pyAirviro.other.utilities import ProgressBar
from pyAirviro.other import logger

try:
    from osgeo import osr
    from osgeo import ogr
    from osgeo import gdal
    from osgeo.gdalconst import *
    __gdal_loaded__=True
except:
    __gdal_loaded__=False
    
#Docstrings for the option parser
usage = "usage: %prog [options] "
version="%prog 1.0"

if __gdal_loaded__:
    gdalDataTypes={"Float32":GDT_Float32,
                   "Int16":GDT_Int16}
    ogrDataTypes={"Real":ogr.OFTReal,
                  "Integer":ogr.OFTInteger}


#-----------Global variables -----------
log = None
# --------------------------------------

    
def block2vector(block,layer,xll,yll,cellsizeX,cellsizeY,nodata,fieldName,filter=None):
    """Write a block to a shape-file"""
    nrows,ncols=block.shape
    #start loop at first row, first col
    yul=yll-nrows*cellsizeY #cellsizeY is negative
    
    for row in np.arange(nrows):
        for col in np.arange(ncols):
            if block[row,col] == nodata:
                continue
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
            feature.SetField(fieldName, float(block[row,col]))
            layer.CreateFeature(feature)
            polygon.Destroy()
            feature.Destroy()


def reclassBlock(block,classDict,errDict):
    """Reclassify values in the data block accoring to a mappings given in a dictionary
    @param block: numpy array with data
    @param classDict: dictionary with oldvalue:newvalue pairs
    """
    try:
        classes=np.unique(block)
    except AttributeError:
        classes=np.unique1d(block)
    outBlock=block[:,:]
    for c in classes:
        try:
            val=classDict[c]
        except KeyError:
            errDesc="No value specified for reclass of category %i" %c
            if errDesc in errDict:
                errDict[errDesc]+=1
            else:
                errDict[errDesc]=1
            continue
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
    Resample the numpy array to a coarser or finer grid
    @param block: a numpy array
    @param cellFactor: a factor to divide the cell side by, if < 1 -> refinement
    @param method: method for resampling to coarser grid (sum, mean or majority)
    @param nodata: nodata value
    """
    orgNrows=block.shape[0]
    orgNcols=block.shape[1]
    
    if cellFactor>1:
        if orgNcols%cellFactor!=0 or orgNrows%cellFactor!=0:
            raise ValueError("Raster dimensions have to be dividable"+ 
                             "with the cellFactor for resampling")
    
    newNcols=int(orgNcols/cellFactor)
    newNrows=int(orgNrows/cellFactor)
    newBlock=np.zeros((newNrows,newNcols))

    if cellFactor>1:
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
            raise IOError("Resampling method %s invalid for coarsening" %method)
    else:
        refinement=int(1/cellFactor)
        newBlock=np.zeros(
            (block.shape[0]*refinement,block.shape[1]*refinement))

        if method == "keepValue":
            for row in range(orgNrows):
                for col in range(orgNcols):
                    minRow=row*refinement
                    minCol=col*refinement
                    newBlock[minRow:minRow+refinement,
                             minCol:minCol+refinement]=block[row,col]

        elif method=="keepTotal":
            nSubcells=pow(refinement,2)
            for row in range(orgNrows):
                for col in range(orgNcols):
                    val=block[row,col]
                    minRow=row*refinement
                    minCol=col*refinement
                    if val!=nodata:
                        newBlock[minRow:minRow+refinement,
                                 minCol:minCol+refinement]=val/float(nSubcells)
                        
        else:
            raise IOError("Resampling method %s invalid for refinement" %method)

    return newBlock

def reprojectBlock(outArray,block,cellFactor,blockDef,outDef,coordTrans):
    for inRow in range(blockDef["nrows"]):
        for inCol in range(blockDef["ncols"]):            
            val=block[inRow,inCol]
            x_in=blockDef["xll"]+(inCol+0.5)*blockDef["cellsize"]
            y_in=blockDef["yul"]-(inRow+0.5)*blockDef["cellsize"]
    
            x_out,y_out,z_out = coordTrans.TransformPoint(x_in, y_in)
            
            if((x_out<outDef["xll"] or y_out<outDef["yll"]) or 
               (x_out > outDef["xur"] or y_out > outDef["yul"])):
                continue
            outCol=int((x_out-outDef["xll"])/outDef["cellsize"])
            outRow=outDef["nrows"]-int(
                np.ceil((y_out-outDef["yll"])/outDef["cellsize"]))
            if outRow == outDef["nrows"]:
                outRow -=1
            if outCol == outDef["ncols"]:
                outRol -=1            
            if val != blockDef["nodata"]:
                outArray[outRow,outCol]+=val
                

def main():
    #-----------Setting up and unsing option parser-----------------------
    parser=OptionParser(usage= usage, version=version)

    parser.add_option("-d","--doc",
                      action="store_true", dest="doc",
                      help="Prints more detailed documentation and exit")
                      
    parser.add_option("-l", "--loglevel",
                      action="store",dest="loglevel",default=2,
                      help="Sets the loglevel (0-3 where 3=full logging)")

    parser.add_option("--no-progress",
                      action="store_const",dest="progressStream",
                      const=None,default=sys.stdout,
                      help="turn off the progress bar")
    
    parser.add_option("-o","--output",
                      action="store",dest="outfileName",default=None,
                      help="Output file")

    parser.add_option("-i","--input",
                      action="store",dest="infileName",
                      help="Input raster")

    parser.add_option("--bbox",
                      action="store",dest="bbox",
                      help="Only read data within bbox, --box <\"x1,y1,x2,y2\"")

    parser.add_option("--reclassify",
                      action="store",dest="classTable",
                      help="Tab-separated table with code "+
                      "and z0 value for each landuse class")

    parser.add_option("--resample",
                      action="store",dest="cellFactor",
                      help="Resample grid by dividing cellsize with an integer factor. Factor <1 results in refinement. Reprojection uses temporary refinement")

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
                      help="Output raster/shape data type",
                      default=None)

    parser.add_option("--toShape",
                      action="store_true",dest="toShape",
                      help="output as shape file",
                      default=None)

    parser.add_option("--fieldName", metavar='FIELD',
                      action="store",dest="fieldName",
                      help="write data in shape file to FIELD, default is 'value'",
                      default=None)

    parser.add_option("--filter",
                      action="store",dest="filter",
                      help="Filter out data equal or below limit in shape output",
                      default=None)

    parser.add_option("-t","--template",
                      action="store",dest="template",
                      help="Header for output raster when reprojecting")
    
    parser.add_option("--fromProj",dest="fromProj",
                      help="Input raster proj4 definition string")
    
    parser.add_option("--toProj",dest="toProj",
                      help="Output raster proj4 definition string")

        
    (options, args) = parser.parse_args()
    
    #------------Setting up logging capabilities -----------
    rootLogger=logger.RootLogger(int(options.loglevel))
    global log
    log=rootLogger.getLogger(sys.argv[0])

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
            log.error("Input raster does not exist")
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

    #Validate fieldName option
    if options.toShape:
        if options.fieldName == "":
            parser.error("fieldName can't be an empty string")
        fieldName = options.fieldName or "value"
    elif options.fieldName is not None:
        parser.error("fieldName option only allowed together with shape output")

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
        classTable=datatable.DataTable()
        classTable.read(options.classTable,defaultType=float)
        oldValueColHeader=classTable.desc[0]["id"]
        newValueColHeader=classTable.desc[0]["id"]
        classTable.setKeys([oldValueColHeader])
        log.debug("Successfully read landuse class table")

        classDict={}
        for row in classTable.data:
<<<<<<< HEAD
            classDict[row[classTable.colIndex["code"]]]=row[classTable.colIndex["z0"]]
=======
            classDict[row[0]]=row[1]
>>>>>>> c92eacf427b26e6386473d4d708db6d83cfa3d30

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
    #Read data from a window defined by option --bbox <"x1,y1,x2,y2">
    if options.bbox is not None:
        xur=xll+ncols*cellsizeX
        try:
            x1,y1,x2,y2=map(float,options.bbox.split(","))
        except:
            log.error("Invalid value for option --bbox <\"x1,y1,x2,y2\">")
            sys.exit(1)


        #Check if totally outside raster extent
        if x2<xll or y2 < yll  or x1>xur or y1>yul:
            log.error("Trying to extract outside grid boundaries")
            sys.exit(1)
        
        #Limit bbox to raster extent
        if x1<xll:
            x1=xll
        if x2>xur:
            x2=xur
        if y1<yll:
            y1=yll
        if y2>yul:
            y2=yul

        #estimate min and max of rows and cols
        colmin=int((x1-xll)/cellsizeX)
        colmaxdec=(x2-xll)/float(cellsizeX)
        rowmin=int((yul-y2)/abs(cellsizeY))
        rowmaxdec=(yul-y1)/float(abs(cellsizeY))

        if (colmaxdec-int(colmaxdec))>0:
            colmax=int(colmaxdec)
        else:
            colmax=int(colmaxdec-1)

        if (rowmaxdec-int(rowmaxdec))>0:
            rowmax=int(rowmaxdec)
        else:
            rowmax=int(rowmaxdec-1)
                
        nrows=rowmax-rowmin+1
        ncols=colmax-colmin+1
        xll= xll+colmin*cellsizeX
        yll=yul+(rowmax+1)*cellsizeY #cellsizeY is negative
        yul=yll-nrows*cellsizeY
    else:
        rowmin=0
        colmin=0

    #process option for resampling
    if options.cellFactor is not None:
        cellFactor=float(options.cellFactor)
    else:
        cellFactor=1
    
    if cellFactor>=1:
        procYBlockSize=int(cellFactor)
    else:
        procYBlockSize=1

    if options.toProj is None:
        #Set output raster dimensions and cellsize
        newNcols=int(ncols/cellFactor)
        newNrows=int(nrows/cellFactor)
        newCellsizeX=cellsizeX*cellFactor
        newCellsizeY=cellsizeY*cellFactor #is a negative value as before
        newXll=xll
        newYll=yll
        newYul=yul
        newNodata=nodata

    else:
        #Create coordinate transform
        if options.fromProj is None:
            src_srs = proj
        else:
            src_srs = osr.SpatialReference()
            src_srs.ImportFromProj4(options.fromProj)
        tgt_srs = osr.SpatialReference()
        tgt_srs.ImportFromProj4(options.toProj)            
        coordTrans= osr.CoordinateTransformation( src_srs, tgt_srs)

        if options.template is None:
            #estimate extent from input
            newNcols=ncols
            newNrows=nrows
            newCellsizeX=cellsizeX
            newCellsizeY=cellsizeY
            newNodata=nodata
            newXll,newYll,z= coordTrans.TransformPoint(xll, yll)
            newXll=int(newXll)
            newYll=int(newYll)
            newYul=newYll-newNrows*newCellsizeY
        else:
            header=open(options.template,"r").read()
            newNcols=int(re.compile("ncols\s*([0-9]*)").search(header).group(1))
            newNrows=int(re.compile("nrows\s*([0-9]*)").search(header).group(1))
            newCellsizeX=float(
                re.compile("cellsize\s*([0-9]*)").search(header).group(1))
            newCellsizeY=-1*newCellsizeX
            newNodata=float(re.compile("NODATA_value\s*(.*?)\n").search(header).group(1))
            newXll=float(re.compile("xllcorner\s*(.*?)\n").search(header).group(1))
            newYll=float(re.compile("yllcorner\s*(.*?)\n").search(header).group(1))
            newYul=newYll-newNrows*newCellsizeY

    #Processing of data is made for blocks of the following size
    #Important - blocks are set to cover all columns for simpler processing
    #This might not always be optimal for read/write speed
    procXBlockSize=ncols 

    #process option for dataType
    if options.toShape:
        dataTypes, defaultDataType = ogrDataTypes, 'Real'
    else:
        dataTypes, defaultDataType = gdalDataTypes, 'Float32'
    try:
        dataType=dataTypes[options.dataType or defaultDataType]
    except KeyError:
        log.error("Unknown datatype choose between: %s" %",".join(dataTypes.keys()))
        sys.exit(1)

    #Create and configure output raster data source
    if not options.toShape and outFilePath is not None:
        #Creates a raster dataset with 1 band

        mem_ds=gdal.GetDriverByName('MEM').Create(outFilePath,newNcols,newNrows,1,dataType)
        
        if mem_ds is None:
            print "Error: could not create output raster"
            sys.exit(1)

        outGeotransform = [newXll, newCellsizeX, 0, newYul, 0, newCellsizeY]
        mem_ds.SetGeoTransform(outGeotransform)
        if options.toProj is not None:
            mem_ds.SetProjection(tgt_srs.ExportToWkt())
        else:
            mem_ds.SetProjection(proj)
        outBand = mem_ds.GetRasterBand(1)
        outBand.SetNoDataValue(newNodata) #Set nodata-value

    if options.fromProj!=options.toProj:
        outArray=np.zeros((newNrows,newNcols))
        outDef={"ncols":newNcols,
                "nrows":newNrows,
                "xll":newXll,
                "yll":newYll,
                "yul":newYll-newNrows*newCellsizeY,
                "xur":newXll+newNcols*newCellsizeX,
                "cellsize":newCellsizeX
                }

    #Create and inititialize output vector data source
    if options.toShape:
        shapeDriver = ogr.GetDriverByName('ESRI Shapefile')
        if path.exists(outFilePath):
            shapeDriver.DeleteDataSource(outFilePath)
        shapeFile = shapeDriver.CreateDataSource(outFilePath)
        if shapeFile is None:
            log.error("Could not open output shapefile %s" %outFilePath)
            sys.exit(1)    
        layer=shapeFile.CreateLayer(outFilePath,geom_type=ogr.wkbPolygon)
        fieldDefn = ogr.FieldDefn(fieldName, dataType)
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
                       "xll":newXll,
                       "yll":newYll,
                       "ncols":newNcols,
                       "nrows":newNrows,
                       "cellsizeX":newCellsizeX,
                       "cellsizeY":newCellsizeY,
                       "nodatavalue":newNodata}
        
    #Loop over block of raster (at least one row in each block)
    rowsOffset=0
    errDict={}
    pg=ProgressBar(nrows,options.progressStream)

    for i in range(0, nrows, int(procYBlockSize)):
        pg.update(i)
        data = band.ReadAsArray(xoff=colmin, yoff=rowmin+i,
                                win_xsize=procXBlockSize, 
                                win_ysize=procYBlockSize)

        if options.summarize:
            gridSummary=updateGridSummary(data,inputGridSummary,nodata)
     
        if options.classTable is not None:
            try:
                data=reclassBlock(data,classDict,errDict)
            except IOError as e:
                log.error(str(e))
                sys.exit(1)

        if options.cellFactor is not None:
            try:
                data=resampleBlock(data[:,:],
                                   cellFactor,
                                   options.resamplingMethod,
                                   nodata)

            except ValueError as err:
                 log.error(err.message)
                 sys.exit(1)
            except IOError as err:
                 log.error(err.message)
                 sys.exit(1)     

        if options.toProj is not None:
            blockDef = {"nrows":int(procYBlockSize/cellFactor),
                        "ncols":int(procXBlockSize/cellFactor),
                        "xll":xll+(colmin)*cellsizeX,
                        "yul":yll-(nrows-rowmin-i)*cellsizeY,
                        "cellsize":cellsizeX*cellFactor,
                        "nodata":nodata}

            reprojectBlock(outArray,
                           data,
                           cellFactor,
                           blockDef,
                           outDef,
                           coordTrans)

        if outFilePath is not None:            
            if options.toShape:
                blockYll=yul+i*newCellsizeY #newCellsizeY is negative
                blockXll=xll
                block2vector(data,layer,blockXll,blockYll,newCellsizeX,newCellsizeY,nodata,fieldName,filter)
            elif options.toProj is None:
                outBand.WriteArray(data,0,rowsOffset) #Write block to raster
                outBand.FlushCache() #Write data to disk

        if options.summarize:
            gridSummary=updateGridSummary(data,outputGridSummary,nodata)

<<<<<<< HEAD
        rowsOffset+=int(procYBlockSize/cellFactor) #Update offset

=======
        rowsOffset+=procYBlockSize/cellFactor #Update offset
                
>>>>>>> c92eacf427b26e6386473d4d708db6d83cfa3d30
    if options.toShape:
        shapeFile.Destroy()

    if not options.toShape and options.toProj is not None:
        outBand.WriteArray(outArray,0,0)
        outBand.FlushCache() #Write data to disk

    if options.outfileName is not None and not options.toShape:
        outputDriver=ds.GetDriver()
        output_ds = outputDriver.CreateCopy(options.outfileName,mem_ds,0)
        output_ds = None
        mem_ds = None

    pg.finished()
    
    if options.summarize:
        print "\nInput raster summary"
        printGridSummary(inputGridSummary)
        print "\nOutput raster summary"
        printGridSummary(outputGridSummary)

    if errDict !={}:
        print "Errors/warnings during processing:"

    for err,nerr in errDict.items():
        print "%s err in %i cells" %(err,nerr)

if __name__=="__main__":
    main()
    # try:
    #     main()
    # except:
    #     import pdb, sys
    #     e, m, tb = sys.exc_info()
    #     pdb.post_mortem(tb)
