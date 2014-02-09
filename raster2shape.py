#!/usr/bin/env python
# -*- coding: latin-1 -*-
"""
***********************************************************************
Name: raster2shape.py
Created: 7 Sep 2013
Author: David Segersson

Description
-----------------------------------------------------------------------
Convert raster into vector grid (fishnet) or points
A filter can be given to produce a smaller 'sparse' raster
"""

#Standard modules
from os import path
import sys
from optparse import OptionParser

import numpy as np

#pyAirviro-modules
from pyAirviro.other import logger
from  pyAirviro.other.utilities import ProgressBar

try:
    from osgeo import ogr
    from osgeo import gdal
    from osgeo.gdalconst import GDT_Float32, GDT_Int16, GA_ReadOnly
    __gdal_loaded__ = True
except:
    __gdal_loaded__ = False

#Docstrings for the option parser
usage = "usage: %prog [options] "
version = "%prog 1.0"

if __gdal_loaded__:
    gdalDataTypes = {"Float32": GDT_Float32,
                     "Int16": GDT_Int16}
    ogrDataTypes = {"Real": ogr.OFTReal,
                    "Integer": ogr.OFTInteger}

#-----------Global variables -----------
log = None
# --------------------------------------


def block2poly(block, layer, xll, yll, cellsizeX, cellsizeY, nodata,
               fieldName, filter=None):
    """
    Write a block to a polygon shape-file
    @param block: raster block
    @param layer: ogr layer object to write to
    @param xll: x-coordinate of lower left block corner
    @param yll: y-coordinate of lower left block corner
    @param cellsizeX: cellsize in X-direction
    @param cellsizeY: cellsize in Y-direction
    @param nodata: nodata value
    @param fieldName: name of field to write
    @param filter: value to be excluded in output shape file
    """

    nrows, ncols = block.shape
    #start loop at first row, first col
    yul = yll - nrows * cellsizeY  # cellsizeY is negative

    for row in np.arange(nrows):
        for col in np.arange(ncols):
            if block[row, col] == nodata:
                continue
            if filter is not None:
                if block[row, col] <= filter:
                    continue
            polygon = ogr.Geometry(ogr.wkbPolygon)
            ring = ogr.Geometry(ogr.wkbLinearRing)

            #Keep in mind that cellsizeY is assumed to be negative
            ring.AddPoint(xll + col * cellsizeX,
                          yul + (row + 1) * cellsizeY)  # lower left corner
            ring.AddPoint(xll + col * cellsizeX,
                          yul + (row + 2) * cellsizeY)  # upper left corner
            ring.AddPoint(xll + col * cellsizeX + cellsizeX,
                          yul + (row + 2) * cellsizeY)  # upper right corner
            ring.AddPoint(xll + col * cellsizeX + cellsizeX,
                          yul + (row + 1) * cellsizeY)  # lower right corner
            ring.AddPoint(xll + col * cellsizeX,
                          yul + (row + 1) * cellsizeY)  # close ring
            #ring.CloseRings()
            polygon.AddGeometry(ring)
            featureDefn = layer.GetLayerDefn()
            feature = ogr.Feature(featureDefn)
            feature.SetGeometry(polygon)
            feature.SetField(fieldName, float(block[row, col]))
            layer.CreateFeature(feature)
            polygon.Destroy()
            feature.Destroy()


def block2point(block, layer, xll, yll, cellsizeX, cellsizeY, nodata,
               fieldName, filter=None):
    """
    Write a block to a polygon shape-file
    @param block: raster block
    @param layer: ogr layer object to write to
    @param xll: x-coordinate of lower left block corner
    @param yll: y-coordinate of lower left block corner
    @param cellsizeX: cellsize in X-direction
    @param cellsizeY: cellsize in Y-direction
    @param nodata: nodata value
    @param fieldName: name of field to write
    @param filter: value to be excluded in output shape file
    """

    nrows, ncols = block.shape
    #start loop at first row, first col
    yul = yll - nrows * cellsizeY  # cellsizeY is negative

    for row in np.arange(nrows):
        for col in np.arange(ncols):
            if block[row, col] == nodata:
                continue
            if filter is not None:
                if block[row, col] <= filter:
                    continue

            point = ogr.Geometry(ogr.wkbPoint)
            #Keep in mind that cellsizeY is assumed to be negative
            point.AddPoint(xll + (col + 0.5) * cellsizeX,
                           yul + (row + 0.5) * cellsizeY)
            featureDefn = layer.GetLayerDefn()
            feature = ogr.Feature(featureDefn)
            feature.SetGeometry(point)
            feature.SetField(fieldName, float(block[row, col]))
            layer.CreateFeature(feature)
            point.Destroy()
            feature.Destroy()


def main():
    #-----------Setting up and unsing option parser-----------------------
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("-d", "--doc",
                      action="store_true", dest="doc",
                      help="Prints more detailed documentation and exit")

    parser.add_option("-l", "--loglevel",
                      action="store", dest="loglevel", default=2,
                      help="Sets the loglevel (0-3 where 3=full logging)")

    parser.add_option("--no-progress",
                      action="store_const", dest="progressStream",
                      const=None, default=sys.stdout,
                      help="turn off the progress bar")

    parser.add_option("-o", "--output",
                      action="store", dest="outfileName", default=None,
                      help="Output file")

    parser.add_option("-i", "--input",
                      action="store", dest="infileName",
                      help="Input raster")

    parser.add_option("--bandIndex",
                      action="store", dest="bandIndex",
                      help="Band index to read from",
                      default=1)

    parser.add_option("--dataType",
                      action="store", dest="dataType",
                      help="Output data type",
                      default=None)

    parser.add_option("--featureType",
                      action="store", dest="featureType",
                      help="Type of feature ('poly' or 'point')," +
                      " default='point'",
                      default='point')

    parser.add_option("--fieldName", metavar='FIELD',
                      action="store", dest="fieldName",
                      help="write data in shape file to FIELD," +
                      " default='value'",
                      default='value')

    parser.add_option("--filter",
                      action="store", dest="filter",
                      help="Filter out data equal or below" +
                      " limit in shape output",
                      default=None)

    (options, args) = parser.parse_args()

    #------------Setting up logging capabilities -----------
    rootLogger = logger.RootLogger(int(options.loglevel))
    global log
    log = rootLogger.getLogger(sys.argv[0])

    #------------Process and validate options---------------
    if options.doc:
        print __doc__
        sys.exit()

    if len(args) > 0:
        parser.error("Incorrect number of arguments")

    #validate infile path
    if options.infileName is not None:
        inFilePath = options.infileName
        if not path.exists(inFilePath):
            log.error("Input raster does not exist")
            sys.exit(1)
    else:
        log.error("No input file specified")
        sys.exit(1)

    #validate outfile path
    if options.outfileName is not None:
        outFilePath = path.abspath(options.outfileName)
        if ".shp" not in outFilePath:
            log.error("Output shape must have .shp extension")
            sys.exit(1)
    else:
        log.error("No output file specified")
        sys.exit(1)

    #Validate fieldName option
    fieldName = options.fieldName
    if fieldName == "":
        parser.error("fieldName can't be an empty string")

    #Validate filter option and convert filter to numeric value if present
    if options.filter is not None:
        filter = float(options.filter)
    else:
        filter = None

    #validate featureType
    if options.featureType not in ('point', 'poly'):
        log.error("Invalid feature type, must be either 'point' or 'poly'")
        sys.exit(1)

    #Assure that gdal is present
    if not __gdal_loaded__:
        raise OSError("Function readGDAL needs GDAL with python bindings")

    # register all of the raster drivers
    gdal.AllRegister()
    ds = gdal.Open(inFilePath, GA_ReadOnly)
    if ds is None:
        log.error('Could not open ' + inFilePath)
        sys.exit(1)

    ncols = ds.RasterXSize
    nrows = ds.RasterYSize
    nbands = ds.RasterCount

    #Info used for georeferencing
    geoTransform = ds.GetGeoTransform()
    xul = geoTransform[0]  # top left x
    cellsizeX = geoTransform[1]  # w-e pixel resolution
    rot1 = geoTransform[2]  # rotation, 0 if image is "north up"
    yul = geoTransform[3]  # top left y
    rot2 = geoTransform[4]  # rotation, 0 if image is "north up"
    cellsizeY = geoTransform[5]  # n-s pixel resolution
    proj = ds.GetProjection()

    #Calculate lower left corner
    xll = xul
    yll = yul + nrows * cellsizeY  # cellsizeY should be a negative value

    #Rotated rasters not handled...yet
    if rot1 != 0 or rot2 != 0:
        log.error('Handling of rotated rasters are not implemented yet')
        sys.exit(1)

    bandIndex = int(options.bandIndex)
    band = ds.GetRasterBand(bandIndex)

    nodata = band.GetNoDataValue()
    #If no nodata value is present in raster, set to -9999 for completeness
    if nodata is None:
        nodata = -9999

    rowmin = 0
    colmin = 0

    #Processing of data is made for blocks of the following size
    #Important - blocks are set to cover all columns for simpler processing
    #This might not always be optimal for read/write speed
    procXBlockSize = ncols
    procYBlockSize = 1

    #process option for dataType
    dataTypes, defaultDataType = ogrDataTypes, 'Real'

    try:
        dataType = dataTypes[options.dataType or defaultDataType]
    except KeyError:
        log.error(
            "Unknown datatype choose between: %s" % ",".join(dataTypes.keys()))
        sys.exit(1)

    #Create and inititialize output vector data source
    shapeDriver = ogr.GetDriverByName('ESRI Shapefile')
    if path.exists(outFilePath):
        shapeDriver.DeleteDataSource(outFilePath)
    shapeFile = shapeDriver.CreateDataSource(outFilePath)
    if shapeFile is None:
        log.error("Could not open output shapefile %s" % outFilePath)
        sys.exit(1)
    if options.featureType == "point":
        layer = shapeFile.CreateLayer(outFilePath, geom_type=ogr.wkbPoint)
    else:
        layer = shapeFile.CreateLayer(outFilePath, geom_type=ogr.wkbPolygon)

    fieldDefn = ogr.FieldDefn(fieldName, dataType)
    layer.CreateField(fieldDefn)

    #Loop over block of raster (at least one row in each block)
    rowsOffset = 0
    pg = ProgressBar(nrows, options.progressStream)

    for i in range(0, nrows, procYBlockSize):
        pg.update(i)
        data = band.ReadAsArray(xoff=colmin, yoff=rowmin + i,
                                win_xsize=procXBlockSize,
                                win_ysize=procYBlockSize)

        blockYll = yul + i * cellsizeY  # newCellsizeY is negative
        blockXll = xll
        if options.featureType == 'point':
            block2point(data, layer, blockXll, blockYll, cellsizeX,
                        cellsizeY, nodata, fieldName, filter)
        else:
            block2poly(data, layer, blockXll, blockYll, cellsizeX,
                       cellsizeY, nodata, fieldName, filter)

        rowsOffset += procYBlockSize  # Update offset

    shapeFile.Destroy()
    pg.finished()

if __name__ == "__main__":
    main()
    # try:
    #     main()
    # except:
    #     import pdb, sys
    #     e, m, tb = sys.exc_info()
    #     pdb.post_mortem(tb)
