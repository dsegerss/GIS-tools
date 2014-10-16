#!/usr/bin/env python
'''
point2raster.py
Created on 13 jan 2013
@author: David Segersson
'''
import sys
import logging
from optparse import OptionParser
from os import path

import numpy


from osgeo import ogr
from osgeo import gdal
from osgeo.gdalconst import GDT_Float32


#Docstrings for the option parser
usage = ("usage: %prog " +
         "--xll <XLL> --yll <YLL> --nx <NX> --ny <NY> -s <SIZE> " +
         "\n" + " " * 23 +
         "[-m <sum|mean> -f <GDAL format id> --field <Field>] " +
         "\n" + " " * 23 +
         "<input shape> <output raster> ")
version = "%prog 1.0"


def main():
    #-----------Setting up and unsing option parser-----------------------
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("--field",
                      action="store", dest="field", metavar='FIELD',
                      help="field to read data from, " +
                      "default is to output number of points per cell",
                      default=None)

    parser.add_option("--xll",
                      action="store", dest="xll",
                      help="X-coord of lower left corner in output raster")

    parser.add_option("--yll",
                      action="store", dest="yll",
                      help="Y-coord of lower left corner in output raster")

    parser.add_option("--nx",
                      action="store", dest="nx",
                      help="Number of columns in output raster")

    parser.add_option("--ny",
                      action="store", dest="ny",
                      help="Number of rows in output raster")

    parser.add_option("-s", "--cellsize",
                      action="store", dest="cellsize",
                      help="Cellsize of output raster")

    parser.add_option("-m", "--method",
                      action="store", dest="method",
                      help="Add to cell as 'sum' or 'mean', default is 'sum'")

    parser.add_option("-f", "--format",
                      action="store", dest="format",
                      help="Output raster format, e.g. 'GTiff' or 'AAIGrid'")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Print statistics when finished")

    (options, args) = parser.parse_args()

    if options.verbose:
        logging.basicConfig(level=logging.INFO)

    if len(args) != 2:
        parser.error("Wrong number of arguments")

    #validate options
    if options.xll is None:
        parser.error("Must specify option --xll")
    if options.yll is None:
        parser.error("Must specify option --yll")
    if options.nx is None:
        parser.error("Must specify option --nx")
    if options.ny is None:
        parser.error("Must specify option --ny")
    if options.cellsize is None:
        parser.error("Must specify option --cellsize")
    if options.method is None:
        method = 'sum'
    elif options.method in ('sum', 'mean'):
        method = options.method
    else:
        parser.error("Invalid method. choose between 'sum' and 'mean'")
    if options.format is None:
        raster_format = "GTiff"
    else:
        raster_format = options.format

    xll = float(options.xll)
    yll = float(options.yll)
    nx = int(options.nx)
    ny = int(options.ny)
    cellsize = float(options.cellsize)
    nodata = -9999

    yul = yll + ny * cellsize
    xur = xll + nx * cellsize

    outArray = numpy.zeros((ny, nx))
    nvals = numpy.zeros((ny, nx))

    #validate infile path
    shapefile = args[0]
    if not path.exists(shapefile):
        logging.error("Input shape file does not exist")
        sys.exit(1)
    elif ".shp" not in shapefile:
        logging.error("Input shape has to to be specified with" +
                  " .shp extension")
        sys.exit(1)

    outfilename = args[1]

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shapefile, update=0)
    if dataSource is None:
        logging.error("Could not open data source: " + shapefile)
        sys.exit(1)

    layer = dataSource.GetLayer()  # Create a layer for the shapefile
    srs = layer.GetSpatialRef()
    npoints = layer.GetFeatureCount()
    feature = layer.GetNextFeature()  # get first polygon in shape file
    noutside = 0
    while feature:
        point = feature.GetGeometryRef()
        x = point.GetX()
        y = point.GetY()
        if options.field is not None:
            try:
                val = float(feature.GetFieldAsString(options.field))
            except:
                logging.error("Could not read field value from shape")
        else:
            val = 1

        #check if outside target extent
        if((x < xll or y < yll) or
           (x > xur or y > yul)):
            feature.Destroy()
            feature = layer.GetNextFeature()
            noutside += 1
            continue

        col = int((x - xll) / cellsize)
        row = ny - int(
            numpy.ceil((y - yll) / cellsize))

        if row == ny:
            row -= 1
        if col == nx:
            col -= 1
        if val != nodata:
            outArray[row, col] += val
            nvals[row, col] += 1

        feature.Destroy()
        feature = layer.GetNextFeature()
    dataSource.Destroy()

    gdal.AllRegister()
    mem_ds = gdal.GetDriverByName('MEM').Create(outfilename,
                                                nx,
                                                ny,
                                                1,
                                                GDT_Float32)
    outGeotransform = [xll, cellsize, 0, yul, 0, -1 * cellsize]
    mem_ds.SetGeoTransform(outGeotransform)
    if srs is not None:
        mem_ds.SetProjection(srs.ExportToWkt())
    outBand = mem_ds.GetRasterBand(1)
    outBand.SetNoDataValue(nodata)

    if method == 'mean':
        outArray = outArray / nvals

    outBand.WriteArray(outArray, 0, 0)
    outBand.FlushCache()
    outputDriver = gdal.GetDriverByName(raster_format)
    output_ds = outputDriver.CreateCopy(outfilename, mem_ds, 0)

    # Once we're done, close properly the dataset
    output_ds = None
    mem_ds = None
    logging.info("Finished, %i points outside raster extent" % noutside)
    logging.info("Output raster sum: %f" % outArray.sum())

if __name__ == "__main__":
    main()
