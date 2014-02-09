#!/usr/bin/env python
# -*- coding: latin-1 -*-
"""
***********************************************************************
Name: raster2STL.py
Created: 28 Oct 2013
Author: David Segersson

Description
-----------------------------------------------------------------------
Convert raster into STL-format
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


class Point:
    def __init__(self, coords):
        self.coords = [coords[0], coords[1], coords[2]]

    def __getitem__(self, i):
        return self.coords[i]

    def __setitem__(self, i, val):
        self.coords[i] = val

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[0]

    @property
    def z(self):
        return self.coords[0]

    def round(self, precision):
        #rounds to specified number of decimals
        for i in range(len(self.coords)):
            self.coords[i] = round(self.coords[i], precision)

    def asvector(self):
        return (self.coords[0], self.coords[1], self.coords[2])

    def __eq__(self, p2):
        try:
            if self.coords == p2.coords:
                return True
            else:
                return False
        except:
            return False

    def __sub__(self, p2):
        p = Point((self.x - p2.x, self.y - p2.y, self.z - p2.z))
        return p

    def __add__(self, p2):
        p = Point((self.x + p2.x, self.y + p2.y, self.z + p2.z))
        return p

    def __mul__(self, p2):
        p = Point((self.x * p2.x, self.y * p2.y, self.z * p2.z))
        return p


class Triangle:

    def __init__(self, p1, p2, p3):
        self.points = [p1, p2, p3]

    def normal(self):
        #vectors that defines the plane
        a = self.p2() - self.p1()
        b = self.p3() - self.p1()
        #cross product p1*p2
        x = a[1] * b[2] - a[2] * b[1]
        y = a[2] * b[0] - a[0] * b[2]
        z = a[0] * b[1] - a[1] * b[0]
        return (x, y, z)

    def p1(self):
        return self.points[0]

    def p2(self):
        return self.points[1]

    def p3(self):
        return self.points[2]

    def round(self, precision):
        for p in self.points:
            p.round(precision)

    def __str__(self):
        t1_str = "  facet normal %.6e %.6e %.6e\n" % self.normal()
        t1_str += "    outer loop\n"
        t1_str += "      vertex %.6e %.6e %.6e\n" % self.p1().asvector()
        t1_str += "      vertex %.6e %.6e %.6e\n" % self.p2().asvector()
        t1_str += "      vertex %.6e %.6e %.6e\n" % self.p3().asvector()
        t1_str += "    endloop\n"
        t1_str += "  endfacet\n"
        return t1_str


def getCentreCoords(row, col, xll, yll, cellsizeX, cellsizeY, nrows, ncols):
    """Return cell centre coordinates"""
    x = xll + 0.5 * cellsizeX + col * cellsizeX
    y = yll - 0.5 * cellsizeY - (nrows - row - 1) * cellsizeY
    return (x, y)


def damp(boundary, avgHeight, dist, bufferDist):
    """Return array of equal size as boundary damped with distance"""
    if np.isscalar(boundary):
        diff = avgHeight - boundary
    else:
        diff = avgHeight * np.ones(boundary.shape) - boundary
    dampingFactor = min(1, (dist / (0.5 * bufferDist)))
    return boundary + diff * dampingFactor


def block2STL(block, stream, xll, yll, cellsizeX, cellsizeY,
              nodata, precision):
    """
    Write a raster block to a stream in STL format
    @param block: raster block
    @param stream: a file-like object
    @param xll: x-coordinate of lower left block corner
    @param yll: y-coordinate of lower left block corner
    @param cellsizeX: cellsize in X-direction
    @param cellsizeY: cellsize in Y-direction
    @param nodata: nodata value
    @param precision: precision of STL-file coordinates
    """
    nrows, ncols = block.shape
    for row in np.arange(1, nrows):
        for col in np.arange(ncols - 1):
            x_c_upper, y_c_upper = getCentreCoords(row - 1, col, xll, yll,
                                                   cellsizeX, cellsizeY,
                                                   nrows, ncols)
            x_c_right, y_c_upper = getCentreCoords(row - 1, col + 1, xll, yll,
                                                   cellsizeX, cellsizeY,
                                                   nrows, ncols)
            x_c, y_c = getCentreCoords(row, col, xll, yll,
                                       cellsizeX, cellsizeY,
                                       nrows, ncols)
            z_c = block[row, col]
            z_c_upper = block[row - 1, col]
            z_c_upper_right = block[row - 1, col + 1]
            z_c_right = block[row, col + 1]

            p = Point((x_c, y_c, z_c))
            p_right = Point((x_c_right, y_c, z_c_right))
            p_upper_right = Point((x_c_right, y_c_upper, z_c_upper_right))
            p_upper = Point((x_c, y_c_upper, z_c_upper))

            t1 = Triangle(p, p_upper, p_upper_right)
            t2 = Triangle(p, p_upper_right, p_right)

            t1.round(precision)
            if (z_c != nodata and z_c_upper != nodata
                and z_c_upper_right != nodata):
                stream.writelines(str(t1))
            if (z_c != nodata and z_c_upper_right != nodata
                and z_c_right != nodata):
                stream.writelines(str(t2))


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

    parser.add_option("-i", "--input",
                      action="store", dest="infileName",
                      help="Input raster")

    parser.add_option("--bandIndex",
                      action="store", dest="bandIndex",
                      help="Band index to read from",
                      default=1)

    parser.add_option("-o", "--output",
                      action="store", dest="output", default=None,
                      help="Output file")

    parser.add_option("--solidName", default="TOPO",
                      action="store", dest="solidName",
                      help="Name of solid in STL-output")

    parser.add_option("--precision", default=6,
                      action="store", dest="precision",
                      help="Precision of coordinates in STL-output")

    parser.add_option("--buffer",
                      action="store", dest="buffer",
                      help="Add buffer distance to topo, " +
                      "damping out differences in elevation")


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
    if options.output is None:
        log.error("No output file specified")
        sys.exit(1)

    if options.buffer is not None:
        bufferDist = float(options.buffer)

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
    #nbands = ds.RasterCount

    #Info used for georeferencing
    geoTransform = ds.GetGeoTransform()
    xul = geoTransform[0]  # top left x
    cellsizeX = geoTransform[1]  # w-e pixel resolution
    rot1 = geoTransform[2]  # rotation, 0 if image is "north up"
    yul = geoTransform[3]  # top left y
    rot2 = geoTransform[4]  # rotation, 0 if image is "north up"
    cellsizeY = geoTransform[5]  # n-s pixel resolution
    #proj = ds.GetProjection()

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

    #Processing of data is made for blocks of the following size
    #two rows are processed at a time since triangles
    # are created between cell centres.
    procXBlockSize = ncols
    procYBlockSize = 2

    with open(options.output, 'w') as stlFile:
        pg = ProgressBar(nrows, options.progressStream)

        stlFile.write("solid " + options.solidName + "\n")

        avgHeight = 0
        #write original cells to STL
        #Loop over blocks of raster

        for i in range(0, nrows - 1, 1):
            pg.update(i)
            data = band.ReadAsArray(xoff=0, yoff=i,
                                    win_xsize=procXBlockSize,
                                    win_ysize=procYBlockSize)
            if nodata in data:
                log.error("Nodata value found in raster," +
                          " this should be interpolated before" +
                          " converting to STL")
                sys.exit(1)

            # newCellsizeY is negative
            blockYll = yul + (procYBlockSize + i) * cellsizeY
            blockXll = xll

            block2STL(data, stlFile, blockXll, blockYll, cellsizeX,
                      cellsizeY, nodata, int(options.precision))

            avgHeight += np.mean(data)

        avgHeight /= (nrows - 1)

        if options.buffer is not None:
            # overlap of one cell since triangles are created between
            # cell centres
            bufferCells = int(bufferDist / cellsizeX) + 1

            #write buffer cells to STL
            leftBoundary = band.ReadAsArray(xoff=0, yoff=0,
                                         win_xsize=1,
                                         win_ysize=nrows)
            leftBuffer = leftBoundary[:]

            for i in range(1, bufferCells):
                dist = i * cellsizeX
                leftBuffer = np.hstack((
                    damp(leftBoundary, avgHeight, dist, bufferDist),
                    leftBuffer))

            block2STL(leftBuffer, stlFile, xll - (bufferCells - 1) * cellsizeX,
                      yll, cellsizeX, cellsizeY, nodata,
                      int(options.precision))

            rightBoundary = band.ReadAsArray(xoff=ncols - 1, yoff=0,
                                         win_xsize=1,
                                         win_ysize=nrows)


            rightBuffer = rightBoundary[:]
            for i in range(1, bufferCells):
                dist = i * cellsizeX
                rightBuffer = np.hstack((
                    rightBuffer,
                    damp(rightBoundary, avgHeight, dist, bufferDist)))

            block2STL(rightBuffer, stlFile, xll + (ncols - 1) * cellsizeX,
                      yll, cellsizeX, cellsizeY, nodata,
                      int(options.precision))

            topBoundary = band.ReadAsArray(xoff=0, yoff=0,
                                         win_xsize=procXBlockSize,
                                         win_ysize=1)
            topBuffer = topBoundary[:]
            for i in range(1, bufferCells):
                dist = abs(i * cellsizeY)
                topBuffer = np.vstack((
                        damp(topBoundary, avgHeight, dist, bufferDist),
                        topBuffer))

            block2STL(topBuffer, stlFile, xll,
                      yul + cellsizeY,
                      cellsizeX, cellsizeY, nodata,
                      int(options.precision))

            bottomBoundary = band.ReadAsArray(xoff=0, yoff=nrows - 1,
                                         win_xsize=procXBlockSize,
                                         win_ysize=1)
            bottomBuffer = bottomBoundary[:]
            for i in range(1, bufferCells):
                dist = abs(i * cellsizeY)
                bottomBuffer = np.vstack((
                        bottomBuffer,
                        damp(bottomBoundary, avgHeight, dist, bufferDist)))

            block2STL(bottomBuffer, stlFile, xll,
                      yll + (bufferCells - 1) * cellsizeY, cellsizeX,
                      cellsizeY, nodata,
                      int(options.precision))

            #Add corner buffers
            cornerBlock = np.ones((bufferCells, bufferCells))

            #upper right
            blockXll = xll + (ncols - 1) * cellsizeX
            blockYll = yul + cellsizeY
            startHeight = rightBoundary[0, 0]
            cornerPoint = (xll + (ncols - 0.5) * cellsizeX,
                           yul + 0.5 * cellsizeY)

            for col in range(bufferCells):
                for row in range(bufferCells):
                    centre = getCentreCoords(row, col, blockXll, blockYll,
                                             cellsizeX, cellsizeY,
                                             bufferCells, bufferCells)

                    dist = (np.sqrt(pow(centre[0] - cornerPoint[0], 2) +
                                    pow(centre[1] - cornerPoint[1], 2)))

                    #dist = min(row, col) * cellsizeX
                    cornerBlock[row, col] = damp(startHeight, avgHeight,
                                                 dist, bufferDist)

            cornerBlock[:, 0] = topBuffer[:, -1]
            cornerBlock[-1, :] = rightBuffer[0, :]

            block2STL(cornerBlock, stlFile, blockXll,
                      blockYll, cellsizeX, cellsizeY, nodata,
                      int(options.precision))

            #lower left
            blockXll = xll - (bufferCells - 1) * cellsizeX
            blockYll = yll + (bufferCells - 1) * cellsizeY
            startHeight = leftBoundary[-1, 0]
            cornerPoint = (xll + 0.5 * cellsizeX,
                           yll - 0.5 * cellsizeY)
            for col in range(bufferCells):
                for row in range(bufferCells):
                    centre = getCentreCoords(row, col, blockXll, blockYll,
                                             cellsizeX, cellsizeY,
                                             bufferCells, bufferCells)

                    dist = (np.sqrt(pow(centre[0] - cornerPoint[0], 2) +
                                    pow(centre[1] - cornerPoint[1], 2)))

                    cornerBlock[row, col] = damp(startHeight, avgHeight,
                                                 dist, bufferDist)
            cornerBlock[:, -1] = bottomBuffer[:, 0]
            cornerBlock[0, :] = leftBuffer[-1, :]
            block2STL(cornerBlock, stlFile, blockXll,
                      blockYll, cellsizeX, cellsizeY, nodata,
                      int(options.precision))

            #upper left
            blockXll = xll - (bufferCells - 1) * cellsizeX
            blockYll = yul + cellsizeY
            startHeight = leftBoundary[0, 0]
            cornerPoint = (xll + 0.5 * cellsizeX,
                           yul + 0.5 * cellsizeY)
            for col in range(bufferCells):
                for row in range(bufferCells):
                    centre = getCentreCoords(row, col, blockXll, blockYll,
                                             cellsizeX, cellsizeY,
                                             bufferCells, bufferCells)

                    dist = (np.sqrt(pow(centre[0] - cornerPoint[0], 2) +
                                    pow(centre[1] - cornerPoint[1], 2)))

                    cornerBlock[row, col] = damp(startHeight, avgHeight,
                                                 dist, bufferDist)
            cornerBlock[:, -1] = topBuffer[:, 0]
            cornerBlock[-1, :] = leftBuffer[0, :]
            block2STL(cornerBlock, stlFile, blockXll,
                      blockYll, cellsizeX, cellsizeY, nodata,
                      int(options.precision))

            #lower right
            blockXll = xll + (ncols - 1) * cellsizeX
            blockYll = yll + (bufferCells - 1) * cellsizeY
            startHeight = rightBoundary[-1, 0]
            cornerPoint = (xll + (ncols - 0.5) * cellsizeX,
                           yll - 0.5 * cellsizeY)
            for col in range(bufferCells):
                for row in range(bufferCells):
                    centre = getCentreCoords(row, col, blockXll, blockYll,
                                             cellsizeX, cellsizeY,
                                             bufferCells, bufferCells)

                    dist = (np.sqrt(pow(centre[0] - cornerPoint[0], 2) +
                                    pow(centre[1] - cornerPoint[1], 2)))

                    cornerBlock[row, col] = damp(startHeight, avgHeight,
                                                 dist, bufferDist)

            cornerBlock[:, 0] = bottomBuffer[:, -1]
            cornerBlock[0, :] = rightBuffer[-1, :]
            block2STL(cornerBlock, stlFile, blockXll,
                      blockYll, cellsizeX, cellsizeY, nodata,
                      int(options.precision))

        stlFile.write("endsolid %s\n" % options.solidName)

    pg.finished()

if __name__ == "__main__":
    main()
    # try:
    #     main()
    # except:
    #     import pdb, sys
    #     e, m, tb = sys.exc_info()
    #     pdb.post_mortem(tb)
