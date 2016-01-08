#!/usr/bin/env python
# -*- coding: latin-1 -*-
"""
Expand raster with a buffer distance by padding with no data.
Edges of expanded raster will be given a value corresponding
to the average of outer edges.
"""

from __future__ import division

# Standard modules
import os
from os import path
import sys
import argparse
import logging

import numpy as np

try:
    from osgeo import gdal
    from osgeo.gdalconst import GDT_Float32, GA_ReadOnly
    __gdal_loaded__ = True
except:
    __gdal_loaded__ = False


log = logging.getLogger(__name__)


def create_terminal_handler(loglevel=logging.INFO, prog=None):
    """Configure a log handler for the terminal."""
    if prog is None:
        prog = os.path.basename(sys.argv[0])
    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(loglevel)
    format = ': '.join((prog, '%(levelname)s', '%(message)s'))
    streamformatter = logging.Formatter(format)
    streamhandler.setFormatter(streamformatter)
    return streamhandler


class VerboseAction(argparse.Action):

    """Argparse action to handle terminal verbosity level."""

    def __init__(self, option_strings, dest,
                 default=logging.WARNING, help=None):
        baselogger = logging.getLogger(__name__)
        baselogger.setLevel(logging.DEBUG)
        self._loghandler = create_terminal_handler(default)
        baselogger.addHandler(self._loghandler)
        super(VerboseAction, self).__init__(
            option_strings, dest,
            nargs=0,
            default=default,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        currentlevel = getattr(namespace, self.dest, logging.WARNING)
        self._loghandler.setLevel(currentlevel - 10)
        setattr(namespace, self.dest, self._loghandler.level)


def get_border_average(data, nodata):
    nrows, ncols = data.shape
    boundary_cells = set()
    for row in range(nrows):
        for col in range(ncols):
            if data[row, col] != nodata:
                if (row, col) not in boundary_cells:
                    boundary_cells.add((row, col))
                break
    for row in range(nrows):
        for col in range(ncols - 1, 0):
            if data[row, col] != nodata:
                if (row, col) not in boundary_cells:
                    boundary_cells.add((row, col))
                break

    for col in range(ncols):
        for row in range(nrows):
            if data[row, col] != nodata:
                if (row, col) not in boundary_cells:
                    boundary_cells.add((row, col))
                break

    for col in range(ncols):
        for row in range(nrows - 1, 0):
            if data[row, col] != nodata:
                if (row, col) not in boundary_cells:
                    boundary_cells.add((row, col))
                break
    boundary_sum = 0
    for row, col in boundary_cells:
        boundary_sum += data[row, col]
    return boundary_sum / len(boundary_cells)


def main():
    # -----------Setting up and unsing option parser----------------
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '-v',
        action=VerboseAction, dest='loglevel',
        default=logging.WARNING,
        help='increase verbosity in terminal',
    )

    parser.add_argument(
        "-i", "--input", required=True,
        action="store", dest="infile",
        help="Input raster"
    )

    parser.add_argument(
        "--bandIndex", type=int,
        action="store", dest="band",
        help="Band index to read from",
        default=1
    )

    parser.add_argument(
        "-o", "--outfile", required=True,
        action="store", dest="outfile", default=None,
        help="Output raster file"
    )

    parser.add_argument(
        "--buffer", required=True, type=float,
        action="store", dest="buffer",
        help="Add buffer distance"
    )

    args = parser.parse_args()

    if not path.exists(args.infile):
        log.error("Input raster does not exist")
        sys.exit(1)

    # Assure that gdal is present
    if not __gdal_loaded__:
        raise OSError("Function readGDAL needs GDAL with python bindings")

    # register all of the raster drivers
    gdal.AllRegister()
    ds = gdal.Open(str(args.infile), GA_ReadOnly)
    if ds is None:
        log.error('Could not open ' + args.infile)
        sys.exit(1)

    ncols = ds.RasterXSize
    nrows = ds.RasterYSize
    # nbands = ds.RasterCount

    # Info used for georeferencing
    geoTransform = ds.GetGeoTransform()
    xul = geoTransform[0]  # top left x
    cellsizex = geoTransform[1]  # w-e pixel resolution
    rot1 = geoTransform[2]  # rotation, 0 if image is "north up"
    yul = geoTransform[3]  # top left y
    rot2 = geoTransform[4]  # rotation, 0 if image is "north up"
    cellsizey = geoTransform[5]  # n-s pixel resolution
    # proj = ds.GetProjection()

    # Calculate lower left corner
    # xll = xul
    # yll = yul + nrows * cellsizey  # cellsizeY should be a negative value

    # Rotated rasters not handled...yet
    if rot1 != 0 or rot2 != 0:
        log.error('Handling of rotated rasters are not implemented yet')
        sys.exit(1)

    band = ds.GetRasterBand(args.band)
    nodata = band.GetNoDataValue()
    # If no nodata value is present in raster, set to -9999 for completeness
    if nodata is None:
        nodata = -9999

    data = band.ReadAsArray(
        0, 0, ds.RasterXSize,
        ds.RasterYSize).astype(np.float)

    border_avg = get_border_average(data, nodata)

    buffercellsx = int(float(args.buffer) / cellsizex)
    buffercellsy = -1 * int(float(args.buffer) / cellsizey)
    new_ncols = ncols + 2 * buffercellsx
    new_nrows = nrows + 2 * buffercellsy
    new_xll = xul - buffercellsx * cellsizex
    new_yul = yul - buffercellsy * cellsizey
    mem_ds = gdal.GetDriverByName(b'MEM').Create(
        str(args.outfile),
        new_ncols,
        new_nrows,
        1,
        GDT_Float32
    )
    outGeotransform = [
        xul - buffercellsx * cellsizex,
        cellsizex,
        0,
        new_yul,
        0,
        cellsizey
    ]
    mem_ds.SetGeoTransform(outGeotransform)
    out_band = mem_ds.GetRasterBand(1)
    out_band.SetNoDataValue(nodata)
    new_data = np.ones((new_nrows, new_ncols)) * nodata

    new_data[buffercellsy:buffercellsy + nrows,
             buffercellsx:buffercellsx + ncols] = data

    new_data[0, :] = border_avg
    new_data[-1, :] = border_avg
    new_data[:, 0] = border_avg
    new_data[:, -1] = border_avg

    out_band.WriteArray(new_data, 0, 0)
    out_band.FlushCache()  # Write data to disk
    out_driver = ds.GetDriver()
    out_ds = out_driver.CreateCopy(str(args.outfile), mem_ds, 0)
    out_ds = None
    mem_ds = None


if __name__ == "__main__":
    main()
