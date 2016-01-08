#!/usr/bin/env python
# -*- coding: us-ascii -*-
"""
Name: grib2gdal.py
Created: 09 Oct 2014
Author: David Segersson

Description
Convert MATCH grib fields to GDAL format of choice.
"""

import sys
import logging
from optparse import OptionParser

import numpy as np

try:
    from osgeo import osr
    from osgeo import gdal
    from osgeo.gdalconst import *
    __gdal_loaded__ = True
except ImportError:
    __gdal_loaded__ = False

try:
    import pygrib as pg
    __pygrib_loaded__ = True
except ImportError:
    __pygrib_loaded__ = False


#Docstrings for the option parser
usage = "usage: %prog [options] "
version = "%prog 1.0"

DEFAULT_RASTER_FORMAT = "GTiff"
DENSITY_AIR = 1.22  # kg/m3
UNIT_CONVERSION_FACTOR = 1.0e6 * DENSITY_AIR


def main():
    #-----------Setting up and unsing option parser-----------------------
    parser = OptionParser(usage=usage, version=version)

    parser.add_option(
        '-v', '--verbose',
        action='store_true', dest='verbose',
        help='produce verbose output'
    )

    parser.add_option(
        '-p', '--param',
        action='store',
        dest='param',
        help='param to extract'
    )

    parser.add_option(
        '-l', '--level',
        action='store',
        dest='level',
        type=int,
        help='Level to extract'
    )

    parser.add_option(
        '-s', '--sort',
        action='store',
        dest='sort',
        help='Sort to extract'
    )

    parser.add_option(
        '-t', '--timeinterp',
        action='store',
        dest='timeinterp',
        help='Timeinterp to extract'
    )

    parser.add_option(
        '-i',
        action='store',
        dest='infile',
        help='Input grib file'
    )

    parser.add_option(
        '-o',
        action='store',
        dest='outfile',
        help='Output gdal file'
    )

    parser.add_option(
        '--format',
        action='store',
        dest='format',
        help='Output format, default is %s' % DEFAULT_RASTER_FORMAT,
        default = DEFAULT_RASTER_FORMAT
    )

    (options, args) = parser.parse_args()
    if options.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.WARNING,

    logging.basicConfig(level=loglevel)
    log = logging.getLogger(__name__)

    #Assure that gdal is present
    if not __gdal_loaded__:
        log.error("Function readGDAL needs GDAL with python bindings")
        sys.exit(1)

    grib = pg.open(options.infile)
    if not grib.has_par(options.param, lev=options.level):
        log.error(
            "Parameter %s and level %i not found in grib: %s" % (
                par, lev, options.infile)
        )
        sys.exit(1)

    grib_field = grib.get(
        options.param,
        lev=options.level,
        sort=options.sort,
        time=options.timeinterp
    )

    gdal_field = np.flipud(
        grib_field.values.transpose()) * UNIT_CONVERSION_FACTOR
    # register all of the raster drivers
    gdal.AllRegister()

    mem_ds = gdal.GetDriverByName('MEM').Create(
        options.outfile,
        gdal_field.shape[1],
        gdal_field.shape[0],
        1,
        GDT_Float32)

    srs = osr.SpatialReference()
    srs.ImportFromProj4(grib_field['proj4'])
    mem_ds.SetProjection(srs.ExportToWkt())
    outGeotransform = [
        1035000,
        grib_field.dlon,
        0,
        6816500,
        0,
        -1 * grib_field.dlat
    ]
    mem_ds.SetGeoTransform(outGeotransform)
    band = mem_ds.GetRasterBand(1)
    band.SetNoDataValue(-9999)
    band.WriteArray(gdal_field)
    band.FlushCache()
    output_driver = gdal.GetDriverByName(options.format)
    output_ds = output_driver.CreateCopy(options.outfile, mem_ds, 0)
    output_ds = None
    mem_ds = None

if __name__ == '__main__':
    main()
