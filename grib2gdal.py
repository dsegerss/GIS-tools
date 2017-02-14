#!/usr/bin/env python
# -*- coding: us-ascii -*-
"""
Convert MATCH grib fields to GDAL format of choice.

Created: 09 Oct 2014
Author: David Segersson
"""

import sys
import logging
import argparse

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

DEFAULT_RASTER_FORMAT = "GTiff"
DENSITY_AIR = 1.22  # kg/m3
UNIT_CONVERSION_FACTOR = 1.0e6 * DENSITY_AIR

log = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '-p', '--param',
        action='store',
        dest='param',
        help='param to extract'
    )

    parser.add_argument(
        '-l', '--level',
        action='store',
        dest='level',
        type=int,
        help='Level to extract'
    )

    parser.add_argument(
        '-s', '--sort',
        action='store',
        dest='sort',
        help='Sort to extract'
    )

    parser.add_argument(
        '-t', '--timeinterp',
        action='store',
        dest='timeinterp',
        help='Timeinterp to extract'
    )

    parser.add_argument(
        '-i',
        action='store',
        dest='infile',
        help='Input grib file'
    )

    parser.add_argument(
        '-o',
        action='store',
        dest='outfile',
        help='Output gdal file'
    )

    parser.add_argument(
        '--format',
        action='store',
        dest='format',
        help='Output format, default is %s' % DEFAULT_RASTER_FORMAT,
        default=DEFAULT_RASTER_FORMAT
    )

    parser.add_argument(
        '--unit-conversion',
        action='store',
        type=float,
        dest='unit_conversion_factor',
        help='Factor for unit conversion',
        default=1.0
    )

    parser.add_argument(
        '-v',
        action=VerboseAction, dest='loglevel', default=logging.WARNING,
        help='increase verbosity in terminal'
    )
    
    args = parser.parse_args()

    # Assure that gdal is present
    if not __gdal_loaded__:
        log.error("Function readGDAL needs GDAL with python bindings")
        sys.exit(1)

    grib = pg.open(args.infile)

    if not grib.has_par(args.param, lev=args.level):
        log.error(
            "Parameter %s and level %i not found in grib: %s" % (
                args.param, args.level, args.infile)
        )
        sys.exit(1)

    grib_field = grib.get(
        args.param,
        lev=args.level,
        sort=args.sort,
        time=args.timeinterp
    )

    gdal_field = np.flipud(
        grib_field.values.transpose()
    )

    gdal_field *= args.unit_conversion_factor

    # to convert to microg/m3, multiply by UNIT_CONVERSION_FACTOR

    # register all of the raster drivers
    gdal.AllRegister()

    mem_ds = gdal.GetDriverByName('MEM').Create(
        args.outfile,
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
    output_driver = gdal.GetDriverByName(args.format)
    output_ds = output_driver.CreateCopy(args.outfile, mem_ds, 0)
    output_ds = None
    mem_ds = None


class VerboseAction(argparse.Action):
    
    """Argparse action to handle terminal verbosity level."""
    
    def __init__(self, option_strings, dest,
                 default=logging.WARNING, help=None):
        baselogger = logging.getLogger('')
        baselogger.setLevel(logging.DEBUG)
        self._loghandler = logging.StreamHandler()
        self._loghandler.setLevel(default)
        format = ': '.join((sys.argv[0], '%(levelname)s', '%(message)s'))
        streamformatter = logging.Formatter(format)
        self._loghandler.setFormatter(streamformatter)
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


if __name__ == '__main__':
    main()
