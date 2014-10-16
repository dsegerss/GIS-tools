#!/usr/bin/env python
"""
Name: gdal_contourf.py
Created: 09 Feb 2014
Author: David Segersson

Description
Produce isolines or isobands in any OGR-supported vector format
from any GDAL supported raster format.
"""

from __future__ import division, unicode_literals
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import sys
import logging
from optparse import OptionParser

try:
    from osgeo import ogr
    from osgeo import gdal
    __gdal_loaded__ = True
except ImportError:
    __gdal_loaded__ = False

try:
    from scipy import interpolate
    import numpy
    import matplotlib.pyplot as plt
    __matplotlib_loaded__ = True
except ImportError:
    __matplotlib_loaded__ = False

#Docstrings for the option parser
usage = "usage: %prog [options] "
version = "%prog 1.0"


def main():
    #-----------Setting up and unsing option parser-----------------------
    parser = OptionParser(usage=usage, version=version)

    parser.add_option(
        '-v', '--verbose',
        action='store_true', dest='verbose',
        help='produce verbose output'
    )

    parser.add_option(
        '--isolines',
        action='store', dest='isolines',
        help='Output vector file for isolines'
    )

    parser.add_option(
        '--isobands',
        action='store', dest='isobands',
        help='Output vector file for isobands (filled contours)'
    )

    parser.add_option(
        '-f', '--format',
        action='store', dest='format', default='ESRI Shapefile',
        help='Output vector format, default is ESRI Shapefile'
    )

    parser.add_option(
        '--degree',
        action='store', dest='degree', default='3',
        help='Interpolate using splines of given order, default is 3 (cubic)'
    )

    parser.add_option(
        '--resample',
        action='store', dest='resample', default='5',
        help='Resampling factor, default is 5 times'
    )

    parser.add_option(
        '--logarithmic',
        action='store_true',
        dest='logarithmic',
        help="Logarithmic auto levels"
    )

    parser.add_option(
        '--interval',
        action='store',
        dest='interval',
        help='Interval between isolones'
    )

    parser.add_option(
        '--levels',
        action='store',
        dest='levels',
        help="Manual levels as a list 'val1 val2 ... val3'"
    )

    parser.add_option(
        '--nlevels',
        action='store',
        dest='nlevels',
        help="Number of levels"
    )

    parser.add_option(
        '--min_lev',
        action='store',
        dest='min_lev',
        help="Lowest level to be created, default is raster minimum"
    )

    parser.add_option(
        '--max_lev',
        action='store',
        dest='max_lev',
        help="Highest level to be created, default is raster maximum"
    )

    parser.add_option(
        '-i',
        '--input',
        action='store',
        dest='raster',
        help='Input raster'
    )

    (options, args) = parser.parse_args()
    if options.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel)
    log = logging.getLogger(__name__)

    if not __matplotlib_loaded__:
        log.error("Must have matplotlib installed")
        sys.exit(1)

    if len(args) != 0:
        log.error("No argument expected")
        sys.exit(1)

    if options.raster is None:
        log.error("Must specify input raster")
        sys.exit(1)

    if options.isolines is None and options.isobands is None:
        log.error("Must specify at least one of options" +
                  " '--isolines <filename>' or '--isobands <filename>'")
        sys.exit(1)

    if not os.path.exists(options.raster):
        log.error("Input raster does not exist")
        sys.exit(1)

    #read input raster
    ds = gdal.Open(options.raster)
    rast = ds.ReadAsArray()
    gt = ds.GetGeoTransform()
    ds = None
    dx = gt[1]
    dy = gt[5]
    xmin = gt[0]
    ymax = gt[3]
    nx = rast.shape[1]
    ny = rast.shape[0]

    spline_degree = int(options.degree)
    resampling_factor = float(options.resample)

    if options.min_lev is not None:
        min_val = float(options.min_lev)
    else:
        min_val = rast.min() + 0.01 * (rast.max() - rast.min())

    if options.max_lev is not None:
        max_val = float(options.max_lev)
    else:
        max_val = rast.max() - 0.01 * (rast.max() - rast.min())

    if options.interval is not None:
        if options.logarithmic:
            log.error(
                "Cannot specify interval and logarithmic at the same time")
            sys.exit(1)
        interval = float(options.interval)
        nlevels = (max_val - min_val) / interval
    elif options.nlevels is not None:
        nlevels = int(options.nlevels)
    else:
        nlevels = 11

    if options.levels is not None:
        levels = map(float, options.levels.split())
    elif options.logarithmic:
        levels = numpy.logspace(numpy.log10(min_val),
                                numpy.log10(max_val), nlevels)
    else:
        levels = numpy.linspace(min_val, max_val, nlevels)

    #prepare output files
    driver = ogr.GetDriverByName(str(options.format))
    if options.isolines is not None:
        if os.path.exists(options.isolines):
            driver.DeleteDataSource(options.isolines)
        isolinefile = driver.CreateDataSource(str(options.isolines))
        if isolinefile is None:
            raise OSError("Could not open data source: " + options.isolines)
        isoline_layer = isolinefile.CreateLayer(str('isolines'),
                                                geom_type=ogr.wkbLineString)

    if options.isobands is not None:
        if os.path.exists(options.isobands):
            driver.DeleteDataSource(options.isobands)
        isobandfile = driver.CreateDataSource(str(options.isobands))
        if isobandfile is None:
            raise OSError("Could not open data source: " + options.isobands)
        isoband_layer = isobandfile.CreateLayer(str('isobands'),
                                                geom_type=ogr.wkbPolygon)

    #create input raster coordinate vectors
    x = numpy.arange(gt[0], xmin + rast.shape[1] * dx, dx)
    y = numpy.arange(ymax, ymax + rast.shape[0] * dy, dy)

    #resample grid using splines of specified order
    x_interp = numpy.linspace(x[0], x[-1], int(nx * resampling_factor))
    y_interp = numpy.linspace(y[-1], y[0], int(ny * resampling_factor))

    kernel = interpolate.RectBivariateSpline(
        y[::-1], x, numpy.flipud(rast),
        kx=spline_degree, ky=spline_degree)
    rast_interp = kernel(y_interp, x_interp)

    #Mesh to be used for contouring
    x_interp, y_interp = numpy.meshgrid(x_interp, y_interp)

    if options.isobands is not None:
        # (nlevels-1) PathCollections
        isobands = plt.contourf(x_interp, y_interp, rast_interp, levels=levels)
        write_isobands(isoband_layer, isobands)
        isobandfile.Destroy()

    if options.isolines is not None:
        # nlevels LineCollections
        isolines = plt.contour(x_interp, y_interp, rast_interp, levels=levels)
        write_isolines(isoline_layer, isolines)
        isolinefile.Destroy()

    #ax = fig.add_subplot(111)
    #patch = patches.PathPatch(paths[1], facecolor='orange', lw=2)
    #ax.add_patch(patch)
    #ax.colorbar()
    #plt.show()


def write_isobands(layer, isobands):
    layer.CreateField(ogr.FieldDefn(str('min'), ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(str('max'), ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(str('interval'), ogr.OFTString))
    bandFeatureDefn = layer.GetLayerDefn()

    for lev, isoband in enumerate(isobands.collections):
        min_level = isobands.levels[lev]
        max_level = isobands.levels[lev + 1]
        paths = isoband.get_paths()
        for path in paths:
            polygon = ogr.Geometry(ogr.wkbPolygon)
            polys = path.to_polygons()
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for seg in polys[0]:
                x = seg[0]
                y = seg[1]
                ring.AddPoint(x, y)
            ring.CloseRings()
            polygon.AddGeometry(ring)
            for poly in polys[1:]:
                ring = ogr.Geometry(ogr.wkbLinearRing)
                for row in range(len(poly) - 1, 0, -1):
                    x = poly[row][0]
                    y = poly[row][1]
                    ring.AddPoint(x, y)
                ring.CloseRings()
                polygon.AddGeometry(ring)
            feature = ogr.Feature(bandFeatureDefn)
            feature.SetGeometry(polygon)
            feature.SetField(str('min'), min_level)
            feature.SetField(str('max'), max_level)
            feature.SetField(str('interval'),
                             str("%f-%f" % (min_level, max_level)))
            layer.CreateFeature(feature)
            polygon.Destroy()
            feature.Destroy()


def write_isolines(layer, isolines):
    layer.CreateField(ogr.FieldDefn(str('value'), ogr.OFTReal))
    lineFeatureDefn = layer.GetLayerDefn()

    for lev, isoline in enumerate(isolines.collections):
        value = isolines.levels[lev]
        lines = isoline.get_paths()
        for line in lines:
            linestring = ogr.Geometry(ogr.wkbLineString)
            for seg in line.vertices:
                x = seg[0]
                y = seg[1]
                linestring.AddPoint(x, y)
            feature = ogr.Feature(lineFeatureDefn)
            feature.SetGeometry(linestring)
            feature.SetField(str('value'), value)
            layer.CreateFeature(feature)
            linestring.Destroy()
            feature.Destroy()

if __name__ == "__main__":
    main()
