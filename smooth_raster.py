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
import logging

import numpy as np
from scipy.ndimage import gaussian_filter


try:
    from osgeo import ogr
#    from osgeo.gdalconst import GDT_Float32, GDT_Int16, GA_ReadOnly
    __gdal_loaded__ = True
except:
    __gdal_loaded__ = False


#Docstrings for the option parser
usage = "usage: %prog [options] "
version = "%prog 1.0"


CONSTANT_FIELDS = ['IX', 'IY', 'ID', 'id_1']
#SMOOTH_FIELDS = ['dMean012', 'dMean013', 'dMean014', 'dMean022', 'dMean023', 'dMean024', 'dMean032', 'dMean033', 'dMean034', 'dMean042', 'dMean043', 'dMean044', 'dMean052', 'dMean053', 'dMean054']
#SMOOTH_FIELDS = ['dMean022', 'dMean023', 'dMean024', 'dMean032', 'dMean033', 'dMean034', 'dMean042', 'dMean043', 'dMean044', 'dMean052', 'dMean053', 'dMean054']
#SMOOTH_FIELDS = ['dMean012', 'dMean013', 'dMean014']
SMOOTH_STDDEV = 5


def smooth_rasters(rasters):
    """Smooth rasters using  a filter."""
    for field_index in range(len(SMOOTH_FIELDS)):
        rasters[:, :, field_index] = gaussian_filter(
            rasters[:, :, field_index],
            sigma=5
        )


def get_cellsize(feature):
    poly = feature.GetGeometryRef()
    ring = poly.GetGeometryRef(0)
    xmin = 9e15
    ymin = 9e15
    xmax = -9e15
    ymax = -9e15
    for point_index in range(ring.GetPointCount()):
        x = ring.GetX(point_index)
        xmin = min(x, xmin)
        xmax = max(x, xmax)
        y = ring.GetY(point_index)
        ymin = min(y, ymin)
        ymax = max(y, ymax)
    return (xmax - xmin, ymax - ymin)


def get_raster_pos(x1, y1, dx, dy, nx, ny, x, y):
    """gets the row and column index for given coordinates."""
    col = int((x - x1) / dx)
    row = ny - int(np.ceil((y - y1) / dy))
    if row == ny:
        row -= 1
    if col == nx:
        col -= 1
    return (row, col)


def main():
    #-----------Setting up and unsing option parser-----------------------
    parser = OptionParser(usage=usage, version=version)

    parser.add_option(
        '-v', '--verbose',
        action='store_true', dest='verbose',
        help='produce verbose output'
    )

    parser.add_option(
        '-i',
        action='store',
        dest='infile',
        help='Input shape file'
    )

    parser.add_option(
        '-o',
        action='store',
        dest='outfile',
        help='Output shape file'
    )

    (options, args) = parser.parse_args()

    if options.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logging.basicConfig(level=log_level)
    log = logging.getLogger(__name__)

    #Assure that gdal is present
    if not __gdal_loaded__:
        log.error("Function readGDAL needs GDAL with python bindings")
        sys.exit(1)

    #Create and inititialize output vector data source
    shape_driver = ogr.GetDriverByName('ESRI Shapefile')
    if path.exists(options.outfile):
        shape_driver.DeleteDataSource(options.outfile)
    out_shape = shape_driver.CreateDataSource(options.outfile)
    if out_shape is None:
        log.error("Could not open output shapefile %s" % options.outfile)
        sys.exit(1)

    out_layer = out_shape.CreateLayer(
        options.outfile,
        geom_type=ogr.wkbPolygon
    )

    in_shape = shape_driver.Open(options.infile, update=0)
    if in_shape is None:
        log.error("Could not open data source: " + options.input)
        sys.exit(1)

    in_layer = in_shape.GetLayer()
    layer_defn = in_layer.GetLayerDefn()

    # get field definitions from input shape
    field_definitions = [
        layer_defn.GetFieldDefn(field_ind) for
        field_ind in range(layer_defn.GetFieldCount())
    ]

    # Get all field names, exept pre-defines constant fields
    # Removes first 4 values in list!
    field_names = [field_def.GetName() for field_def in field_definitions]
    for emil in range(0,4):  # Loops 0:3
        field_names.pop(0) # Remove first item in list 
    nr_fields = len(field_names)
    print "Found "+str(nr_fields)+" fields to smooth: "
    print field_names
    global SMOOTH_FIELDS
    SMOOTH_FIELDS = field_names   

    # create fields on output shape
    for field_def in field_definitions:
        out_layer.CreateField(field_def)

    feature = in_layer.GetNextFeature()
    dx, dy = get_cellsize(feature)
    x1, y1, x2, y2 = in_layer.GetExtent()
    nx = int((x2 - x1) / dx)
    ny = int((y2 - y1) / dy)

    rasters = np.zeros((ny, nx, len(SMOOTH_FIELDS)))

    while feature:
        polygon = feature.GetGeometryRef()
        # ring = polygon.GetGeometryRef(0)
        centroid = polygon.Centroid()
        x = centroid.GetX()
        y = centroid.GetY()
        row, col = get_raster_pos(x1, y1, dx, dy, nx, ny, x, y)
        for field_ind, field_name in enumerate(SMOOTH_FIELDS):
            val = float(feature.GetFieldAsString(field_name.encode('ascii')))
            rasters[row, col, field_ind] = val
        feature = in_layer.GetNextFeature()

    smooth_rasters(rasters)

    in_layer.ResetReading()
    for feature in in_layer:
        polygon = feature.GetGeometryRef()
        #ring = polygon.GetGeometryRef(0)
        centroid = polygon.Centroid()
        x = centroid.GetX()
        y = centroid.GetY()
        row, col = get_raster_pos(x1, y1, dx, dy, nx, ny, x, y)

        for field_ind, field_name in enumerate(SMOOTH_FIELDS):
            val = rasters[row, col, field_ind]
            feature.SetField(field_name.encode('ascii'), val)

        if out_layer.CreateFeature(feature) != 0:
            log.error("Failed to create feature in shapefile.\n")
            sys.exit(1)

        # feature = None
    # out_layer = None
    # out_shape = None
    # in_shape = None


if __name__ == '__main__':
    main()
