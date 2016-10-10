#!/usr/bin/env python
# coding=utf-8

"""Tools for post processing of CFD results."""

from __future__ import unicode_literals
from __future__ import division

from os import path
import glob
import sys
import argparse
import logging

import numpy as np
from osgeo import ogr

log = logging.getLogger(__name__)


def create_ogr_geometry(points, indices):
    """Build ogr polygon from points and polygon indices."""

    poly_points = points[indices, :]
    poly = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for p in poly_points:
        ring.AddPoint(*p)
    ring.CloseRings()
    poly.AddGeometry(ring)
    return poly


def write_features(
        layer, points, polygons,
        field_data, fieldnames, layer_definition):
    """Write polygon features to layer."""

    for i in range(len(polygons)):
        poly = create_ogr_geometry(points, polygons[i])
        feature = ogr.Feature(layer_definition)
        feature.SetGeometry(poly)
        for j, fieldname in enumerate(fieldnames):
            polygon_average = field_data[j][polygons[i]].mean()
            feature.SetField(
                fieldnames[j].encode('ascii'),
                polygon_average
            )
        layer.CreateFeature(feature)


def read_vtk_geometry(vtkfile):
    """Read points and polygon point indices from vtk-file object."""

    log.debug('Reading VTK-file')
    line = vtkfile.readline()
    while 'POINTS' not in line:
        line = vtkfile.readline()
    npoints = int(line.split()[1])
    points = np.fromfile(vtkfile, count=npoints * 3, sep=' ')
    points = np.reshape(points, (npoints, 3))
    log.debug('Found %i points' % npoints)

    while 'POLYGONS' not in line:
        line = vtkfile.readline()
    line = vtkfile.readline()

    polygons = []
    while 'POINT_DATA' not in line:
        point_indices = np.array(map(int, line.split()[1:]))
        polygons.append(point_indices)
        line = vtkfile.readline()
    log.debug('Found %i polygons' % len(polygons))
    return points, polygons
    

def read_vtk_point_attribute(vtkfile):
    """Convert vtk to ogr supported vector format."""
    line = vtkfile.readline()
    while 'FIELD attributes' not in line:
        line = vtkfile.readline()
    # currently only handles one float attribute
    # nattributes = int(line.split()[-1])
    line = vtkfile.readline()
    fieldname, _, nvals, type = line.split()
    nvals = int(nvals)
    values = np.fromfile(vtkfile, count=nvals, sep=' ')
    log.debug('Found %i field values' % values.size)
    return (fieldname, values)


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '-i',
        dest='infiles', action='store', metavar='INFILES', nargs='*',
        help='Input vtk-files (geometry only read from first file)'
    )

    parser.add_argument(
        '-o',
        dest='outfile', action='store', metavar='OUTFILE',
        help='Output file'
    )

    parser.add_argument(
        '--format',
        dest='format', action='store', metavar='format',
        default='ESRI Shapefile',
        help='Output file'
    )

    parser.add_argument(
        '--overwrite',
        dest='overwrite', action='store_true',
        help='Overwrite output file if it already exists'
    )

    parser.add_argument(
        '--move', type=float,
        dest='move', action='store', nargs=3, metavar=('X', 'Y', 'Z'),
        help='Move coordinates in vtk-file'
    )

    parser.add_argument(
        '-v',
        action=VerboseAction, dest='loglevel', default=logging.WARNING,
        help='increase verbosity in terminal'
    )

    args = parser.parse_args()

    driver = ogr.GetDriverByName(args.format.encode('ascii'))

    fieldnames = []
    field_data = []
    with open(args.infiles[0], 'r') as vtkfile:
        points, polygons = read_vtk_geometry(vtkfile)
        fieldname, values = read_vtk_point_attribute(vtkfile)
        fieldnames.append(fieldname)
        field_data.append(values)

    if args.move is not None:
        points = points + np.array(args.move)

    for infile in args.infiles[1:]:
        with open(infile, 'r') as vtkfile:
            fieldname, values = read_vtk_point_attribute(vtkfile)
        if values.size != points.shape[0]:
            log.error(
                'Inconsistency between number of points and number ' +
                'of attribute values in %s' % infile
            )
            sys.exit(1)
        fieldnames.append(fieldname)
        field_data.append(values)

    # create file
    if path.exists(args.outfile):
        if args.overwrite:
            driver.DeleteDataSource(args.outfile)
        else:
            raise IOError('File %s already exist' % args.outfile)
    out = driver.CreateDataSource(args.outfile)
    layer = out.CreateLayer(
        args.outfile,
        geom_type=ogr.wkbPolygon
    )

    # define layer
    for i in range(len(fieldnames)):
        fieldnames[i] = fieldnames[i][:10]
        field_definition = ogr.FieldDefn(
            fieldnames[i].encode('ascii'),
            ogr.OFTReal
        )
        layer.CreateField(field_definition)
    layer_definition = layer.GetLayerDefn()

    write_features(
        layer, points, polygons, field_data, fieldnames, layer_definition
    )
    out.Destroy()
    log.debug('finished')


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
