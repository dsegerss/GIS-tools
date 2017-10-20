#!/usr/bin/env python
# -*- coding: latin-1 -*-
"""Converts rasters supported by GDAL to FoUl-GRIB"""

from __future__ import division
from __future__ import unicode_literals

import datetime
import sys
import re
import argparse
import logging
import csv
from os import path
from collections import OrderedDict


import numpy as np

from osgeo import osr
from osgeo import gdal
from osgeo.gdalconst import (
    GA_ReadOnly
)

import pygrib as pg

log = logging.getLogger(__name__)

NODATA_VALUE = -9999

SUBSTANCE2PAR = {
    'NOx': 'NOX',
    'BC': 'EC',
    'SOx': 'SO2',
    'PM10': 'PM10',
    'PM2.5': 'PM2.5',
    'NH3': 'NH3',
    'CO': 'CO',
    'benzo_a_pyrene': 'BaP',
    'NMVOC': 'NMVOC',
    'PAH-4': 'PAH_4',
    'CH4': 'CH4',
    'OC': 'OC'
}

SECTORS = OrderedDict(
    [
        ('01', {'long': 'Combustion in energy and transformation industries',
                'short': 'Energy', 'lev': 1}),
        ('02', {'long': 'Non-inudustrial combustion',
                'short': 'Non-ind combustion', 'lev': 2}),
        ('0201', {'long': 'Commercial and institutional plants',
                  'short': 'Commercial', 'lev': 201}),
        ('0202', {'long': 'Residential plants',
                  'short': 'Residential', 'lev': 202}),
        ('0203', {'long': 'Agriculture, forestry and aquaculture',
                  'short': 'Other', 'lev': 203}),
        ('03', {'long': 'Combustion in manufacturing industry',
                'short': 'Industry', 'lev': 3}),
        ('04', {'long': 'Production processes', 'short': 'Process', 'lev': 4}),
        ('05', {'long': 'Extraction and distribution of fossil fuels and geothermal energy',
                'short': 'Extraction', 'lev': 5}),
        ('06', {'long': 'Solvent and other product use',
                'short': 'Solvent', 'lev': 6}),
        ('07', {'long': 'Road import pdb; pdb.set_trace()ansport',
                'short': 'Road transport', 'lev': 7}),
        ('0701-0705', {'long': 'Road transport - exhaust',
                       'short': 'Road exhaust', 'lev': 7015}),
        ('0706-0708', {'long': 'Road transport - non-exhaust',
                       'short': 'Road non-exhaust', 'lev': 768}),
        ('08', {'long': 'Other mobile sources and machinery',
                'short': 'Off-road', 'lev': 8}),
        ('0801', {'long': 'Military', 'short': 'Military', 'lev': 801}),
        ('0802', {'long': 'Railways', 'short': 'Railways', 'lev': 802}),
        ('080501-080502', {'long': 'Aviation - LTO',
                           'short': 'Aviation', 'lev': 8512}),
        ('0806-0807', {'long': 'Agriculture and forestry',
                       'short': 'Agriculture and forestry', 'lev': 867}),
        ('0808', {'long': 'Industry', 'short': 'Industry', 'lev': 808}),
        ('0809', {'long': 'Houshold and gardening',
                  'short': 'Household', 'lev': 809}),
        ('0810', {'long': 'Other off-road',
                  'short': 'Other off-road', 'lev': 810}),
        ('09', {'long': 'Waste treatment and disposal', 'short': 'Waste', 'lev': 9}),
        ('10', {'long': 'Agriculture', 'short': 'Agriculture', 'lev': 10}),
        ('shipping', {'long': 'Shipping', 'short': 'Shipping', 'lev': 11})
    ]
)


# def xy2lonlat(proj4, x, y):
#     source_srs = osr.SpatialReference()
#     source_srs.ImportFromProj4(proj4.encode('ascii'))
#     target_srs = osr.SpatialReference()
#     target_srs.ImportFromEPSG(4326)  # wgs84 lat lon
#     coord_trans = osr.CoordinateTransformation(source_srs, target_srs)
#     lon, lat, _ = coord_trans.TransformPoint(x, y)
#     return (lon, lat)


def add_point_source(data, extent, x, y, value):
    ny, nx = data.shape
    x1, y1, x2, y2 = extent
    dx = (x2 - x1) / nx
    dy = (y2 - y1) / ny

    if((x < x1 or y < y1) or (x > x2 or y > y2)):
        log.debug(
            'Coord %f %f is outside %f %f, %f %f' % (
                x, y, x1, y1, x2, y2
            )
        )
    else:
        col = int((x - x1) / dx)
        row = int((y - y2) / (-1 * dy))  # dy must be positive
        if row == ny:
            row -= 1
        if col == nx:
            col -= 1
        data[row, col] += value


def create_grib(extent, nx, ny, proj4):

    x1, y1, x2, y2 = extent
    dx = (x2 - x1) / nx
    dy = (y2 - y1) / ny
    lonw, lats = pg.proj4_inv(x1, y1, proj4)

    proj = {}
    proj['proj4'] = proj4
    proj["lonw"] = lonw
    proj["lats"] = lats
    proj["dlon"] = dx  # may be in meters
    proj["dlat"] = dy  # may be in meters
    proj["nx"] = nx
    proj["ny"] = ny
    proj["projid"] = 99

    g = pg.seed()
    g.setprojection(proj)

    return g


def extract_by_bbox(extent, nx, ny, bbox, strict_bounds=True):

    xll, yll, xur, yur = extent
    dx = (xur - xll) / nx
    dy = (yur - yll) / ny

    x1, y1, x2, y2 = bbox

    if strict_bounds:
        if x1 < xll or y1 < yll or x2 > xur or y2 > yur:
            raise ValueError(
                'Trying to extract outside grid boundaries\n' +
                'input extent %f %f %f %f\n' % (xll, yll, xur, yur) +
                'requested extent %f %f %f %f' % (x1, y1, x2, y2)
            )
    else:
        if x2 < xll or y2 < yll or x1 > xur or y1 > yur:
            raise ValueError(
                'Trying to read data completely outside of grid boundaries\n' +
                'input extent %f %f %f %f\n' % (xll, yll, xur, yur) +
                'requested extent %f %f %f %f' % (x1, y1, x2, y2)
            )

    colmin = int(round((x1 - xll) / dx))
    colmax = int(round((x2 - xll) / dx))
    rowmin = int(round((yur - y2) / dy))
    rowmax = int(round((yur - y1) / dy))

    colmin_from_file = max(colmin, 0)
    colmax_from_file = min(colmax, nx - 1)
    rowmin_from_file = max(rowmin, 0)
    rowmax_from_file = min(rowmax, ny - 1)

    nrows = rowmax - rowmin
    ncols = colmax - colmin

    nrows_from_file = rowmax_from_file - rowmin_from_file
    ncols_from_file = colmax_from_file - colmin_from_file

    # size of padding following numpy.pad doc
    # [(before_y, after_y), (before_x, after_x)]

    padding = (
        (abs(min(0, rowmin)), rowmax - rowmax_from_file),
        (abs(min(0, colmin)), colmax - colmax_from_file)
    )

    extent = (
        xll + colmin * dx,
        yur - rowmax * dy,
        xll + colmax * dx,
        yur - rowmin * dy
    )

    return (
        extent,
        rowmin_from_file,
        colmin_from_file,
        nrows_from_file,
        ncols_from_file,
        padding
    )


def read_raster(filename, band=1, bbox=None,
                strict_bounds=False, pad=True, padding_values=None):
    """Read raster,
    bbox (x1, y1, x2, y2) specifies window to read
    strict_bounds=True means that bbox must be within raster bounds
    or ValueError will be raised
    pad means that padding with nodata will be carried out for
    parts of bbox outside of raster
    padding_values is a tuple of values to add in
    1st and second dimension, default is (nodata, nodata)
    """

    ds = gdal.Open(filename, GA_ReadOnly)
    if ds is None:
        raise ValueError('Could not open %s' % filename)

    meta = {}
    # dimensions
    ny = ds.RasterYSize
    nx = ds.RasterXSize

    # Info used for georeferencing
    geoTransform = ds.GetGeoTransform()
    xul = geoTransform[0]  # top left x
    dx = geoTransform[1]  # w-e pixel resolution
    rot1 = geoTransform[2]  # rotation, 0 if image is "north up"
    yul = geoTransform[3]  # top left y
    rot2 = geoTransform[4]  # rotation, 0 if image is "north up"
    dy = geoTransform[5]  # n-s pixel resolution

    # calculate lower left corner
    xll = xul
    yll = yul + ny * dy  # cellsizeY should be a negative value
    xur = xll + nx * dx
    yur = yll - ny * dy
    extent = (xll, yll, xur, yur)
    extent = map(lambda v: round(v, 5), extent)
    if band <= ds.RasterCount:
        band = ds.GetRasterBand(band)
    else:
        raise ValueError('Band %i not available' % band)

    if bbox is not None:
        extent, rowmin, colmin, ny, nx, padding = extract_by_bbox(
            extent, nx, ny, bbox, strict_bounds=strict_bounds
        )
    else:
        rowmin = colmin = 0
        padding = [(0, 0), (0, 0)]

    data = band.ReadAsArray(
        xoff=colmin, yoff=rowmin,
        win_xsize=nx,
        win_ysize=ny
    )
    nodata = band.GetNoDataValue() or NODATA_VALUE

    if pad:
        data = np.pad(
            data, padding, mode=b'constant',
            constant_values=padding_values or (nodata, nodata)
        )
    meta = {
        'extent': extent,
        'proj': ds.GetProjection(),
        'nodata': nodata,
        'padding': padding
    }

    return data, meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        action="store", dest="infile",
        help="Input raster file"
    )

    parser.add_argument(
        action="store", dest="outfile",
        help="Output GRIB file"
    )

    parser.add_argument(
        '--epsg',
        action='store', type=int,
        dest='epsg',
        help='EPSG code of projection used \
        (default is to read projection from dataset)'
    )

    parser.add_argument(
        "--par",
        action="store", dest="par",
        help="Par in output grib"
    )

    parser.add_argument(
        "--lev", type=int,
        action="store", dest="lev", default=0,
        help="Lev in target grib file"
    )

    parser.add_argument(
        '--sort',
        action='store', dest="sort", default='Q(TON)',
        help='Sort in target grib file'
    )

    parser.add_argument(
        '--surf',
        action='store', dest="surf", default='ALL',
        help='Surface in target grib file'
    )

    parser.add_argument(
        '--time',
        action='store', dest="time", default='YEARACCUM',
        help='Time repr. in target grib file'
    )

    parser.add_argument(
        "--overwrite",
        action="store_true", dest="overwrite",
        help="Overwrite grib file"
    )

    parser.add_argument(
        '--timestamp',
        action='store', dest='timestamp',
        type=lambda x: datetime.datetime.strptime(x, '%y%m%d%H'),
        help='Timestamp as YYMMDDHH'
    )

    parser.add_argument(
        '--name-pattern',
        action="store", dest='name_pattern',
        help='Specify reg-exp to extract named groups for \
        <substance> and <sector>'
    )

    parser.add_argument(
        '--extent',
        action='store', dest='extent', nargs=4,
        metavar=('X1', 'Y1', 'X2', 'Y2'),
        help='Extent of target GRIB file (default is extent of raster)'
    )

    parser.add_argument(
        '--pointsources',
        action='store', dest='pointsources',
        metavar='pointsource-table',
        help='csv-file with point sources to include in grib'
    )

    parser.add_argument(
        '-v',
        action=VerboseAction, dest='loglevel', default=logging.WARNING,
        help='increase verbosity in terminal'
    )

    args = parser.parse_args()

    if args.name_pattern is None:
        par = args.par
        lev = args.lev
    else:
        match = re.compile(args.name_pattern).search(args.infile)
        if match is None:
            log.error('Specified name-pattern provides no match')
            sys.exit(1)
        groups = match.groupdict()
        try:
            substance = groups.get('substance', None)
            par = SUBSTANCE2PAR[substance]
        except KeyError:
            log.error('No par-mapping for substance %s, skipping' % substance)
            sys.exit(1)
        try:
            sector = groups.get('sector', None)
            lev = SECTORS[sector]['lev']
        except KeyError:
            log.error('No lev-mapping for sector %s, skipping' % sector)
            sys.exit(1)

    if args.timestamp is not None:
        date = pg.setdate(
            args.timestamp.year,
            args.timestamp.month,
            args.timestamp.day,
            args.timestamp.minute,
            len=12
        )
    else:
        date = pg.setdate(1971, 01, 01, 00, len=12)

    outfile = pg.getfilename(args.outfile, "GRIB1", date)

    if not path.exists(outfile):
        grbs = pg.open(args.outfile, "w", date=date, suffix='GRIB1')
    else:
        # grb_old = get_field_if_exist(args.outfile, date)
        grbs = pg.open(args.outfile, "a", suffix='GRIB1', date=date)

    # register all of the raster drivers
    gdal.AllRegister()

    data, meta = read_raster(
        args.infile, bbox=args.extent, padding_values=NODATA_VALUE
    )
    data = np.where(data == meta['nodata'], 0, data)
    nx = data.shape[1]
    ny = data.shape[0]
    extent = meta['extent']
    proj4 = meta['proj']

    if args.epsg is not None:
        source_srs = osr.SpatialReference()
        source_srs.ImportFromEPSG(args.epsg)
        proj4 = source_srs.ExportToProj4()
    elif proj4 is None:
        log.error(
            'No projection info in dataset, please specify using --epsg'
        )
        sys.exit(1)

    if args.pointsources is not None:
        pstable = open(args.pointsources, 'r')
        psreader = csv.DictReader(pstable,  delimiter=b';')
        ps_emis_tot = 0
        for rownr, row in enumerate(psreader):
            if args.timestamp is not None:
                year = int(row['Year'])
                if year != args.timestamp.year:
                    continue

            code = row['SNAP_sector_code']
            if code not in SECTORS:
                log.error('Sector code of point source %s not defined' % code)
                sys.exit(1)
            if sector == code:
                emis = row.get(substance, 0)
                if emis != 0 and emis.strip() != '':
                    try:
                        emis = float(emis)
                    except ValueError:
                        log.error(
                            'invalid point-source emission value ' +
                            '%s on row %i in %s' %
                            (str(emis), rownr, args.pointsources)
                        )
                        sys.exit(1)

                    ps_x = float(row['X_ETRS89'])
                    ps_y = float(row['Y_ETRS89'])
                    add_point_source(data, extent, ps_x, ps_y, emis)
                    ps_emis_tot += emis
        if ps_emis_tot > 0:
            log.debug(
                'Added %f ton/year %s from point-sources' % (
                    ps_emis_tot, substance
                )
            )

    grb = create_grib(extent, nx, ny, proj4)

    # fortran data ordering used in pygrib
    grb.setvalues(np.flipud(data).T)

    grb.setpar(
        par,
        lev=lev,
        sort=args.sort,
        surf=args.surf,
        time=args.time
    )

    grbs.put(grb)
    log.info(
        'Wrote par %s, lev %i, sort %s surf %s to %s' % (
            grb.par, grb.lev, grb.sort, grb.surf, args.outfile)
    )

    grbs.close()


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


if __name__ == "__main__":
    main()
