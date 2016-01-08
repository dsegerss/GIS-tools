#!/usr/bin/env python
'''
Created on 4 okt 2010
@author: David Segersson
'''

import sys
from os import path
try:
    from osgeo import ogr
except:
    import ogr


integerCoordinates = True


def main():

    outfilename = path.abspath(sys.argv[2])
    driver = ogr.GetDriverByName('ESRI Shapefile')

    if path.exists(outfilename):
        driver.DeleteDataSource(outfilename)
    datasource = driver.CreateDataSource(outfilename)

    if datasource is None:
        raise OSError("Could not open data source: " + outfilename)

    layer = datasource.CreateLayer('roads', geom_type=ogr.wkbLineString)

    column_names = sys.argv[2].split(',')
    fieldDefn = ogr.FieldDefn(b'FLOW', ogr.OFTInteger)
    layer.CreateField(fieldDefn)
    fieldDefn = ogr.FieldDefn(b'NAME', ogr.OFTString)
    layer.CreateField(fieldDefn)
    fieldDefn = ogr.FieldDefn(b'INFO', ogr.OFTInteger)
    layer.CreateField(fieldDefn)

    featureDefn = layer.GetLayerDefn()

    lines = open(sys.argv[1], 'r').readlines()

    for i, line in enumerate(lines[1:]):
        print (float(i) / float(len(lines) - 1) * 100)
        feature = ogr.Feature(featureDefn)
        roadid, name, flow = line.split('\t')[:3]
        flow = int(flow)
        coords = line.split('\t')[3:]
        nodei = 0
        road = ogr.Geometry(ogr.wkbLineString)
        while 1:
            try:
                x = int(coords[nodei])
                y = int(coords[nodei + 1])
                road.AddPoint(x, y)
                nodei += 2
            except:
                feature.SetGeometry(road)
                break

#         for colId in inputTable.listIds():
        feature.SetField(b'FLOW', flow)
        feature.SetField(b'NAME', name)
        feature.SetField(b'INFO', roadid)
        layer.CreateFeature(feature)
        road.Destroy()
        feature.Destroy()
    datasource.Destroy()


if __name__ == "__main__":
    main()

