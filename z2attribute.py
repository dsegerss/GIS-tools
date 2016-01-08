#!/usr/bin/env python
'''
Created on 4 okt 2010
@author: David Segersson
'''
import sys
try:
    from osgeo import ogr
except:
    print('Could not import ogr from osgeo')
    sys.exit(1)


def main():
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(sys.argv[1], update=1)
    if dataSource is None:
        raise OSError("Could not open data source: " + sys.argv[1])
    layer = dataSource.GetLayer()
    layer_defn = layer.GetLayerDefn()
    field_names = [
        layer_defn.GetFieldDefn(i).GetName()
        for i in range(layer_defn.GetFieldCount())
    ]
    if 'ELEVATION' not in field_names:
        new_field = ogr.FieldDefn(b'ELEVATION', ogr.OFTReal)
        layer.CreateField(new_field)

    feature = layer.GetNextFeature()
    while feature:
        geom = feature.GetGeometryRef()
        if geom is not None and geom.GetPointCount() > 0:
            z = geom.GetZ(0)
            feature.SetField(b'ELEVATION', z)
            layer.SetFeature(feature)

        feature.Destroy()
        feature = layer.GetNextFeature()
    dataSource.Destroy()

if __name__ == '__main__':
    main()
