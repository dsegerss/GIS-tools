#!/usr/bin/env python
'''
Created on 4 okt 2010
@author: David Segersson
'''
import sys,os
from math import floor
from os import path
try:
    from osgeo import ogr
except:
    import ogr
from pyAirviro_dev.other import dataTable
import pdb

integerCoordinates=True

def main():

    inputTablePath=path.abspath(sys.argv[1])
    inputTable=dataTable.DataTable()
    inputTable.read(inputTablePath)
    
        
    outputPath=path.abspath(sys.argv[2])

    driver = ogr.GetDriverByName('ESRI Shapefile')

    if path.exists(outputPath):
        driver.DeleteDataSource(outputPath)
    dataSource = driver.CreateDataSource(outputPath)
    
    if dataSource is None:
        raise OSError("Could not open data source: " + outputPath)


    layer = dataSource.CreateLayer('roads', geom_type=ogr.wkbLineString)

    # add an id field to the output
    for colId in inputTable.listIds():
        fieldDefn = ogr.FieldDefn(str(colId),ogr.OFTString)
        layer.CreateField(fieldDefn)
    
    feature = layer.GetNextFeature() #get first polygon in shape file    
    



    featureDefn = layer.GetLayerDefn()
    t=inputTable

    for i,road in enumerate(inputTable.data):
        print float(i)/float(inputTable.nrows())*100
        
        feature = ogr.Feature(featureDefn)


        x1=int(float(road[t.colIndex("X1")]))
        x2=int(float(road[t.colIndex("X2")]))
        y1=int(float(road[t.colIndex("Y1")]))
        y2=int(float(road[t.colIndex("Y2")]))

        seg=ogr.Geometry(ogr.wkbLineString)
        seg.AddPoint(x1,y1)
        seg.AddPoint(x2,y2)
        feature.SetGeometry(seg)
        
        for colId in inputTable.listIds():
            feature.SetField(str(colId),road[t.colIndex(colId)]) 
        
        layer.CreateFeature(feature)
        
        seg.Destroy()
        feature.Destroy()
        
    dataSource.Destroy()

if __name__=="__main__":
    main()

