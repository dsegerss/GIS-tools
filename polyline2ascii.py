#!/usr/bin/env python
'''
Created on 4 okt 2010
@author: David Segersson
'''
from math import floor
try:
    from osgeo import ogr
except:
    import ogr


shapeFile="/usr/airviro/data/sverige/prod/projects/karlstad_2011/karlstad_rt90.shp"
ouputPath="/usr/airviro/data/sverige/prod/projects/karlstad_2011/karlstad.asc"

#Type is "numeric" or "text"
fields=[
        {"name":"klassning","type":"numeric"}
        ]
integerCoordinates=True

def main():
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource =  driver.Open(shapeFile, update=0)
    if dataSource is None:
        raise OSError("Could not open data source: " + shapeFile)
    outfile=open(ouputPath,"w")
    layer = dataSource.GetLayer() #Create a layer for the shapefile (fil_mall)        
    feature = layer.GetNextFeature() #get first polygon in shape file    
    
    print "Reading lines..."   
    

    while feature:
        
        for f in fields[:-1]:
            outfile.write(f["name"]+"\t")
        outfile.write(fields[-1]["name"]+"\n") 

        for f in fields:
            if f["type"]=="text":
                value=feature.GetFieldAsString(f["name"])
            if f["type"]=="numeric":
                value=feature.GetFieldAsDouble(f["name"])
            outfile.write(str(value)+"\t")
        outfile.write("\n")
        #Reading outer ring and stores as a list of points
        centreline = feature.GetGeometryRef()
        outerList=[]
        #print "pointcount: ", centreline.GetPointCount()
        for i in range(centreline.GetPointCount()):
            x=centreline.GetX(i)
            y=centreline.GetY(i)
            #z=centreline.GetZ(i)-offset[2]
            outfile.write("X%i   %i Y%i  %i\n" % (i,floor(x),i,floor(y)))
        outfile.write("\n")
        feature.Destroy()
        feature = layer.GetNextFeature()            
    outfile.close()
    dataSource.Destroy()
    print "finished"

if __name__=="__main__":
    main()
