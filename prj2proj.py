#!/usr/bin/env python
"""Export ESRI projection definition to proj4-definition"""

from osgeo import osr
import sys

def main(prj_file):
    prj_text = open(prj_file, 'r').read()
    srs = osr.SpatialReference()
    if srs.ImportFromWkt(prj_text):
        raise ValueError("Error importing PRJ information from: %s" % prj_file)

    print srs.ExportToProj4()
         
    if __name__=="__main__":
        main(*sys.argv[1])
