#!/usr/bin/env python
# -*- coding: latin-1 -*-

"""
Created on 27 Oct 2010
author: David Segersson

Script exports ESRI ascii rasters to grib

"""
import os, sys, glob,re
from os import path
from optparse import OptionParser

import numpy as np
import pygrib as pg
import raster

import pdb

#Docstrings for the option parser
usage = "usage: %prog <rasterFile> <gribFile>  -[-s 'substances' -c 'sectors']"
version="%prog 1.0"

def main():
    #-----------Setting up and unsing option parser-----------------------
    parser=OptionParser(usage= usage, version=version)
        
    parser.add_option("-l", "--lev",
                      action="store",dest="lev",
                      help="Lev in target grib file")   
 
    parser.add_option("-o","--overwrite",
                      action="store_true",dest="overwrite",
                      help="Overwrite grib file, [default: %default]")

    parser.add_option("--levfilter",
                      action="store",dest="levfilter",
                      help="Regular expression to extract level from\
                      raster path, level should be\
                      given as 1:st group in expression")

#   parser.add_option("-y","--year",
#                      action="store",dest=year,
#                      help="Year for emission grib filename")


    (options, args) = parser.parse_args()


    if len(args)!=2:
        parser.error("Wrong number of arguments")
        
    rastFileName=path.abspath(args[0])
    rast=raster.raster()
    try:
        rast.read(rastFileName)
    except:
        print "Could not read raster file"
        sys.exit(1)

#    if options.year is not None:
#        date=pg.setdate(int(year),01,01,00,len=12)
#    else:
#        date=pg.setdate(1971,01,01,00,len=12)
    
    grbName=path.abspath(args[1])

    #Read template grib
    grbs=pg.open(grbName,"r",suffix='FIX')
    
    templateGrb=grbs.next()
    grbs.close()      
    
    if options.lev is not None:
        templateGrb.lev=int(options.lev)       
    elif options.levfilter is not None:
        pattern=re.compile(options.levfilter)
        match = pattern.match(rastFileName)
        if match is None:
            parser.error("levfilter does not match raster path")
            sys.exit(1)
        templateGrb.lev=int(match.groups()[0])
            
    
    #Add raster data to template
    rast.flip()
    templateGrb["values"]=rast.data.transpose()

    if options.overwrite:
        try:
            os.rename(grbName,grbName+".tmp")
            grbs=pg.open(grbName,"w",suffix="FIX")
        except:
            os.rename(grbName+".tmp",grbName)
            print "Could not overwrite Grib: %s" %gribName
            sys.exit(1)
        os.remove(grbName+".tmp")        
    else:
        grbs=pg.open(grbName,"a",suffix="FIX")        
    grbs.put(templateGrb)
    grbs.close()

    print "Finished!"

if __name__=="__main__":
    main()
