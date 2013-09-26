#!/usr/bin/env python
# -*- coding: latin-1 -*-
import os, sys, subprocess, urllib, time
from os import path
from optparse import OptionParser

import pdb

#Docstrings for the option parser
usage = "usage: %prog --ll '(lon lat') --ur '(lon lat')"
version="%prog 1.0"


def main():
    parser=OptionParser(usage= usage, version=version)
    
    parser.add_option("--ll",action="store",dest="ll",
                      help="Lower left corner of extraction area (whole degrees)")
    
    parser.add_option("--ur",action="store",dest="ur",
                      help="Upper right corner of extraction area (whole degrees)")

    parser.add_option("--url", action="store",dest="url",
                      help="Base URL for SRTM-files",
                      default="http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/South_America")
    parser.add_option("-o","--outputDir",action="store",dest="outputDir",
                      help="Output directory")

    parser.add_option("-f","--force",action="store_true",dest="force",
                      help="To overwrite already downloaded tiles")

     
    (options, args) = parser.parse_args()

    if options.outputDir is None:
        parser.error("Has to specify output directory")

    outputDir=path.abspath(options.outputDir)
    if not path.exists(outputDir):
        parser.error("Outut dir does not exist")
        
    if options.ll is None:
        parser.error("Has to specify lower left corner")
    if options.ur is None:
        parser.error("Has to specify upper right corner")
        
    try:
        ll=eval(options.ll)
        ll=map(int,ll)
        ur=eval(options.ur)
        ur=map(int,ur)
    except:
        print "Could not read coordinate pairs, make sure they are in format '(lon lat)'"
        sys.exit()

    url=options.url

    ntiles=len(range(ll[0],ur[0]))*len(range(ll[1],ur[1]))
    
    itile=1
    for lon in range(ll[0],ur[0]):
        for lat in range(ll[1],ur[1]):
            tile=""
            if lat<0:
                tile+="S%02i" %(abs(lat))
            else:
                tile+="N%02i" %(lat)
            if lon<0:
                tile+="W%03i" %(abs(lon))
            else:
                tile+="E%03i" %(lon)
            tif=tile+".tif"
            tile+=".hgt.zip"
            tilePath=path.join(outputDir,tile)
            tifPath=path.join(outputDir,tif)
            if path.exists(tifPath):
                print "Keeping tile: %s" %tile
                continue
            sys.stdout.write("Downloading tile %i out of %i: %s..." %(itile,ntiles,tile))
            #subprocess.call("wget "+path.join(url,tile))
            urllib.urlcleanup()
#            print "url is: %s" %
            pdb.set_trace()
            urllib.urlretrieve(path.join(url,tile),tilePath)
            time.sleep(5)
            sys.stdout.write("done, converting to geotiff...")
            os.chdir(outputDir)
            print "calling: "+"srtm_generate_hdr.sh "+tile
            returnCode = subprocess.call(["srtm_generate_hdr.sh",tile])
            if returnCode!=0:

                print "Could not download tile from url: %s" %url
                sys.exit(1)
                
            sys.stdout.write("done\n")
            itile+=1

    
if __name__=="__main__":
    main()
