#!/bin/bash

###################################################################
### C. Asker                                         2014-03-13 ###
### Stupid script to get elevation contours based on SRTM data. ###
### Files are downloaded from http://www.opendem.info           ###
### All files are merged into one shapefile                     ###
### input should be in WHOLE degrees!                           ###
###################################################################
### FIXME: user-friendly CLI interface
### FIXME: what happens for negative lat/lons...?

WEST=$1
EAST=$2
NORTH=$3
SOUTH=$4

baseurl="http://opendemdata.info/data/srtm_contour/"
#http://opendemdata.info/data/srtm_contour/N48E034.zip
outfile="contour.shp"

dir=`pwd`

mkdir -p contour.tmp
cd contour.tmp

### Loop over all "cells"
for lat in `seq --format="%02.0f" $4 $3`; do

    for lon in `seq --format="%03.0f" $1 $2`; do
        echo ${baseurl}N${lat}E${lon}.zip
        wget -q ${baseurl}N${lat}E${lon}.zip && unzip -q N${lat}E${lon}.zip && mv N${lat}E${lon}/* .
    done

done

### Create empty shapefile
echo "\nCreating shapefiles"

first=1
for fil in N*E*.shp; do
    echo $fil
    if [ $first == 1 ] ; then 
        ogr2ogr -a_srs EPSG:4326 $outfile $fil
        first=0
    elif [ $first != 1 ] ; then
        echo ogr2ogr -update -append $outfile $fil -f "esri shapefile" -nln contour
        ogr2ogr -update -append $outfile $fil -f "esri shapefile" -nln contour
    fi
done
rm -rf N*E*

cd $dir
