#!/usr/bin/env python
# -*- coding: latin-1 -*-
doc="""
********************************************************************************
Name: setBuildingHeights.py
Created: 22 Mar 2013
Author: David Segersson

Description
--------------------------------------------------------------------------------
Set building height as attribute to 3D line features
Height is read from z-coordinate of lines, and from a raster DEM
Difference between building z and DEM z is calculated for all vertices
Average height difference is added as attribute BUILDHGT to result polyline feature class
"""
#Standard modules
from os import path
import sys,os
from optparse import OptionParser
import pdb
from math import ceil
import shutil
import glob

#pyAirviro-modules
from pyAirviro.other import logger
from  pyAirviro.other.utilities import ProgressBar
from pyAirviro.geo.raster import readGDAL

#help module for this script
from geometry import Segment,perp2D,intersect2D_2Segments

#plot modules
from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as Patches

import numpy as np

#GDAL bindings
try:
    from osgeo import ogr
    from osgeo import gdal
    from osgeo.gdalconst import *
    __gdal_loaded__=True
except:
    __gdal_loaded__=False


#Docstrings for the option parser
usage = "usage: %prog [options] "
version="%prog 1.0"

#Global constants
CELLSIZE=75 #Cellsize of spatial index for building contours
MAXDIST=100  #Maximal distance from road within which buildings are processed
CSDIST=40   #Max distance between cross-sections
HEIGHTCORR=1 #Added height to buildings to account for non-flat roofs
PLOTIND=-1 #Index To plot for debugging
CSIND=1 #cross section index to plot for debugging
#Creating global log
log=None



def splitRoad(inRoadFeature,splitLimit,layerDefn):
    """
    Split road in sections with same direction
    Returns list of features with same field values as original feature
    @param inRoadFeature: road feature to split
    @param splitLimit: threshold for splitting the road
    @param outRoadLayer: layer to create new features in    
    """
    outRoadFeatures=[]
    inRoadGeom = inRoadFeature.GetGeometryRef()
    nPoints=inRoadGeom.GetPointCount()
    if splitLimit is None or nPoints<=2:
        outRoadFeatures.append(ogr.Feature(layerDefn))
        outRoadFeatures[0].SetGeometry(inRoadGeom)
        outRoadFeatures[0].SetFrom(inRoadFeature)
        return outRoadFeatures

    x1=inRoadGeom.GetX(0)
    y1=inRoadGeom.GetY(0)
    z1=inRoadGeom.GetZ(0)

    x2=inRoadGeom.GetX(1)
    y2=inRoadGeom.GetY(1)
    z2=inRoadGeom.GetZ(1)

    #Direction vector of first segment
    Ulast=np.array([x2-x1,y2-y1])
    Ulast=Ulast/np.linalg.norm(Ulast)
    #Convert to polar
    #lastAng=np.degrees(np.arctan2(U[1],U[0]))

    #Create new feature
    outRoadFeatures.append(ogr.Feature(layerDefn))

    #create new geometry and add first point
    newGeom=ogr.Geometry(type=ogr.wkbLineString)
    newGeom.AddPoint(x1,y1,z1)
    
    for pInd in range(1,nPoints):
        x1=inRoadGeom.GetX(pInd-1)
        y1=inRoadGeom.GetY(pInd-1)
        z1=inRoadGeom.GetZ(pInd-1)

        x2=inRoadGeom.GetX(pInd)
        y2=inRoadGeom.GetY(pInd)
        z2=inRoadGeom.GetZ(pInd)
                
        #Estimate direction vector
        U=np.array([x2-x1,y2-y1])
        U=U/np.linalg.norm(U)        
        diff=np.degrees(np.arccos(np.dot(U,Ulast)))
        if np.isnan(diff):
            diff=0

        if diff>180:
            diff-=180
        
        #If direction change is acceptable,
        #point is added to existing features geometry
        if abs(diff)<=splitLimit:
            newGeom.AddPoint(x2,y2,z2)
        #If to large change, a new feature and geometry is created
        else:
            outRoadFeatures[-1].SetFrom(inRoadFeature)
            outRoadFeatures[-1].SetGeometryDirectly(newGeom)
            outRoadFeatures.append(ogr.Feature(layerDefn))
            newGeom=ogr.Geometry(type=ogr.wkbLineString)
            newGeom.AddPoint(x1,y1,z1)
            newGeom.AddPoint(x2,y2,z2)
            Ulast=U

        #if last point, set geometry
        if pInd==nPoints-1:
            outRoadFeatures[-1].SetFrom(inRoadFeature)                
            outRoadFeatures[-1].SetGeometryDirectly(newGeom)

        
    return outRoadFeatures
        
    

def bheight2sect(hgt1,hgt2,angle1):
    """
    Returns string list of buildingsheights in 12 sectors relative to north
    @param: building height on side 1
    @param: building height in side 2
    @param: angle relative north to buildings on side 1
    """
    
    nSec=12

    if hgt1 is None or hgt2 is None:
        return None
        
    elif hgt1==hgt2:
        sectHgt=[hgt1]*nSec
        sectStr=" ".join(map(str,sectHgt))
        return sectStr

    elif angle1 is None:
        return None
      
    sectHgt=[None]*nSec
    secSize=360/nSec
    sects=range(0,360,secSize)
    i=0
    indicesSide1=[]

    if hgt1 is None or hgt2 is None:
        return None
    
    i=int(angle1/secSize)
    rest=angle1/float(secSize)-i
    if rest>0.5: 
        i+=1
    if i>=nSec:
        i-=nSec
                    
    for j in range(3):
        if i+j<=11:
            sectHgt[i+j]=hgt1
            indicesSide1.append(i+j)
        else:
            sectHgt[i+j-nSec]=hgt1
            indicesSide1.append(i+j-nSec)

    for j in range(-1,-4,-1):
        if i+j>=0:
            sectHgt[i+j]=hgt1
            indicesSide1.append(i+j)
        else:
            sectHgt[i+j+nSec]=hgt1
            indicesSide1.append(i+j+nSec)

    for ind in range(nSec):
        if ind not in indicesSide1:
            sectHgt[ind]=hgt2

    sectStr=" ".join(map(str,sectHgt))
    return sectStr
 
    

def plotSegments(ax,segments,color="black",style="-",width=0.5):
    for seg in segments:
        X=[seg.P0[0],seg.P1[0]]
        Y=[seg.P0[1],seg.P1[1]]
        ax.plot(X,Y,color=color,linestyle=style,linewidth=width)
        
def getBuildingSegments(building):
    """
    Return list of segments in building contour
    @param building: ogr geometry reference for building contour
    """
    nPoints=building.GetPointCount()
    segments=[]
    for i in range(1,nPoints):
        x1=building.GetX(i-1)
        y1=building.GetY(i-1)
        z1=building.GetZ(i-1)
        x2=building.GetX(i)
        y2=building.GetY(i)
        z2=building.GetZ(i)
        P1=np.array([x1,y1,z1])
        P2=np.array([x2,y2,z2])
        segments.append(Segment(P1,P2))
    return segments


def getIntersectingFacade(ax,cs,segments,plot=False):
    """
    Get nearest instersecting facade for cs
    @param cs: cross section expressed as a Segment instance
    """
    minDist=MAXDIST
    nearestIntersection=None
    log.debug("Found %i segments" %len(segments))
    for segInd,seg in enumerate(segments):
        if plot and segInd==9:
            ax.text(0.5*(seg.P0[0]+seg.P1[0]),
                    0.5*(seg.P0[1]+seg.P1[1]),
                    '%i' %segInd, fontsize=10)
        Pintersect,Poverlap = intersect2D_2Segments(cs,seg)
        if Pintersect is not None:
            #add segement z-level to intersection
            Pintersect=np.append(Pintersect,seg.P1[2])
            #evaluate horizontal distance    
            dist=np.linalg.norm(Pintersect[:-1] - cs.P0)
            if dist<=minDist:
                minDist=dist
                nearestIntersection=Pintersect

    return (minDist,nearestIntersection)
    

class SpatialIndex:
    def __init__(self,xmin,ymin,nrows,ncols,cellsize):
        """
        @param xmin: min x of extent
        @param ymin: min y of extent
        @param ncols: number of index columns
        @param nrows: number of index rows
        @param cellsize: resolution of spatial index
        """
        self.xmin=xmin
        self.xmax=xmin+ncols*cellsize
        self.ymin=ymin
        self.ymax=ymin+nrows*cellsize
        self.nrows=nrows
        self.ncols=ncols
        self.cellsize=cellsize

        self.ind={}

    def inside(self,road):
        Pll,Pur=road.getBoundingBox()
        try:
            row,col=self.getInd(Pll[0]-100,Pll[1]-100)
            row,col=self.getInd(Pur[0]+100,Pur[1]+100)
        except ValueError:
            return False
        else:
            return True

    def plotCell(self,ax,row,col,color="grey",style="--",width=0.5):
        #plot index grid cell
        #indexing clockwise from lower left corner

        x1=x2=x5=self.xmin+col*self.cellsize
        x3=x4=self.xmin+(col+1)*self.cellsize
        y1=y4=y5=self.ymax-(row+1)*self.cellsize
        y2=y3=self.ymax-row*self.cellsize
        S1=Segment(np.array([x1,y1]),(np.array([x2,y2])))
        S2=Segment(np.array([x2,y2]),(np.array([x3,y3])))
        S3=Segment(np.array([x3,y3]),(np.array([x4,y4])))
        S4=Segment(np.array([x4,y4]),(np.array([x5,y5])))
        gridSegments=[S1,S2,S3,S4]
        plotSegments(ax,gridSegments,color=color,style=style,width=width)

        
    def plot(self,ax):
        gridSegments=[]
        buildSegments=[]        
        for row in self.ind:
            for col in self.ind[row]:
                buildSegments+=self.ind[row][col]
                self.plotCell(ax,row,col)
                ax.text(self.xmin+(col+0.5)*self.cellsize,
                        self.ymax-(row+0.5)*self.cellsize,
                        '%i' %len(self.ind[row][col]), fontsize=10)

        plotSegments(ax,buildSegments,color="grey",style="-",width=1.0)
        
                
                
    def getBuildingSegments(self,x,y):
        """
        Building segments from index cells close to given point are returned
        Also segments in neighbouring index cells are included
        @param x: x coord of point to find segments for
        @param y: y coord of point to find segments for
        """
        row,col=self.getInd(x,y)
        try:
            segments=self.ind[row][col][:]
        except KeyError:
            segments=[]
            
        if col>0:
            try:
                segments+=self.ind[row][col-1]
            except KeyError:
                pass
        if col<self.ncols-1:
            try:
                segments+=self.ind[row][col+1]
            except KeyError:
                pass
        if row>0:
            try:
                segments+=self.ind[row-1][col]
            except KeyError:
                pass
            if col>0:
                try:
                    segments+=self.ind[row-1][col-1]
                except KeyError:
                    pass
            if col<self.ncols-1:
                try:
                    segments+=self.ind[row-1][col+1]
                except KeyError:
                    pass
                
        if row<self.nrows-1:
            try:
                segments+=self.ind[row+1][col]
            except KeyError:
                pass
            if col>0:
                try:
                    segments+=self.ind[row+1][col-1]
                except KeyError:
                    pass
            if col<self.ncols-1:
                try:
                    segments+=self.ind[row+1][col+1]
                except KeyError:
                    pass

        return segments    

    def getInd(self,x,y):
        """
        Function to estimate row and col of spatial index
        @param x: x coordinate to retrieve index for
        @param y: y coordinate to retrieve index for
        """
        if x<self.xmin or y<self.ymin or x > self.xmax or y > self.ymax:    
            raise ValueError("Outside extent")
        col=int((x-self.xmin)/self.cellsize)
        row=self.nrows-int(ceil((y-self.ymin)/self.cellsize))
        #If on boundary of extent
        if row == self.nrows:
            row -=1
        if col == self.ncols:
            col -=1            
        return (row,col)

    
    def indexBuildingContours(self,buildingFile):
        """
        Read building contours from shape file and create a spatial indexing.
        The indexing is implemented as a nested dict
        cols={row1,row2,row3}, where row1..n = {cell1,cell2,...,celln},
        and cell1...n = [contour_1,contour_2,...,contour_n]
        Extent of index is set to extent of the used DEM (topo)
        @param buildingFile: shapefile with 3D building contours
        """
        #open building contours shape-file
        log.info("Reading building contours")
        #Init reading of shape-files using OGR
        shapeDriver = ogr.GetDriverByName('ESRI Shapefile')
        buildings = shapeDriver.Open(buildingFile, update=0)
        if buildings is None:
            log.error("Could not open file with building contours")
            return 1
        buildingsLayer = buildings.GetLayer()
    
        self.ind={}

        #set up progress indicator
        pg=ProgressBar(buildingsLayer.GetFeatureCount(),sys.stdout)

        contourInd=0    #Counter for building contour
        buildHeights=[] #List to store building height for plotting
        noGeom=0        #Counter for invalid features without geometry reference

        #Read and process building contours
        contour = buildingsLayer.GetNextFeature() #get first feature
        while contour:      
            if log.level==0:
                pg.update(contourInd)
            log.debug("Contour %i" %contourInd)

            contourGeom = contour.GetGeometryRef()

            if contourGeom is None or contourGeom.GetPointCount()==0:
                noGeom+=1
            else:
                for i in range(1,contourGeom.GetPointCount()):
                    x1=contourGeom.GetX(i-1)
                    y1=contourGeom.GetY(i-1)                
                    z1=contourGeom.GetZ(i-1)

                    x2=contourGeom.GetX(i)
                    y2=contourGeom.GetY(i)
                    z2=contourGeom.GetZ(i)
                    
                    P1=np.array([x1,y1,z1])
                    P2=np.array([x2,y2,z2])
                    try:
                        row,col=self.getInd((x2+x1)/2.0,(y1+y2)/2.0)
                    except ValueError: #outside extent
                        row=col=None
                    if row is not None and col is not None:
                        if row not in self.ind:
                            self.ind[row]={}
                            self.ind[row][col]=[]
                        elif col not in self.ind[row]:
                            self.ind[row][col]=[]

                        self.ind[row][col].append(Segment(P1,P2))

            #Should the contour be destroyed here?
            contour.Destroy()
            contour = buildingsLayer.GetNextFeature()            
            contourInd+=1
        #close datasource for building contours
        buildings.Destroy()
        pg.finished()
        if noGeom>0:
            log.warning("Found %i building contours without geometry" %noGeom)
        

    def getBuildingHeight(self,x,y,z,topo):
        try:
            groundLevel=topo.getVal(x,y)
        except ValueError:
            buildhgt=-9999
        else:
            buildhgt=z-groundLevel
            if buildhgt<0:
                buildhgt=0
            
        return buildhgt


    def plotBuildHgtDistr(self):
        fig = plt.figure(1)
        fig.clf()
        ax = fig.add_subplot(111)
        binLims=np.arange(-20,60,1)
        n, bins, patches = ax.hist(heights,bins=binLims, normed=0)
        plt.setp(patches, 'facecolor', 'black', 'alpha', 0.75)
        plt.title("Building height distribution")
        plt.xlabel('Height [m]')
        plt.axis([-20,60,0,max(n)*1.1])
        plt.show()

class Road:
    def __init__(self,road):
        self.nPoints=road.GetPointCount()
        self.points=[]
        for i in range(self.nPoints):
            self.points.append(np.array([road.GetX(i),road.GetY(i)]))

    def normalAngles(self):
        x1=self.points[-1][0]
        y1=self.points[0][1]
        x2=self.points[-1][0]
        y2=self.points[-1][1]
        roadSegLen=np.sqrt((x2-x1)**2 + (y2-y1)**2)
        U=self.points[-1]-self.points[0]
        roadDir=U/roadSegLen
        
        roadNormal1=np.array([-1*roadDir[1],roadDir[0]])
        roadNormal2=np.array([roadDir[1],-1*roadDir[0]])
        ang1=np.degrees(np.arctan2(roadNormal1[1],roadNormal1[0]))
        ang2=np.degrees(np.arctan2(roadNormal2[1],roadNormal2[0]))
        
        ang1-=90 #relative to north
        if ang1<0:
            ang1+=360 #Make sure angle is given in positive direction
        ang2-=90 #relative to north
        if ang2<0:
            ang2+=360 #Make sure angle is given in positive direction
        if np.isnan(ang1) or np.isnan(ang2):
            return (None,None)
        return (ang1,ang2)

    

    def getBoundingBox(self):
        """
        Get bounding box for road 
        returns (minPoint,maxPoint)
        """
        P1=np.array([min(self.points[0][0],self.points[-1][0]),
                     min(self.points[0][1],self.points[-1][1]),])
        P2=np.array([max(self.points[0][0],self.points[-1][0]),
                     max(self.points[0][1],self.points[-1][1]),])
        return (P1,P2)

    def getSegments(self):
        segments=[]
        for i in range(1,self.nPoints):
            segments.append(Segment(self.points[i-1],self.points[i]))
        return segments
        
    def defineCrossSections(self):
        """
        Define cross-sections as two lists of segments, one for each side of the road
        returns ("segment list side1","segment list side 2")
        """
        cs1=[]
        cs2=[]

        for i in range(1,self.nPoints):
            x1=self.points[i-1][0]
            y1=self.points[i-1][1]
            x2=self.points[i][0]
            y2=self.points[i][1]

            roadSegLen=np.sqrt((x2-x1)**2 + (y2-y1)**2)

            nCS=int(np.ceil(roadSegLen/float(CSDIST)))

            segFrac=1/float(nCS)

            P0=self.points[i-1]
            P1=self.points[i]
            U=P1-P0
            roadDir=U/roadSegLen
            PM=P0+U*0.5 #middle of segment
            roadNormal1=np.array([-1*roadDir[1],roadDir[0]])
            roadNormal2=np.array([roadDir[1],-1*roadDir[0]])

            for csInd in range(0,nCS):
                Pi=P0+U*(csInd+0.5)*segFrac
                cs1.append(Segment(Pi,Pi+roadNormal1*MAXDIST))
                cs2.append(Segment(Pi,Pi+roadNormal2*MAXDIST))
                
        return (cs1,cs2)    

def main():
    #-----------Setting up and unsing option parser-----------------------
    parser=OptionParser(usage= usage, version=version)

    parser.add_option("-d","--doc",
                      action="store_true", dest="doc",
                      help="Prints more detailed documentation and exit")
                      
    parser.add_option("-l", "--loglevel",
                      action="store",dest="loglevel",default=2,
                      help="Sets the loglevel (0-3 where 3=full logging)")

    parser.add_option("-f","--format",
                      action="store",dest="format",default="ESRI Shapefile",
                      help="Format of road network")
    
    parser.add_option("--inRoads",
                      action="store",dest="inRoads",
                      help="Input road network")

    parser.add_option("--outRoads",
                      action="store",dest="outRoads",
                      help="Output road network")

    parser.add_option("--buildings",
                      action="store",dest="buildings",
                      help="Input building roof contours (3D)")
    
    parser.add_option("--topo",
                      action="store",dest="topo",
                      help="Input raster DEM")

    parser.add_option("--split",
                      action="store",dest="split",
                      help="Threshold in changed road direction (degrees) for when to split road")
            
    (options, args) = parser.parse_args()
    
    #------------Setting up logging capabilities -----------
    global log
    rootLogger=logger.RootLogger(int(options.loglevel))
    log=rootLogger.getLogger(sys.argv[0])

    #------------Process and validate options---------------
    if options.doc:
        print doc
        sys.exit()
    
    if len(args) > 0:
        parser.error("No arguments needed")

    #validate topo path
    if options.topo is None:
        log.error("No topo specified")
        return 1        
    if not path.exists(options.topo):
        log.error("Input raster does not exist")
        return 1

    #validate road input network path
#    if options.roads is None:
#        log.error("No input road network is specified")
#        return 1        
#    if not path.exists(options.roads):
#        log.error("Input road network does not exist")
#        return 1

    #validate building contour path
    if options.buildings is None:
        log.error("No input building contours are specified")
        return 1        
    if not path.exists(options.buildings):
        log.error("Input building contours does not exist")
        return 1


    log.info("Reading DEM")
    topo= readGDAL(options.topo,bandIndex=1)[0]

    #Init reading of shape-files using OGR
    shapeDriver = ogr.GetDriverByName('ESRI Shapefile')

    #Opening driver for road networks
    try:
        driver = ogr.GetDriverByName(options.format)
    except:
        log.error("Invalid format for road networks, check ogr documentation")
        return 1

    if options.split is None:
        log.info("Do not split roads")
        splitLimit=None
    else:
        splitLimit=float(options.split)
        log.info("Split roads that change direction"+
                 " more than %f" %splitLimit)

        
    #extract extent from topo raster
    xmin=topo.xll
    ymin=topo.yll
    xmax=topo.xur()
    ymax=topo.yur()

    #Calculate dimensions of spatial index
    ncols=int((xmax-xmin)/CELLSIZE)
    nrows=int((ymax-ymin)/CELLSIZE)
    
    #Init spatial index of building contours
    spatInd=SpatialIndex(xmin,ymin,nrows,ncols,CELLSIZE)

    log.info("Reading and indexing building contours")
    #Read buildings and store using spatial indices
    spatInd.indexBuildingContours(options.buildings)    

    #open road network shape-file
    log.info("Reading road network")
    inRoadFile = driver.Open(options.inRoads, update=0)

    if path.exists(options.outRoads):
        #outFilePath, filename = path.split(path.abspath(options.outRoads))
        name, ext = path.splitext(options.outRoads)
        for ext in [".shp",".shx",".dbf",".prj",".sbn",".qpj",".sbx"]:
            filename=name +ext
            if path.exists(filename):
                os.remove(filename)
        
    outRoadFile = driver.CreateDataSource(options.outRoads)
    if inRoadFile is None:
        log.error("Could not open file with input road network")
        return 1

    if outRoadFile is None:
        log.error("Could not open file with input road network")
        return 1

    
        
    #Get layer definition and first feature of input road network
    inRoadLayer = inRoadFile.GetLayer()        
    inRoadLayerDefn = inRoadLayer.GetLayerDefn()

    outRoadLayer=outRoadFile.CreateLayer("first_layer",geom_type=inRoadLayer.GetGeomType())
    for fieldInd  in range(inRoadLayerDefn.GetFieldCount()):
        fieldDefn = inRoadLayerDefn.GetFieldDefn(fieldInd)
        outRoadLayer.CreateField(fieldDefn)

    outRoadLayerDefn=outRoadLayer.GetLayerDefn()
    fieldNames= [outRoadLayerDefn.GetFieldDefn(i).GetName()
                 for i in range(outRoadLayerDefn.GetFieldCount())]

                                     
    log.info("Adding attributes to road feature (if missing)")
    #Add attributes for street canyon geometry
    if "BHGT1" not in fieldNames:
        outRoadLayer.CreateField(ogr.FieldDefn("BHGT1",ogr.OFTInteger))
    if "BHGT2" not in fieldNames:
        outRoadLayer.CreateField(ogr.FieldDefn("BHGT2",ogr.OFTInteger))
    if "BHGT1W" not in fieldNames:
        outRoadLayer.CreateField(ogr.FieldDefn("BHGT1W",ogr.OFTInteger))
    if "BHGT2W" not in fieldNames:
        outRoadLayer.CreateField(ogr.FieldDefn("BHGT2W",ogr.OFTInteger))
    if "BANG1" not in fieldNames:
        outRoadLayer.CreateField(ogr.FieldDefn("BANG1",ogr.OFTInteger))
    if "BANG2" not in fieldNames:
        outRoadLayer.CreateField(ogr.FieldDefn("BANG2",ogr.OFTInteger))
    if "BDIST" not in fieldNames:
        outRoadLayer.CreateField(ogr.FieldDefn("BDIST",ogr.OFTInteger))
    if "BSECT" not in fieldNames:
        fieldDefn= ogr.FieldDefn("BSECT",ogr.OFTString)
        fieldDefn.SetWidth(40)
        outRoadLayer.CreateField(fieldDefn)
    if "BSECTW" not in fieldNames:
        fieldDefn= ogr.FieldDefn("BSECTW",ogr.OFTString)
        fieldDefn.SetWidth(40)
        outRoadLayer.CreateField(fieldDefn)

        
    fig1=plt.figure(1)
    ax1=plt.subplot(111)
    ax1.axis('equal')
    if PLOTIND>0:
        spatInd.plot(ax1)
        
    roadInd=0
    noGeom=0
    nsplit=0
    inRoadFeature = inRoadLayer.GetNextFeature() #get first road feature
    #Loop over all roads
    log.info("Finding nearest facades, setting heights...")
    pg=ProgressBar(inRoadLayer.GetFeatureCount(),sys.stdout)
    while inRoadFeature:

        #name=inRoadFeature.GetFieldAsString("Namn_130")
        #adt=inRoadFeature.GetFieldAsString("ADT_f_117")        
        #roadLen=int(float(inRoadFeature.GetFieldAsString("ExtLen")))
        #print "\n"+name[:6]
        #print roadLen
        #if name[:6]=="Landsv" and roadLen ==121:
        #    global PLOTIND
        #    PLOTIND=roadInd
            
        pg.update(roadInd)

        outRoadFeatures=splitRoad(inRoadFeature,splitLimit,outRoadLayer.GetLayerDefn())
        if len(outRoadFeatures)>1:            
            log.debug("Raod split into %s parts" %len(outRoadFeatures))
        nsplit+=len(outRoadFeatures)-1
        for outRoadFeature in outRoadFeatures:                    
            intersections=[]
            outRoadGeom=outRoadFeature.GetGeometryRef()
            road=Road(outRoadGeom)
            if outRoadGeom is None or outRoadGeom.GetPointCount()==0 or not spatInd.inside(road):
                noGeom+=1
                maxHeight1=None
                maxHeight2=None
                avgDist=None
                bAngle1=None
                bAngle2=None
                avgHeight1=None
                avgHeight2=None
            else:
                sumHeight1=0
                sumHeight2=0
                maxHeight1=0
                maxHeight2=0
                sumDist=0
                #Define crossections along the road,
                #defined by start and endpoints at both side of the road
                cs1List,cs2List=road.defineCrossSections()
                nCS=len(cs1List)
                log.debug("Defined %i cross sections" %nCS)
                #Check intersections with building contours for all cross-sections
                for csInd in range(nCS):
                    cs1=cs1List[csInd]
                    cs2=cs2List[csInd]
                    cs1MidPoint=cs1.P0+0.5*(cs1.P1-cs1.P0)
                    buildingSegments=spatInd.getBuildingSegments(
                        cs1MidPoint[0],cs1MidPoint[1])

                    log.debug("Calculating intersection")
                    if PLOTIND==roadInd and CSIND==csInd:
                        dist1,Pint1=getIntersectingFacade(ax1,cs1,buildingSegments,True)
                    else:
                        dist1,Pint1=getIntersectingFacade(ax1,cs1,buildingSegments,False)
                    if Pint1 is None:
                        log.debug("No intersection on side 1")
                        height1=0
                        dist1=MAXDIST
                    else:
                        log.debug("Intersection1 in (%f,%f,%f)" %(Pint1[0],Pint1[1],Pint1[2]))
                        height1= spatInd.getBuildingHeight(
                            Pint1[0],Pint1[1],Pint1[2],topo)+HEIGHTCORR
                        intersections.append(Pint1[:2])

                    if PLOTIND==roadInd and csInd == CSIND:
                        plotSegments(ax1,buildingSegments,color='red',width=2.0)
                        row,col=spatInd.getInd(cs1MidPoint[0],cs1MidPoint[1])             
                        spatInd.plotCell(ax1,row,col,color="purple",width=2.0)
                        plotSegments(ax1,[cs1List[csInd]],color="pink",style="-",width=1.0)
                        plt.draw()



                    cs2MidPoint=cs2.P0+0.5*(cs2.P1-cs2.P0)
                    buildingSegments=spatInd.getBuildingSegments(
                        cs2MidPoint[0],cs2MidPoint[1])

                    if PLOTIND==roadInd and csInd == CSIND:
                        plotSegments(ax1,buildingSegments,color='red',width=2.0)
                        row,col=spatInd.getInd(cs2MidPoint[0],cs2MidPoint[1])
                        spatInd.plotCell(ax1,row,col,color="brown",width=2.0)
                        plotSegments(ax1,[cs2List[csInd]],color="red",style="-",width=1.0)
                        plt.draw()

                    log.debug("Calculating intersection")
                    if PLOTIND==roadInd and CSIND==csInd:
                        dist2,Pint2=getIntersectingFacade(ax1,cs2,buildingSegments,True)
                    else:
                        dist2,Pint2=getIntersectingFacade(ax1,cs2,buildingSegments,False)



                    if Pint2 is None:
                        log.debug("No intersection on side 2")
                        height2=0
                    else:
                        log.debug("Intersection2 in (%f,%f,%f)" %(Pint2[0],Pint2[1],Pint2[2]))
                        height2= spatInd.getBuildingHeight(
                            Pint2[0],Pint2[1],Pint2[2],topo)+HEIGHTCORR
                        intersections.append(Pint2[:2])

                    sumHeight1+=height1
                    sumHeight2+=height2
                    sumDist+=dist1+dist2
                    maxHeight1=int(max(height1,maxHeight1))
                    maxHeight2=int(max(height2,maxHeight2))
                    if PLOTIND==roadInd and CSIND==csInd:
                        if Pint1 is not None:
                            ax1.text(Pint1[0],Pint1[1],"Distance=%f" %dist1)
                        if Pint2 is not None:
                            ax1.text(Pint2[0],Pint2[1],"Distance=%f" %dist2)

                avgHeight1=int(sumHeight1/float(nCS))
                avgHeight2=int(sumHeight2/float(nCS))
                #averaging over both sides of street
                #distance refers to between facades on opposite sides
                avgDist=int(round(sumDist/float(nCS)))
                bAngle1,bAngle2=road.normalAngles()


                if PLOTIND>0:
                    plotSegments(ax1,road.getSegments(),color='grey',width=0.3)
                if PLOTIND==roadInd:
                    plotSegments(ax1,road.getSegments(),color='black',width=2.0)
                    plotSegments(ax1,cs1List,color="green",style="--",width=0.5)
                    plotSegments(ax1,cs2List,color="green",style="--",width=0.5)

                    X=[intersect[0] for intersect in intersections]
                    Y=[intersect[1] for intersect in intersections]
                    if len(X)>0:
                        ax1.plot(X,Y,"*")
                    plt.title("Road %i, cross-section %i" %(PLOTIND,CSIND))
                    plt.draw()

            #building height as list of sectors
            bsect=bheight2sect(avgHeight1,avgHeight2,bAngle1)
            bsectw=bheight2sect(maxHeight1,maxHeight2,bAngle1)

            outRoadFeature.SetField("BSECT",bsect)
            outRoadFeature.SetField("BSECTW",bsectw)
            outRoadFeature.SetField("BHGT1",avgHeight1)
            outRoadFeature.SetField("BHGT2",avgHeight2)
            outRoadFeature.SetField("BHGT1W",maxHeight1)
            outRoadFeature.SetField("BHGT2W",maxHeight2)
            outRoadFeature.SetField("BANG1",bAngle1)
            outRoadFeature.SetField("BANG2",bAngle2)
            outRoadFeature.SetField("BDIST",avgDist)
            
            outRoadLayer.CreateFeature(outRoadFeature)
            outRoadFeature.Destroy()
        inRoadFeature.Destroy()
        inRoadFeature = inRoadLayer.GetNextFeature()            
        roadInd+=1

    inRoads=inRoadLayer.GetFeatureCount()
    outRoads=outRoadLayer.GetFeatureCount()
    #close datasource for building contours
    inRoadFile.Destroy()
    outRoadFile.Destroy()
    pg.finished()
    if PLOTIND>0:
        plt.show()

    log.info("Read %i roads, wrote %i roads (created %i by splitting)" %(
        inRoads,outRoads,nsplit))

    if noGeom>0:
        log.warning("Found %i roads without geometry" %noGeom)

    log.info("Finished")

if __name__=="__main__":
    sys.exit(main())
    # try:
    #     main()
    # except:
    #     import pdb, sys
    #     e, m, tb = sys.exc_info()
    #     pdb.post_mortem(tb)
