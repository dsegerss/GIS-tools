""" This script can be used to weight a set of emission rasters by regional statistical data. The resulting emission rasters will have the same total sum as the original.
The script creates an output raster sum of the emissions within a region equals the regional total
The emissions within a given each region is distributed using a specified distribution key
The script is hard coded for the SMED-raster extents and resolution
David Segersson, 080212
"""
import os,future,numpy,sys,pyEmissions,fileInput,PyRaster
from os import path
#-----------------------------------------
#If the parameter uniformDistribution is set to True no keyRaster is needed. The data is distributed uniformly within the regions
uniformDistribution=False 
keyRasterPath="//Winfs/data/prod/Smeddata/Geofordelning/2007/Energi/"
keyRastPrefix="Energi"
keyRastPostfix="tonperyear2007.asc"

substances=["DIOX","Hg",
           "N2O","NH3","Ni","NMVOC","NOx",
           "PAH","Pb","PM10","PM25","Se",
           "SO2","TSP","Zn"] 
#klara: "As","BAP","Cd","CH4","CO2", "CO","Cr","Cu"

regionalTotalsPath="//Winfs/data/prod/Smeddata/Geofordelning/2007/regionalTotals/Cat_01_excl_1A1b_and_1A2fi_lan_2007.txt"
regionRasterPath="//Winfs/data/prod/Smeddata/Geofordelning/geodata/ascii/ln07.asc"
resultRasterPath="//Winfs/data/prod/Smeddata/Geofordelning/2007/Energi"
#--------------------------------

def distributeRegion(argList):
    regionArray=argList[0]
    data=argList[1]
    code=argList[2]
    result=argList[3]
    regionalTotal=argList[4]

    mask=numpy.array(regionArray==code)
    regKey=mask*data
    regSum=numpy.sum(numpy.sum(regKey))
    #These two are put here to ensure that all threads are executed equally fast.
    onesRast=mask*(numpy.ones(mask.shape)*1.0)
    regMaskSum=numpy.sum(numpy.sum(onesRast))
    if regSum>0:
        result=result+regKey/regSum
    elif regionalTotal>0:
        print "Warning, for code: ",code," the distribution key is zero but the regional total is: ",regionalTotal
        print "The regional total is distributed uniformly over the whole region!"        
        result=result+onesRast/regMaskSum      
    return result

keyRast= PyRaster.raster() # skapar ett initialt tomt raster
regionRast= PyRaster.raster()# skapar ett initialt tomt raster

#Laser in indata efter att ha kontrollerat att filerna finns
if os.path.exists(regionRasterPath):
    print "Reading region raster"
    regionRast.readAscii(regionRasterPath)  
else:
    sys.exit("regionRaster does not exist")

regionalTotals=pyEmissions.RegTotalEmis("dummy")
if os.path.exists(regionalTotalsPath):
    print "Reading regional totals"
    regionalTotals.readAscii(regionalTotalsPath)
else:
    sys.exit("Regional totals file does not exist")


#Skapar en lista med de koder som finns representerade i regionsrastret
codesInRegionRast=regionRast.unique()
if regionRast.nodata in codesInRegionRast:
    codesInRegionRast.remove(regionRast.nodata)
if 0.0 in codesInRegionRast:
    codesInRegionRast.remove(0.0)

codesInRegionalTotals=regionalTotals.getIDs()
if "National" in codesInRegionalTotals:
    codesInRegionalTotals.remove("National")


for subst in substances:
    totFractionList=[]
    codeNumList=[]
    for code in codesInRegionalTotals:
        totalEmis=regionalTotals.totalDict[code]
        emission=totalEmis.emisDict[subst]
        if emission==-9999:
            emission=0
        emissionFraction=emission/regionalTotals.regSum(subst)
        totFractionList.append((int(code),emissionFraction))
        #creates a copy of codesInRegionalTotals with numercal values instead of strings
        codeNumList.append(int(code))
        #Assuring that the regional total codes are present in the regional raster
        if int(code) not in codesInRegionRast:
            sys.exit("Input Error: code:"+ str(int(code))+"in regional totals is not represented in regional raster") 

    #Assuring that the regional raster IDs are present in the regional totals                
    for code in codesInRegionRast:
        if code not in codeNumList:
            sys.exit("Input Error: ID:"+str(code)+" in region raster is not represented in regional totals!")
        
    keyRastName=path.join(keyRasterPath,keyRastPrefix+subst+keyRastPostfix)
    if os.path.exists(keyRasterPath):
        print "Reading key raster for substance: ", subst
        keyRast.readAscii(keyRastName)
    else:
        sys.exit("keyRaster does not exist")
    
    keyRast.nodataToZero()
    regionRast.nodataToZero()
    total=keyRast.sum()

    #Creates a raster with all cells within a region given the regions fraction of the national total
    totFractionRast=regionRast.replace(totFractionList)
    totFractionRast.nodataToZero()
    #For all regions, calculate the regional sum, set the key values to key/regSum
    resKeyArray=numpy.zeros(keyRast.data.shape)
    maskArray=numpy.zeros(regionRast.data.shape)
    print "Regionalizing key raster"
    print "Original sum: ",total
    for code in codesInRegionRast:
        totalEmis=regionalTotals.totalDict[str(int(code))]
        emission=totalEmis.emisDict[subst]
        print "Regionalizing region: " + str(code)+" with the regional total: ",emission     
        regDistThread=future.Future(distributeRegion,[regionRast.data,keyRast.data,code,resKeyArray,emission])         
        resKeyArray=regDistThread()
    keyRast.data=resKeyArray
    
    
    resultRast=(totFractionRast*keyRast)*regionalTotals.regSum(subst)
    
    #Assure that the total is unchanged
    resultRast=(resultRast/resultRast.sum())*total
    
    resultRasterName=path.join(resultRasterPath,"redist_"+keyRastPrefix+subst+keyRastPostfix)        
    resultRast.write(resultRasterName)
    print "wrote result for substance: ",subst

