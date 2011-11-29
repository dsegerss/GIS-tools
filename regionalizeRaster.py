""" This script can be used to weight a key raster by regional statistical data.
The script creates an output raster sum of the emissions within a region equals the regional total
The emissions within a given each region is distributed using a specified distribution key
The script is hard coded for the SMED-raster extents and resolution
David Segersson, 080212
"""
import os,future,numpy,sys,pyEmissions,fileInput,PyRaster

#-----------------------------------------
#If the parameter uniformDistribution is set to True no keyRaster is needed. The data is distributed uniformly within the regions
#uniformDistribution=False 
#keyRasterPath="Y:/PROJEKT/2010/Fredrikstad_Norge/Konvertering_griddade_emissioner/road_buffer_zone.asc"
#regionRasterPath="Y:/PROJEKT/2010/Fredrikstad_Norge/Konvertering_griddade_emissioner/grunnkretser.asc"
#regionalTotalsPath="Y:/PROJEKT/2010/Fredrikstad_Norge/Konvertering_griddade_emissioner/PM10_boliger.txt"
#resultRasterPath="Y:/PROJEKT/2010/Fredrikstad_Norge/Konvertering_griddade_emissioner/residential_PM10.asc"
#substance="PM10" 
#distributeEmis=True
uniformDistribution=False 
keyRasterPath="//Winfs/data/prod/Smeddata/Geofordelning/keys/industrimark_v1.asc"
regionRasterPath="//Winfs/data/prod/Smeddata/Geofordelning/geodata/ascii/ln07_sr99tm.asc"
regionalTotalsPath="//Winfs/data/prod/Smeddata/Geofordelning/originalstatistik/skogsbruk/avverkning_2006_08.txt"
resultRasterPath="//Winfs/data/prod/Smeddata/Geofordelning/keys/avverkning_v5.asc"
substance="all" 
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

if not uniformDistribution:
    if os.path.exists(keyRasterPath):
        print "Reading key raster"
        keyRast.readAscii(keyRasterPath)
    else:
        sys.exit("keyRaster does not exist")
else:
    keyRast.nodataToZero()
    keyRast=keyRast+1
        
    
regionalTotals=pyEmissions.RegTotalEmis("dummy")
if os.path.exists(regionalTotalsPath):
    print "Reading regional totals"
    regionalTotals.readAscii(regionalTotalsPath)
else:
    sys.exit("Regional totals file does not exist")
    
print "Finished reading data!"

keyRast.nodataToZero()
regionRast.nodataToZero()

print "Ensuring consistency and completeness in codes in region raster and regional totals"
#Skapar en lista med de koder som finns representerade i regionsrastret
codesInRegionRast=regionRast.unique()
if regionRast.nodata in codesInRegionRast:
    codesInRegionRast.remove(regionRast.nodata)
if 0.0 in codesInRegionRast:
    codesInRegionRast.remove(0.0)

codesInRegionalTotals=regionalTotals.getIDs()
if "National" in codesInRegionalTotals:
    codesInRegionalTotals.remove("National")

totFractionList=[]
codeNumList=[]
for code in codesInRegionalTotals:
    totalEmis=regionalTotals.totalDict[code]
    emission=totalEmis.emisDict[substance]
    if emission==-9999:
        emission=0
    emissionFraction=emission/regionalTotals.regSum(substance)
    totFractionList.append((int(code),emissionFraction))
    #creates a copy of codesInRegionalTotals with numerical values instead of strings
    codeNumList.append(int(code))
    #Assuring that the regional total codes are present in the regional raster
    if int(code) not in codesInRegionRast:
        sys.exit("Input Error: code:"+ str(int(code))+"in regional totals is not represented in regional raster") 

#Assuring that the regional raster IDs are present in the regional totals                
for code in codesInRegionRast:
    if code not in codeNumList:
        sys.exit("Input Error: ID:"+str(code)+" in region raster is not represented in regional totals!")

#Creates a raster with all cells within a region given the regions fraction of the national total
totFractionRast=regionRast.replace(totFractionList)
totFractionRast.nodataToZero()
#For all regions, calculate the regional sum, set the key values to key/regSum
resKeyArray=numpy.zeros(keyRast.data.shape)
maskArray=numpy.zeros(regionRast.data.shape)
print "Regionalizing key raster"
for code in codesInRegionRast:
    totalEmis=regionalTotals.totalDict[str(int(code))]
    emission=totalEmis.emisDict[substance]
    print "Regionalizing region: " + str(code)+" with the regional total: ",emission     
    regDistThread=future.Future(distributeRegion,[regionRast.data,keyRast.data,code,resKeyArray,emission])         
    resKeyArray=regDistThread()
keyRast.data=resKeyArray

resultRast=(totFractionRast*keyRast)*regionalTotals.regSum(substance)        
resultRast.write(resultRasterPath)
print "Result with sum: ",resultRast.sum()," written to disk!"

