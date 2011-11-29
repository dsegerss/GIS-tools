import PyRaster
import os,sys,numpy

#maskFile="//Winfs/data/proj/Luftkvalitet/Masker/alla_masker.asc"
maskFile="//Winfs/data/prod/Smeddata/Geofordelning/2007/regionRasters/Cat_01_2007_kommungrupper.asc"
emisFile="//Winfs/data/prod/Smeddata/Geofordelning/2007/kvalitetskontroller/ElVarmeverk_u/ElVarmeverk_uCO2tonperyear2007.asc"
flippa=0
#--------------------
outputDir=os.path.dirname(emisFile)
summaryFile="sumByMask_"+os.path.basename(emisFile)
summaryPath=os.path.join(outputDir,summaryFile)
output=open(summaryPath,'w')

maskRast=PyRaster.raster()
maskRast.readAscii(maskFile)
maskRast.nodataToZero()

emisRast=PyRaster.raster()
emisRast.readAscii(emisFile)
emisRast.nodataToZero()
if flippa:
    emisRast.flip()
print "Finding all mask codes"
maskCodes=maskRast.unique()

if maskRast.nodata in maskCodes:
    maskCodes.remove(maskRast.nodata)
if 0.0 in maskCodes:
    maskCodes.remove(0.0)

ortRast=PyRaster.raster()
output.write("Summary by mask code for file "+os.path.basename(emisFile)+"\n")
output.write("code\tsum\n")

#print "Doing summary for code:"
#for code in maskCodes:
#    print code
    
for code in maskCodes:
    print "summing for mask code: ", code
    ortRast.data=numpy.where(maskRast.data==code,emisRast.data,0.0)
    sum=ortRast.sum()
    output.write(str(code)+"\t"+str(sum)+"\n")
output.close()
print "finished!"
    
    